#' Estimate publication bias-corrected meta-analysis
#'
#' For a chosen ratio of publication probabilities, \code{eta}, estimates a
#' publication bias-corrected pooled point estimate and confidence interval per
#' Mathur & VanderWeele (2020). Model options include fixed-effects (a.k.a.
#' "common-effect"), robust independent, and robust clustered specifications.
#' @export
#'
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances for the point estimates.
#' @param sei A vector of estimated standard errors for the point estimates
#'   (only relevant when not using \code{vi}).
#' @param eta The number of times more likely an affirmative study is to be
#'   published than a nonaffirmative study; see Details.
#' @param clustervar A character, factor, or numeric vector with the same length
#'   as yi. Unique values should indicate unique clusters of point estimates. By
#'   default, assumes all point estimates are independent.
#' @param model "fixed" for fixed-effects (a.k.a. "common-effect") or "robust"
#'   for robust random-effects.
#' @param selection_tails 1 (for one-tailed selection, recommended for its
#'   conservatism) or 2 (for two-tailed selection).
#' @param favor_positive \code{TRUE} if publication bias is assumed to favor
#'   positive estimates; \code{FALSE} if assumed to favor negative estimates;
#'   see Details.
#' @param alpha_select Alpha-level at which publication probability is assumed
#'   to change.
#' @param ci_level Confidence interval level (as proportion) for the corrected
#'   point estimate. (The alpha level for inference on the corrected point
#'   estimate will be calculated from \code{ci_level}.)
#' @param small Should inference allow for a small meta-analysis? We recommend
#'   always using \code{TRUE}.
#'
#' @details The ratio \code{eta} represents the number of times more likely
#'   affirmative studies (i.e., those with a "statistically significant" and
#'   positive estimate) are to be published than nonaffirmative studies (i.e.,
#'   those with a "nonsignificant" or negative estimate).
#'
#'   If \code{favor_positive == FALSE}, such that publication bias is assumed to
#'   favor negative rather than positive estimates, the signs of \code{yi} will
#'   be reversed prior to performing analyses. The corrected estimate will be
#'   reported based on the recoded signs rather than the original sign
#'   convention.
#'
#' @return A list with three elements, \code{values}, \code{stats} and
#'   \code{fit}. Stats is a list that contains the bias-corrected pooled point
#'   estimate (\code{estimate}) and inference on the bias-corrected estimate
#'   (\code{se}, \code{ci_lower}, \code{ci_upper}, \code{p_value}). Values is a
#'   list that contains the user's specified \code{eta}, the number of
#'   affirmative and nonaffirmative studies (\code{k_affirmative} and
#'   \code{k_nonaffirmative}), and a dataframe combining \code{yi}, \code{vi},
#'   \code{clustervar}.
#'
#' @references Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for
#'   publication bias in meta-analyses. \emph{Journal of the Royal Statistical
#'   Society, Series C.} Preprint available at https://osf.io/s9dp6/.
#'
#' @examples
#'  # calculate effect sizes from example dataset in metafor
#'  require(metafor)
#'  dat = metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
#'                        di = cneg, data = dat.bcg)
#'
#'  # first fit fixed-effects model without any bias correction
#'  # since the point estimate is negative here, we'll assume publication bias
#'  # favors negative log-RRs rather than positive ones
#'  metafor::rma( yi, vi, data = dat, method = "FE" )
#'
#'  # warmup
#'  # note that passing eta = 1 (no publication bias) yields the naive point
#'  # estimate from rma above, which makes sense
#'  pubbias_eta_corrected( yi = dat$yi,
#'                         vi = dat$vi,
#'                         eta = 1,
#'                         model = "fixed",
#'                         favor_positive = FALSE )
#'
#'  # assume a known selection ratio of 5
#'  # i.e., affirmative results are 5x more likely to be published than
#'  # nonaffirmative ones
#'  pubbias_eta_corrected( yi = dat$yi,
#'                         vi = dat$vi,
#'                         eta = 5,
#'                         favor_positive = FALSE,
#'                         model = "fixed" )
#'
#'  # same selection ratio, but now account for heterogeneity and clustering via
#'  # robust specification
#'  pubbias_eta_corrected( yi = dat$yi,
#'                         vi = dat$vi,
#'                         eta = 5,
#'                         favor_positive = FALSE,
#'                         clustervar = dat$author,
#'                         model = "robust" )
#'
#'  ##### Make sensitivity plot as in Mathur & VanderWeele (2020) #####
#'  # range of parameters to try (more dense at the very small ones)
#'  etas = c( 200, 150, 100, 50, 40, 30, 20, seq(15, 1) )
#'
#'  # compute estimate for each value of eta
#'  estimates = lapply(etas, function(e) {
#'    pubbias_eta_corrected( yi = dat$yi, vi = dat$vi, eta = e, model = "robust",
#'                           clustervar = dat$author, favor_positive = FALSE )$stats
#'  })
#'  estimates = dplyr::bind_rows(estimates)
#'  estimates$eta = etas
#'
#'  require(ggplot2)
#'  ggplot( estimates, aes( x = eta, y = estimate ) ) +
#'    geom_ribbon( aes( ymin = ci_lower, ymax = ci_upper ), fill = "gray" ) +
#'    geom_line( lwd = 1.2 ) +
#'    labs( x = bquote( eta ), y = bquote( hat(mu)[eta] ) ) +
#'    theme_classic()

pubbias_eta_corrected = function( yi,
                                  vi,
                                  sei,
                                  eta,
                                  clustervar = 1:length(yi),
                                  model,
                                  selection_tails = 1,
                                  favor_positive,
                                  alpha_select = 0.05,
                                  ci_level = 0.95,
                                  small = TRUE ) {

  # stop if eta doesn't make sense
  if ( eta < 1 ) stop( "Eta must be at least 1.")

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  # number of point estimates
  k = length(yi)

  # calculate alpha for inference on point estimate
  alpha = 1 - ci_level

  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k & model == "fixed" ) {
    warning( "Clusters exist, but will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }

  ##### Flip Estimate Signs If Needed #####
  # if favor_positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor_positive == TRUE ) {
    # keep track of whether we flipped for reporting at the end
    flipped = FALSE
    yif = yi
  } else {
    flipped = TRUE
    yif = -yi
  }

  # OLD VERSION: decides whether to flip signs based on naive meta-analysis
  # # check and flip if naive point estimate is negative
  # # do standard meta
  # m0 = rma.uni(yi, vi)
  #
  # # reverse signs if needed to have pooled point estimate > 0
  # if ( m0$b < 0 ) {
  #   # keep track so that we can flip back at the end
  #   flipped = TRUE
  #   yif = -yi
  # } else {
  #   flipped = FALSE
  #   yif = yi
  # }

  # 2-sided p-values for each study even if 1-tailed selection
  pvals = 2 * ( 1 - pnorm( abs(yif) / sqrt(vi) ) )

  # affirmative indicator based on selection tails
  if ( selection_tails == 1 ) A = (pvals < alpha_select) & (yif > 0)
  if ( selection_tails == 2 ) A = (pvals < alpha_select)

  k_affirmative = sum(A)
  k_nonaffirmative = k - sum(A)

  if ( k_affirmative == 0 | k_nonaffirmative == 0 ) {
    stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
  }

  dat = data.frame( yi, yif, vi, A, clustervar )


  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {

    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat = dat %>% dplyr::group_by(A) %>%
      dplyr::summarise( nu = sum( 1 / vi ),
                        ybar = sum( yi / vi ) )

    # components of bias-corrected estimate by affirmative status
    ybarN = strat$ybar[ strat$A == 0 ]
    ybarS = strat$ybar[ strat$A == 1 ]
    nuN = strat$nu[ strat$A == 0 ]
    nuS = strat$nu[ strat$A == 1 ]

    # corrected pooled point estimate
    est = ( eta * ybarN + ybarS ) / ( eta * nuN + nuS )

    # inference
    var = ( eta^2 * nuN + nuS ) / ( eta * nuN + nuS )^2
    se = sqrt(var)

    # z-based inference
    if ( small == FALSE ) {
      lo = est - qnorm( 1 - (alpha/2) ) * sqrt(var)
      hi = est + qnorm( 1 - (alpha/2) ) * sqrt(var)
      z =  abs( est / sqrt(var) )
      pval_est = 2 * ( 1 - pnorm( z ) )
    }

    # t-based inference
    if ( small == TRUE ) {
      df = k - 1
      lo = est - qt( 1 - (alpha/2), df = df ) * sqrt(var)
      hi = est + qt( 1 - (alpha/2), df = df ) * sqrt(var)
      t =  abs( est / sqrt(var) )
      pval_est = 2 * ( 1 - pt( t, df = df ) )
    }
  } # end fixed = TRUE

  ##### Robust Independent and Robust Clustered #####
  if ( model == "robust" ) {

    # weight for model
    weights = rep( 1, length(pvals) )
    weights[ A == FALSE ] = eta

    # initialize a dumb (unclustered and uncorrected) version of tau^2
    # which is only used for constructing weights
    meta_re = metafor::rma.uni( yi = yi,
                                vi = vi)
    t2hat_naive = meta_re$tau2

    # fit weighted robust model
    meta_robu = robumeta::robu( yi ~ 1,
                                studynum = clustervar,
                                data = dat,
                                userweights = weights / (vi + t2hat_naive),
                                var.eff.size = vi,
                                small = small )

    est = as.numeric(meta_robu$b.r)
    se = meta_robu$reg_table$SE
    lo = meta_robu$reg_table$CI.L
    hi = meta_robu$reg_table$CI.U
    pval_est = meta_robu$reg_table$prob
    eta = eta
  } # end robust = TRUE

  values = list(eta = eta,
                k = k,
                k_affirmative = k_affirmative,
                k_nonaffirmative = k_nonaffirmative,
                data = dplyr::rename(dat, affirm = .data$A))

  stats = list(estimate = est,
               se = se,
               ci_lower = lo,
               ci_upper = hi,
               p_value = pval_est)

  fit = list()
  if ( exists("meta_robu") ) {
    fit$robust = meta_robu
  }

  return(list(values = values, stats = stats, fit = fit))

}


#' @rdname pubbias_eta_corrected
#' @param selection.tails (deprecated) see selection_tails
#' @param favor.positive (deprecated) see favor_positive
#' @param alpha.select (deprecated) see alpha_select
#' @param CI.level (deprecated) see ci_level
#' @export
corrected_meta <- function( yi,
                            vi,
                            eta,
                            clustervar = 1:length(yi),
                            model,
                            selection.tails = 1,
                            favor.positive,
                            alpha.select = 0.05,
                            CI.level = 0.95,
                            small = TRUE ) {
  .Deprecated("pubbias_eta_corrected")
  pubbias_eta_corrected(yi = yi,
                        vi = vi,
                        eta = eta,
                        clustervar = clustervar,
                        model = model,
                        selection_tails = selection.tails,
                        favor_positive = favor.positive,
                        alpha_select = alpha.select,
                        ci_level = CI.level,
                        small = small)
}
