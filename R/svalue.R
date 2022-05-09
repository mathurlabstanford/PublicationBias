#' Severity of publication bias needed to "explain away" results
#'
#' Estimates the S-value, defined as the severity of publication bias (i.e., the ratio
#' by which affirmative studies are more likely to be published than nonaffirmative studies)
#' that would be required to shift the pooled point estimate or its confidence interval limit
#' to the value \code{q}.
#' @param yi A vector of point estimates to be meta-analyzed. Their signs should be coded such that publication bias is
#' assumed to favor positive, rather than negative, estimates.
#' @param vi A vector of estimated variances for the point estimates
#' @param sei A vector of estimated standard errors for the point estimates (only relevant when not using vi)
#' @param q The attenuated value to which to shift the point estimate or CI. Should be specified on the same scale as \code{yi}
#' (e.g., if \code{yi} is on the log-RR scale, then \code{q} should be as well).
#' @param clustervar A character, factor, or numeric vector with the same length as \code{yi}. Unique values should indicate
#' unique clusters of point estimates. If left unspecified, assumes studies are independent.
#' @param model \code{"fixed"} for fixed-effects (a.k.a. "common-effect") or \code{"robust"} for robust random-effects
#' @param alpha_select Alpha-level at which publication probability is assumed to change
#' @param eta_grid_hi The largest value of \code{eta} that should be included in the grid search. This argument is only needed when \code{model = "robust"}.
#' @param favor_positive \code{TRUE} if publication bias is assumed to favor positive estimates; \code{FALSE} if assumed to favor negative estimates.
#' See Details.
#' @param ci_level Confidence interval level (as a proportion) for the corrected point estimate
#' @param small Should inference allow for a small meta-analysis? We recommend using always using \code{TRUE}.
#' @param return_worst_meta Should the worst-case meta-analysis of only the nonaffirmative studies be returned?
#' @import
#' metafor
#' stats
#' robumeta
#' ggplot2
#' @importFrom
#' dplyr %>% group_by summarise
#' @export
#' @details
#' To illustrate interpretation of the S-value, if the S-value for the point estimate is 30 with \code{q=0}, this indicates that affirmative studies
#' (i.e., those with a "statistically significant" and positive estimate) would need to be 30-fold more likely to be published
#' than nonaffirmative studies (i.e., those with a "nonsignificant" or negative estimate) to attenuate the pooled point estimate to
#' \code{q}.
#'
#' If \code{favor_positive == FALSE}, such that publication bias is assumed to favor negative rather than positive estimates, the signs of \code{yi} will be reversed prior to
#' performing analyses. The returned number of affirmative and nonaffirmative studies will reflect the recoded signs, and accordingly the returned value \code{signs.recoded} will be \code{TRUE}.
#' @return
#' The function returns: the amount of publication bias required to attenuate the pooled point estimate to \code{q} (\code{sval.est}),
#' the amount of publication bias required to attenuate the confidence interval limit of the pooled point estimate to \code{q} (\code{sval.ci}),
#' the number of affirmative and nonaffirmative studies after any needed recoding of signs (\code{k.affirmative} and \code{k.nonaffirmative}),
#' and an indicator for whether the point estimates' signs were recoded (\code{signs.recoded}).
#'
#' If \code{return_worst_meta = TRUE}, also returns the worst-case meta-analysis of only the nonaffirmative studies. If \code{model = "fixed"}, the worst-case meta-analysis is fit by \code{metafor::rma.uni}. If \code{model = "robust"}, it is fit by \code{robumeta::robu}. Note that in the latter case, custom inverse-variance weights are used, which are the inverse of the sum of the study's variance and a heterogeneity estimate from a naive random-effects meta-analysis (Mathur & VanderWeele, 2020). This is done for consistency with the results of \code{pubbias_eta_corrected}, which is used to determine \code{sval.est} and \code{sval.ci}. Therefore, the worst-case meta-analysis results may differ slightly from what you would obtain if you simply fit \code{robumeta::robu} on the nonaffirmative studies with the default weights.
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. \emph{Journal of the Royal Statistical Society, Series C.} Preprint available at https://osf.io/s9dp6/.
#' @examples
#'  # calculate effect sizes from example dataset in metafor
#'  require(metafor)
#'  dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#'  ##### Fixed-Effects Specification #####
#'  # S-values and worst-case meta-analysis under fixed-effects specification
#'  svals.FE.0 = pubbias_svalue( yi = dat$yi,
#'                               vi = dat$vi,
#'                               q = 0,
#'                               favor_positive = FALSE,
#'                               model = "fixed" )
#'
#'  # publication bias required to shift point estimate to 0
#'  svals.FE.0$sval.est
#'
#'  # and to shift CI to include 0
#'  svals.FE.0$sval.ci
#'
#'  # now try shifting to a nonzero value (RR = 0.90)
#'  svals.FE.q = pubbias_svalue( yi = dat$yi,
#'                               vi = dat$vi,
#'                               q = log(.9),
#'                               favor_positive = FALSE,
#'                               model = "fixed" )
#'
#'  # publication bias required to shift point estimate to RR = 0.90
#'  svals.FE.q$sval.est
#'
#'  # and to shift CI to RR = 0.90
#'  svals.FE.q$sval.ci
#'
#'  ##### Robust Clustered Specification #####
#'  pubbias_svalue( yi = dat$yi,
#'                  vi = dat$vi,
#'                  q = 0,
#'                  favor_positive = FALSE,
#'                  model = "robust" )

pubbias_svalue = function( yi,
                           vi,
                           sei,
                           q,
                           clustervar = 1:length(yi),
                           model,
                           alpha_select = 0.05,
                           eta_grid_hi = 200,
                           favor_positive,
                           ci_level = 0.95,
                           small = TRUE,
                           return_worst_meta = FALSE ) {

  # stop if eta doesn't make sense
  if ( eta_grid_hi < 1 ) stop( "eta_grid_hi must be at least 1.")

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  # number of point estimates
  k.studies = length(yi)

  alpha = 1 - ci_level

  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k.studies & model == "fixed" ) {
    warning( "You indicated there are clusters, but these will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }

  # fit uncorrected model
  m0 = pubbias_eta_corrected( yi = yi,
                              vi = vi,
                              sei = sei,
                              eta = 1,
                              model = model,
                              clustervar = clustervar,
                              selection_tails = 1,
                              favor_positive = favor_positive,
                              ci_level = ci_level,
                              small = small )

  # stop if q is on wrong side of null
  if ( m0$est > 0 & q > m0$est ) stop( paste( "The uncorrected pooled point estimate is ", round2(m0$est),
                                              ". q must be less than this value (i.e., closer to zero).",
                                              sep = "" ) )
  if ( m0$est < 0 & q < m0$est ) stop( paste( "The uncorrected pooled point estimate is ", round2(m0$est),
                                              ". q must be greater than this value (i.e., closer to zero).",
                                              sep = "" ) )

  # # reverse signs if needed to have pooled point estimate > 0
  # if ( m0$est < 0 ) {
  #   # keep track so that we can flip back at the end
  #   flipped = TRUE
  #   yi = -yi
  #   q = -q
  # } else {
  #   flipped = FALSE
  # }
  ##### Flip Estimate Signs If Needed #####

  # if favor_positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor_positive == TRUE ) {
    # keep track of whether we flipped for reporting at the end
    flipped = FALSE
  } else {
    flipped = TRUE
    yi = -yi
    q = -q
  }

  # 2-sided p-values for each study even if 1-tailed selection
  pvals = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )

  # affirmative indicator under 1-tailed selection
  A = (pvals < alpha_select) & (yi > 0)

  k.affirmative = sum(A)
  k.nonaffirmative = k.studies - sum(A)

  if ( k.affirmative == 0 | k.nonaffirmative == 0 ) {
    stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
  }

  dat = data.frame( yi, vi, A, clustervar )


  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {

    if (k.nonaffirmative > 1){
      # first fit worst-case meta
      meta.worst = rma.uni( yi = yi,
                            vi = vi,
                            data = dat[ A == FALSE, ],
                            method = "FE" )


      est.worst = as.numeric(meta.worst$b)
      lo.worst = meta.worst$ci.lb
    }

    if (k.nonaffirmative == 1) {
      est.worst = dat$yi[ A == FALSE ]
      lo.worst = dat$yi[ A == FALSE ] - qnorm(0.975) * sqrt(dat$vi[ A == FALSE ])
    }

    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat = dat %>% group_by(A) %>%
      summarise( nu = sum( 1 / vi ),
                 ybar = sum( yi / vi ) )

    # components of bias-corrected estimate by affirmative status
    ybarN = strat$ybar[ strat$A == 0 ]
    ybarA = strat$ybar[ strat$A == 1 ]
    nuN = strat$nu[ strat$A == 0 ]
    nuA = strat$nu[ strat$A == 1 ]

    # S-value for point estimate
    sval.est = ( nuA * q - ybarA ) / ( ybarN - nuN * q )

    # S-value for CI (to shift it to q)
    # match term names used in Wolfram Alpha
    a = ybarN
    b = ybarA
    c = nuN
    d = nuA

    if ( small == FALSE ) k = qnorm( 1 - (alpha/2) )
    if ( small == TRUE ) {
      df = k.studies - 1
      k = qt( 1 - (alpha/2), df = df )
    }

    # # version directly from Wolfram
    # termA = a^2 * d * k^2 - (2 * a * c * d * k^2 * q) +
    #           b^2 * c * k^2 -
    #           (2 * b * c * d * k^2 * q) +
    #           c^2 * d * k^2 * q^2 +
    #           c * d^2 * k^2 * q^2 -
    #           c * d * k^4

    # manually simplied version
    termA = k^2 * ( a^2 * d -
                      (2 * c * d * q) * (a + b) +
                      b^2 * c +
                      q^2 * (c^2 * d + d^2 * c) -
                      c * d * k^2 )

    termB = -a*b + a*d*q + b*c*q - c*d*q^2

    termC = a^2 - 2*a*c*q + c^2*q^2 - c*k^2

    sval.ci = ( -sqrt(termA) + termB ) / termC
    if ( sval.ci < 0 ) sval.ci = ( sqrt(termA) + termB ) / termC

    # # sanity check by inversion
    # # corrected CI limit
    # eta = sval.ci
    # termD = (eta * a + b) / (eta * c + d)
    # termE = k * sqrt( (eta^2 * c + d) / (eta * c + d)^2 )
    # expect_equal( termD - termE,
    #               q )
    # # WORKS!!!

  } # end fixed = TRUE


  ##### Robust Independent and Robust Clustered #####
  if ( model == "robust" ) {

    ##### Worst-Case Meta to See if We Should Search at All

    if (k.nonaffirmative > 1){
      # first fit worst-case meta to see if we should even attempt grid search
      # initialize a dumb (unclustered and uncorrected) version of tau^2
      # which is only used for constructing weights
      meta.re = rma.uni( yi = yi,
                         vi = vi)
      t2hat.naive = meta.re$tau2

      # fit model exactly as in pubbias_eta_corrected
      meta.worst =  robu( yi ~ 1,
                          studynum = clustervar,
                          data = dat[ A == FALSE, ],
                          userweights = 1 / (vi + t2hat.naive),
                          var.eff.size = vi,
                          small = small )

      est.worst = as.numeric(meta.worst$b.r)
      lo.worst = meta.worst$reg_table$CI.L
    }

    # robumeta above can't handle meta-analyzing only 1 nonaffirmative study
    if (k.nonaffirmative == 1) {
      est.worst = dat$yi[ A == FALSE ]
      lo.worst = dat$yi[ A == FALSE ] - qnorm(0.975) * sqrt(dat$vi[ A == FALSE ])
    }

    ##### Get S-value for estimate
    if ( est.worst > q ) {
      sval.est = "Not possible"
    } else {

      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func = function(.eta) {
        est.corr = pubbias_eta_corrected( yi = yi,
                                          vi = vi,
                                          sei = sei,
                                          eta = .eta,
                                          model = model,
                                          clustervar = clustervar,
                                          selection_tails = 1,
                                          favor_positive = TRUE,  # always TRUE because we've already flipped signs if needed
                                          ci_level = ci_level,
                                          small = small )$est
        return( abs(est.corr - q))
      }

      opt = optimize( f = func,
                      interval = c(1, eta_grid_hi),
                      maximum = FALSE )
      sval.est = opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta_grid_hi
      if ( abs(sval.est - eta_grid_hi) < 0.0001 & diff > 0.0001 ) sval.est = paste(">", eta_grid_hi)
    }

    # do something similar for CI
    if ( lo.worst > q ) {
      sval.ci = "Not possible"

    } else {
      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func = function(.eta) {
        lo.corr = pubbias_eta_corrected( yi = yi,
                                         vi = vi,
                                         sei = sei,
                                         eta = .eta,
                                         model = model,
                                         clustervar = clustervar,
                                         selection_tails = 1,
                                         favor_positive = TRUE, # always TRUE because we've already flipped signs if needed
                                         ci_level = ci_level,
                                         small = small )$lo
        return( abs(lo.corr - q))
      }

      opt = optimize( f = func,
                      interval = c(1, eta_grid_hi),
                      maximum = FALSE )
      sval.ci = opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta_grid_hi
      if ( abs(sval.ci - eta_grid_hi) < 0.0001 & diff > 0.0001 ) sval.ci = paste(">", eta_grid_hi)
    }

  }

  # s-values less than 1 indicate complete robustness
  # is.numeric is in case we have a "< XXX" string instead of a number
  if ( is.numeric(sval.est) & !is.na(sval.est) & sval.est < 1) sval.est = "Not possible"
  if ( is.numeric(sval.ci) & !is.na(sval.ci) & sval.ci < 1) sval.ci = "Not possible"

  # m0 was fit BEFORE flipping signs
  # but q has now been flipped in the latter case in "or" statement below
  if ( (m0$est > 0 & m0$lo < q) | (m0$est < 0 & m0$hi > -q) ) {
    # important: Shiny website assumes that this exact string ("--") for CI can be interpreted as
    #  the naive CI's already containing q
    sval.ci = "--"
    message("sval.ci is not applicable because the naive confidence interval already contains q")
  }

  # meta.worst might not exist if, for example, there is only 1 nonaffirmative study
  if ( return_worst_meta == TRUE & exists("meta.worst") ) {
    return( list( stats = data.frame( sval.est,
                                      sval.ci = sval.ci,
                                      k.affirmative,
                                      k.nonaffirmative,
                                      signs.recoded = flipped ),
                  meta.worst = meta.worst ) )
  } else {
    return( data.frame( sval.est,
                        sval.ci = sval.ci,
                        k.affirmative,
                        k.nonaffirmative,
                        signs.recoded = flipped ) )
  }


}


#' @rdname pubbias_svalue
#' @param alpha.select (deprecated) see alpha_select
#' @param eta.grid.hi (deprecated) see eta_grid_hi
#' @param favor.positive (deprecated) see favor_positive
#' @param CI.level (deprecated) see ci_level
#' @param return.worst.meta (deprecated) see return_worst_meta
#' @export
svalue <- function( yi,
                    vi,
                    q,
                    clustervar = 1:length(yi),
                    model,
                    alpha.select = 0.05,
                    eta.grid.hi = 200,
                    favor.positive,
                    CI.level = 0.95,
                    small = TRUE,
                    return.worst.meta = FALSE ) {
  .Deprecated("pubbias_svalue")
  pubbias_svalue(yi = yi,
                 vi = vi,
                 q = q,
                 clustervar = clustervar,
                 model = model,
                 alpha_select = alpha.select,
                 eta_grid_hi = eta.grid.hi,
                 favor_positive = favor.positive,
                 ci_level = CI.level,
                 small = small,
                 return_worst_meta = return.worst.meta)
}
