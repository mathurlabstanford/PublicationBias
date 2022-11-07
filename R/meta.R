#' Estimate publication bias-corrected meta-analysis
#'
#' For a chosen ratio of publication probabilities, `eta`, estimates a
#' publication bias-corrected pooled point estimate and confidence interval per
#' Mathur & VanderWeele (2020). Model options include fixed-effects (a.k.a.
#' "common-effect"), robust independent, and robust clustered specifications.
#'
#' @inheritParams metabias::params
#' @param selection_tails 1 (for one-tailed selection, recommended for its
#'   conservatism) or 2 (for two-tailed selection).
#' @param model_type "fixed" for fixed-effects (a.k.a. "common-effect") or
#'   "robust" for robust random-effects.
#'
#' @details The `selection_ratio` represents the number of times more likely
#'   affirmative studies (i.e., those with a "statistically significant" and
#'   positive estimate) are to be published than nonaffirmative studies (i.e.,
#'   those with a "nonsignificant" or negative estimate).
#'
#'   If `favor_positive` is `FALSE`, such that publication bias is assumed to
#'   favor negative rather than positive estimates, the signs of `yi` will be
#'   reversed prior to performing analyses. The corrected estimate will be
#'   reported based on the recoded signs rather than the original sign
#'   convention.
#'
#' @return A list with three elements, `values`, `stats` and `fit`. Stats is a
#'   list that contains the bias-corrected pooled point estimate (`estimate`)
#'   and inference on the bias-corrected estimate (`se`, `ci_lower`, `ci_upper`,
#'   `p_value`). Values is a list that contains the user's specified
#'   `selection_ratio`, the number of affirmative and nonaffirmative studies
#'   (`k_affirmative` and `k_nonaffirmative`), and a dataframe combining `yi`,
#'   `vi`, `cluster`.
#'
#' @references Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for
#'   publication bias in meta-analyses. *Journal of the Royal Statistical
#'   Society, Series C.* Preprint available at https://osf.io/s9dp6/.
#'
#' @export
#' @example inst/examples/pubbias_meta.R
pubbias_meta = function(yi, # data
                        vi,
                        sei,
                        cluster = 1:length(yi),

                        selection_ratio, # params

                        model_type = "robust", # opts
                        favor_positive = TRUE,
                        selection_tails = 1,
                        alpha_select = 0.05,
                        ci_level = 0.95,
                        small = TRUE) {

  # stop if selection_ratio doesn't make sense
  if ( selection_ratio < 1 ) stop( "selection_ratio must be at least 1.")

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
  nclusters = length( unique( cluster ) )
  if ( nclusters < k & model_type == "fixed" ) {
    warning( "Clusters exist, but will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model_type = robust.")
  }

  # warn if naive estimate is in opposite direction than favor_positive
  naive_pos <- metafor::rma(yi, vi, method = "FE")$beta > 0
  if (naive_pos != favor_positive)
    warning("Favored direction is opposite of the pooled estimate.")

  ##### Flip Estimate Signs If Needed #####
  if ( favor_positive ) yif = yi else yif = -yi

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

  dat = data.frame( yi, yif, vi, A, cluster )


  ##### Fixed-Effects Model #####
  if ( model_type == "fixed" ) {

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
    est = ( selection_ratio * ybarN + ybarS ) / ( selection_ratio * nuN + nuS )

    # inference
    var = ( selection_ratio^2 * nuN + nuS ) / ( selection_ratio * nuN + nuS )^2
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
  if ( model_type == "robust" ) {

    # weight for model
    weights = rep( 1, length(pvals) )
    weights[ A == FALSE ] = selection_ratio

    # initialize a dumb (unclustered and uncorrected) version of tau^2
    # which is only used for constructing weights
    meta_re = metafor::rma.uni( yi = yi,
                                vi = vi)
    t2hat_naive = meta_re$tau2

    # fit weighted robust model
    meta_robu = robumeta::robu( yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = weights / (vi + t2hat_naive),
                                var.eff.size = vi,
                                small = small )

    est = as.numeric(meta_robu$b.r)
    se = meta_robu$reg_table$SE
    reg_table <- meta_robu$reg_table
    lo <- est - qt(1 - alpha / 2, reg_table$dfs) * reg_table$SE
    hi <- est + qt(1 - alpha / 2, reg_table$dfs) * reg_table$SE
    pval_est = meta_robu$reg_table$prob
  } # end robust = TRUE

  data = dat %>% dplyr::rename(affirm = .data$A)
  values = list(selection_ratio = selection_ratio,
                model_type = model_type,
                selection_tails = selection_tails,
                favor_positive = favor_positive,
                alpha_select = alpha_select,
                ci_level = ci_level,
                small = small,
                k = k,
                k_affirmative = k_affirmative,
                k_nonaffirmative = k_nonaffirmative)

  stats = list(estimate = est,
               se = se,
               ci_lower = lo,
               ci_upper = hi,
               p_value = pval_est)

  fit = list()
  if ( exists("meta_robu") ) {
    fit$robust = meta_robu
  }

  results <- list(data = data, values = values, stats = stats, fit = fit)
  class(results) <- "metabias"
  return(results)

}


#' @rdname pubbias_meta
#' @param eta (deprecated) see selection_ratio
#' @param clustervar (deprecated) see cluster
#' @param model (deprecated) see model_type
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
  .Deprecated("pubbias_meta")
  pubbias_meta(yi = yi,
                        vi = vi,
                        cluster = clustervar,
                        selection_ratio = eta,
                        selection_tails = selection.tails,
                        model_type = model,
                        favor_positive = favor.positive,
                        alpha_select = alpha.select,
                        ci_level = CI.level,
                        small = small)
}
