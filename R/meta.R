#' Estimate publication bias-corrected meta-analysis
#'
#' For a chosen ratio of publication probabilities, `selection_ratio`, estimates
#' a publication bias-corrected pooled point estimate and confidence interval
#' per \insertCite{mathur2020;textual}{metabias}. Model options include
#' fixed-effects (a.k.a. "common-effect"), robust independent, and robust
#' clustered specifications.
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
#' @return An object of class [metabias::metabias()], a list containing:
#' \describe{
#'   \item{data}{A tibble with one row per study and the columns
#'               `r meta_names_str("data")`.}
#'   \item{values}{A list with the elements `r meta_names_str("values")`.}
#'   \item{stats}{A tibble with the columns `r meta_names_str("stats")`.}
#'   \item{fit}{A list of fitted models, if any.}
#' }
#'
#' @references
#' \insertRef{mathur2020}{metabias}
#'
#' @export
#' @example inst/examples/pubbias_meta.R
pubbias_meta <- function(yi, # data
                         vi,
                         sei,
                         cluster = 1:length(yi),

                         selection_ratio, # params

                         selection_tails = 1, # opts
                         model_type = "robust",
                         favor_positive = TRUE,
                         alpha_select = 0.05,
                         ci_level = 0.95,
                         small = TRUE) {

  # stop if selection_ratio doesn't make sense
  if (selection_ratio < 1) stop("selection_ratio must be at least 1.")

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  # number of point estimates
  k <- length(yi)

  # calculate alpha for inference on point estimate
  alpha <- 1 - ci_level

  # warn if clusters but user said fixed
  nclusters <- length(unique(cluster))
  if (nclusters < k && model_type == "fixed") {
    warning('Clusters exist, but will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model_type = "robust".')
  }

  # warn if naive estimate is in opposite direction than favor_positive
  naive_pos <- metafor::rma(yi, vi, method = "FE")$beta > 0
  if (naive_pos != favor_positive)
    warning("Favored direction is opposite of the pooled estimate.")

  ##### Flip Estimate Signs If Needed #####
  yif <- if (favor_positive)  yi else -yi

  # 2-sided p-values for each study even if 1-tailed selection
  pvals <- 2 * (1 - pnorm(abs(yif) / sqrt(vi)))

  # affirmative indicator based on selection tails
  if (selection_tails == 1) affirm <- (pvals < alpha_select) & (yif > 0)
  if (selection_tails == 2) affirm <- (pvals < alpha_select)

  k_affirmative <- sum(affirm)
  k_nonaffirmative <- k - k_affirmative

  k_zero_msg <- \(dir) glue("There are zero {dir} studies. Model estimation cannot proceed.")
  if (k_affirmative == 0) stop(k_zero_msg("affirmative"))
  if (k_nonaffirmative == 0) stop(k_zero_msg("nonaffirmative"))

  dat <- tibble(yi, yif, vi, affirm, cluster)


  ##### Fixed-Effects Model #####
  if (model_type == "fixed") {

    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat <- dat |>
      group_by(.data$affirm) |>
      summarise(nu = sum(1 / .data$vi), ybar = sum(.data$yi / .data$vi))

    # components of bias-corrected estimate by affirmative status
    ybar_n <- strat$ybar[!strat$affirm]
    ybar_s <- strat$ybar[strat$affirm]
    nu_n <- strat$nu[!strat$affirm]
    nu_s <- strat$nu[strat$affirm]

    # corrected pooled point estimate
    est <- (selection_ratio * ybar_n + ybar_s) / (selection_ratio * nu_n + nu_s)

    # inference
    var <- (selection_ratio ^ 2 * nu_n + nu_s) / (selection_ratio * nu_n + nu_s) ^ 2
    se <- sqrt(var)

    if (!small) {
      # z-based inference
      lo <- est - qnorm(1 - alpha / 2) * sqrt(var)
      hi <- est + qnorm(1 - alpha / 2) * sqrt(var)
      z <- abs(est / sqrt(var))
      pval_est <- 2 * (1 - pnorm(z))
    } else {
      # t-based inference
      df <- k - 1
      lo <- est - qt(1 - alpha / 2, df = df) * sqrt(var)
      hi <- est + qt(1 - alpha / 2, df = df) * sqrt(var)
      t <- abs(est / sqrt(var))
      pval_est <- 2 * (1 - pt(t, df = df))
    }

    stats <- tibble(estimate = est,
                    se = se,
                    ci_lower = lo,
                    ci_upper = hi,
                    p_value = pval_est)
    fits <- list()

  } # end fixed = TRUE

  ##### Robust Independent and Robust Clustered #####
  if (model_type == "robust") {

    # weight for model
    weights <- rep(1, length(pvals))
    weights[!affirm] <- selection_ratio

    # initialize a naive (unclustered and uncorrected) version of tau^2
    # which is only used for constructing weights
    meta_re <- metafor::rma.uni(yi = yi, vi = vi)
    t2hat_naive <- meta_re$tau2

    # fit weighted robust model
    meta_robu <- robumeta::robu(yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = weights / (vi + t2hat_naive),
                                var.eff.size = vi,
                                small = small)

    stats <- metabias::robu_ci(meta_robu)
    fits <- list("robust" = meta_robu)
  } # end robust = TRUE

  values <- list(selection_ratio = selection_ratio,
                 selection_tails = selection_tails,
                 model_type = model_type,
                 favor_positive = favor_positive,
                 alpha_select = alpha_select,
                 ci_level = ci_level,
                 small = small,
                 k = k,
                 k_affirmative = k_affirmative,
                 k_nonaffirmative = k_nonaffirmative)

  metabias::metabias(data = dat, values = values, stats = stats, fits = fits)

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
corrected_meta <- function(yi,
                           vi,
                           eta,
                           clustervar = 1:length(yi),
                           model,
                           selection.tails = 1,
                           favor.positive,
                           alpha.select = 0.05,
                           CI.level = 0.95,
                           small = TRUE) {
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
