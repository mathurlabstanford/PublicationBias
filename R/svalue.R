#' Severity of publication bias needed to "explain away" results
#'
#' Estimates the S-value, defined as the severity of publication bias (i.e., the
#' ratio by which affirmative studies are more likely to be published than
#' nonaffirmative studies) that would be required to shift the pooled point
#' estimate or its confidence interval limit to the value `q`.
#'
#' @inheritParams metabias::params
#' @inheritParams pubbias_meta
#' @param selection_ratio_max The largest value of `selection_ratio` that should
#'   be included in the grid search. This argument is only needed when
#'   `model_type = "robust"`.
#'
#' @details To illustrate interpretation of the S-value, if the S-value for the
#'   point estimate is 30 with `q=0`, this indicates that affirmative studies
#'   (i.e., those with a "statistically significant" and positive estimate)
#'   would need to be 30-fold more likely to be published than nonaffirmative
#'   studies (i.e., those with a "nonsignificant" or negative estimate) to
#'   attenuate the pooled point estimate to `q`.
#'
#'   If `favor_positive = FALSE`, such that publication bias is assumed to favor
#'   negative rather than positive estimates, the signs of `yi` will be reversed
#'   prior to performing analyses. The returned number of affirmative and
#'   nonaffirmative studies will reflect the recoded signs.
#'
#'   If `return_worst_meta = TRUE`, also returns the worst-case meta-analysis of
#'   only the nonaffirmative studies. If `model_type = "fixed"`, the worst-case
#'   meta-analysis is fit by `metafor::rma.uni()`. If `model_type = "robust"`,
#'   it is fit by `robumeta::robu()`. Note that in the latter case, custom
#'   inverse-variance weights are used, which are the inverse of the sum of the
#'   study's variance and a heterogeneity estimate from a naive random-effects
#'   meta-analysis (Mathur & VanderWeele, 2020). This is done for consistency
#'   with the results of `pubbias_meta()`, which is used to determine `sval_est`
#'   and `sval_ci`. Therefore, the worst-case meta-analysis results may differ
#'   slightly from what you would obtain if you simply fit `robumeta::robu()` on
#'   the nonaffirmative studies with the default weights.
#'
#' @return An object of class [metabias::metabias()], a list containing:
#' \describe{
#'   \item{data}{A tibble with one row per study and the columns
#'               `r meta_names_str("data")`.}
#'   \item{values}{A list with the elements `r meta_names_str("values")`.}
#'   \item{stats}{A tibble with the columns `r meta_names_str("stats")`.
#'                `sval_est` represents the amount of publication bias required
#'                to attenuate the pooled point estimate to `q`; `sval_ci`
#'                represents the amount of publication bias required to
#'                attenuate the confidence interval limit of the pooled point
#'                estimate to `q`.}
#'   \item{fit}{A list of fitted models, if any.}
#' }
#'
#' @references
#' \insertRef{mathur2020}{metabias}
#'
#' @export
#' @example inst/examples/pubbias_svalue.R
pubbias_svalue <- function(yi, # data
                           vi,
                           sei,
                           cluster = 1:length(yi),

                           q = 0, # params

                           model_type = "robust", # opts
                           favor_positive = TRUE,
                           alpha_select = 0.05,
                           ci_level = 0.95,
                           small = TRUE,
                           selection_ratio_max = 200,
                           return_worst_meta = FALSE) {

  # error if selection_ratio_max doesn't make sense
  if (selection_ratio_max < 1) stop("selection_ratio_max must be at least 1.")

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

  # warn if clusters are present but model_type == "fixed"
  nclusters <- length(unique(cluster))
  if (nclusters < k && model_type == "fixed") {
    warning("You indicated there are clusters, but these will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model_type = robust.")
  }

  # fit uncorrected model
  m0 <- pubbias_meta(yi = yi,
                     vi = vi,
                     sei = sei,
                     cluster = cluster,
                     selection_ratio = 1,
                     selection_tails = 1,
                     model_type = model_type,
                     favor_positive = favor_positive,
                     ci_level = ci_level,
                     small = small)

  # error if q is on wrong side of null
  est0 <- m0$stats$estimate
  q_error <- function(dir) {
    glue("The uncorrected pooled point estimate is {est0}. q must be {dir} than this value (i.e., closer to zero).")
  }
  if (est0 > 0 && q > est0) stop(q_error("less"))
  if (est0 < 0 && q < est0) stop(q_error("greater"))

  ##### Flip Estimate Signs If Needed #####
  if (!favor_positive) {
    yi <- -yi
    q <- -q
  }

  # 2-sided p-values for each study even if 1-tailed selection
  pvals <- 2 * (1 - pnorm(abs(yi) / sqrt(vi)))

  # affirmative indicator under 1-tailed selection
  affirm <- (pvals < alpha_select) & (yi > 0)

  k_affirmative <- sum(affirm)
  k_nonaffirmative <- k - k_affirmative

  # error if there are zero affirmative or zero nonaffirmative studies
  k_zero_msg <- \(dir) glue("There are zero {dir} studies. Model estimation cannot proceed.")
  if (k_affirmative == 0) stop(k_zero_msg("affirmative"))
  if (k_nonaffirmative == 0) stop(k_zero_msg("nonaffirmative"))

  dat <- data.frame(yi, vi, affirm, cluster)
  fits <- list()

  # fit worst-case meta-analysis of only nonaffirmative studies
  # always need model for svalue search initilization
  meta_worst <- fit_meta_worst(dat, model_type = model_type,
                               ci_level = ci_level, small = small)
  if (return_worst_meta) fits$meta_worst <- meta_worst$meta

  ##### Fixed-Effects Model #####
  if (model_type == "fixed") {

    strat <- dat |>
      group_by(.data$affirm) |>
      summarise(nu = sum(1 / .data$vi), ybar = sum(.data$yi / .data$vi))

    # components of bias-corrected estimate by affirmative status
    ybar_n <- strat$ybar[!strat$affirm]
    ybar_a <- strat$ybar[strat$affirm]
    nu_n <- strat$nu[!strat$affirm]
    nu_a <- strat$nu[strat$affirm]

    # S-value for point estimate
    sval_est <- (nu_a * q - ybar_a) / (ybar_n - nu_n * q)

    # S-value for CI (to shift it to q)
    # match term names used in Wolfram Alpha
    a <- ybar_n
    b <- ybar_a
    c <- nu_n
    d <- nu_a

    k <- if (!small) qnorm(1 - (alpha / 2)) else qt(1 - (alpha / 2),
                                                    df = k - 1)
    term_a <- k ^ 2 * (a ^ 2 * d -
                         (2 * c * d * q) * (a + b) +
                         b ^ 2 * c +
                         q ^ 2 * (c ^ 2 * d + d ^ 2 * c) -
                         c * d * k ^ 2)

    term_b <- -a * b + a * d * q + b * c * q - c * d * q ^ 2

    term_c <- a ^ 2 - 2 * a * c * q + c ^ 2 * q ^ 2 - c * k ^ 2

    sval_ci <- (-sqrt(term_a) + term_b) / term_c
    if (sval_ci < 0) sval_ci <- (sqrt(term_a) + term_b) / term_c

  } # end fixed = TRUE


  ##### Robust Independent and Robust Clustered #####
  if (model_type == "robust") {

    .meta_fun <- function(.selection_ratio) {
      pubbias_meta(yi = yi,
                   vi = vi,
                   sei = sei,
                   cluster = cluster,
                   selection_ratio = .selection_ratio,
                   selection_tails = 1,
                   model_type = model_type,
                   favor_positive = TRUE,  # always TRUE because we've already flipped signs if needed
                   ci_level = ci_level,
                   small = small)
    }

    .svalue_fun <- function(.param) {
      find_svalue(param = .param, meta_fun = .meta_fun,
                  meta_worst = meta_worst, q = q,
                  selection_ratio_max = selection_ratio_max)
    }

    sval_est <- .svalue_fun("estimate")
    sval_ci <- .svalue_fun("ci_lower")
  }

  # check if naive confidence interval contains q
  # m0 was fit BEFORE flipping signs
  # but q has now been flipped in the latter case in "or" statement below
  if ((est0 > 0 && m0$stats$ci_lower < q) ||
      (est0 < 0 && m0$stats$ci_upper > -q)) {
    # important: Shiny website assumes that this exact string ("--") for CI can
    # be interpreted as the naive CI's already containing q
    sval_ci <- "--"
    message("sval_ci is not applicable because the naive confidence interval already contains q")
  }

  values <- list(
    q = q,
    model_type = model_type,
    favor_positive = favor_positive,
    alpha_select = alpha_select,
    ci_level = ci_level,
    small = small,
    selection_ratio_max = selection_ratio_max,
    k = k,
    k_affirmative = k_affirmative,
    k_nonaffirmative = k_nonaffirmative
  )

  stats <- tibble(sval_est = sval_est, sval_ci = sval_ci)

  metabias::metabias(data = dat, values = values, stats = stats, fits = fits)

}

#' @keywords internal
find_svalue <- function(param, meta_fun, meta_worst, q, selection_ratio_max,
                        tol = 0.0001) { # param is estimate or ci_lower

  if (meta_worst$stats[[param]] > q) return("Not possible")

  # define the function we need to minimize
  # i.e., distance between corrected estimate and the target value of q
  func <- function(.selection_ratio) {
    corrected <- suppressWarnings(meta_fun(.selection_ratio))
    return(abs(corrected$stats[[param]] - q))
  }

  opt <- optimize(f = func,
                  interval = c(1, selection_ratio_max),
                  maximum = FALSE)
  sval <- opt$minimum

  # discrepancy between the corrected estimate and the s-value
  diff <- opt$objective

  # s-values less than 1 indicate complete robustness
  if (!is.na(sval) && sval < 1) return("Not possible")

  # if the optimal value is very close to the upper range of grid search
  #  AND we're still not very close to the target q,
  #  that means the optimal value was above selection_ratio_max
  if (abs(sval - selection_ratio_max) < tol && diff > tol)
    sval <- paste(">", selection_ratio_max)

  return(sval)
}


#' @rdname pubbias_svalue
#' @param clustervar (deprecated) see cluster
#' @param model (deprecated) see model_type
#' @param alpha.select (deprecated) see alpha_select
#' @param eta.grid.hi (deprecated) see selection_ratio_max
#' @param favor.positive (deprecated) see favor_positive
#' @param CI.level (deprecated) see ci_level
#' @param return.worst.meta (deprecated) see return_worst_meta
#' @export
svalue <- function(yi,
                   vi,
                   q,
                   clustervar = 1:length(yi),
                   model,
                   alpha.select = 0.05,
                   eta.grid.hi = 200,
                   favor.positive,
                   CI.level = 0.95,
                   small = TRUE,
                   return.worst.meta = FALSE) {
  .Deprecated("pubbias_svalue")
  pubbias_svalue(yi = yi,
                 vi = vi,
                 cluster = clustervar,
                 q = q,
                 model_type = model,
                 favor_positive = favor.positive,
                 alpha_select = alpha.select,
                 ci_level = CI.level,
                 small = small,
                 selection_ratio_max = eta.grid.hi,
                 return_worst_meta = return.worst.meta)
}
