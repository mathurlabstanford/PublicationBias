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
#' @param return_worst_meta Should the worst-case meta-analysis of only the
#'   nonaffirmative studies be returned?
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
#' @return An object of class [metabias::metabias()].
#'   `stats` includes: the amount of publication bias required to
#'   attenuate the pooled point estimate to `q` (`sval_est`), the amount of
#'   publication bias required to attenuate the confidence interval limit of the
#'   pooled point estimate to `q` (`sval_ci`), the number of affirmative and
#'   nonaffirmative studies after any needed recoding of signs (`k_affirmative`
#'   and `k_nonaffirmative`), and an indicator for whether the point estimates'
#'   signs were recoded (`signs_recoded`).
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

  # stop if selection_ratio doesn't make sense
  if (selection_ratio_max < 1) stop("selection_ratio_max must be at least 1.")

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  # number of point estimates
  k_studies <- length(yi)

  alpha <- 1 - ci_level

  # warn if clusters but user said fixed
  nclusters <- length(unique(cluster))
  if (nclusters < k_studies && model_type == "fixed") {
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

  est0 <- m0$stats$estimate
  q_error <- function(dir) {
    glue("The uncorrected pooled point estimate is {est0}. q must be {dir} than this value (i.e., closer to zero).")
  }
  # stop if q is on wrong side of null
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
  k_nonaffirmative <- k_studies - affirm

  k_zero_msg <- \(dir) glue("There are zero {dir} studies. Model estimation cannot proceed.")
  if (k_affirmative == 0) stop(k_zero_msg("affirmative"))
  if (k_nonaffirmative == 0) stop(k_zero_msg("nonaffirmative"))

  dat <- data.frame(yi, vi, affirm, cluster)
  dat_naff <- dat |> dplyr::filter(!.data$affirm)

  ##### Fixed-Effects Model #####
  if (model_type == "fixed") {

    # special case worst-case meta for 1 non-affirmative study
    if (k_nonaffirmative == 1) {
      est_worst <- dat_naff$yi
      lo_worst <- dat_naff$yi - qnorm(1 - alpha / 2) * sqrt(dat_naff$vi)
    }

    # first fit worst-case meta
    meta_worst <- metafor::rma.uni(yi = dat_naff$yi, vi = dat_naff$vi,
                                   method = "FE", level = ci_level)
    est_worst <- as.numeric(meta_worst$b)
    lo_worst <- meta_worst$ci.lb

    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
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

    if (!small) k <- qnorm(1 - (alpha / 2))
    k <- qt(1 - (alpha / 2), df = k_studies - 1)
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

    ##### Worst-Case Meta to See if We Should Search at All

    if (k_nonaffirmative > 1) {
      # first fit worst-case meta to see if we should even attempt grid search
      # initialize a dumb (unclustered and uncorrected) version of tau^2
      # which is only used for constructing weights
      meta_re <- metafor::rma.uni(yi = yi,
                                  vi = vi)
      t2hat_naive <- meta_re$tau2

      # fit model exactly as in pubbias_meta
      meta_worst <-  robumeta::robu(yi ~ 1,
                                    studynum = cluster,
                                    data = dat[!affirm, ],
                                    userweights = 1 / (vi + t2hat_naive),
                                    var.eff.size = vi,
                                    small = small)

      est_worst <- as.numeric(meta_worst$b.r)
      table_worst <- meta_worst$reg_table
      lo_worst <- est_worst -
        qt(1 - alpha / 2, table_worst$dfs) * table_worst$SE
    }

    # robumeta above can't handle meta-analyzing only 1 nonaffirmative study
    if (k_nonaffirmative == 1) {
      est_worst <- dat$yi[!affirm]
      lo_worst <- dat$yi[!affirm] - qnorm(1 - alpha / 2) * sqrt(dat$vi[!affirm])
    }

    ##### Get S-value for estimate
    if (est_worst > q) {
      sval_est <- "Not possible"
    } else {

      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func <- function(.selection_ratio) {
        corrected <- suppressWarnings(
          pubbias_meta(yi = yi,
                       vi = vi,
                       sei = sei,
                       cluster = cluster,
                       selection_ratio = .selection_ratio,
                       selection_tails = 1,
                       model_type = model_type,
                       favor_positive = TRUE,  # always TRUE because we've already flipped signs if needed
                       ci_level = ci_level,
                       small = small))
        return(abs(corrected$stats$estimate - q))
      }

      opt <- optimize(f = func,
                      interval = c(1, selection_ratio_max),
                      maximum = FALSE)
      sval_est <- opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff <- opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above selection_ratio_max
      if (abs(sval_est - selection_ratio_max) < 0.0001 && diff > 0.0001)
        sval_est <- paste(">", selection_ratio_max)
    }

    # do something similar for CI
    if (lo_worst > q) {
      sval_ci <- "Not possible"

    } else {
      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func <- function(.selection_ratio) {
        corrected <- suppressWarnings(
          pubbias_meta(yi = yi,
                       vi = vi,
                       sei = sei,
                       cluster = cluster,
                       selection_ratio = .selection_ratio,
                       selection_tails = 1,
                       model_type = model_type,
                       favor_positive = TRUE, # always TRUE because we've already flipped signs if needed
                       ci_level = ci_level,
                       small = small))
        return(abs(corrected$stats$ci_lower - q))
      }

      opt <- optimize(f = func,
                      interval = c(1, selection_ratio_max),
                      maximum = FALSE)
      sval_ci <- opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff <- opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above selection_ratio_max
      if (abs(sval_ci - selection_ratio_max) < 0.0001 && diff > 0.0001)
        sval_ci <- paste(">", selection_ratio_max)
    }

  }

  # s-values less than 1 indicate complete robustness
  # is.numeric is in case we have a "< XXX" string instead of a number
  if (is.numeric(sval_est) && !is.na(sval_est) && sval_est < 1)
    sval_est <- "Not possible"
  if (is.numeric(sval_ci) && !is.na(sval_ci) && sval_ci < 1)
    sval_ci <- "Not possible"

  # m0 was fit BEFORE flipping signs
  # but q has now been flipped in the latter case in "or" statement below
  if ((est0 > 0 && m0$stats$ci_lower < q) ||
      (est0 < 0 && m0$stats$ci_upper > -q)) {
    # important: Shiny website assumes that this exact string ("--") for CI can
    # be interpreted as the naive CI's already containing q
    sval_ci <- "--"
    message("sval_ci is not applicable because the naive confidence interval already contains q")
  }

  data <- dat |> dplyr::rename(affirm = .data$A)

  values <- list(
    q = q,
    model_type = model_type,
    alpha_select = alpha_select,
    selection_ratio_max = selection_ratio_max,
    favor_positive = favor_positive,
    ci_level = ci_level,
    small = small,
    k = k_studies,
    k_affirmative = k_affirmative,
    k_nonaffirmative = k_nonaffirmative
  )

  stats <- list(sval_est = sval_est,
                sval_ci = sval_ci)

  fit <- list()
  # meta_worst might not exist, e.g. if there is only 1 nonaffirmative study
  if (return_worst_meta && exists("meta_worst")) fit$meta_worst <- meta_worst

  metabias::metabias(data = data, values = values, stats = stats, fit = fit)

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
