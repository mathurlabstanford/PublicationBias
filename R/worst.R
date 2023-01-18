#' @keywords internal
fit_meta_worst <- function(dat, model_type, ci_level, small) {

  dat_naff <- dat |> dplyr::filter(!.data$affirm)
  alpha <- 1 - ci_level
  k_nonaffirmative <- nrow(dat_naff)

  # special case worst-case meta for 1 non-affirmative study
  if (k_nonaffirmative == 1) {
    meta_worst <- NULL

    stats_worst <- tibble(estimate = dat_naff$yi,
                          se = sqrt(dat_naff$vi),
                          ci_lower = .data$estimate - qnorm(1 - alpha / 2) * .data$se,
                          ci_upper = .data$estimate + qnorm(1 - alpha / 2) * .data$se,
                          p_value = 2 * (1 - pnorm(abs(.data$estimate / .data$se))))

  } else {

    if (model_type == "fixed") {
      meta_worst <- metafor::rma.uni(yi = dat_naff$yi, vi = dat_naff$vi,
                                     method = "FE", level = ci_level)

      stats_worst <- tibble(estimate = as.numeric(meta_worst$b),
                            se = meta_worst$se,
                            ci_lower = meta_worst$ci.lb,
                            ci_upper = meta_worst$ci.ub,
                            p_value = meta_worst$pval)
    }

    if (model_type == "robust") {

      # initialize a naive (unclustered and uncorrected) version of tau^2
      # which is only used for constructing weights
      t2hat_naive <- metafor::rma.uni(yi = dat$yi, vi = dat$vi)$tau2

      # fit model exactly as in pubbias_meta()
      meta_worst <- robumeta::robu(yi ~ 1,
                                   studynum = dat_naff$cluster,
                                   data = dat_naff,
                                   userweights = 1 / (dat_naff$vi + t2hat_naive),
                                   var.eff.size = dat_naff$vi,
                                   small = small)

      stats_worst <- metabias::robu_ci(meta_worst) |> select(-.data$param)
    }
  }

  stats_worst <- stats_worst |>
    mutate(model ="worst_case", .before = everything())

  list(meta = meta_worst, stats = stats_worst)
}
