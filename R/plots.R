#' Make significance funnel plot
#'
#' Creates a modified funnel plot that distinguishes between affirmative and
#' nonaffirmative studies, helping to detect the extent to which the
#' nonaffirmative studies' point estimates are systematically smaller than the
#' entire set of point estimates. The estimate among only nonaffirmative studies
#' (gray diamond) represents a corrected estimate under worst-case publication
#' bias. If the gray diamond represents a negligible effect size or if it is
#' much smaller than the pooled estimate among all studies (black diamond), this
#' suggests that the meta-analysis may not be robust to extreme publication
#' bias. Numerical sensitivity analyses (via [pubbias_svalue()]) should still be
#' carried out for more precise quantitative conclusions.
#'
#' @export
#'
#' @inheritParams metabias::params
#' @param xmin x-axis (point estimate) lower limit for plot.
#' @param xmax x-axis (point estimate) upper limit for plot.
#' @param ymin y-axis (standard error) lower limit for plot.
#' @param ymax y-axis (standard error) upper limit for plot.
#' @param xlab Label for x-axis (point estimate).
#' @param ylab Label for y-axis (standard error).
#' @param est_all Regular meta-analytic estimate among all studies (optional).
#' @param est_worst Worst-case meta-analytic estimate among only nonaffirmative
#'   studies (optional).
#' @param plot_pooled Should the pooled estimates within all studies and within
#'   only the nonaffirmative studies be plotted as well?
#'
#' @details By default (`plot_pooled = TRUE`), also plots the pooled point
#'   estimate within all studies, supplied by the user as `est_all` (black
#'   diamond), and within only the nonaffirmative studies, supplied by the user
#'   as `est_worst` (gray diamond). The user can calculate `est_all` and
#'   `est_worst` using their choice of meta-analysis model. If instead these
#'   are not supplied but `plot_pooled = TRUE`, these pooled estimates will
#'   be automatically calculated using a fixed-effects (a.k.a. "common-effect")
#'   model.
#'
#' @references
#' \insertRef{mathur2020}{metabias}
#'
#' @example inst/examples/significance_funnel.R
significance_funnel <- function(yi,
                                vi,
                                sei,
                                favor_positive = TRUE,
                                alpha_select = 0.05,
                                plot_pooled = TRUE,
                                est_all = NA,
                                est_worst = NA,
                                xmin = min(yi),
                                xmax = max(yi),
                                ymin = 0,  # so that pooled points are shown
                                ymax = max(sqrt(vi)),
                                xlab = "Point estimate",
                                ylab = "Estimated standard error") {

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  d <- tibble(yi, vi, sei = sqrt(vi),
              pval = 2 * (1 - pnorm(abs(yi) / sqrt(vi))))

  if (!is.na(est_all) && (est_all > 0) != favor_positive)
    warning("Favored direction is opposite of the pooled estimate.")

  if (favor_positive) d <- d |>
    mutate(affirm = .data$yi > 0 & .data$pval < alpha_select)
  if (!favor_positive) d <- d |>
    mutate(affirm = .data$yi < 0 & .data$pval < alpha_select)

  # stop if no studies in either group
  if (all(d$affirm)) stop("There are no non-affirmative studies.")
  if (!any(d$affirm)) stop("There are no affirmative studies.")

  # pooled fixed-effects estimates
  # if not supplied, gets them from common-effect model
  if (is.na(est_worst) && is.na(est_all)) {
    d_naff <- d |> dplyr::filter(!.data$affirm)
    est_worst <- metafor::rma.uni(yi = d_naff$yi, vi = d_naff$vi,
                                  method = "FE")$b
    est_all <- metafor::rma.uni(yi = d$yi, vi = d$vi, method = "FE")$b
  }
  d <- d |>
    mutate(affirm = factor(.data$affirm,
                           labels = c("Non-affirmative", "Affirmative")))

  # set up pooled estimates for plotting
  pooled_pts <- data.frame(yi = c(est_worst, est_all), sei = c(0, 0))

  # calculate slope and intercept of the "just affirmative" line
  # i.e., 1.96 = (just affirmative estimate) / se
  sl <- 1 / qnorm(1 - alpha_select / 2)
  if (!favor_positive) sl <- -sl

  ##### Make the Plot #####
  colors <- c("darkgray", "orange")

  p_funnel <- ggplot(d, aes(x = .data$yi, y = .data$sei, color = .data$affirm))

  if (plot_pooled) {

    # plot the pooled points
    p_funnel <- p_funnel +
      geom_point(aes(x = .data$yi, y = .data$sei), data = pooled_pts,
                 size = 4, shape = 5, fill = NA,
                 color = c(colors[1], "black")) +

      geom_point(aes(x = .data$yi, y = .data$sei), data = pooled_pts,
                 size = 4, shape = 18, alpha = 1,
                 color = c(colors[1], "black")) +

      # just for visual separation of pooled ests
      geom_hline(yintercept = 0) +

      # diagonal "just significant" line
      geom_abline(slope = sl, intercept = 0, color = "gray")
  }

  p_funnel <- p_funnel +

    # semi-transparent points with solid circles around them
    geom_point(size = 3, alpha = .3) +
    geom_point(size = 3, shape = 1) +

    scale_color_manual(values = colors) +
    scale_x_continuous(name = xlab, limits = c(xmin, xmax)) +
    scale_y_continuous(name = ylab, limits = c(ymin, ymax)) +

    theme_classic() +
    theme(legend.title = element_blank())

  return(p_funnel)
}



#' Plot one-tailed p-values
#'
#' Plots the one-tailed p-values. The leftmost red line indicates the cutoff for
#' one-tailed p-values less than 0.025 (corresponding to "affirmative" studies;
#' i.e., those with a positive point estimate and a two-tailed p-value less than
#' 0.05). The rightmost red line indicates one-tailed p-values greater than
#' 0.975 (i.e., studies with a negative point estimate and a two-tailed p-value
#' less than 0.05). If there is a substantial point mass of p-values to the
#' right of the rightmost red line, this suggests that selection may be
#' two-tailed rather than one-tailed.
#' @export
#'
#' @inheritParams metabias::params
#' @param yi A vector of point estimates to be meta-analyzed. The signs of the
#'   estimates should be chosen such that publication bias is assumed to operate
#'   in favor of positive estimates.
#'
#' @references
#' \insertRef{mathur2020}{metabias}
#'
#' @examples
#' # compute meta-analytic effect sizes
#' require(metafor)
#' dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
#'                        di = cneg, data = dat.bcg)
#'
#' # flip signs since we think publication bias favors negative effects
#' dat$yi <- -dat$yi
#'
#' pval_plot(yi = dat$yi, vi = dat$vi)
pval_plot <- function(yi,
                      vi,
                      sei,
                      alpha_select = 0.05) {

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  # calculate 1-tailed p-values
  pval <- 1 - pnorm(yi / sqrt(vi))

  ggplot(data.frame(pval = pval), aes(x = .data$pval)) +
    geom_vline(xintercept = alpha_select / 2, color = "red", lwd = 1) +
    geom_vline(xintercept = 1 - (alpha_select / 2), color = "red", lwd = 1) +
    geom_histogram(binwidth = 0.025) +
    xlab("One-tailed p-value") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16, face = "bold"))
}
