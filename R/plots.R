#' Make significance funnel plot
#'
#' Creates a modified funnel plot that distinguishes between affirmative and nonaffirmative studies, helping to detect the extent to which
#' the nonaffirmative studies' point estimates are systematically smaller than the entire set of point estimates. The estimate among only nonaffirmative studies (gray diamond)
#' represents a corrected estimate under worst-case publication bias. If the gray diamond represents a negligible effect size or if it is much smaller than
#' the pooled estimate among all studies (black diamond), this suggests that the meta-analysis may not be robust to extreme publication bias.
#' Numerical sensitivity analyses (via \code{PublicationBias::pubbias_svalue}) should still be carried out for more precise quantitative conclusions.
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances for the point estimates
#' @param sei A vector of estimated standard errors for the point estimates (only relevant when not using vi)
#' @param xmin x-axis (point estimate) lower limit for plot
#' @param xmax x-axis (point estimate) upper limit for plot
#' @param ymin y-axis (standard error) lower limit for plot
#' @param ymax y-axis (standard error) upper limit for plot
#' @param xlab Label for x-axis (point estimate)
#' @param ylab Label for y-axis (standard error)
#' @param est_all Regular meta-analytic estimate among all studies (optional)
#' @param est_N Worst-case meta-analytic estimate among only nonaffirmative studies (optional)
#' @param favor_positive \code{TRUE} if publication bias is assumed to favor positive estimates; \code{FALSE} if assumed to favor negative estimates.
#' @param alpha_select Alpha-level at which publication probability is assumed to change
#' @param plot_pooled Should the pooled estimates within all studies and within only the nonaffirmative
#' studies be plotted as well?
#' @import
#' metafor
#' stats
#' ggplot2
#' graphics
#' robumeta
#' @details
#' By default (\code{plot_pooled = TRUE}), also plots the pooled point
#' estimate within all studies, supplied by the user as \code{est_all} (black diamond), and within only the nonaffirmative studies, supplied
#' by the user as \code{est_N} (grey diamond). The user can calculate \code{est_all} and \code{est_N} using their choice of meta-analysis model. If instead
#' these are not supplied but \code{plot_pooled = TRUE}, these pooled estimates will be automatically calculated using a fixed-effects (a.k.a. "common-effect") model.
#' @export
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. \emph{Journal of the Royal Statistical Society, Series C.} Preprint available at https://osf.io/s9dp6/.
#' @examples
#'
#' ##### Make Significance Funnel with User-Specified Pooled Estimates #####
#'
#' # compute meta-analytic effect sizes for an example dataset
#' require(metafor)
#' dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#' # flip signs since we think publication bias operates in favor of negative effects
#' # alternatively, if not flipping signs, could pass favor_positive = FALSE to
#' #  significance_funnel
#' dat$yi = -dat$yi
#'
#' # optional: regular meta-analysis of all studies (for the black diamond)
#' # for flexibility, you can use any choice of meta-analysis model here
#' # in this case, we'll use the robust independent specification since the point estimates
#' #  seem to be from unique papers
#' # thus, each study gets its own studynum
#' require(robumeta)
#' meta.all =  robu( yi ~ 1,
#'                   studynum = 1:nrow(dat),
#'                   data = dat,
#'                   var.eff.size = vi,
#'                   small = TRUE )
#'
#' # optional: calculate worst-case estimate (for the gray diamond)
#' #  by analyzing only the nonaffirmative studies
#' dat$pval = 2 * ( 1 - pnorm( abs( dat$yi / sqrt(dat$vi) ) ) )  # two-tailed p-value
#' dat$affirm = (dat$yi > 0) & (dat$pval < 0.05)  # is study affirmative?
#' meta.worst =  robu( yi ~ 1,
#'                     studynum = 1:nrow( dat[ dat$affirm == TRUE, ] ),
#'                     data = dat[ dat$affirm == TRUE, ],
#'                     var.eff.size = vi,
#'                     small = TRUE )
#'
#' ##### Make Significance Funnel with Alpha = 0.50 and Default Pooled Estimates #####
#' # change alpha to 0.50 just for illustration
#' # now the pooled estimates are from the fixed-effect specification because they are
#' #  not provided by the user
#' significance_funnel( yi = dat$yi,
#'                      vi = dat$vi,
#'                      favor_positive = TRUE,
#'                      alpha_select = 0.50,
#'                      plot_pooled = TRUE )

significance_funnel = function( yi,
                                vi,
                                sei,
                                xmin = min(yi),
                                xmax = max(yi),
                                ymin = 0,  # so that pooled points are shown
                                ymax = max( sqrt(vi) ),
                                xlab = "Point estimate",
                                ylab = "Estimated standard error",
                                favor_positive = NA,
                                est_all = NA,
                                est_N = NA,
                                alpha_select = 0.05,
                                plot_pooled = TRUE ) {

  # resolve vi and sei
  if (missing(vi)) {
    if (missing(sei)) {
      stop("Must specify 'vi' or 'sei' argument.")
    }
    vi <- sei ^ 2
  }

  d = data.frame(yi, vi)
  d$sei = sqrt(vi)

  # calculate p-values
  d$pval = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )

  # which direction of effects are favored?
  # if we have the pooled point estimate, but not the favored direction,
  #  assume favored direction matches sign of pooled estimate (but issue warning)
  if ( !is.na(est_all) & is.na(favor_positive) ) {
    favor_positive = (est_all > 0)
    warning("favor_positive not provided, so assuming publication bias favors estimates whose sign matches est_all")
  }
  if ( is.na(est_all) & is.na(favor_positive) ) {
    stop("Need to specify favor_positive")
  }

  # affirmative vs. nonaffirmative indicator
  d$affirm = rep(NA, nrow(d))

  if ( favor_positive == TRUE ) {
    d$affirm[ (d$yi > 0) & (d$pval < alpha_select) ] = "Affirmative"
    d$affirm[ (d$yi < 0) | (d$pval >= alpha_select) ] = "Non-affirmative"
  }
  if ( favor_positive == FALSE ) {
    d$affirm[ (d$yi < 0) & (d$pval < alpha_select) ] = "Affirmative"
    d$affirm[ (d$yi > 0) | (d$pval >= alpha_select) ] = "Non-affirmative"
  }

  # reorder levels for plotting joy
  d$affirm = factor( d$affirm, c("Non-affirmative", "Affirmative") )

  # stop if no studies in either group
  if ( sum( d$affirm == "Non-affirmative" ) == 0 ) {
    stop("There are no non-affirmative studies. The plot would look silly.")
  }

  if ( sum( d$affirm == "Affirmative" ) == 0 ) {
    stop("There are no affirmative studies. The plot would look silly.")
  }

  # pooled fixed-effects estimates
  # if not supplied, gets them from common-effect model
  if ( is.na(est_N) & is.na(est_all) ) {
    est_N = rma.uni(yi = d$yi[ d$affirm == "Non-affirmative" ],
                    vi = d$vi[ d$affirm == "Non-affirmative" ],
                    method="FE")$b

    est_all = rma.uni(yi = d$yi,
                      vi = d$vi,
                      method="FE")$b
  }

  # set up pooled estimates for plotting
  pooled.pts = data.frame( yi = c(est_N, est_all),
                           sei = c(0,0) )

  # for a given SE (y-value), return the "just significant" point estimate value (x-value)
  just_signif_est = function( .sei ) .sei * qnorm(1 - alpha_select/2)

  # calculate slope and intercept of the "just affirmative" line
  # i.e., 1.96 = (just affirmative estimate) / se
  if (favor_positive == TRUE) sl = 1/qnorm(1 - alpha_select/2)
  if (favor_positive == FALSE) sl = -1/qnorm(1 - alpha_select/2)
  int = 0
  # # sanity check: should be exactly alpha_select
  # 2 * ( 1 - pnorm( abs(1) / sl ) )


  ##### Make the Plot #####
  colors = c("darkgray", "orange")

  p.funnel = ggplot( data = d, aes( x = d$yi,
                                    y = d$sei,
                                    color = d$affirm ) )

  if ( plot_pooled == TRUE ) {

    # plot the pooled points
    p.funnel = p.funnel + geom_point(
      data = pooled.pts,
      aes( x = pooled.pts$yi, y = pooled.pts$sei ),
      size = 4,
      shape = 5,
      fill = NA,
      color = c(colors[1], "black")
    ) +

      geom_point(
        data = pooled.pts,
        aes( x = pooled.pts$yi, y = pooled.pts$sei ),
        size = 4,
        shape = 18,
        color = c(colors[1], "black"),
        alpha = 1
      ) +

      # just for visual separation of pooled ests
      geom_hline( yintercept = 0 ) +

      # diagonal "just significant" line
      geom_abline(slope=sl,intercept = int, color = "gray")
  }

  p.funnel = p.funnel +

    # semi-transparent points with solid circles around them
    geom_point( size = 3, alpha=.3) +
    geom_point( size = 3, shape = 1) +

    scale_color_manual(values = colors) +

    xlab(xlab) +
    ylab(ylab) +

    scale_x_continuous( limits = c(xmin, xmax) ) +
    scale_y_continuous( limits = c(ymin, ymax) ) +

    theme_classic() +
    theme(legend.title=element_blank())

  plot(p.funnel)
  return(p.funnel)
}



#' Plot one-tailed p-values
#'
#' Plots the one-tailed p-values. The leftmost red line indicates the cutoff for one-tailed p-values less than 0.025
#' (corresponding to "affirmative" studies; i.e., those with a positive point estimate and a two-tailed p-value
#' less than 0.05). The rightmost red line indicates one-tailed p-values greater than 0.975 (i.e., studies with a
#' negative point estimate and a two-tailed p-value less than 0.05). If there is a substantial point mass of p-values
#' to the right of the rightmost red line, this suggests that selection may be two-tailed rather than one-tailed.
#' @param yi A vector of point estimates to be meta-analyzed. The signs of the estimates should be chosen
#' such that publication bias is assumed to operate in favor of positive estimates.
#' @param vi A vector of estimated variances for the point estimates
#' @param sei A vector of estimated standard errors for the point estimates (only relevant when not using vi)
#' @param alpha_select Alpha-level at which publication probability is assumed to change
#' @import
#' stats
#' ggplot2
#' @export
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. \emph{Journal of the Royal Statistical Society, Series C.} Preprint available at https://osf.io/s9dp6/.
#' @examples
#'
#'  # compute meta-analytic effect sizes
#'  require(metafor)
#'  dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#'  # flip signs since we think publication bias operates in favor of negative effects
#'  dat$yi = -dat$yi
#'
#'  pval_plot( yi = dat$yi,
#'             vi = dat$vi )

pval_plot = function( yi,
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
  pval = 1 - pnorm( yi / sqrt(vi) )

  ggplot( data = data.frame(pval = pval),
          aes( x = pval ) ) +
    geom_vline(xintercept = alpha_select/2, color = "red", lwd = 1) +
    geom_vline(xintercept = 1 - (alpha_select/2), color = "red", lwd = 1) +
    geom_histogram( binwidth = 0.025 ) +
    xlab("One-tailed p-value") +
    theme_classic() +
    theme( panel.grid = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text=element_text(size=16),
           axis.title=element_text(size=16, face = "bold") )
}
