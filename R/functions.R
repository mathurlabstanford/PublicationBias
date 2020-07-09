

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             FNs FOR MODEL 1 (FIXED EFFECTS)                         #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


################################ FN: CORRECTED ESTIMATE ################################

#' Estimate publication bias-corrected meta-analysis
#'
#' For a chosen ratio of publication probabilities, \code{eta}, estimates a publication bias-corrected pooled point
#' estimate and confidence interval per Mathur & VanderWeele (2020). Model options include fixed-effects (a.k.a. "common-effects"), robust independent, and robust
#' clustered specifications.
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances for the point estimates
#' @param eta The number of times more likely an affirmative study is to be published than a nonaffirmative study; see Details
#' @param clustervar A character, factor, or numeric vector with the same length as yi. Unique values should indicate
#' unique clusters of point estimates. By default, assumes all point estimates are independent.
#' @param model "fixed" for fixed-effects (a.k.a. "common-effects") or "robust" for robust random-effects
#' @param selection.tails 1 (for one-tailed selection, recommended for its conservatism) or 2 (for two-tailed selection)
#' @param favor.positive \code{TRUE} if publication bias is assumed to favor positive estimates; \code{FALSE} if assumed to favor negative estimates.
#' See Details.
#' @param alpha.select Alpha-level at which publication probability is assumed to change
#' @param CI.level Confidence interval level (as proportion) for the corrected point estimate. (The alpha level for inference on the corrected
#' point estimate will be calculated from \code{CI.level}.)
#' @param small Should inference allow for a small meta-analysis? We recommend always using TRUE.
#' @import
#' metafor
#' stats
#' robumeta
#' ggplot2
#' graphics
#' MetaUtility
#' @importFrom
#' dplyr %>% group_by summarise
#' @export
#' @details
#' The ratio \code{eta} represents the number of times more likely affirmative studies (i.e., those with a "statistically significant" and positive estimate)
#' are to be published than nonaffirmative studies (i.e., those with a "nonsignificant" or negative estimate).
#'
#' If \code{favor.positive == FALSE}, such that publication bias is assumed to favor negative rather than positive estimates, the signs of \code{yi} will be reversed prior to
#' performing analyses. The corrected estimate will be reported based on the recoded signs rather than the original sign convention, and accordingly the returned value \code{signs.recoded} will be \code{TRUE}.
#' @return
#' The function returns: the corrected pooled point estimate (\code{est}) potentially with its sign recoded as indicated by \code{signs.recoded},
#' inference on the bias-corrected estimate (\code{se}, \code{lo}, \code{hi}, \code{pval}), the user's
#' specified \code{eta}, the number of affirmative and nonaffirmative studies after any needed recoding of signs (\code{k.affirmative} and \code{k.nonaffirmative}),
#' and an indicator for whether the point estimates' signs were recoded (\code{signs.recoded}).
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. Preprint available at https://osf.io/s9dp6/.
#' @examples
#'  # calculate effect sizes from example dataset in metafor
#'  require(metafor)
#'  dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#'  # first fit fixed-effects model without any bias correction
#'  # since the point estimate is negative here, we'll assume publication bias favors negative
#'  #  log-RRs rather than positive ones
#'  rma( yi, vi, data = dat, method = "FE" )
#'
#'  # warmup
#'  # note that passing eta = 1 (no publication bias) yields the naive point estimate
#'  #  from rma above, which makes sense
#'  corrected_meta( yi = dat$yi,
#'                  vi = dat$vi,
#'                  eta = 1,
#'                  model = "fixed",
#'                  favor.positive = FALSE )
#'
#'  # assume a known selection ratio of 5
#'  # i.e., affirmative results are 5x more likely to be published
#'  #  than nonaffirmative
#'  corrected_meta( yi = dat$yi,
#'                  vi = dat$vi,
#'                  eta = 5,
#'                  favor.positive = FALSE,
#'                  model = "fixed" )
#'
#'  # same selection ratio, but now account for heterogeneity
#'  # and clustering via robust specification
#'  corrected_meta( yi = dat$yi,
#'                  vi = dat$vi,
#'                  eta = 5,
#'                  favor.positive = FALSE,
#'                  clustervar = dat$author,
#'                  model = "robust" )
#'
#'  ##### Make sensitivity plot as in Mathur & VanderWeele (2020) #####
#'  # range of parameters to try (more dense at the very small ones)
#'  eta.list = as.list( c( 200, 150, 100, 50, 40, 30, 20, rev( seq(1,15,1) ) ) )
#'  res.list = lapply( eta.list, function(x) {
#'                      cat("\n Working on eta = ", x)
#'                      return( corrected_meta( yi = dat$yi,
#'                                              vi = dat$vi,
#'                                              eta = x,
#'                                              model = "robust",
#'                                              favor.positive = FALSE,
#'                                              clustervar = dat$author ) )
#'                                          }
#'                        )
#'
#'  # put results for each eta in a dataframe
#'  res.df = as.data.frame( do.call( "rbind", res.list ) )
#'
#'  require(ggplot2)
#'  ggplot( data = res.df, aes( x = eta, y = est ) ) +
#'
#'    geom_ribbon( data = res.df, aes( x = eta, ymin = lo, ymax = hi ), fill = "gray" ) +
#'
#'    geom_line( lwd = 1.2 ) +
#'    xlab( bquote( eta ) ) +
#'    ylab( bquote( hat(mu)[eta] ) ) +
#'
#'    theme_classic()


corrected_meta = function( yi,
                           vi,
                           eta,
                           clustervar = 1:length(yi),
                           model,
                           selection.tails = 1,
                           favor.positive,
                           alpha.select = 0.05,
                           CI.level = 0.95,
                           small = TRUE ) {

  # stop if eta doesn't make sense
  if ( eta < 1 ) stop( "Eta must be at least 1.")

  # number of point estimates
  k = length(yi)

  # calculate alpha for inference on point estimate
  alpha = 1 - CI.level

  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k & model == "fixed" ) {
    warning( "Clusters exist, but will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }

  ##### Flip Estimate Signs If Needed #####
  # if favor.positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor.positive == TRUE ) {
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
  if ( selection.tails == 1 ) A = (pvals < alpha.select) & (yif > 0)
  if ( selection.tails == 2 ) A = (pvals < alpha.select)

  k.affirmative = sum(A)
  k.nonaffirmative = k - sum(A)

  if ( k.affirmative == 0 | k.nonaffirmative == 0 ) {
    stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
  }

  dat = data.frame( yi, yif, vi, A, clustervar )


  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {

    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat = dat %>% group_by(A) %>%
      summarise( nu = sum( 1 / vi ),
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
      pval.est = 2 * ( 1 - pnorm( z ) )
    }

    # t-based inference
    if ( small == TRUE ) {
      df = k - 1
      lo = est - qt( 1 - (alpha/2), df = df ) * sqrt(var)
      hi = est + qt( 1 - (alpha/2), df = df ) * sqrt(var)
      t =  abs( est / sqrt(var) )
      pval.est = 2 * ( 1 - pt( t, df = df ) )
    }
  } # end fixed = TRUE

  ##### Robust Independent and Robust Clustered #####
  if ( model == "robust" ) {

    # weight for model
    weights = rep( 1, length(pvals) )
    weights[ A == FALSE ] = eta

    # initialize a dumb (unclustered and uncorrected) version of tau^2
    # which is only used for constructing weights
    meta.re = rma.uni( yi = yi,
                       vi = vi)
    t2hat.naive = meta.re$tau2

    # fit weighted robust model
    meta.robu = robu( yi ~ 1,
                      studynum = clustervar,
                      data = dat,
                      userweights = weights / (vi + t2hat.naive),
                      var.eff.size = vi,
                      small = small )

    est = as.numeric(meta.robu$b.r)
    se = meta.robu$reg_table$SE
    lo = meta.robu$reg_table$CI.L
    hi = meta.robu$reg_table$CI.U
    pval.est = meta.robu$reg_table$prob
    eta = eta
  } # end robust = TRUE

  return( data.frame( est,
                      se,
                      lo,
                      hi,
                      pval = pval.est,
                      eta = eta,
                      k.affirmative,
                      k.nonaffirmative,
                      signs.recoded = flipped ) )
}




###############################

#' Severity of publication bias needed to "explain away" results
#'
#' Estimates the S-value, defined as the severity of publication bias (i.e., the ratio
#' by which affirmative studies are more likely to be published than nonaffirmative studies)
#' that would be required to shift the pooled point estimate or its confidence interval limit
#' to the value \code{q}.
#' @param yi A vector of point estimates to be meta-analyzed. Their signs should be coded such that publication bias is
#' assumed to favor positive, rather than negative, estimates.
#' @param vi A vector of estimated variances for the point estimates
#' @param q The attenuated value to which to shift the point estimate or CI. Should be specified on the same scale as \code{yi}
#' (e.g., if \code{yi} is on the log-RR scale, then \code{q} should be as well).
#' @param clustervar A character, factor, or numeric vector with the same length as \code{yi}. Unique values should indicate
#' unique clusters of point estimates. If left unspecified, assumes studies are independent.
#' @param model \code{"fixed"} for fixed-effects (a.k.a. "common effects") or \code{"robust"} for robust random-effects
#' @param alpha.select Alpha-level at which publication probability is assumed to change
#' @param eta.grid.hi The largest value of \code{eta} that should be included in the grid search. This argument is only needed when \code{model = "robust"}.
#' @param favor.positive \code{TRUE} if publication bias is assumed to favor positive estimates; \code{FALSE} if assumed to favor negative estimates.
#' See Details.
#' @param CI.level Confidence interval level (as a proportion) for the corrected point estimate
#' @param small Should inference allow for a small meta-analysis? We recommend using always using \code{TRUE}.
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
#' If \code{favor.positive == FALSE}, such that publication bias is assumed to favor negative rather than positive estimates, the signs of \code{yi} will be reversed prior to
#' performing analyses. The returned number of affirmative and nonaffirmative studies will reflect the recoded signs, and accordingly the returned value \code{signs.recoded} will be \code{TRUE}.
#' @return
#' The function returns: the amount of publication bias required to attenutate the pooled point estimate to \code{q} (\code{sval.est}),
#' the amount of publication bias required to attenutate the confidence interval limit of the pooled point estimate to \code{q} (\code{sval.ci}),
#' the number of affirmative and nonaffirmative studies after any needed recoding of signs (\code{k.affirmative} and \code{k.nonaffirmative}),
#' and an indicator for whether the point estimates' signs were recoded (\code{signs.recoded}).
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. Preprint available at https://osf.io/s9dp6/.
#' @examples
#'  # calculate effect sizes from example dataset in metafor
#'  require(metafor)
#'  dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#'  ##### Fixed-Effects Specification #####
#'  # S-values and worst-case meta-analysis under fixed-effects specification
#'  svals.FE.0 = svalue( yi = dat$yi,
#'                     vi = dat$vi,
#'                     q = 0,
#'                     favor.positive = FALSE,
#'                     model = "fixed" )
#'
#'  # publication bias required to shift point estimate to 0
#'  svals.FE.0$sval.est
#'
#'  # and to shift CI to include 0
#'  svals.FE.0$sval.ci
#'
#'  # now try shifting to a nonzero value (RR = 0.90)
#'  svals.FE.q = svalue( yi = dat$yi,
#'                       vi = dat$vi,
#'                       q = log(.9),
#'                       favor.positive = FALSE,
#'                       model = "fixed" )
#'
#'  # publication bias required to shift point estimate to RR = 0.90
#'  svals.FE.q$sval.est
#'
#'  # and to shift CI to RR = 0.90
#'  svals.FE.q$sval.ci
#'
#'  ##### Robust Clustered Specification #####
#'  svalue( yi = dat$yi,
#'          vi = dat$vi,
#'          q = 0,
#'          favor.positive = FALSE,
#'          model = "robust" )


svalue = function( yi,
                   vi,
                   q,
                   clustervar = 1:length(yi),
                   model,
                   alpha.select = 0.05,
                   eta.grid.hi = 200,
                   favor.positive,
                   CI.level = 0.95,
                   small = TRUE ) {

  # # # ~~~ TEST ONLY
  # # require(metafor)
  # # dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  # dat = sim_data( data.frame( k = 50,
  #                             per.cluster = 5,
  #                             mu = -0.5,
  #                             V = 0.25,
  #                             V.gam = 0,
  #                             sei.min = 0.1,
  #                             sei.max = 1,
  #                             eta = 2 ) )
  # yi = dat$yi
  # vi = dat$vi
  # q = -.1
  # clustervar = 1:length(yi)
  # CI.level = 0.95
  # small = TRUE
  # model = "robust"
  # # # ~~ end test

  # stop if eta doesn't make sense
  if ( eta.grid.hi < 1 ) stop( "eta.grid.hi must be at least 1.")

  # number of point estimates
  k.studies = length(yi)

  alpha = 1 - CI.level

  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k.studies & model == "fixed" ) {
    warning( "You indicated there are clusters, but these will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }

  # fit uncorrected model
  m0 = corrected_meta( yi = yi,
                       vi = vi,
                       eta = 1,
                       model = model,
                       clustervar = clustervar,
                       selection.tails = 1,
                       favor.positive = favor.positive,
                       CI.level = CI.level,
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
  # if favor.positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor.positive == TRUE ) {
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
  A = (pvals < alpha.select) & (yi > 0)

  k.affirmative = sum(A)
  k.nonaffirmative = k.studies - sum(A)

  if ( k.affirmative == 0 | k.nonaffirmative == 0 ) {
    stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
  }

  dat = data.frame( yi, vi, A, clustervar )


  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {

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

      # fit model exactly as in corrected_meta
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
        est.corr = corrected_meta( yi = yi,
                                   vi = vi,
                                   eta = .eta,
                                   model = model,
                                   clustervar = clustervar,
                                   selection.tails = 1,
                                   favor.positive = TRUE,  # always TRUE because we've already flipped signs if needed
                                   CI.level = CI.level,
                                   small = small )$est
        return( abs(est.corr - q))
      }

      opt = optimize( f = func,
                      interval = c(1, eta.grid.hi),
                      maximum = FALSE )
      sval.est = opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta.grid.hi
      if ( abs(sval.est - eta.grid.hi) < 0.0001 & diff > 0.0001 ) sval.est = paste(">", eta.grid.hi)
    }

    # do something similar for CI
    if ( lo.worst > q ) {
      sval.ci = "Not possible"
    } else {
      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func = function(.eta) {
        lo.corr = corrected_meta( yi = yi,
                                  vi = vi,
                                  eta = .eta,
                                  model = model,
                                  clustervar = clustervar,
                                  selection.tails = 1,
                                  favor.positive = TRUE, # always TRUE because we've already flipped signs if needed
                                  CI.level = CI.level,
                                  small = small )$lo
        return( abs(lo.corr - q))
      }

      opt = optimize( f = func,
                      interval = c(1, eta.grid.hi),
                      maximum = FALSE )
      sval.ci = opt$minimum

      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective

      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta.grid.hi
      if ( abs(sval.ci - eta.grid.hi) < 0.0001 & diff > 0.0001 ) sval.ci = paste(">", eta.grid.hi)
    }


    # ##### Worst-case bound #####
    # # as eta -> infinity
    # # not supposed to use KNHA for FE model
    #
    # # flip signs back to original direction if needed
    # if ( flipped == TRUE ) {
    #   yi = -yi
    # }
    #
    # meta.bd = corrected_meta( yi = yi[ A == 0 ],
    #                 vi = vi[ A == 0 ],
    #                 eta = 1,
    #                 model = model,
    #                 CI.level = CI.level,
    #                 small = small )
  }

  # s-values less than 1 indicate complete robustness
  # is.numeric is in case we have a "< XXX" string instead of a number
  if ( is.numeric(sval.est) & !is.na(sval.est) & sval.est < 1) sval.est = "Not possible"
  if ( is.numeric(sval.ci) & !is.na(sval.ci) & sval.ci < 1) sval.ci = "Not possible"

  return( data.frame( sval.est,
                      sval.ci = sval.ci,
                      k.affirmative,
                      k.nonaffirmative,
                      signs.recoded = flipped ) )

}



############################# FN: SIGNIFICANCE FUNNEL PLOT #############################

#' Make significance funnel plot
#'
#' Creates a modified funnel plot that distinguishes between affirmative and nonaffirmative studies, helping to detect the extent to which
#' the nonaffirmative studies' point estimates are systematically smaller than the entire set of point estimates. The estimate among only nonaffirmative studies (gray diamond)
#' represents a corrected estimate under worst-case publication bias. If the gray diamond represents a negligible effect size or if it is much smaller than
#' the pooled estimate among all studies (black diamond), this suggests that the meta-analysis may not be robust to extreme publication bias.
#' Numerical sensitivity analyses (via \code{PublicationBias::svalue}) should still be carried out for more precise quantitative conclusions.
#' @param yi A vector of point estimates to be meta-analyzed.
#' @param vi A vector of estimated variances for the point estimates
#' @param xmin x-axis (point estimate) lower limit for plot
#' @param xmax x-axis (point estimate) upper limit for plot
#' @param ymin y-axis (standard error) lower limit for plot
#' @param ymax y-axis (standard error) upper limit for plot
#' @param xlab Label for x-axis (point estimate)
#' @param ylab Label for y-axis (standard error)
#' @param est.all Regular meta-analytic estimate among all studies (optional)
#' @param est.N Worst-case meta-analytic estimate among only nonaffirmative studies (optional)
#' @param favor.positive \code{TRUE} if publication bias is assumed to favor positive estimates; \code{FALSE} if assumed to favor negative estimates.
#' @param alpha.select Alpha-level at which publication probability is assumed to change
#' @param plot.pooled Should the pooled estimates within all studies and within only the nonaffirmative
#' studies be plotted as well?
#' @import
#' metafor
#' stats
#' ggplot2
#' graphics
#' robumeta
#' @details
#' By default (\code{plot.pooled = TRUE}), also plots the pooled point
#' estimate within all studies, supplied by the user as \code{est.all} (black diamond), and within only the nonaffirmative studies, supplied
#' by the user as \code{est.N} (grey diamond). The user can calculate \code{est.all} and \code{est.N} using their choice of meta-analysis model. If instead
#' these are not supplied but \code{plot.pooled = TRUE}, these pooled estimates will be automatically calculated using a fixed-effects (a.k.a. "common-effects") model.
#' @export
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. Preprint available at https://osf.io/s9dp6/.
#' @examples
#'
#' ##### Make Significance Funnel with User-Specified Pooled Estimates #####
#'
#' # compute meta-analytic effect sizes for an example dataset
#' require(metafor)
#' dat = metafor::escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#'
#' # flip signs since we think publication bias operates in favor of negative effects
#' # alternatively, if not flipping signs, could pass favor.positive = FALSE to
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
#'                      favor.positive = TRUE,
#'                      alpha.select = 0.50,
#'                      plot.pooled = TRUE )



significance_funnel = function( yi,
                                vi,
                                xmin = min(yi),
                                xmax = max(yi),
                                ymin = 0,  # so that pooled points are shown
                                ymax = max( sqrt(vi) ),
                                xlab = "Point estimate",
                                ylab = "Estimated standard error",
                                favor.positive = NA,
                                est.all = NA,
                                est.N = NA,
                                alpha.select = 0.05,
                                plot.pooled = TRUE ) {

  d = data.frame(yi, vi)
  d$sei = sqrt(vi)

  # calculate p-values
  d$pval = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )

  # which direction of effects are favored?
  # if we have the pooled point estimate, but not the favored direction,
  #  assume favored direction matches sign of pooled estimate (but issue warning)
  if ( !is.na(est.all) & is.na(favor.positive) ) {
    favor.positive = (est.all > 0)
    warning("favor.positive not provided, so assuming publication bias favors estimates whose sign matches est.all")
  }
  if ( is.na(est.all) & is.na(favor.positive) ) {
    stop("Need to specify favor.positive")
  }

  # affirmative vs. nonaffirmative indicator
  d$affirm = rep(NA, nrow(d))

  if ( favor.positive == TRUE ) {
    d$affirm[ (d$yi > 0) & (d$pval < alpha.select) ] = "Affirmative"
    d$affirm[ (d$yi < 0) | (d$pval >= alpha.select) ] = "Non-affirmative"
  }
  if ( favor.positive == FALSE ) {
    d$affirm[ (d$yi < 0) & (d$pval < alpha.select) ] = "Affirmative"
    d$affirm[ (d$yi > 0) | (d$pval >= alpha.select) ] = "Non-affirmative"
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
  # if not supplied, gets them from common-effects model
  if ( is.na(est.N) & is.na(est.all) ) {
    est.N = rma.uni(yi = d$yi[ d$affirm == "Non-affirmative" ],
                    vi = d$vi[ d$affirm == "Non-affirmative" ],
                    method="FE")$b

    est.all = rma.uni(yi = d$yi,
                      vi = d$vi,
                      method="FE")$b
  }

  # set up pooled estimates for plotting
  pooled.pts = data.frame( yi = c(est.N, est.all),
                           sei = c(0,0) )

  # for a given SE (y-value), return the "just significant" point estimate value (x-value)
  just_signif_est = function( .sei ) .sei * qnorm(1 - alpha.select/2)

  # calculate slope and intercept of the "just affirmative" line
  # i.e., 1.96 = (just affirmative estimate) / se
  if (favor.positive == TRUE) sl = 1/qnorm(1 - alpha.select/2)
  if (favor.positive == FALSE) sl = -1/qnorm(1 - alpha.select/2)
  int = 0
  # # sanity check: should be exactly alpha.select
  # 2 * ( 1 - pnorm( abs(1) / sl ) )


  ##### Make the Plot #####
  colors = c("darkgray", "orange")

  p.funnel = ggplot( data = d, aes( x = d$yi,
                                    y = d$sei,
                                    color = d$affirm ) )

  if ( plot.pooled == TRUE ) {

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



############################# FN: P-VALUE PLOT #############################

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
#' @param alpha.select Alpha-level at which publication probability is assumed to change
#' @import
#' stats
#' ggplot2
#' @export
#' @references
#' 1. Mathur MB & VanderWeele TJ (2020). Sensitivity analysis for publication bias in meta-analyses. Preprint available at https://osf.io/s9dp6/.
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
                      alpha.select = 0.05) {

  # calculate 1-tailed p-values
  pval = 1 - pnorm( yi / sqrt(vi) )

  ggplot( data = data.frame(pval = pval),
          aes( x = pval ) ) +
    geom_vline(xintercept = alpha.select/2, color = "red", lwd = 1) +
    geom_vline(xintercept = 1 - (alpha.select/2), color = "red", lwd = 1) +
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

