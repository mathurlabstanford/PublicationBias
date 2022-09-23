library(testthat)
library(devtools)
library(metafor)


# helper fn for simulating meta-analysis with publication bias

# # p: row of parameters dataframe
sim_data <- function(p) {
  withr::with_seed(42, {

    N = p$k * p$per.cluster

    # generate cluster random intercepts
    gam1 = rnorm( n = p$k, mean = 0, sd = sqrt( p$V.gam ) )
    gam1i = rep( gam1, each = p$per.cluster )

    # generate individual-study random intercepts
    gam2i = rnorm( n = N, mean = 0, sd = sqrt( p$V - p$V.gam ) )

    # individual study means
    mui = p$mu + gam1i + gam2i
    sei = runif( n = N, min = p$sei.min, max = p$sei.max )
    yi = rnorm( n = N, mean = mui, sd = sei )

    d = data.frame( cluster = rep(1:p$k, each = p$per.cluster),
                    Study.name = 1:N,
                    yi = yi,
                    sei = sei,
                    vi = sei^2,
                    pval = 2 * ( 1 - pnorm( abs(yi) / sei ) ) )

    # 1-tailed publication bias
    signif = d$pval < 0.05 & d$yi > 0
    publish = rep( 1, nrow(d) )
    publish[ signif == FALSE ] = rbinom( n = sum(signif == FALSE), size = 1, prob = 1/p$selection_ratio )

    d$weight = 1
    d$weight[ signif == 0 ] = p$selection_ratio
    d = d[ publish == 1, ]

    return(d)
  })
}

# ##### Sanity Check #####
# d = sim_data( data.frame( k = 5,
#                           per.cluster = 20,
#                           mu = .5,
#                           V = 1,
#                           V.gam = .5,
#                           sei.min = 0.01,
#                           sei.max = 0.01,
#                           selection_ratio = 1 ) )


# does pubbias_meta agree with regular meta-analysis fns when there
#  is no selection?
test_that("pubbias_meta #1", {

  ##### Recover Regular FE model With selection_ratio = 1 #####
  # when using z-based inference and selection_ratio = 1,
  # should match metafor
  dat = escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  FE.plain = rma( yi, vi, data = dat, method = "FE" )
  # should match corrected with selection_ratio = 1

  FE.adj = pubbias_meta( yi = dat$yi,
                             vi = dat$vi,
                             selection_ratio = 1,
                             model = "fixed",
                             selection_tails = 1,
                             ci_level = 0.95,
                             small = FALSE,
                           favor_positive = FALSE )$stats

  expect_equal( FE.adj$estimate, as.numeric(FE.plain$b), tol = 0.001 )
  expect_equal( FE.adj$se, as.numeric(FE.plain$se), tol = 0.001 )
  expect_equal( FE.adj$ci_lower, as.numeric(FE.plain$ci.lb), tol = 0.001 )
  expect_equal( FE.adj$ci_upper, as.numeric(FE.plain$ci.ub), tol = 0.001 )
  expect_equal( FE.adj$p_value, as.numeric(FE.plain$pval), tol = 0.001 )

  ##### Recover Regular Robust Indepenent Model With selection_ratio = 1 #####
  for ( .small in c(TRUE, FALSE) ) {
    cluster = 1:length(dat$yi)

    meta.re = rma.uni( yi = dat$yi,
                       vi = dat$vi)
    t2hat.naive = meta.re$tau2

    RI.robust = robumeta::robu( yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = 1 / (vi + t2hat.naive),
                                var.eff.size = vi,
                                small = .small )


    RI.adj = pubbias_meta( yi = dat$yi,
                             vi = dat$vi,
                             selection_ratio = 1,
                             model = "robust",
                             selection_tails = 1,
                             ci_level = 0.95,
                             small = .small,
                             favor_positive = FALSE)$stats

    expect_equal( RI.adj$estimate, as.numeric( as.numeric(RI.robust$b.r) ), tol = 0.001 )
    expect_equal( RI.adj$se, as.numeric( RI.robust$reg_table$SE ), tol = 0.001 )
    expect_equal( RI.adj$ci_lower, as.numeric( RI.robust$reg_table$CI.L ), tol = 0.001 )
    expect_equal( RI.adj$ci_upper, as.numeric( RI.robust$reg_table$CI.U ), tol = 0.001 )
    expect_equal( RI.adj$p_value, as.numeric( RI.robust$reg_table$prob ), tol = 0.001 )
  }


  ##### Recover Regular Robust Clustered Model With selection_ratio = 1 #####
  for ( .small in c(TRUE, FALSE) ) {
    cluster = dat$author

    meta.re = rma.uni( yi = dat$yi,
                       vi = dat$vi)
    t2hat.naive = meta.re$tau2

    RI.robust = robumeta::robu( yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = 1 / (vi + t2hat.naive),
                                var.eff.size = vi,
                                small = .small )


    RI.adj = pubbias_meta( yi = dat$yi,
                             vi = dat$vi,
                             selection_ratio = 1,
                             model = "robust",
                             cluster = cluster,
                             selection_tails = 1,
                             ci_level = 0.95,
                             small = .small,
                             favor_positive = FALSE)$stats

    expect_equal( RI.adj$estimate, as.numeric( as.numeric(RI.robust$b.r) ), tol = 0.001 )
    expect_equal( RI.adj$se, as.numeric( RI.robust$reg_table$SE ), tol = 0.001 )
    expect_equal( RI.adj$ci_lower, as.numeric( RI.robust$reg_table$CI.L ), tol = 0.001 )
    expect_equal( RI.adj$ci_upper, as.numeric( RI.robust$reg_table$CI.U ), tol = 0.001 )
    expect_equal( RI.adj$p_value, as.numeric( RI.robust$reg_table$prob ), tol = 0.001 )
  }

})


# do pubbias_svalue and pubbias_meta agree?
test_that("pubbias_svalue #1", {

  dat = escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

  for ( .ci_level in c(.8, .95) ) {
    for ( .small in c(TRUE, FALSE) ) {

      svals = pubbias_svalue( yi = dat$yi,
                      vi = dat$vi,
                      q = 0,
                      model = "fixed",
                      ci_level = 0.95,
                      small = .small,
                      favor_positive = FALSE )

      # CI upper limit should be exactly 0 when selection_ratio = sval_ci
      meta = pubbias_meta( yi = dat$yi,
                             vi = dat$vi,
                             selection_ratio = as.numeric(svals$stats$sval_ci),
                             model = "fixed",
                             selection_tails = 1,
                             ci_level = 0.95,
                             small = .small,
                             favor_positive = FALSE )$stats

      expect_equal( meta$ci_upper, 0 )

    }
  }
} )




# # is worst-case meta correct for both 1-tailed and 2-tailed selection?
# test_that("pubbias_svalue #2", {
#
#   dat = escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
#   # doctor one data point so that 2-tailed will be different from 1-tailed
#   dat$yi[2] = -dat$yi[2]
#
#   # first get worst-case meta-analysis using package
#   m2 = pubbias_svalue( yi = dat$yi,
#                vi = dat$vi,
#                q = 0,
#                model = "fixed",
#                ci_level = 0.95,
#                small = FALSE )$meta.bd
#
#   # flip signs
#   rma.uni( dat$yi,
#            dat$vi,
#            method = "FE" )$b
#   dat$yi = -dat$yi
#
#   # calculate p-values
#   Z = abs( dat$yi / sqrt(dat$vi) )
#   pval = 2 * ( 1 - pnorm(Z) )
#
#   # affirmative indicators
#   A.1tail = (pval < .05) & (dat$yi > 0)
#   A.2tail = (pval < 0.05)
#
#   table(A.1tail)
#   table(A.2tail)
#
#   # note sign flip back to original
#   m1 = rma.uni( -dat$yi[ A.1tail == FALSE ],
#                  dat$vi[ A.1tail == FALSE ],
#                  method = "FE" )
#
#   expect_equal( as.numeric( m1$b ),
#                 as.numeric( m2$b ) )
# } )



# does it reject bad choices of q?
# i.e., the point estimate is already closer to the null than q
test_that( "pubbias_svalue #3", {
  dat = escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

  expect_error( pubbias_svalue( yi = dat$yi,
          vi = dat$vi,
          q = -3,
          model = "fixed",
          ci_level = 0.95,
          small = FALSE,
          favor_positive = FALSE) )

  # reverse signs
  expect_error( pubbias_svalue( yi = -dat$yi,
                        vi = dat$vi,
                        q = .8,
                        model = "fixed",
                        ci_level = 0.95,
                        small = FALSE,
                        favor_positive = FALSE) )
} )


# does it produce appropriate warnings when q is already inside the CI?
test_that("pubbias_svalue #3.5", {
  dat = escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

  # fixed case
  # naive: -0.4302852 [-0.5096613, -0.3509091]
  # should give message about sval.ci not applying
  expect_message( pubbias_svalue( yi = dat$yi,
                          vi = dat$vi,
                          q = -0.4,  # closer to null than naive estimate, but within CI (from being in browser and looking at m0)
                          model = "fixed",
                          ci_level = 0.95,
                          small = FALSE,
                          favor_positive = FALSE) )

  # should run without message/error
  pubbias_svalue( yi = dat$yi,
          vi = dat$vi,
          q = -0.3,  # closer to null than naive estimate, but within CI (from being in browser and looking at m0)
          model = "fixed",
          ci_level = 0.95,
          small = FALSE,
          favor_positive = FALSE)

  # robust case and flipped signs
  # naive: 0.7145323 [0.3241296, 1.104935]
  expect_message( pubbias_svalue( yi = -dat$yi,
                          vi = dat$vi,
                          q = 0.6,  # closer to null than naive estimate, but within CI
                          model = "robust",
                          ci_level = 0.95,
                          small = FALSE,
                          favor_positive = TRUE) )

  # should run without message or error
  pubbias_svalue( yi = -dat$yi,
                          vi = dat$vi,
                          q = 0.3,  # closer to null than naive estimate, but within CI
                          model = "robust",
                          ci_level = 0.95,
                          small = FALSE,
                          favor_positive = TRUE)
})



# does pubbias_svalue give correct results when the s-value is greater than the highest
#  value in selection_ratio grid?
test_that( "pubbias_svalue #4", {

  selection_ratio = 3

  d = sim_data( data.frame( k = 50,
                            per.cluster = 5,
                            mu = 0.5,
                            V = 0,
                            V.gam = 0,
                            sei.min = 0.1,
                            sei.max = 1,
                            selection_ratio = selection_ratio ) )

  svals = pubbias_svalue( yi = d$yi,
          vi = d$vi,
          q = 0,
          model = "robust",
          selection_ratio_grid_hi = 10,
          # selection_ratio_grid = seq(1,10,1),
          ci_level = 0.95,
          small = FALSE,
          favor_positive = FALSE)


} )


# does significance_funnel give error if no non-affirmative studies?
test_that( "significance_funnel #1", {
  selection_ratio = 50

  d = sim_data( data.frame( k = 50,
                            per.cluster = 1,
                            mu = 0.5,
                            V = 0.1,
                            V.gam = 0,
                            # sei.min = 0.2,
                            # sei.max = .9,
                            sei.min = .1,
                            sei.max=.1,
                            selection_ratio = selection_ratio ) )


  # also see if significance_funnel works
  expect_error( significance_funnel( yi = d$yi,
                       vi = d$vi,
                       favor_positive = FALSE) )
} )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                     MANUAL TESTS                                    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# ##### Visually Test Significance Funnel on Simulated Data #####
# selection_ratio = 5
#
# d = sim_data( data.frame( k = 500,
#                           per.cluster = 1,
#                           mu = -0.1,
#                           V = 0.1,
#                           V.gam = 0,
#                           sei.min = 0.01,
#                           sei.max = 0.05,
#                           selection_ratio = selection_ratio ) )
#
#
# # also see if significance_funnel works
# significance_funnel( yi = d$yi,
#                      vi = d$vi )


# ##### Reproduce Manuscript Stats for Boehm Data #####
#
# setwd("~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Applied examples/Data/Boehm data")
#
# d = read.csv("boehm_prepped.csv")
#
# significance_funnel( yi = d$yi,
#                      vi = d$vi )
#
# ### fixed effects
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = 0,
#         model = "fixed",
#         small = TRUE)
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = r_to_z(0.1),
#         model = "fixed",
#         small = TRUE)
#
#
# ### robust independent
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = 0,
#         model = "robust" )
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = r_to_z(0.1),
#         model = "robust" )
#
# ### robust clusters
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         cluster = d$study,
#         q = 0,
#         model = "robust" )
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         cluster = d$study,
#         q = r_to_z(0.1),
#         model = "robust" )
#
#
#
# ##### Reproduce Manuscript Stats for Anderson Data #####
#
# setwd("~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Private data component/Anderson data")
#
# d = read.csv("anderson_prepped.csv")
#
# significance_funnel( yi = d$yi,
#                      vi = d$vi )
#
# ### fixed effects
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = 0,
#         model = "fixed",
#         small = TRUE)
# # matches :)
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = r_to_z(0.1),
#         model = "fixed",
#         small = TRUE)
#
#
# ### robust independent
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = 0,
#         model = "robust" )
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         q = r_to_z(0.1),
#         model = "robust" )
# # should move CI to 0.10
# pubbias_meta( yi = d$yi,
#                 vi = d$vi,
#                 selection_ratio = 5.25,
#                 model = "robust" )
# # yes :)
#
# ### robust clusters
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         cluster = d$cluster,
#         q = 0,
#         model = "robust" )
#
# pubbias_svalue( yi = d$yi,
#         vi = d$vi,
#         cluster = d$cluster,
#         q = r_to_z(0.1),
#         model = "robust" )
# # should move CI to 0.10
# pubbias_meta( yi = d$yi,
#                 vi = d$vi,
#                 selection_ratio = 3.5,
#                 model = "robust",
#                 cluster = d$cluster )






