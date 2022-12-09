# does pubbias_meta agree with regular meta-analysis fns when there
#  is no selection?
test_that("pubbias_meta no selection agreement", {

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

  expect_equal( FE.adj$estimate, as.numeric(FE.plain$b), tolerance = 0.001 )
  expect_equal( FE.adj$se, as.numeric(FE.plain$se), tolerance = 0.001 )
  expect_equal( FE.adj$ci_lower, as.numeric(FE.plain$ci.lb), tolerance = 0.001 )
  expect_equal( FE.adj$ci_upper, as.numeric(FE.plain$ci.ub), tolerance = 0.001 )
  expect_equal( FE.adj$p_value, as.numeric(FE.plain$pval), tolerance = 0.001 )

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

    expect_equal( RI.adj$estimate, as.numeric( as.numeric(RI.robust$b.r) ), tolerance = 0.001 )
    expect_equal( RI.adj$se, as.numeric( RI.robust$reg_table$SE ), tolerance = 0.001 )
    expect_equal( RI.adj$ci_lower, as.numeric( RI.robust$reg_table$CI.L ), tolerance = 0.001 )
    expect_equal( RI.adj$ci_upper, as.numeric( RI.robust$reg_table$CI.U ), tolerance = 0.001 )
    expect_equal( RI.adj$p_value, as.numeric( RI.robust$reg_table$prob ), tolerance = 0.001 )
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

    expect_equal( RI.adj$estimate, as.numeric( as.numeric(RI.robust$b.r) ), tolerance = 0.001 )
    expect_equal( RI.adj$se, as.numeric( RI.robust$reg_table$SE ), tolerance = 0.001 )
    expect_equal( RI.adj$ci_lower, as.numeric( RI.robust$reg_table$CI.L ), tolerance = 0.001 )
    expect_equal( RI.adj$ci_upper, as.numeric( RI.robust$reg_table$CI.U ), tolerance = 0.001 )
    expect_equal( RI.adj$p_value, as.numeric( RI.robust$reg_table$prob ), tolerance = 0.001 )
  }

})
