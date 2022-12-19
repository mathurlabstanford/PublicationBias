# does pubbias_meta agree with regular meta-analysis fns when there
#  is no selection?
test_that("pubbias_meta no selection agreement", {
  tol <- 0.001
  ##### Recover Regular FE model With selection_ratio = 1 #####
  # when using z-based inference and selection_ratio = 1,
  # should match metafor
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                data = dat.bcg)
  fe_plain <- rma(yi, vi, data = dat, method = "FE")
  # should match corrected with selection_ratio = 1

  fe_adj <- pubbias_meta(yi = dat$yi,
                         vi = dat$vi,
                         selection_ratio = 1,
                         model = "fixed",
                         selection_tails = 1,
                         ci_level = 0.95,
                         small = FALSE,
                         favor_positive = FALSE)$stats

  expect_equal(fe_adj$estimate, as.numeric(fe_plain$b), tolerance = tol)
  expect_equal(fe_adj$se, as.numeric(fe_plain$se), tolerance = tol)
  expect_equal(fe_adj$ci_lower, as.numeric(fe_plain$ci.lb), tolerance = tol)
  expect_equal(fe_adj$ci_upper, as.numeric(fe_plain$ci.ub), tolerance = tol)
  expect_equal(fe_adj$p_value, as.numeric(fe_plain$pval), tolerance = tol)

  ##### Recover Regular Robust Indepenent Model With selection_ratio = 1 #####
  for (.small in c(TRUE, FALSE)) {
    cluster <- 1:length(dat$yi)

    meta_re <- rma.uni(yi = dat$yi, vi = dat$vi)
    t2hat_naive <- meta_re$tau2

    ri_robust <- robumeta::robu(yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = 1 / (vi + t2hat_naive),
                                var.eff.size = vi,
                                small = .small)


    ri_adj <- pubbias_meta(yi = dat$yi,
                           vi = dat$vi,
                           selection_ratio = 1,
                           model = "robust",
                           selection_tails = 1,
                           ci_level = 0.95,
                           small = .small,
                           favor_positive = FALSE)$stats

    rt <- ri_robust$reg_table
    expect_equal(ri_adj$estimate, as.numeric(ri_robust$b.r), tolerance = tol)
    expect_equal(ri_adj$se, as.numeric(rt$SE), tolerance = tol)
    expect_equal(ri_adj$ci_lower, as.numeric(rt$CI.L), tolerance = tol)
    expect_equal(ri_adj$ci_upper, as.numeric(rt$CI.U), tolerance = tol)
    expect_equal(ri_adj$p_value, as.numeric(rt$prob), tolerance = tol)
  }


  ##### Recover Regular Robust Clustered Model With selection_ratio = 1 #####
  for (.small in c(TRUE, FALSE)) {
    cluster <- dat$author

    meta_re <- rma.uni(yi = dat$yi, vi = dat$vi)
    t2hat_naive <- meta_re$tau2

    ri_robust <- robumeta::robu(yi ~ 1,
                                studynum = cluster,
                                data = dat,
                                userweights = 1 / (vi + t2hat_naive),
                                var.eff.size = vi,
                                small = .small)

    ri_adj <- pubbias_meta(yi = dat$yi,
                           vi = dat$vi,
                           selection_ratio = 1,
                           model = "robust",
                           cluster = cluster,
                           selection_tails = 1,
                           ci_level = 0.95,
                           small = .small,
                           favor_positive = FALSE)$stats

    rt <- ri_robust$reg_table
    expect_equal(ri_adj$estimate, as.numeric(ri_robust$b.r), tolerance = tol)
    expect_equal(ri_adj$se, as.numeric(rt$SE), tolerance = tol)
    expect_equal(ri_adj$ci_lower, as.numeric(rt$CI.L), tolerance = tol)
    expect_equal(ri_adj$ci_upper, as.numeric(rt$CI.U), tolerance = tol)
    expect_equal(ri_adj$p_value, as.numeric(rt$prob), tolerance = tol)
  }

})
