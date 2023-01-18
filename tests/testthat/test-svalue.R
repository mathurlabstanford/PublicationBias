dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg)

# do pubbias_svalue and pubbias_meta agree?
test_that("pubbias_svalue / pubbias_meta agreement", {

  for (.ci_level in c(.8, .95)) {
    for (.small in c(TRUE, FALSE)) {

      svals <- pubbias_svalue(yi = dat$yi,
                              vi = dat$vi,
                              q = 0,
                              model = "fixed",
                              ci_level = .ci_level,
                              small = .small,
                              favor_positive = FALSE)

      # CI upper limit should be exactly 0 when selection_ratio = sval_ci
      meta <- pubbias_meta(yi = dat$yi,
                           vi = dat$vi,
                           selection_ratio = as.numeric(svals$stats$sval_ci),
                           model = "fixed",
                           selection_tails = 1,
                           ci_level = .ci_level,
                           small = .small,
                           favor_positive = FALSE)$stats

      expect_equal(meta$ci_upper, 0)

    }
  }
})


# does it reject bad choices of q?
# i.e., the point estimate is already closer to the null than q
test_that("pubbias_svalue q outside point estimate", {
  expect_error(
    regexp = "q must be greater",
    pubbias_svalue(yi = dat$yi,
                   vi = dat$vi,
                   q = -3,
                   model = "fixed",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = FALSE)
  )

  # reverse signs
  expect_error(
    regexp = "q must be less",
    pubbias_svalue(yi = -dat$yi,
                   vi = dat$vi,
                   q = 0.8,
                   model = "fixed",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = TRUE)
  )
})


# does it produce appropriate warnings when q is already inside the CI?
test_that("pubbias_svalue q inside CI", {
  q_inside_msg <- "naive confidence interval already contains q"

  # fixed case
  # calculated naive: -0.4302852 [-0.5096613, -0.3509091]
  # should give message about sval.ci not applying
  expect_message(
    regexp = q_inside_msg,
    pubbias_svalue(yi = dat$yi,
                   vi = dat$vi,
                   q = -0.4,  # closer to null than naive estimate, within CI
                   model = "fixed",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = FALSE)
  )

  # should run without message/error
  expect_message(
    regexp = NA,
    pubbias_svalue(yi = dat$yi,
                   vi = dat$vi,
                   q = -0.3,  # closer to null than naive estimate, within CI
                   model = "fixed",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = FALSE)
  )

  # robust case and flipped signs
  # calculated naive: 0.7145323 [0.3241296, 1.104935]
  expect_message(
    regexp = q_inside_msg,
    pubbias_svalue(yi = -dat$yi,
                   vi = dat$vi,
                   q = 0.6,  # closer to null than naive estimate, within CI
                   model = "robust",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = TRUE)
  )

  # should run without message or error
  expect_message(
    regexp = NA,
    pubbias_svalue(yi = -dat$yi,
                   vi = dat$vi,
                   q = 0.3,  # closer to null than naive estimate, within CI
                   model = "robust",
                   ci_level = 0.95,
                   small = FALSE,
                   favor_positive = TRUE))
})



# does pubbias_svalue give correct results when the s-value is greater than the
# highest value in selection_ratio grid?
test_that("pubbias_svalue s-value outside grid", {

  selection_ratio <- 3
  svals <- pubbias_svalue(yi = dat$yi,
                          vi = dat$vi,
                          q = 0,
                          selection_ratio_max = 2,
                          model = "robust",
                          ci_level = 0.95,
                          small = TRUE,
                          favor_positive = FALSE)
  expect_equal(svals$stats$sval_ci, "> 2")
})
