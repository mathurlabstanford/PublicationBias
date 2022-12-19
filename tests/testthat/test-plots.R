# does significance_funnel give error if no non-affirmative studies?
test_that("significance_funnel no non-affirmative", {

  d <- sim_data(data.frame(k = 50,
                           per.cluster = 1,
                           mu = 0.5,
                           V = 0.1,
                           V.gam = 0,
                           sei.min = .1,
                           sei.max = .1,
                           selection_ratio = 50))


  # also see if significance_funnel works
  expect_error(
    significance_funnel(yi = d$yi, vi = d$vi, favor_positive = FALSE)
  )
})
