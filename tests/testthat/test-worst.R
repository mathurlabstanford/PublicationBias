test_that("worst-case meta fixed", {
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                data = dat.bcg)
  pubbias <- pubbias_meta(dat$yi, dat$vi, model_type = "fixed",
                          selection_ratio = 1, favor_positive = FALSE,
                          return_worst_meta = TRUE)
  worst <- fit_meta_worst(pubbias$data, model_type = "fixed",
                          ci_level = 0.95, small = TRUE)

  expect_equal(worst$stats, pubbias$stats |> filter(model == "worst_case"))
})

test_that("worst-case meta robust", {
  dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                data = dat.bcg)
  pubbias <- pubbias_meta(dat$yi, dat$vi, cluster = , dat$author,
                          model_type = "robust", selection_ratio = 1,
                          favor_positive = FALSE, return_worst_meta = TRUE)
  worst <- fit_meta_worst(pubbias$data, model_type = "robust",
                          ci_level = 0.95, small = TRUE)

  expect_equal(worst$stats, pubbias$stats |> filter(model == "worst_case"))
})
