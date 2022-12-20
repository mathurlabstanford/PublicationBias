test_that("correct output structure for pubbias_meta", {

  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                         di = cneg, data = dat.bcg)

  # change values of all params and opts
  params <- list(selection_ratio = 1, model_type = "robust",
                 selection_tails = 2, favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE)

  meta <- rlang::exec(pubbias_meta, yi = dat$yi, vi = dat$vi,
                      cluster = dat$author, !!!params)

  expect_s3_class(meta, "metabias")
  expect_named(meta, c("data", "values", "stats", "fits"))

  expect_s3_class(meta$data, "data.frame")
  expect_equal(nrow(meta$data), nrow(dat))
  expect_named(meta$data, meta_names("data"))
  expect_equal(meta$data$yi, dat$yi)
  expect_equal(meta$data$yif, -dat$yi)
  expect_equal(meta$data$vi, dat$vi)
  expect_equal(meta$data$cluster, dat$author)

  expect_type(meta$values, "list")
  expect_named(meta$values, meta_names("values"))
  purrr::walk(names(params), \(p) expect_equal(params[[p]], meta$values[[p]]))

  expect_s3_class(meta$stats, "data.frame")
  expect_equal(nrow(meta$stats), 1)
  expect_named(meta$stats, meta_names("stats"))
  expect_true(all(meta$stats$se < 1))
  expect_true(all(meta$stats$ci_lower < meta$stats$estimate))
  expect_true(all(meta$stats$ci_upper > meta$stats$estimate))
  expect_true(all(meta$stats$p_value < 1))

  expect_length(meta$fit, 1)

})

test_that("correct output structure for pubbias_svalue", {

  dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                         di = cneg, data = dat.bcg)

  # change values of all params and opts
  params <- list(q = 1, model_type = "robust", favor_positive = FALSE,
                 alpha_select = 0.01, ci_level = 0.90, small = FALSE)

  svalue <- rlang::exec(pubbias_svalue, yi = dat$yi, vi = dat$vi,
                        cluster = dat$author, !!!params,
                        return_worst_meta = TRUE)

  expect_s3_class(svalue, "metabias")
  expect_named(svalue, c("data", "values", "stats", "fits"))

  expect_s3_class(svalue$data, "data.frame")
  expect_equal(nrow(svalue$data), nrow(dat))
  expect_named(svalue$data, svalue_names("data"))
  expect_equal(svalue$data$yi, -dat$yi)
  expect_equal(svalue$data$vi, dat$vi)
  expect_equal(svalue$data$cluster, dat$author)

  expect_type(svalue$values, "list")
  expect_named(svalue$values, svalue_names("values"))
  params$q <- -params$q
  purrr::walk(names(params), \(p) expect_equal(params[[p]], svalue$values[[p]]))

  expect_s3_class(svalue$stats, "data.frame")
  expect_named(svalue$stats, svalue_names("stats"))

  expect_length(svalue$fit, 1)

})
