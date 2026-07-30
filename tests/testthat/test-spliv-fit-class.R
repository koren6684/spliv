test_that("spliv fits use only the new class", {
  d <- make_synth_panel(n_gid = 8, n_t = 8, seed = 701)
  prior <- sp_prior_ltz(y ~ x + w1 + w2 | z + w1 + w2, d,
                        inst_vary = "z", mean = 0, sd = 0.1)
  expect_silent(
    fit <- spliv(y ~ x + w1 + w2 | z + w1 + w2, d, method = "ltz",
                 prior = prior, vcov = "hc1")
  )

  expect_s3_class(fit, "spliv_fit")
  expect_identical(class(fit), "spliv_fit")

  expect_output(print(fit), "spliv fit")

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  expect_silent(plot(fit))
  expect_identical(
    utils::getS3method("print", "spliv_fit"),
    getFromNamespace("print.spliv_fit", "spliv")
  )
  expect_identical(
    utils::getS3method("plot", "spliv_fit"),
    getFromNamespace("plot.spliv_fit", "spliv")
  )
})
