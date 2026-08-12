test_that("fit_baygmst validates inputs before touching CmdStan", {
  # These should all fail on argument validation, before fit_baygmst() ever
  # calls instantiate::stan_package_model() -- so they run even in
  # environments without CmdStan installed.
  expect_error(
    fit_baygmst(
      proxy = 1:5, instrumental_T = 1:4, forcing_G = 1:5,
      forcing_V = 1:5, forcing_S = 1:5, years = 1:5
    ),
    regexp = NULL
  )
  expect_error(
    fit_baygmst(
      proxy = 1:5, instrumental_T = c(1, NA, NA, NA, 5),
      forcing_G = 1:5, forcing_V = 1:5, forcing_S = 1:5, years = 1:5,
      iter_sampling = 10
    ),
    "iter_sampling"
  )
  expect_error(
    fit_baygmst(
      proxy = rep(NA_real_, 5), instrumental_T = c(1, NA, NA, NA, 5),
      forcing_G = 1:5, forcing_V = 1:5, forcing_S = 1:5, years = 1:5
    )
  )
})

test_that("fit_baygmst fits and reconstruct()/summary() work end to end", {
  skip_if_not_installed("instantiate")
  skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not available")

  set.seed(1)
  n <- 60
  years <- seq_len(n)
  true_T <- cumsum(rnorm(n, 0, 0.05))
  instrumental_T <- c(rep(NA_real_, 40), true_T[41:60])
  proxy <- true_T + rnorm(n, 0, 0.1)
  forcing <- transform_forcings(
    G = 280 + cumsum(rgamma(n, 0.05, 1)),
    V = abs(rnorm(n, 0, 0.05)),
    S = 1361 + rnorm(n, 0, 0.3)
  )

  fit <- fit_baygmst(
    proxy = proxy, instrumental_T = instrumental_T,
    forcing_G = forcing$G, forcing_V = forcing$V, forcing_S = forcing$S,
    years = years, chains = 1, parallel_chains = 1,
    iter_warmup = 200, iter_sampling = 1000, seed = 1
  )
  expect_s3_class(fit, "baygmst_fit")

  recon <- reconstruct(fit)
  expect_equal(nrow(recon), n)
  expect_setequal(unique(recon$type), c("instrumental", "reconstruction"))

  s <- summary(fit)
  expect_s3_class(s, "summary.baygmst_fit")
  expect_true(all(c("RMSE", "MAE", "Bias", "Correlation", "R2") %in% names(s$performance)))
})
