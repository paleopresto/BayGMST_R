# CmdStan-free tests of the post-fitting pipeline (reconstruct, summary,
# print, and the plot functions), using a mocked "baygmst_fit" whose $fit
# element imitates the small slice of the cmdstanr::CmdStanMCMC interface
# these functions actually use.

make_mock_fit <- function(n_obs = 20, n_mis = 30, n_draws = 40, seed = 42) {
  set.seed(seed)
  years   <- seq_len(n_obs + n_mis)
  idx_mis <- seq_len(n_mis)          # contiguous missing block first,
  idx_obs <- n_mis + seq_len(n_obs)  # observed block after (as in the data)
  y_obs   <- rnorm(n_obs)

  y_mis_draws <- matrix(
    rnorm(n_draws * n_mis), nrow = n_draws,
    dimnames = list(NULL, paste0("y_mis[", seq_len(n_mis), "]"))
  )
  y_ins_draws <- matrix(
    rnorm(n_draws * n_obs), nrow = n_draws,
    dimnames = list(NULL, paste0("y_ins_fitted[", seq_len(n_obs), "]"))
  )
  all_draws <- cbind(y_mis_draws, y_ins_draws)

  mock_cmdstan <- list(
    # reconstruct(): fit$summary(variables = "y_ins_fitted", ~..., "mean")
    # summary.baygmst_fit(): fit$summary(variables = <parameters>)
    summary = function(variables, ...) {
      if (identical(variables, "y_ins_fitted")) {
        data.frame(
          variable = paste0("y_ins_fitted[", seq_len(n_obs), "]"),
          lo_inner = y_obs - 0.5, hi_inner = y_obs + 0.5,
          lo_outer = y_obs - 1, hi_outer = y_obs + 1,
          mean = y_obs
        )
      } else {
        data.frame(
          variable = variables,
          mean = rnorm(length(variables)),
          median = rnorm(length(variables)),
          sd = rexp(length(variables)),
          rhat = rep(1.0, length(variables)),
          ess_bulk = rep(800, length(variables))
        )
      }
    },
    # reconstruct(): fit$draws("y_mis"); plot_trace()/plot_posterior_
    # densities(): fit$draws(variables = ..., format = "df")
    draws = function(variables = "y_mis", format = NULL, ...) {
      if (identical(format, "df")) {
        out <- as.data.frame(
          matrix(rnorm(n_draws * length(variables)), nrow = n_draws,
                 dimnames = list(NULL, variables))
        )
        out$.chain <- rep(1:2, each = n_draws / 2)
        out$.iteration <- rep(seq_len(n_draws / 2), times = 2)
        out$.draw <- seq_len(n_draws)
        out
      } else {
        # matrix draws of whole vector variables, e.g. "y_mis"
        keep <- sub("\\[.*", "", colnames(all_draws)) %in% variables
        posterior::as_draws_matrix(all_draws[, keep, drop = FALSE])
      }
    },
    num_chains = function() 2L,
    metadata = function() list(iter_sampling = n_draws / 2)
  )

  structure(
    list(
      fit     = mock_cmdstan,
      years   = years,
      idx_obs = idx_obs,
      idx_mis = idx_mis,
      data    = list(y_obs = y_obs),
      draws_mock = all_draws,
      call    = quote(fit_baygmst(mocked = TRUE))
    ),
    class = "baygmst_fit"
  )
}

test_that("reconstruct() builds a tidy, year-sorted data.frame", {
  fit <- make_mock_fit()
  recon <- reconstruct(fit)
  expect_s3_class(recon, "data.frame")
  expect_named(
    recon,
    c("year", "type", "T_obs", "T_mean",
      "T_lo_inner", "T_hi_inner", "T_lo_outer", "T_hi_outer")
  )
  expect_equal(nrow(recon), length(fit$years))
  expect_true(!is.unsorted(recon$year))
  expect_setequal(unique(recon$type), c("instrumental", "reconstruction"))
  expect_true(all(is.na(recon$T_obs[recon$type == "reconstruction"])))
  expect_false(anyNA(recon$T_mean))
})

test_that("reconstruction_ensemble() returns coherent thinned members", {
  fit <- make_mock_fit(n_obs = 20, n_mis = 30, n_draws = 40)
  ens <- reconstruction_ensemble(fit, n_members = 10)
  expect_s3_class(ens, "data.frame")
  expect_equal(nrow(ens), length(fit$years))
  expect_named(ens, c("year", "type", "T_obs", sprintf("m%04d", 1:10)))
  keep <- attr(ens, "draws")
  expect_equal(keep, round(seq(1, 40, length.out = 10)))
  # member j is draw keep[j] of y_mis and of y_ins_fitted (same draw)
  expect_equal(ens$m0003[fit$idx_mis],
               unname(fit$draws_mock[keep[3], paste0("y_mis[", 1:30, "]")]))
  expect_equal(ens$m0003[fit$idx_obs],
               unname(fit$draws_mock[keep[3],
                                     paste0("y_ins_fitted[", 1:20, "]")]))
  expect_equal(ens$T_obs[fit$idx_obs], fit$data$y_obs)
  expect_true(all(is.na(ens$T_obs[fit$idx_mis])))
  expect_false(anyNA(ens[, -(1:3)]))

  obs <- reconstruction_ensemble(fit, n_members = 5, instrumental = "observed")
  expect_equal(obs$m0005[fit$idx_obs], fit$data$y_obs)
  expect_equal(obs$m0001[fit$idx_mis], ens$m0001[fit$idx_mis])

  expect_equal(ncol(reconstruction_ensemble(fit, n_members = NULL)), 3 + 40)
  expect_error(reconstruction_ensemble(fit, n_members = 41), "exceeds")
  expect_error(reconstruction_ensemble(fit, n_members = 0), "positive")
})

test_that("summary() and print methods work on a mocked fit", {
  fit <- make_mock_fit()
  s <- summary(fit)
  expect_s3_class(s, "summary.baygmst_fit")
  expect_named(s, c("posterior", "performance", "convergence"))
  expect_named(
    s$performance,
    c("RMSE", "MAE", "Bias", "Correlation", "R2", "R2_detrended")
  )
  expect_output(print(s), "Posterior parameter summary")
  expect_output(print(fit), "baygmst_fit")

  # a converged mock (rhat = 1.0) must not print the convergence warning
  expect_false(any(grepl("not converged", capture.output(print(s)))))
})

test_that("print.summary.baygmst_fit warns on non-convergence", {
  s <- structure(
    list(
      posterior = data.frame(variable = "alpha1", mean = 0, rhat = 1.9,
                             ess_bulk = 3),
      performance = data.frame(RMSE = 0.1, MAE = 0.1, Bias = 0,
                               Correlation = 0.9, R2 = 0.8,
                               R2_detrended = 0.8),
      convergence = data.frame(max_rhat = 1.9, min_ess_bulk = 3)
    ),
    class = "summary.baygmst_fit"
  )
  expect_output(print(s), "not converged")
})

test_that("plot functions return ggplot objects", {
  fit <- make_mock_fit()
  expect_s3_class(plot_reconstruction(fit), "ggplot")
  expect_s3_class(plot_trace(fit), "ggplot")
  expect_s3_class(plot_posterior_densities(fit), "ggplot")
})

test_that("plot() on a fit builds the combined summary figure", {
  skip_if_not_installed("patchwork")
  fit <- make_mock_fit()
  p <- plot(fit)
  expect_s3_class(p, "patchwork")
  expect_match(p$patches$annotation$subtitle, "Instrumental period: \\(31, 50\\)")
  expect_null(plot(fit, subtitle = NULL)$patches$annotation$subtitle)
})

test_that("plot_cv() and print.baygmst_cv work on a mocked cv object", {
  set.seed(7)
  folds <- data.frame(
    year = 1951:1980,
    fold = rep(1:3, each = 10),
    T_true = rnorm(30),
    T_mean = rnorm(30),
    T_lo_inner = rnorm(30) - 0.5, T_hi_inner = rnorm(30) + 0.5,
    T_lo_outer = rnorm(30) - 1, T_hi_outer = rnorm(30) + 1
  )
  cv <- structure(
    list(folds = folds, r2 = c(fold1 = 0.8, fold2 = 0.7, fold3 = 0.75),
         mse = c(fold1 = 0.1, fold2 = 0.2, fold3 = 0.15),
         years = 1951:1980, call = quote(cv_baygmst(mocked = TRUE))),
    class = "baygmst_cv"
  )
  expect_s3_class(plot_cv(cv), "ggplot")
  expect_output(print(cv), "baygmst_cv")
})
