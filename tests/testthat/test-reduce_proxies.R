make_synthetic_proxies <- function(n = 100, seed = 1) {
  set.seed(seed)
  years <- seq_len(n)
  true_signal <- cumsum(rnorm(n, 0, 0.1))
  proxy_matrix <- cbind(
    p1 = true_signal + rnorm(n, 0, 0.05),
    p2 = true_signal + rnorm(n, 0, 0.05),
    p3 = true_signal + rnorm(n, 0, 0.05)
  )
  list(years = years, proxy_matrix = proxy_matrix, true_signal = true_signal)
}

test_that("reduce_proxies(method = 'PCR') runs on a single, fully-covered segment", {
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  out <- reduce_proxies(
    proxy_matrix = d$proxy_matrix,
    years        = d$years,
    temp_calib   = d$true_signal[calib_years],
    calib_years  = calib_years,
    method       = "PCR",
    chunk        = 100 # single segment, fully contains calib_years
  )
  expect_s3_class(out$composite, "data.frame")
  expect_named(out$composite, c("year", "RP1"))
  expect_equal(nrow(out$composite), 100)
  expect_false(anyNA(out$composite$RP1))
  # a well-specified PCR reconstruction should correlate strongly with the
  # (unobservable, in real life) true signal it was built to recover
  expect_gt(stats::cor(out$composite$RP1, d$true_signal), 0.8)
})

test_that("reduce_proxies rejects an unknown method", {
  d <- make_synthetic_proxies()
  expect_error(
    reduce_proxies(d$proxy_matrix, d$years, d$true_signal[51:100], 51:100, method = "bogus")
  )
})

test_that("reduce_proxies errors clearly when a segment doesn't fully contain calib_years", {
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  expect_error(
    reduce_proxies(
      proxy_matrix = d$proxy_matrix,
      years        = d$years,
      temp_calib   = d$true_signal[calib_years],
      calib_years  = calib_years,
      method       = "PCR",
      chunk        = 10 # too short: later segments won't span all of calib_years
    ),
    "does not fully contain the calibration"
  )
})

test_that("reduce_proxies validates temp_calib/calib_years length match", {
  d <- make_synthetic_proxies()
  expect_error(
    reduce_proxies(d$proxy_matrix, d$years, temp_calib = 1:5, calib_years = 51:100)
  )
})

test_that("reduce_proxies returns a classed object that fit_baygmst accepts", {
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  out <- reduce_proxies(
    d$proxy_matrix, d$years, d$true_signal[calib_years], calib_years,
    method = "PCR", chunk = 100
  )
  expect_s3_class(out, "baygmst_rp")
  p <- as_baygmst_proxy(out)
  expect_s3_class(p, "baygmst_proxy")
  expect_equal(p$proxy, out$composite$RP1)
})

test_that("reduce_proxies(method = 'LASSO') runs", {
  skip_if_not_installed("glmnet")
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  out <- reduce_proxies(
    d$proxy_matrix, d$years, d$true_signal[calib_years], calib_years,
    method = "LASSO", chunk = 100
  )
  expect_false(anyNA(out$composite$RP1))
  expect_gt(stats::cor(out$composite$RP1, d$true_signal), 0.7)
})

test_that("reduce_proxies(method = 'sPCR') runs", {
  skip_if_not_installed("superpc")
  # superpc needs more features than make_synthetic_proxies()'s 3
  set.seed(3)
  n <- 100
  true_signal <- cumsum(rnorm(n, 0, 0.1))
  proxy_matrix <- sapply(1:12, function(i) true_signal + rnorm(n, 0, 0.05))
  calib_years <- 51:100
  out <- suppressWarnings(reduce_proxies(
    proxy_matrix, seq_len(n), true_signal[calib_years], calib_years,
    method = "sPCR", chunk = 100, n_components = 1
  ))
  expect_false(anyNA(out$composite$RP1))
})

test_that("reduce_proxies(method = 'SPLS') runs", {
  skip_if_not_installed("spls")
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  # cv.spls prints progress; keep the grid tiny for speed
  out <- suppressWarnings(reduce_proxies(
    d$proxy_matrix, d$years, d$true_signal[calib_years], calib_years,
    method = "SPLS", chunk = 100, spls_eta = c(0.3, 0.6)
  ))
  expect_false(anyNA(out$composite$RP1))
})

test_that("reduce_proxies(method = 'SIR') warns that it is experimental", {
  skip_if_not_installed("dr")
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  expect_warning(
    reduce_proxies(
      d$proxy_matrix, d$years, d$true_signal[calib_years], calib_years,
      method = "SIR", chunk = 100
    ),
    "experimental"
  )
})

test_that("reduce_proxies screens missingness before zero-filling gaps", {
  set.seed(4)
  n <- 300
  years <- seq_len(n)
  true_signal <- cumsum(rnorm(n, 0, 0.1))
  noisy <- function() true_signal + rnorm(n, 0, 0.05)
  decadal <- noisy()
  decadal[years %% 10 != 0] <- NA # 90% missing after its first value
  ends_early <- noisy()
  ends_early[years > 150] <- NA # no data in the calibration window
  late_sparse <- noisy()
  late_sparse[1:120] <- NA # starts after row 101; only seen by segment 2
  late_sparse[seq(130, n, by = 4)] <- NA # 25% internal gaps
  proxy_matrix <- cbind(
    p1 = noisy(), p2 = noisy(), p3 = noisy(),
    decadal = decadal, ends_early = ends_early, late_sparse = late_sparse
  )
  calib_years <- 201:300
  out <- reduce_proxies(
    proxy_matrix, years, true_signal[calib_years], calib_years,
    method = "PCR", chunk = 100
  )
  expect_named(out$proxies, c("RP1", "RP2", "RP3"))
  for (k in seq_along(out$proxies)) {
    expect_setequal(out$proxies[[k]], c("p1", "p2", "p3"))
  }
  expect_false(anyNA(out$composite$RP1))
})

test_that("reduce_proxies records proxy indices when columns are unnamed", {
  d <- make_synthetic_proxies()
  calib_years <- 51:100
  out <- reduce_proxies(
    unname(d$proxy_matrix), d$years, d$true_signal[calib_years], calib_years,
    method = "PCR", chunk = 100
  )
  expect_equal(out$proxies, list(RP1 = c("1", "2", "3")))
})
