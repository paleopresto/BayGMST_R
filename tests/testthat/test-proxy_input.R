test_that("as_baygmst_proxy.numeric requires years and matching lengths", {
  expect_error(as_baygmst_proxy(rnorm(10)), "`years` is required")
  expect_error(
    as_baygmst_proxy(rnorm(10), years = 1:5),
    "must match length"
  )
  p <- as_baygmst_proxy(rnorm(10), years = 2001:2010)
  expect_s3_class(p, "baygmst_proxy")
  expect_named(p, c("year", "proxy"))
  expect_equal(p$year, 2001:2010)
  expect_identical(attr(p, "n_ensemble"), 1L)
})

test_that("as_baygmst_proxy.data.frame detects the year column", {
  df <- data.frame(Year = 1901:1910, RP1 = rnorm(10))
  p <- as_baygmst_proxy(df)
  expect_equal(p$year, 1901:1910)
  expect_equal(p$proxy, df$RP1)

  expect_error(
    as_baygmst_proxy(data.frame(a = 1:5, b = rnorm(5))),
    "must have a `year`"
  )
  expect_error(
    as_baygmst_proxy(data.frame(year = 1:5, a = rnorm(5), b = rnorm(5))),
    "exactly one value column"
  )
})

test_that("as_baygmst_proxy.matrix collapses ensembles", {
  set.seed(1)
  m <- matrix(rnorm(50), nrow = 10, ncol = 5)
  expect_error(as_baygmst_proxy(m), "`years` is required")
  expect_error(as_baygmst_proxy(m, years = 1:3), "must match nrow")

  p <- as_baygmst_proxy(m, years = 1901:1910)
  expect_identical(attr(p, "n_ensemble"), 5L)
  expect_equal(p$proxy, apply(m, 1, stats::median))

  p_mean <- as_baygmst_proxy(m, years = 1901:1910, collapse = "mean")
  expect_equal(p_mean$proxy, rowMeans(m))
})

test_that("as_baygmst_proxy.list handles composite (paleoComposite-shaped) input", {
  set.seed(2)
  # ages in yr BP, descending, with an all-NA bin, tibble-like composite
  comp <- list(
    ages = seq(950, 0, by = -50),
    composite = as.data.frame(matrix(rnorm(20 * 8), nrow = 20, ncol = 8))
  )
  comp$composite[3, ] <- NA

  expect_message(p <- as_baygmst_proxy(comp), "dropped")
  expect_s3_class(p, "baygmst_proxy")
  expect_equal(nrow(p), 19)
  expect_true(all(diff(p$year) > 0))
  expect_equal(range(p$year), c(1000, 1950)) # 1950 - c(950, 0)
  expect_identical(attr(p, "n_ensemble"), 8L)

  # ages already in CE
  comp_ce <- list(ages = 1001:1010, composite = matrix(rnorm(10), ncol = 1))
  p_ce <- as_baygmst_proxy(comp_ce, age_units = "CE")
  expect_equal(p_ce$year, 1001:1010)

  # shape mismatch and non-composite lists error informatively
  expect_error(
    as_baygmst_proxy(list(ages = 1:5, composite = matrix(rnorm(12), 6, 2))),
    "must match length"
  )
  expect_error(as_baygmst_proxy(list(a = 1)), "composite object")
})

test_that("as_baygmst_proxy accepts reduce_proxies()-style output", {
  rp_like <- list(
    segments = NULL,
    composite = data.frame(year = 1:50, RP1 = rnorm(50)),
    method = "PCR"
  )
  # unclassed (e.g. saved from an older version)
  p <- as_baygmst_proxy(rp_like)
  expect_equal(p$proxy, rp_like$composite$RP1)

  # classed, as reduce_proxies() now returns
  class(rp_like) <- "baygmst_rp"
  p2 <- as_baygmst_proxy(rp_like)
  expect_equal(p2$proxy, rp_like$composite$RP1)
  expect_match(attr(p2, "source"), "PCR")
})

test_that("as_baygmst_proxy rejects duplicated years and all-missing input", {
  expect_error(
    as_baygmst_proxy(rnorm(4), years = c(1, 1, 2, 3)),
    "duplicated years"
  )
  expect_error(
    suppressMessages(
      as_baygmst_proxy(rep(NA_real_, 3), years = 1:3)
    ),
    "No finite proxy values"
  )
})

test_that("print.baygmst_proxy and print.baygmst_rp print summaries", {
  p <- as_baygmst_proxy(rnorm(10), years = 2001:2010)
  expect_output(print(p), "baygmst_proxy")

  rp <- structure(
    list(
      segments = data.frame(year = 1:10, RP1 = rnorm(10)),
      composite = data.frame(year = 1:10, RP1 = rnorm(10)),
      method = "PCR"
    ),
    class = "baygmst_rp"
  )
  expect_output(print(rp), "baygmst_rp")
  expect_output(print(rp), "PCR")
})

test_that("resolve_proxy passes bare vectors through positionally", {
  expect_identical(resolve_proxy(1:5, years = 2001:2005), 1:5)
})

test_that("resolve_proxy aligns year-aware objects and errors on gaps", {
  df <- data.frame(year = 1901:1950, RP1 = seq_len(50))
  # aligned subset, in order
  expect_equal(resolve_proxy(df, years = 1910:1919), 10:19)
  # requesting uncovered years errors informatively
  expect_error(
    resolve_proxy(df, years = 1899:1905),
    "does not cover"
  )
})

test_that("fit_baygmst validates inputs informatively (no CmdStan needed)", {
  years <- 1:50
  ok <- rnorm(50)
  expect_error(
    fit_baygmst(ok[1:10], ok, ok, ok, ok, years),
    "one value per year"
  )
  expect_error(
    fit_baygmst(c(NA, ok[-1]), ok, ok, ok, ok, years),
    "Missing values are not allowed"
  )
  expect_error(
    fit_baygmst(ok, rep(NA_real_, 50), ok, ok, ok, years),
    "no observed"
  )
  expect_error(
    fit_baygmst(ok, ok, ok, ok, ok, years, iter_sampling = 10),
    "at least 1000"
  )
})
