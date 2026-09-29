test_that("transform_forcings applies the documented transforms", {
  out <- transform_forcings(
    G = c(280, 560),
    V = c(0, 1),
    S = c(1360, 1362)
  )
  expect_equal(out$G, 5.35 * log(c(280, 560) / 280))
  expect_equal(out$V, -abs(25) * (1 - exp(-c(0, 1))))
  expect_equal(out$S, c(1360, 1362) - mean(c(1360, 1362)))
})

test_that("transform_forcings honors custom coefficients", {
  out <- transform_forcings(G = 300, V = 0.5, S = 1361, co2_c0 = 300, co2_coef = 1, vol_coef = 10)
  expect_equal(out$G, 1 * log(300 / 300))
  expect_equal(out$V, -10 * (1 - exp(-0.5)))
})

test_that("transform_forcings validates input lengths", {
  expect_error(transform_forcings(G = c(1, 2), V = 1, S = c(1, 2)))
})
