# This file is part of the standard testthat setup; see
# https://testthat.r-lib.org/articles/special-files.html
#
# NOTE (BayGMST packaging): instantiate documents that pkgload::load_all()
# (which devtools::test() uses) is not compatible with packages built around
# instantiate::stan_package_model(). Run tests via `R CMD check` (which
# installs the package first) rather than `devtools::test()`. See
# CRAN-READINESS.md.

library(testthat)
library(BayGMST)

test_check("BayGMST")
