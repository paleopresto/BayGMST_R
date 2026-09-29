# Reduce a multi-proxy network to a single representative proxy series

Collapses a matrix of paleoclimate proxy records into one (or a small
set of) "representative proxy" series suitable for
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md),
following the segmented-composite approach of Barboza et al. (2019) as
ported from `utils/PAGES2k_reducedProxy_UNSC.R` (see
`system.file("legacy-scripts", package = "BayGMST")` for the original).
The proxy network is broken into overlapping time segments of `chunk`
years; within each segment, only proxies that already have data are
used, and are regressed against instrumental temperature during
`calib_years` using the chosen `method`. The final composite series
averages across all segments that cover a given year.

## Usage

``` r
reduce_proxies(
  proxy_matrix,
  years,
  temp_calib,
  calib_years,
  method = c("PCR", "LASSO", "sPCR", "SPLS", "SIR"),
  chunk = 250,
  target_adj_r2 = 0.7,
  max_na_frac = 0.05,
  n_components = 3,
  spls_eta = seq(0.1, 0.9, 0.1)
)
```

## Arguments

- proxy_matrix:

  Numeric matrix or data.frame of proxy values, one row per year (in the
  same order as `years`), one column per proxy record. Do not include a
  year column here – pass it separately via `years`.

- years:

  Integer vector of calendar years, one per row of `proxy_matrix`.

- temp_calib:

  Numeric vector of instrumental temperatures during the calibration
  period, one value per year in `calib_years`.

- calib_years:

  Integer vector of calendar years used for calibration (must be a
  subset of `years`, and `length(calib_years) == length(temp_calib)`).

- method:

  One of `"PCR"` (principal component regression, the default),
  `"LASSO"`, `"sPCR"` (supervised PCR), `"SPLS"` (sparse partial least
  squares), or `"SIR"`. **`"SIR"` is experimental** – see Details.

- chunk:

  Segment length, in years. Default `250`. Every segment (each of which
  extends from its own start year through to the *end* of the record,
  not just `chunk` years) must fully contain the calibration window,
  i.e. `chunk` must be small enough, relative to `length(years)`, that
  even the last, shortest segment still spans all of `calib_years` –
  otherwise `reduce_proxies()` errors rather than silently comparing
  mismatched vectors (this mirrors an explicit check present in the
  original script).

- target_adj_r2:

  For `method = "PCR"`: the calibration-period adjusted R-squared at
  which to stop adding principal components (the smallest number of PCs
  reaching this threshold is used; if none reach it, the number
  maximizing adjusted R-squared is used). Default `0.70`.

- max_na_frac:

  Proxies missing more than this fraction of observations within a
  segment (evaluated both over the full segment and over its calibration
  window) are dropped from that segment. Default `0.05`.

- n_components:

  For `method = "sPCR"` only: number of supervised principal components
  to use for prediction. A single value (recycled across segments) or a
  vector with one value per segment. Default `3`.

- spls_eta:

  For `method = "SPLS"` only: grid of `eta` sparsity values passed to
  [`spls::cv.spls()`](https://rdrr.io/pkg/spls/man/cv.spls.html).
  Default `seq(0.1, 0.9, 0.1)`.

## Value

An object of class `"baygmst_rp"`, a list with components:

- segments:

  A data.frame with columns `year`, `RP1`, ..., `RPns` (one column per
  time segment).

- composite:

  A data.frame with columns `year` and `RP1`, the across-segment average
  – the single representative proxy series to pass to
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
  (which accepts the whole object directly, via
  [`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md)).

- method:

  The method used, echoed back.

## Method support

`"PCR"`, `"LASSO"`, and `"SPLS"` are fully adaptive (all tuning
parameters are chosen by calibration-period cross-validation within each
segment) and are considered supported. `"sPCR"` is supported but
requires the caller to choose `n_components` (it is not cross-validated
in the underlying superpc workflow, matching the original script).

`"SIR"` is a simplified, **unvalidated** reimplementation using a single
sliced-inverse-regression direction from dr, with no per-segment
variable pre-selection. It intentionally does *not* reproduce the
original script's SIR branch, which called an undocumented `edrSelec()`
step from the non-CRAN package 'edrGraphicalTools' that the original
author flagged as never fully working ("NEED TO FIGURE OUT HOW TO
INSTALL edrGraphicalTools!!"). `config.yml` in the original repository
likewise documented SIR as "still under construction." Treat
`method = "SIR"` results with corresponding caution.

## Required packages

`"LASSO"` requires glmnet; `"sPCR"` requires superpc; `"SPLS"` requires
spls; `"SIR"` requires dr. These are `Suggests`, not hard dependencies
of BayGMST, and are only required for the method you actually call.

## Reproducibility

`"LASSO"`, `"sPCR"`, and `"SPLS"` use randomized cross-validation
internally. Call [`set.seed()`](https://rdrr.io/r/base/Random.html)
before calling this function if you need reproducible output;
`reduce_proxies()` itself never calls
[`set.seed()`](https://rdrr.io/r/base/Random.html), so it never perturbs
the caller's random number stream as a side effect.

## Examples

``` r
# \donttest{
set.seed(1)
years <- 1:500
proxy_matrix <- cbind(
  p1 = c(rep(NA_real_, 50), cumsum(rnorm(450))),
  p2 = cumsum(rnorm(500))
)
calib_years <- 401:500
temp_calib <- cumsum(rnorm(100)) / 10
rp <- reduce_proxies(
  proxy_matrix, years, temp_calib, calib_years,
  method = "PCR", chunk = 250
)
head(rp$composite)
#>   year       RP1
#> 1    1 -2.177415
#> 2    2 -2.135523
#> 3    3 -2.083341
#> 4    4 -1.975523
#> 5    5 -1.874755
#> 6    6 -1.878659
# }
```
