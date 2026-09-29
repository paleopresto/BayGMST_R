# K-fold cross-validate the BayGMST model over the instrumental period

Partitions the instrumental (non-`NA`) years of `instrumental_T` into
`nfold` contiguous, roughly equal blocks; for each block in turn, refits
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)'s
underlying model with that block's temperatures hidden (treated as
missing/to-be-reconstructed), and reports how well the model recovers
them. Ported from the fold loop in the original `cv_v0.1.R` script (see
`system.file("legacy-scripts", package = "BayGMST")`).

## Usage

``` r
cv_baygmst(
  proxy,
  instrumental_T,
  forcing_G,
  forcing_V,
  forcing_S,
  years,
  nfold = 3,
  chains = 4,
  parallel_chains = 1,
  iter_warmup = 500,
  iter_sampling = 1500,
  seed = NULL,
  ...
)
```

## Arguments

- proxy:

  The representative/reduced proxy series. Either a bare numeric vector
  with one value per year in `years` (aligned by position, no missing
  values), or any object accepted by
  [`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md)
  – the output of
  [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md),
  a `data.frame` with year and value columns, or a composite object from
  another package (e.g. a 'compositeR' `paleoComposite`) – in which case
  it is aligned to `years` by calendar year, with an informative error
  if any requested year is not covered.

- instrumental_T:

  Numeric vector, instrumental temperature, one value per year in
  `years`. `NA` for years without instrumental coverage (the
  pre-instrumental years the model reconstructs).

- forcing_G, forcing_V, forcing_S:

  Numeric vectors of *already transformed* greenhouse-gas, volcanic, and
  solar forcing – see
  [`transform_forcings()`](https://paleopresto.github.io/BayGMST_R/reference/transform_forcings.md)
  – one value per year in `years`. No missing values are allowed.

- years:

  Integer vector of calendar years, defining the row order for all of
  the above.

- nfold:

  Number of folds. Default `3`, matching the original script.

- chains, parallel_chains, iter_warmup, iter_sampling, seed:

  Passed to `cmdstanr`'s `$sample()` method. `iter_sampling` must be at
  least 1000 (the original script enforced the same minimum).

- ...:

  Additional arguments passed on to `$sample()`.

## Value

An object of class `"baygmst_cv"`, a list with elements:

- folds:

  A single combined data.frame (one row per held-out year, across all
  folds) with columns `year`, `fold`, `T_true`, `T_mean`, `T_lo_inner`,
  `T_hi_inner`, `T_lo_outer`, `T_hi_outer`.

- r2:

  Named numeric vector, the posterior mean of the Stan-computed `r2_cv`
  generated quantity for each fold (fully Bayesian; this differs
  slightly from the plug-in point-estimate R-squared computed in the
  original R script – see Details in the package vignette).

- mse:

  As `r2`, for `mse_cv`.

- years, call:

  Echoed back.

## Details

Two blocks present in the original script were **not** ported: a
`Forcings.projections`/RCP-scenario block that was built but never
referenced again (dead code, and it hardcoded a path to a specific
contributor's machine), and a large commented-out "OLD" duplicate of the
fold loop at the end of the file. See `NEWS.md`.

## CmdStan required

See
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

## See also

[`plot_cv()`](https://paleopresto.github.io/BayGMST_R/reference/plot_cv.md)
