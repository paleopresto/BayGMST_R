# Summarize a fitted BayGMST model

Reports posterior summaries for the structural parameters (matching the
variable set used in the original `BayGMST_v1.0.R` script) together with
simple instrumental-period fit performance statistics (RMSE, MAE, bias,
correlation, and R-squared, both raw and linearly detrended).

## Usage

``` r
# S3 method for class 'baygmst_fit'
summary(object, ...)
```

## Arguments

- object:

  A `"baygmst_fit"` object, as returned by
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

- ...:

  Unused; present for S3 consistency.

## Value

An object of class `"summary.baygmst_fit"`, a list with elements
`posterior` (a data.frame of parameter posterior summaries),
`performance` (a one-row data.frame of instrumental-period fit
statistics), and `convergence` (a one-row data.frame with `max_rhat` and
`min_ess_bulk` across the structural parameters; the print method warns
when `max_rhat` exceeds 1.05).
