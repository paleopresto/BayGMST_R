# Plot k-fold cross-validation reconstructions

Ported from the `p_ts_cv` plotting block of the original `cv_v0.1.R`
script.

## Usage

``` r
plot_cv(cv)
```

## Arguments

- cv:

  A `"baygmst_cv"` object, as returned by
  [`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md).

## Value

A `ggplot` object showing each fold's held-out reconstruction against
its true instrumental value, faceted by fold.
