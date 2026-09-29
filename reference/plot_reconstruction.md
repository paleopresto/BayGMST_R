# Plot a BayGMST reconstruction

Ported from the `p_ts` plotting block of the original `BayGMST_v1.0.R`
script. Unlike the original script, this returns a `ggplot` object
instead of calling `ggsave()` as a side effect – save it yourself with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
if you want a file.

## Usage

``` r
plot_reconstruction(fit, title = "GMST Reconstruction")
```

## Arguments

- fit:

  A `"baygmst_fit"` object, as returned by
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

- title:

  Plot title. Default `"GMST Reconstruction"`.

## Value

A `ggplot` object.
