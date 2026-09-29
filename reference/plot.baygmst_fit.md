# Plot a fitted BayGMST model as a single summary figure

Reproduces the combined figure of the original `BayGMST_v1.0.R` script
(and the PReSto manuscript): the reconstruction from
[`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md)
on top, and below it three posterior-density panels grouping the proxy
coefficient (`alpha1`), the forcing sensitivities (`betaG`, `betaV`,
`betaS`), and the autoregressive coefficients (`phi_R`, `phi_T`). The
individual pieces remain available through
[`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md)
and
[`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md).

## Usage

``` r
# S3 method for class 'baygmst_fit'
plot(
  x,
  title = "GMST Reconstruction using a Reduced Proxy",
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A `"baygmst_fit"` object, as returned by
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

- title:

  Figure title. Default `"GMST Reconstruction using a Reduced Proxy"`.

- subtitle:

  Figure subtitle. By default, states the instrumental period and the
  AR(1) model structure; pass `NULL` to omit it.

- ...:

  Unused; for compatibility with the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.

## Value

A 'patchwork' object (which is also a `ggplot`); save it with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
if a file is wanted.

## Details

Requires the 'patchwork' package.

## See also

[`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md),
[`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md),
[`plot_trace()`](https://paleopresto.github.io/BayGMST_R/reference/plot_trace.md)
