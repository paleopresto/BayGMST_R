# Trace plots for a fitted BayGMST model

Ported from the trace-plot block of the original `BayGMST_v1.0.R`
script.

## Usage

``` r
plot_trace(
  fit,
  parameters = c("alpha1", "betaG", "betaV", "betaS", "phi_R", "phi_T")
)
```

## Arguments

- fit:

  A `"baygmst_fit"` object, as returned by
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

- parameters:

  Character vector of Stan parameter names to plot. Default
  `c("alpha1", "betaG", "betaV", "betaS", "phi_R", "phi_T")`.

## Value

A `ggplot` object, faceted by parameter.
