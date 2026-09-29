# Extract a tidy reconstruction from a fitted BayGMST model

Builds a tidy `data.frame` combining the observed instrumental
temperatures with the posterior reconstruction of the pre-instrumental
years, equivalent to the `gmst_reconstruction_data.csv` produced as a
side effect by the original `BayGMST_v1.0.R` script – but returned as an
object rather than written to disk. Save it yourself
(`write.csv(reconstruct(fit), "path/to/file.csv")`) if you want a file.

## Usage

``` r
reconstruct(fit, probs_inner = c(0.16, 0.84), probs_outer = c(0.025, 0.975))
```

## Arguments

- fit:

  A `"baygmst_fit"` object, as returned by
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md).

- probs_inner, probs_outer:

  Numeric vectors of length 2, the lower/ upper quantile levels for the
  "inner" and "outer" credible bands. Defaults `c(0.16, 0.84)` (68%) and
  `c(0.025, 0.975)` (95%), matching the original script.

## Value

A `data.frame` with one row per year in `fit$years`, and columns:

- year:

  Calendar year.

- type:

  `"instrumental"` or `"reconstruction"`.

- T_obs:

  Observed instrumental temperature (`NA` for reconstructed years).

- T_mean:

  Posterior mean (instrumental years: posterior predictive mean of the
  fitted value; reconstructed years: posterior mean of the latent
  temperature).

- T_lo_inner, T_hi_inner:

  Inner credible band (`probs_inner`).

- T_lo_outer, T_hi_outer:

  Outer credible band (`probs_outer`).

## See also

[`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md)
