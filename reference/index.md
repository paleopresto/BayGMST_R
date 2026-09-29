# Package index

## Prepare proxy inputs

Reduce a multi-proxy network to a representative series, or bring a
reduced/composited series from elsewhere (including compositeR-style
composite ensembles).

- [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
  : Reduce a multi-proxy network to a single representative proxy series
- [`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md)
  : Coerce a reduced-proxy series into BayGMST's input format

## Transform forcings

- [`transform_forcings()`](https://paleopresto.github.io/BayGMST_R/reference/transform_forcings.md)
  :

  Transform raw radiative forcings for use with
  [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)

## Fit and cross-validate the model

Requires a working CmdStan installation.

- [`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
  : Fit the BayGMST Bayesian hierarchical AR(1) reconstruction model
- [`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md)
  : K-fold cross-validate the BayGMST model over the instrumental period

## Extract and visualize results

- [`reconstruct()`](https://paleopresto.github.io/BayGMST_R/reference/reconstruct.md)
  : Extract a tidy reconstruction from a fitted BayGMST model
- [`summary(`*`<baygmst_fit>`*`)`](https://paleopresto.github.io/BayGMST_R/reference/summary.baygmst_fit.md)
  : Summarize a fitted BayGMST model
- [`print(`*`<baygmst_fit>`*`)`](https://paleopresto.github.io/BayGMST_R/reference/print.baygmst_fit.md)
  : Print a fitted BayGMST model
- [`plot(`*`<baygmst_fit>`*`)`](https://paleopresto.github.io/BayGMST_R/reference/plot.baygmst_fit.md)
  : Plot a fitted BayGMST model as a single summary figure
- [`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md)
  : Plot a BayGMST reconstruction
- [`plot_trace()`](https://paleopresto.github.io/BayGMST_R/reference/plot_trace.md)
  : Trace plots for a fitted BayGMST model
- [`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md)
  : Posterior density plots for a fitted BayGMST model
- [`plot_cv()`](https://paleopresto.github.io/BayGMST_R/reference/plot_cv.md)
  : Plot k-fold cross-validation reconstructions

## Package overview

- [`BayGMST`](https://paleopresto.github.io/BayGMST_R/reference/BayGMST-package.md)
  [`BayGMST-package`](https://paleopresto.github.io/BayGMST_R/reference/BayGMST-package.md)
  : BayGMST: Bayesian Reconstruction of Global Mean Surface Temperature
