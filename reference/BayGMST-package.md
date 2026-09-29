# BayGMST: Bayesian Reconstruction of Global Mean Surface Temperature

Fits a Bayesian hierarchical AR(1) state-space model relating global
mean surface temperature (GMST) to a reduced paleoclimate proxy network,
instrumental temperature observations, and radiative forcings
(greenhouse gases, volcanic aerosols, and solar irradiance). See
[`vignette("baygmst-intro", package = "BayGMST")`](https://paleopresto.github.io/BayGMST_R/articles/baygmst-intro.md)
for a worked end-to-end example, and
[`vignette("baygmst-model", package = "BayGMST")`](https://paleopresto.github.io/BayGMST_R/articles/baygmst-model.md)
for the full statistical model description, which builds on the
reduced-proxy Bayesian framework of Barboza et al. (2014, 2019) and Wang
(2020).

## Model fitting requires CmdStan

[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
and
[`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md)
compile and sample a 'Stan' model via the instantiate and cmdstanr
packages, which in turn require a working CmdStan installation on the
user's machine. Use
[`instantiate::stan_cmdstan_exists()`](https://wlandau.github.io/instantiate/reference/stan_cmdstan_exists.html)
to check whether one is available, and see
[`vignette("cmdstanr", package = "cmdstanr")`](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
for setup instructions.
[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
and
[`transform_forcings()`](https://paleopresto.github.io/BayGMST_R/reference/transform_forcings.md)
do not require CmdStan.

## Provenance

This package was restructured from a collection of top-level analysis
scripts; see `NEWS.md` for exactly what changed (fixed bugs, dropped
dead code, and modeling decisions reviewed by the authors) and
`system.file("legacy-scripts", package = "BayGMST")` for the original,
unmodified scripts.

## See also

Useful links:

- <https://github.com/paleopresto/BayGMST_R>

- Report bugs at <https://github.com/paleopresto/BayGMST_R/issues>

## Author

**Maintainer**: Nick McKay <nick@nau.edu>

Authors:

- Tyler Bagwell <teb6@rice.edu>

- Julien Emile-Geay <julieneg@usc.edu>

- Frederi Viens
