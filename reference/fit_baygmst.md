# Fit the BayGMST Bayesian hierarchical AR(1) reconstruction model

Prepares the data list expected by the `baygmst` 'Stan' model and
samples from it with cmdstanr, via instantiate. Ported from the
data-preparation and `mod$sample()` block of the original
`BayGMST_v1.0.R` script (see
`system.file("legacy-scripts", package = "BayGMST")`), restructured to
take its inputs as explicit arguments (rather than reading `config.yml`
and hardcoded relative paths) and to return an object rather than
leaving variables in the global environment.

## Usage

``` r
fit_baygmst(
  proxy,
  instrumental_T,
  forcing_G,
  forcing_V,
  forcing_S,
  years,
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

- chains, parallel_chains, iter_warmup, iter_sampling, seed:

  Passed to `cmdstanr`'s `$sample()` method. `iter_sampling` must be at
  least 1000 (the original script enforced the same minimum).

- ...:

  Additional arguments passed on to `$sample()`.

## Value

An object of class `"baygmst_fit"`, a list with elements:

- fit:

  The underlying `CmdStanMCMC` object (see
  [`cmdstanr::CmdStanMCMC`](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html)).

- years:

  The `years` argument, echoed back.

- idx_obs, idx_mis:

  Integer indices (into `years`) of the instrumental and
  pre-instrumental years, respectively.

- data:

  The data list passed to Stan.

- cmdstan_version:

  The CmdStan version used, for provenance.

- call:

  The matched call.

Use
[`reconstruct()`](https://paleopresto.github.io/BayGMST_R/reference/reconstruct.md)
to extract a tidy reconstruction `data.frame`,
[plot()](https://paleopresto.github.io/BayGMST_R/reference/plot.baygmst_fit.md)
or
[`plot_reconstruction()`](https://paleopresto.github.io/BayGMST_R/reference/plot_reconstruction.md)
/
[`plot_trace()`](https://paleopresto.github.io/BayGMST_R/reference/plot_trace.md)
/
[`plot_posterior_densities()`](https://paleopresto.github.io/BayGMST_R/reference/plot_posterior_densities.md)
to visualize it, and [`summary()`](https://rdrr.io/r/base/summary.html)
for posterior parameter summaries.

## CmdStan required

This function requires a working CmdStan installation (via cmdstanr /
instantiate). Check availability first with
[`instantiate::stan_cmdstan_exists()`](https://wlandau.github.io/instantiate/reference/stan_cmdstan_exists.html).

## Examples

``` r
# \donttest{
if (instantiate::stan_cmdstan_exists()) {
  set.seed(1)
  years <- 1:200
  forcing <- transform_forcings(
    G = 280 + cumsum(rgamma(200, 0.05, 1)),
    V = abs(rnorm(200, 0, 0.05)),
    S = 1361 + rnorm(200, 0, 0.3)
  )
  true_T <- cumsum(rnorm(200, 0, 0.05))
  instrumental_T <- c(rep(NA_real_, 150), true_T[151:200])
  proxy <- true_T + rnorm(200, 0, 0.1)

  fit <- fit_baygmst(
    proxy = proxy,
    instrumental_T = instrumental_T,
    forcing_G = forcing$G,
    forcing_V = forcing$V,
    forcing_S = forcing$S,
    years = years,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 1000
  )
  summary(fit)
}
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 1 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 1 Iteration:  201 / 1200 [ 16%]  (Sampling) 
#> Chain 1 Iteration:  300 / 1200 [ 25%]  (Sampling) 
#> Chain 1 Iteration:  400 / 1200 [ 33%]  (Sampling) 
#> Chain 1 Iteration:  500 / 1200 [ 41%]  (Sampling) 
#> Chain 1 Iteration:  600 / 1200 [ 50%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1200 [ 58%]  (Sampling) 
#> Chain 1 Iteration:  800 / 1200 [ 66%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1200 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1200 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 2.1 seconds.
#> Posterior parameter summary:
#>    variable         mean       median          sd         mad          q5
#> 1    alpha0 -0.016950119 -0.013973503 0.050045327 0.051853721 -0.10217389
#> 2    alpha1  1.075244685  1.075303650 0.126091868 0.130335366  0.87270268
#> 3     phi_R -0.067253481 -0.071341387 0.076960773 0.074319327 -0.18873288
#> 4     phi_T  0.908247713  0.911418210 0.034134300 0.032896752  0.84624421
#> 5     beta0 -0.006005973 -0.005524942 0.013614786 0.013834332 -0.02972632
#> 6     betaG  0.310139779  0.305632435 0.170446204 0.157340651  0.04823441
#> 7     betaS -0.021733753 -0.021592369 0.022038289 0.022350523 -0.05818163
#> 8     betaV -0.006413797 -0.006511145 0.008950821 0.008947334 -0.02162157
#> 9   sigma_y  0.059351832  0.058771752 0.005904955 0.005725961  0.05048344
#> 10  sigma_z  0.101165445  0.100895755 0.007381874 0.007371398  0.08996498
#>           q95      rhat   ess_bulk ess_tail
#> 1  0.06298867 1.0094340   87.82032 226.9748
#> 2  1.28255167 1.0100986  131.06270 331.4295
#> 3  0.06174513 0.9991444  532.52555 743.7651
#> 4  0.96143386 1.0098840  313.46586 393.9780
#> 5  0.01537595 0.9997764  395.33910 468.7778
#> 6  0.59477989 1.0028570  503.71110 409.5329
#> 7  0.01334192 1.0024985 1204.23694 753.8090
#> 8  0.00838420 0.9991603 1077.54430 734.3666
#> 9  0.06947748 0.9990908  205.36513 305.5818
#> 10 0.11394000 1.0004887  514.20078 642.2759
#> 
#> Instrumental-period fit performance:
#>        RMSE        MAE        Bias Correlation        R2 R2_detrended
#>  0.06164432 0.05021732 -0.01488405   0.7210831 0.5199609    0.4878794
# }
```
