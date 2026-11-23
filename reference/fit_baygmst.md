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
#> 1    alpha0 -0.013662596 -0.013184454 0.039141454 0.040973870 -0.07799226
#> 2    alpha1  1.072471662  1.075100200 0.110571398 0.110049691  0.88886279
#> 3     phi_R -0.069508411 -0.071189997 0.080238314 0.078607244 -0.19965416
#> 4     phi_T  0.908668982  0.909504425 0.034599827 0.035535372  0.85206721
#> 5     beta0 -0.006073826 -0.005847249 0.013036583 0.012968865 -0.02847861
#> 6     betaG  0.305529045  0.310311970 0.175751155 0.188486934  0.02328199
#> 7     betaS -0.022286779 -0.021661999 0.022609129 0.023069567 -0.06138355
#> 8     betaV -0.006373070 -0.006214704 0.008777555 0.008640606 -0.02030879
#> 9   sigma_y  0.059026317  0.058858932 0.005693137 0.005509446  0.04988436
#> 10  sigma_z  0.100791573  0.100549830 0.007499653 0.007700372  0.08885259
#>            q95     rhat   ess_bulk ess_tail
#> 1  0.047943400 1.014372   94.78349 247.2896
#> 2  1.249779740 1.000515  182.67420 400.3008
#> 3  0.062509676 1.003780  433.94397 488.5911
#> 4  0.961824242 1.004517  309.41784 579.6001
#> 5  0.015554666 1.005957  596.08729 709.5762
#> 6  0.600663912 1.005744  677.40961 669.5255
#> 7  0.014419137 1.002659  817.30643 715.0619
#> 8  0.008365304 1.002921 1075.91131 858.2240
#> 9  0.068976836 1.038388   45.11408 512.5929
#> 10 0.112859806 1.010026  365.11416 425.0081
#> 
#> Instrumental-period fit performance:
#>        RMSE        MAE        Bias Correlation        R2 R2_detrended
#>  0.06241168 0.05122751 -0.01745454    0.720217 0.5187126    0.4871168
# }
```
