# Introduction to BayGMST

This vignette walks through the full BayGMST workflow end to end, using
small example datasets bundled with the package (trimmed subsets of the
data used in Bagwell et al., in prep.; see
`system.file("extdata", package = "BayGMST")`). It replaces the
`config.yml`-driven scripts the package was originally built from – see
`NEWS.md` and `system.file("legacy-scripts", package = "BayGMST")`.

``` r

library(BayGMST)
```

## 1. Load example data

``` r

proxy_matrix <- read.csv(
  system.file("extdata", "example_proxy_matrix.csv", package = "BayGMST")
)
forcing_raw <- read.csv(
  system.file("extdata", "example_forcing.csv", package = "BayGMST")
)
instrumental_raw <- read.csv(
  system.file("extdata", "example_instrumental_temperature.csv", package = "BayGMST")
)
# read.csv() mangles "Time"/"Anomaly (deg C)" headers via make.names(); give
# them stable names instead of relying on the mangled spelling.
colnames(instrumental_raw)[1:2] <- c("Time", "Anomaly")

years <- forcing_raw$year
head(proxy_matrix)
#>   year SAm_030  SAm_026 SAm_025 Ocn_095 NAm_186 NAm_104
#> 1 1700  9.2943 -20.9760   0.463 -4.5324   0.652   1.021
#> 2 1701  9.2947 -18.8329   0.474 -4.6960   0.577   0.941
#> 3 1702  9.2951 -18.7600   0.660 -4.6556   0.831   0.982
#> 4 1703  9.2955 -19.9733   0.980 -4.3724   1.045   0.994
#> 5 1704  9.2959 -19.9486   0.735 -4.4295   0.876   0.973
#> 6 1705  9.2963 -21.1983  -0.357 -4.6110   0.672   0.990
```

## 2. Reduce the proxy network to a single representative series

``` r

calib_years <- instrumental_raw$Time

rp <- reduce_proxies(
  proxy_matrix = proxy_matrix[, -1],
  years        = proxy_matrix$year,
  temp_calib   = instrumental_raw$Anomaly,
  calib_years  = calib_years,
  method       = "PCR",
  # `chunk` must be large enough that every segment (including the last,
  # shortest one) still fully contains the calibration window (1850-2000,
  # 151 years). This toy example only has 301 years total, so we use a
  # single segment (chunk >= length(years)) rather than the package
  # default of 250, which would create a too-short trailing segment here.
  chunk        = nrow(proxy_matrix)
)
rp
#> <baygmst_rp>
#>   Method:   PCR
#>   Years:    1700 to 2000 (n = 301)
#>   Segments: 1
#>   Proxies per segment: 6
#>   Pass to fit_baygmst() directly, or use x$composite$RP1.
```

`rp` can be passed to
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
directly (it is aligned to the model’s year grid by calendar year, via
[`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md));
the raw series lives in `rp$composite$RP1` if you need it. You do not
have to use
[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
at all – any reduced or composited proxy series works here; see “Using
composites from other packages” below.

## 3. Transform the radiative forcings

``` r

forcing <- transform_forcings(
  G = forcing_raw$CO2,
  V = forcing_raw$volcanic,
  S = forcing_raw$solar
)
```

## 4. Align everything onto the common year grid and fit the model

The proxy series is aligned automatically (by calendar year) when you
pass the
[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
output itself; the instrumental series still needs aligning onto `years`
by hand, with `NA` marking the pre-instrumental years the model will
reconstruct:

``` r

instrumental_T <- instrumental_raw$Anomaly[match(years, instrumental_raw$Time)]
```

Fitting requires a working CmdStan installation (via / ); see
[`vignette("cmdstanr", package = "cmdstanr")`](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
if you don’t have one yet. The chunk below is skipped if none is found.

``` r

fit <- fit_baygmst(
  proxy          = rp,
  instrumental_T = instrumental_T,
  forcing_G      = forcing$G,
  forcing_V      = forcing$V,
  forcing_S      = forcing$S,
  years          = years,
  chains         = 2,
  parallel_chains = 2,
  iter_warmup    = 300,
  iter_sampling  = 1000,
  seed           = 1
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 1300 [  0%]  (Warmup) 
#> Chain 2 Iteration:    1 / 1300 [  0%]  (Warmup) 
#> Chain 1 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 1 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 2 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 1 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 1 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 2 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 2 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 2 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 2 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 2 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 2 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 2 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 1 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 2 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 finished in 1.4 seconds.
#> Chain 2 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 2 finished in 1.9 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.7 seconds.
#> Total execution time: 2.1 seconds.
summary(fit)
#> Posterior parameter summary:
#>    variable        mean      median          sd         mad          q5
#> 1    alpha0 -0.06483890 -0.06479368 0.012864394 0.012919098 -0.08616641
#> 2    alpha1  0.26192061  0.26065676 0.040014450 0.038979978  0.19558742
#> 3     phi_R  0.42065623  0.42223359 0.056968593 0.058957976  0.32676857
#> 4     phi_T  0.55946381  0.55955643 0.067630249 0.067263212  0.45444761
#> 5     beta0 -0.20119179 -0.20130989 0.034409116 0.034346371 -0.25773441
#> 6     betaG  0.27374461  0.27326430 0.051265366 0.052224778  0.18955831
#> 7     betaS -0.01146663 -0.01275857 0.050571557 0.051397436 -0.09430322
#> 8     betaV  0.06126381  0.06103615 0.017320465 0.017011913  0.03283503
#> 9   sigma_y  0.10387777  0.10362326 0.005999494 0.005874133  0.09449950
#> 10  sigma_z  0.11275083  0.11263524 0.004739298 0.004792831  0.10534086
#>            q95      rhat  ess_bulk  ess_tail
#> 1  -0.04345803 1.0003871 2324.5065 1430.0226
#> 2   0.32906312 1.0011778 1010.4032 1413.9274
#> 3   0.51241216 1.0005227 1327.0289 1628.7055
#> 4   0.67551736 1.0035327  595.4946  647.7400
#> 5  -0.14381323 1.0018334  624.9322  996.4503
#> 6   0.35909683 1.0012401  849.9260 1200.6283
#> 7   0.07454485 0.9995202  799.0300 1010.7448
#> 8   0.09076808 0.9996033 1016.0437 1378.6688
#> 9   0.11424406 1.0069576  728.4398 1349.8923
#> 10  0.12087903 0.9998866 2135.5942 1236.4933
#> 
#> Instrumental-period fit performance:
#>      RMSE       MAE        Bias Correlation        R2 R2_detrended
#>  0.130057 0.1061169 -0.02977442   0.8540961 0.7294801    0.3461799
```

## 5. Visualize the reconstruction

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fitted
model draws the summary figure used in the PReSto manuscript: the
reconstruction on top and the posterior densities of the model
parameters below (this needs the ‘patchwork’ package).

``` r

plot(fit)
```

![](baygmst-intro_files/figure-html/plot-fit-1.png)

The panels are also available individually:

``` r

plot_reconstruction(fit)
```

![](baygmst-intro_files/figure-html/plot-reconstruction-1.png)

``` r

plot_trace(fit)
```

![](baygmst-intro_files/figure-html/plot-trace-1.png)

``` r

plot_posterior_densities(fit)
```

![](baygmst-intro_files/figure-html/plot-posteriors-1.png)

## 6. Cross-validate

``` r

cv <- cv_baygmst(
  proxy          = rp,
  instrumental_T = instrumental_T,
  forcing_G      = forcing$G,
  forcing_V      = forcing$V,
  forcing_S      = forcing$S,
  years          = years,
  nfold          = 3,
  chains         = 2,
  parallel_chains = 2,
  iter_warmup    = 300,
  iter_sampling  = 1000,
  seed           = 1
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 1300 [  0%]  (Warmup) 
#> Chain 2 Iteration:    1 / 1300 [  0%]  (Warmup)
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: normal_lpdf: Scale parameter is 0, but must be positive! (in '/tmp/Rtmpj67fZw/model-1c1842aaadfe.stan', line 75, column 2 to column 40)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2
#> Chain 1 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 2 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 1 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 2 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 2 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 1 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 2 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 2 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 2 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 2 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 2 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 2 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 2 finished in 2.2 seconds.
#> Chain 1 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 1 finished in 2.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.3 seconds.
#> Total execution time: 2.4 seconds.
#> 
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 1300 [  0%]  (Warmup) 
#> Chain 2 Iteration:    1 / 1300 [  0%]  (Warmup)
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: normal_lpdf: Scale parameter is 0, but must be positive! (in '/tmp/Rtmpj67fZw/model-1c1842aaadfe.stan', line 75, column 2 to column 40)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2
#> Chain 2 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 2 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 1 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 2 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 2 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 2 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 2 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 2 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 2 finished in 1.1 seconds.
#> Chain 1 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 1 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 1 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 1 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 1 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 1 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 1 finished in 2.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.9 seconds.
#> Total execution time: 2.7 seconds.
#> 
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:    1 / 1300 [  0%]  (Warmup) 
#> Chain 2 Iteration:    1 / 1300 [  0%]  (Warmup)
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: normal_lpdf: Scale parameter is 0, but must be positive! (in '/tmp/Rtmpj67fZw/model-1c1842aaadfe.stan', line 75, column 2 to column 40)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2
#> Chain 1 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 2 Iteration:  100 / 1300 [  7%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  200 / 1300 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 1 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 2 Iteration:  300 / 1300 [ 23%]  (Warmup) 
#> Chain 2 Iteration:  301 / 1300 [ 23%]  (Sampling) 
#> Chain 1 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 2 Iteration:  400 / 1300 [ 30%]  (Sampling) 
#> Chain 1 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 2 Iteration:  500 / 1300 [ 38%]  (Sampling) 
#> Chain 1 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 2 Iteration:  600 / 1300 [ 46%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 2 Iteration:  700 / 1300 [ 53%]  (Sampling) 
#> Chain 1 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 2 Iteration:  800 / 1300 [ 61%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 2 Iteration:  900 / 1300 [ 69%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1300 [ 76%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 1300 [ 84%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 1300 [ 92%]  (Sampling) 
#> Chain 1 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 1300 [100%]  (Sampling) 
#> Chain 1 finished in 2.2 seconds.
#> Chain 2 finished in 2.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.2 seconds.
#> Total execution time: 2.3 seconds.
cv
#> <baygmst_cv>
#>   Folds:  3
#>   R2 (posterior mean, per fold): -3.73, 0.25, -0.29
#>   Median R2: -0.29
plot_cv(cv)
```

![](baygmst-intro_files/figure-html/cv-1.png)

## 7. Using composites from other packages

[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)
is one way to produce the model’s proxy input, but not the only one.
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)’s
proxy equation estimates its own intercept and slope against
temperature, so any reduced or composited proxy series on a reasonable
scale works – a standardized (z-score) composite or a
temperature-anomaly-scale series is appropriate; a raw-unit series with
large numerical magnitude should be standardized first.

[`as_baygmst_proxy()`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md)
normalizes the common shapes such a series arrives in. In particular, it
accepts composite objects of the shape produced by compositing packages
in the LiPD ecosystem (such as the in-development
[compositeR](https://github.com/nickmckay/compositeR) package, whose
`compositeEnsembles2()` is the intended primary pathway for building
BayGMST inputs from paleoclimate data): a list carrying `$ages` (in
years BP, i.e. before 1950 CE) and `$composite` (rows = time steps,
columns = ensemble members). Here is a synthetic stand-in:

``` r

set.seed(1)
composite_ens <- list(
  ages      = seq(250, -50, by = -1),               # 1700-2000 CE, in yr BP
  composite = matrix(rnorm(301 * 50, sd = 0.2), 301, 50) +
    seq(-0.2, 0.4, length.out = 301)                # 50 ensemble members
)

proxy_input <- as_baygmst_proxy(composite_ens)
proxy_input
#> <baygmst_proxy>
#>   Source:  composite (50 members, median; ages BP)
#>   Years:   1700 to 2000 CE (n = 301)
#>   Ensemble: collapsed from 50 members
#>   Range:   [-0.247, 0.487]
```

Ages are converted to calendar years and the ensemble is collapsed to
its median series (compositing uncertainty is therefore not yet
propagated into the posterior; the object records `n_ensemble` so this
can become an ensemble-aware pathway in a future release). The result
can be passed straight to `fit_baygmst(proxy = proxy_input, ...)`, which
aligns it to `years` by calendar year and errors informatively if any
requested year is not covered. Bare numeric vectors, `data.frame`s with
year and value columns, and time-by-ensemble matrices work the same way
– see
[`?as_baygmst_proxy`](https://paleopresto.github.io/BayGMST_R/reference/as_baygmst_proxy.md).

## Known limitations

See `NEWS.md` for the modeling decisions and known limitations of this
release (for example, `reduce_proxies(method = "SIR")` is experimental).
The data in this vignette are synthetic, so its reconstructions should
not be over-interpreted.
