# Coerce a reduced-proxy series into BayGMST's input format

Normalizes the many shapes a reduced/composited proxy series can arrive
in – the output of
[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md),
a plain `data.frame`, a bare vector, a time-by-ensemble matrix, or a
composite object from another package (such as a `paleoComposite` from
the in-development 'compositeR' package) – into a single, year-aware
`"baygmst_proxy"` object that
[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md)
and
[`cv_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/cv_baygmst.md)
accept directly and align by calendar year (rather than by position).

## Usage

``` r
as_baygmst_proxy(x, ...)

# S3 method for class 'baygmst_proxy'
as_baygmst_proxy(x, ...)

# S3 method for class 'baygmst_rp'
as_baygmst_proxy(x, ...)

# S3 method for class 'numeric'
as_baygmst_proxy(x, years, ...)

# S3 method for class 'matrix'
as_baygmst_proxy(x, years, collapse = c("median", "mean"), ...)

# S3 method for class 'data.frame'
as_baygmst_proxy(x, ...)

# S3 method for class 'list'
as_baygmst_proxy(
  x,
  age_units = c("BP", "CE"),
  collapse = c("median", "mean"),
  ...
)
```

## Arguments

- x:

  The object to coerce. Supported inputs:

  `baygmst_rp`

  :   The output of
      [`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md);
      uses its `$composite` series.

  `data.frame`

  :   Must contain a year column (named `year` or `Year`) and exactly
      one other numeric column (e.g. the `example_proxy.csv` shipped in
      `inst/extdata`, with columns `Year, RP1`).

  `numeric` vector

  :   One value per element of `years`, which must also be supplied.

  `matrix` (or `data.frame` of ensemble columns)

  :   Rows are time steps (matching `years`), columns are ensemble
      members; collapsed to a single series with `collapse`.

  composite `list`

  :   Any list with elements `$ages` and `$composite` (rows = time steps
      of `$ages`, columns = ensemble members), the shape produced by
      compositing packages in the LiPD ecosystem (e.g.
      `compositeR::compositeEnsembles2()`). Ages are assumed to be in
      years BP (before present, i.e. before 1950 CE) unless
      `age_units = "CE"`. The ensemble is collapsed with `collapse`.

- ...:

  Passed between methods.

- years:

  Numeric vector of calendar years (CE). Required when `x` carries no
  time axis of its own (bare vectors and matrices); ignored otherwise.

- collapse:

  How to collapse an ensemble (matrix or composite list) to a single
  series: `"median"` (default) or `"mean"`, taken across ensemble
  members within each time step, ignoring `NA`s. Time steps whose
  collapsed value is not finite (e.g. bins no record covers) are dropped
  with a message.

- age_units:

  For composite lists only: `"BP"` (default; ages are years before 1950
  CE and are converted via `year = 1950 - age`) or `"CE"` (ages are
  already calendar years).

## Value

An object of class `"baygmst_proxy"`: a `data.frame` with columns `year`
(ascending) and `proxy`, plus attributes `source` (a short label for
where the series came from) and `n_ensemble` (the number of ensemble
members collapsed; `1L` for deterministic inputs).

## Scale expectations

The BayGMST 'Stan' model regresses the proxy series on temperature with
an estimated intercept and slope, so the proxy does not need to be
pre-calibrated to temperature units. The (weakly informative) priors do,
however, assume the series is on roughly unit scale: a standardized
(z-score) composite or a temperature-anomaly-scale series is
appropriate; a raw-unit series (per mil, mm, ...) with large numerical
magnitude should be standardized first.

## Uncertainty

In this version the ensemble dimension is collapsed before fitting, so
compositing/ensemble uncertainty is not propagated into the posterior. A
future release may add an ensemble-aware fitting pathway; the
`"baygmst_proxy"` object records `n_ensemble` so existing code will not
need to change shape when it does.

## See also

[`fit_baygmst()`](https://paleopresto.github.io/BayGMST_R/reference/fit_baygmst.md),
[`reduce_proxies()`](https://paleopresto.github.io/BayGMST_R/reference/reduce_proxies.md)

## Examples

``` r
# from a bare vector
p <- as_baygmst_proxy(rnorm(100), years = 1901:2000)

# from a data.frame like inst/extdata/example_proxy.csv
df <- data.frame(Year = 1901:2000, RP1 = rnorm(100))
p <- as_baygmst_proxy(df)

# from a composite ensemble with ages in yr BP (compositeR-style)
comp <- list(
  ages = seq(950, 0, by = -10),               # 1000-1950 CE
  composite = matrix(rnorm(96 * 20), 96, 20)  # 20 ensemble members
)
p <- as_baygmst_proxy(comp)
head(p)
#> <baygmst_proxy>
#>   Source:  composite (20 members, median; ages BP)
#>   Years:   1000 to 1050 CE (n = 6)
#>   Ensemble: collapsed from 20 members
#>   Range:   [-0.422, 0.665]
```
