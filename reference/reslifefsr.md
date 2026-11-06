# Calculate Residual Life Values Using a flexsurvreg() Object

Calculate residual life values using a fitted `flexsurvreg` object.
Allows specifying covariate profiles via `newdata`.

## Usage

``` r
reslifefsr(obj, life, p = 0.5, type = "mean", newdata = data.frame())
```

## Arguments

- obj:

  A fitted `flexsurvreg` object from package flexsurv.

- life:

  Scalar numeric value of elapsed time at which to calculate residual
  life.

- p:

  Percentile for quantile residual life. Default is 0.5.

- type:

  Type of residual life to return: `"mean"`, `"median"`, `"quantile"`,
  `"percentile"`, or `"all"`. Default is `"mean"`.

- newdata:

  A data frame containing covariate values to evaluate. Default is empty
  (evaluates baseline).

## Value

A numeric vector (if `type` is `"mean"`, `"median"`, or `"quantile"`) or
a data frame (if `type == "all"`), with one row per covariate profile in
`newdata`.

## References

Jackson CH (2016). “flexsurv: a platform for parametric survival
modeling in R.” Journal of statistical software, 70.

Poynor V (2010). “Bayesian inference for mean residual life functions in
survival analysis.” Masters diss., University of California, Santa Cruz.

Prentice RL (1975). “Discrimination among some parametric models.”
Biometrika, 62(3), 607–614.

Stacy EW (1962). “A generalization of the gamma distribution.” The
Annals of mathematical statistics, pp. 1187–1192.

## Examples

``` r
library(flexsurv)
fitg <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "gengamma")
reslifefsr(obj = fitg, life = 6, p = 0.75, type = "all")
#>       mean   median percentile
#> 1 913.9498 903.1636   2856.237

fitg2 <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
df_new <- data.frame(age = c(50, 65))
reslifefsr(obj = fitg2, life = 3, type = "median", newdata = df_new)
#> [1] 1368.6166  366.1353
```
