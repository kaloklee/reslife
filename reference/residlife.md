# Calculating Residual Life Values

Calculates residual life values over a range of values. Allows the user
to specify the distribution and the parameters.

## Usage

``` r
residlife(values, distribution, parameters, p = 0.5, type = "mean")
```

## Arguments

- values:

  Range of values over which residual life is calculated.Usually given
  as a vector.

- distribution:

  Name of the distribution. Needs to be one of 'weibull', 'gompertz',
  'gamma', 'gengamma.orig', 'exponential', 'lnorm', or 'llogis', 'genf',
  or 'genf.orig'.

- parameters:

  Parameters of the survival function. Needs to be inputted in order as
  a vector, with the name of the parameter included.

- p:

  Percentile at which to calculate residual life. Default is .5.

- type:

  Type of residual life outputted. Must be "mean", "median",
  "percentile", or "all". Default is "mean".

## Value

The residual life for a specified sequence of values.

## References

Jackson CH (2016). “flexsurv: a platform for parametric survival
modeling in R.” Journal of statistical software, 70.

Poynor V (2010). “Bayesian inference for mean residual life functions in
survival analysis.” Masters diss., University of California, Santa Cruz.

Prentice RL (1975). “Discrimination among some parametric models.”
Biometrika, 62(3), 607–614.

Stacy EW (1962). “A generalization of the gamma distribution.” The
Annals of mathematical statistics, pp. 1187–1192

## Examples

``` r
residlife(values = 0:60, distribution= 'weibull', parameters = c(shape = 1.2, scale = 3))
#>  [1] 2.821968 2.526989 2.375535 2.268927 2.186629 2.119796 2.063703 2.015509
#>  [9] 1.973364 1.935994 1.902487 1.872166 1.844513 1.819126 1.795687 1.773938
#> [17] 1.753668 1.734702 1.716895 1.700123 1.684281 1.669280 1.655040 1.641494
#> [25] 1.628583 1.616253 1.604460 1.593162 1.582321 1.571906 1.561886 1.552236
#> [33] 1.542930 1.533948 1.525268 1.516874 1.508748 1.500874 1.493240 1.485832
#> [41] 1.478638 1.471646 1.464848 1.458233 1.451792 1.445518 1.439402 1.433437
#> [49] 1.427617 1.421935 1.416386 1.410964 1.405663 1.400479 1.395408 1.390444
#> [57] 1.385584 1.380824 1.376159 1.371588 1.367105
residlife(values = 15:35, distribution= 'gamma', parameters =  c(shape = 1.2, rate =  1.7),
p = .25, type ='all')
#>         mean    median percentile
#> 1  0.5925534 0.4107970  0.1705060
#> 2  0.5922994 0.4106130  0.1704284
#> 3  0.5920736 0.4104498  0.1703597
#> 4  0.5918716 0.4103042  0.1702985
#> 5  0.5916899 0.4101733  0.1702435
#> 6  0.5915254 0.4100552  0.1701939
#> 7  0.5913760 0.4099480  0.1701489
#> 8  0.5912395 0.4098502  0.1701079
#> 9  0.5911144 0.4097607  0.1700703
#> 10 0.5909993 0.4096785  0.1700359
#> 11 0.5908931 0.4096027  0.1700041
#> 12 0.5907947 0.4095326  0.1699747
#> 13 0.5907034 0.4094675  0.1699475
#> 14 0.5906183 0.4094070  0.1699222
#> 15 0.5905389 0.4093506  0.1698986
#> 16 0.5904647 0.4092979  0.1698765
#> 17 0.5903951 0.4092485  0.1698558
#> 18 0.5903297 0.4092021  0.1698364
#> 19 0.5902681 0.4091584  0.1698182
#> 20 0.5902101 0.4091173  0.1698010
#> 21 0.5901553 0.4090785  0.1697848
```
