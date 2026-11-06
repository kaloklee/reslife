# reslife

The repository for the R package [reslife](https://cran.r-project.org/web/packages/reslife/index.html).

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/reslife)](https://cran.r-project.org/web/packages/reslife/index.html)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

`reslife` computes residual life functions (mean residual life, median residual life, and quantile residual life) for parametric survival distributions and fitted `flexsurvreg` survival models.

## Installation

### Stable release from CRAN

```r
install.packages("reslife")
```

### Development version from GitHub

```r
# install.packages("devtools")
devtools::install_github("kaloklee/reslife")
```

> **System Dependency Note:** If you encounter an installation failure related to the dependency package `gsl` on Debian, Ubuntu, or LMDE, install `libgsl-dev` via your terminal:
> ```bash
> sudo apt install libgsl-dev
> ```

---

## Documentation & Website

* **Package Website & Vignette:** [https://kaloklee.github.io/reslife/](https://kaloklee.github.io/reslife/)
* **Methodology Paper:** Wang, Z., Crawford, A., Lee, K. L., Lu, L., & Jaganathan, S. (2025). *reslife: Residual Lifetime Analysis Tool in R*. arXiv preprint [arXiv:2308.07410](https://arxiv.org/abs/2308.07410).

---

## Quick Example

```r
library(reslife)
library(flexsurv)
library(survival)

# 1. Evaluate standalone distribution parameters
pars_weibull <- c(shape = 1.27, scale = 6.19)
residlife(values = c(1, 3, 5), distribution = "weibull", parameters = pars_weibull, type = "all")

# 2. Evaluate from a fitted flexsurvreg model
fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "weibull")
reslifefsr(obj = fit, life = 2, type = "all")
```