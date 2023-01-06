# PublicationBias

<!-- badges: start -->
  [![R-CMD-check](https://github.com/mayamathur/PublicationBias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mayamathur/PublicationBias/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`PublicationBias` provides sensitivity analyses for publication bias in meta-analyses (per [Mathur & VanderWeele, 2020](https://osf.io/s9dp6)).

## Installation

You can install PublicationBias from CRAN with:
```
install.packages("PublicationBias")
```

You can install the development version of PublicationBias from [GitHub](https://github.com/) with:
``` r
# install.packages("devtools")
devtools::install_github("mayamathur/PublicationBias")
```

## Example

Start by generating some example data from the `metafor` package.

``` r
library(PublicationBias)
dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                       di = cneg, data = dat.bcg)
```

Calculate the meta-analytic effect size estimate, correcting for an assumed
selection ratio of 5 (i.e., affirmative results are 5x more likely to be
published than nonaffirmative ones).

``` r
pubbias_meta(yi = dat$yi, vi = dat$vi, selection_ratio = 5,
             model_type = "fixed", favor_positive = FALSE)
```

Calculate how high the selection ratio would need to be to attenuate the effect
size estimate to the null.

``` r
pubbias_svalue(yi = dat$yi, vi = dat$vi, q = 0,
               model_type = "fixed", favor_positive = FALSE)
```
