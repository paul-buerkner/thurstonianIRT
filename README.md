thurstonianIRT
==============

[![Build
Status](https://travis-ci.org/paul-buerkner/thurstonianIRT.svg?branch=master)](https://travis-ci.org/paul-buerkner/thurstonianIRT)
[![CRAN
Version](http://www.r-pkg.org/badges/version/thurstonianIRT)](https://cran.r-project.org/package=thurstonianIRT)

Overview
--------

The **thurstonianIRT** package allows to fit various models from Item
Response Theory (IRT) for forced-choice questionnaires, most notably the
Thurstonian IRT model originally proposed by (Brown & Maydeu-Olivares,
2011). The key characteristic of forced-choice questionnaires is that
participants cannot endorse all items at the same time and instead have
to make a comparative judgement between two or more items. Such a format
comes with the hope of providing more valid inference in situation where
participants have motivation to not answer honestly (e.g., personnel
selection), but instead respond in a way that appears favorable in the
given situation. Whether forced-choice questionnaires and the
corresponding IRT models live up to this hope remains a topic of debate
(Bürkner, Schulte, & Holling, 2019) but it is in any case necessary to
provide software for fitting these statistical models both for practical
and research purposes.

In the original formulation, the Thurstonian IRT model works on
dichotomous pairwise comparisons and models the probability of endorsing
one versus the other item. This probability depends on parameters
related to the items under comparison as well as on parameters related
to the participants’ latent traits which are assumed to be measured by
the items. For more details see Brown and Maydeu-Olivares (2011), Brown
and Maydeu-Olivares (2012), and Bürkner et al. (2019).

How to use thurstonianIRT
-------------------------

``` r
library(thurstonianIRT)
```

As a simple example consider a data set of 4 blocks of items each
containing 3 items (i.e., triplets) for which we have data of 100
participants:

``` r
data("triplets")
head(triplets)
#>   i1i2 i1i3 i2i3 i4i5 i4i6 i5i6 i7i8 i7i9 i8i9 i10i11 i10i12 i11i12
#> 1    1    0    0    1    0    0    1    1    1      0      1      1
#> 2    0    0    1    0    0    0    0    0    1      0      0      0
#> 3    0    0    1    0    0    1    0    1    1      0      0      0
#> 4    0    0    1    1    1    0    1    1    0      0      0      0
#> 5    1    1    1    0    0    1    1    1    0      1      0      0
#> 6    1    1    1    0    0    1    1    0    0      0      1      1
```

In the data set, a `1` indicates that the first item has been selected
over the second item while a `0` indicates that the second items has
been selected over the first item. In order to fit a Thurstonian IRT
model on this data, we have to tell **thurstonianIRT** about the block
structure of the items, the traits on which the items load, and the sign
of these loadings, that is, whether items have been inverted. For the
present data, we specify this as follows:

``` r
blocks <-
  set_block(c("i1", "i2", "i3"), traits = c("t1", "t2", "t3"),
            signs = c(1, 1, 1)) +
  set_block(c("i4", "i5", "i6"), traits = c("t1", "t2", "t3"),
            signs = c(-1, 1, 1)) +
  set_block(c("i7", "i8", "i9"), traits = c("t1", "t2", "t3"),
            signs = c(1, 1, -1)) +
  set_block(c("i10", "i11", "i12"), traits = c("t1", "t2", "t3"),
            signs = c(1, -1, 1))
```

Next, we transform the data into a format that **thurstonianIRT**
understands.

``` r
triplets_long <- make_TIRT_data(
  data = triplets, blocks = blocks, direction = "larger",
  format = "pairwise", family = "bernoulli", range = c(0, 1)
)
head(triplets_long)
#> # A tibble: 6 x 11
#>   person block comparison itemC trait1 trait2 item1 item2 sign1 sign2 response
#>    <int> <int>      <int> <dbl> <fct>  <fct>  <fct> <fct> <dbl> <dbl>    <dbl>
#> 1      1     1          1     1 t1     t2     i1    i2        1     1        1
#> 2      2     1          1     1 t1     t2     i1    i2        1     1        0
#> 3      3     1          1     1 t1     t2     i1    i2        1     1        0
#> 4      4     1          1     1 t1     t2     i1    i2        1     1        0
#> 5      5     1          1     1 t1     t2     i1    i2        1     1        1
#> 6      6     1          1     1 t1     t2     i1    i2        1     1        1
```

Finally, we can fit the model using several model fitting engines.
Currently supported are [Stan](https://mc-stan.org/),
[lavaan](http://lavaan.ugent.be/), and
[Mplus](http://www.statmodel.com/). Here, we choose Stan to fit the
Thurstonian IRT model in a Bayesian framework.

``` r
fit <- fit_TIRT_stan(triplets_long, chains = 1)
```

As basic summary and convergence checks can be obtained via

``` r
print(fit)
```

Finally, we obtain predictions of participants’ trait scores in a tidy
data format via

``` r
pr <- predict(fit)
head(pr)
#> # A tibble: 6 x 6
#>      id trait estimate    se lower_ci upper_ci
#>   <int> <chr>    <dbl> <dbl>    <dbl>    <dbl>
#> 1     1 t1      0.360  0.544   -0.623    1.50 
#> 2     1 t2     -1.49   0.567   -2.64    -0.454
#> 3     1 t3      0.0196 0.603   -1.06     1.29 
#> 4     2 t1     -0.868  0.536   -1.98     0.165
#> 5     2 t2      0.654  0.583   -0.409    1.91 
#> 6     2 t3      0.383  0.607   -0.723    1.62
```

The thurstonianIRT package not only comes with model fitting functions
but also with the possibility to simulate data from Thurstonian IRT
models. Below we simulate data with a very similar structure to the
`triplets` data set we have used above.

``` r
sim_data <- sim_TIRT_data(
  npersons = 100,
  ntraits = 3,
  nblocks_per_trait = 4,
  gamma = 0,
  lambda = runif(12, 0.5, 1),
  Phi = diag(3)
)
#> Computing standardized psi as 1 - lambda^2
head(sim_data)
#> # A tibble: 6 x 19
#>   person block comparison itemC trait1 trait2 item1 item2 sign1 sign2 gamma lambda1
#>    <int> <int>      <int> <dbl>  <int>  <int> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1      1     1          1     1      3      1     1     2     1     1     0   0.970
#> 2      2     1          1     1      3      1     1     2     1     1     0   0.970
#> 3      3     1          1     1      3      1     1     2     1     1     0   0.970
#> 4      4     1          1     1      3      1     1     2     1     1     0   0.970
#> 5      5     1          1     1      3      1     1     2     1     1     0   0.970
#> 6      6     1          1     1      3      1     1     2     1     1     0   0.970
#> # ... with 7 more variables: lambda2 <dbl>, psi1 <dbl>, psi2 <dbl>, eta1 <dbl>,
#> #   eta2 <dbl>, mu <dbl>, response <int>
```

The structure of the data is the same as what we have obtained via the
`make_TIRT_data` function above and can readily be passed to the model
fitting functions.

How to install thurstonianIRT
-----------------------------

The current developmental version can be downloaded from github via

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("paul-buerkner/thurstonianIRT")
```

References
----------

Brown, A., & Maydeu-Olivares, A. (2011). Item response modeling of
forced-choice questionnaires. *Educational and Psychological
Measurement*, 71(3), 460-502.

Brown, A., & Maydeu-Olivares, A. (2012). Fitting a Thurstonian IRT model
to forced-choice data using Mplus. *Behavior Research Methods*, 44(4),
1135-1147.

Bürkner P. C., Schulte N., & Holling H. (2019). On the Statistical and
Practical Limitations of Thurstonian IRT Models. *Educational and
Psychological Measurement.*
<a href="doi:10.1177/0013164419832063" class="uri">doi:10.1177/0013164419832063</a>
