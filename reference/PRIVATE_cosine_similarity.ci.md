# Confidence interval of cosine similarity

Computes the confidence interval of a cosine similarity coefficient by
bootstraping. Adapted from the implementation of spearman.ci in
RVAidemoire Version 0.9-83-7.

## Usage

``` r
PRIVATE_cosine_similarity.ci(var1, var2, nrep = 1000, conf.level = 0.95)
```

## Arguments

- var1:

  numeric vector (first variable).

- var2:

  nuermic verctor (second variable).

- nrep:

  number of replicates for bootstrapping.

- conf.level:

  confidence level of the interval.

## Value

description method name of the test.

data.name a character string giving the name(s) of the data.

conf.level confidence level.

rep number of replicates.

estimate Spearman's rank correlation coefficient.

conf.int confidence interval.
