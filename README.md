
<!-- README.md is generated from README.Rmd. Please edit that file -->

# covcpd

<!-- badges: start -->
<!-- badges: end -->

For single-cell RNA sequencing data with continuous or ambiguous cell
states, we develop a covariance-based change point detection procedure
to infer the discrete subgroups by utilizing the continuous pseudotime
of single cells. Little research suggests whether and how well the
existing multivariate change point detection methods work for scRNA-seq
data. Hence, popular existing methods are first benchmarked and
evaluated in the simulation study and are shown to be powered for
detecting mean but not covariance changes. We then propose the algorithm
covcpd (Covariance Structure Change Point Detection), which
hierarchically partitions single-cell samples into homogeneous network
groups by utilizing a covariance equality testing statistic.

<figure>
<img src="vignettes/covcpd_F1.PNG" style="width:80.0%" alt="Overview" />
<figcaption aria-hidden="true">Overview</figcaption>
</figure>

## Installation

You can install the development version of covcpd like so:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("meichendong/covcpd")
```

## Example

The input data matrix can be raw data or other formats. Samples should
be ordered by pseudotime or other prior information. Pseudotime can be
calculated by trajectory inference methods. Please refer to the website
[dynverse](https://dynverse.org/) to calculate pseudotime for samples.

The matrix should be genes by sample. `k` is the user-defined minimum
number of cells per segment. `maxseg` is the user-defined maximum number
segments for the dataset. `search_by` is the number of cells/points
apart set per candidate point. `siglvl` significance level. `nperm`
number of permutations when learning the Gumbel distribution parameters.

<!--The following example will take around 4 mins.-->

``` r
library(covcpd)
## basic example code
set.seed(1)

# simulate data matrix with one change point.
# example: the first 50 genes have different structure, the rest of genes have the same structures.
nsample = 100

nivec.list.diff <- list(nivec= c(40,60,100),nivec2 = c(100,100))
diffblk = list( 1:2,1)
sigma.list.1 <- generateSigmaList(nivec.list = nivec.list.diff, structure = "Diff S, Identical W", diffblk = diffblk)
count1 <- CountMap5(sigma = sigma.list.1[[1]], ngene = 200, n = nsample)
count2 <- CountMap6(sigma = sigma.list.1[[2]], ngene = 200, n = nsample, scmeanvec = count1$thetaj, pijvec = count1$pij)
y = rbind(count1$count,
          count2$count)

# Run covcpd
t1 = proc.time()
result = covcpd(X = t(y), k=50, maxseg = 3, search_by = NULL, siglvl = 0.05, nperm = 1000)
t2 = proc.time()

# result
result$cps

t2-t1

# > result$cps
# [1]   1  47 100
# > t2-t1
# user  system elapsed 
# 1405.11  102.25 1568.05 
```

Simulation and Real data analysis R code can be found at the R folder on
github.
