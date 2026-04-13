# kBET - k-nearest neighbour batch effect test

`kBET` runs a chi square test to evaluate the probability of a batch
effect.

## Usage

``` r
.kBET_fct(
  df,
  batch,
  k0 = NULL,
  knn = NULL,
  testSize = NULL,
  do.pca = TRUE,
  dim.pca = 50,
  heuristic = TRUE,
  n_repeat = 100,
  alpha = 0.05,
  addTest = FALSE,
  verbose = FALSE,
  plot = TRUE,
  adapt = TRUE
)
```

## Arguments

- df:

  dataset (rows: cells, columns: features)

- batch:

  batch id for each cell or a data frame with both condition and
  replicates

- k0:

  number of nearest neighbours to test on (neighbourhood size)

- knn:

  an n x k matrix of nearest neighbours for each cell (optional)

- testSize:

  number of data points to test, (10 percent sample size default, but at
  least 25)

- do.pca:

  perform a pca prior to knn search? (defaults to TRUE)

- dim.pca:

  if do.pca=TRUE, choose the number of dimensions to consider (defaults
  to 50)

- heuristic:

  compute an optimal neighbourhood size k (defaults to TRUE)

- n_repeat:

  to create a statistics on batch estimates, evaluate 'n_repeat' subsets

- alpha:

  significance level

- addTest:

  perform an LRT-approximation to the multinomial test AND a multinomial
  exact test (if appropriate)

- verbose:

  displays stages of current computation (defaults to FALSE)

- plot:

  if stats \> 10, then a boxplot of the resulting rejection rates is
  created

- adapt:

  In some cases, a number of cells do not contribute to any
  neighbourhood and this may cause an imbalance in observed and expected
  batch label frequencies. Frequencies will be adapted if adapt=TRUE
  (default).

## Value

list object

1.  `summary` - a rejection rate for the data, an expected rejection
    rate for random labeling and the significance for the observed
    result

2.  `results` - detailed list for each tested cells; p-values for
    expected and observed label distribution

3.  `average.pval` - significance level over the averaged batch label
    distribution in all neighbourhoods

4.  `stats` - extended test summary for every sample

5.  `params` - list of input parameters and adapted parameters,
    respectively

6.  `outsider` - only shown if `adapt=TRUE`. List of samples without
    mutual nearest neighbour:

    - `index` - index of each outsider sample)

    - `categories` - tabularised labels of outsiders

    - `p.val` - Significance level of outsider batch label distribution
      vs expected frequencies. If the significance level is lower than
      `alpha`, expected frequencies will be adapted

If the optimal neighborhood size (k0) is smaller than 10, NA is
returned.

## Author

Mariano Ruz Jurado (edited from: Maren Buettner)
