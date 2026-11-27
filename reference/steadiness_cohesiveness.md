# Steadiness and Cohesiveness of Groups

Steadiness checks whether groups sampled in the embedding are coherent
in the reference. Cohesiveness checks the reverse direction. Larger
scores indicate better preservation. Class labels are unnecessary.

## Usage

``` r
steadiness_cohesiveness(Xin, Xout, k = NULL, n_iter = 150L, ret_extra = FALSE)
```

## Arguments

- Xin:

  Reference coordinates, as a dense numeric matrix or data frame with
  observations in rows.

- Xout:

  Embedding coordinates for the same observations in the same order.
  Nonnumeric data-frame columns are ignored in both inputs.

- k:

  Number of other Euclidean neighbors. The default is the integer part
  of the square root of the number of rows. Distance ties follow row
  order.

- n_iter:

  Number of independently grown groups per direction. Use
  [`set.seed()`](https://rdrr.io/r/base/Random.html) to reproduce a
  draw.

- ret_extra:

  If true, also return scores for the reliability map.

## Value

A named vector of `steadiness` and `cohesiveness`. With
`ret_extra = TRUE`, returns `scores` and a row-aligned data frame
`pointwise` with both map values. Each map accumulates weighted losses
by partner, averages repeated events with that partner, then rescales
across observations. Map values depend on the sampled groups. A point
with no recorded loss receives one. Map means do not recover the global
scores. The other direction's global loss can amplify a map value's
distance from one, so map values can fall below zero.

## Details

The walk follows links with probability given by normalized shared
neighbor similarity, using rank weights from k + 1 down to 2. A walk
accepts about 0.3 times the number of observations in links, or stops
when its queue is empty. The opposite space partitions the sampled group
with
[`dbscan::hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.html)
(minimum points 5). Each noise point forms a singleton. The
dissimilarity between two subgroups is the reciprocal of the sum of
their mean similarity and 0.1.

Each score subtracts the weighted mean positive directional distortion
from one. Weights are the products of subgroup sizes. A direction with
no observed positive distortion receives one. Both `dbscan` and `Matrix`
are required.

## Scaling

The calculation builds dense `N` by `N` distances and shared-neighbor
similarities internally, even for global scores. A single double matrix
at 5,000 rows uses about 191 MiB. Several coexist, and
`ret_extra = TRUE` adds square matrices of partner losses and visit
counts. For large datasets, evaluate the same representative subset in
both spaces and report its size. Reducing `n_iter` saves repeated
clustering work. The dense graphs still cover every supplied
observation.

## References

Jeon et al. (2022). Measuring and Explaining the Inter-Cluster
Reliability of Multidimensional Projections.
[doi:10.1109/TVCG.2021.3114833](https://doi.org/10.1109/TVCG.2021.3114833)
The [authors' code](https://github.com/hj-n/steadiness-cohesiveness)
accompanies the paper. For shared-neighbor similarities and pointwise
map weights, Quadra follows [ZADU 0.5.4](https://github.com/hj-n/zadu).
It uses R's random draws and the dbscan package's HDBSCAN partition.

## Examples

``` r
x <- as.matrix(iris[, -5])
set.seed(42)
steadiness_cohesiveness(x, prcomp(x)$x[, 1:2], n_iter = 15)
#>   steadiness cohesiveness 
#>    0.9068442    0.8811552 
```
