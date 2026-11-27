# Exact Rank Profiles from Euclidean Coordinates

Compare several embeddings against the same candidate neighbors. For
each query, compute its reference ranks once, score every embedding,
then release the ranks before moving to the next query.

## Usage

``` r
exact_rank_profiles(Xin, embeddings, k, mrre = FALSE, query_rows = NULL)
```

## Arguments

- Xin:

  Dense reference coordinates with observations in rows.

- embeddings:

  A nonempty list of dense coordinate matrices or data frames with the
  same observations in the same order as `Xin`. Pass `list(y)` for a
  single embedding. Nonnumeric data-frame columns are ignored. All
  coordinate values must be finite.

- k:

  Unique neighborhood sizes, each positive and less than half the number
  of observations. This range supports both rank penalties.

- mrre:

  Logical. Include mean relative rank errors (MRRE) for input neighbors
  and output neighbors? Defaults to `FALSE`.

- query_rows:

  Optional vector of unique, one-based row positions in `Xin` to score,
  in the requested order. The default (`NULL`) scores all rows. Each
  query still ranks all other observations, including rows outside the
  selection.

## Value

A list in embedding order, retaining the names of `embeddings`. Each
entry contains pointwise `trustworthiness` and `continuity` matrices,
the RNX `auc` vector, and the RNX `curve` matrix at the requested sizes.
The integer matrix `retained` counts neighbors shared by the reference
and embedding among the first `k` other observations. Divide a column by
its `k` to get the unadjusted fraction retained. Rows retain the
selected reference row names. Column means summarize the selected
queries. With `mrre = TRUE`, each entry also contains `mrre_input` and
`mrre_output` matrices at the requested sizes. Smaller errors indicate
better rank preservation.

## Details

Each row's candidate pool contains all other observations. Tied
Euclidean distances follow observation order, as with
[`pointwise_rank_preservation()`](https://jlmelville.github.io/quadra/reference/pointwise_rank_preservation.md)
and
[`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)
on matrices from [`stats::dist()`](https://rdrr.io/r/stats/dist.html).
RNX AUC uses every size from one through `n - 2`.

For MRRE, select the first `k` neighbors in one space. For each of these
neighbors, divide the absolute change in rank by its rank in the
selecting space. Sum and divide by `sum((n - 2 * r) / r)` over ranks 1
through `k`. Identical rankings give zero and reversed rankings give
one. For large neighborhoods, other permutations can exceed one under
this normalization.

## Scaling

Keep all `N` rows of the coordinates. Each of `q` selected queries
computes and sorts `N - 1` distances per space. Scoring all rows
requires quadratic distance work, and wide inputs make each distance
more expensive. Temporary rank workspace grows linearly with `N`. Output
size grows with the number of queries, requested neighborhood sizes, and
embeddings.

For a large-data audit, sample `query_rows` and reuse them for every
embedding. Each score is exact for its query against the full pool. The
column means summarize the sampled queries. Slicing `Xin` and the
embeddings before the call changes the candidate pool and can change the
scores.

## References

Venna and Kaski (2001) give the rank penalties for trustworthiness and
continuity.
[doi:10.1007/3-540-44668-0_68](https://doi.org/10.1007/3-540-44668-0_68)

Lee et al. (2013) define RNX. Lee, Peluffo-Ordóñez and Verleysen (2015)
summarize the curve with weights proportional to `1 / k`.
[doi:10.1016/j.neucom.2012.12.036](https://doi.org/10.1016/j.neucom.2012.12.036)
and
[doi:10.1016/j.neucom.2014.12.095](https://doi.org/10.1016/j.neucom.2014.12.095)

Lee and Verleysen (2008) give the two mean relative rank errors in
equations 5 and 6. Quadra reports the errors directly, with zero for
agreement. <https://proceedings.mlr.press/v4/lee08a.html> [ZADU
0.5.4](https://github.com/hj-n/zadu/tree/v0.5.4) instead reports one
minus the normalized error.

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x)$x[, 1:2]
profile <- exact_rank_profiles(x, list(PCA = y), k = c(5, 15), mrre = TRUE)
colMeans(profile$PCA$mrre_input)
#>         k5        k15 
#> 0.01600334 0.02117424 
set.seed(42)
queries <- sample.int(nrow(x), 30)
spot <- exact_rank_profiles(x, list(PCA = y), k = 5, query_rows = queries)
colMeans(spot$PCA$curve) # 30 queries, 149 candidates per query
#>      rnx5 
#> 0.6412963 
```
