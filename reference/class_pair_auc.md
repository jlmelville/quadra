# Retrieval AUC Within Each Pair of Classes

Restrict both queries and candidates to two classes at a time. Each
defined query receives equal weight within its class pair.

## Usage

``` r
class_pair_auc(dm, labels, metrics = c("roc", "pr"))
```

## Arguments

- dm:

  Finite square numeric distance matrix. The query itself is excluded.

- labels:

  A label vector without missing values, with one label per row of `dm`.

- metrics:

  Which AUCs to report. Defaults to both `"roc"` and `"pr"`.

## Value

A data frame with one row per unordered class pair and columns `label_a`
and `label_b` in the same first-appearance order as
[`label_trustworthiness_continuity()`](https://jlmelville.github.io/quadra/reference/label_trustworthiness_continuity.md).
Each requested metric adds its mean (`roc` or `pr`) and defined-query
count (`n_defined_roc` or `n_defined_pr`). Queries from singleton
classes are undefined. A pair of singletons has count zero and mean
`NA_real_`.

## Note

Quadra applies
[`roc_auc()`](https://jlmelville.github.io/quadra/reference/roc_auc.md)
and
[`pr_auc()`](https://jlmelville.github.io/quadra/reference/pr_auc.md) to
each two-class subset. Both use PRROC. The PR score uses Davis and
Goadrich's area.

## Scaling

The input is a dense `N` by `N` distance matrix. Each class pair also
requires a square submatrix for its members. For a particular pair,
select its two classes before computing distances. Keeping all members
of those classes gives the same pair score as extracting it from the
full analysis.

## References

Davis and Goadrich (2006).
[doi:10.1145/1143844.1143874](https://doi.org/10.1145/1143844.1143874)

Grau, Grosse and Keilwagen (2015). PRROC.
[doi:10.1093/bioinformatics/btv153](https://doi.org/10.1093/bioinformatics/btv153)

## Examples

``` r
class_pair_auc(as.matrix(dist(iris[, -5])), iris$Species)
#>      label_a    label_b       roc n_defined_roc        pr n_defined_pr
#> 1     setosa versicolor 0.9937673           100 0.9957411          100
#> 2     setosa  virginica 0.9992245           100 0.9995351          100
#> 3 versicolor  virginica 0.8031510           100 0.8217076          100
```
