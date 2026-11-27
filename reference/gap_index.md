# Gap Index for Two-Dimensional Embeddings

Measure how the relative areas of empty regions change in an embedding.
The Delaunay triangles of the embedding define the regions. Each
triangle is compared with the same three observations in the reference
coordinates.

## Usage

``` r
gap_index(Xin, Xout, n_threads = 0, ret_extra = FALSE)
```

## Arguments

- Xin:

  Reference coordinates, as a dense numeric matrix or data frame with
  observations in rows.

- Xout:

  Two-dimensional embedding coordinates for the same observations in the
  same row order. Use at least three non-collinear points with distinct
  coordinates. Nonnumeric data-frame columns are ignored in both inputs.

- n_threads:

  Maximum number of threads to use. `0` or `1` runs serially.

- ret_extra:

  Whether to return the triangles and their deformations.

## Value

The scalar Gap Index. With `ret_extra = TRUE`, a list with `gap_index`,
a matrix `triangles` of one-based row indices, and the corresponding
vectors `relative_area_in`, `relative_area_out`, and signed
`deformation`.

## Details

Each area's relative share is computed by dividing by the sum of the
areas of the same triangles in its own space. Both spaces use Euclidean
geometry. Deformation is the output share minus the input share, divided
by the larger share. Positive values indicate expansion, and negative
values indicate compression. A triangle with zero area in both spaces
has zero deformation. Positive total area is needed in each space.

The scalar is the mean absolute deformation weighted by each triangle's
larger share across the two spaces. It ranges from 0 (no change in
relative area) to 1. Relative areas are invariant to uniform rescaling
of either space. A shear can preserve them while distorting shape. The
triangles fill the convex hull, so regions at the boundary participate
in the result. The triangulation can change with the embedding or the
rows included. Requires the `geometry` package.

## Scaling

The 2D Delaunay triangulation has a number of triangles proportional to
`N`. Only their edge distances are computed in the two spaces. Use all
eligible rows when the gaps of the full map matter. Subsetting the rows
changes the triangulation and can change which gaps are visible.

## Numerical failures

A nonfinite computed distance causes an error, even with finite
coordinates. Finite distances alone do not guarantee accuracy at extreme
coordinate scales.

## References

Ros, J., Arleo, A., & Paulovich, F. (2026). Measuring Distortion in the
Empty Regions of Dimensionality Reduction Scatterplots with the Gap
Index. <https://arxiv.org/abs/2607.28324> The authors provide [reference
code](https://codeberg.org/jros/gap-index).

## See also

[`local_procrustes()`](https://jlmelville.github.io/quadra/reference/local_procrustes.md),
[`pointwise_rnx()`](https://jlmelville.github.io/quadra/reference/pointwise_rnx.md)

## Examples

``` r
x <- as.matrix(iris[, -5])
y <- prcomp(x, rank. = 2)$x
# Use one observation per distinct location in the embedding.
keep <- !duplicated(y)
gap_index(x[keep, ], y[keep, ])
#> [1] 0.4750517
```
