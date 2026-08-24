# Neighborhood Preservation Between Distance Matrices

Returns the overlap between the `k`-neighborhoods defined by two
distance matrices, separately for each observation.

## Usage

``` r
nbr_pres(din, dout, k)
```

## Arguments

- din:

  Reference distance matrix.

- dout:

  Distance matrix to compare with `din`.

- k:

  Neighborhood size.

## Value

Per-observation neighborhood overlap in `[0, 1]`.

## Details

Diagonal entries are excluded; off-diagonal entries must be finite. Ties
at the `k`th distance are included, with each score capped at 1.
