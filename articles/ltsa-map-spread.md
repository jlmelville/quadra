# When two digits set the scale of an LTSA map

Two of the 5,000 handwritten digits below account for 99.72% of the sum
of squared distances from the center of the LTSA map. The other 4,998
are there if you zoom in. How much would we have noticed from the
fraction of pixel neighbors the map retained?

## What changed between the two maps?

Local Tangent Space Alignment (LTSA) fits a plane to each patch of
neighboring images and finds global coordinates that align those local
fits. We can keep the patches and the alignment matrix fixed and change
one choice. Ordinary LTSA gives each image equal weight when it
normalizes the coordinate columns. Jacobi normalization takes the
weights from the diagonal of the alignment matrix. Both maps below use
the same 15 images in each patch, including its center. [Zhang and Zha’s
LTSA paper](https://arxiv.org/abs/cs/0212008) · [Jacobi normalization in
flotsam](https://jlmelville.github.io/flotsam/articles/ltsa-normalization.html)

[![Ordinary LTSA has a tight clump and two isolated points. The
normalized map has several arms. Two images account for nearly all the
ordinary map's spread, while the normalized map takes 1,848 images to
reach 90
percent.](figures/ltsa-mnist-spread.png)](https://jlmelville.github.io/quadra/articles/figures/ltsa-mnist-spread.png)

The rings mark the same 1 and 2 in both maps. Each panel shows every
digit at its own map’s full range. To compare concentration
independently of the axis scales, measure each image’s squared distance
from its map’s centroid, then add up those distances starting with the
farthest image. In ordinary LTSA the first two images take us almost to
100%. With Jacobi normalization the top two account for 0.28%, and it
takes 1,848 images to reach 90%.

Fetch the digits, fit LTSA and draw the comparison

This is the balanced training sample from the [MNIST
workflow](https://jlmelville.github.io/quadra/articles/mnist-pca-umap.html#load-the-digits).
The pixel graph is exact. We build it once with 30 entries per row so we
can also try a larger patch. Both 15-image fits had a connected
neighborhood graph and the eigensolver reported `ok`.

``` r

# pak::pak(c("jlmelville/quadra", "jlmelville/snedata", "jlmelville/flotsam",
#            "jlmelville/vizier", "uwot"))
library(quadra)
digits <- snedata::download_mnist(as = "list")
train <- which(digits$meta$split == "training")
set.seed(413)
id <- sort(unlist(
  lapply(split(train, digits$meta$label[train]), sample, size = 500),
  use.names = FALSE
))
x <- digits$data[id, , drop = FALSE] / 255
label <- factor(digits$meta$label[id], levels = 0:9)
rm(digits)
graph <- rnndescent::brute_force_knn(x, k = 30, n_threads = 2)$idx
fit_ltsa <- function(k, normalize = FALSE) {
  set.seed(439)
  fit <- flotsam::ltsa(
    x,
    n_neighbors = k,
    nn_method = graph[, seq_len(k)],
    normalize = normalize,
    n_threads = 1,
    output = "result"
  )
  stopifnot(fit$eigen$status == "ok", fit$assembly$component_count == 1)
  fit$embedding
}
fits <- list(Ordinary = fit_ltsa(15), Jacobi = fit_ltsa(15, TRUE))
spread <- function(y) {
  d2 <- rowSums(scale(y, scale = FALSE)^2)
  d2 / sum(d2)
}
shares <- lapply(fits, spread)
extreme <- head(order(shares$Ordinary, decreasing = TRUE), 2)
data.frame(training_row = id[extreme], digit = label[extreme])
```

``` r

par(mfrow = c(1, 3), mar = c(1, 0.4, 2.2, 0.4))
for (name in names(fits)) {
  y <- fits[[name]]
  vizier::embed_plot(
    y,
    colors = grDevices::hcl.colors(10, "Set 2")[label],
    title = name,
    fixed_aspect = TRUE,
    show_axes = FALSE,
    cex = 0.45
  )
  points(y[extreme, ], pch = 1, cex = 1.4, lwd = 1.5)
}
par(mar = c(3.5, 3.6, 2.2, 0.7), mgp = c(2, 0.6, 0))
curves <- vapply(
  shares,
  function(v) {
    cumsum(sort(v, decreasing = TRUE))
  },
  numeric(nrow(x))
)
cols <- c("#17638A", "#AD4E25")
matplot(
  seq_len(nrow(x)),
  100 * curves,
  type = "l",
  log = "x",
  lty = 1:2,
  lwd = 2,
  col = cols,
  xlab = "Images, farthest first",
  ylab = "Cumulative spread (%)",
  main = "Who sets the scale?"
)
legend("left", names(fits), lty = 1:2, col = cols, lwd = 2, bty = "n")
```

## Can the neighbor score see the difference?

Look at the overlap column first. Both 15-image fits keep a little over
one in ten of the nearest pixel neighbors. Now look at fitted-scale
stress. It measures the error in the map distances after fitting the
single multiplier that best matches them to the pixel distances. Zero
would be a perfect fit. Even after rescaling, the ordinary map has
almost as much error as setting every map distance to zero.

| Fit | Spread in the two farthest images | Overlap at 15 | Fitted-scale stress |
|----|---:|---:|---:|
| Ordinary, 15 per patch | 99.72% | 0.112 | 0.999 |
| Jacobi, 15 per patch | 0.28% | 0.119 | 0.622 |
| Ordinary, 30 per patch | 0.84% | 0.118 | 0.617 |

The last row changes only the size of the patches used to fit ordinary
LTSA. The two-point concentration all but disappears, even though we
still score overlap at 15 and get much the same answer. Neighbors can
keep their order in a tiny clump. Their distances tell us how much of
the original geometry still fits inside it.

Do the two farthest digits also account for most of the error? Only 21
of the 20,000 sampled pairs contain either digit. Together those pairs
account for 0.025% of the squared error after fitting one scale to all
20,000 pairs. Almost all the error comes from distances between images
in the clump.

Measure spread, overlap and stress on the same queries and pairs

For overlap, use 100 queries per digit and all 4,999 other images as
candidates. Divide the number of retained neighbors by 15 to get the
fraction shared. For stress, use the same 20,000 pairs in every fit and
include all their distances. [Chen and Buja on neighbor
overlap](https://doi.org/10.1198/jasa.2009.0111) · [Lee and colleagues
on RNX](https://doi.org/10.1016/j.neucom.2012.12.036) · [Smelser, Miller
and Kobourov on scale](https://arxiv.org/abs/2408.07724) · [Takaeda’s
distance-band
benchmark](https://github.com/toorpia/dr-fidelity-benchmark)

``` r

fits30 <- c(fits, list(Ordinary30 = fit_ltsa(30)))
set.seed(417)
audit <- sort(unlist(lapply(
  split(seq_along(label), label),
  sample,
  size = 100
)))
set.seed(415)
pairs <- sample_pairs(nrow(x), 20000)
rank15 <- exact_rank_profiles(x, fits30, k = 15, query_rows = audit)
stress <- function(ref, y, pairs) {
  random_pair_distance_band_profile(
    ref,
    y,
    bands = 1,
    pairs = pairs
  )$normalized_stress
}
scores <- lapply(names(fits30), function(name) {
  y <- fits30[[name]]
  share <- sort(spread(y), decreasing = TRUE)
  c(
    top_two_pct = 100 * sum(head(share, 2)),
    overlap = mean(rank15[[name]]$retained) / 15,
    stress = stress(x, y, pairs)
  )
})
round(do.call(rbind, setNames(scores, names(fits30))), 3)
```

``` r

detail <- random_pair_distance_band_profile(
  x,
  fits$Ordinary,
  bands = 1,
  pairs = pairs,
  ret_extra = TRUE
)
squared_error <- (detail$distance_in -
  detail$fitted_scale * detail$distance_out)^2
touches <- pairs[, 1] %in% extreme | pairs[, 2] %in% extreme
c(
  pairs = sum(touches),
  error_pct = 100 * sum(squared_error[touches]) / sum(squared_error)
)
```

## Which distances should a Swiss roll preserve?

MNIST gave us one reference, Euclidean distance between pixel vectors.
The Swiss roll gives us two. Imagine two points directly across a fold
from each other. The straight line between them can be short even when
the route along the sheet is long. Fit ordinary LTSA and UMAP to the
same 612 points in three dimensions, then score each map against both
sets of distances. Color follows the length of the sheet from the inner
turn (purple) to the outer turn (yellow).

[![The flat Swiss roll is a strip colored from purple to yellow along
its length. LTSA bunches together parts of the strip. UMAP follows the
color progression along an elongated
sheet.](figures/ltsa-roll.png)](https://jlmelville.github.io/quadra/articles/figures/ltsa-roll.png)

Make the roll, fit both maps and measure stress on every pair

The [Swiss-roll
workflow](https://jlmelville.github.io/quadra/articles/reference-geometry.html#make-the-roll-and-its-flat-reference)
gives the arc-length formula. Keeping arc length and height in the same
units gives us distances on the unrolled sheet. Both fits see only the
rolled coordinates.

``` r

set.seed(63)
angle <- seq(1.5 * pi, 4.5 * pi, length.out = 2001)
arc <- (angle * sqrt(1 + angle^2) + asinh(angle)) / 2
grid <- expand.grid(
  s = seq(min(arc) + 1.5, max(arc) - 1.5, length.out = 34),
  h = seq(1, 29, length.out = 18)
)
n <- nrow(grid)
s <- grid$s + runif(n, -0.65, 0.65)
h <- grid$h + runif(n, -0.55, 0.55)
t <- approx(arc, angle, xout = s)$y
rolled <- cbind(t * cos(t), h, t * sin(t))
flat <- cbind(s, h)
ltsa <- flotsam::ltsa(
  rolled,
  n_neighbors = 15,
  nn_method = "exact",
  eig_method = "eig"
)
set.seed(414)
umap <- uwot::umap(
  rolled,
  n_neighbors = 15,
  min_dist = 0.1,
  n_threads = 1,
  n_sgd_threads = 1
)
roll_maps <- list(LTSA = ltsa, UMAP = umap)
all_pairs <- t(combn(n, 2))
round(
  sapply(list(Rolled = rolled, Unrolled = flat), function(ref) {
    vapply(
      roll_maps,
      function(y) {
        stress(ref, y, all_pairs)
      },
      numeric(1)
    )
  }),
  3
)
```

``` r

par(mfrow = c(1, 3), mar = c(0.4, 0.4, 2.3, 0.4))
views <- c(list("Unrolled sheet" = flat), roll_maps)
for (name in names(views)) {
  vizier::embed_plot(
    views[[name]],
    x = s,
    color_scheme = "grDevices::Viridis",
    title = name,
    fixed_aspect = TRUE,
    show_axes = FALSE
  )
}
```

| Reference for fitted-scale stress |  LTSA |  UMAP |
|-----------------------------------|------:|------:|
| Straight through the roll (3D)    | 0.335 | 0.530 |
| Along the unrolled sheet (2D)     | 0.453 | 0.064 |

If we want to read straight-line distances through the roll, LTSA has
the lower stress. If we want a map of the sheet we would get by
unrolling it, UMAP has the lower stress by a much wider margin. We
changed neither the maps nor the pairs. We changed which distance we
wanted the map to reproduce.
