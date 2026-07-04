test_that("README quick-start examples run", {
  iris_x <- as.matrix(iris[, -5])
  iris_pca2 <- stats::prcomp(iris_x, rank. = 2, scale. = FALSE, retx = TRUE)$x

  set.seed(42)
  expect_length(
    random_triplet_accuracy(iris_x, iris_pca2, n_triplets = 1),
    1
  )
  expect_length(
    random_pair_distance_correlation(iris_x, iris_pca2, n_pairs = 20),
    1
  )
  expect_length(
    random_pair_distance_correlation(
      iris_x,
      iris_pca2,
      n_pairs = 20,
      method = "spearman"
    ),
    1
  )
  expect_length(
    random_pair_distance_emd(iris_x, iris_pca2, n_pairs = 20),
    1
  )
  expect_length(
    random_pair_distance_stress(iris_x, iris_pca2, n_pairs = 20),
    1
  )

  nnp <- nn_preservation(
    iris_x,
    iris_pca2,
    k = c(5, 15),
    nn_method_in = "brute",
    nn_method_out = "brute"
  )
  expect_named(nnp, c("nnp5", "nnp15"))

  mnc <- mutual_neighbor_correlation(
    iris_x,
    iris_pca2,
    k = c(5, 15),
    nn_method_in = "brute",
    nn_method_out = "brute"
  )
  expect_named(mnc, c("mnc5", "mnc15"))

  din <- as.matrix(stats::dist(iris_x))
  dout <- as.matrix(stats::dist(iris_pca2))
  expect_length(nbr_pres(din, dout, k = 5), nrow(iris_x))
  expect_length(rnx_auc(din, dout), 1)
})
