test_that("deterministic overlap preserves existing and absent RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  })
  idx <- matrix(c(2, 1), ncol = 1)
  calls <- list(
    function(threads) nbr_pres_knn(idx, idx, k = 1, n_threads = threads),
    function(threads) {
      nn_preservation(
        list(idx = idx),
        list(idx = idx),
        k = 1,
        n_threads = threads
      )
    }
  )
  for (fun in calls) {
    for (threads in c(0, 2)) {
      set.seed(42)
      seed <- .Random.seed
      invisible(fun(threads))
      expect_identical(.Random.seed, seed)
      rm(list = ".Random.seed", envir = .GlobalEnv)
      invisible(fun(threads))
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  }
})
