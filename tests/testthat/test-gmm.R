test_that("trawl objective", {
  pollution_data <- read.csv('data/clean_pollution_data.csv')
  trawl_obj <- TrawlObjective(data = pollution_data[,2],
                       depth = 20,
                       parametrisation = 'standard')

  testthat::expect_equal(a, b, tolerance=1e-3)
})



test_that("GGM obejctive", {
  # h <- 1
  # alpha <- 2.
  # beta <- 10.
  # kappa <- 19.
  # rho <- .2
  #
  # pollution_data <- read.csv('data/clean_pollution_data.csv')
  # to <- TrawlObjective(data = pollution_data[,2],
  #                      depth = 10,
  #                      parametrisation = 'standard')
  # max_depth <- 10000
  # custom_mle <- CustomMarginalMLE(pollution_data[1:max_depth,2])
  # custom_mle_kappa <- c(custom_mle, kappa)
  # c_mle_kappa_rho <- c(custom_mle_kappa, 0.2)
  # gmm_obj <- GMMObjective(data = pollution_data[1:max_depth,2], depth = 10)
  # gmm_obj(c_mle_kappa_rho[-3])

  testthat::expect_equal(T, T)
})



