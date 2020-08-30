# test_that("trawl objective", {
#   pollution_data <- read.csv('../../data/clean_pollution_data.csv')
#   test_column <- 2
#   pars_gpd <- evir::gpd(pollution_data[, test_column], threshold = 0)$par.ses
#   p_plus <- mean(pollution_data[, test_column] > 0)
#   kappa <- 1/p_plus - 1.
#   pars <- c(pars_gpd, kappa)
#   trawl_obj <- TrawlObjective(
#       data = pollution_data[, test_column],
#       depth = 10,
#       parametrisation = 'standard')
#   trawl_obj_as_trawl_params <- trawl_obj(pars)
#
#   trawl_values <- vapply(1:10/10, function(x){trawl_obj_as_trawl_params(x)}, .1)
#   testthat::expect_equal(which.min(trawl_values), 2)
# })
#
#
# test_that("trawl objective - grad", {
#   pollution_data <- read.csv('../../data/clean_pollution_data.csv')
#   test_column <- 2
#   pars_gpd <- evir::gpd(pollution_data[, test_column], threshold = 0)$par.ses
#   p_plus <- mean(pollution_data[, test_column] > 0)
#   kappa <- 1/p_plus - 1.
#   pars <- c(pars_gpd, kappa)
#   trawl_obj <- TrawlObjective(
#     data = pollution_data[, test_column],
#     depth = 10,
#     parametrisation = 'standard')
#   trawl_obj_as_trawl_params <- trawl_obj(pars)
#
#   trawl_grad_values <- lapply(as.list(1:5/5*.3), function(x){pracma::grad(trawl_obj_as_trawl_params, x0 = x)})
#   trawl_grad_values <- unlist(trawl_grad_values)
#   testthat::expect_equal(which.min(abs(trawl_grad_values)), 3)
# })
#
# test_that("GMM objective - positive", {
#   max_length <- 30000
#   pollution_data <- read.csv('../../data/clean_pollution_data.csv')
#   test_column <- 2
#   init_guess_bds <- GetInitialGuessAndBounds(
#     data=pollution_data[1:max_length, test_column]
#   )
#
#   gmm_obj <- FullGMMObjective(
#     data = pollution_data[, test_column],
#     depth = 10)
#
#   testthat::expect_true(gmm_obj(c(init_guess_bds$init_guess, .2)) > 0.0)
# })

# test_that("GMM objective", {
#   max_length <- 30000
#   pollution_data <- read.csv('../../data/clean_pollution_data.csv')
#   test_column <- 2
#   init_guess_bds <- GetInitialGuessAndBounds(
#       data=pollution_data[1:max_length, test_column]
#   )
#
#   gmm_obj <- FullGMMObjective(
#       data = pollution_data[, test_column],
#       depth = 3
#   )
#
#   print(microbenchmark::microbenchmark(gmm_obj(c(init_guess_bds$init_guess, .2)), times=5))
#   print(optim(
#     par=c(init_guess_bds$init_guess, .2),
#     fn = gmm_obj, method = 'L-BFGS-B',
#     lower = c(init_guess_bds$lower, .001),
#     upper =c(init_guess_bds$upper, 2),
#     control = list(trace=3)
#   )$par)
#
#   testthat::expect_equal(T, T)
# })

# test_that("Two-stage GMM objective", {
#   max_length <- 30000
#   pollution_data <- read.csv('../../data/clean_pollution_data.csv')
#   test_column <- 2
#   init_guess_bds <- GetInitialGuessAndBounds(
#     data=pollution_data[1:max_length, test_column]
#   )
#
#   two_step_gmm_obj <- TwoStageGMMObjective(
#     data = pollution_data[, test_column],
#     depth = 10
#   )
#
#   start_time <- Sys.time()
#   trawl_param_value <- optim(
#     par=c(.2),
#     fn = two_step_gmm_obj, method = 'L-BFGS-B',
#     lower = c(.001),
#     upper =c(2),
#   )$par
#   time_delta <- Sys.time() - start_time
#
#   testthat::expect_equal(trawl_param_value, .15, tolerance=1e-2)
#   testthat::expect_true(time_delta < 30)
# })



