implemented_names <- c(
  "exponential_trawl", "sum_exponential_trawl", "sup_ig_trawl", "gamma_trawl"
)
type_names <- c("exp", "sum_exp", "sup_ig", "gamma")

test_that("trawl_utils__get_trawl_env", {
  vals <- get_trawl_envs_list()
  testthat::expect_true(
    all(vapply(implemented_names, function(x) x %in% vals, T))
  )
})

test_that("trawl_utils__check_attibutes", {
  for (nm in c(type_names, "tmp")) {
    if (nm %in% type_names) {
      trwl_env <- get_trawl_env(nm)
      testthat::expect_true("trawl_b_one" %in% names(trwl_env))
      testthat::expect_true("trawl_b_two" %in% names(trwl_env))
      testthat::expect_true("trawl_b_three" %in% names(trwl_env))
    } else {
      testthat::expect_null(get_trawl_env(nm))
    }
  }
})

test_that("trawl_utils__params_cfg", {
  for (nm in type_names) {
    cfg <- get_trawl_params_config(nm)
    testthat::expect_true(is.list(cfg))
    testthat::expect_false(any(is.na(cfg)))
    testthat::expect_true("n_params" %in% names(cfg))
  }
})
