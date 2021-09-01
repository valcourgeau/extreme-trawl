
#' Fitting the LTME or EV Trawl model.
#' @param data Data to fit.
#' @param depth Depth to fit.
#' @param method Selecting `"GMM"` or `"PL"`.
#' @param mode Either one-stage/full (`"full"`) or two-stage (`"two-stage"`).
#' @param type Trawl type (e.g. `"exp"`, `"sum_exp"`, etc.)
#' @param bounds Parameters bounds. `"config"` uses trawl cfg, `"multiplier"`
#'     uses scaled initial guess.
#' @param cl Parallel cluster, defaults to `NULL`.
#' @return Vector of fitted parameters.
#' @examples
#' n <- 2000
#' test_column <- 2
#' data <- pollution_data[seq_len(n), test_column]
#' depth <- 5
#' ev_trawl_fit(data, depth, method = "GMM")
#' @importFrom stats runif
#' @export
ev_trawl_fit <- function(data, depth, method, mode = "two-stage",
                         type = "exp", bounds = "config", cl = NULL) {
  stopifnot(mode %in% c("two-stage", "full"))
  stopifnot(method %in% c("GMM", "PL"))
  stopifnot(bounds %in% c("config", "multiplier"))

  # method 'PL' or 'GMM' either 'full' or 'two-stage' modes
  trawl_cfg <- get_trawl_params_config(type)
  init_guess_lower_upper <- get_initial_guess_and_bounds(data)
  init_guess_model <- init_guess_lower_upper$init_guess
  lower_model <- init_guess_lower_upper$lower
  upper_model <- init_guess_lower_upper$upper

  init_trawl <- vapply(1:trawl_cfg$n_params, function(i) {
    stats::runif(n = 1, min = trawl_cfg$lower[i], max = trawl_cfg$upper[i])
  }, 1.)

  if (mode == "two-stage") {
    marginal_params <- composite_marginal_mle(data)
    # choose function
    if (method == "PL") {
      optim_fn <- pl_two_stage_trawl(
        data = data, depth = depth,
        type = type, cl = cl
      )
    } else {
      # method GMM
      optim_fn <- trawl_gmm$two_stage_gmm_objective(
        data = data, depth = depth, type = type
      )
    }

    # choose bounds
    lower <- trawl_cfg$lower
    upper <- trawl_cfg$upper
    if (bounds == "multiplier") {
      init_trawl <- pl_init_guess(
        data = data, depth = depth, n_trials = 40, type = type
      )
      lower <- init_trawl * 0.5
      upper <- init_trawl * 1.5
    }

    init_guess <- init_trawl
  } else {
    # mode full below
    # choose function
    if (method == "PL") {
      optim_fn <- pl_trawl(data = data, depth = depth, type = type, cl = cl)
    } else {
      # Method GMM
      optim_fn <- trawl_gmm$full_gmm_objective(
        data = data, depth = depth, type = type
      )
    }

    # choose bounds
    lower <- trawl_cfg$lower
    upper <- trawl_cfg$upper
    if (bounds == "multiplier") {
      init_trawl <- pl_init_guess(
        data = data, depth = depth, n_trials = 40, type = type
      )
      lower <- init_trawl * 0.5
      upper <- init_trawl * 1.5
    }
    lower <- c(lower_model, lower)
    upper <- c(upper_model, upper)

    # start values
    init_guess <- c(init_guess_model, init_trawl)
  }

  trawl_inference <- stats::optim(
    fn = optim_fn, par = init_guess,
    lower = lower, upper = upper, method = "L-BFGS-B"
  )

  if (mode == "two-stage") {
    return(c(marginal_params, trawl_inference$par))
  } else {
    return(trawl_inference$par)
  }
}

#' Fitting the LTME or EV Trawl model on block subsamples.
#' @param data Data to fit.
#' @param sample_length Length of block subsamples.
#' @param depth Depth to fit.
#' @param method Selecting `"GMM"` or `"PL"`.
#' @param mode Either one-stage/full (`"full"`) or two-stage (`"two-stage"`).
#' @param type Trawl type (e.g. `"exp"`, `"sum_exp"`, etc.)
#' @param bounds Parameters bounds. `"config"` uses trawl cfg, `"multiplier"`
#'     uses scaled initial guess.
#' @param trials Number of blocks on which to perform the inference.
#' @param parallel Parallel or not, defaults to `FALSE`.
#' @param seed Seed for block selection, not used by default `NULL`.
#' @return Vector of fitted parameters.
#' @examples
#' n <- 1000
#' test_column <- 2
#' data <- pollution_data[seq_len(n), test_column]
#' depth <- 5
#' ev_trawl_fit(data, depth, method = "GMM")
#' @importFrom stats runif
#' @export
sub_sample_fit <- function(data, sample_length, depth,
                           method, mode, type, bounds,
                           trials, parallel = F, seed = NULL) {
  # method 'PL' or 'GMM'
  # depth for PL is the length of blocks, GMM depth is the ACF depth
  n <- length(data)
  stopifnot(n >= sample_length)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  start_points <- sample(1:(n - sample_length), size = trials, replace = F)

  results_list <- list()
  results <- matrix(
    0,
    nrow = trials, ncol = 3 + get_trawl_params_config(type)$n_params
  )

  if (parallel) {
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(max(cores - 1, 1))
    parallel::clusterExport(
      cl, c(
        "transformation_map_inverse", "transformation_map",
        "transformation_jacobian", "parametrisation_translator",
        "cpp_acf_trawl", "ev_trawl_fit", "cpp_case_separator"
      )
    )

    # TODO check to include sample_length, etc in clusterExport when testing

    sub_sample_time <- Sys.time()
    if (method == "GMM") {
      # we perform parallel estimations
      results <- parallel::parLapply(
        X = start_points,
        cl = cl,
        fun = function(start_pt) {
          res <- ev_trawl_fit(
            data = data[start_pt:(start_pt + sample_length)], depth = depth,
            mode = mode, method = method, type = type, bounds = bounds,
            cl = NULL
          )
          return(res)
        }
      )
    } else {
      # PL Method
      # we perform PL computation in parallel
      results <- lapply(
        X = start_points,
        FUN = function(start_pt) {
          res <- ev_trawl_fit(
            data = data[start_pt:(start_pt + sample_length)],
            depth = depth, mode = mode, method = method,
            type = type, bounds = bounds, cl = cl
          )
          return(res)
        }
      )
    }

    print(Sys.time() - sub_sample_time)
    parallel::stopCluster(cl)
    results <- matrix(unlist(results), ncol = length(results[[1]]), byrow = T)
  } else {
    print("No parallel trials.")
    sub_sample_time <- Sys.time()
    results <- t(vapply(start_points,
      FUN = function(start_pt) {
        ev_trawl_fit(
          data = data[start_pt:(start_pt + sample_length)], depth = depth,
          type = type, mode = mode, method = method, bounds = bounds, cl = NULL
        )
      },
      FUN.VALUE = rep(0, ncol(results))
    ))
    print(Sys.time() - sub_sample_time)
  }

  results_list$estimators <- results
  results_list$sample_length <- sample_length
  results_list$depth <- depth
  results_list$mode <- mode
  results_list$method <- method
  results_list$bounds <- bounds
  results_list$type <- type
  results_list$data_length <- n
  results_list$trials <- trials
  results_list$start_indices <- start_points

  return(results_list)
}


#' Fitting the LTME or EV Trawl model.
#' @param data Data to fit.
#' @param depth Depth to fit.
#' @param method Selecting `"GMM"` or `"PL"`.
#' @param mode Either one-stage/full (`"full"`) or two-stage (`"two-stage"`).
#' @param type Trawl type (e.g. `"exp"`, `"sum_exp"`, etc.)
#' @param bounds Parameters bounds. `"config"` uses trawl `cfg`, `"multiplier"`
#'     uses scaled initial guess.
#' @param parallel Boolean, defaults to `FALSE`.
#' @return List of fitted parameters, hyperparameters, etc.
#' @examples
#' n <- 2000
#' test_column <- 2
#' data <- pollution_data[seq_len(n), test_column]
#' depth <- 5
#' ev_trawl_fit(data, depth, method = "GMM")
#' @importFrom stats runif
#' @export
fit <- function(data, depth, method, mode, type, bounds, parallel = F) {
  return(sub_sample_fit(
    data = data, depth = depth, sample_length = length(data) - 1,
    method = method, mode = mode, type = type, bounds = bounds,
    parallel = parallel, trials = 1
  ))
}
