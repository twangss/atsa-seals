#' fit_stan is the primary function which calls pre-written stan scripts for time series data.
#'
#' @param y The response variable (numeric)
#' @param x The predictors, either a vector or matrix
#' @param model_name The specific name of the model to be fitted. Currently supported are 'regression', 'ar', 'rw', 'ma', 'ss_ar' (state space univariate AR), or 'ss_rw' (state space univariate random walk).
#' @param est_drift Whether or not to estimate a drift parameter (default = FALSE). Only applicable to the rw and ar models.
#' @param est_mean Whether to estimate a mean or not (for state space autoregressive model only)
#' @param P The order of the ar model, with minimum value = 1 (default).
#' @param Q The order of the ma model, with minimum value = 1 (default).
#' @param mcmc_list A list of MCMC control parameters. These include the number of 'iterations' (default = 1000), burn in or warmup (default = 500), chains (default = 3), and thinning (default = 1)
#' @param family A named distribution for the observation model, defaults to gaussian
#' @param est_nu Boolean, whether to model process deviations as Student-t or not (default).
#' @param marss A named list containing the following elements for specifying marss models: (states=NULL, obsVariances=NULL, proVariances=NULL, trends=NULL
#' @param map_estimation Whether to do maximum a posteriori estimation via [rstan::optimizing()] (defualts to FALSE)
#' @param hessian Whether to return hessian if map_estimation is TRUE via [rstan::optimizing()]
#' @param ... Any other arguments passed to [rstan::sampling()].
#' @return an object of class 'rstan'
#' @importFrom rstan sampling
#' @export
#'
fit_stan <- function(y, x = NA,
                     mcmc_list = list(n_mcmc = 1000, n_burn = 500, n_chain = 3, n_thin = 1),
                     family = "gaussian",
                     est_nu = FALSE,
                     est_trend = FALSE,
                     # est_sigma_process_prior = FALSE,
                     marss = list(states = NULL, obsVariances = NULL, proVariances = NULL, trends = NULL,
                                  # sigma_process_prior=NULL
                                  ),
                     map_estimation = FALSE,
                     hessian = FALSE, ...) {
  dist <- c("gaussian", "binomial", "poisson", "gamma", "lognormal")
  family <- which(dist == family)

  # process potential NAs in data
  if (!is.matrix(y)) {
    N <- length(y)
    pos_indx <- which(!is.na(y))
    y <- y[pos_indx]
    n_pos <- length(pos_indx)
    # catch case where -- needs to be 2 elements for stan vec
    if (length(pos_indx) == 0) {
      pos_indx <- rep(0, 2)
    } else {
      pos_indx <- c(pos_indx, 0, 0)
    }

  }

  data <- NA
  
    object <- stanmodels$marss
    if (is.null(marss$states)) marss$states <- rep(1, nrow(y))
    if(length(marss$states) != nrow(y)) stop("Error: state vector must be same length as number of time series in y")
    if (is.null(marss$obsVariances)) marss$obsVariances <- rep(1, nrow(y))
    if(length(marss$obsVariances) != nrow(y)) stop("Error: vector of observation error variances must be same length as number of time series in y")
    if (is.null(marss$proVariances)) marss$proVariances <- rep(1, max(marss$states))
    if(length(marss$proVariances) < max(marss$states)) stop("Error: vector of process error variances is fewer than the number of states")
    if(length(marss$proVariances) > max(marss$states)) stop("Error: vector of process error variances is larger than the number of states")
    if (is.null(marss$trends)) marss$trends <- rep(1, max(marss$states))
    if(length(marss$trends) < max(marss$states)) stop("Error: vector of trends is fewer than the number of states")
    if(length(marss$trends) > max(marss$states)) stop("Error: vector of trends is larger than the number of states")

    proVariances <- c(marss$proVariances, 0) # to keep types in stan constant
    trends <- c(marss$trends, 0) # to keep types in stan constant
    N <- ncol(y)
    M <- nrow(y)
    row_indx_pos <- matrix((rep(1:M, N)), M, N)[which(!is.na(y))]
    col_indx_pos <- matrix(sort(rep(1:N, M)), M, N)[which(!is.na(y))]
    n_pos <- length(row_indx_pos)
    y <- y[which(!is.na(y))]
    
    est_A <- rep(1, M)
    for(i in 1:max(marss$states)) {
      indx <- which(marss$states==i)
      est_A[indx[1]] <- 0
    }
    est_A <- which(est_A > 0)
    est_A <- c(est_A, 0, 0)
    n_A <- length(est_A) - 2
    
    data = list("N"=N,"M"=M, "y"=y,
                     "states"=marss$states, "S" = max(marss$states), "obsVariances"=marss$obsVariances,
                     "n_obsvar" = max(marss$obsVariances), "proVariances" = proVariances,
                     "n_provar" = max(proVariances),
                     "trends"=trends, "n_trends" = max(trends),
                     "n_pos" = n_pos,
                     "col_indx_pos" = col_indx_pos,
                     "row_indx_pos" = row_indx_pos,
                     "est_A" = est_A,
                     "n_A" = n_A,
                     "est_nu" = est_nu,
                     "est_trend" = est_trend,
                     # "est_sigma_process_prior" = est_sigma_process_prior,
                     # "sigma_process_prior" = marss$sigma_process_prior,
                     "family"=1)
              
    pars = c("pred", "log_lik","sigma_process","sigma_obs","x0")
    # pars = c("pred", "sigma_process","sigma_obs","x0", "log_lik")
    #if(marss$est_B) pars = c(pars, "B")
    if(est_trend) pars = c(pars, "U")
    if(n_A > 0) pars = c(pars,"A")
    if(est_nu) pars = c(pars,"nu")
    out <- rstan::sampling(
      object = object,
      data = data,
      pars = pars,
      control = list(adapt_delta = 0.95, max_treedepth = 13),
      warmup = mcmc_list$n_burn,
      iter = mcmc_list$n_mcmc,
      thin = mcmc_list$n_thin,
      chains = mcmc_list$n_chain, ...
    )
  return(out)
  #return(list(model = out, data = data, pars = pars))
}
