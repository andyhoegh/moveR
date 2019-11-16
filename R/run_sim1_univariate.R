#' A Function to fit univariate ABM for simulation 1
#'
#' @param num_mcmc The number of MCMC iterations
#' @param num_particles The number of particles to use in particle filter
#' @param data Time_points x 2 matrix of observed locations
#' @param m_mu Prior values for step size mean: # mu ~ N(m_mu, sigmasq_m)
#' @param sigmasq_m Prior values for step size mean: # mu ~ N(m_mu, sigmasq_m)
#' @param nu Prior values for step size standard deviation: sigmasq_u ~ IG(nu/2,sigmasq0 * nu / 2)
#' @param sigmasq0 Prior values for step size standard deviation: sigmasq_u ~ IG(nu/2,sigmasq0 * nu / 2)
#' @param nu_eps Prior values for observation error: # sigmasq_eps ~ IG(nu_eps/2,nu_eps * sigmasq0_eps / 2)
#' @param sigmasq0_eps Prior values for observation error: # sigmasq_eps ~ IG(nu_eps/2,nu_eps * sigmasq0_eps / 2)
#' @param mu_theta_mean Prior values for projected normal: mu_theta ~ N(mu_theta_mean, var_theta)
#' @param var_theta Prior values for projected normal: mu_theta ~ N(mu_theta_mean, var_theta)
#' @return A list containing the path and the log probability of the path
#' @export
run_sim1_univariate <- function(num_mcmc, num_particles, data, m_mu, sigmasq_m,
                                nu, sigmasq0, nu_eps, sigmasq0_eps, mu_theta_mean, var_theta){
  accept_count <- rep(0, num_mcmc)
  log_pi <- rep(0, num_mcmc)
  time_points <- nrow(data)

  # priors
  # mu ~ N(m_mu, sigmasq_m)
  a_lower <- 0
  b_upper <- Inf
  # sigmasq_u ~ IG(a,b)
  a_sigmau <- nu / 2
  b_sigmau <- nu *sigmasq0 /2
  # sigmasq_eps ~ IG(a,b)
  a_sigmasq_eps <- nu_eps / 2
  b_sigmasq_eps <- nu_eps *sigmasq0_eps /2
  # mu_theta ~ N(mu_theta_mean, var_theta)
  var_theta_inv <- solve(var_theta)
  m <- 50
  r_grid <- seq(.1,6, length.out = m)
  r <- rep(0,time_points)
  sigmasq_eta <- 0

  mu_samples <- rep(m_mu, num_mcmc)
  sigma_samples <- rep(sqrt(sigmasq0), num_mcmc)
  mu_theta_samples <- matrix(mu_theta_mean, num_mcmc,2, byrow = T)
  sigmasq_eps_samples <- rep(sigmasq0_eps, num_mcmc)
  state_variables <- array(0, dim=c(time_points, num_mcmc, 6))

  pb <- progress::progress_bar$new(
    format = " MCMC iterations [:bar] :percent complete in :elapsed",
    total = 100, clear = FALSE, width= 60)

  smc <- moveR::run_smc(num_particles,data, mu_samples[1], sigma_samples[1], mu_theta_samples[1,],
                        sigmasq_eta = sigmasq_eta, sigmasq_eps = sigmasq_eps_samples[1])

  state_variables[,1,] <- smc$path
  log_pi[1] <- smc$log_pi


  for (iter in 2:num_mcmc){
    if (iter %% (num_mcmc / 100) == 0) pb$tick()

    ## inverse CDF step
    U <- (stats::pnorm((state_variables[, iter - 1, 5] - mu_samples[iter-1]) / sigma_samples[iter-1]) - stats::pnorm((a_lower - mu_samples[iter-1]) / sigma_samples[iter-1])) /
      (stats::pnorm((b_upper - mu_samples[iter-1]) / sigma_samples[iter-1]) - stats::pnorm((a_lower - mu_samples[iter-1]) / sigma_samples[iter-1]))
    y <- mu_samples[iter-1] + sigma_samples[iter-1] * stats::qnorm(U)

    # sample mu
    cov_mu <- solve(1/sigmasq_m + time_points / sigma_samples[iter - 1]^2)
    mean_mu <- cov_mu * (m_mu / sigmasq_m + sum(y) / sigma_samples[iter - 1]^2)
    mu_samples[iter] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_mu, sd = sqrt(cov_mu))

    # sample sd
    sigma_samples[iter] <- sqrt(LearnBayes::rigamma(1,time_points / 2 + a_sigmau, b_sigmau +
                                                      .5 * sum(((y - mu_samples[iter])^2))))

    ## PMMH Step
    smc <- moveR::run_smc(num_particles, data, mu_samples[iter], sigma_samples[iter], mu_theta_samples[iter - 1, ], sigmasq_eta = sigmasq_eta, sigmasq_eps = sigmasq_eps_samples[iter-1])

    log.ratio <- (smc$log_pi) - (log_pi[iter - 1])
    if (log(stats::runif(1)) < (log.ratio)){
      # accept
      log_pi[iter] <- smc$log_pi
      state_variables[,iter,] <- smc$path
      accept_count[iter] <- 1
    } else{
      log_pi[iter] <- log_pi[iter - 1]
      state_variables[,iter,] <- state_variables[,iter - 1,]
    }

    # sample sigmasq eps
    sigmasq_eps_samples[iter] <- LearnBayes::rigamma(1, a_sigmasq_eps + time_points * 2 / 2,
                                                     b_sigmasq_eps + .5 * sum((state_variables[,iter, 1:2] - data)^2))

    # sample mu theta
    state_variables[,iter,6]
    w_theta <- matrix(r_grid, time_points,m, byrow  = T) * exp(matrix(-.5 * r_grid^2, time_points,m, byrow  = T) +
                                                                 matrix(cbind(cos(state_variables[,iter,6]), sin(state_variables[,iter,6]))%*% mu_theta_samples[iter - 1, ],time_points,1) %*% matrix(r_grid, 1, m))
    for (pts in 1:time_points){
      r[pts] <- sample(r_grid,1, prob = w_theta[pts,])
    }
    x <- r * cbind(cos(state_variables[,iter,6]), sin(state_variables[,iter,6]))

    cov_mu_theta <- base::solve(time_points * diag(2) + var_theta_inv)
    exp_mu_theta <- cov_mu_theta %*% (colSums(x) + var_theta_inv %*% mu_theta_mean)
    mu_theta_samples[iter,] <- mnormt::rmnorm(1, exp_mu_theta, cov_mu_theta)
  }
  return(list(state_variables = state_variables, accept_count = mean(accept_count), mu_theta_samples = mu_theta_samples,sigmasq_eps_samples=sigmasq_eps_samples,
              sigma_samples=sigma_samples, mu_samples=mu_samples))
}
