#' A Function to fit multistate ABM for simulation 2
#'
#' @param num_mcmc The number of MCMC iterations
#' @param num_particles The number of particles to use in particle filter
#' @param data num_agents X time_points x 2 array of observed locations
#' @param m_mu1 Prior values for step size mean for state 1: # mu ~ N(m_mu, sigmasq_m)
#' @param m_mu2 Prior values for step size mean for state 2: # mu ~ N(m_mu, sigmasq_m)
#' @param sigmasq_m1 Prior values for step size variance for state 1: # mu ~ N(m_mu, sigmasq_m)
#' @param sigmasq_m2 Prior values for step size variance for state 2: # mu ~ N(m_mu, sigmasq_m)
#' @param nu Prior values for step size standard deviation: sigmasq_u ~ IG(nu/2,sigmasq0 * nu / 2)
#' @param sigmasq01 Prior values for step size standard deviation for state 1: sigmasq_u ~ IG(nu/2,sigmasq0 * nu / 2)
#' @param sigmasq02 Prior values for step size standard deviation for state 2: sigmasq_u ~ IG(nu/2,sigmasq0 * nu / 2)
#' @param nu_eps Prior values for observation error: # sigmasq_eps ~ IG(nu_eps/2,nu_eps * sigmasq0_eps / 2)
#' @param sigmasq0_eps Prior values for observation error: # sigmasq_eps ~ IG(nu_eps/2,nu_eps * sigmasq0_eps / 2)
#' @param mu_theta_mean Prior values for projected normal: mu_theta ~ N(mu_theta_mean, var_theta)
#' @param var_theta Prior values for projected normal: mu_theta ~ N(mu_theta_mean, var_theta)
#' @return A list containing the path and the log probability of the path
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @export
run_sim2 <- function(num_mcmc, num_particles, data, m_mu1, m_mu2, sigmasq_m1,sigmasq_m2,
                     nu, sigmasq01,sigmasq02, nu_eps, sigmasq0_eps, mu_theta_mean, var_theta){
  time_points <- dim(data)[2]
  num_agents <- dim(data)[1]
  log_pi <- matrix(0, num_mcmc, num_agents)
  accept_count <- matrix(0, num_mcmc,num_agents)
  # priors
  # mu ~ N(m_mu, sigmasq_m)
  a_lower <- 0
  b_upper <- Inf
  # sigmasq_u ~ IG(a,b)
  a_sigmau1 <- a_sigmau2 <- nu / 2
  b_sigmau1 <- nu *sigmasq01 / 2
  b_sigmau2 <- nu *sigmasq02 /2
  # sigmasq_eps ~ IG(a,b)
  a_sigmasq_eps <- nu_eps / 2
  b_sigmasq_eps <- nu_eps *sigmasq0_eps /2
  # mu_theta ~ N(mu_theta_mean, var_theta)
  var_theta_inv <- solve(var_theta)
  m <- 50 # number of grid points for griddy Gibbs
  r_grid <- seq(.1,6, length.out = m)
  r <- rep(0,time_points*num_agents)
  sigmasq_eta <- 0

  mu_samples <- matrix(c(m_mu1,m_mu2),byrow = T, nrow = num_mcmc, ncol = 2)
  sigma_samples <- matrix(c(sqrt(sigmasq01), sqrt(sigmasq02)),  byrow = T, nrow = num_mcmc, ncol = 2)
  mu_theta_samples <- matrix(mu_theta_mean, num_mcmc,2, byrow = T)
  sigmasq_eps_samples <- rep(sigmasq0_eps, num_mcmc)
  pi_values <- matrix(c(.9,.1,.1,.9), byrow = T, nrow = num_mcmc, ncol = 4)
  state_variables <- array(0, dim=c(time_points, num_mcmc, 7, num_agents))

  pb <- progress::progress_bar$new(
    format = " MCMC iterations [:bar] :percent complete in :elapsed",
    total = 100, clear = FALSE, width= 60)
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  smc_init <- foreach::foreach(agent = 1:num_agents) %dopar% {
    moveR::run_smc_msm(num_particles,data[agent,,], mu_samples[1,], sigma_samples[1,], mu_theta_samples[1,],
                       sigmasq_eta = sigmasq_eta, sigmasq_eps = sigmasq_eps_samples[1], pi_vals = rbind(pi_values[1,1:2], pi_values[1,3:4]))
  }

  state_variables[,1,,] <- array(unlist(purrr::map(smc_init, 'path')), dim = c(time_points,7, num_agents))
  log_pi[1,] <- unlist(purrr::map(smc_init, 'log_pi'))

  for (iter in 2:num_mcmc){
    print(iter)
    if (iter %% (num_mcmc / 100) == 0) pb$tick()

    ### pi_values
    state_changes <- state_variables[,iter-1,7,] %>% reshape2::melt() %>% dplyr::group_by(Var2) %>% dplyr::mutate(value_prev = dplyr::lag(value)) %>% dplyr::ungroup() %>% dplyr::filter(!is.na(value_prev)) %>% dplyr::group_by(value, value_prev) %>% dplyr::summarize(count = sum(value_prev, na.rm = T)) %>% dpylr::ungroup()
    pi_11 <- stats::rbeta(1, state_changes %>% dplyr::filter(value ==1, value_prev == 1) %>% dplyr::select(count) + 1, state_changes %>% dplyr::filter(value ==2, value_prev == 1) %>% dplyr::select(count) )
    pi_values[iter, 1] <- pi_11
    pi_values[iter, 2] <- 1 - pi_11
    pi_22 <- stats::rbeta(1, state_changes %>% dplyr::filter(value ==2, value_prev == 2) %>% dplyr::select(count),state_changes %>% dplyr::filter(value ==1, value_prev == 2) %>% dplyr::select(count) )
    pi_values[iter, 3] <- 1 - pi_22
    pi_values[iter, 4] <- pi_22

    ## inverse CDF step
    steps <- apply(state_variables[,iter - 1,c(5,7),], 2, rbind)
    steps1 <- steps[steps[,2] == 1,1]
    steps2 <- steps[steps[,2] == 2,1]

    U1 <- (stats::pnorm((steps1 - mu_samples[iter-1,1]) / sigma_samples[iter-1,1]) - stats::pnorm((a_lower - mu_samples[iter-1,1]) / sigma_samples[iter-1,1])) /
      (stats::pnorm((b_upper - mu_samples[iter-1,1]) / sigma_samples[iter-1,1]) - stats::pnorm((a_lower - mu_samples[iter-1,1]) / sigma_samples[iter-1,1]))
    y1 <- mu_samples[iter-1,1] + sigma_samples[iter-1,1] * stats::qnorm(U1)

    U2 <- (stats::pnorm((steps2 - mu_samples[iter-1,2]) / sigma_samples[iter-1,2]) - stats::pnorm((a_lower - mu_samples[iter-1,2]) / sigma_samples[iter-1,2])) /
      (stats::pnorm((b_upper - mu_samples[iter-1,2]) / sigma_samples[iter-1,2]) - stats::pnorm((a_lower - mu_samples[iter-1,2]) / sigma_samples[iter-1,2]))
    y2 <- mu_samples[iter-1,2] + sigma_samples[iter-1,2] * stats::qnorm(U2)

    # sample mu
    # Need to update time points number
    cov_mu1 <- solve(1/sigmasq_m1 + length(y1) / sigma_samples[iter - 1,1]^2)
    mean_mu1 <- cov_mu1 * (m_mu1 / sigmasq_m1 + sum(y1) / sigma_samples[iter - 1,1]^2)
    mu_samples[iter,1] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_mu1, sd = sqrt(cov_mu1))

    cov_mu2 <- solve(1/sigmasq_m2 + length(y2) / sigma_samples[iter - 1,2]^2)
    mean_mu2 <- cov_mu2 * (m_mu2 / sigmasq_m2 + sum(y2) / sigma_samples[iter - 1,2]^2)
    mu_samples[iter,2] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_mu2, sd = sqrt(cov_mu2))

    # sample sd
    sigma_samples[iter,1] <- sqrt(LearnBayes::rigamma(1,length(y1) / 2 + a_sigmau1, b_sigmau1 +
                                                        .5 * sum(((y1 - mu_samples[iter,1])^2))))
    sigma_samples[iter,2] <- sqrt(LearnBayes::rigamma(1,length(y2) / 2 + a_sigmau1, b_sigmau1 +
                                                        .5 * sum(((y2 - mu_samples[iter,2])^2))))

    ## PMMH Step
    smc <- foreach::foreach(agent = 1:num_agents) %dopar% {
      moveR::run_smc_msm(num_particles,data[agent,,], mu_samples[iter,], sigma_samples[iter,], mu_theta_samples[iter - 1,],
                         sigmasq_eta = sigmasq_eta, sigmasq_eps = sigmasq_eps_samples[iter - 1], pi_vals = rbind(pi_values[iter,1:2], pi_values[iter,3:4]))
    }

    state_variables_star <- array(unlist(purrr::map(smc, 'path')), dim = c(time_points,7, num_agents))
    log_pi_star <- unlist(purrr::map(smc, 'log_pi'))


    for (agent in 1:num_agents){
      log.ratio <- (log_pi_star[agent]) - (log_pi[iter - 1,agent])
      if (log(stats::runif(1)) < (log.ratio)){
        # accept
        log_pi[iter,agent] <- log_pi_star[agent]
        state_variables[,iter,,agent] <- state_variables_star[,,agent]
        accept_count[iter,agent] <- 1
      } else{
        log_pi[iter,agent] <- log_pi[iter - 1,agent]
        state_variables[,iter,,agent] <- state_variables[,iter - 1,,agent]
      }
    }


    # sample sigmasq eps
    ss <- sum((state_variables[,iter, 1,] - t(data[,,1]))^2 +
                (state_variables[,iter, 2,] - t(data[,,2]))^2)
    sigmasq_eps_samples[iter] <- LearnBayes::rigamma(1, a_sigmasq_eps + time_points * 2 * num_agents / 2,
                                                     b_sigmasq_eps + .5 * ss)

    # sample mu theta
    w_theta <- matrix(r_grid, time_points * num_agents,m, byrow  = T) * exp(matrix(-.5 * r_grid^2, time_points* num_agents,m, byrow  = T) +
                                                                              matrix(cbind(cos(as.numeric(state_variables[,iter,6,])), sin(as.numeric(state_variables[,iter,6,])))%*% mu_theta_samples[iter - 1, ],time_points* num_agents,1) %*% matrix(r_grid, 1, m))
    for (pts in 1:(time_points*num_agents)){
      r[pts] <- sample(r_grid,1, prob = w_theta[pts,])
    }
    x <- r * cbind(cos(as.numeric(state_variables[,iter,6,])), sin(as.numeric(state_variables[,iter,6,])))

    cov_mu_theta <- base::solve(time_points * num_agents * diag(2) + var_theta_inv)
    exp_mu_theta <- cov_mu_theta %*% (colSums(x) + var_theta_inv %*% mu_theta_mean)
    mu_theta_samples[iter,] <- mnormt::rmnorm(1, exp_mu_theta, cov_mu_theta)
  }
  return(list(state_variables = state_variables, accept_count = mean(accept_count), mu_theta_samples = mu_theta_samples,sigmasq_eps_samples=sigmasq_eps_samples,
              sigma_samples=sigma_samples, mu_samples=mu_samples))
}
