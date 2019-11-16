#' A Function to simulate animal movement with an ABM
#'
#' @param num_agents The number of agents to simulate
#' @param time_points The number of discrete time steps
#' @param mu_true Value for mean step size in truncated normal distribution
#' @param sigma_true Value for the standard deviation of mean step size in truncated normal distribution
#' @param mu_theta_true Vector for mean of multivariate normal (for projected normal distribution)
#' @param sigmasq_eta Value for eta
#' @param sigmasq_eps Value for eps
#' @param x_range Range of x-values to simulate starting points of agents
#' @param y_range Range of y-values to simulate starting points of agents
#' @return An array of dimension (num_agents X time_points x 2) containing the agent locations
#' @export
sim_agents <- function(num_agents, time_points, mu_true, sigma_true, mu_theta_true,
                       sigmasq_eps, sigmasq_eta, x_range, y_range){
  z <- s <- delta <- array(0, dim=c(num_agents, time_points, 2))

  u <- theta <- matrix(0, num_agents, time_points)
  theta_tmp <- mnormt::rmnorm(num_agents, mean = mu_theta_true, varcov = diag(2))

  theta[,1] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta

  ## starting locations

  s[,1,] <- cbind(stats::runif(num_agents, min = x_range[1], max = x_range[2]),
                  stats::runif(num_agents, min = y_range[1], max = y_range[2]))
  z[,1,] <- s[,1,] + matrix(stats::rnorm(num_agents * 2, sd = sqrt(sigmasq_eps)), nrow = num_agents, ncol = 2)

  for (t in 2:time_points){
    # update angle
    home_path <- tibble::tibble(x =z[,1,1]  - z[,t-1, 1], y = z[,1,2]  - z[,t-1, 2])
    home_angle <-  useful::cart2pol(home_path$x, home_path$y)$theta
    theta_tmp <- mnormt::rmnorm(num_agents, mean = mu_theta_true, varcov = diag(2))

    theta[,t] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta
    delta[,t,1] <- cos(theta[,t] + home_angle)
    delta[,t,2] <- sin(theta[,t] + home_angle)

    # update speed
    u[,t] <- truncnorm::rtruncnorm(num_agents, a = 0, b= Inf, mean = mu_true, sd = sigma_true)

    # update latent agent locations
    s[,t,] <- s[,t-1,] + u[,t] * delta[,t,] + matrix(stats::rnorm(num_agents * 2, mean = 0 , sd = sqrt(sigmasq_eta)), nrow = num_agents, ncol = 2)

    # update observed locations
    z[,t,] <- s[,t,] + matrix(stats::rnorm(num_agents * 2, sd = sqrt(sigmasq_eps)), nrow = num_agents, ncol = 2)
  }
  return(z)
}
