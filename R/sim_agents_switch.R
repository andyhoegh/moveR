#' A Function to simulate animal movement with an ABM with a Markov switching process
#'
#' @param num_agents The number of agents to simulate
#' @param time_points The number of discrete time steps
#' @param mu_true Vector with mean step size in truncated normal distribution
#' @param sigma_true vector with the standard deviation of mean step size in truncated normal distribution
#' @param mu_theta_true Vector for mean of multivariate normal (for projected normal distribution)
#' @param sigmasq_eta Value for eta
#' @param sigmasq_eps Value for eps
#' @param x_range Range of x-values to simulate starting points of agents
#' @param y_range Range of y-values to simulate starting points of agents
#' @param pi1 Vector with probability of moving from state 1 to states 1 and 2
#' @param pi2 Vector with probability of moving from state 2 to states 1 and 2
#' @return An array of dimension (num_agents X time_points x 2) containing the agent locations
#' @export
sim_agents_switch <- function(num_agents, time_points, mu_true, sigma_true, mu_theta_true,
                              sigmasq_eps, sigmasq_eta, x_range, y_range, pi1, pi2){
  pi <- cbind(pi1, pi2)
  z <- s <- delta <- array(0, dim=c(num_agents, time_points, 2))

  states <- u <- theta <- matrix(0, num_agents, time_points)
  theta_tmp <- mnormt::rmnorm(num_agents, mean = mu_theta_true, varcov = diag(2))

  theta[,1] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta

  ## starting locations (T = 1)

  s[,1,] <- cbind(stats::runif(num_agents, min = x_range[1], max = x_range[2]),
                  stats::runif(num_agents, min = y_range[1], max = y_range[2]))
  z[,1,] <- s[,1,] + matrix(stats::rnorm(num_agents * 2, sd = sqrt(sigmasq_eps)), nrow = num_agents, ncol = 2)

  states[,1] <- sample(2, 10,replace = T)

  #  First Step (T = 2)
  for (agent in 1:num_agents){
    states[agent,2] <- sample(2,1, prob = pi[,states[agent,1]])
  }

  # update angle
  theta_tmp <- mnormt::rmnorm(num_agents, mean = c(0,0), varcov = diag(2))

  theta[,2] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta
  delta[,2,1] <- cos(theta[,2])
  delta[,2,2] <- sin(theta[,2])


  for (agent in 1:num_agents){
    # update speed
    u[agent,2] <- truncnorm::rtruncnorm(1, a = 0, b= Inf, mean = mu_true[states[agent,2]], sd = sigma_true[states[agent,2]])
  }

  # update latent agent locations
  s[,2,] <- s[,1,] + u[,2] * delta[,2,] + matrix(stats::rnorm(num_agents * 2, mean = 0 , sd = sqrt(sigmasq_eta)), nrow = num_agents, ncol = 2)

  # update observed locations
  z[,2,] <- s[,2,] + matrix(stats::rnorm(num_agents * 2, sd = sqrt(sigmasq_eps)), nrow = num_agents, ncol = 2)

  for (t in 3:time_points){
    for (agent in 1:num_agents){
      states[agent,t] <- sample(2,1, prob = pi[,states[agent,t-1]])
    }
    states[,t]

    # update angle
    home_path <- tibble::tibble(x =z[,1,1]  - z[,t-1, 1], y = z[,1,2]  - z[,t-1, 2])
    home_angle <-  useful::cart2pol(home_path$x, home_path$y)$theta
    last_path <- tibble::tibble(x =z[,t-1,1]  - z[,t-2, 1], y = z[,t-1,2]  - z[,t-2, 2])
    last_angle <-  useful::cart2pol(last_path$x, last_path$y)$theta

    theta_tmp <- mnormt::rmnorm(num_agents, mean = mu_theta_true, varcov = diag(2))

    theta[,t] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta
    angles <- cbind(home_angle,last_angle)
    for (agent in 1:num_agents){
      delta[agent,t,1] <- cos(theta[agent,t] + angles[agent,states[agent,t]])
      delta[agent,t,2] <- sin(theta[agent,t] + angles[agent,states[agent,t]])

      # update speed
      u[agent,t] <- truncnorm::rtruncnorm(1, a = 0, b= Inf, mean = mu_true[states[agent,t]], sd = sigma_true[states[agent,2]])
    }

    # update latent agent locations
    s[,t,] <- s[,t-1,] + u[,t] * delta[,t,] + matrix(stats::rnorm(num_agents * 2, mean = 0 , sd = sqrt(sigmasq_eta)), nrow = num_agents, ncol = 2)

    # update observed locations
    z[,t,] <- s[,t,] + matrix(stats::rnorm(num_agents * 2, sd = sqrt(sigmasq_eps)), nrow = num_agents, ncol = 2)
  }
  return(z)
}
