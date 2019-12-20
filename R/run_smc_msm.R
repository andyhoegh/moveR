#' A Function to run a particle filter with state switching behavior
#'
#' @param num_particles The number of particles for the particle filter for animal movement model
#' @param data A data frame with num_times rows and 2 columns with location data
#' @param mu_val Value for the mean step size
#' @param sigma_val Value for the standard deviation of step size
#' @param mu_theta_val Vector for 2-D Normal mean for projected Normal Distribution
#' @param sigmasq_eta Value for eta
#' @param sigmasq_eps Value for eps
#' @param pi_vals 2 X 2 matrix with markov transition probabilities
#' @return A list containing the path and the log probability of the path
#' @export
run_smc_msm <- function(num_particles, data, mu_val, sigma_val, mu_theta_val, sigmasq_eta, sigmasq_eps, pi_vals){
  samp_state <- function(pi){
    return(base::sample(1:2, 1, prob= pi))
  }
  time_points <- nrow(data)
  particle_values <- array(0, dim=c(time_points, num_particles, 7))
  w <- matrix(0, nrow = num_particles, ncol = time_points)

  # Time 1
  particle_values[1,,1:2] <- LearnBayes::rmnorm(num_particles, mean = c(data[1,1], data[1,2]), varcov = diag(2)*.1)
  theta_tmp <- mnormt::rmnorm(num_particles, mean = mu_theta_val, varcov = diag(2))
  theta <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta
  particle_values[1,,3] <- sin(theta)
  particle_values[1,,3] <- cos(theta)
  particle_values[1,,5] <- truncnorm::rtruncnorm(num_particles, a = 0, b = Inf, mean = mu_val, sd = sigma_val)
  particle_values[1,,7] <- sample(1:2, num_particles, replace = T)
  descendents <- array(0, dim=c(num_particles, time_points))

  #calculate weights
  log_w <- LearnBayes::dmnorm(particle_values[1,,1:2], mean = c(data[1,1],data[1,2]), varcov = diag(2) * sigmasq_eps, log = T) -
    LearnBayes::dmnorm(particle_values[1,,1:2], mean = c(data[1,1],data[1,2]), varcov = diag(2) * .1, log = T)
  w[,1] <- exp(log_w)
  log_w <- smcUtils::renormalize(log_w, log = T)
  descendents[,1] <- sample(num_particles, replace = T, prob = log_w)
  particle_values[1,,] <- particle_values[1,descendents[,1] ,]

  # Time 2
  # update state
  old_state <- particle_values[1,,7]
  pi_vec <- pi_vals[old_state,]
  particle_values[2,,7] <- apply(pi_vec,1,samp_state)

  # propose angles
  home_path <- tibble::tibble(x =data[1,1]  - particle_values[1, ,1],
                              y = data[1, 2] - particle_values[1,,2]  )
  home_angle <-  useful::cart2pol(home_path$x, home_path$y)$theta

  last_path <- tibble::tibble(x = particle_values[1,,1]  - data[1,1],
                              y = particle_values[1,,2]  - data[1,2])
  last_angle <-  useful::cart2pol(last_path$x, last_path$y)$theta
  angles <- cbind(home_angle,last_angle)

  particle_values[2,,6] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta

  for (particle in 1:num_particles){
    particle_values[2,particle,3] <- cos(particle_values[2,particle,6] + angles[particle,particle_values[2,particle,7]])
    particle_values[2,particle,4] <- sin(particle_values[2,particle,6] + angles[particle,particle_values[2,particle,7]])

    # update speed
    particle_values[2,particle,5]  <- truncnorm::rtruncnorm(1, a = 0, b= Inf, mean = mu_val[particle_values[2,particle,7]], sd = sigma_val[particle_values[2,particle,7]])
  }

  # update particle locations
  particle_values[2,,1:2] <- particle_values[1,,1:2] + particle_values[2,,5] * particle_values[2,,3:4] +
    stats::rnorm(num_particles * 2, mean = 0, sd = sqrt(sigmasq_eta))

  # calculate weights
  log_w <- LearnBayes::dmnorm(particle_values[2,,1:2], mean = c(data[2,1], data[2,2]), varcov = diag(2) * sigmasq_eps, log = T)
  w[,2] <- exp(log_w)
  log_w <- smcUtils::renormalize(log_w, log = T)
  descendents[,2] <- sample(num_particles, replace = T, prob = log_w)
  particle_values[2,,] <- particle_values[2, descendents[,2],]

  # Time 3:T
  for (t in 3:time_points){
    # update state
    old_state <- particle_values[t-1,,7]
    pi_vec <- pi_vals[old_state,]
    particle_values[t,,7] <- apply(pi_vec,1,samp_state)

    # propose angles
    home_path <- tibble::tibble(x =data[1,1]  - particle_values[t-1, ,1],
                                y = data[1, 2] - particle_values[t-1,,2]  )
    home_angle <-  useful::cart2pol(home_path$x, home_path$y)$theta

    last_path <- tibble::tibble(x = particle_values[t-1,,1]  - particle_values[t-2,,1],
                                y = particle_values[t-1,,2]  - particle_values[t-2,,2])
    last_angle <-  useful::cart2pol(last_path$x, last_path$y)$theta
    angles <- cbind(home_angle,last_angle)

    particle_values[t,,6] <- useful::cart2pol(theta_tmp[,1], theta_tmp[,2])$theta

    for (particle in 1:num_particles){
      particle_values[t,particle,3] <- cos(particle_values[t,particle,6] + angles[particle,particle_values[t,particle,7]])
      particle_values[t,particle,4] <- sin(particle_values[t,particle,6] + angles[particle,particle_values[t,particle,7]])

      # update speed
      particle_values[t,particle,5]  <- truncnorm::rtruncnorm(1, a = 0, b= Inf, mean = mu_val[particle_values[t,particle,7]], sd = sigma_val[particle_values[t,particle,7]])
    }

    # update particle locations
    particle_values[t,,1:2] <- particle_values[t-1,,1:2] + particle_values[t,,5] * particle_values[t,,3:4] +
      stats::rnorm(num_particles * 2, mean = 0, sd = sqrt(sigmasq_eta))

    # calculate weights
    log_w <- LearnBayes::dmnorm(particle_values[t,,1:2], mean = c(data[t,1], data[t,2]), varcov = diag(2) * sigmasq_eps, log = T)
    w[,t] <- exp(log_w)
    log_w <- smcUtils::renormalize(log_w, log = T)
    descendents[,t] <- sample(num_particles, replace = T, prob = log_w)
    particle_values[t,,] <- particle_values[t, descendents[,t],]
  }
  index <- rep(0, time_points)
  index[time_points] <- sample(num_particles, 1)
  path <- array(0, dim=c(time_points,7))
  path[time_points,] <- particle_values[time_points,index[time_points],]
  for (iter in time_points:2){
    index[iter- 1] <- descendents[index[iter], iter]
    path[iter - 1,] <- particle_values[iter - 1, index[iter-1],]
  }
  return(list(path = path, log_pi = sum(log(colMeans(w)))))
}
