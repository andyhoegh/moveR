#' A Function to run a multivariate particle filter
#'
#' @param num_particles The number of particles for the particle filter for animal movement model
#' @param data A array with dimension num_bears x num_times X 2 columns with location data
#' @param mu_val Value for the mean step size
#' @param sigma_val Value for the standard deviation of step size
#' @param kappa_val Value for kappa concentration parameter in VonMises distribution
#' @param sigmasq_eta Value for eta
#' @param sigmasq_eps Value for eps
#' @return A list containing the path and the log probability of the path
#' @export
run_mvsmc <- function(num_particles, data, mu_val, sigma_val, kappa_val, sigmasq_eta, sigmasq_eps){
  num_bears <- dim(data)[1]
  time_points <- dim(data)[2]
  particle_values <- array(0, dim=c(time_points, num_particles, 6, num_bears))
  w <- matrix(0, nrow = num_particles, ncol = time_points)

  # Time 1
  for (i in 1:num_bears){
    particle_values[1,,1:2,i] <- LearnBayes::rmnorm(num_particles, mean = c(data[i,1,1], data[i,1,2]), varcov = diag(2)*.2)
  }

  theta <- circular::rvonmises(num_particles * num_bears, mu = circular::circular(0), kappa = .1)
  particle_values[1,,3,] <- sin(theta)
  particle_values[1,,3,] <- cos(theta)
  particle_values[1,,5,] <- truncnorm::rtruncnorm(num_particles, a = 0, b = Inf, mean = mu_val, sd = sigma_val)
  descendents <- array(0, dim=c(num_particles, time_points))

  #calculate weights
  log_w <- LearnBayes::dmnorm(matrix(particle_values[1,,1:2,], nrow = num_particles, ncol = 2 * num_bears),
                              mean = c(t(data[,1,])), varcov = diag(num_bears * 2) * sigmasq_eps, log = T) -
    LearnBayes::dmnorm(matrix(particle_values[1,,1:2,], nrow = num_particles, ncol = 2 * num_bears),
                       mean = c(t(data[,1,])), varcov = diag(num_bears * 2) * .2, log = T)
  w[,1] <- exp(log_w)
  log_w <- smcUtils::renormalize(log_w, log = T)
  descendents[,1] <- sample(num_particles, replace = T, prob = log_w)
  particle_values[1,,,] <- particle_values[1,descendents[,1] ,,]

  # Time 2:T
  for (t in 2:time_points){
    # propose angles
    x_vals <-
      home_path <- tibble::tibble(x = rep(data[,1,1], each = num_particles) -
                                    as.numeric(particle_values[t-1,,1,])  ,
                                  y = rep(data[,1, 2], each = num_particles)
                                  - as.numeric(particle_values[t-1,,2,]))
    home_angle <-  useful::cart2pol(home_path$x, home_path$y)$theta
    particle_values[t,,6,] <- circular::rvonmises(num_particles * num_bears, mu = circular::circular(0), kappa = kappa_val)
    particle_values[t,,3,] <- cos(particle_values[t,,6,] + home_angle) # x coord
    particle_values[t,,4,] <- sin(particle_values[t,,6,] + home_angle) # y coord

    # propose distance
    particle_values[t,,5,] <- truncnorm::rtruncnorm(num_particles * num_bears, a = 0, b= Inf, mean = mu_val, sd = sigma_val)
    u_array <- array(0, dim = c(num_particles, 2, num_bears) )
    u_array[,1,] <- u_array[,2,] <- particle_values[t,,5,]

    # update particle locations
    particle_values[t,,1:2,] <- particle_values[t-1,,1:2,] + u_array * particle_values[t,,3:4,] +
      array(stats::rnorm(num_particles * 2 * num_bears, mean = 0, sd = sqrt(sigmasq_eta)), dim = c(num_particles, 2, num_bears))

    # calculate weights
    log_w <- LearnBayes::dmnorm(matrix(particle_values[t,,1:2,], nrow = num_particles, ncol = 2 * num_bears),
                                mean = c(t(data[,t,])), varcov = diag(num_bears * 2) * sigmasq_eps, log = T)
    w[,t] <- exp(log_w)
    log_w <- smcUtils::renormalize(log_w, log = T)
    descendents[,t] <- sample(num_particles, replace = T, prob = log_w)
    particle_values[t,,,] <- particle_values[t, descendents[,t],,]
  }
  index <- rep(0, time_points)
  index[time_points] <- sample(num_particles, 1)
  path <- array(0, dim=c(time_points,6, num_bears))
  path[time_points,,] <- particle_values[time_points,index[time_points],,]
  for (iter in time_points:2){
    index[iter- 1] <- descendents[index[iter], iter]
    path[iter - 1,,] <- particle_values[iter - 1, index[iter-1],,]
  }
  return(list(path = path, log_pi = sum(log(colMeans(w)))))
}
