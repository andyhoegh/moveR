% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_smc.R
\name{run_smc}
\alias{run_smc}
\title{A Function to run a particle filter}
\usage{
run_smc(num_particles, data, mu_val, sigma_val, mu_theta_val, sigmasq_eta,
  sigmasq_eps)
}
\arguments{
\item{num_particles}{The number of particles for the particle filter for animal movement model}

\item{data}{A data frame with num_times rows and 2 columns with location data}

\item{mu_val}{Value for the mean step size}

\item{sigma_val}{Value for the standard deviation of step size}

\item{mu_theta_val}{Vector for 2-D Normal mean for projected Normal Distribution}

\item{sigmasq_eta}{Value for eta}

\item{sigmasq_eps}{Value for eps}
}
\value{
A list containing the path and the log probability of the path
}
\description{
A Function to run a particle filter
}
