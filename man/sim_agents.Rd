% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_agents.R
\name{sim_agents}
\alias{sim_agents}
\title{A Function to simulate animal movement with an ABM}
\usage{
sim_agents(num_agents, time_points, mu_true, sigma_true, mu_theta_true,
  sigmasq_eps, sigmasq_eta, x_range, y_range)
}
\arguments{
\item{num_agents}{The number of agents to simulate}

\item{time_points}{The number of discrete time steps}

\item{mu_true}{Value for mean step size in truncated normal distribution}

\item{sigma_true}{Value for the standard deviation of mean step size in truncated normal distribution}

\item{mu_theta_true}{Vector for mean of multivariate normal (for projected normal distribution)}

\item{sigmasq_eps}{Value for eps}

\item{sigmasq_eta}{Value for eta}

\item{x_range}{Range of x-values to simulate starting points of agents}

\item{y_range}{Range of y-values to simulate starting points of agents}
}
\value{
An array of dimension (num_agents X time_points x 2) containing the agent locations
}
\description{
A Function to simulate animal movement with an ABM
}
