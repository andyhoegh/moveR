% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_agents_switch.R
\name{sim_agents_switch}
\alias{sim_agents_switch}
\title{A Function to simulate animal movement with an ABM with a Markov switching process}
\usage{
sim_agents_switch(num_agents, time_points, mu_true, sigma_true,
  mu_theta_true, sigmasq_eps, sigmasq_eta, x_range, y_range, pi1, pi2)
}
\arguments{
\item{num_agents}{The number of agents to simulate}

\item{time_points}{The number of discrete time steps}

\item{mu_true}{Vector with mean step size in truncated normal distribution}

\item{sigma_true}{vector with the standard deviation of mean step size in truncated normal distribution}

\item{mu_theta_true}{Vector for mean of multivariate normal (for projected normal distribution)}

\item{sigmasq_eps}{Value for eps}

\item{sigmasq_eta}{Value for eta}

\item{x_range}{Range of x-values to simulate starting points of agents}

\item{y_range}{Range of y-values to simulate starting points of agents}

\item{pi1}{Vector with probability of moving from state 1 to states 1 and 2}

\item{pi2}{Vector with probability of moving from state 2 to states 1 and 2}
}
\value{
list with z as an array of dimension (num_agents X time_points x 2) containing the agent locations and state a vector that denotes state
}
\description{
A Function to simulate animal movement with an ABM with a Markov switching process
}
