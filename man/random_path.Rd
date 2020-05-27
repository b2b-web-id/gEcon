\name{random_path}
\alias{random_path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the model with a random path of shocks
}
\description{
The function generates random shock paths based on the shock covariance matrix
specified by the user and simulates the behaviour of the system.
}
\usage{
random_path(model, variables = NULL, sim_length = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or indices of variables whose paths are to be simulated. By default all variables are selected.
}
  \item{sim_length}{
the length of simulation path, the default value is 40.
}
}
\details{
Cf. gEcon users' guide, chapter `Model analysis'.
}
\value{
An object of \code{\link{gecon_simulation}} class storing simulated paths of the model's variables.
}
\seealso{
The \code{\link{simulate_model}} function allows the user to specify her own shock paths and simulate
their impact on the model's variables. This function also returns an object of \code{gecon_simulation} class.
Generic functions such as \code{print} and \code{summary} provide information about simulations' results.
The \code{plot_simulation} function allows to visualize the impact of simulated shock paths on the model's variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# perform simulation and plot the results
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc_rndpath <- random_path(rbc, variables = c("K_s", "C", "Z", "I", "Y"),
                           sim_length = 100)
plot_simulation(rbc_rndpath)
}
