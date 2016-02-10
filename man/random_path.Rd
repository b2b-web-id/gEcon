\name{random_path}
\alias{random_path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the model using a random path of shocks
}
\description{
This function draws random shocks from distribution with user specified covariance matrix
and then simulates the behaviour of the system.
}
\usage{
random_path(model, shock_list = NULL, var_list = NULL, path_length = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{shock_list}{
a list of shock names that should be taken into account. If not specified, the system of all the shocks is simulated.
}
  \item{var_list}{
a list of variables on which the impact of shocks is to be computed. By default, the impact of random path is evaluated
for the state variables only.
}
  \item{path_length}{
the length of stochastic path, default value = 100.
}
}
\details{
Cf. gEcon manual, chapter "Model analysis".
}
\value{
An object of the \code{\link{gecon_simulation}} class with simulated paths of variables.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
the \code{\link{simulate_model}} function enables user to specify her own path of shocks and simulate the impact.
The function returns an object of the \code{gecon_simulation} class. Generic functions such as \code{print}
and \code{summary} provide information about the simulations.
The \code{plot_simulation} function allows to visualize the impact on variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# perform simulation and plot the results
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
irf_rbc <- random_path(rbc, path_length = 100,
                          var_list = c('K_s', 'C', 'Z', 'I', 'Y'))
plot_simulation(irf_rbc)
}
