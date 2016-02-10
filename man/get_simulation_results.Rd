\name{get_simulation_results}
\alias{get_simulation_results}
\title{
Retrieve series of simulated variables
}
\description{
The \code{get_simulation_results} function retrieves the series of simulated variables
from an object of the \code{gecon_simulation} class.
}
\usage{
get_simulation_results(sim_obj)
}
\arguments{
  \item{sim_obj}{
An object of the \code{\link{gecon_simulation-class}} class.
}
}
\value{
The results are returned as one element list when the simulation has been invoked by
the \code{random_path} or \code{simulate_model} functions or a list of more elements corresponding to the number
of shocks when the simulation has been performed with the \code{compute_irf} function.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
For details, see \code{\link{gecon_simulation-class}}.
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

# compute and print the IRFs
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
irf_rbc <- compute_irf(rbc, cholesky = TRUE,
                       var_list = c('K_s', 'C', 'Z', 'I', 'Y'))
get_simulation_results(irf_rbc)
}
