\name{simulate_model}
\alias{simulate_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the model
}
\description{
The \code{simulate_model} function simulates model based on realisations of shock values given by the user.
In particular it enables to compute the impact of negative shocks.
}
\usage{
simulate_model(model, shock_list = NULL, var_list = NULL, shock_m = NULL,
               periods = NULL, path_length = 40,
               sim_type = NULL, model_name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{shock_list}{
the shock names for the rows \code{shock_m} specified by the user.
The default names are the names of the first shocks from the list of shocks up to the number of \code{shock_m} matrix rows.
}
  \item{var_list}{
the list of variables for which the impact has to be computed.  By default, the impact of shocks is evaluated
for the state variables only.
}
  \item{shock_m}{
a matrix or vector of shocks given by the user. Values for different shocks should be stored in rows and values for periods in columns.
}
  \item{periods}{
the number of periods for which, shocks in the \code{shock_m} function have been specified.
The default values are from 1 to the number of columns of the shock matrix.
}
  \item{path_length}{
the number of periods for which the model is simulated. The default number is 40.
}
\item{sim_type}{
the type of simulation performed on model. It does not have to be specified when user invokes this function directly.
}
\item{model_name}{
the name of the \code{gecon_model-class} object based on which simulations are created.
The user does not have to specify the name explicitly (by default, the variable is deparsed and name is retrieved automatically).
It does not have to be specified when user invokes this function directly.
}
}
\value{
An object of the \code{\link{gecon_simulation}} class with simulated paths of variables.
}
\details{
The \code{\link{random_path}} and \code{\link{compute_irf}} functions are wrappers for this function. They generate
a path of shock(s) values and pass it on to \code{simulate_model} function, which performs computations and returns relevant results.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
the \code{\link{random_path}} function to use random path of shocks for the simulation of the model.
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

# compute and print the IRFs
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
irf_rbc <- simulate_model(rbc, shock_m = matrix(c(-0.05, -0.05), 1, 2),
                             periods=c(1, 4), var_list = c('K_s', 'C', 'Z', 'I', 'Y'))
plot_simulation(irf_rbc)
}
