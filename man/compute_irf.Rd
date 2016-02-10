\name{compute_irf}
\alias{compute_irf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute impulse response functions (IRFs)
}
\description{
The \code{compute_irf} function computes the impulse response functions for the specified set of variables and shocks and returns an object of the \code{gecon_simulation} class.
}
\usage{
compute_irf(model, shock_list = NULL, var_list = NULL,
            path_length = 40, cholesky = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{shock_list}{
a list of shocks for which the IRFs are to be computed. If this argument is missing, the IRFs are computed for all the shocks in the model.
By default, the impulse response functions are created for all the shocks in the model.
}
  \item{var_list}{
a list of variables, for which the impact of shocks has to be computed. By default, the impulse response functions are created
for the state variables only.
}
  \item{path_length}{
the number of periods for which the IRFs are to be computed.
}
  \item{cholesky}{
logical. If this option is set to FALSE, the function computes the IRFs based on uncorrelated shocks, otherwise the variance-covariance matrix is orthogonalized by using the Cholesky decomposition and the IRFs are computed by using this matrix.
}
}
\details{
Cf. gEcon manual, chapter "Model analysis".
}
\value{
The function returns an object of \code{gecon_simulation} class. Generic functions such as \code{print} and \code{summary} provide information about the impulse response functions. The \code{plot_simulation} function allows to visualize the IRFs.
}
\seealso{
For details, see \code{\link{gecon_simulation-class}}.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
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
summary(irf_rbc)
plot_simulation(irf_rbc)
}