\name{compute_moments}
\alias{compute_moments}
\title{
Computation of correlations
}
\description{
This function computes the statistics of the model by using spectral and simulation methods.
}
\usage{
compute_moments(model, ngrid = (64 * 16), filter = TRUE,
             sim = FALSE, nrun = 1e+05, lambda = 1600,
             ref_var = NULL, n_leadlags = 6)
}
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{ngrid}{
density of grid used by the Fourier transform. It has to be a multiplicity of 8.
}
  \item{filter}{
logical. If TRUE, the statistics are computed for the HP-filtered series, otherwise non-filtered series are used for the statistics computation.
}
  \item{sim}{
logical. If TRUE, simulation methods are used for correlations computations, otherwise spectral methods are used.
}
  \item{nrun}{
the number of MC simulation runs.
}
  \item{lambda}{
the lambda parameter for the HP filter.
}
  \item{ref_var}{
the name or the number of a variable in relation to which correlations are computed. When not specified, the first variable in variables list is treated as the reference value.
}
  \item{n_leadlags}{
the number of leads/lags for computing correlation tables.
}
}
\details{
Cf. gEcon manual, chapter "Model analysis".
}
\value{
An object of \code{gecon_model} class representing a model. Generic functions such as \code{print} and \code{summary} allow to show model elements. The \code{\link{get_moments}} function returns various statistics of the model (both absolute and relative).
}
\references{
Hamilton. J.D. (1994), "Time Series Analysis", \emph{Princeton University Press}
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\note{
The grid has to be large enough (at least 64 * 8)  for spectral methods to converge to simulation means.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load model
rbc <- make_model('rbc.gcn')

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# compute and print correlations
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
rbc <- compute_moments(rbc, ref_var = 'Y')
get_moments(rbc,
            moments = TRUE,
            correlations=TRUE,
            autocorrelations=TRUE,
            var_dec=TRUE)
}