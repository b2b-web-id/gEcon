\name{get_moments}
\alias{get_moments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistics of the model
}
\description{
The \code{get_moments} function prints and returns the statistics of the model (absolute and relative to the reference variable).
}
\usage{
get_moments(model, var_names, relative_to = FALSE, moments = TRUE,
            correlations = TRUE, autocorrelations = TRUE, var_dec = TRUE,
            to_tex = FALSE)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{var_names}{
the names of the variables of interest.
}
  \item{relative_to}{
logical. If TRUE, the function returns moments relative to one of the variables in accordance with relevant options chosen in the \code{\link{compute_moments}} (then only 'moments' and 'correalations' are active options). The default value is FALSE.
}
  \item{moments}{
logical. If TRUE, the moments of variables: steady state values, standard deviations and variances are returned with the information about which variables have been log-linearised.
  If 'relative_to' is set to TRUE then the moments and steady-state values relative to the reference variable are returned.
}
  \item{correlations}{
logical. If TRUE, a correlation matrix is returned. If relative_to is set to TRUE, then the corraletions of variables with lagged and leading values of a chosen variable are returned. The default value is TRUE.}
  \item{autocorrelations}{
logical. If TRUE then the autocorrelations of variables are returned. If the \code{relative_to} is set to TRUE, this option is inactive. The default value is TRUE.
}
  \item{var_dec}{
logical. If TRUE then the variance decomposition (of shocks) is returned. If the \code{relative_to} is set to TRUE, the option is inactive. The default value is TRUE.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
}
\value{
The function returns a list of absolute or relative moments of variables depending on the value of the \code{relative_to} argument.

When the \code{relative_to} is set to FALSE, the list may consist of the following elements:
\itemize{
\item \code{moments} - means, standard deviations, and variances of variables,
\item \code{correlation_matrix} - a matrix of correlation of variables,
\item \code{autocorrelations} - a matrix of correlation of variables with their own lagged values (autocorrelations),
\item \code{variance_decomposition} - the variance decomposition, describing the amount of variable variability that can be ascribed to each of shocks.
}
When the \code{relative_to} is set to TRUE, the list may consist of two elements:
\itemize{
\item \code{relative_moments} - means, standard deviations, and variance of variables with respect to reference variable specified in the \code{compute_moments} function,
\item \code{correlations_variable} - a matrix of correlation of variables with lead and lagged values of a reference variable (usually GDP).
}
}
\seealso{
the \code{\link{compute_moments}} function to see how the statistics are computed.
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
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')

# compute and retrieve correlations
rbc <- compute_moments(rbc,  ref_var = 'Y')
get_moments(model = rbc, relative_to = FALSE, moments = TRUE, correlations = TRUE,
            autocorrelations = TRUE, var_dec = TRUE)
get_moments(model = rbc, relative_to = TRUE, moments = TRUE, correlations = TRUE)
}

