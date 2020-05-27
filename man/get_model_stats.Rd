\name{get_model_stats}
\alias{get_model_stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistics of the model
}
\description{
The \code{get_model_stats} function prints and returns the statistics of the model (absolute and relative to the reference variable).
}
\usage{
get_model_stats(model, variables = NULL,
                basic_stats = TRUE,
                corr = TRUE,
                autocorr = TRUE,
                var_dec = TRUE,
                to_tex = FALSE, silent = FALSE)
}
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or the indices of the variables of interest. All variables are selected by default.
}
  \item{basic_stats}{
a logical value. If TRUE, the following information is returned for selected variables:
steady-state value, standard deviation, variance, and information whether a variable
has been log-linearised. The default value is TRUE.
}
  \item{corr}{
a logical value. If TRUE, a correlation matrix is returned. If a reference variable was not NULL
while invoking the 'compute_model_stats' function, then correlations of selected variables with leads and lags of
the reference variable are also returned. The default value is TRUE.
}
  \item{autocorr}{
a logical value. If TRUE, autocorrelations of selected variables are returned. The default value is TRUE.
}
  \item{var_dec}{
a logical value. If TRUE, variance decomposition (contributions of shocks to the variables' variances) is returned.
The default value is TRUE.
}
\item{to_tex}{
a logical value. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
a logical value. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
The function returns a list of model statistics, which may contain the following fields:
\itemize{
\item \code{basic_stats} - a data frame with steady-state values, standard deviations, variances, and information whether variables have been log-linearised,
\item \code{corr} - a correlation matrix,
\item \code{corr_refvar} - a matrix of correlations of selected variables with the reference variable (its lags and leads),
\item \code{autocorr} - a matrix of autocorrelations,
\item \code{var_dec} - a matrix of variance decomposition, i.e. contributions of shocks to the variables' variances.
}
}
\seealso{
The \code{\link{compute_model_stats}} function computes statistics of the model using spectral (FFT) or simulation methods.
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
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")

# compute and retrieve model statistics
rbc <- compute_model_stats(rbc,  ref_var = "Y")
get_model_stats(model = rbc)
}

