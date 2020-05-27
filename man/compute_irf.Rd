\name{compute_irf}
\alias{compute_irf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute impulse response functions (IRFs)
}
\description{
The \code{compute_irf} function computes the impulse response functions for selected variables
and shocks and returns an object of \code{gecon_simulation} class.
}
\usage{
compute_irf(model, variables = NULL, shocks = NULL,
            sim_length = 40, cholesky = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or indices of variables whose responses are to be simulated.
}
  \item{shocks}{
the names or indices of shocks for which IRFs are to be computed.
If missing, IRFs are computed for all shocks with non-zero variance.
}
  \item{sim_length}{
the length of simulation path, the default value is 40.
}
  \item{cholesky}{
a logical value. If set to FALSE, IRFs are computed with initial values of all shocks equal to 1,
otherwise the Cholesky decomposition of shock covariance matrix is used (the default).
}
}
\details{
Cf. gEcon users' guide, chapter `Model analysis'.
}
\value{
The function returns an object of \code{gecon_simulation} class.
}
\seealso{
For details, see \code{\link{gecon_simulation-class}}.
Generic functions such as \code{print} and \code{summary} provide information
about the impulse response functions.
The \code{plot_simulation} function allows to plot the IRFs.
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

# compute and print the IRFs
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc_irf <- compute_irf(rbc, variables = c("K_s", "C", "Z", "I", "Y"))
summary(rbc_irf)
plot_simulation(rbc_irf)
}