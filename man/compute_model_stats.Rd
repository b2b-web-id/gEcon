\name{compute_model_stats}
\alias{compute_model_stats}
\title{
Computation of model statistics
}
\description{
This function computes statistics of the model using spectral (FFT) or simulation methods.
}
\usage{
compute_model_stats(model, n_leadlags = 5, ref_var = NULL, lambda = 1600,
                    ngrid = 64 * 16,
                    sim = FALSE, sim_length = 1e5)
}
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{n_leadlags}{
the number of leads and lags of the model's variables for which correlations are to be computed.
}
  \item{ref_var}{
the name or the index of the reference variable with respect to which correlations are to be computed.
}
  \item{lambda}{
HP filter parameter, if it is set to 0 no filtering is performed, 1600 is the default value (quarterly data).
}
\item{ngrid}{
the density of grid used by the Fast Fourier transform (used only if the \code{sim}
option is set to FALSE). It has to be a multiple of 8 and has to be large enough
to guarantee unbiased results.
}
\item{sim}{
a logical value. If TRUE simulation method is used for computing correlations, if FALSE,
the Fast Fourier transform is used.
}
\item{sim_length}{
the length of simulation path (used only if the \code{sim} option is set to TRUE).
}
}
\details{
Cf. gEcon users' guide, chapter `Model analysis'.
}
\value{
An object of \code{gecon_model} class representing the model.
}
\seealso{
Generic functions such as \code{print} and \code{summary} allow to show the model's components.
The \code{\link{get_model_stats}} function returns various statistics of the model (both absolute and relative).
}
\references{
Hamilton. J.D. (1994), Time Series Analysis, \emph{Princeton University Press}
}
\note{
The density of grid used by the FFT has to be large enough (at least 64 * 8) for spectral method to give accurate results.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load model
rbc <- make_model("rbc.gcn")

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# compute and print correlations
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc <- compute_model_stats(rbc, ref_var = "Y")
get_model_stats(rbc)
}