\name{gecon_simulation-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_simulation-class}


\title{\code{gecon_simulation} class}
\description{
The class storing simulation results.
}
\section{Slots}{
  \describe{
    \item{\code{sim}:}{a three-dimensional array with impulse response functions (the dimensions: variables, time, shocks)
                       or two-dimensional array with simulation results.}
    \item{\code{shocks}:}{a character vector with the names of the shocks for which the simulations have been performed.}
    \item{\code{shocks_tex}:}{a character vector with LaTeX names of the shocks for which the simulations have been performed.}
    \item{\code{variables}:}{a character vector with the names of the variables for which the simulations have been performed.}
    \item{\code{variables_tex}:}{a character vector with LaTeX names of the variables for which the simulations have been performed.}
    \item{\code{sim_name}:}{a character string, the simulation name.}
    \item{\code{model_info}:}{a character vector of length 3 containing information about the model: the input file name, the input file path and the date of creation.}
    \item{\code{r_object_name}:}{a character string, the name of an R object storing the model for which the simulations have been performed.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "gecon_simulation")}: prints short information about the simulation.}
    \item{print}{\code{signature(x = "gecon_simulation")}: prints information about the simulation.}
    \item{summary}{\code{signature(object = "gecon_simulation")}: prints the simulation results.}
    }
}
\seealso{
\code{\link{get_simulation_results}} to retrieve the simulated series from the \code{sim} slot.
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
irf_rbc <- compute_irf(rbc)
summary(irf_rbc)
class(irf_rbc)
}
