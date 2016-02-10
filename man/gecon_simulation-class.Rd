\name{gecon_simulation-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_simulation-class}


\title{Class \code{"gecon_simulation"}}
\description{
The class storing simulation results.
}
\section{Objects from the Class}{
Objects can be created by calling the \code{new("gecon_simulation", ...)} form or (preferably) using \code{\link{gecon_simulation}} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{sim}:}{a three-dimensional array with impulse response functions (the dimensions are variables, time, shocks) or two-dimensional array
   when storing the results of user-specified path of shocks or random path of shocks.}
    \item{\code{shock_list}:}{a vector of shocks for which simulations have been computed.}
    \item{\code{var_list}:}{a vector of names of simulated variables.}
    \item{\code{sim_type}:}{a type of simulation.}
    \item{\code{time_n}:}{the number of simulation periods.}
     \item{\code{model_info}:}{a character vector of length 3, containing information about the model: the input file name, the input file path and the date of creation.}
    \item{\code{model_variable_name}:}{a character string denoting the name of a variable storing the model for which the simulations have been performed.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_simulation")}: prints diagnostic information about the simulation performed. }
    \item{show}{\code{signature(object = "gecon_simulation")}: prints short information about the simulation. }
    \item{summary}{\code{signature(object = "gecon_simulation")}: prints and returns the simulation results in the form of list.}
	 }
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{get_simulation_results}} to retrieve the simulated series from \code{sim} slot.
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
irf_rbc <- compute_irf(rbc)
summary(irf_rbc)
class(irf_rbc)
}
