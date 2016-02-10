\name{gecon_simulation}
\alias{gecon_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructor for "gecon_simulation" class object
}
\description{
This function creates an object of \code{gecon_simulation} class.
}
\usage{
gecon_simulation(sim, shock_list, var_list, sim_type, time_n, model_info, model_variable_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
the array of simulation results (three dimensional when the impulse response functions have been computed for more than one shock).
}
  \item{shock_list}{
a list of shocks for which the simulations have been computed. 
}
  \item{var_list}{
a list of variables used.
}
  \item{sim_type}{
a type of simulation.
}
  \item{time_n}{
the number of periods for which the simulation has been performed.
}
 \item{model_info}{
a character vector of length 3, containing information about the model: the input file name, the input file path, and the date of creation.
}
  \item{model_variable_name}{
a character string denoting the name of the model for which the simulation has been performed.
}
}
\value{
An object of the \code{gecon_simulation} class.
}
\seealso{
Generic functions such as \code{print} and \code{summary} can provide information about the simulations. 
The \code{\link{plot_simulation}} function allows to visualize the impact on variables.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
