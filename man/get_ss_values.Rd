\name{get_ss_values}
\alias{get_ss_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the steady-state (equilibrium) values
}
\description{
The \code{get_ss_values} function returns (and prints) the steady state of the model for dynamic models (equilibrium for static models).
}
\usage{
get_ss_values(model, var_names = NULL, to_tex = FALSE, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{var_names}{
the names or the indices of the variables, whose steady-state values (equilibrium values) are to be returned. The default option is a vector containing all the variable names.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
A numeric vector of the steady-state (equilibrium for static models) values.
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

# find and print the steady state values
rbc <- steady_state(rbc)
get_ss_values(rbc)
get_ss_values(rbc, var_names=c('K_s', 'L_s'))
}