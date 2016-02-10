\name{initval_calibr_par}
\alias{initval_calibr_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting initial values of calibrated parameters
}
\description{
The \code{initval_calibr_par} function enables setting the initial values of calibrated parameters for the nonlinear solver searching for the steady state of dynamic models (equilibrium for static models) and the values of calibrated parameters. If not set by this function, the default values of parameters are assumed to be 0.5.
}
\usage{
initval_calibr_par(model, calibr_par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{calibr_par}{
a named list or a vector of parameters.
}
}
\details{
The values of parameters passed to the \code{gecon_model} are treated as initial values for the steady-state solver when the user specifies calibrating equations in a .gcn file and requests that \code{steady_state} function shall use it. If the calibration is omitted, the initial values of calibrated parameters are treated as their final values, so one has to specify the right set of calibrated parameters values when decides to omit the calibrating equations.
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model's elements. The \code{\link{get_par_values}} function return parameter values.
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

# set initial values for the calibrated parameters
rbc <- initval_calibr_par(rbc, calibr_par=list(alpha = 0.4))

# find the steady state and values of the calibrated parameters
rbc <- steady_state(rbc)
get_par_values(rbc, c('alpha'))
}
