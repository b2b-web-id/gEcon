\name{initval_calibr_par}
\alias{initval_calibr_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting initial values of calibrated parameters
}
\description{
The \code{initval_calibr_par} function enables setting the initial values of calibrated parameters 
for the non-linear steady-state (equilibrium) solver or their expected values if calibration is omitted. 
If not set by this function, parameters' values are set to 0.5 by default.
}
\usage{
initval_calibr_par(model, calibr_par, warnings = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{calibr_par}{
a named list or vector of calibrated parameters' initial values.
}
  \item{warnings}{logical, should warnings be displayed?}
}
\details{
The values of calibrated parameters passed to the object of \code{gecon_model} class are treated as initial values for the steady-state solver 
when the user specifies calibrating equations in a .gcn file and requests that the \code{steady_state} function shall use it. 
If the calibration is omitted, the initial values of calibrated parameters are treated as their final values. 
Calibrated parameters have to be set to correct values upon the decision to make calibrating equations inactive.
}
\value{
An updated object of \code{gecon_model} class representing the model. 
Generic functions such as \code{print} and \code{summary} allow to show the model's elements. 
The \code{\link{get_par_values}} function returns parameters' values.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# set initial values for calibrated parameters
rbc <- initval_calibr_par(rbc, calibr_par=list(alpha = 0.4))

# find the steady state and values of calibrated parameters
rbc <- steady_state(rbc)
get_par_values(rbc, c("alpha"))
}
