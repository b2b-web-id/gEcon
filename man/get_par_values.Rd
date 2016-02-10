\name{get_par_values}
\alias{get_par_values}
\title{
Parameters of the model
}
\description{
The \code{get_par_values} function prints and returns the values of parameters.
}
\usage{
get_par_values(model, par_names, to_tex)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{par_names}{
a list of requested parameters names.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
}
\value{
This function returns both free and calibrated parameter values.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
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

# find the steady state
rbc <- steady_state(rbc)

# get parameter values
get_par_values(rbc)
}
