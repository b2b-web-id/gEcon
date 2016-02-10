\name{list_calibr_eq}
\alias{list_calibr_eq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find calibrating equations
}
\description{
The \code{list_calibr_eq} function returns calibrating equations with given indices.
}
\usage{
list_calibr_eq(model, no_eq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon model} class.
}
  \item{no_eq}{
a numeric variable, specifies the indices of requested equations.
}
}
\value{
A character vector of requested equations.
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

# list calibrating equations
list_calibr_eq(rbc, c(1))
}