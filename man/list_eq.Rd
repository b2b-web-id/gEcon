\name{list_eq}
\alias{list_eq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find model equations
}
\description{
The \code{list_eq} function returns equations with the specified indices.
}
\usage{
list_eq(model, no_eq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
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

# get the 3rd and the 5th model equation
list_eq(rbc, c(3, 5))
}
