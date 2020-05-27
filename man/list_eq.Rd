\name{list_eq}
\alias{list_eq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List model equations
}
\description{
The \code{list_eq} function returns equations according to the specified indices.
}
\usage{
list_eq(model, eq_idx = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{eq_idx}{
an integer/numeric value/vector specifying the indices of requested equations.
}
}
\value{
A character (one-column) matrix with requested equations.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# get the 3rd and the 5th model equation
list_eq(rbc, c(3, 5))
}
