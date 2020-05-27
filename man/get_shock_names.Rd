\name{get_shock_names}
\alias{get_shock_names}
\title{Accessing shock names used by objects of \code{gecon_model} class}
\description{
The \code{get_shock_names} function allows to
retrieve the names of shocks from objects of \code{gecon_model} class.
}
\usage{
get_shock_names(model)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
}
\value{
The function returns a character vector of shock names stored in a given object of \code{gecon_model} class.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# get shock names
get_shock_names(rbc)
}

