\name{get_par_names}
\alias{get_par_names}
\title{Accessing parameter names used by objects of \code{gecon_model} class}
\description{
The \code{get_par_names} function allows to
retrieve the names of parameters from objects of \code{gecon_model} class.
}
\usage{
get_par_names(model, free_par = TRUE, calibr_par = TRUE)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
\item{free_par}{logical. If TRUE (default), free parameters' names are
                added to the vector of parameters' names.}
\item{calibr_par}{logical. If TRUE (default), calibrated parameters' names are
                  added to the vector of parameters' names.}
}
\value{
The function returns a character vector of parameters' names stored in a given object of \code{gecon_model} class.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load model
rbc <- make_model("rbc.gcn")

# get names of all parameters
get_par_names(rbc)

# get names of free parameters
get_par_names(rbc, calibr_par = FALSE)

# get names of calibrated parameters
get_par_names(rbc, free_par = FALSE)
}
