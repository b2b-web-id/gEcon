\name{get_var_names}
\alias{get_var_names}
\title{Accessing variable names used by \code{gecon_model} class objects}
\description{
The \code{get_var_names} function allows to
retrieve the names of variables from \code{gecon_model} class objects.
}
\usage{
get_var_names(model)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
}
\value{
The function returns a character vector of variable names, stored by the chosen object of \code{gecon_model} class.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# get variable names
get_var_names(rbc)
}

