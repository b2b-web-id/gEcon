\name{par_info}
\alias{par_info}
\title{
Information about parameters
}
\description{
The \code{par_info} function allows to create an object of \code{gecon_par_info} class,
containing information about the model's parameters. It allows to check types and values of a set of parameters
as well as their incidence matrix.
}
\usage{
par_info(model, parameters = NULL, all = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{parameters}{
the names or indices of the parameters of interest.
}
  \item{all}{
logical value. If TRUE, information about all model parameters is generated (FALSE by default).
}
}
\details{
If the function's result is not assigned to any variable, the information about the requested parameters
is printed to the console.
}
\value{
An object of \code{gecon_par_info} class.
}
\seealso{
\code{\link{shock_info}} for information about the shocks and \code{\link{var_info}} for information about the variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# find the steady state and the values of calibrated parameters
rbc <- steady_state(rbc)

# get information about all parameters
par_info(rbc, all = TRUE)
}
