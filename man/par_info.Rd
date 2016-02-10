\name{par_info}
\alias{par_info}
\title{
Information about parameters
}
\description{
The \code{par_info} function allows to create an object of class \code{gecon_par_info},
containing information about given parameters of the chosen model. It allows to check type and value of a given set of parameters,
and the incidence matrix.
}
\usage{
par_info(model, par_names = NULL, all_parameters = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{par_names}{
the names of the parameters of interest.
}
  \item{all_parameters}{
the logical value. If set to TRUE, the par_names argument is overwritten with a vector of all parameters
appearing in the model. The default value is FALSE.
}
}
\details{
If the function result is not assigned to any variable, the information about the requested parameters
is printed in the console.
}
\value{
An object of \code{gecon_par_info} class.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{shock_info}} for information about the shocks and \code{\link{var_info}} for information about the variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# find the steady state and the values of the calibrated parameters
rbc <- steady_state(rbc)

# get information about all parameters
par_info(rbc, all_parameters = TRUE)
}
