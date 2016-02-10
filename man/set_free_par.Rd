\name{set_free_par}
\alias{set_free_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting free parameters of model
}
\description{
The \code{set_free_par} function specifies parameters of a \code{gecon_model} class object.
}
\usage{
set_free_par(model, free_par, reset = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class \code{gecon_model}.
}
  \item{free_par}{
a named list or a vector of parameters.
}
  \item{reset}{
logical value. If TRUE, the function allows to reset free parameters to values specified in the .gcn file.
}
}
\value{
An object of the \code{gecon_model} class representing the model. If the \code{reset} option has been set to TRUE, the model's parameters will be set back to values from the .gcn file. Generic functions such as \code{print} and \code{summary} allow to show model elements. The \code{\link{get_par_values}} function returns parameter values currently in use.
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

# set different free parameters values from those in the .gcn file
rbc <- set_free_par(rbc, free_par = list(beta = 0.98, delta = 0.01))
rbc <- steady_state(rbc, options = list(method = 'Broyden',
                                        global = 'gline'))
get_ss_values(rbc)

# reset values to .gcn file values
rbc <- set_free_par(rbc, reset = TRUE)
rbc <- steady_state(rbc)
get_ss_values(rbc)
}
