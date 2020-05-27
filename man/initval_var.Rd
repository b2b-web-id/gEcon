\name{initval_var}
\alias{initval_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting initial values of variables.
}
\description{
The \code{initval_var} function sets the initial values of the model's variables to values specified by the user. The initial values close to solution will help the nonlinear equations solver to find the solution.
}
\usage{
initval_var(model, init_var, warnings = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{init_var}{
a named list or vector of the initial values of variables.
}
  \item{warnings}{logical, should warnings be displayed?}
}
\value{
An object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show model elements. The \code{\link{get_ss_values}} function returns the steady-state (equilibrium) values of the model variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# set initial values of variables
rbc <- initval_var(rbc, list(K_s = 10, C = 2, I = 0.5, Y = 2.5))
rbc <- initval_var(rbc, c(pi=0))
}