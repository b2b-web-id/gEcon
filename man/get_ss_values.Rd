\name{get_ss_values}
\alias{get_ss_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variables' steady-state (equilibrium) values
}
\description{
The \code{get_ss_values} function returns (and prints) the steady state/equilibrium of the model.
}
\usage{
get_ss_values(model, variables = NULL, to_tex = FALSE, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or the indices of the variables whose steady-state (equilibrium) values are to be returned. All variables are listed by default.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
A numeric vector of the steady-state (equilibrium) values of variables.
}
\details{
In case the steady-state solver has been started but has not converged, the function returns a vector of variables' values from the last solver iteration.
}    
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# find and print the steady-state values
rbc <- steady_state(rbc)
get_ss_values(rbc)
get_ss_values(rbc, variables = c("K_s", "L_s"))
}