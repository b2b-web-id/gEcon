\name{get_residuals}
\alias{get_residuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieving equations' residuals
}
\description{
The \code{get_residuals} function allows to check the residuals of the steady-state (equilibrium) equations of a dynamic (static) model 
and identify equations with the largest errors. This may help to assign initial values to the model's variables more accurately
when the solver cannot find the steady state (equilibrium).
}
\usage{
get_residuals(model, largest = 5, calibration = TRUE, last_solver_iter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{largest}{
the number of equations with the largest errors which are to be printed.
}
  \item{calibration}{
if FALSE, calibrating equations will not be taken into account when computing equations' residuals.
Initial values of calibrated parameters will be then treated as their values.
}
  \item{last_solver_iter}{
logical. if TRUE, solver searches for a solution starting from the values of last saved iteration of the former search process instead 
of initial values specified by the user. FALSE by default.
}
}
\value{
The function returns a list of two elements: \code{initial} and \code{final}. 
Initial residuals are steady-state (equilibrium) equations' residuals computed using the variables' initial values 
or the values of last saved iteration of the former search process. 
Final residuals are residuals computed after the solver has exited. 
The function prints the indices of equations with the largest initial and final errors. 
The equations can be further investigated using the \code{\link{list_eq}} function.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
                           "home_production_templ.gcn"), to = getwd())

# make and load the model
home_prod_templ <- make_model("home_production_templ.gcn")

# for the purpose of the example, initial values are set very far from the solution
home_prod_templ <- initval_var(home_prod_templ, c(N = 0.02,
                                                  N__H = 0.01,
                                                  N__M = 0.01))
home_prod_templ <- steady_state(home_prod_templ)
get_residuals(home_prod_templ)

# after setting more reasonable values the steady state is found
home_prod_templ <- initval_var(home_prod_templ, c(N = 0.5,
                                                  N__H = 0.25,
                                                  N__M = 0.25))
home_prod_templ <- steady_state(home_prod_templ)
get_residuals(home_prod_templ)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
