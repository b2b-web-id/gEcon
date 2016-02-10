\name{get_residuals}
\alias{get_residuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retriving residuals
}
\description{
The \code{get_residuals} function allows to check the residuals of the steady-state equations (equations characterising equilibrium
in case of static models) and identify equations with the highest errors. This information may help to assign better initial values to variables when the solver cannot find the steady state (equilibrium).
}
\usage{
get_residuals(model, highest = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{highest}{
the number of equations with the highest error to be printed.
}
}
\value{
This function returns a list with the \code{initial} and \code{final} elements. The initial residuals are residuals computed using the initial values. The final residuals are residuals computed after the solver has exited. The function prints the indices of equations with the highest initial and final errors. The equations can be investigated by using the \code{\link{list_eq}} function.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
                           'home_prod_templ.gcn'), to = getwd())

# make and load the model
home_prod_templ <- make_model('home_prod_templ.gcn')

# for the purpose of the example, we set the initial values extremely far from the solution
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
