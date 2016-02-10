\name{get_pert_solution}
\alias{get_pert_solution}
\title{
Recursive laws of motion for the model variables
}
\description{
The \code{get_pert_solution} function prints and returns the recursive laws of motion for the model's variables.
}
\usage{
# getting recursive laws of motion
get_pert_solution(model, to_tex = FALSE, silent = FALSE)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
A list with \code{P}, \code{Q}, \code{R}, \code{S} elements. \code{P} and \code{Q} matrices denote the impact of lagged state variables and current values of shocks variables on current values of state variables. \code{R} and \code{S} matrices denote the impact of lagged state variables and current values of shocks variables on current values of non-state variables.
}
\seealso{
\code{\link{solve_pert}} for the description of solution procedure and description of output.
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

# find the steady state
rbc <- steady_state(rbc)

# solve the model in log-linearised form and get the results
rbc <- solve_pert(rbc)
get_pert_solution(rbc)
}
