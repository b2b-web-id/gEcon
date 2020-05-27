\name{get_pert_solution}
\alias{get_pert_solution}
\title{
Recursive laws of motion of the model's variables
}
\description{
The \code{get_pert_solution} function prints and returns the recursive laws of motion of the model's variables.
}
\usage{
# getting recursive laws of motion
get_pert_solution(model, to_tex = FALSE, silent = FALSE)
}
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
\item{to_tex}{
a logical value. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
a logical value. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
A list of \code{P}, \code{Q}, \code{R}, \code{S} matrices
(\code{P} and \code{R} only in case of a deterministic model),
a vector of variables' steady-state values (\code{ss_val}),
a vector of flags indicating which variables were log-linearised (\code{loglin}),
and a vector of indices of state variables (\code{state_ind}).
\code{P} and \code{Q} matrices describe the impact of lagged state variables
and current values of shocks on current values of state variables.
\code{R} and \code{S} matrices describe the impact of lagged state variables
and current values of shocks on current values of non-state (jump) variables.
}
\seealso{
\code{\link{solve_pert}} for the details of finding-a-solution procedure and the description of output which is returned.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# find the steady state
rbc <- steady_state(rbc)

# solve the model in log-linearised form and get the results
rbc <- solve_pert(rbc)
get_pert_solution(rbc)
}
