\name{check_bk}
\alias{check_bk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blanchard-Kahn conditions and eigenvalues
}
\description{
The \code{check_bk} function checks the Blanchard-Kahn conditions and prints information about eigenvalues.
}
\usage{
check_bk(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
}
\details{
The function checks if the Blanchard-Kahn conditions have been satisfied and prints info about eigenvalues
larger than 1 in modulus and the number of forward looking variables.
Eigenvalues are computed when \code{gEcon} attempts to solve the 1st order perturbation
(solver uses the Lapack \code{zgges} function  to compute eigenvalues).
The \code{\link{solve_pert}} function must be called before eigenvalues can be retrieved.
}
\references{
Blanchard, O., Kahn C. M. (1980), The Solution of Linear Difference Models
under Rational Expectations, \emph{Econometrica}
}
\examples{
# copy the example to the current working directory, process it
# and solve for the steady state
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())
rbc <- make_model("rbc.gcn")
rbc <- steady_state(rbc)

# solve the model in log-linearised form
rbc <- solve_pert(rbc)

# check eigenvalues
check_bk(rbc)
}