\name{solve_pert}
\alias{solve_pert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve the model in a linearised form (1st order perturbation)
}
\description{
This function solves the model in a linearised form using Christopher Sims' \code{gensys} solver.
}
\usage{
solve_pert(model, loglin = TRUE, not_loglin_var = NULL,
           norm_tol = 1e-08, solver = "sims_solver")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{loglin}{
an option to log-linearise the perturbation. If FALSE, the model is only linearised.
}
  \item{not_loglin_var}{
a vector of variables that will not be log-linearised.
}
  \item{norm_tol}{
the tolerance for residuals of model (default 1e-08).
}
  \item{solver}{
the name linear RE solver. The default solver is Christopher Sims' solver. Currently no other solvers are available.
}
}
\details{
Cf. gEcon manual, chapter "Solving the model in linearised form".
}
\value{
an object of the \code{gecon_model} class representing the model. Generic functions such as \code{print} and \code{summary} allow to show the model elements. The \code{\link{get_pert_solution}}  function returns computed recursive laws of motion for the model's variables. The \code{\link{check_bk}} function displays the eigenvalues of the system and checks the Blanchard-Kahn conditions.
}
\references{
Sims, Ch. A. (2002), "Solving Linear Rational Expectations Models.", \emph{Computational Economics}
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

# solve in log-linearised form
rbc <- solve_pert(rbc)
get_pert_solution(rbc)

# solve in linearised form
rbc <- solve_pert(rbc, loglin = FALSE)
get_pert_solution(rbc)

# solve with all variables except L_s log-linearised
rbc <- solve_pert(rbc, not_loglin_var = c('L_s'))
get_pert_solution(rbc)
}
