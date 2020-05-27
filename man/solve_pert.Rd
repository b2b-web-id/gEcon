\name{solve_pert}
\alias{solve_pert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve the model in a (log-)linearised form (the first order perturbation)
}
\description{
The function solves the first order perturbation of a model in a (log-)linearised form
using Christopher Sims' \code{gensys} solver.
}
\usage{
solve_pert(model, loglin = TRUE,
           loglin_var = NULL, not_loglin_var = NULL,
           tol = 1e-6, solver = "gensys")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{loglin}{
a logical value. If TRUE, all variables are selected for log-linearisation.
}
  \item{loglin_var}{
a vector of variables that are to be log-linearised (effective only if the \code{loglin} argument is set to FALSE).
}
  \item{not_loglin_var}{
a vector of variables that are not to be log-linearised (overrides previous settings).
}
  \item{tol}{
a numeric value. Tolerance level of a solution (1e-6 by default).
}
  \item{solver}{
the name of the first order perturbation solver. The default solver is Christopher Sims' \code{gensys} solver.
}
}
\details{
Cf. gEcon users' guide, chapter `Solving the model in linearised form'.
}
\value{
an object of \code{gecon_model} class representing the model.
Generic functions such as \code{print} and \code{summary} allow to show the model's components.
The \code{\link{get_pert_solution}} function returns computed recursive laws of motion of the model's variables.
The \code{\link{check_bk}} function displays the eigenvalues of the system and checks the Blanchard-Kahn conditions.
}
\references{
Sims, Ch. A. (2002), Solving Linear Rational Expectations Models, \emph{Computational Economics}
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())
# make and load the model
rbc <- make_model("rbc.gcn")
# find the steady state
rbc <- steady_state(rbc)

# solve in log-linearised form
rbc <- solve_pert(rbc)
get_pert_solution(rbc)

# solve in linearised form
rbc <- solve_pert(rbc, loglin = FALSE)
get_pert_solution(rbc)

# solve with all variables except L_s log-linearised
rbc <- solve_pert(rbc, not_loglin_var = c("L_s"))
get_pert_solution(rbc)
}
