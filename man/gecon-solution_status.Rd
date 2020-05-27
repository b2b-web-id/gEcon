\name{gecon-solution_status}
\alias{re_solved}
\alias{ss_solved}
\title{Model solution status}
\usage{
ss_solved(model)

re_solved(model)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
}
\description{
Functions allowing to check the solution status of \code{gecon_model} objects.
}
\section{Value}{
 \subsection{\code{ss_solved}}{
             TRUE, if the steady state (equilibrium) of the model has been found. FALSE otherwise.
     }


 \subsection{\code{re_solved}}{
             TRUE, if the perturbation has been solved. FALSE otherwise.
     }
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# before the model is solved both return FALSE
ss_solved(rbc)
re_solved(rbc)

# find the steady state
rbc <- steady_state(rbc)
# solve the model in log-linearised form
rbc <- solve_pert(rbc)
get_pert_solution(rbc)

# after the model has been solved both return TRUE
ss_solved(rbc)
re_solved(rbc)
}
