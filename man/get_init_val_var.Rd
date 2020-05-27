\name{get_init_val_var}
\alias{get_init_val_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variables' initial values
}
\description{
The \code{get_init_val_var} function returns (and prints) initial values of model variables.
}
\usage{
get_init_val_var(model, variables = NULL, to_tex = FALSE, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or the indices of the variables whose initial values are to be returned. All variables are listed by default.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
A numeric vector of initial values of model variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# find the steady state and print variables' initial values
rbc <- steady_state(rbc)
get_init_val_var(rbc)
get_init_val_var(rbc, variables = c("K_s", "L_s"))
}