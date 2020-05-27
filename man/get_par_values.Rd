\name{get_par_values}
\alias{get_par_values}
\title{
Parameters' values
}
\description{
The \code{get_par_values} function prints and returns parameters' values.
}
\usage{
get_par_values(model, parameters = NULL, to_tex = FALSE, silent = FALSE)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{parameters}{
the names or the indices of the parameters whose values are to be returned. All parameters are listed by default.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
This function returns both free and calibrated parameters' values.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# find the steady state
rbc <- steady_state(rbc)

# get parameters' values
get_par_values(rbc)
}
