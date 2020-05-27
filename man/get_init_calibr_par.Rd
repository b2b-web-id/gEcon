\name{get_init_calibr_par}
\alias{get_init_calibr_par}
\title{
Calibrated parameters' initial values
}
\description{
The \code{get_init_calibr_par} function prints and returns calibrated parameters' initial values.
}
\usage{
get_init_calibr_par(model, calibr_par = NULL, to_tex = FALSE, silent = FALSE)
}
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{calibr_par}{
the names or the indices of calibrated parameters whose initial values are to be returned. All calibrated parameters are listed by default.
}
\item{to_tex}{
logical. If TRUE, the output is written to a .tex file. The default value is FALSE.
}
\item{silent}{
logical. If TRUE, console output is suppressed. The default value is FALSE.
}
}
\value{
This function returns calibrated parameters' initial values.
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

# get calibrated parameters' initial values
get_init_calibr_par(rbc)
}
