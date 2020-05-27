\name{set_free_par}
\alias{set_free_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting free parameters' values
}
\description{
The \code{set_free_par} function allows to set values of free parameters occurring in a \code{gecon_model} class object.
}
\usage{
set_free_par(model, free_par = NULL, reset = FALSE, warnings = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{free_par}{
a named list or a vector of parameters.
}
  \item{reset}{
a logical value, if TRUE, the function allows to reset free parameters to values from the .gcn file.
}
  \item{warnings}{logical, if TRUE, a warning is displayed whenever the default parameter value (specified in the .gcn file)
        is overwritten.}
}
\value{
An updated object of \code{gecon_model} class representing the model. If the \code{reset} option is set to TRUE, 
the model's parameters will be set back to values from the .gcn file. Generic functions such as \code{print} and \code{summary} 
allow to show model elements. The \code{\link{get_par_values}} function returns parameters' values currently in use.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# set free parameters' values different from those in the .gcn file
rbc <- set_free_par(rbc, free_par = list(beta = 0.98, delta = 0.01))
rbc <- steady_state(rbc, options = list(method = "Broyden",
                                        global = "gline"))
get_ss_values(rbc)

# reset values to .gcn file values
rbc <- set_free_par(rbc, reset = TRUE)
rbc <- steady_state(rbc)
get_ss_values(rbc)
}
