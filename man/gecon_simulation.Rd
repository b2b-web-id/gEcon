\name{gecon_simulation}
\alias{gecon_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create objects of \code{gecon_simulation} class
}
\description{
The function creates an object of \code{gecon_simulation} class.
}
\usage{
gecon_simulation(sim,
                 shocks, shocks_tex, variables, variables_tex,
                 sim_name, model_info, r_object_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
an array with simulation results (3-dimensional if the IRFs have been performed,
2-dimensional otherwise).
}
  \item{shocks}{
a character vector with the names of the shocks for which the simulations have been performed.
}
  \item{shocks_tex}{
a character vector with LaTeX names of the shocks for which the simulations have been performed.
}
  \item{variables}{
a character vector with the names of the variables for which the simulations have been performed.
}
  \item{variables_tex}{
a character vector with LaTeX names of the variables for which the simulations have been performed.
}
  \item{sim_name}{
a character string, the simulation name.
}
 \item{model_info}{
a character vector of length 3 containing information about the model: the input file name, the input file path, and the date of creation.
}
  \item{r_object_name}{
a character string with the name of an R object storing the model for which the simulations have been performed.
}
}
\value{
An object of \code{gecon_simulation} class.
}
\seealso{
Generic functions such as \code{print} and \code{summary} provide information about the simulations.
The \code{\link{plot_simulation}} function allows to visualize the impact of simulations on the model's variables.
}
\note{
The \code{gecon_simulation} constructor is invoked by the \code{random_path},
\code{simulate_model}, and \code{compute_irf} functions.
In general, users should not call this function explicitly.
}
