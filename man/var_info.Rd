\name{var_info}
\alias{var_info}
\title{
Information about variables
}
\description{
The \code{var_info} function allows to create an object of \code{gecon_var_info} class,
which contains information about the chosen variables. It allows to check which equations
given variables appear in. In addition, this function prints the already computed statistics
of the given set of variables.
}
\usage{
var_info(model, var_names = NULL, all_variables = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{var_names}{
the names of the variables of interest.
}
  \item{all_variables}{
the logical value. If set to TRUE, the \code{var_names} argument is overwritten with a vector of all variables
appearing in the model. The default value is FALSE.
}
}
\details{
The \code{var_info} function may be useful in debugging model and quick retrieval of information when the model is large.
If the object returned by the function is not assigned to any variable, the information about the requested parameters
is printed to the console.
One or more of the following elements may be printed:
\itemize{
\item incidence information,
\item steady-state (equilibrium) values,
\item variables info (which ones are log-linearised and which are state variables),
\item state variables impact on the chosen variables,
\item shocks impact on the chosen variables,
\item moments,
\item correlations,} depending on which operations have been performed on \code{gecon\_model} class object.
}
\value{
An object of \code{\link{gecon_var_info-class}} class.
}
\note{
The function only displays the elements of a model that have been already set or computed. Eg. if the model has been solved but the statistics have not been computed, the correlations will not be passed to the \code{gecon_var_info} class.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{shock_info}} for information about the shocks. \cr \code{\link{get_ss_values}}, \code{\link{get_pert_solution}},
\code{\link{get_moments}} to extract the steady-state (equilibrium) values, the solution, and various moments and statistics of the model.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# get information about variables
rbc_var_info <- var_info(rbc, var_names = c('K_s'))
print(rbc_var_info)
show(rbc_var_info)
summary(rbc_var_info)
}
