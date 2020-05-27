\name{var_info}
\alias{var_info}
\title{
Information about variables
}
\description{
The \code{var_info} function allows to create an object of \code{gecon_var_info} class,
containing information about selected model variables. It allows to check variables'
equation incidence matrix as well as the already computed statistics for them.
}
\usage{
var_info(model, variables = NULL, all = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{variables}{
the names or indices of the variables of interest.
}
  \item{all}{
logical value. If TRUE, information about all model variables is generated (FALSE by default).
}
}
\details{
The \code{var_info} function may be useful while debugging a model. It allows also to retrieve information quickly when a model is large.
If the R command is not assigned to any R object, the information about the requested variables is printed to the console.
Depending on which operations have been performed on an object of \code{gecon\_model} class, the information printed consists of a subset of the following elements:
\itemize{
\item incidence information,
\item steady-state (equilibrium) values,
\item initial values,                     
\item variables info (which variables are log-linearised and which are state variables),
\item state variables impact on the selected variables,
\item shocks' impact on the selected variables,
\item basic statistics,
\item correlations.}
}
\value{
An object of \code{\link{gecon_var_info-class}} class.
}
\note{
The function only displays and returns the model's characteristics that have been already set or computed.
Eg. if the model has been solved but the statistics have not been computed, the correlations will not be passed to the object of \code{gecon_var_info} class.
}
\seealso{
\code{\link{shock_info}} for information about the shocks. \cr \code{\link{get_ss_values}}, \code{\link{get_pert_solution}},
\code{\link{get_model_stats}} to extract statistics of the model variables.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# solve the model
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)

# get information about variables
rbc_var_info <- var_info(rbc, variables = c("K_s"))
print(rbc_var_info)
show(rbc_var_info)
summary(rbc_var_info)
}
