\name{set_shock_cov_mat}
\alias{set_shock_cov_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting a covariance matrix of stochastic shocks.
}
\description{
The \code{set_shock_cov_mat} function allows to set a covariance matrix of model shocks.
}
\usage{
set_shock_cov_mat(model, cov_matrix, shock_order = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{cov_matrix}{
a numeric matrix. Shock covariance matrix.
}
  \item{shock_order}{
a character vector specifying the order of shocks in the cov_matrix.
If not specified, shocks will be sorted according to their internal ordering (as reported by the \code{get_shock_names} function).
}
}
\value{
An (updated) object of \code{gecon_model} class representing the model.
}
\details{
The order of rows/columns of shock covariance matrix must agree with the internal order of shocks in a corresponding
\code{gecon_model-class} object, unless the \code{shock_order} argument is supplied.
Shocks' internal order can be checked with the \code{\link{shock_info}} function as well as the generic function \code{print}.
}
\seealso{
Generic functions such as \code{print} and \code{summary} allow to show the model's components.
The \code{\link{shock_info}} function returns the names of shocks, information about equations which they appear in, and their current covariance matrix.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")
shock_info(rbc, all = TRUE)

# set the shock covariance matrix
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
shock_info(rbc, all = TRUE)
}
