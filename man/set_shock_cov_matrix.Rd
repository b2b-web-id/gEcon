\name{set_shock_cov_mat}
\alias{set_shock_cov_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting a variance-covariance matrix of stochastic shocks.
}
\description{
The \code{set_shock_cov_mat} function allows to set a variance-covariance matrix for the model shocks.
}
\usage{
set_shock_cov_mat(model, shock_matrix, shock_order = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of the \code{gecon_model} class.
}
  \item{shock_matrix}{
a symmetric, positive definite matrix with the dimensions (n * n), where n is the number of shocks in the model.
}
  \item{shock_order}{
a character vector specifying the order of shocks in the shock_matrix.
If not specified, it is assumed that the order is in accordance
with the internal order of the model. The default order
can be displayed by using the \code{\link{shock_info}} function with
the \code{all_shocks} argument set to TRUE.
}
}
\value{
An object of the \code{gecon_model} class, which is representing the model.
Generic functions such as \code{print} and \code{summary} allow to show the model elements. The \code{\link{shock_info}} function returns names of shocks, information about which equations they appear in and the current variance-covariance matrix.
}
\details{
The rows and columns of shock matrix must agree with the order of shocks stored in a \code{gecon_model-class} object if the \code{shock_order} argument is not supplied. This order can be checked by using the \code{\link{shock_info}} function and the generic function \code{print}.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')
shock_info(rbc, all_shocks = TRUE)

# set the shock covariance matrix
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
shock_info(rbc, all_shocks = TRUE)
}
