\name{shock_info}
\alias{shock_info}
\title{
Information about shocks
}
\description{
The \code{shock_info} function allows to create an object of the \code{gecon_shock_info} class,
which contains the information about the model's shocks (occurence in equations and the variance-covariance matrix).
}
\usage{
shock_info(model, shock_names = NULL, all_shocks = FALSE)
}
\arguments{
 \item{model}{
an object of the \code{gecon_model} class.
}
  \item{shock_names}{
the names of shocks of interest.
}
  \item{all_shocks}{
the logical value. If set to TRUE, the \code{shock_names} argument is overwritten with the vector of all shocks
appearing in the model. The default value is FALSE.
}
}
\value{
An object of the \code{\link{gecon_shock_info-class}} class.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load model
rbc <- make_model('rbc.gcn')

# set the shock covariance matrix
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')

# get information about shocks
shock_info(rbc, all_shocks = TRUE)
}
