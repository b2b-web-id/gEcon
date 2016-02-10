\name{get_cov_mat}
\alias{get_cov_mat}
\title{Accessing a variance-covariance matrix of model shocks.}
\description{
The \code{get_cov_mat} function returns a variance-covariance
matrix of model shocks.
}
\usage{
get_cov_mat(model)
}
\arguments{
\item{model}{an object of \code{gecon_model} class.}
}
\value{
The function returns a variance-covariance matrix of model shocks.
}
\seealso{
For details, see \code{\link{gecon_model-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'home_production.gcn'), to = getwd())

# make and load the model
home_production <- make_model('home_production.gcn')
shock_info(home_production, all_shocks = TRUE)

# set the shock distribution parameters
home_production <- set_shock_distr_par(home_production,
                                       distr_par = list('sd(epsilon_h)' = 0.7,
                                                        'var(epsilon_m)' = 0.49,
                                                        'cor(epsilon_m,
                                                             epsilon_h)' = 2/3))

# retrieve and show the variance-covariance matrix of model shocks
cov_mat <- get_cov_mat(home_production)
cov_mat
}
