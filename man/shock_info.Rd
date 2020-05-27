\name{shock_info}
\alias{shock_info}
\title{
Information about shocks
}
\description{
The \code{shock_info} function prints information about the model's shocks (occurrence in equations, covariance matrix).
It also allows to create an object of the \code{gecon_shock_info} class, which stores this information.
}
\usage{
shock_info(model, shocks = NULL, all = FALSE)
}
\arguments{
 \item{model}{
an object of \code{gecon_model} class.
}
  \item{shocks}{
the names or indices of shocks of interest.
}
  \item{all}{
logical value. If TRUE, information about all model shocks is generated (FALSE by default).
}
}
\value{
An object of the \code{\link{gecon_shock_info-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load model
rbc <- make_model("rbc.gcn")

# set the shock covariance matrix
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")

# get information about shocks
shock_info(rbc, all = TRUE)
}
