\name{set_shock_distr_par}
\alias{set_shock_distr_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting distribution parameters of model shocks
}
\description{
The \code{set_shock_distr_par} function sets distribution parameters
(standard deviations, correlations, etc.) of shocks in an object of \code{gecon_model} class.
}
\usage{
set_shock_distr_par(model, distr_par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of \code{gecon_model} class.
}
  \item{distr_par}{
a list or a vector of distribution parameters with named elements.
}
}
\value{
An (updated) object of the \code{gecon_model} class representing the model.
}
\details{
By default, \code{gEcon} uses an identity matrix as the covariance matrix of shocks.
Valid parameters' names should match any of the following patterns:

\code{"sd( SHOCK_NAME )"}

\code{"var( SHOCK_NAME )"}

\code{"cov( SHOCK_NAME_1, SHOCK_NAME_2 )"}

\code{"cor( SHOCK_NAME_1, SHOCK_NAME_2 )"}

There are two issues which the user should be careful about while
using the \code{set_shock_distr_par} function. First, in contrast to other parameters,
shock distribution parameters require quotation marks to be assigned properly.
If quotation marks are omitted, R parser treats elements of the \code{distr_par} list
or vector as functions and attempts to evaluate them, producing errors.
Second, parameters passed to the \code{distr_par} argument should not be specified twice.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "home_production.gcn"), to = getwd())

# make and load the model
home_production <- make_model("home_production.gcn")
shock_info(home_production, all = TRUE)

# set the shock distribution parameters
home_production <- set_shock_distr_par(home_production,
                                       distr_par = list("sd(epsilon_h)" = 0.7,
                                                        "var(epsilon_m)" = 0.49,
                                                        "cor(epsilon_m,
                                                             epsilon_h)" = 2/3))
# get information about shocks in the model
shock_info(home_production, all = TRUE)
}
