\name{gecon_shock_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_shock_info-class}

\title{\code{gecon_shock_info} class}
\description{
The class storing information about shocks selected by the user.
}
\section{Slots}{
\describe{
    \item{\code{r_object_name}:}{a character string. The name of an R object of \code{gecon_model} class storing the model for which simulations have been performed.}
    \item{\code{shocks}:}{a character vector of shocks' names.}
    \item{\code{cov_matrix}:}{a numeric matrix containing columns of shock covariance matrix corresponding to selected shocks.}
    \item{\code{cov_matrix_flag}:}{a logical value. Set to TRUE when the user enters non-default data into a covariance matrix of shocks.}
    \item{\code{incid_mat}:}{an object of \code{Matrix} class representing the mapping of shocks to equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_shock_info")}: Prints all the available information (the incidence matrix, the covariance matrix) about the shocks which were selected when creating a \code{gecon_shock_info-class} object.}
    \item{show}{\code{signature(object = "gecon_shock_info")}: Prints all the available information (the incidence matrix, the covariance matrix) about the shocks which were selected when creating a \code{gecon_shock_info-class} object.}
    \item{summary}{\code{signature(object = "gecon_shock_info")}: Prints all the available information (the incidence matrix, the covariance matrix) about the shocks which were selected when creating a \code{gecon_shock_info-class} object.}
 }
}
\seealso{
\code{\link{shock_info}} creates an object of \code{gecon_shock_info} class.
\code{\link{gecon_var_info-class}} and \code{\link{gecon_par_info-class}} are
similar classes storing the information about variables and parameters, respectively.
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
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(0.01, 1, 1),
                         shock_order = "epsilon_Z")
rbc <- compute_model_stats(rbc,  ref_var = "Y")

# create an object of gecon_shock_info class and display info about the shocks
rbc_shock_info <- shock_info(rbc, all = TRUE)
print(rbc_shock_info)
summary(rbc_shock_info)
show(rbc_shock_info)
}
