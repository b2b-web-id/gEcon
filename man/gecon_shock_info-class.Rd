\name{gecon_shock_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_shock_info-class}

\title{Class \code{"gecon_shock_info"}}
\description{
The class storing information about shocks chosen by the user.
}
\section{Objects from the Class}{
Objects of this class are created when calling the \code{\link{shock_info}} function.
}
\section{Slots}{
\describe{
    \item{\code{model_info}:}{a character vector of length 3, containing information about the model: the input file name, the input file path, and the date of creation.}
    \item{\code{model_variable_name}:}{a character string denoting the name of a variable storing the model for which simulations have been performed.}
    \item{\code{shock_names}:}{a character vector of the shock names.}
    \item{\code{shock_matrix}:}{a \code{matrix} object containing columns of the variance-covariance matrix corresponding to given shocks.}
    \item{\code{shock_matrix_flag}:}{logical. Set to TRUE when the user specifies non-default entries in a variance-covariance matrix of shocks.}
    \item{\code{incid_mat}:}{a \code{Matrix} object representing the mapping of shocks to equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_shock_info")}: Prints all the available information (the incidence matrix, the variance-covariance matrix) about the shocks.}
    \item{show}{\code{signature(object = "gecon_shock_info")}: Prints the incidence matrix and the variance-covariance matrix of shocks specified when creating a \code{gecon_shock_info} object. }
    \item{summary}{\code{signature(object = "gecon_shock_info")}: Prints all the available information (the incidence matrix, the variance-covariance matrix) about the shocks.}
	 }
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{shock_info}} to create a \code{gecon_shock_info} object. The analogous classes storing the information about
variables and parameters are \code{\link{gecon_var_info-class}} and \code{\link{gecon_par_info-class}}.
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
rbc <- set_shock_cov_mat(rbc, shock_matrix = matrix(0.01, 1, 1), shock_order = 'epsilon_Z')
rbc <- compute_moments(rbc,  ref_var = 'Y')

# display info about the shocks
rbc_shock_info <- shock_info(rbc, all_shocks = TRUE)
print(rbc_shock_info)
summary(rbc_shock_info)
show(rbc_shock_info)
}
