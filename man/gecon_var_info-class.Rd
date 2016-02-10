\name{gecon_var_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_var_info-class}


\title{Class \code{"gecon_var_info"}}
\description{
The class storing information about variables chosen by the user.
}
\section{Objects from the Class}{
Objects of this class are created when calling the \code{\link{var_info}} function.
}
\section{Slots}{
  \describe{
    \item{\code{model_info}:}{a character vector of length 3, containing information about the model: the input file name, the input file path, and the date of creation.}
    \item{\code{model_variable_name}:}{a character string, the name of a variable storing the model for which the simulations have been performed.}
    \item{\code{var_names}:}{a character vector of the variable names.}
    \item{\code{is_stochastic}:}{logical. If TRUE, the model, based on which the info was generated, has stochastic shocks.}
    \item{\code{is_dynamic}:}{logical. If TRUE, the model, based on which the info was generated, has any lead or lagged variables.}
    \item{\code{ss_solved}:}{logical. If TRUE, the steady state (equilibrium for static models) for the model has been found.}
    \item{\code{re_solved}:}{logical. It is set to TRUE if the model, based on which the info was generated, has been solved. The default value is FALSE.}
    \item{\code{corr_computed}:}{logical. If TRUE, it indicates that the correlations and other statistics of variables have been computed. The default value is FALSE.}
    \item{\code{ss_val}:}{a vector of the steady-state values of variables (dynamic models) or equilibrium (static models). If the steady state has not been computed,
                         this slot contains initial values of variables.}
    \item{\code{state}:}{a logical vector of the length equal to the number of the variables. The TRUE entries denote that a corresponding variable is a state variable.}
    \item{\code{state_var_impact}:}{the rows of the matrices P and R of state space representation corresponding to the chosen variables.}
    \item{\code{shock_impact}:}{the rows of the matrices Q and S of state space representation corresponding to the chosen variables.}
    \item{\code{std_dev_val}:}{a numeric vector of standard deviations of chosen variables.}
    \item{\code{loglin_flag}:}{a logical vector of the length equal to the number of the variables. The TRUE entries denote that a corresponding variable has been loglinearised before solving the model.}
    \item{\code{cr}:}{a matrix containing the correlations of the chosen variables with all the model variables.}
    \item{\code{incid_mat}:}{a \code{Matrix} object representing the mapping of variables to equations and calibrating equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_var_info")}: Prints all the available information (short listing, values, statistics, incidence) about the variables, retrieved from the model when creating a \code{gecon_var_info-class} object. }
    \item{show}{\code{signature(object = "gecon_var_info")}: Prints information about the variables incidence and the results already obtained for the variables. }
    \item{summary}{\code{signature(object = "gecon_var_info")}: Prints all the available information (short listing, values, statistics, incidence) about the variables, retrieved from the model when creating a \code{gecon_var_info-class} object.}
 }
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{var_info}} to create a \code{gecon_var_info} object. The analogous classes storing the information about
shocks and parameters are \code{\link{gecon_shock_info-class}} and \code{\link{gecon_par_info-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# compute the steady state
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)
rbc_var_info <- var_info(rbc, all_variables = TRUE)
print(rbc_var_info)
summary(rbc_var_info)
show(rbc_var_info)
}
