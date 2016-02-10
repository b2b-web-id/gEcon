\name{gecon_model-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_model-class}


\title{Class \code{"gecon_model"}}
\description{
The class for storing models.
}
\section{Objects from the Class}{
Objects can be created by calling \code{new("gecon_model", ...)} or preferably using \code{\link{gecon_model}} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{model_info}:}{a character vector of length 3, containing information about the model: the input file name, the input file path, and the date of creation.}
    \item{\code{index_sets}:}{a list containing information about index sets. The names of the list components correspond to the set names. Each component contains character vector of the names of the relevant set elements.}
    \item{\code{parameters}:}{a character vector of all parameter names.}
    \item{\code{parameters_free}:}{a character vector of free parameter names.}
    \item{\code{map_free_into_params}:}{the mapping of free parameters to parameters indices.}
    \item{\code{parameters_free_mod_flag}:}{a logical vector of flags denoting if particular free parameter has been modified with respect to the .gcn file calibration}
    \item{\code{parameters_calibr}:}{a character vector of the names of calibrated parameters.}
    \item{\code{parameters_calibr_mod_flag}:}{a logical vector of the flags denoting if any non-default value has been given for the calibrated parameter.}
    \item{\code{map_calibr_into_params}:}{a numeric vector of the mapping of free parameters to parameters indices.}
    \item{\code{variables}:}{a character vector of the names of variables.}
    \item{\code{shocks}:}{a character vector of the names of shocks.}
    \item{\code{equations}:}{a character vector of the names of equations.}
    \item{\code{calibr_equations}:}{a character vector of the names of calibrating equations.}
    \item{\code{var_eq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of variables to equations.}
    \item{\code{shock_eq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of shocks to equations.}
    \item{\code{var_ceq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of variables to calibrating equations.}
    \item{\code{cpar_eq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of calibrated parameters to equations.}
    \item{\code{cpar_ceq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of calibrated parameters to calibrating equations.}
    \item{\code{fpar_eq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of free parameters to equations.}
    \item{\code{fpar_ceq_map}:}{a sparse matrix (a \code{Matrix} class object) representing the mapping of free parameters to calibrating equations.}
    \item{\code{init_calib_pars_supplied}:}{logical. It informs if calibrated parameters values have been supplied.}
    \item{\code{init_vals_supplied}:}{logical. It informs if values of variables have been supplied.}
    \item{\code{loglin}:}{logical. It informs if the model has to be log-linearised. The default value is TRUE.}
    \item{\code{loglin_var}:}{logical. Flags are set to TRUE for log-linearised variables.}
    \item{\code{re_solved}:}{logical. It is set to TRUE if the model has been solved. The default value is FALSE.}
    \item{\code{corr_computed}:}{logical. If TRUE, indicates that the correlations and other statistics of variables have been computed. The default value is FALSE.}
    \item{\code{is_stochastic}:}{logical. If TRUE, the model has stochastic shocks.}
    \item{\code{is_dynamic}:}{logical. If TRUE, the model has any lead or lagged variables.}
    \item{\code{is_calibrated}:}{logical. If TRUE, the model takes into account calibrating equations when solving for the steady state for a dynamic model (equilibrium in case of static model).}
    \item{\code{parameters_free_init_val}:}{a vector of free parameters values which have been declared in .gcn file.}
    \item{\code{parameters_free_val}:}{a vector of current free parameter values.}
    \item{\code{parameters_calibr_val}:}{a vector of current values of calibrated parameters.}
    \item{\code{params}:}{a vector of the model parameters.}
    \item{\code{steady}:}{a vector of the steady-state values of variables for dynamic models or equilibrium for static models.}
    \item{\code{ss_function}:}{a function returning the steady state (equilibrium for static models) equations residuals.}
    \item{\code{ss_function_jac}:}{a function computing the Jacobian of steady-state function (equilibrium for static models).}
    \item{\code{ss_calibr_function_jac}:}{a function computing the Jacobian of both steady-state (equilibrium) and calibrating functions.}
    \item{\code{calibr_function}:}{calibrating functions}
    \item{\code{init_residual_vector}:}{a numeric vector of residuals of the steady-state (equilibrium) function computed for intial values and calibrating parameters.}
    \item{\code{residual_vector}:}{a numeric vector of residuals of the steady-state (equilibrium) function  computed for the values of variables and calibrating parameters after the nonlinear solver has exited.}
    \item{\code{solver_status}:}{a character string describing the steady-state (equilibrium) solver status.}
    \item{\code{ss_solved}:}{logical. If TRUE, steady state (equilibrium for static models) has been found.}
    \item{\code{pert}:}{functions returning the perturbation of first order (returning a list of matrices).}
    \item{\code{eig_vals}:}{a matrix of system eigenvalues.}
    \item{\code{solution}:}{a list with elements P, Q, R, S storing solution of the model.}
    \item{\code{state_var_indices}:}{a numeric vector containing the indices of state variables.}
    \item{\code{solver_exit_info}:}{a character string containing information about perturbation solver exit information.}
    \item{\code{solution_resid}:}{residuals of checking equations, verifying if the model has been solved.}
    \item{\code{active_shocks}:}{a logical vector of the length equal to the number of shocks. If entry is set to FALSE, the shock is not taken into account while performing stochastic simulations of the model.}
    \item{\code{cov_mat}:}{a variance-covariance matrix of model shocks.}
    \item{\code{shock_mat_flag}:}{logical. Set to TRUE when the user specifies non-default entries in variance-covariance matrix for shocks.}
    \item{\code{corr_mat}:}{a matrix of the model variables correlations.}
    \item{\code{autocorr_mat}:}{a matrix of the model variables autocorrelations.}
    \item{\code{corr_variable_mat}:}{a matrix of correlations of variables with the reference variable lead and lagged values.}
    \item{\code{var_position}:}{a numeric value indicating position of reference variable for the computation of statistics.}
    \item{\code{var_dec}:}{a matrix of variance decomposition of shocks.}
    \item{\code{sdev}:}{a vector of standard deviations of variables.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_model")}: prints information about the model solution status.}
    \item{show}{\code{signature(object = "gecon_model")}: prints short information about the model solution status.}
    \item{summary}{\code{signature(object = "gecon_model")}: prints detailed results of the model.}
	 }
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
print(rbc)
class(rbc)
}
