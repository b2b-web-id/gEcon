\name{gecon_model-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_model-class}


\title{\code{gecon_model} class}
\description{
The class for storing models.
}
\section{Slots}{
  \describe{
    \item{\code{model_info}:}{a character vector of length 4, containing information about the model: the input file name, the input file path, the date of creation, and a logical flag informing if the model requires compiling.}

    \item{\code{index_sets}:}{a list containing information about index sets. The names of its elements correspond to sets' names. Each element contains a character vector of the names of the relevant set's components.}

    \item{\code{parameters}:}{a character vector of all parameters' names.}
    \item{\code{parameters_tex}:}{a character vector of all parameters' LaTeX names.}
    \item{\code{parameters_free}:}{a character vector of free parameters' names.}
    \item{\code{map_free_into_params}:}{an integer vector of free parameters' indices.}
    \item{\code{parameters_calibr}:}{a character vector of calibrated parameters' names.}
    \item{\code{map_calibr_into_params}:}{an integer vector of calibrated parameters' indices.}
    \item{\code{variables}:}{a character vector of variables' names.}
    \item{\code{variables_tex}:}{a character vector of variables' LaTeX names.}
    \item{\code{shocks}:}{a character vector of shocks' names.}
    \item{\code{shocks_tex}:}{a character vector of shocks' LaTeX names.}

    \item{\code{equations}:}{a character vector of model equations (equilibrium relationships).}
    \item{\code{calibr_equations}:}{a character vector of calibrating equations.}
    \item{\code{var_eq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of variables to equations.}
    \item{\code{shock_eq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of shocks to equations.}
    \item{\code{var_ceq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of variables to calibrating equations.}
    \item{\code{cpar_eq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of calibrated parameters to equations.}
    \item{\code{cpar_ceq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of calibrated parameters to calibrating equations.}
    \item{\code{fpar_eq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of free parameters to equations.}
    \item{\code{fpar_ceq_map}:}{a sparse matrix (an object of \code{Matrix} class) representing the mapping of free parameters to calibrating equations.}

    \item{\code{is_stochastic}:}{logical. If TRUE, the model has stochastic shocks.}
    \item{\code{is_dynamic}:}{logical. If TRUE, the model has at least one lead or lagged variable.}
    \item{\code{is_calibrated}:}{logical. If TRUE, calibrating equations are taken into account when solving for the steady state of a dynamic model (the equilibrium in case of a static model).}

    \item{\code{ss_function}:}{a function returning the steady-state/equilibrium equations' residuals.}
    \item{\code{ss_calibr_jac_function}:}{a function computing the Jacobian of both steady-state (equilibrium) and calibrating functions or NULL.}
    \item{\code{calibr_function}:}{a function returning calibrating equations' residuals.}
    \item{\code{parameters_free_init_val}:}{a vector of free parameters' values declared in the .gcn file.}
    \item{\code{parameters_free_val}:}{a vector of current free parameters' values.}
    \item{\code{parameters_free_mod_flag}:}{a logical vector indicating which free parameters' values have been modified relative to the .gcn file setting.}
    \item{\code{parameters_calibr_val}:}{a vector of current calibrated parameters' values.}

    \item{\code{variables_init_val}:}{a numeric vector of variables' initial values.}
	\item{\code{parameters_calibr_init_val}:}{a numeric vector of calibrated parameters' initial values.}
    \item{\code{init_residual_vector}:}{a numeric vector of residuals of the steady-state (equilibrium) function computed for initial values and calibrated parameters.}
    \item{\code{residual_vector}:}{a numeric vector of residuals of the steady-state (equilibrium) function computed for variables' values and calibrated parameters after the non-linear solver has exited.}
    \item{\code{solver_status}:}{a character string describing the steady-state (equilibrium) solver status.}
    \item{\code{parameters_val}:}{a vector of model parameters' values.}
    \item{\code{variables_ss_val}:}{a vector of variables' steady-state/equilibrium values.}
    \item{\code{ss_solved}:}{logical. If TRUE, the steady state (equilibrium in case of static models) has been found.}

    \item{\code{pert}:}{a function of the first order perturbation (returning a list of matrices).}
    \item{\code{loglin_var}:}{logical. Flags are set to TRUE for log-linearised variables.}

    \item{\code{eig_vals}:}{a matrix of system eigenvalues.}
    \item{\code{solution}:}{a list with elements P, Q, R, S storing solution of the model.}
    \item{\code{state_var_indices}:}{a numeric vector containing the indices of state variables.}
    \item{\code{solver_exit_info}:}{a character string, solver exit information.}
    \item{\code{solution_resid}:}{a list of residuals of perturbation equations, verifying if the model has been solved.}
    \item{\code{re_solved}:}{logical. It is set to TRUE if the model has been solved. The default value is FALSE.}

    \item{\code{active_shocks}:}{a logical vector of the length equal to the number of shocks. If an entry is set to FALSE, the variance of the corresponding shock is zero (the shock is not taken into account during stochastic simulations of the model).}
    \item{\code{shock_cov_mat}:}{a covariance matrix of model shocks.}
    \item{\code{shock_cov_mat_flag}:}{logical. Set to TRUE when the user specifies non-default entries in a covariance matrix of model shocks.}
    \item{\code{corr_mat}:}{a matrix of the model variables' correlations.}
    \item{\code{autocorr_mat}:}{a matrix of the model variables' autocorrelations.}
    \item{\code{ref_var_corr_mat}:}{a matrix of correlations of model variables with the reference variable's lead and lagged values.}
    \item{\code{ref_var_idx}:}{an integer value, the index of the reference variable used in statistics computation.}
    \item{\code{var_dec}:}{a matrix of the variance decomposition of shocks.}
    \item{\code{sdev}:}{a vector of variables' standard deviations.}
    \item{\code{corr_computed}:}{logical. If TRUE, it indicates that the correlations and other variables' statistics have been computed. The default value is FALSE.}
    \item{\code{ext}:}{a list containing additional information about the class instance and potential extensions.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_model")}: prints short diagnostic information about the model.}
    \item{show}{\code{signature(object = "gecon_model")}: prints general information about the model.}
    \item{summary}{\code{signature(object = "gecon_model")}: prints detailed information about the model's computation results.}
    }
}
\seealso{
\code{\link{gecon_model}} is a constructor of objects of \code{gecon_model} class.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")
print(rbc)
class(rbc)
}
