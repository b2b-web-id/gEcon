\name{gecon_model}
\alias{gecon_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create objects of \code{gecon_model} class
}
\description{
The \code{gecon_model} function is a constructor of objects of \code{gecon_model} class.
}
\usage{
gecon_model(model_info, index_sets,
            variables, variables_tex,
            shocks, shocks_tex,
            parameters, parameters_tex,
            parameters_free, parameters_free_val,
            equations, calibr_equations,
            var_eq_map, shock_eq_map, var_ceq_map, cpar_eq_map,
            cpar_ceq_map, fpar_eq_map, fpar_ceq_map,
            ss_function, calibr_function, ss_calibr_jac_function,
            pert, ext)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_info}{
a character vector of length 4, containing information about the model: the input file name, the input file path, the date of creation, and a logical flag informing if the model requires compiling.
}
  \item{index_sets}{
a list containing information about index sets. The names of its elements correspond to sets' names.
Each element contains a character vector of the names of the relevant set's components.
}
  \item{variables}{
a character vector of variables' names.
}
  \item{variables_tex}{
a character vector of variables' LaTeX names.
}
  \item{shocks}{
a character vector of shocks' names.
}
  \item{shocks_tex}{
a character vector of shocks' LaTeX names.
}
  \item{parameters}{
a character vector of all parameters' names.
}
  \item{parameters_tex}{
a character vector of all parameters' LaTeX names.
}
  \item{parameters_free}{
a character vector of free parameters' names.
}
  \item{parameters_free_val}{
a vector of free parameters' values.
}
  \item{equations}{
a character vector of model equations.
}
  \item{calibr_equations}{
a character vector of calibrating equations.
}
  \item{var_eq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of variables to equations.
}
  \item{shock_eq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of shocks to equations.
}
  \item{var_ceq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of variables to calibrating equations.
}
  \item{cpar_eq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of calibrated parameters to equations.
}
  \item{cpar_ceq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of calibrated parameters to calibrating equations.
}
  \item{fpar_eq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of free parameters to equations.
}
  \item{fpar_ceq_map}{
a sparse matrix (an object of \code{Matrix} class) representing the mapping of free parameters to calibrating equations.
}
  \item{ss_function}{
a function returning the steady-state/equilibrium equations' residuals.
}
  \item{calibr_function}{
a function used for the calibration of parameters.
}
  \item{ss_calibr_jac_function}{
a function computing the Jacobian of both steady-state (equilibrium) and calibrating functions or NULL.
}
  \item{pert}{
a function returning matrices representing the first order perturbation of the model
(in the canonical form).}
  \item{ext}{
a list containing additional information about the class instance and potential extensions.}
}
\value{
An object of \code{gecon_model} class.
}
\note{
The \code{gecon_model} constructor is used in .R files created by \code{gEcon}. In general, users should not call this function explicitly.
}
