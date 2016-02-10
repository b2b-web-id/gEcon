\name{gecon_model}
\alias{gecon_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructor of "gecon_model" class
}
\description{
The \code{gecon_model} function is a constructor of \code{gecon_model} class objects.
}
\usage{
gecon_model(model_info, index_sets, variables, shocks, parameters, parameters_free,
            parameters_free_val, equations, calibr_equations, 
            var_eq_map, shock_eq_map, var_ceq_map, cpar_eq_map,
            cpar_ceq_map, fpar_eq_map, fpar_ceq_map, ss_function, 
            calibr_function, ss_calibr_function_jac, pert)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_info}{
a character vector of length 3,  containing information about the model: the input file name, the input file path and the date of creation.
}
  \item{index_sets}{
a list containing information about index sets. The names of the list components correspond to the set names.
Each component contains character vector of the names of the relevant set elements.
}
  \item{variables}{
a character vector of all the variable names.
}
  \item{shocks}{
a character vector of all the shock names.
}
  \item{parameters}{
a character vector of all the parameter names.
}
  \item{parameters_free}{
a character vector of all the free parameter names.
}
  \item{parameters_free_val}{
a numeric vector of the values of all the free parameters.
}
  \item{equations}{
a character vector of model equations.
}
  \item{calibr_equations}{
a character vector of model calibrating equations.
}
  \item{var_eq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of variables to equations.
}
  \item{shock_eq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of shocks to equations.
}
  \item{var_ceq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of variables to calibrating equations.
}
  \item{cpar_eq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of calibrated parameters to equations.
}
  \item{cpar_ceq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of calibrated parameters to calibrating equations.
}
  \item{fpar_eq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of free parameters to equations.
}
  \item{fpar_ceq_map}{
a sparse matrix (a \code{Matrix} class object) representing the mapping of free parameters to calibrating equations.
}
  \item{ss_function}{
a function returning residuals from the steady-state (equilibrium for the static models) equations.
}
  \item{calibr_function}{
a function used for the calibration of parameters.
}
  \item{ss_calibr_function_jac}{
a function returning a Jacobian of functions returning residuals from the steady-state (equilibrium for the static models) equations.
}
  \item{pert}{
a function returning a list with the matrices representing canonical form of the model.}
}
\value{
An object of the \code{gecon_model} class.
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}

\note{
The \code{gecon_model} constructor is used in R files created by \code{gEcon}.
}
