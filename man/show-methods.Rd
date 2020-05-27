\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,gecon_simulation-method}
\alias{show,gecon_model-method}
\alias{show,gecon_par_info-method}
\alias{show,gecon_shock_info-method}
\alias{show,gecon_var_info-method}

\title{Show methods for classes in the gEcon package}
\description{
This method prints general information about objects of the \code{gecon_simulation},
 \code{gecon_model}, \code{gecon_var_info}, \code{gecon_shock_info},
and \code{gecon_par_info} classes.}
\section{Methods}{
\describe{

\item{\code{signature(object = "gecon_simulation")}}{
Prints the name of the model for which the simulations have been performed, information about time span, shocks, and variables used.
}

\item{\code{signature(object = "gecon_model")}}{
Prints a short description of the model, its creation date, model's solution status, and the information about the numbers of model variables and parameters.
}

\item{\code{signature(object = "gecon_var_info")}}{Prints the incidence matrix and the results that have been already obtained in terms of model variables. }

\item{\code{signature(object = "gecon_shock_info")}}{Prints the incidence matrix and the covariance matrix of shocks.}

\item{\code{signature(object = "gecon_par_info")}}{Prints the incidence matrix and information about parameters' types and values.}

}}
\keyword{methods}
