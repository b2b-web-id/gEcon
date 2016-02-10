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
This method shows objects of the \code{gecon_simulation}, 
 \code{gecon_model}, \code{gecon_var_info}, \code{gecon_shock_info},
and \code{gecon_par_info} classes.}
\section{Methods}{
\describe{

\item{\code{signature(object = "gecon_simulation")}}{
Shows the name of the model for which the simulations have been created, information about time span, shocks, and variables used.
}

\item{\code{signature(object = "gecon_model")}}{
Shows the type of the model, the date of creation, the solving status, and the information about number of variables and parameters.
}

\item{\code{signature(object = "gecon_var_info")}}{Prints information about variables' incidence and the results that have been already obtained for the variables. }

\item{\code{signature(object = "gecon_shock_info")}}{Prints the incidence matrix and the variance-covariance matrix of shocks, retrieved from the model when creating \code{gecon_shock_info} object. }

\item{\code{signature(object = "gecon_par_info")}}{Prints information about parameters' type, value, and the incidence matrix. }

}}
\keyword{methods}
