\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,gecon_simulation-method}
\alias{print,gecon_model-method}
\alias{print,gecon_var_info-method}
\alias{print,gecon_shock_info-method}
\alias{print,gecon_par_info-method}
\title{Print methods for classes in the gEcon package}
\description{
Prints information about objects of the \code{gecon_simulation}, \code{gecon_model},
\code{gecon_var_info}, \code{gecon_shock_info}, and \code{gecon_par_info} classes.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "gecon_simulation")}}{
Prints the name of the model for which the simulations have been performed, information about the time span, shocks, and variables used.
}

\item{\code{signature(x = "gecon_model")}}{
Prints a short description of the model, its creation date, model's solution status, and more detailed information concerning model variables and parameters then the \code{show} generic function.
}

\item{\code{signature(object = "gecon_var_info")}}{Prints the incidence matrix and the results that have been already obtained in terms of model variables. }

\item{\code{signature(object = "gecon_shock_info")}}{Prints the incidence matrix and the covariance matrix of shocks.}

\item{\code{signature(object = "gecon_par_info")}}{Prints the incidence matrix and information about parameters' types and values.}

}}
\keyword{methods}
