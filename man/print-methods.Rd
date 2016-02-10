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
Prints the name of the model for which the simulations have been generated, information about the time span, shock, and variables used.
}

\item{\code{signature(x = "gecon_model")}}{
Shows the type of the model, the date of creation, the solving status, and more detailed information about the number of variables and parameters then the \code{show} generic.
}

\item{\code{signature(x = "gecon_var_info")}}{Prints all the available information (a short listing, the incidence matrix, the statistics) about the variables, retrieved from the model when creating a \code{gecon_var_info-class} object.}

\item{\code{signature(x = "gecon_shock_info")}}{Prints all the available information (a short listing, the incidence matrix , the variance-covariance matrix) about the shocks, retrieved from the model when creating a \code{gecon_shock_info-class} object.}

\item{\code{signature(x = "gecon_par_info")}}{Prints all the available information (a short listing, the values, the type, and the incidence) about the parameters, retrieved from the model when creating a \code{gecon_par_info-class} object.}

}}
\keyword{methods}
