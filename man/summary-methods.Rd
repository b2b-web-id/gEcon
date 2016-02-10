\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,gecon_simulation-method}
\alias{summary,gecon_model-method}
\alias{summary,gecon_par_info-method}
\alias{summary,gecon_shock_info-method}
\alias{summary,gecon_var_info-method}

\title{Summary methods for classes in the gEcon package} 
\description{
This method summarizes the objects of the \code{gecon_simulation}, 
\code{gecon_model}, \code{gecon_var_info}, \code{gecon_shock_info}, and
\code{gecon_par_info} classes.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "gecon_simulation")}}{
Prints a summary of a \code{gecon_simulation} class object consisting of a shock matrix and the simulation for each shock.
}

\item{\code{signature(object = "gecon_model")}}{
Prints a summary of a \code{gecon_model} class object consisting of all the computed statistics and values.
}

\item{\code{signature(object = "gecon_var_info")}}{Prints all the available information (short listing, the incidence matrix, the statistics) about the variables, retrieved from the model when creating a \code{gecon_var_info-class} object.}

\item{\code{signature(object = "gecon_shock_info")}}{Prints all the available information (short listing, the incidence matrix , the variance-covariance matrix) about the shocks, retrieved from the model when creating a \code{gecon_shock_info-class} object.}

\item{\code{signature(object = "gecon_par_info")}}{Prints all the available information (short listing, the values, the type, and the incidence) about the parameters, retrieved from the model when creating a \code{gecon_par_info-class} object.}

}}
\keyword{methods}
