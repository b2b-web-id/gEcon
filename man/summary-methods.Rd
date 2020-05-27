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
This method summarizes information about objects of \code{gecon_simulation},
\code{gecon_model}, \code{gecon_var_info}, \code{gecon_shock_info}, and
\code{gecon_par_info} classes.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "gecon_simulation")}}{
Prints a summary of an object of \code{gecon_simulation} class consisting of shock covariance matrix and a simulation for each shock.
}

\item{\code{signature(object = "gecon_model")}}{
Prints a summary of an object of \code{gecon_model} class consisting of computed statistics and values.
}

\item{\code{signature(object = "gecon_var_info")}}{Prints the incidence matrix and the results that have been already obtained in terms of model variables. }

\item{\code{signature(object = "gecon_shock_info")}}{Prints the incidence matrix and the covariance matrix of shocks.}

\item{\code{signature(object = "gecon_par_info")}}{Prints the incidence matrix and information about parameters' types and values.}

}}
\keyword{methods}
