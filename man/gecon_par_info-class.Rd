\name{gecon_par_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_par_info-class}

\title{\code{gecon_par_info} class}
\description{
The class storing information about parameters selected by the user.
}
\section{Slots}{
  \describe{
    \item{\code{r_object_name}:}{a character string. The name of an R object of \code{gecon_model} class storing the model from which the information about parameters comes from.}
    \item{\code{parameters}:}{a character vector of parameters' names.}
    \item{\code{gcn_values}:}{a numeric vector of free parameters' values specified in the .gcn file.}
    \item{\code{current_values}:}{a numeric vector of parameters' current values.}
    \item{\code{calibr_flag}:}{a logical vector of the length equal to the number of the parameters. If TRUE, a corresponding parameter is a calibrated parameter.}
    \item{\code{calibr_init_val}:}{a numeric vector of calibrated parameters' initial values.}                                                                                              
    \item{\code{incid_mat}:}{a object of \code{Matrix} class representing the mapping of parameters to equations and calibrating equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_par_info")}: Prints all the available information (values, types, incidence matrix) about the parameters which were selected when creating a \code{gecon_par_info-class} object.}
    \item{show}{\code{signature(object = "gecon_par_info")}: Prints all the available information (values, types, incidence matrix) about the parameters which were selected when creating a \code{gecon_par_info-class} object.}
    \item{summary}{\code{signature(object = "gecon_par_info")}: Prints all the available information (values, types, incidence matrix) about the parameters which were selected when creating a \code{gecon_par_info-class} object.}
  }
}
\seealso{
\code{\link{par_info}} creates an object of \code{gecon_par_info} class.
\code{\link{gecon_shock_info-class}} and \code{\link{gecon_var_info-class}} are
similar classes storing the information about shocks and variables, respectively.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# compute the steady state and calibrate alpha
rbc <- steady_state(rbc)
rbc_par_info <- par_info(rbc, all = TRUE)
print(rbc_par_info)
summary(rbc_par_info)
show(rbc_par_info)
}
