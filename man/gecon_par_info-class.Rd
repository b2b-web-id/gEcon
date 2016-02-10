\name{gecon_par_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_par_info-class}

\title{Class \code{"gecon_par_info"}}
\description{
The class storing information about parameters chosen by the user.
}
\section{Objects from the Class}{
Objects of this class are created when calling the \code{\link{par_info}} function.
}
\section{Slots}{
  \describe{
    \item{\code{model_info}:}{a character vector of length 3, containing information about the model: the input file name, the input file path, and the date of creation.}
    \item{\code{model_variable_name}:}{a character string denoting the name of a variable storing the model for which the information about parameters has been created.}
    \item{\code{par_names}:}{a character vector of parameter names.}
    \item{\code{gcn_values}:}{a numeric vector of the values of free parameters specified in the .gcn file.}
    \item{\code{current_values}:}{a numeric vector of the most recent values of the parameters.}
    \item{\code{calibr_flag}:}{a logical vector of the length equal to the number of the parameters. The TRUE entries denote that a corresponding parameter is a calibrated parameter.}
    \item{\code{incid_mat}:}{a \code{Matrix} object representing the mapping of parameters to equations and calibrating equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_par_info")}: Prints all the available information (short listing, values, type, incidence) about the parameters retrieved from the model when creating a \code{gecon_par_info-class} object. }
    \item{show}{\code{signature(object = "gecon_par_info")}: Prints information about parameters' types, values, and the incidence matrix.}
    \item{summary}{\code{signature(object = "gecon_par_info")}: Prints all the available information (short listing, values, type, incidence) about the parameters, retrieved from the model when creating a \code{gecon_par_info-class} object. }
  }
}
\author{
Karol Podemski \email{karol.podemski@gmail.com}, \cr
Kaja Retkiewicz-Wijtiwiak \email{kaja.retkiewicz@gmail.com}
}
\seealso{
\code{\link{par_info}} to create a \code{gecon_par_info} object. The analogous classes storing the information about
shocks and variables are \code{\link{gecon_shock_info-class}} and \code{\link{gecon_var_info-class}}.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = 'gEcon'),
          'rbc.gcn'), to = getwd())

# make and load the model
rbc <- make_model('rbc.gcn')

# compute the steady state and calibrate alpha
rbc <- steady_state(rbc)
rbc_par_info <- par_info(rbc, all_parameters = TRUE)
print(rbc_par_info)
summary(rbc_par_info)
show(rbc_par_info)
}
