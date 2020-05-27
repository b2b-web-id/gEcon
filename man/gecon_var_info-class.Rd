\name{gecon_var_info-class}
\Rdversion{1.1}
\docType{class}
\alias{gecon_var_info-class}

\title{\code{gecon_var_info} class}
\description{
The class storing information about variables selected by the user.
}
\section{Slots}{
  \describe{
    \item{\code{r_object_name}:}{a character string. The name of an R object of \code{gecon_model} class storing the model from which the information about variables comes from.}
    \item{\code{variables}:}{a character vector of selected variables' names.}
    \item{\code{is_stochastic}:}{logical. Set to TRUE, if the model, based on which the info was generated, has stochastic shocks.}
    \item{\code{is_dynamic}:}{logical. Set to TRUE, if the model, based on which the info was generated, has any lead or lagged variables.}
    \item{\code{ss_solved}:}{logical. Set to TRUE, if the steady state/equilibrium of the model has been found.}
    \item{\code{re_solved}:}{logical. Set to TRUE, if the model, based on which the info was generated, has been solved.}
    \item{\code{corr_computed}:}{logical. Set to TRUE, if correlations and other statistics of variables have been computed.}
 	\item{\code{init_val}:}{a vector of selected variables' initial values.}                                                                        
    \item{\code{ss_val}:}{a vector of the steady-state (dynamic models) or equilibrium (static models) values of selected variables. If the steady state has not been computed, the vector contains last solver iteration values of variables.}
    \item{\code{state}:}{a logical vector of length equal to the number of selected variables. A TRUE entry denotes that a corresponding variable is a state variable.}
    \item{\code{state_var_impact}:}{the rows of the matrices P and R of state space representation corresponding to selected variables.}
    \item{\code{shock_impact}:}{the rows of the matrices Q and S of state space representation corresponding to selected variables.}
    \item{\code{std_dev_val}:}{a numeric vector of standard deviations of selected variables.}
    \item{\code{loglin_flag}:}{a logical vector of length equal to the number of selected variables. A TRUE entry denotes that a corresponding variable was log-linearised before solving the model.}
    \item{\code{cr}:}{a matrix containing the correlations of selected variables with all model variables.}
    \item{\code{incid_mat}:}{an object of \code{Matrix} class representing the mapping of variables to model equations and calibrating equations.}
}
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "gecon_var_info")}: Prints all the available information (values, statistics, incidence matrix, etc.) about the variables, which were selected when creating a \code{gecon_var_info-class} object.}
    \item{show}{\code{signature(object = "gecon_var_info")}: Prints all the available information (values, statistics, incidence matrix, etc.) about the variables, which were selected when creating a \code{gecon_var_info-class} object.}
    \item{summary}{\code{signature(object = "gecon_var_info")}: Prints all the available information (values, statistics, incidence matrix, etc.) about the variables, which were selected when creating a \code{gecon_var_info-class} object.}
 }
}
\seealso{
\code{\link{var_info}} creates an object of \code{gecon_var_info} class.
\code{\link{gecon_shock_info-class}} and \code{\link{gecon_par_info-class}} are
similar classes storing the information about shocks and parameters, respectively.
}
\examples{
# copy the example to the current working directory
file.copy(from = file.path(system.file("examples", package = "gEcon"),
          "rbc.gcn"), to = getwd())

# make and load the model
rbc <- make_model("rbc.gcn")

# compute the steady state
rbc <- steady_state(rbc)
rbc <- solve_pert(rbc)
rbc_var_info <- var_info(rbc, all = TRUE)
print(rbc_var_info)
summary(rbc_var_info)
show(rbc_var_info)
}
