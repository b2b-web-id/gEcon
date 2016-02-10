# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# The nleqslv_wrap function is a wrapper for the solvers from
# the nleqslv package. It solves the system of nonlinear equations.
# ###################################################################
# Input
#   start_val - (numeric vector) starting values
#   fn - (function) a system of equations to be solved
#   jac - (function) the Jacobian of the system of equations
#   options_list - solver settings. Following options can be set
#                    in the solver:
#                   -> method - a character specifying method of solving
#                               system can be set to "Newton" or "Broyden".
#                               Default option is "Newton"
#                   -> global - a character specifying search strategy 
#                               can be set to "dbldog", "pwldog", 
#                               "qline", "gline", "none".
#                               Default option is "qline"
#                   -> xscalm - a character specifying method of scaling x
#                               can be set to "fixed", "auto".
#                               Default option is "fixed"
#                   -> max_iter - a numeric, max. number of 
#                               iterations default option is 150
#                   -> tol - a numeric, the tolerance of solution,
#                            default is 1e-6
# Output
#   a list with the following fields:
#       value - the solution
#       function_value - value of function evaluated at the solution
#       exit_code - termination code
# ###################################################################
nleqslv_wrap <- function(start_val, 
                         fn, 
                         jac = NULL,
                         options_list = NULL)
{
    # Default options
    def_options = list(method = "Newton",
                       global = "qline",
                       xscalm = "fixed",
                       max_iter =  150,
                       xtol = 1e-6,
                       tol = 1e-6)
    
    # Checking if specified options are correct
    if (!is.null(options_list)) {
        fitting_options <- which(names(options_list) %in% names(def_options))
        if (length(fitting_options) < length(options_list)) {
            warning('Unrecognised options were specified for nleqslv')
        }
        
        # Overwriting default options
        mo <- match(names(options_list[fitting_options]), names(def_options))
        def_options[mo] <- options_list[fitting_options]
    }

    
    # Calling nleqslv
    err_info <- tryCatch(output <- nleqslv(x = start_val, 
                                           fn = fn, 
                                           jac = jac, 
                                           method = def_options$method, 
                                           global = def_options$global, 
                                           xscalm = def_options$xscalm, 
                                           jacobian = def_options$jacobian, 
                                           control = list(xtol = def_options$xtol,
                                                          ftol = def_options$tol,
                                                          maxit = def_options$max_iter)),
                                           warning = function(w) w,
                                           error = function(w) w)
        

    if(inherits(err_info, 'warning') || inherits(err_info, 'error')) {
        err_info <- gsub(pattern=".* in ", replacement='',
                           x= as.character(err_info), perl = TRUE)
        return(list(value = start_val, 
                    blank_field = NULL, 
                    exit_code = err_info))
    }
    
    return(output)
}
