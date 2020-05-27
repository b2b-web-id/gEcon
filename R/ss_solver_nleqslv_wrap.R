# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################

# ############################################################################
# The nleqslv_wrap function is a wrapper for non-linear equations solvers
# from the nleqslv package.
# ############################################################################
# Input
#   init_val - (numeric vector) initial values
#   fn - (function) a system of equations to be solved
#   jac - (function) the Jacobian of the steady-state equations
#   options_list - nleqslv solver specific settings. The following options can be set:
#                  -> method - a character string specifying the method used
#                              for finding a solution; it can be set to "Newton"
#                              or "Broyden"; default option is "Newton"
#                  -> global - a character string specifying global search strategy
#                              to use; it can be set to "dbldog", "pwldog", "qline",
#                              "gline", "none"; the default option is "qline"
#                  -> xscalm - a character string specifying the type of init_val
#                              scaling to be applied; it can be set to: "fixed", "auto";
#                              the default option is "fixed"
#                  -> max_iter - a numeric value, the maximal number of iterations;
#                                the default value is 150
#                  -> tol - a numeric value setting a numeric tolerance for a solution,
#                           the default value is 1e-6
# Output
#   a list containing the following components:
#       value - the solution
#       function_value - value of function evaluated at the solution
#       exit_code - termination code as integer
# ############################################################################
nleqslv_wrap <- function(init_val, fn, jac = NULL, options_list = NULL)
{
    # Default options
    def_options <- list(method = "Newton",
                        global = "qline",
                        xscalm = "fixed",
                        max_iter = 150,
                        xtol = 1e-6,
                        tol = 1e-6)

    # Checking if specified options are correct
    if (!is.null(options_list)) {
        fitting_options <- which(names(options_list) %in% names(def_options))
        if (length(fitting_options) < length(options_list)) {
            warning("unrecognised options were specified for \'nleqslv\'")
        }
        # Overwriting default options
        mo <- match(names(options_list[fitting_options]), names(def_options))
        def_options[mo] <- options_list[fitting_options]
    }

    # nleqslv expects dense matrices
    if (!is.null(jac)) {
        jacob <- function (x) return (as.matrix(jac(x)))
    } else {
        jacob <- jac
    }

    # Calling nleqslv function
    err_info <- tryCatch(output <- nleqslv(x = init_val,
                                           fn = fn,
                                           jac = jacob,
                                           method = def_options$method,
                                           global = def_options$global,
                                           xscalm = def_options$xscalm,
                                           jacobian = def_options$jacobian,
                                           control = list(xtol = def_options$xtol,
                                                          ftol = def_options$tol,
                                                          maxit = def_options$max_iter)),
                                           warning = function(w) w,
                                           error = function(w) w)

    if(inherits(err_info, "warning") || inherits(err_info, "error")) {
        err_info <- gsub(pattern = ".* in ", replacement = "",
                         x = as.character(err_info), perl = TRUE)
        return (list(value = init_val,
                    blank_field = NULL,
                    exit_code = err_info))
    }

    return (output)
}
