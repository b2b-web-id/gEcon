# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################

# ############################################################################
# The ss_find function returns steady-state (equilibrium) values of variables
# and calibrated parameters in a dynamic (static) model
# ############################################################################
# Input
#   ss_guess - a vector of estimated steady-state (equilibrium) values of
#              variables and calibrated parameters
#   ss_eval - a function returning residuals for the conditions characterizing
#             the steady state (equilibrium) of a model
#   ss_jac - the Jacobian of the steady state (equilibrium) function
#   solver - the name of non-linear equations solver (by default "nleqslv")
#   options_list - a list of nleqslv solver specific settings; a detailed
#                  description of options can be found in the documentation
#                  of the nleqslv_wrap function
# Output
#   a list containing:
#       a vector of steady-state (equilibrium) values of variables and parameters
#       solver status' description
# ############################################################################
ss_find <- function(ss_guess, ss_eval, ss_jac, solver, options_list)
{
    if (!is.character(solver) || length(solver) != 1) {
        stop("invalid solver argument")
    }
    solver_list <- list("nleqslv")

    pos <- as.numeric(match(solver, solver_list, nomatch = 0))

    # solver not specified in solver list
    if (pos == 0) {
        stop(paste0("solver called \'", solver, "\' is not supported"))
    }
    # nleqslv
    if (pos == 1) {
        solutions <- nleqslv_wrap(init_val = ss_guess,
                                  fn = ss_eval,
                                  jac = ss_jac,
                                  options_list = options_list)
        solution <- solutions[[1]]
        if (is.character(solutions[[3]])) {
            descr <- solutions[[3]]
        } else {
            result <- as.numeric(solutions[[3]])
            descr <- switch(result,
                            "Convergence of function values has been achieved.",
                            "x-values within tolerance.",
                            "The algorithm has stalled and cannot find an acceptable new point.",
                            "Iteration limit max_iter exceeded.",
                            "Jacobian is too ill-conditioned.",
                            "Jacobian is singular.",
                            "Analytical Jacobian is most likely incorrect.")
        }
        return(list(solution = solution, solver_status = descr))
    }
}
