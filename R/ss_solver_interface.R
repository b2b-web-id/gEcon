# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# The ss_par_compute function returns steady-state values of variables  
# in case of dynamic model (equilibrium in case of static model)
# and calibrated parameters
# ###################################################################
#
# Input  
#   f_ssequations <- a function returning residuals for the    
#                conditions characterizing the equilibrium of the model
#   ss_guess <- a vector of estimated steady-state
#               (equilibrium) values of variables 
#   param_guess <- a vector of estimated parameter values
#   calib_par_l <- the number of calibrated parameters
#   calibration <- if TRUE, the calibrating equations are taken
#                  into account while solving for steady state
#                  (equilibrium) and parameters are calibrated. Otherwise
#                  calibrating equations are dropped.
#   ss_function_jac <- the Jacobian of the steady state
#                      (equilibrium) function 
#   use_jac <- a logical, defines whether Jacobian computed 
#               by symbolic library or numerical 
#               approximation should be used
#   solver <- a character string naming a chosen solver of a system 
#             of non linear equations defined by ssequations (by 
#             default "slv1_nleqslv"); solvers available: 
#             "slv1_nleqslv" - solver from nleqslv package
#   options_list - a list of settings. Solver specific
#                   check in nleqslv_wrap and model class
#                   for detailed description of options
# Output
#   a list containing:
#       a vector of steady-state (equilibrium) values of variables, 
#       a vector of parameter values, 
#       a vector of residuals,
#   (optionally) a list with computation's results
# ###################################################################
ss_par_compute <- function(f_ssequations, 
                           ss_guess, 
                           param_guess, 
                           calib_par_l,
                           calibration,
                           ss_function_jac,
                           use_jac = FALSE,
                           solver = "slv1_nleqslv",
                           options_list = NULL)
{
    nvars <- length(ss_guess)
    npar <- length(param_guess)
    if (calibration) {
        guess <- c(ss_guess, param_guess)
    } else {
        guess <- ss_guess
    }
    # wrapping function
    if (calibration) {
        steady_state <- function(yy) 
        {
            ss <- yy[1:nvars]
            if (calib_par_l) {
                ppp <- yy[(nvars + 1):(nvars + npar)]
            } else ppp <- NULL
            f_ssequations(ss, ppp)
        }
    } else steady_state <- f_ssequations
    
    
    # wrapping Jacobian if uses Jacobian supplied by user
    if (use_jac) {
        if (calibration) {
            jac <- function(yy) 
            {
                ss <- yy[1:nvars]
                pp <- yy[(nvars + 1):(nvars + npar)]
                if (npar) {
                    a <- as.matrix(ss_function_jac(ss, pp))
                } else a <- as.matrix(ss_function_jac(ss, NULL))
                return(a)
            }
        } else jac <- function(yy) 
            {
                as.matrix(ss_function_jac(yy, NULL))
            }   
    } else jac <- NULL
    
    solver_list <- list("slv1_nleqslv")       
    pos <- as.numeric(match(solver, solver_list, nomatch = 0))
    
    # solver not specified in solver list
    if (pos == 0) {
        stop(paste(paste("A solver named", "'", sep=" "), 
                solver, paste("'", "is not available.", sep=" "), sep=""))
        # using nleqslv solver
        print(options_list)
    } else if (pos == 1) {
        if (use_jac) {
            solutions = nleqslv_wrap(start_val = guess, 
                                     fn = steady_state, 
                                     jac = jac, 
                                     options_list = options_list)
        } else {
            solutions = nleqslv_wrap(start_val = guess, 
                                     fn = steady_state, 
                                     jac = NULL, 
                                     options_list = options_list)
        }

        solution <- solutions[[1]]
        if (is.character(solutions[[3]])) {
            descr <- solutions[[3]]
        } else {
            result <- as.numeric(solutions[[3]])
            descr <- switch(result,
                            "Convergence of function values has been achieved.",
                            "x-values within tolerance.",
                            "The algorithm has stalled and cannot find an acceptable new point.",
                            "Iteration limit maxit exceeded.",
                            "Jacobian is too ill-conditioned.",
                            "Jacobian is singular.",
                            "Analytical Jacobian is most likely incorrect.")
        }
        if (npar) {
            results <- list(ss_values = solution[1:nvars], 
                        par_values = solution[(nvars + 1):(nvars + npar)], 
                        solver_status = descr)
        } else {
            results <- list(ss_values = solution[1:nvars], 
                        par_values = NULL, 
                        solver_status = descr)            
        }         
        return(results)
    }
}
