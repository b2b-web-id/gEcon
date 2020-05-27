# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# (c) Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2018-2019                    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################
# Class for storing information about model variables
# ############################################################################

# ############################################################################
# Class definition
# ############################################################################

setClass(
    Class = "gecon_var_info",
    representation = representation(r_object_name = "character",
                                    variables = "character",
                                    is_stochastic = "logical",
                                    is_dynamic = "logical",
                                    ss_solved = "logical",
                                    re_solved = "logical",
                                    corr_computed = "logical",
                                    init_val = "matrix",                             
                                    ss_val = "matrix",
                                    state = "logical",
                                    state_var_impact = "matrix",
                                    shock_impact = "matrix",
                                    std_dev_val = "numeric",
                                    loglin_flag = "logical",
                                    cr = "matrix",
                                    incid_mat = "Matrix")
)

# ############################################################################
# The function gecon_var_info is
# a gecon_var_info class object constructor
# ############################################################################
# Input
#   r_object_name - a character string, the name of an R object of gecon_model
#                   class storing the model from which the information about
#                   variables comes from.
#   variables - a character vector of selected variables' names.
#   is_stochastic - logical. Set to TRUE, if the model has stochastic shocks.
#   is_dynamic - logical. Set to TRUE, if the model has any lead
#                or lagged variables.
#   ss_solved - logical. Set to TRUE, if the steady state/equilibrium of the model
#               has been found.
#   re_solved - logical. Set to TRUE, if the model has been solved.
#   corr_computed - logical. Set to TRUE, if the correlations and other
#                   statistics of variables have been computed.
#   init_val - 	a vector of initial values of variables.                                                 
#   ss_val - a vector of the steady-state (equilibrium) values of selected variables.
#            If the steady state has not been computed, the vector contains
#            last solver iteration values of variables.
#   state - logical. A TRUE entry denotes that a corresponding variable is
#           a state variable.
#   state_var_impact - a matrix with the impact of state variables' values
#                      in the previous period on the variables' values.
#   shock_impact - a matrix with the impact of exogenous variables' values
#                  on the variables' values.
#   std_dev_val - a numeric vector of standard deviations of selected variables.
#   loglin_flag - a logical vector of length equal to the number of selected
#                 variables. A TRUE entry denotes that a corresponding variable
#                 has been log-linearised before solving the model.
#   cr - a matrix with correlations of selected variables with all model variables.
#   incid_mat - a Matrix object with the mapping of variables into model equations
#               and calibrating equations.
# Output
#   An object of the gecon_var_info class.
# ############################################################################
gecon_var_info <- function(r_object_name,
                           variables, is_stochastic,
                           is_dynamic, ss_solved, re_solved,
                           corr_computed, init_val, ss_val,
                           state, state_var_impact, shock_impact,
                           std_dev_val, loglin_flag,
                           cr, incid_mat)
{
    var_info_object <- new("gecon_var_info")

    if (!is.character(r_object_name)) {
        stop("r_object_name should be a character string")
    } else var_info_object@r_object_name <- r_object_name

    if (!is.character(variables)) {
        stop("variables should be of character type")
    } else var_info_object@variables <- variables

    if (!is.logical(is_stochastic)) {
        stop("is_stochastic should be of logical type")
    } else var_info_object@is_stochastic <- is_stochastic

    if (!is.logical(is_dynamic)) {
        stop("is_dynamic should be of logical type")
    } else var_info_object@is_dynamic <- is_dynamic

    if (!is.logical(ss_solved)) {
        stop("ss_solved should be of logical type")
    } else var_info_object@ss_solved <- ss_solved

    if (!is.logical(re_solved)) {
        stop("re_solved should be of logical type")
    } else var_info_object@re_solved <- re_solved

    if (!is.logical(corr_computed)) {
        stop("corr_computed should be of logical type")
    } else var_info_object@corr_computed <- corr_computed

    if (!is.matrix(init_val)) {
        stop("init_val should be of matrix type")
    } else var_info_object@init_val <- init_val

    if (!is.matrix(ss_val)) {
        stop("ss_val should be of matrix type")
    } else var_info_object@ss_val <- ss_val

    if (!is.logical(state)) {
        stop("state should be of logical type")
    } else var_info_object@state <- state

    if (!is.matrix(state_var_impact)) {
        stop("state_var_impact should be of matrix type")
    } else var_info_object@state_var_impact <- state_var_impact

    if (!is.matrix(shock_impact)) {
        stop("shock_impact should be of matrix type")
    } else var_info_object@shock_impact <- shock_impact

    if (!is.numeric(std_dev_val)) {
        stop("std_dev_val should be of numeric type")
    } else var_info_object@std_dev_val <- std_dev_val

    if (!is.logical(loglin_flag)) {
        stop("loglin_flag should be of logical type")
    } else var_info_object@loglin_flag <- loglin_flag

    if (!is.matrix(cr)) {
        stop("cr should be of matrix type")
    } else var_info_object@cr <- cr

    if (!inherits(incid_mat, "Matrix")) {
        stop("incid_mat should be of Matrix class")
    } else var_info_object@incid_mat <- incid_mat

    return (var_info_object)
}


# ############################################################################
# Showing information stored in an object of gecon_var_info class.
# ############################################################################
# Input
#   object - object of gecon_var_info class
# Output
#   Information about variables' incidence and already computed
#   statistics for the variables
# ############################################################################
setMethod("show", signature(object = "gecon_var_info"),
function(object) {
    # Print incidence matrix
    ch_incid_mat <- matrix("      .",
                           nrow = nrow(object@incid_mat),
                           ncol = ncol(object@incid_mat))
    if (object@is_dynamic) {
        ch_incid_mat[which(object@incid_mat == 1)] <-  "    t-1"
        ch_incid_mat[which(object@incid_mat == 2)] <-  "      t"
        ch_incid_mat[which(object@incid_mat == 3)] <-  " t-1, t"
        ch_incid_mat[which(object@incid_mat == 4)] <-  "    t+1"
        ch_incid_mat[which(object@incid_mat == 5)] <-  "t-1, t+1"
        ch_incid_mat[which(object@incid_mat == 6)] <-  " t, t+1"
        ch_incid_mat[which(object@incid_mat == 7)] <-  "t-1, t, t+1"
        ch_incid_mat[which(object@incid_mat == 8)] <-  "     ss"
        ch_incid_mat[which(object@incid_mat == 9)] <-  "ss, t-1"
        ch_incid_mat[which(object@incid_mat == 10)] <- "  ss, t"
        ch_incid_mat[which(object@incid_mat == 11)] <- "ss, t-1, t"
        ch_incid_mat[which(object@incid_mat == 12)] <- "ss, t+1"
        ch_incid_mat[which(object@incid_mat == 13)] <- "ss, t-1, t+1"
        ch_incid_mat[which(object@incid_mat == 14)] <- "ss, t, t+1"
        ch_incid_mat[which(object@incid_mat == 15)] <- "ss, t-1, t, t+1"
    } else {
        ch_incid_mat[which(object@incid_mat != 0)] <-  "      X"
    }
    colnames(ch_incid_mat) <- colnames(object@incid_mat)
    rownames(ch_incid_mat) <- rownames(object@incid_mat)
    cat("Incidence info:\n\n")
    print(as.data.frame(ch_incid_mat))

    # Print info about value
    cat("\n----------------------------------------------------------", "\n")
    if (object@ss_solved) {
        if (object@is_dynamic) {
            cat("\nSteady-state values:\n\n")
        } else {
            cat("\nEquilibrium values:\n\n")
        }
    } else if (!all(is.na(object@ss_val))) {
        if (object@is_dynamic) {
            cat("Last solver iteration values for steady state computation:\n\n")
        } else {
            cat("Last solver iteration values for equilibrium computation:\n\n")
        }
    }
    print(round(object@ss_val, digits = 4))

    cat("\n----------------------------------------------------------", "\n")
    if (object@is_dynamic) {
        cat("\nInitial values for steady state computation:\n\n")
    } else {
        cat("\nInitial values for equilibrium computation:\n\n")
    }
    print(round(object@init_val, digits = 4))
    # Print info about loglinearized and state variables
    if (object@re_solved) {
        ch_state <- matrix("", length(object@variables), 2)
        colnames(ch_state) <- c("Is a state variable?", "Is log-linearised?")
        rownames(ch_state) <- object@variables
        ch_state[which(object@state), 1] <- "Y          "
        ch_state[, 2] <- "N          "
        ch_state[which(object@loglin_flag), 2] <- "Y          "
        cat("\n----------------------------------------------------------", "\n")
        cat("\nVariable info:\n\n")
        print(as.data.frame(ch_state))
    }

    # Print info about the recursive laws of motion
    if (object@re_solved) {
        cat("\n----------------------------------------------------------", "\n")
        cat("\nRecursive laws of motion for the variables\n")

        cat("\nState variables impact:\n\n")
        print(round(object@state_var_impact, digits = 4))

        if (object@is_stochastic) {
            cat("\nShocks impact:\n\n")
            print(round(object@shock_impact, digits = 4))
        }
    }

    # Print info about moments and correlations
    if (object@corr_computed) {
        loglin_indic <- rep("N   ", length(object@variables))
        loglin_indic[object@loglin_flag] <- "Y   "
        mom <- round(cbind(object@ss_val,
                            object@std_dev_val,
                            object@std_dev_val ^ 2), 4)

        rownames(mom) <- object@variables
        mom <- as.data.frame(cbind(mom, loglin_indic))
        colnames(mom) <- c("Steady-state value", " Std. dev.", " Variance", "Loglin")
        cat("\n----------------------------------------------------------", "\n\n")
        cat("Basic statistics:\n\n")
        print(mom)
        cat("\n")

        # Correlations
        if (dim(object@cr)[1]) {
            cat("Correlations:\n\n")
            print(round(object@cr, digits = 4))
            cat("\n")
        }
    }
})

# ############################################################################
# Showing information stored in an object of gecon_var_info class.
# ############################################################################
# Input
#   object - object of class gecon_var_info
# Output
#   Information about variables incidence and already computed
#   statistics for the variables
# ############################################################################
setMethod("print", signature(x = "gecon_var_info"),
function(x, ...)
{
    show(x)
})

# ############################################################################
# Showing information stored in an object of gecon_var_info class.
# ############################################################################
# Input
#   object - object of class gecon_var_info
# Output
#   Information about variables incidence and already computed
#   statistics for the variables
# ############################################################################
setMethod("summary", signature(object = "gecon_var_info"),
function(object)
{
    show(object)
})

