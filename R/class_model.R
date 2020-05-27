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
# Class representing general equilibrium model
# ############################################################################


# ############################################################################
# Class to store Jacobian function or NULL value
# ############################################################################
setClassUnion("function_or_null", c("function", "NULL"))


# ############################################################################
# Class definition
# ############################################################################
setClass(
    Class = "gecon_model",
    representation = representation(
        # Info about the model
        model_info = "character",

        # Index sets
        index_sets = "list",

        # Parameters, variables, shocks
        parameters = "character",
        parameters_tex = "character",
        parameters_free = "character",
        map_free_into_params = "numeric",
        parameters_calibr = "character",
        map_calibr_into_params = "numeric",
        variables = "character",
        variables_tex = "character",
        shocks = "character",
        shocks_tex = "character",

        # Equations
        equations = "character",
        calibr_equations = "character",
        var_eq_map = "Matrix",
        shock_eq_map = "Matrix",
        var_ceq_map = "Matrix",
        cpar_eq_map = "Matrix",
        cpar_ceq_map = "Matrix",
        fpar_eq_map = "Matrix",
        fpar_ceq_map = "Matrix",

        # Model type
        is_stochastic = "logical",
        is_dynamic = "logical",
        is_calibrated = "logical",

        # Steady-state / equilibrium
        ss_function = "function",
        calibr_function = "function",
        ss_calibr_jac_function = "function_or_null",
        parameters_free_init_val = "vector",
        parameters_free_val = "vector",
        parameters_free_mod_flag = "logical",
        parameters_calibr_val = "numeric",

        # Steady-state / equilibrium solution
        variables_init_val = "numeric",
        parameters_calibr_init_val = "numeric",
        init_residual_vector = "numeric",
        residual_vector = "numeric",
        solver_status = "character",
        parameters_val = "numeric",
        variables_ss_val = "numeric",
        ss_solved = "logical",

        # 1st order perturbation
        pert = "function",
        loglin_var = "logical",

        # Solution to 1st order perturbation
        eig_vals = "matrix",
        solution = "list",
        state_var_indices = "numeric",
        solver_exit_info = "character",
        solution_resid = "list",
        re_solved = "logical",

        # Stochastic simulation
        active_shocks = "logical",
        shock_cov_mat = "matrix",
        shock_cov_mat_flag = "logical",
        corr_mat = "matrix",
        autocorr_mat = "matrix",
        ref_var_corr_mat = "matrix",
        ref_var_idx = "integer",
        var_dec = "matrix",
        sdev = "matrix",
        corr_computed = "logical",
        
        # Potential extensions
        ext = "list"
    ) ,
    prototype = prototype(
        parameters = character(0),
        variables = character(0),
        shocks = character(0),
        ref_var_idx = 0L,
        solution_resid = list(),
        ss_solved = FALSE,
        init_residual_vector = numeric(0),
        residual_vector = numeric(0),
        pert = function(x) NULL,
        corr_computed = FALSE,
        loglin_var = logical(0),
        is_calibrated = TRUE,
        re_solved = FALSE,
        is_dynamic = FALSE,
        is_stochastic = TRUE,
        shock_cov_mat_flag = FALSE,
        solution = list(P = NULL,
                        Q = NULL,
                        R = NULL,
                        S = NULL),
        ext = list()
    )
)

# ############################################################################
# The gecon_model function is a constructor of the class gecon_model
# ############################################################################
# Input
#   model_info - [character vector, length = 4] information about
#                   model: input file name, file path, the date of creation
#                   and a flag informing if the model contains compiled C++ code
#   index_sets - [list] each of the list components corresponds to one
#                 set specified in gecon model class. Each components stores
#                 all the elements in the set.
#   variables - [character vector] of all variable names
#   variables_tex - [character vector] of all variable LaTeX names
#   shocks - [character vector]  of all shock names
#   shocks_tex - [character vector]  of all shock LaTeX names
#   parameters - [character vector] of all parameter names
#   parameters_tex - [character vector] of all parameter LaTeX names
#   parameters_free - [character vector] of all free parameter names
#   parameters_free_val - [numeric vector] values of free parameters
#   equations - [character vector] of all the model equations
#   calibr_equations - [character vector] of all calibrating equations
#   var_eq_map - [(sparse) Matrix class] the mapping of variables
#                to equations
#   shock_eq_map - [(sparse) Matrix class] the mapping of shocks
#                  to equations
#   var_ceq_map - [sparse Matrix class] the mapping of variables
#                  to calibrating equations.
#   cpar_eq_map - [sparse Matrix class] the mapping of calibrated
#                  parameters to equations.
#   cpar_ceq_map - [sparse Matrix class] the mapping of calibrated
#                   parameters to calibrating equations.
#   fpar_eq_map - [sparse Matrix class] the mapping of free
#                  parameters to equations.
#   fpar_ceq_map - [sparse Matrix class] the mapping of free
#                   parameters to calibrating equations.
#   ss_function - [function] function returning residuals of steady state equations
#   calibr_function - [function] a function used for calibration of variables
#   ss_calibr_jac_function - [function or NULL] Jacobian of system
#                             of functions defining steady state
#                             a dynamic model or equilibrium
#                             in a static model
#                             and calibration equations
#   pert - [function] the function returning perturbation matrices.
#   ext - [list] a list containing additional information about the class instance
#          and potential extensions. 
#
# Output
#   An object of class "gecon_model"
# ############################################################################
gecon_model <- function(model_info,
                        index_sets,
                        variables,
                        variables_tex,
                        shocks,
                        shocks_tex,
                        parameters,
                        parameters_tex,
                        parameters_free,
                        parameters_free_val,
                        equations,
                        calibr_equations,
                        var_eq_map,
                        shock_eq_map,
                        var_ceq_map,
                        cpar_eq_map,
                        cpar_ceq_map,
                        fpar_eq_map,
                        fpar_ceq_map,
                        ss_function,
                        calibr_function,
                        ss_calibr_jac_function,
                        pert,
                        ext)
{
    mod <- new("gecon_model")

    if (!is.character(model_info)) {
        stop("model_info should be of character type")
    } else {
        mod@model_info <- model_info
        names(mod@model_info) <- c("Model name", "Source .gcn file",
                                   "Generated", "Is compiled")
    }

    if (!is.list(index_sets)) {
        stop("index_sets must be a list of sets (character vectors)")
    } else mod@index_sets <- index_sets

    if (!is.character(variables)) {
        stop("variables should be of character type")
    } else mod@variables <- variables

    if (!is.character(variables_tex)) {
        stop("variables_tex should be of character type")
    } else mod@variables_tex <- variables_tex

    if (!is.character(shocks)) {
        stop("shocks should be of character type")
    } else mod@shocks <- shocks

    if (!is.character(shocks_tex)) {
        stop("shocks_tex should be of character type")
    } else mod@shocks_tex <- shocks_tex

    if (!is.character(parameters)) {
        stop("parameters should be of character type")
    } else mod@parameters <- parameters

    if (!is.character(parameters_tex)) {
        stop("parameters_tex should be of character type")
    } else mod@parameters_tex <- parameters_tex

    if (!is.character(parameters_free)) {
        stop("parameters_free should be of character type")
    } else mod@parameters_free <- parameters_free

    mod@parameters_calibr <-
        parameters[-which(parameters %in% parameters_free)]

    if (!is.logical(parameters_free_val) & !is.numeric(parameters_free_val)) {
        stop("parameters_free_val should be a numeric or logical vector")
    } else {
        mod@parameters_free_val <- as.numeric(parameters_free_val)
        mod@parameters_free_init_val <- as.numeric(parameters_free_val)
    }

    if (!is.character(equations)) {
        stop("equations should be of character type")
    } else mod@equations <- equations

    if (!is.character(calibr_equations)) {
        stop("calibr_equations should be of character type")
    } else mod@calibr_equations <- calibr_equations

    if (!inherits(var_eq_map, "Matrix")) {
        stop("var_eq_map should be of Matrix class")
    } else mod@var_eq_map <- var_eq_map

    if (!inherits(shock_eq_map, "Matrix")) {
        stop("shock_eq_map should be of Matrix class")
    } else mod@shock_eq_map <- shock_eq_map

    if (!inherits(var_ceq_map, "Matrix")) {
        stop("var_ceq_map should be of Matrix class")
    } else mod@var_ceq_map <- var_ceq_map

    if (!inherits(cpar_eq_map, "Matrix")) {
        stop("cpar_eq_map should be of Matrix class")
    } else mod@cpar_eq_map <- cpar_eq_map

    if (!inherits(cpar_ceq_map, "Matrix")) {
        stop("cpar_ceq_map should be of Matrix class")
    } else mod@cpar_ceq_map <- cpar_ceq_map

    if (!inherits(fpar_eq_map, "Matrix")) {
        stop("fpar_eq_map should be of Matrix class")
    } else mod@fpar_eq_map <- fpar_eq_map

    if (!inherits(fpar_ceq_map, "Matrix")) {
        stop("fpar_ceq_map should be of Matrix class")
    } else mod@fpar_ceq_map <- fpar_ceq_map

    if (!is.function(ss_function)) {
        stop("ss_function should be a function")
    } else mod@ss_function <- ss_function

    if (!is.function(calibr_function)) {
        stop("calibr_function should be a function")
    } else mod@calibr_function <- calibr_function

    if (!is.function(ss_calibr_jac_function) & !is.null(ss_calibr_jac_function)) {
        stop("ss_calibr_jac_function should be a function or a NULL value")
    } else mod@ss_calibr_jac_function <- ss_calibr_jac_function

    if (!is.function(pert)) {
        stop("pert should be a function")
    } else mod@pert <- pert

    if (!is.list(ext)) {
        stop("ext should be a list")
    } else mod@ext <- ext
    
    if (length(parameters) == length(parameters_free)) {
        mod@is_calibrated <- FALSE
    }

    if (!length(shock_eq_map))
        mod@is_stochastic <- FALSE

    if (mod@is_stochastic) {
        mod@active_shocks <- rep(TRUE, length(mod@shocks))
        mod@shock_cov_mat <- diag(length(mod@shocks))
    }

    if (any((var_eq_map != 2) & (var_eq_map != 0))) {
        mod@is_dynamic <- TRUE
    }

    mod@variables_ss_val <- as.numeric(rep(NA, length(mod@variables)))
    mod@parameters_calibr_val <- as.numeric(rep(NA, length(mod@parameters_calibr)))
    mod@variables_init_val <- as.numeric(rep(NA, length(mod@variables)))
    mod@parameters_calibr_init_val <- as.numeric(rep(NA, length(mod@parameters_calibr)))
    mod@parameters_free_mod_flag <- rep(FALSE, length(mod@parameters_free))
    mod@map_free_into_params <- match(mod@parameters_free, mod@parameters)
    mod@map_calibr_into_params <- match(mod@parameters_calibr, mod@parameters)
    mod@parameters_val[mod@map_free_into_params] <- mod@parameters_free_val

    return (mod)
}

# ############################################################################
# The is.gecon_model function checks if given object
# is of class "gecon_model"
# ############################################################################
# Input
#   x - any R object
# Output
#   Logical value indicating if object is of class "gecon_model"
# ############################################################################
is.gecon_model <- function(x)
{
    if (is(x, "gecon_model")) return (TRUE)
    return (FALSE)
}



# ############################################################################
# Function list_eq returns equations with specified indices
# ############################################################################
# Input
#   model - object of gecon_model class
#   eq_idx -  indices of equations to be returned
# Output
#   character matrix - equations with given indices
# ############################################################################
list_eq <- function(model, eq_idx = NULL)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }

    if (is.null(eq_idx)) {
        eq_idx <- (1:length(model@equations))
    }

    if (any(eq_idx < 1 | eq_idx > length(model@equations))) {
        stop("invalid equation index")
    }

    eq <- matrix(model@equations[eq_idx], length(eq_idx), 1)
    rownames(eq) <- paste0("Eq. ", eq_idx, ": ")
    colnames(eq) <- ""

    return (eq)
}

# ############################################################################
# Function list_calibr_eq returns calibrating equations with specified indices
# ############################################################################
# Input
#   model - object of gecon_model class
#   eq_idx -  indices of equations to be returned
# Output
#   character matrix - calibrating equations with given indices
# ############################################################################
list_calibr_eq <- function(model, eq_idx = NULL)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }

    if (is.null(eq_idx)) {
        if (!length(model@calibr_equations))
            stop("the model has no calibrating equations")
        eq_idx <- (1:length(model@calibr_equations))
    }

    if (any(eq_idx < 1 | eq_idx > length(model@calibr_equations))) {
        stop("invalid calibrating equation index")
    }

    eq <- matrix(model@calibr_equations[eq_idx], length(eq_idx), 1)
    rownames(eq) <- paste0("Calibr. Eq. ", eq_idx, ": ")
    colnames(eq) <- ""

    return (eq)
}



# ############################################################################
# The get_index_sets function allows to retrieve a list with
# the index sets from an object of gecon_model class
# ############################################################################
# Input
#   model - object of the gecon_model class.
# Output
#   list with index sets specified for the model.
# ############################################################################
get_index_sets <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    if (!length(model@index_sets))
        stop("no index sets have been specified")

    return (model@index_sets)
}


# ############################################################################
# The re_solved function returns a logical value indicating
# if the perturbation has been solved.
# ############################################################################
# Input
#   model - object of the gecon_model class.
# Output
#   Logical value. If TRUE, the model has been solved.
# ############################################################################
re_solved <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    return (model@re_solved)
}


# ############################################################################
# The ss_solved function returns a logical value indicating
# if the steady state (equilibrium) for the model has been found.
# ############################################################################
# Input
#   model - object of the gecon_model class.
# Output
#   Logical value. If TRUE, the steady state (equilibrium)
#   for the model has been found.
# ############################################################################
ss_solved <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    return (model@ss_solved)
}


# ############################################################################
# The get_model_info function returns a character vector
# with information about the model.
# ############################################################################
# Input
#   model - object of the gecon_model class.
# Output
#   character vector of length 3, containing information
#   about the model: the input file name, the input file path,
#   and the date of creation.
# ############################################################################
get_model_info <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    return (model@model_info)
}




# ########################ja####################################################
# The show method controls how object of class gecon_model is printed
# ############################################################################
# Input
#   object - an object of class gecon_model
# Output
#   None
# ############################################################################
setMethod("show", signature(object = "gecon_model"),
function(object)
{
    if (!is.gecon_model(object)) {
        stop("argument should be of gecon_model class")
    }
    model_name <- object@model_info[1]
    no_var <- length(object@variables)
    no_shocks <- length(object@shocks)
    no_par <- length(object@parameters)
    no_par_mod <- sum(object@parameters_free_mod_flag)
    no_calibr_par <- length(object@parameters_calibr)
    no_state <- length(object@state_var_indices)

#     cat("\n")
    if (!object@is_dynamic) {
        cat("\'", model_name, "\' is a static model.\n\n", sep = "")
    } else if (object@is_stochastic) {
        cat("\'", model_name, "\' is a dynamic, stochastic model.\n\n", sep = "")
    } else {
        cat("\'", model_name, "\' is a dynamic, deterministic model.\n\n", sep = "")
    }

    cat("Generated ", object@model_info[3], " from \'", add_gcn_ext(object@model_info[2]),
        "\'.\n\n", sep = "")

    if (object@is_dynamic)
        cat("Steady-state values of variables have ")
    else
        cat("Equilibrium has ")
    if (object@ss_solved) {
        cat("been FOUND.\n")
    } else {
        cat("NOT been FOUND.\n")
    }

    if (length(object@calibr_equations)) {
        if (object@is_calibrated && object@ss_solved) {
            if (object@is_dynamic) {
                cat("Calibrating equations have been taken",
                    "into account when solving for the steady state.\n")
            } else {
                cat("Calibrating equations have been taken",
                    "into account when solving for the equilibrium.\n")
            }
        } else if (!object@is_calibrated && object@ss_solved) {
            if (object@is_dynamic) {
                cat("Calibrating equations have not been taken",
                    "into account when solving for the steady state.\n")
            } else {
                cat("Calibrating equations have not been taken",
                    "into account when solving for the equilibrium.\n")
            }
        }
    }

    if (object@re_solved) {
        cat("The (log-)linearised version of the model has been SOLVED.\n")
    } else if (object@is_dynamic) {
        cat("The model in its (log-)linearised form has NOT been SOLVED.\n")
    } else {
        cat("\n")
    }

    if (object@corr_computed) {
        cat("Correlations have been COMPUTED.\n\n")
    } else if(object@is_stochastic) {
            cat("Correlations have NOT been COMPUTED.\n\n")
    }

    cat("Number of variables:", no_var, "\n")
    if (no_par) {
        cat("Number of parameters:", no_par, "\n")
    }
    if (no_calibr_par) {
        cat("Number of calibrated parameters:", no_calibr_par, "\n")
    }
    if (object@is_stochastic) {
        cat("Number of stochastic shocks:", no_shocks, "\n")
    }
    if (object@is_dynamic) {
        cat("Number of state variables: ")
        if (object@re_solved) {
            cat(no_state, "\n")
        } else {
            cat("?\n")
        }
    }
#     cat("\n")
})


# ############################################################################
# The print method prints basic information and diagnostics of
# gecon_model objects
# ############################################################################
# Input
#   x - an object of class gecon_model
# Output
#   None
# ############################################################################
setMethod("print", signature(x = "gecon_model"),
function(x)
{
    if (!is.gecon_model(x)) {
        stop("argument should be of gecon_model class")
    }
    object <- x

    # model_name <- deparse(substitute(x))
    model_name <- object@model_info[1]
    no_var <- length(object@variables)
    no_shocks <- length(object@shocks)
    no_par <- length(object@parameters)
    no_par_mod <- sum(object@parameters_free_mod_flag)
    no_calibr_par <- length(object@parameters_calibr)
    no_state <- length(object@state_var_indices)
    variables <- object@variables
    calibr_parameters <- object@parameters_calibr
    shocks <- object@shocks
    parameters <- object@parameters
    parameters_modified <-
        object@parameters_free[object@parameters_free_mod_flag]
    state_variables <-
        object@variables[object@state_var_indices]

#     cat("\n")
    if (!object@is_dynamic) {
        cat("\'", model_name, "\' is a static model.\n\n", sep = "")
    } else if (object@is_stochastic) {
        cat("\'", model_name, "\' is a dynamic, stochastic model.\n\n", sep = "")
    } else {
        cat("\'", model_name, "\' is a dynamic, deterministic model.\n\n", sep = "")
    }

    cat("Generated ", object@model_info[3], " from \'", add_gcn_ext(object@model_info[2]),
        "\'.\n\n", sep = "")

    if (object@is_dynamic) {
        cat("Steady-state values of variables have ")
    } else {
        cat("Equilibrium has ")
    }
    if (object@ss_solved) {
        cat("been FOUND.\n")
    } else {
        cat("NOT been FOUND.\n")
    }

    if (length(object@calibr_equations)) {
        if (object@is_calibrated && object@ss_solved) {
            if (object@is_dynamic) {
                cat("Calibrating equations have been taken \n")
                cat("into account when solving for the steady state. \n")
            } else {
                cat("Calibrating equations have been taken \n")
                cat("into account when solving for the equilibrium. \n")
            }
        } else if (!object@is_calibrated && object@ss_solved) {
            if (object@is_dynamic) {
                cat("Calibrating equations have not been taken \n")
                cat("into account when solving for the steady state. \n")
            } else {
                cat("Calibrating equations have not been taken \n")
                cat("into account when solving for the equilibrium. \n")
            }
        }
    }

    if (object@re_solved) {
        cat("The (log-)linearised version of the model has been SOLVED.\n")
    } else if (object@is_dynamic) {
            cat("The model in its (log-)linearised form has NOT been SOLVED.\n")
    } else {
        cat("\n")
    }

    if (object@corr_computed) {
        if (object@ref_var_idx) {
            cat(paste0("Correlations have been COMPUTED (with ",
                       object@variables[object@ref_var_idx],
                       " as the reference variable).\n"))
        } else {
            cat("Correlations have been COMPUTED.\n")
        }
    } else {
        if (object@is_stochastic) {
            cat("Correlations have NOT been COMPUTED.\n")
        }
    }

    cat("\n")
    cat(paste0("Variables (", no_var, "):\n"))
    cat(list2str2(variables))
    cat("\n")

    if (no_par) {
        cat(paste0("Parameters (", no_par, "):\n"))
        cat(list2str2(parameters))
        cat("\n")
    }

    if(no_par_mod) {
        cat(paste0("Parameters whose values have been modified at R level (", no_par_mod, "):\n"))
        cat(list2str2(parameters_modified))
        cat("\n")
    }

    if (no_calibr_par) {
        cat(paste0("Calibrated parameters (", no_calibr_par, "):\n"))
        cat(list2str2(calibr_parameters))
        cat("\n")
    }

    if (object@re_solved) {
        if (no_state) {
            cat(paste0("State variables (", no_state, "):\n"))
            cat(list2str2(state_variables))
            cat("\n")
        } else {
            cat("The are NO STATE VARIABLES in the model.",
                "Solution matrices are empty.\n")
        }
    } else if (object@is_dynamic) {
            cat("In order to identify state variables you have to solve",
                "the model in its (log-)linearised form first.\n")
    }

    if (no_shocks) {
        cat(paste0("Stochastic shocks (", no_shocks, "):\n"))
        cat(list2str2(shocks))
        cat("\n")
    }

    cat("\n")
})


# ############################################################################
# The summary method prints and returns the summary of model solution
# ############################################################################
# Input
#   object - an object of class gecon_model
# Output
#   Summary of model solution
# ############################################################################
setMethod("summary", signature(object = "gecon_model"),
function (object)
{
    if (!is.gecon_model(object)) {
        stop("argument should be of gecon_model class");
    }
    if (!object@ss_solved) {
        cat("Model has NOT been SOLVED.\n")
        return (invisible())
    }
    if (object@is_dynamic) {
        cat("\nSteady state:\n")
    } else {
        cat("\nEquilibrium:\n")
    }
    cat_ss <- as.data.frame(object@variables_ss_val)
    colnames(cat_ss) <- c("")
    rownames(cat_ss) <- object@variables
    print(round(cat_ss, digits = 6))

    if (length(object@parameters_val)) {
        cat("\n----------------------------------------------------------\n")
        cat("\nParameter values:\n")
        cat_par <- as.data.frame(object@parameters_val)
        colnames(cat_par) <- c("")
        rownames(cat_par) <- object@parameters
        print(cat_par)
    }

    if (object@is_dynamic) {
        cat("\n----------------------------------------------------------\n")
    } else {
        return (invisible())
    }

    if (!object@re_solved) {
        cat("\nThe model in its (log-)linearised form has NOT been SOLVED.\n\n")
        if (length(object@solver_exit_info)) {
            cat("Solver exit infomation:\n", object@solver_exit_info, "\n", sep = "")
        }
        return (invisible())
    }

    if (!length(object@state_var_indices)) {
        cat("The are NO STATE VARIABLES in the model.",
            "Solution matrices are empty.\n")
        return (invisible())
    }

    cat("\nLinearisation:\n")
    cat("x_t = P x_{t-1} + Q epsilon_t\n")
    cat("y_t = R x_{t-1} + S epsilon_t\n")

    cat("\nP:\n")
    print(round(object@solution$P, digits = 6))

    if (object@is_stochastic) {
        cat("\nQ:\n")
        print(round(object@solution$Q, digits = 6))
    }

    cat("\nR:\n")
    print(round(object@solution$R, digits = 6))

    if (object@is_stochastic) {
        cat("\nS:\n")
        print(round(object@solution$S, digits = 6))
    }

    if (object@is_stochastic) {
        cat("\n----------------------------------------------------------\n")
    } else {
        return (invisible())
    }

    cat("\nShock covariance matrix:\n\n")
    rownames(object@shock_cov_mat) <- object@shocks
    colnames(object@shock_cov_mat) <- object@shocks
    print(object@shock_cov_mat)

    if (!object@corr_computed) {
        cat("\nCorrelations have NOT been COMPUTED.\n\n")
        return (invisible())
    }

    cat("\nBasic statistics:\n\n")
    loglin_indic <- rep("Y   ", length(object@loglin_var))
    loglin_indic[which(!object@loglin_var)] <- "N   "
    cat_mom <- round(cbind(object@variables_ss_val,  object@sdev, object@sdev ^ 2), digits = 4)
    cat_mom <- as.data.frame(cbind(cat_mom, loglin_indic))
    rownames(cat_mom) <- object@variables
    colnames(cat_mom) <- c("Steady-state value", "Std. dev.", " Variance", "Loglin")
    print(cat_mom)

    nzind <- which(object@sdev != 0)

    cat("\nCorrelation matrix:\n\n")
    cat_corrm <- object@corr_mat[nzind, nzind]
    cat_corrm <- as.data.frame(cat_corrm)
    N <- length(nzind)
    for (j in (1:N)) {
        cat_corrm[, j] <- as.character(round(cat_corrm[, j], digits = 3))
        if (j < N) {
            for (i in ((j + 1):N)) {
                cat_corrm[i, j] <- ""
            }
        }
    }
    print(cat_corrm)

    cat("\nAutocorrelations:\n\n")
    print(round(object@autocorr_mat[nzind, ], digits = 3))

    if (object@ref_var_idx) {
        cat(paste0("\nCross correlations with reference variable (",
                   object@variables[object@ref_var_idx], "):\n\n"))
        print(round(object@ref_var_corr_mat[nzind, ], digits = 3))
    }
    cat("\n")

    if (dim(object@var_dec)[2]) {
        cat("Variance decomposition:\n\n")
        cat_vd <- object@var_dec[nzind, ]
        dim(cat_vd) <- c(length(nzind), dim(object@var_dec)[2])
        colnames(cat_vd) <- colnames(object@var_dec)
        rownames(cat_vd) <- rownames(object@var_dec)[nzind]
        print(round(cat_vd, digits = 3))
        cat("\n")
    }
})
