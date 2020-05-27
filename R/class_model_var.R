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
# Model variables
# ############################################################################


# ############################################################################
# The get_var_names function returns names of variables
# stored in an object of gecon_model class
# ############################################################################
# Input
#   model - an object of the gecon_model class
# Output
#   returns a vector of variables' names
# ############################################################################
get_var_names <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    return (model@variables)
}



# ############################################################################
# The get_var_names_by_index function returns
# the names of variables indexed by the indices specified
# in the index_names argument
# ############################################################################
# Input
#   model - an object of the gecon_model class
#   index_names - a character vector of indices
# Output
#   A character vector with relevant variables' names
# ############################################################################
get_var_names_by_index <- function(model, index_names)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    if (!is.character(index_names))
        stop("index_names argument has to be of character type")

    if (!length(model@index_sets))
        stop("no index sets have been specified in this model")

    # Check if index sets are correctly specified
    ind_sets <- unlist(model@index_sets)

    not_indices <- which(!(index_names %in% ind_sets))

    if (length(not_indices))
        stop(paste0("the following are not valid indices: ",
                    list2str(not_indices, "\'")))

    # Finding relevant variables
    index_names_alternative <- paste("(", paste(index_names,  collapse = "|"), ")", 
                                     collapse = "", sep = "")
    assoc_names  <- grep(pattern = paste("__", index_names_alternative,
                                         "(__)*(?!(_?[a-zA-Z0-9]))", sep = ""),
                                x = model@variables,
                                perl = T)

    if (!length(assoc_names))
        stop("no variables with specified indices found")

    return (model@variables[assoc_names])
}





# ############################################################################
# The initval_var function sets the initial values of variables
# to values specified by the user.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   init_var -  a (named) list or vector of initial values of variables.
# Output
#   An updated object of gecon_model class.
# ############################################################################
initval_var <- function(model, init_var, warnings = TRUE)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    init_var_names <- names(init_var)
    if (is.null(init_var_names))
        stop("init_var argument has to be a list or a vector with named elements")

    is_proper_numeric <- tryCatch(expr = {init_var <- as.numeric(init_var)
                                            TRUE
                                       },
                               warning = function(w) FALSE,
                               error = function(w) FALSE)
                               
    if (!is_proper_numeric) {
        stop(paste("the following list elements cannot be coerced into initial variables' values:",
                    paste(init_var_names[union(which(as.logical(lapply(init_var, is.character))),
                          which(as.numeric(lapply(init_var, length)) > 1))], 
                          collapse = ", ")))
    }
    
    if (any(is.nan(as.numeric(init_var)) |
        is.infinite(as.numeric(init_var)))) {

        # storing names before retrieving values
        init_var <- as.numeric(init_var)

        # NaN and Inf indices
        nan_ind <- which(is.nan(init_var))
        inf_ind <- which(is.infinite(init_var))

        # names of misspecified parameters
        misspec_ind <- c(nan_ind, inf_ind)
        misspec_names <- init_var_names[misspec_ind]

        # labels for the print
        labs <- rep("Inf", length(misspec_ind))
        labs[which(is.nan(init_var[misspec_ind]))] <- "NaN"

        stop("initial values of the following variables have been specified incorrectly: ",
             paste(paste0("\"", misspec_names, "\"", paste0(" (=", labs, ")")),
                   collapse = ", "))
    }

    if (any(is.na(as.numeric(init_var)))) {
        na_ind <- which(is.na(as.numeric(init_var)))
        if (warnings) {
            warning("initial values of the following variables have not been modified ",
                    "(NA values have been supplied): ",
                    list2str(init_var_names[na_ind], "\""))
        }
        init_var <- init_var[-na_ind]
        if (!length(init_var)) {
            if (warnings) {
                warning("no initial values of model variables have been modified")
            }
            return (model)
        }
    }

    not_variables <- which(!(init_var_names %in% model@variables))
    if (length((not_variables)) != 0) {
        stop(paste0("the following are not model variables: ",
                   list2str(init_var_names[not_variables], "\"")))
    }

    if (model@ss_solved) {
        if (warnings) {
            warning("new initial values have been set, the existing solution has been erased")
        }
    }

    matches <- which(model@variables %in% init_var_names)
    ord <- match(model@variables[matches], init_var_names)
    model@variables_init_val[matches] <- as.numeric(init_var[ord])

    if (length(matches) != length(model@variables)) {
        if (warnings) {
            warning(paste0("initial values of the following variables have not been set: ",
                           list2str(model@variables[-matches], "\"")))
        }
    }

    # clearing slots
    model@is_calibrated <- TRUE
    if (length(model@parameters) == length(model@parameters_free)) {
        model@is_calibrated <- FALSE
    }
    model@loglin_var <- logical(length = 0)
    model@re_solved <- FALSE
    model@corr_computed <- FALSE
    model@residual_vector <- numeric(0)
    model@solver_status <- character(0)
    model@ss_solved <- FALSE
    model@eig_vals <- matrix(nrow = 0, ncol = 0)
    model@solution <- list(P = NULL,
                          Q = NULL,
                          R = NULL,
                          S = NULL)
    model@state_var_indices <- numeric(0)
    model@solver_exit_info <- character(0)
    model@solution_resid <- list(NULL)
    model@corr_mat <- matrix(nrow = 0, ncol = 0)
    model@autocorr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_corr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_idx <- 0L
    model@var_dec <- matrix(nrow = 0, ncol = 0)
    model@sdev <- matrix(nrow = 0, ncol = 0)
    return (model)
}


# ############################################################################
# The get_init_val_var function prints and returns initial values
# of model variables
# ############################################################################
# Input
#   model - an object of gecon_model class.
#   variables - names or indices of the variables, whose initial values 
#               are to be returned, default option is a vector of all variables.
#   to_tex - logical, if TRUE results are written to .tex file
#   silent - logical, if TRUE, console output is suppressed (FALSE by default).
# Output
#   returns a vector of initial values.
# ############################################################################
get_init_val_var <- function(model, variables = NULL,
                             to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")
    
    if (is.null(variables)) {
        var_ind <- 1:length(model@variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
    }

    init_vector <- model@variables_init_val[var_ind]
    init_names <- model@variables[var_ind]
    
    def_ind <- which(init_vector %in% 0.9)
    if (!length(def_ind) == 0) {
        warning(paste0("In case of the following variables the default value of 0.9 ",
                       "is used as their initial value: ", 
                       paste0(init_names[def_ind], collapse = ", "), "."))
    }
    
    if (to_tex || !silent) {
        init_val <- matrix(round(init_vector, digits = 4), length(var_ind), 1)
        colnames(init_val) <- "Initial value"
    }

    if (to_tex) {
        texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
        tex <- "\\section{Initial values}\n\n"
        rownames(init_val) <- paste0("$", model@variables_tex[var_ind], "$")
        tex <- paste0(tex, df2textable(init_val))
    }
    if (!silent) {
        cat("Initial values:\n\n")
        rownames(init_val) <- init_names
        print(round(init_val, digits = 4))
    }
    if (to_tex) {
        write_tex(texf, tex, !silent)
    }

    names(init_vector) <- init_names
    return (invisible(init_vector))
}


# ############################################################################
# The get_ss_values function prints and returns steady-state
# (equilibrium) values of model variables
# ############################################################################
# Input
#   model - an object of gecon_model class.
#   variables - names or indices of the variables, whose steady state
#               (or equilibrium) values are to be returned,
#               default option is a vector of all variables.
#   to_tex - logical, if TRUE results are written to .tex file
#   silent - logical, if TRUE, console output is suppressed (FALSE by default).
# Output
#   returns a vector of steady-state (equilibrium) values.
# ############################################################################
get_ss_values <- function(model, variables = NULL,
                          to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")
    if (!model@ss_solved) {
        if (all(is.na(model@variables_ss_val))) {
            if (model@is_dynamic)
                stop(paste0("Compute steady state first using \'steady_state\' function. ",
                            "The solver has not been started or has not returned any values."))
            else
                stop(paste0("compute equilibrium first using \'steady_state\' function. ",
                            "The solver has not been started or has not returned any values."))
        } else {
            warning(paste0("The steady state solver has not converged at given tolerance level. ",
                           "The solver status was:\n", model@solver_status, "\n",
                           "Presented results are not the solution to the steady-state/equlibrium system ",    
                           "at given tolerance level. They come from the last solver iteration."))    
        }
    }

    if (is.null(variables)) {
        var_ind <- 1:length(model@variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
    }

    ss_vector <- model@variables_ss_val[var_ind]
    ss_names <- model@variables[var_ind]

    if (to_tex || !silent) {
        ss_val <- matrix(round(ss_vector, digits = 4), length(var_ind), 1)
        if (model@is_dynamic) {
            colnames(ss_val) <- "Steady-state value"
        } else {
            colnames(ss_val) <- "Equilibrium value"
        }
    }

    if (to_tex) {
        texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
        if (model@is_dynamic) {
            tex <- "\\section{Steady-state values}\n\n"
        } else {
            tex <- "\\section{Equilibrium values}\n\n"
        }
        rownames(ss_val) <- paste0("$", model@variables_tex[var_ind], "$")
        tex <- paste0(tex, df2textable(ss_val))
    }
    if (!silent) {
        if (model@is_dynamic) {
            cat("Steady-state values:\n\n")
        } else {
            cat("Equilibrium values:\n\n")
        }
        rownames(ss_val) <- ss_names
        print(round(ss_val, digits = 4))
    }
    if (to_tex) {
        write_tex(texf, tex, !silent)
    }

    names(ss_vector) <- ss_names
    return (invisible(ss_vector))
}




# ############################################################################
# var_info function returns and prints information about model variables
# ############################################################################
# Input
#   model - an object of gecon_model class.
#   variables - the names or indices of the variables of interest.
#   all - logical value. If TRUE, information about all model variables 
#         is generated (default is FALSE).                                        
# Output
#   An object of gecon_var_info class.
# ############################################################################
var_info <- function(model, variables = NULL, all = FALSE)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    if (is.null(variables) && !all)
        stop("no variables have been specified")

    if (all) {
        if (!is.null(variables)) {
            warning("ignoring variables argument when all option is set to TRUE")
        }
        variables <- model@variables
        var_ind <- 1:length(variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
        variables <- model@variables[var_ind]
    }

    # Creating map of parameters
    if (length(model@var_ceq_map)) {
        var_ceq_map <- model@var_ceq_map
        var_ceq_map[which(var_ceq_map == 1)] <- 8
    } else var_ceq_map <-
        as(model@var_ceq_map, Class = "dgCMatrix")

    var_eq_map <- rbind(model@var_eq_map,
                         var_ceq_map)

    end_eq_block <- nrow(var_eq_map) - nrow(model@cpar_ceq_map)

    var_eq_map <- var_eq_map[, var_ind, drop = FALSE]
    incid_eq <- sort(unique(var_eq_map@i)) + 1
    var_eq_map <- var_eq_map[incid_eq, , drop = FALSE]

    colnames(var_eq_map) <- model@variables[var_ind]

    if (any(incid_eq > end_eq_block)) {
        calibr_eq_numb <- which(incid_eq > end_eq_block)
        incid_eq[calibr_eq_numb]  <- incid_eq[calibr_eq_numb] - end_eq_block
        rnames <- vector(length = length(incid_eq))
        rnames[-c(calibr_eq_numb)] <-
            paste("Equation    ", incid_eq[-c(calibr_eq_numb)])
        rnames[calibr_eq_numb] <-
            paste("Calibr. Eq. ", incid_eq[calibr_eq_numb])
        rownames(var_eq_map) <- rnames
    } else {
        rownames(var_eq_map) <- paste("Equation    ", incid_eq)
    }
    
    #Current initial values
    initval <- matrix(model@variables_init_val[var_ind], length(var_ind), 1)
    rownames(initval) <- model@variables[var_ind]
    colnames(initval) <- "Initial values"

    # SS information
    if (model@ss_solved) {
        # Steady state values
        ss_val <- matrix(model@variables_ss_val[var_ind], length(var_ind), 1)
        rownames(ss_val) <- model@variables[var_ind]
        if (model@is_dynamic) {
            colnames(ss_val) <- "Steady state"
        } else {
            colnames(ss_val) <- "Equilibrium"
        }
    } else if (!all(is.na(model@variables_ss_val))) {
        ss_val <- matrix(model@variables_ss_val[var_ind], length(var_ind), 1)
        rownames(ss_val) <- model@variables[var_ind]
        colnames(ss_val) <- "Last solver iteration"
    } else {
        # Initial values
        cat("The steady-state solver has not been started or has not returned any values.\n")
        ss_val <- initval
    }

    # Block obtained when perturbation is solved
    if (model@re_solved) {
        # Lag info
        state <- logical(length = length(var_ind))
        state[which(var_ind %in% model@state_var_indices)] <- TRUE

        # Get impact of state and exog variables
        state_var_impact <- matrix(0, length(var_ind),
                        length(model@state_var_indices))
        shock_impact <- matrix(0, length(var_ind), length(model@shocks))


        # state variables
        svi <- model@state_var_indices
        if (any(var_ind %in% svi)) {
            ind <- which(model@variables[model@state_var_indices] %in%
                    model@variables[svi[which(svi %in% var_ind)]])
            state_var_impact[var_ind %in% svi, ] <- model@solution$P[ind, ]
            shock_impact[var_ind %in% svi, ] <- model@solution$Q[ind, ]
        }

        non_svi <- c(1:length(model@variables))
        non_svi <- non_svi[-model@state_var_indices]

        # non state variables
        no_s_ind <- which(model@variables[non_svi] %in%
                            model@variables[non_svi[which(non_svi %in% var_ind)]])
        state_var_impact[var_ind %in% non_svi, ] <- model@solution$R[no_s_ind, ]
        shock_impact[var_ind %in% non_svi, ] <- model@solution$S[no_s_ind, ]

        # adding names
        rownames(state_var_impact) <- model@variables[var_ind]
        rownames(shock_impact) <- model@variables[var_ind]
        colnames(state_var_impact) <- paste(model@variables[model@state_var_indices], "[-1]", sep = "")
        colnames(shock_impact) <- model@shocks
        loglin_indic <- model@loglin_var[var_ind]
    } else {
        loglin_indic <- logical(length = 0)
        state <- vector(length = 0)
        state_var_impact <- matrix(nrow = 0, ncol = 0)
        shock_impact <- matrix(nrow = 0, ncol = 0)
    }

    # Block obtained when correlations are computed
    if (model@corr_computed) {
        # Moments
        steady_state_val <- model@variables_ss_val[var_ind]
        std_dev_val <- model@sdev[var_ind]

        # Correlations
        nzind <- which(model@sdev != 0)
        nzinds <- intersect(var_ind, nzind)
        cr <- matrix(model@corr_mat[nzinds, nzind], length(nzinds), length(nzind))
        rownames(cr) <- model@variables[nzinds]
        colnames(cr) <- model@variables[nzind]
    } else {
        std_dev_val <- numeric(length = 0)
        cr <- matrix(nrow = 0, ncol = 0)
    }

    new_var_info <- gecon_var_info(r_object_name = deparse(substitute(model)),
                                   variables = model@variables[var_ind],
                                   is_stochastic = model@is_stochastic,
                                   is_dynamic = model@is_dynamic,
                                   ss_solved = model@ss_solved,
                                   re_solved = model@re_solved,
                                   corr_computed = model@corr_computed,
                                   ss_val = ss_val,
                                   init_val = initval,
                                   state = state,
                                   state_var_impact = state_var_impact,
                                   shock_impact = shock_impact,
                                   std_dev_val = std_dev_val,
                                   loglin_flag = loglin_indic,
                                   cr = cr,
                                   incid_mat = var_eq_map)
    return (new_var_info)
}
