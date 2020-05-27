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
# Model parameters
# ############################################################################


# ############################################################################
# The get_par_names function returns the names of
# parameters stored in an object of gecon_model class.
# ############################################################################
# Input
#   model - an object of gecon_model class
#   free_par - logical, controls if free parameters should be
#              added to the vector of parameter names.
#   calibr_par - logical, controls if calibrated parameters should be
#                added to the vector of parameter names.
# Output
#   returns character vector of parameter names.
# ############################################################################
get_par_names <- function(model, free_par = TRUE, calibr_par = TRUE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    if (!free_par & !calibr_par)
        stop("at least one of the free_par and calibr_par options has to be set to TRUE.")

    ret_list <- character(0)

    if (calibr_par)
        ret_list <- model@parameters_calibr

    if (free_par)
        ret_list <- c(ret_list, model@parameters_free)

    return (ret_list)
}

# ############################################################################
# The get_par_names_by_index function returns
# the names of parameters indexed by the indices specified
# in the index_names argument
# ############################################################################
# Input
#   model - an object of the gecon_model class
#   index_names - a character vector of indices
# Output
#   A character vector with relevant parameters' names
# ############################################################################
get_par_names_by_index <- function(model, index_names)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    if (!is.character(index_names))
        stop("expected index_names argument should be a character vector")

    if (!length(model@index_sets))
        stop("no index sets have been specified in this model")

    # Check if index sets are correctly specified
    ind_sets <- unlist(model@index_sets)

    not_indices <- which(!(index_names %in% ind_sets))

    if (length(not_indices))
        stop(paste0("the following are not valid indices: ",
                    list2str(not_indices, "\'")))

    # Finding relevant parameters
    index_names_alternative <- paste("(", paste(index_names,  collapse = "|"), ")",
                                     collapse = "", sep = "")
    assoc_names  <- grep(pattern = paste("__", index_names_alternative,
                                         "(__)*(?!(_?[a-zA-Z0-9]))", sep = ""),
                                x = model@parameters,
                                perl = T)

    if (!length(assoc_names))
        stop("no parameters with specified indices found")

    return (model@parameters[assoc_names])
}



# ############################################################################
# The set_free_par function assigns parameters' values to gecon_model object
# ############################################################################
# Input
#   model - an object of the gecon_model class, to which
#           we want to assign parameters
#   free_par - a (named) list or vector of parameters
#   reset - logical, allows to reset the values parameters to values
#           specified in .gcn file (if set to TRUE)
#   warnings - logical, if TRUE, a warning is displayed whenever default
#            (specified in .gcn file) parameters' values are overwritten
# Output
#   A gecon_model object with added/updated free parameters' values
# ############################################################################
set_free_par <- function(model, free_par = NULL, reset = FALSE, warnings = TRUE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    free_par_names <- names(free_par)
    if (is.null(free_par_names) & (reset == FALSE))
        stop("free_par argument has to be a list or a vector with named elements")

    is_proper_numeric <- tryCatch(expr = {free_par <- as.numeric(free_par)
                                                TRUE
                                         },
                                   warning = function(w) FALSE,
                                   error = function(w) FALSE)

    if (!is_proper_numeric) {
        stop(paste("the following elements of the list cannot be coerced into parameters' values:",
                    paste(free_par_names[union(which(as.logical(lapply(free_par, is.character))),
                                         which(as.numeric(lapply(free_par, length)) > 1))], collapse = ", ")))
    }

    if (any(is.nan(as.numeric(free_par)) |
            is.infinite(as.numeric(free_par)))) {
        # NaN and Inf indices
        nan_ind <- which(is.nan(free_par))
        inf_ind <- which(is.infinite(free_par))

        # names of misspecified parameters
        misspec_ind <- c(nan_ind, inf_ind)
        misspec_names <- free_par_names[misspec_ind]

        # labels for the print
        labs <- rep("Inf", length(misspec_ind))
        labs[which(is.nan(free_par[misspec_ind]))] <- "NaN"

        stop("values of the following free parameters have not been correctly specified: ",
             paste(paste0("\"", misspec_names, "\"", paste0(" (=", labs, ")")),
                   collapse = ", "))
    }

    if (!is.list(free_par) && !(is.numeric(free_par)) &&
        !(is.null(free_par) && reset))  {
        stop("parameters' values have to be passed as a numeric vector or a hash list")
    }

    # resetting to the default values
    if (reset) {
        if (all(is.na(model@parameters_free_init_val)))
            stop("no free parameters' values have been specified in the .gcn file")
        matches <- 1:length(model@parameters_free_val)
        free_par <- model@parameters_free_init_val
        mod_flag <- FALSE
    } else {
        # validating inputs
        not_free_parameters <-
            which(free_par_names %in% model@parameters_calibr)

        if (length(not_free_parameters)) {
            stop(paste0("the following are calibrated parameters ",
                        "(their values can be specified through a call to 'initval_calibr_par'): ",
                        list2str(free_par_names[not_free_parameters], "\"")))
        }

        not_parameters <-
            which(!(free_par_names %in% model@parameters_free))

        if (length(not_parameters)) {
            stop(paste0("the following are not model parameters: ",
                        list2str(free_par_names[not_parameters], "\"")))
        }

        # initializing
        matches <- which(model@parameters_free %in% free_par_names)
        ord <- match(model@parameters_free[matches], free_par_names)
        if (warnings) {
            ow <- which(matches %in%
                            which(!is.na(model@parameters_free_val)))
            ow <- matches[ow]
            if (length(ow)) {
                warning(paste0("the values of the following parameters have been overwritten: ",
                                list2str(model@parameters_free[ow], "\"")))
            }
        }
        free_par <- as.numeric(free_par[ord])
        mod_flag <- TRUE
    }

    # vector of parameters
    model@parameters_free_val[matches] <- free_par
    model@parameters_free_mod_flag[matches] <- mod_flag
    model@parameters_val[model@map_free_into_params[matches]] <- free_par

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

    ss_indic <- FALSE
    if (model@ss_solved) {
        if (length(model@parameters_calibr_val)) {
            ss_indic_res <- abs(model@ss_function(model@variables_ss_val,
                                                  model@parameters_calibr_val,
                                                  model@parameters_free_val))
            ss_indic <- (max(ss_indic_res[which(!is.na(ss_indic_res))]) > 1e-6) |
                         any(is.na(ss_indic_res))
        } else {
            ss_indic_res <- abs(model@ss_function(model@variables_ss_val,
                                                  numeric(0),
                                                  model@parameters_free_val))
            ss_indic <- (max(ss_indic_res[which(!is.na(ss_indic_res))]) > 1e-6) |
                         any(is.na(ss_indic_res))
        }
    }

    if (ss_indic)
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
# The initval_calibr_par function assigns initial values of
# calibrated parameters to gecon_model object.
# Calibrated parameters may be treated as free parameters
# when the calibration option is set to FALSE while solving for the steady state
# or equilibrium.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   calibr_par - a (named) list or vector of parameters.
#   warnings - logical, should warnings be displayed?
# Output
#   A gecon_model object with calibrated parameters set.
# ############################################################################
initval_calibr_par <- function(model, calibr_par, warnings = TRUE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    calibr_par_names <- names(calibr_par)
    if (is.null(calibr_par_names))
        stop("calibr_par argument has to be a list or a vector with named elements")

    is_proper_numeric <- tryCatch(expr = {distr_par <- as.numeric(calibr_par)
                                            TRUE
                                         },
                               warning = function(w) FALSE,
                               error = function(w) FALSE)
    if (!is_proper_numeric) {
        stop(paste("the following elements of the list cannot be coerced into parameters' values:",
                    paste(calibr_par_names[union(which(as.logical(lapply(calibr_par, is.character))),
                           which(as.numeric(lapply(calibr_par, length)) > 1))], 
                           collapse = ", ")))
    }


    if (any(is.nan(as.numeric(calibr_par)) |
        is.infinite(as.numeric(calibr_par)))) {

        # storing names before retrieving values
        calibr_par <- as.numeric(calibr_par)

        # NaN and Inf indices
        nan_ind <- which(is.nan(calibr_par))
        inf_ind <- which(is.infinite(calibr_par))

        # names of misspecified parameters
        misspec_ind <- c(nan_ind, inf_ind)
        misspec_names <- calibr_par_names[misspec_ind]

        # labels for the print
        labs <- rep("Inf", length(misspec_ind))
        labs[which(is.nan(calibr_par[misspec_ind]))] <- "NaN"

        stop("initial values of the following calibrated parameters have not been correctly specified: ",
             paste(paste0("\"", misspec_names, "\"", paste0(" (=", labs, ")")),
                   collapse = ", "))
    }

    if (any(is.na(as.numeric(calibr_par)))) {
        na_ind <- which(is.na(as.numeric(calibr_par)))
        if (warnings) {
            warning(paste0("initial values of the following calibrated parameters have not been modified ",
                           "(NA values have been supplied): ",
                           list2str(calibr_par_names[na_ind], "\"")))
        }
        calibr_par <- calibr_par[-na_ind]
        if (!length(calibr_par)) {
            if (warnings) {
                warning("none of the initial values of calibrated parameters has been modified")
            }
            return (model)
        }
    }

    # validating inputs
    not_free_parameters <-
        which(calibr_par_names %in% model@parameters_free)
    if (length((not_free_parameters)) != 0) {
        stop(paste0("the following are free parameters ",
                    "(their values can be specified through a call to 'set_free_par'): ",
                    list2str(calibr_par_names[not_free_parameters], "\"")))
    }

    # validating inputs
    not_parameters <- which(!calibr_par_names %in% model@parameters_calibr)
    if (length((not_parameters)) != 0) {
        stop(paste0("the following are not model parameters: ",
                    list2str(calibr_par_names[not_parameters], "\"")))
    }

    matches <- which(model@parameters_calibr %in% calibr_par_names)
    ord <- match(model@parameters_calibr[matches], calibr_par_names)
    model@parameters_calibr_init_val[matches] <- as.numeric(calibr_par[ord])
    model@parameters_val[model@map_calibr_into_params[matches]] <-
        as.numeric(calibr_par[ord])

    if (length(matches) != length(model@parameters_calibr)) {
        if (warnings) {
            warning(paste0("initial values of the following calibrated parameters have not been set: ",
                           list2str(model@parameters_calibr[-matches], "\"")))
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
# The get_init_calibr_par function retrieves calibrated parameters' initial values
# ############################################################################
# Input
#   model - an object of the gecon_model class
#    calibr_par - the names or indices of calibrated parameters whose initial  
#                values one wants to check, default option is to print all values
#   to_tex - logical, if TRUE parameter values are written to .tex file
#   silent - logical, if TRUE, console output is suppressed (FALSE by default).
# Output
#   returns vector of calibrated parameters' initial values
# ############################################################################
get_init_calibr_par <- function(model, calibr_par = NULL,
                                to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    if (is.null(calibr_par)) {
        par_ind <- which(model@parameters %in% model@parameters_calibr)
        calibr_par_names <- model@parameters_calibr
    } else {
        par_ind <- list2ind(calibr_par, model@parameters, "parameter")
        calibr_par_names <- model@parameters[par_ind]
    }

    # validating inputs
    not_free_parameters <- which(calibr_par_names %in% model@parameters_free)
    if (length((not_free_parameters)) != 0) {
        stop(paste0("the following are free parameters ",
                    "(their values can be checked through a call to 'get_par_values'): ",
                    list2str(calibr_par_names[not_free_parameters], "\"")))
    }

    matches <- which(model@parameters_calibr %in% calibr_par_names)
    par_vector <- model@parameters_calibr_init_val[matches]
    
    def_ind <- which(par_vector %in% 0.5)
    if (!length(def_ind) == 0) {
        warning(paste0("In case of the following parameters the default value of 0.5 ",
                       "is used as their initial value: ", 
                       paste0(calibr_par_names[def_ind], collapse = ", "), "."))
    }

    if (to_tex || !silent) {
        par_val <- matrix(round(par_vector, digits = 4), length(matches), 1)
        colnames(par_val) <- "Initial value"
    }

    if (to_tex) {
        texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
        tex <- "\\section{Calibrated parameters}\n\n"
        rownames(par_val) <- paste0("$", model@parameters_tex[par_ind], "$")
        tex <- paste0(tex, df2textable(par_val))
    }
    if (!silent) {
        rownames(par_val) <- calibr_par_names #parameters
        cat("Initial values of calibrated parameters:\n\n")
        print(round(par_val, digits = 3))
    }
    if (to_tex) {
        write_tex(texf, tex, !silent)
    }

    names(par_vector) <- calibr_par_names #parameters
    return (invisible(par_vector))
}

# ############################################################################
# The get_par_values function retrieves parameters' values
# ############################################################################
# Input
#   model - an object of the gecon_model class
#   parameters - the names (or indices) of parameters whose values one
#            wants to check, default option is to print all values
#   to_tex - logical, if TRUE parameter values are written to .tex file
#   silent - logical, if TRUE, console output is suppressed (FALSE by default).
# Output
#   returns vector of parameters' values
# ############################################################################
get_par_values <- function(model, parameters = NULL,
                           to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    if (is.null(parameters)) {
        par_ind <- 1:length(model@parameters)
    } else {
        par_ind <- list2ind(parameters, model@parameters, "parameter")
    }

    par_vector <- model@parameters_val[par_ind]
    parameters <- model@parameters[par_ind]
    
    calibr_ind <- which(parameters %in% model@parameters_calibr)
    
    if(length(calibr_ind) != 0) {
        if (!model@ss_solved && !all(is.na(model@parameters_calibr_val))) {
            warning(paste0("Presented results for calibrated parameters: ",
                           paste0(parameters[calibr_ind], collapse = ", "),
                           " are not the solution to the steady-state/equlibrium system at given tolerance level. ", 
                           "They come from the last solver iteration. ", 
                           "The steady state solver has not converged at given tolerance level. ",
                           "The solver status was:\n", model@solver_status, "\n"))    
        }
    }    

    if (to_tex || !silent) {
        par_val <- matrix(round(par_vector, digits = 4), length(par_ind), 1)
        colnames(par_val) <- "Value"
    }

    if (to_tex) {
        texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
        tex <- "\\section{Model parameters}\n\n"
        rownames(par_val) <- paste0("$", model@parameters_tex[par_ind], "$")
        tex <- paste0(tex, df2textable(par_val))
    }
    if (!silent) {
        rownames(par_val) <- parameters
        cat("Model parameters:\n\n")
        print(round(par_val, digits = 3))
    }
    if (to_tex) {
        write_tex(texf, tex, !silent)
    }

    names(par_vector) <- parameters
    return (invisible(par_vector))
}


# ############################################################################
# par_info function prints information about model parameters
# ############################################################################
# Input
#   model - an object of the gecon_model class
#   parameters - the names or indices of the parameters of interest.
#   all - logical value. If TRUE, information about all model parameters is generated
#         (default is FALSE).
# Output
#   An object of gecon_par_info class.
# ############################################################################
par_info <- function(model, parameters = NULL, all = FALSE)
{
    if (!is.gecon_model(model))
        stop("expected model argument should be of gecon_model class")

    if (is.null(parameters) && !all)
        stop("no parameters have been specified")

    if (all) {
        if (!is.null(parameters)) {
            warning("ignoring parameters argument when 'all' option is set to TRUE")
        }
        parameters <- model@parameters
        par_ind <- 1:length(parameters)
    } else {
        par_ind <- list2ind(parameters, model@parameters, "parameter")
        parameters <- model@parameters[par_ind]
    }

    calibr_par_ind <- match(parameters, model@parameters_calibr)
    calibr_par_ind <- calibr_par_ind[which(!is.na(calibr_par_ind))]

    free_par_ind <- match(parameters, model@parameters_free)
    free_par_ind <- free_par_ind[which(!is.na(free_par_ind))]

    # Creating map of parameters
    cpar_eq_map <- rbind(model@cpar_eq_map,
                         model@cpar_ceq_map)

    fpar_eq_map <- rbind(model@fpar_eq_map,
                         model@fpar_ceq_map)

    par_eq_map <- sparseMatrix(i = NULL, j = NULL,
                               dims = c(nrow(cpar_eq_map),
                                        length(par_ind)))

    par_eq_map <- as(object = par_eq_map, Class = "dgCMatrix")
    par_eq_map[, which(parameters %in% model@parameters_calibr)] <-
                cpar_eq_map[, calibr_par_ind]
    par_eq_map[, which(parameters %in% model@parameters_free)] <-
                fpar_eq_map[, free_par_ind]
    if (is.null(ncol(par_eq_map)))
        par_eq_map <- as(Matrix(par_eq_map), "dgCMatrix")
    incid_eq <- sort(unique(par_eq_map@i)) + 1
    par_eq_map <- par_eq_map[incid_eq, , drop = FALSE]

    colnames(par_eq_map) <- model@parameters[par_ind]
    end_eq_block <- nrow(model@var_eq_map)

    if (any(incid_eq > end_eq_block)) {
        calibr_eq_numb <- which(incid_eq > end_eq_block)
        incid_eq[calibr_eq_numb]  <- incid_eq[calibr_eq_numb] - end_eq_block
        rnames <- vector(length = length(incid_eq))
        rnames[-c(calibr_eq_numb)] <-
            paste("Equation    ", incid_eq[-c(calibr_eq_numb)])
        rnames[calibr_eq_numb] <-
            paste("Calibr. Eq. ", incid_eq[calibr_eq_numb])
        rownames(par_eq_map) <- rnames
    } else {
        rownames(par_eq_map) <- paste("Equation    ", incid_eq)
    }

    #Finding initial values specified in gcn file
    initial <- rep(NA, length(parameters))
    initial[which(parameters %in% model@parameters_free)] <-
                model@parameters_free_init_val[free_par_ind]

    calibr_flag <- logical(length(par_ind))
    calibr_flag[which(parameters %in% model@parameters_calibr)] <- TRUE

    
    calibr_initial <- rep(NA, length(parameters))
    calibr_initial[which(parameters %in% model@parameters_calibr)] <-
                    model@parameters_calibr_init_val[calibr_par_ind]                                            
    info <- gecon_par_info(r_object_name = deparse(substitute(model)),
                           parameters = model@parameters[par_ind],
                           gcn_values = initial,
                           current_values = model@parameters_val[par_ind],
                           calibr_flag = calibr_flag,
                           calibr_init_val = calibr_initial,
                           incid_mat = par_eq_map)

    return (info)
}
