# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Class for storing equations, variables, solution and statistics
# ###################################################################

# ###################################################################
# Class definition
# ###################################################################
setClass(
    Class = "gecon_model",
    representation = representation(
        # Info about the model
        model_info = "character",

        # Mapping of variables
        parameters = "character",
        parameters_free = "character",
        map_free_into_params = 'numeric',
        parameters_free_mod_flag = "numeric",
        parameters_calibr = "character",
        parameters_calibr_mod_flag = "numeric",
        map_calibr_into_params = 'numeric',
        variables = "character",
        shocks = "character",
        equations = "character",
        calibr_equations = "character",
        var_eq_map = "Matrix",
        shock_eq_map = "Matrix",
        var_ceq_map = "Matrix",
        cpar_eq_map = "Matrix",
        cpar_ceq_map = "Matrix",
        fpar_eq_map = "Matrix",
        fpar_ceq_map = "Matrix",
        
        # Controls
        init_calib_pars_supplied = "logical",
        init_vals_supplied = "logical",
        loglin = 'logical',
        loglin_var = 'logical',
        re_solved = 'logical',
        corr_computed = 'logical',
        is_stochastic = 'logical',
        is_dynamic = 'logical',
        is_calibrated = 'logical',

        # Steady-state parameters
        parameters_free_init_val = "vector",
        parameters_free_val = "vector",
        parameters_calibr_val = "numeric",
        params = "numeric",
        steady = "numeric",

        # Equilibrium relationships
        ss_function = "function",
        ss_function_jac = 'function',
        ss_calibr_function_jac = 'function',
        calibr_function = 'function',

        # Steady-state solution
        init_residual_vector = "numeric",
        residual_vector = "numeric",
        solver_status = "character",
        ss_solved = "logical",

        # Necessary to solve RE
        pert = 'function',
        eig_vals = 'matrix',

        # Solution
        solution = "list",
        state_var_indices = "numeric",
        solver_exit_info = "character",
        solution_resid = "list",

        # Stochastic simulation
        active_shocks = "logical",
        cov_mat = "matrix",
        shock_mat_flag = "logical",
        corr_mat = "matrix",
        autocorr_mat = "matrix",
        corr_variable_mat = "matrix",
        var_position = "numeric",
        var_dec = "matrix",
        sdev = "matrix",
        
        # Index sets
        index_sets = "list"
    ) ,
    prototype = prototype(
        parameters = character(0),
        variables = character(0),
        shocks = character(0),
        init_vals_supplied = FALSE,
        init_calib_pars_supplied = FALSE,
        var_position = numeric(0),
        solution_resid = list(),
        ss_solved = FALSE,
        init_residual_vector = numeric(0),
        residual_vector = numeric(0),
        pert = function(x) NULL,
        corr_computed = FALSE,
        loglin = TRUE,
        is_calibrated = TRUE,
        re_solved = FALSE,
        is_dynamic = FALSE,
        is_stochastic = TRUE,
        shock_mat_flag = FALSE,
        solution = list(P = NULL,
                        Q = NULL,
                        R = NULL,
                        S = NULL)
    )
)

# ###################################################################
# The gecon_model function is a constructor of the class gecon_model
# ###################################################################
# Input
#   model_info - [character vector, length = 3] information about
#                   model: input file name, file path and the date of creation
#   index_sets - [list] each of the list components corresponds to one
#                 set specified in gecon model class. Each components stores
#                 all the elements in the set.
#   variables - [character vector] of all variable names
#   shocks - [character vector]  of all shock names
#   parameters - [character vector] of all parameter names
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
#   ss_calibr_function_jac - [function] Jacobian of system
#                             of functions defining steady state
#                             a dynamic model or equilibrium
#                             in a static model
#                             and calibration equations
#   pert - [function] the function returning perturbation matrices.
#
# Output
#   An object of class "gecon_model"
# ###################################################################
gecon_model <- function(model_info,
                        index_sets,
                        variables,
                        shocks,
                        parameters,
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
                        ss_calibr_function_jac,
                        pert)
{
    mod <- new('gecon_model')

    if (!is.character(model_info)) {
        stop('model_info should be of character type')
    } else mod@model_info = model_info

    if (!is.list(index_sets)) {
        stop('index_sets must be a list of sets and their elements')
    } else mod@index_sets = index_sets
    
    if (!is.character(variables)) {
        stop('variables should be of character type')
    } else mod@variables = variables

    if (!is.character(shocks)) {
        stop('shocks should be of character type')
    } else mod@shocks = shocks

    if (!is.character(parameters)) {
        stop('parameters should be of character type')
    } else mod@parameters = parameters

    if (!is.character(parameters_free)) {
        stop('parameters_free should be of character type')
    } else mod@parameters_free = parameters_free

    mod@parameters_calibr <-
        parameters[-which(parameters %in% parameters_free)]

    if (!is.logical(parameters_free_val) & !is.numeric(parameters_free_val)) {
        stop('parameters_free_val should be a numeric or logical vector')
    } else {
        mod@parameters_free_val = as.numeric(parameters_free_val)
        mod@parameters_free_init_val = as.numeric(parameters_free_val)
    }

    if (!is.character(equations)) {
        stop('equations should be of character type')
    } else mod@equations = equations

    if (!is.character(calibr_equations)) {
        stop('equations should be of character type')
    } else mod@calibr_equations = calibr_equations

    if (!inherits(var_eq_map, 'Matrix')) {
        stop('var_eq_map should be of Matrix type')
    } else mod@var_eq_map = var_eq_map

    if (!inherits(shock_eq_map, 'Matrix')) {
        stop('shock_eq_map should be of Matrix type')
    } else mod@shock_eq_map = shock_eq_map

    if (!inherits(var_ceq_map, 'Matrix')) {
        stop('var_ceq_map should be of Matrix type')
    } else mod@var_ceq_map = var_ceq_map    

    if (!inherits(cpar_eq_map, 'Matrix')) {
        stop('cpar_eq_map should be of Matrix type')
    } else mod@cpar_eq_map = cpar_eq_map  

    if (!inherits(cpar_ceq_map, 'Matrix')) {
        stop('cpar_ceq_map should be of Matrix type')
    } else mod@cpar_ceq_map = cpar_ceq_map  

    if (!inherits(fpar_eq_map, 'Matrix')) {
        stop('fpar_eq_map should be of Matrix type')
    } else mod@fpar_eq_map = fpar_eq_map  

    if (!inherits(fpar_ceq_map, 'Matrix')) {
        stop('fpar_ceq_map should be of Matrix type')
    } else mod@fpar_ceq_map = fpar_ceq_map  
    
    if (!is.function(ss_function)) {
        stop('ss_function should be a function')
    } else mod@ss_function = ss_function

    if (!is.function(calibr_function)) {
        stop('calibr_function should be a function')
    } else mod@calibr_function = calibr_function

    if (!is.function(ss_calibr_function_jac)) {
        stop('ss_calibr_function_jac should be a function')
    } else mod@ss_calibr_function_jac = ss_calibr_function_jac

    if (!is.function(pert)) {
        stop('pert should be a function')
    } else mod@pert = pert

    if (length(parameters) == length(parameters_free)) {
        mod@is_calibrated <- FALSE
    }
    
    if (!length(shock_eq_map)) 
        mod@is_stochastic = FALSE

    if (mod@is_stochastic)
        mod@active_shocks <- rep(TRUE, length(mod@shocks))
    
    if (any((var_eq_map != 2) & (var_eq_map != 0))) {
        mod@is_dynamic <- TRUE
    }
    
    mod@steady <- rep(0.9, length(mod@variables))
    mod@parameters_calibr_val <-
                rep(0.5, length(mod@parameters_calibr))

    mod@parameters_free_mod_flag <- rep(0,
                        length(mod@parameters_free))
    mod@parameters_calibr_mod_flag <- rep(0,
                        length(mod@parameters_calibr))

    mod@map_free_into_params <- match(mod@parameters_free,
                                        mod@parameters)

    mod@map_calibr_into_params <- match(mod@parameters_calibr,
                                          mod@parameters)

    mod@params[mod@map_free_into_params] <-
        mod@parameters_free_val

    mod@cov_mat <- diag(length(mod@shocks))
    
    return(mod)
}

# ###################################################################
# The is.gecon_model function checks if given object
# is of class "gecon_model"
# ###################################################################
# Input
#   x - any R object
# Output
#   A logical value indicating if object is of class "gecon_model"
# ###################################################################
is.gecon_model <- function(x)
{
    if (is(x, "gecon_model"))
        return(TRUE)
    else return(FALSE)
}

# ###################################################################
# The set_free_par function assigns parameters to gecon_model object
# ###################################################################
# Input
#   model - an object of the gecon_model class, to which
#           we want to assign parameters
#   free_par - a list or vector of parameters with names
#   reset - logical, allows to reset the values parameters to values
#           specified in .gcn file (if set to TRUE)
# Output
#   A gecon_model object with added free parameter values
# ###################################################################
set_free_par <- function(model, free_par = NULL, reset = FALSE)
{     
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    free_par_names <- names(free_par)
    if (is.null(free_par_names) & (reset == FALSE)) 
        stop('free_par argument has to be a list or a vector with named elements')  

    is_proper_numeric <- tryCatch(expr = {free_par <- as.numeric(free_par)
                                                TRUE
                                           },
                                   warning = function(w) FALSE, 
                                   error = function(w) FALSE)
                       
    if (!is_proper_numeric) {
        stop(paste("The following list elements either do NOT contain numeric values or",
                   "or contain more than one entry and cannot be coerced into parameter values:",
                    paste(free_par_names[union(which(as.logical(lapply(free_par, is.character))),
                                         which(as.numeric(lapply(free_par, length)) > 1))], collapse= ', ')))
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
        labs <- rep('Inf', length(misspec_ind))
        labs[which(is.nan(free_par[misspec_ind]))] <- 'NaN'
             
        stop('Values of the following free parameters have been incorrectly specified: \n',
             paste(paste(misspec_names,
                          paste('(= ', labs, ')', sep = '')), 
                   collapse = ', '))
    }        
            
    if (!is.list(free_par) && !(is.numeric(free_par)) && 
        !(is.null(free_par) && reset))  {
        stop(paste("Parameter values have to be passed",
                   "as a numeric vector or a hash list."))
    }
            
    # resetting to the default values
    if (reset) {
        if (all(is.na(model@parameters_free_init_val)))
            stop("Free parameter values have NOT been specified in the .gcn file.")
        matches <- c(1:length(model@parameters_free_val))
        free_par <- model@parameters_free_init_val
        mod_flag <- rep(0, length(model@parameters_free_init_val))
    } else {
        # validating inputs
        not_free_parameters <-
            which(free_par_names %in% model@parameters_calibr)

        if (length(not_free_parameters)) {
            stop(paste('Following parameters',
                       'are calibrated parameters:',
                       paste(free_par_names[not_free_parameters],
                             collapse=' '),
                       'Their values should be specified',
                       'using initval_calibr_par'))
        }
        
        not_parameters <-
            which(!(free_par_names %in% model@parameters_free))
            
        if (length(not_parameters)) {
            stop(paste('Following names of parameters',
                       'have been misspelled or are NOT model parameters:',
                       paste(free_par_names[not_parameters],
                             collapse=' ')))
        }
        
        # initializing
        matches <- which(model@parameters_free %in% free_par_names)
        ord <- match(model@parameters_free[matches], free_par_names)
        ow <- which(matches %in%
                        which(!is.na(model@parameters_free_val)))
        ow <- matches[ow]
        if (length(ow)) {
            warning('Following parameter values have been overwritten:\n',
                    paste(model@parameters_free[ow], collapse=' '))
        }
        free_par <- as.numeric(free_par[ord])
        mod_flag <- 1
    }

    # vector of parameters
    model@parameters_free_val[matches] <- free_par
    model@parameters_free_mod_flag[matches] <- mod_flag
    model@params[model@map_free_into_params[matches]] <- free_par
        
    # clearing slots
    model@is_calibrated = TRUE
    if (length(model@parameters) == length(model@parameters_free)) {
        model@is_calibrated <- FALSE
    }
    model@loglin_var = logical(length = 0)
    model@loglin = TRUE
    model@re_solved = FALSE
    model@corr_computed = FALSE
    model@residual_vector = numeric(0)
    model@solver_status = character(0)

    ss_indic <- FALSE
    if (model@ss_solved) {
        if (length(model@parameters_calibr_val)) {
            ss_indic_res <- abs(model@ss_function(model@steady, 
                                                  model@parameters_calibr_val,
                                                  model@parameters_free_val))
            ss_indic <- (max(ss_indic_res[which(!is.na(ss_indic_res))]) > 1e-6) |
                         any(is.na(ss_indic_res))   
        } else {
            ss_indic_res <- abs(model@ss_function(model@steady, 
                                                  NULL,
                                                  model@parameters_free_val))
            ss_indic <- (max(ss_indic_res[which(!is.na(ss_indic_res))]) > 1e-6) |
                         any(is.na(ss_indic_res))           
        }
    }
    
    if (ss_indic)
        model@ss_solved = FALSE
    
    model@eig_vals = matrix(nrow=0, ncol=0)
    model@solution = list(P = NULL,
                          Q = NULL,
                          R = NULL,
                          S = NULL)
    model@state_var_indices = numeric(0)
    model@solver_exit_info = character(0)
    model@solution_resid = list(NULL)
    model@corr_mat = matrix()
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
            
    return(model)
}

# ###################################################################
# The set_shock_distr_par function assigns distribution parameters
# (standard deviations, correlations of shocks etc)
# to gecon_model object
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   distr_par - a list or vector of distribution parameters.
# Output
#   A gecon_model object with shock distribution parameters set 
#   to values given by the user.
# ###################################################################
set_shock_distr_par <- function(model, distr_par = NULL)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")

    distr_par_names <- names(distr_par)
    if(is.null(distr_par_names)) 
        stop('distr_par argument has to be a list or a vector with named elements')    
    
    is_proper_numeric <- tryCatch(expr = {distr_par <- as.numeric(distr_par)
                                            TRUE
                                       },
                               warning = function(w) FALSE, 
                               error = function(w) FALSE)
                       
    if (!is_proper_numeric) {
        stop(paste("The following list elements either do NOT contain numeric values or",
                   "or contain more than one entry and cannot be coerced into parameter values:",
                    paste(distr_par_names[union(which(as.logical(lapply(distr_par, is.character))),
                                         which(as.numeric(lapply(distr_par, length)) > 1))], collapse= ', ')))
    }

    if (any(is.nan(distr_par) | is.infinite(distr_par) | is.na(distr_par))) {
    
        # NaN and Inf indices
        nan_ind <- which(is.nan(distr_par))
        inf_ind <- which(is.infinite(distr_par))
        na_ind <- which(is.infinite(distr_par))
        
        # names of misspecified parameters
        misspec_ind <- c(nan_ind, inf_ind, na_ind)
        misspec_names <- distr_par_names[misspec_ind]        
        
        # labels for the print
        labs <- rep('Inf', length(misspec_ind))
        labs[which(is.nan(distr_par[misspec_ind]))] <- 'NaN'
        labs[which(is.na(distr_par[misspec_ind]))] <- 'NA'
        
        stop('Values of the following free parameters have been incorrectly specified: \n',
             paste(paste(misspec_names,
                          paste('(= ', labs, ')', sep = '')), 
                   collapse = ', '))
    }
    
    par_type <- rep(0, length(distr_par))
    
    par_type <- grepl("^sd\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) +
                grepl("^var\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 2 +
                grepl("^cor\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*,\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 3 +
                grepl("^cov\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*,\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 4
    
    
    if (any(par_type == 0))
        stop(paste("The following shock distribution parameter names have been misspelled:\n",
                   paste(distr_par_names[which(par_type == 0)], collapse = ', ')))
    
    single_var <- which(par_type %in% c(1, 2))               
    mult_var <- which(par_type %in% c(3, 4))
    
    # validating values
    if (any(distr_par[single_var] < 0))
        stop('Variances and standard deviations cannot be negative')
    
    if (any(distr_par[which(par_type == 3)] > 1 | distr_par[which(par_type == 3)] < -1))
        stop('Correlations have to be numbers from the interval [-1, 1]')    
    
    # retrieving shock names
    shock_names_matrix <- matrix(data = character(), length(distr_par_names), 2)
    
    shock_names <- regexpr('(?<=\\()(_?[\\s,_a-zA-Z0-9])*(?=\\))', distr_par_names, perl = T)
    shock_names <- regmatches(distr_par_names, shock_names)
    shock_names <- gsub('\\s', '', shock_names)
    
    shock_names_matrix[single_var, ] <- shock_names[single_var]
    
    name_vec <- regexpr('(?<=,)(_?[\\s,_a-zA-Z0-9])*', shock_names[mult_var], perl = T)
    shock_names_matrix[mult_var, 1] <- regmatches(shock_names[mult_var], name_vec)

    name_vec <- regexpr('(_?[\\s,_a-zA-Z0-9])*(?=,)', shock_names[mult_var], perl = T)
    shock_names_matrix[mult_var, 2] <- regmatches(shock_names[mult_var], name_vec)
    
    shock_names_matrix <- matrix(match(shock_names_matrix, model@shocks), 
                                 length(distr_par_names), 2)

    # validation of inputs
    wrong_names <- unique(which(is.na(shock_names_matrix), arr.ind=1)[, 1])
    
    if(length(wrong_names))
        stop(paste('The following entries do NOT correspond to any shock distribution parameters:',
                   paste(distr_par_names[wrong_names], collapse=', ')))
    
    complying_ind <- shock_names_matrix[, 1] == shock_names_matrix[, 2]   
    
    corr_err <- which(par_type == 3 & complying_ind)
    var_spec <- which(par_type == 4 & complying_ind)

    if (length(corr_err))
        stop(paste("Correlations have to be specified for different shocks\n",
                   "The specification of:", paste(distr_par_names[corr_err],  collapse = ', '), 
                   "is NOT acceptable"))
    
    if (length(var_spec)) {
        par_type[var_spec] <- 2
        single_var[var_spec] <- TRUE
    }
    
    cov_matrix_flags <- sparseMatrix(i = shock_names_matrix[, 1], 
                                     j = shock_names_matrix[, 2], 
                                     dims = c(length(model@shocks), 
                                              length(model@shocks)))  
                                              
    cov_matrix_flags <- as(cov_matrix_flags, "dgCMatrix")
    check <- tril(cov_matrix_flags, k = -1) + t(triu(cov_matrix_flags, k = 1))
    double_entries <- which(check > 1, arr.ind = T)
    
    if (length(double_entries))
        stop(paste('For the following pairs of shocks the covariances/correlations have been specified more than once:\n',
                    paste(model@shocks[double_entries[, 1]], ' and ', model@shocks[double_entries[, 2]], '\n', collapse = '')))
    duplicates <- duplicated(shock_names_matrix[which(par_type < 3), 1])
    duplicates <- model@shocks[unique(shock_names_matrix[duplicates, 1])]
    if (length(duplicates))
        stop(paste("Variance or standard deviation has been specified more than once for the following shocks:",
                   paste(duplicates, collapse = ', ')))
    
    # updating shock matrix
    var_cov_matrix <- model@cov_mat
    sd_init <- sqrt(diag(var_cov_matrix))
    correlation_matrix <- matrix(0, dim(model@cov_mat)[1], dim(model@cov_mat)[1])
    correlation_matrix[model@active_shocks, model@active_shocks] <- var_cov_matrix[model@active_shocks, model@active_shocks] / 
                            kronecker(sd_init, t(sd_init))[model@active_shocks, model@active_shocks]
    diag(correlation_matrix) <- 1
    
    # updating diagonal of the variance-covariance matrix
    distr_par[which(par_type == 1)] <- distr_par[which(par_type == 1)] ^ 2   
    var_cov_matrix[cbind(shock_names_matrix[which(par_type < 3)], shock_names_matrix[par_type < 3])] <- 
                                                distr_par[which(par_type < 3)]
    sd_fin <- sqrt(diag(var_cov_matrix))
    
    # updating correlation matrix
    col_ind <-  c(shock_names_matrix[which(par_type == 3), 1], 
                      shock_names_matrix[which(par_type == 3), 2])
    row_ind <-  c(shock_names_matrix[which(par_type == 3), 2],
                      shock_names_matrix[which(par_type == 3), 1])                 
    correlation_matrix[(col_ind - 1) * dim(model@cov_mat)[1] + row_ind] <- 
                                               c(distr_par[which(par_type == 3)],
                                                 distr_par[which(par_type == 3)])
    var_cov_matrix <- kronecker(sd_fin, t(sd_fin)) *  correlation_matrix 

                                   
    # updating covariances
    pure_covariances <- which(par_type == 4)

    col_ind <-  shock_names_matrix[pure_covariances, 1]
    row_ind <-  shock_names_matrix[pure_covariances, 2]                    
    
    var_cov_matrix[(col_ind - 1) * dim(model@cov_mat)[1] + row_ind] <- distr_par[pure_covariances]
    var_cov_matrix[(row_ind - 1) * dim(model@cov_mat)[1] + col_ind] <- distr_par[pure_covariances]    
        
    active_shocks <- rep(TRUE, length(model@shocks))
    zero_shocks <- which(diag(var_cov_matrix) == 0)
    if (length(zero_shocks) == length(model@shocks))
        stop("At least one shock has to have variance greater than zero.")
    
    if (length(zero_shocks)) {  
        active_shocks[zero_shocks] <- FALSE 
        var_cov_matrix[zero_shocks, zero_shocks] <- 0
        warning(paste("The following shocks will NOT be taken into account when simulating model: ",
                      paste(model@shocks[zero_shocks], collapse = ', '),
                      "\nAll the covariances/correlations with these shocks have been set to zero."))
    }
    # check if matrix is positive definite
    is_pos_def <- tryCatch(expr = {chol(var_cov_matrix[active_shocks, 
                                                       active_shocks])
                                    TRUE
                                   },
                           warning = function(w) FALSE, 
                           error = function(w) FALSE)
                           
    if (!is_pos_def) {
        stop(paste('Shock matrix is NOT positive definite.',
                   'One or more parameter values have to be modified.'))
    }
    
    model@active_shocks <- active_shocks
    model@cov_mat <- var_cov_matrix
    model@shock_mat_flag <- TRUE
    
    # clearing slots
    model@corr_computed = FALSE
    model@corr_mat = matrix()
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
    
    return(model)
}

# ###################################################################
# The initval_calibr_par function assigns initial values of 
# the calibrating parameters to gecon_model object. 
# Calibrated parameters may be treated as free parameters 
# when the calibration option while solving for the steady state 
# or equilibrium is set to FALSE.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   calibr_par - a list or vector of parameters with names.
# Output
#   A gecon_model object with calibrated parameters set.
# ###################################################################
initval_calibr_par <- function(model, calibr_par)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    calibr_par_names <- names(calibr_par)
    
        
    is_proper_numeric <- tryCatch(expr = {distr_par <- as.numeric(calibr_par)
                                            TRUE
                                       },
                               warning = function(w) FALSE, 
                               error = function(w) FALSE)
                       
    if (!is_proper_numeric) {
        stop(paste("The following list elements either do NOT contain numeric values or",
                   "or contain more than one entry and cannot be coerced into parameter values:",
                    paste(calibr_par_names[union(which(as.logical(lapply(calibr_par, is.character))),
                                         which(as.numeric(lapply(calibr_par, length)) > 1))], collapse= ', ')))
    }
    
    if (is.null(calibr_par_names)) 
        stop('calibr_par argument has to be a list or a vector with named elements')  
        
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
        labs <- rep('Inf', length(misspec_ind))
        labs[which(is.nan(calibr_par[misspec_ind]))] <- 'NaN'
             
        stop('Initial values of the following calibrated parameters have been incorrectly specified: \n',
             paste(paste(misspec_names,
                          paste('(= ', labs, ')', sep = '')), 
                   collapse = ', '))
    }
    
    if (any(is.na(as.numeric(calibr_par)))) {
        na_ind <- which(is.na(as.numeric(calibr_par)))        
        warning('The initial values for the calibrated parameters:\n',
                  paste(calibr_par_names[na_ind], collapse = ' '),
                  '\n  have NOT been modified (NA values have been supplied by the user).')
        calibr_par <- calibr_par[-na_ind]
        if (!length(calibr_par)) {
            warning('NONE of the initial values for calibrated parameters has been modified.')
            return(model)   
        }
    }
    
    # validating inputs
    not_free_parameters <-
        which(calibr_par_names %in% model@parameters_free)
    if (length((not_free_parameters)) != 0) {
        stop(paste('Following parameters',
                   'are free parameters.: \n',
                   paste(calibr_par_names[not_free_parameters],
                         collapse=' '),
                   '\n Their values should be specified',
                   'using set_free_par'))
    }
    
    # validating inputs
    not_parameters <- which(!calibr_par_names %in% model@parameters_calibr)
    if (length((not_parameters)) != 0) {
        stop(paste('Following names of parameters',
                   'have been misspelled or are NOT model parameters:',
                   paste(calibr_par_names[not_parameters],
                         collapse=' ')))
    }
    
    matches <- which(model@parameters_calibr %in% calibr_par_names)
    ord <- match(model@parameters_calibr[matches], calibr_par_names)
    model@parameters_calibr_val[matches] <- as.numeric(calibr_par[ord])
    model@params[model@map_calibr_into_params[matches]] <-
        as.numeric(calibr_par[ord])
    model@parameters_calibr_mod_flag[matches] <- 1
    if (length(matches) == length(model@parameters_calibr)) {
        model@init_calib_pars_supplied <- TRUE
    }

    if (length(matches) != length(model@parameters_calibr)) {
        warning(paste('The following calibrated parameters:\n',
                       paste(model@parameters_calibr[-matches], collapse= ', '),
                      'have retained their previous values',
                      'or have been set to the default value of 0.5'))
    }
    
    
    # clearing slots
    model@is_calibrated = TRUE
    if (length(model@parameters) == length(model@parameters_free)) {
        model@is_calibrated <- FALSE
    }
    model@loglin_var = logical(length = 0)
    model@loglin = TRUE
    model@re_solved = FALSE
    model@corr_computed = FALSE
    model@residual_vector = numeric(0)
    model@solver_status = character(0)
    model@ss_solved = FALSE
    model@eig_vals = matrix(nrow=0, ncol=0)
    model@solution = list(P = NULL,
                          Q = NULL,
                          R = NULL,
                          S = NULL)
    model@state_var_indices = numeric(0)
    model@solver_exit_info = character(0)
    model@solution_resid = list(NULL)
    model@corr_mat = matrix()
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
    
    return(model)

}

# ###################################################################
# The initval_var function sets the initial values of variables
# to values specified by the user.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   init_var -  a list or vector of initial values of variables.
# Output
#   A gecon_model object with nonempty init_var slot.
# ###################################################################
initval_var <- function(model, init_var)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")

    init_var_names <- names(init_var)
    if (is.null(init_var_names)) 
        stop('init_var argument has to be a list or a vector with named elements')  
          
    is_proper_numeric <- tryCatch(expr = {init_var <- as.numeric(init_var)
                                            TRUE
                                       },
                               warning = function(w) FALSE, 
                               error = function(w) FALSE)
                       
    if (!is_proper_numeric) {
        stop(paste("The following list elements either do NOT contain numeric values or",
                   "or contain more than one entry and cannot be coerced into variable values:",
                    paste(init_var_names[union(which(as.logical(lapply(init_var, is.character))),
                                         which(as.numeric(lapply(init_var, length)) > 1))], collapse= ', ')))
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
        labs <- rep('Inf', length(misspec_ind))
        labs[which(is.nan(init_var[misspec_ind]))] <- 'NaN'
             
        stop('Initial values of the following variables have been incorrectly specified: \n',
             paste(paste(misspec_names,
                          paste('(= ', labs, ')', sep = '')), 
                   collapse = ', '))
    }

    if (any(is.na(as.numeric(init_var)))) {
        na_ind <- which(is.na(as.numeric(init_var)))        
        warning('The initial values for the variables:\n',
                  paste(init_var_names[na_ind], collapse = ' '),
                  '\n  have NOT been modified (NA values have been supplied by the user).')
        init_var <- init_var[-na_ind]
        if (!length(init_var)) {
            warning('NONE of the initial values for the model variables has been modified.')
            return(model)   
        }
    }

    if (model@ss_solved)
        warning("New initial values have been set. The existing solution has been erased.")
        
    not_variables <- which(!(init_var_names %in% model@variables))
    if (length((not_variables)) != 0) {
        stop(paste('Following names of variables have been ',
                   'misspelled or are NOT model variables:',
                   paste(init_var_names[not_variables], collapse=' ')))
    }
    
    matches <- which(model@variables %in% init_var_names)
    ord <- match(model@variables[matches], init_var_names)
    model@steady[matches] <- as.numeric(init_var[ord])
    
    if (length(matches) == length(model@variables)) {
        model@init_vals_supplied <- TRUE
    } else {
        warning(paste('The following variables:\n',
                       paste(model@variables[-matches], collapse= ', '),
                      'have retained their previous values',
                      'or have been set to the default value of 0.9'))
        model@init_vals_supplied <- TRUE
    }

    # clearing slots
    model@is_calibrated = TRUE
    if (length(model@parameters) == length(model@parameters_free)) {
        model@is_calibrated <- FALSE
    }
    model@loglin_var = logical(length = 0)
    model@loglin = TRUE
    model@re_solved = FALSE
    model@corr_computed = FALSE
    model@residual_vector = numeric(0)
    model@solver_status = character(0)
    model@ss_solved = FALSE
    model@eig_vals = matrix(nrow=0, ncol=0)
    model@solution = list(P = NULL,
                          Q = NULL,
                          R = NULL,
                          S = NULL)
    model@state_var_indices = numeric(0)
    model@solver_exit_info = character(0)
    model@solution_resid = list(NULL)
    model@corr_mat = matrix()
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
    return(model)
}

# ###################################################################
# The set_shock_cov_mat function sets the shock matrix to matrix specified
# by the user.
# ###################################################################
# Input
#   model - an object of the gecon_model class, for which the 
#              variance covariance matrix of shocks is to be created.
#   shock_matrix -  a positive definite matrix with dimensions 
#                   (n * n) where n is the number of shocks 
#                   in the model.
#   shock_order - a character vector declaring the ordering 
#                 of shocks in the shock_matrix. If not specified,
#                 it is assumed that the ordering is in accordance 
#                 with the internal ordering of model. 
#                 The default ordering can be displayed 
#                 by using the shock_info function with 
#                 the all_shocks argument set to TRUE.
# Output
#   A gecon_model object with the variance covariance matrix set to
#   the values given by the user.
# ###################################################################
set_shock_cov_mat <- function(model, 
                       shock_matrix, 
                       shock_order = NULL)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")    
    
    if ((!is.matrix(shock_matrix) |  !is.numeric(shock_matrix))) {
        stop('Shock matrix has to be passed as a numeric matrix.')
    }

    smd <- dim(shock_matrix)
    if (smd[1] != smd[2]) {
        stop('Dimensions of shock matrix do NOT agree.')
    } else if (smd[1] != length(model@shocks)) {
        stop(paste('Shock matrix dimensions do NOT',
                   'agree with the number of shocks in the model.'))
    }
    
    if(!isSymmetric(shock_matrix))
        stop('Shock matrix has to be symmetric.')
    
    if (length(as.numeric(shock_matrix)) == 1 & (length(model@shocks) == 1)) {
        shock_order <- model@shocks
    }
    
    if (is.null(shock_order)) {
        shock_order <- model@shocks
        warning(paste('The argument shock_order has NOT been specified.',
                'It has been assumed that the ordering of shocks is', 
                'in accordance with internal ordering of model.',
                'This ordering can be displayed using shock_info function with the argument', 
                'all_shocks set to TRUE.'))
    }
    

    if (!is.character(shock_order) | length(shock_order) != length(model@shocks)) {
            stop(paste('The shock_order argument has to be a character vector',
                        'declaring the ordering of ALL the shocks in the shock_matrix.'))
    }    
    
    shock_perm <- match(model@shocks, shock_order)
    
    if (any(is.na(shock_perm)))
        stop('One or more names in shock_order argument is NOT a proper name of shock.')
    
    shock_matrix <- shock_matrix[shock_perm, shock_perm, drop = FALSE]
    
    # shocks with zero variance
    active_shocks <- rep(TRUE, length(model@shocks))
    zero_shocks <- which(diag(shock_matrix) == 0)
    
    if (length(zero_shocks) == length(model@shocks))
        stop("At least one shock has to have variance greater than zero.")
        
    if (length(zero_shocks)) {  
        active_shocks[zero_shocks] <- FALSE 
        shock_matrix[zero_shocks, zero_shocks] <- 0
        warning(paste("The following shocks will NOT be taken into account when simulating model: ",
                      paste(model@shocks[zero_shocks], collapse = ', '),
                      "\nAll the covariances/correlations with these shocks have been set to zero."))
    }
    

    
    is_pos_def <- tryCatch(expr = {chol(shock_matrix[active_shocks, active_shocks])
                                    TRUE
                                  },
                            warning = function(w) FALSE, 
                            error = function(w) FALSE)

    
    if (!is_pos_def) {
        stop(paste('Shock matrix is NOT positive definite.',
                   'Try with another matrix.'))
    }
    
    rownames(shock_matrix) <- model@shocks
    colnames(shock_matrix) <- model@shocks
    model@cov_mat <- shock_matrix
    model@shock_mat_flag <- TRUE
    model@active_shocks <- active_shocks
    
    # clearing slots
    model@corr_computed = FALSE
    model@corr_mat = matrix()
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
    
    return(model)
}

# ###################################################################
# A var_info function prints information about model variables
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   var_names - the names of the variables of interest.
#   all_variables - the logical value. If set to TRUE,
#                   the var_names argument is overwritten 
#                   with all the variables appearing in the model. 
#                   The default value is FALSE.
# Output
#   An object of gecon_var_info class.
# ###################################################################
var_info <- function(model, var_names = NULL,
                     all_variables = FALSE)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    if (is.null(var_names) & !all_variables)
        stop(paste("Variables names have to be specified",
                   "using var_names argument or option",
                   "all_variables has to be set to TRUE",
                   "in order to print the information about variables."))
        
    if (all_variables)
        var_names <- model@variables  
    
    # Validation and getting indices of variables
    if(is(var_names, 'numeric')) {
        if (any(var_names > length(model@variables))) {
            stop(paste('Supplied index is greater than the number of variables.',
                       'Model has:', length(model@variables),
                       'variables.'))
        } else var_ind <- var_names
    } else {
        # Checking if all supplied variables are model variables
        not_in_model <-
            (var_names[which(!(var_names %in% model@variables))])
        if (length(not_in_model))
            stop(paste('Following variables are NOT model variables:',
                       not_in_model, collapse =' '))
        var_ind <- match(var_names, model@variables)
        var_ind <- var_ind[which(!is.na(var_ind))]
    }

    # Creating map of parameters
    if (length(model@var_ceq_map)) {
        var_ceq_map <- model@var_ceq_map
        var_ceq_map[which(var_ceq_map == 1)] <- 8
    } else var_ceq_map <- 
        as(model@var_ceq_map, Class = 'dgCMatrix')
    
    var_eq_map <- rBind(model@var_eq_map,
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
            paste('Equation    ', incid_eq[-c(calibr_eq_numb)])
        rnames[calibr_eq_numb] <- 
            paste('Calibr. Eq. ', incid_eq[calibr_eq_numb])
        rownames(var_eq_map) <- rnames
    } else {
        rownames(var_eq_map) <- paste('Equation    ', incid_eq)
    }

    # SS information
    if (model@ss_solved) {
        # Steady state values
        ss_val <- matrix(model@steady[var_ind], length(var_ind), 1)
        rownames(ss_val) <- model@variables[var_ind]
        if (model@is_dynamic) {
            colnames(ss_val) <- 'Steady state'
        } else {
            colnames(ss_val) <- 'Equilibrium'
        }
    } else {
        # Initial values
        ss_val <- matrix(model@steady[var_ind], length(var_ind), 1)
        rownames(ss_val) <- model@variables[var_ind]
        if (model@is_dynamic) {
            colnames(ss_val) <- 'Initial values for SS computation'
        } else {
            colnames(ss_val) <- 'Initial values for equilibrium computation'
        }
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
        colnames(state_var_impact) <- paste(model@variables[model@state_var_indices], '[-1]', sep='')
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
        steady_state_val <- model@steady[var_ind]
        std_dev_val <- model@sdev[var_ind]

        # Correlations
        cr <- matrix(model@corr_mat[var_ind, ], length(var_ind),
                      dim(model@corr_mat)[2])
        rownames(cr) <- model@variables[var_ind]
        colnames(cr) <- model@variables
    } else {
        std_dev_val <- numeric(length = 0)
        cr <- matrix(nrow = 0, ncol = 0)       
    }
    
    new_var_info <- gecon_var_info(model_info = model@model_info,
                                   model_variable_name = deparse(substitute(model)),
                                   var_names = model@variables[var_ind], 
                                   is_stochastic = model@is_stochastic,
                                   is_dynamic = model@is_dynamic, 
                                   ss_solved = model@ss_solved,
                                   re_solved = model@re_solved,
                                   corr_computed = model@corr_computed, 
                                   ss_val = ss_val, 
                                   state = state,
                                   state_var_impact = state_var_impact, 
                                   shock_impact = shock_impact, 
                                   std_dev_val = std_dev_val,
                                   loglin_flag = loglin_indic,
                                   cr = cr, 
                                   incid_mat = var_eq_map)
    return(new_var_info)
}

# ###################################################################
# A par_info function prints information about model parameters
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   par_names - the names of the parameters of interest.
#   all_parameters - the logical value. If set to TRUE, the par_names
#                    argument is overwritten with all the parameters
#                    appearing in the model. 
#                    The default value is FALSE.      
# Output
#   An object of gecon_par_info class.
# ###################################################################
par_info <- function(model, 
                           par_names = NULL, 
                           all_parameters = FALSE)
{ 
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    if (is.null(par_names) & !all_parameters)
        stop(paste("Parameter names have to be specified",
                   "using par_names argument or option",
                   "all_parameters has to be set to TRUE",
                   "in order to print the information about parameters."))
        
    if (all_parameters)
        par_names <- model@parameters
    
    
    # Checking if all supplied variables are model variables
    not_in_model <-
        par_names[which(!(par_names %in% model@parameters))]
    if (length(not_in_model))
        stop(paste('Following parameters are NOT model parameters:',
                   not_in_model, collapse = ' '))
    
    params_in_map_ind <- match(par_names, model@parameters)
    params_in_map_ind <- params_in_map_ind[which(!is.na(params_in_map_ind))]
    
    calibr_params_in_map_ind <- match(par_names, model@parameters_calibr) 
    calibr_params_in_map_ind <- calibr_params_in_map_ind[which(!is.na(calibr_params_in_map_ind))]
    
    free_params_in_map_ind <- match(par_names, model@parameters_free) 
    free_params_in_map_ind <- free_params_in_map_ind[which(!is.na(free_params_in_map_ind))]
    
    # Creating map of parameters
    cpar_eq_map <- rBind(model@cpar_eq_map,
                         model@cpar_ceq_map)
    
    fpar_eq_map <- rBind(model@fpar_eq_map,
                         model@fpar_ceq_map)
    
    par_eq_map <- sparseMatrix(i = NULL, j = NULL,
                               dims = c(nrow(cpar_eq_map), 
                                        length(params_in_map_ind)))
                                        
    par_eq_map <- as(object = par_eq_map, Class = 'dgCMatrix')
    par_eq_map[, which(par_names %in% model@parameters_calibr)] <- 
                cpar_eq_map[, calibr_params_in_map_ind]
    par_eq_map[, which(par_names %in% model@parameters_free)] <- 
                fpar_eq_map[, free_params_in_map_ind]
    if (is.null(ncol(par_eq_map)))
        par_eq_map <- as(Matrix(par_eq_map), 'dgCMatrix')
    incid_eq <- sort(unique(par_eq_map@i)) + 1
    par_eq_map <- par_eq_map[incid_eq, , drop = FALSE] 

    colnames(par_eq_map) <- model@parameters[params_in_map_ind]
    end_eq_block <- nrow(model@var_eq_map)

    if (any(incid_eq > end_eq_block)) {
        calibr_eq_numb <- which(incid_eq > end_eq_block)
        incid_eq[calibr_eq_numb]  <- incid_eq[calibr_eq_numb] - end_eq_block
        rnames <- vector(length = length(incid_eq))
        rnames[-c(calibr_eq_numb)] <- 
            paste('Equation    ', incid_eq[-c(calibr_eq_numb)])
        rnames[calibr_eq_numb] <- 
            paste('Calibr. Eq. ', incid_eq[calibr_eq_numb])
        rownames(par_eq_map) <- rnames
    } else {
        rownames(par_eq_map) <- paste('Equation    ', incid_eq)
    }
    
    #Finding initial values specified in gcn file
    initial <- rep(NA, length(par_names))
    free_params_ind_in_free <- match(par_names, model@parameters_free)
    free_params_ind_in_free <- 
            free_params_ind_in_free[which(!is.na(free_params_ind_in_free))]
    initial[which(par_names %in% model@parameters_free)]  <- 
                model@parameters_free_init_val[free_params_ind_in_free]    

    calibr_flag <- logical(length(params_in_map_ind))
    calibr_flag[which(par_names %in% model@parameters_calibr)] <- TRUE
    
    info <- gecon_par_info(model_info = model@model_info,
                             model_variable_name = deparse(substitute(model)),
                             par_names = model@parameters[params_in_map_ind],
                             gcn_values = initial,
                             current_values = model@params[params_in_map_ind],
                             calibr_flag = calibr_flag,
                             incid_mat = par_eq_map)
    
    return(info)
}

# ###################################################################
# A shock_info function prints information about model shocks
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   shock_names - the names of the shocks of interest.
#   all_shocks- the logical value. If set to TRUE, the shock_names
#                    argument is overwritten with all the shocks
#                    appearing in the model. 
#                    The default value is FALSE. 
# Output
#   An object of gecon_shock_info class.
# ###################################################################    
shock_info <- function(model, 
                       shock_names = NULL,
                       all_shocks = FALSE)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
        
    if (!model@is_stochastic)
        stop('Model is deterministic. It does not have any shocks.')

    if (is.null(shock_names) & !all_shocks)
        stop(paste("Shocks names have to be specified",
                   "using shock_names argument or option",
                   "all_shocks has to be set to TRUE",
                   "in order to print the information about shocks."))
        
    if (all_shocks)
        shock_names <- model@shocks
    
    
    # Validation and getting indices of shocks
    if(is(shock_names, 'numeric')) {
        if (any(shock_names > length(model@shocks))) {
            stop(paste('Supplied index is greater than the number of variables.',
                       'Model has:', length(model@shocks), 'shocks.'))
        } else shock_ind <- shock_names
    } else {
        # Checking if all supplied shocks are model shocks
        not_in_model <-
            (shock_names[which(!(shock_names %in% model@shocks))])
        if (length(not_in_model))
            stop(paste('Following shocks are NOT model shocks:',
                       not_in_model))
        shock_ind <- match(shock_names, model@shocks)
    }

    shocks_in_map_ind <- match(shock_names, model@shocks)
    shocks_in_map_ind <- shocks_in_map_ind[which(!is.na(shocks_in_map_ind))]
    
    
    # Information about variable incidence
    shock_eq_map <- model@shock_eq_map[, shock_ind, drop = FALSE] 
    incid_eq <- sort(unique(shock_eq_map@i)) + 1
    shock_eq_map <- shock_eq_map[incid_eq, , drop = FALSE]  

    
    rownames(shock_eq_map) <- paste('Eq. ', incid_eq, sep='')
    colnames(shock_eq_map) <- model@shocks[shock_ind]

    # Information about shocks correlation matrix
    cov_mat <-  model@cov_mat[ , shock_ind, drop = FALSE]
    rownames(cov_mat) <- model@shocks
    colnames(cov_mat) <- model@shocks[shock_ind]

    
    info <- gecon_shock_info(model_info = model@model_info,
                             model_variable_name = deparse(substitute(model)),
                             shock_names = model@shocks[shocks_in_map_ind],
                             shock_matrix = cov_mat,
                             shock_matrix_flag = model@shock_mat_flag,
                             incid_mat = shock_eq_map)
                             
    return(info)
}

# ###################################################################
# Function list_eq prints the equations with specified numbers
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   no_eq -  the number of equation(s) which we want to print
# Output
#   prints the equation(s) with given number(s)
# ###################################################################
list_eq <- function(model, no_eq = NULL)
{
    if (!is.gecon_model(model)) {
        stop("The model argument has to be of class gecon_model")
    }

    if (is.null(no_eq)) {
        no_eq = (1:length(model@equations))
    }

    if (any(no_eq < 1 | no_eq > length(model@equations))) {
        stop('Invalid equation index')
    }

    eq <- matrix(model@equations[no_eq], length(no_eq), 1)
    rownames(eq) <- paste('Eq. ', no_eq, ': ', sep = '')
    colnames(eq) <- ''

    print(eq)
    cat('\n')
}

# ###################################################################
# Function list_calibr_eq prints calibration equation(s) with 
# the given number(s)
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   no_eq -  the number of equation(s) which we want to print
# Output
#   prints calibration(s) equation with given number(s)
# ###################################################################
list_calibr_eq <- function(model, no_eq = NULL)
{
    if (!is.gecon_model(model)) {
        stop("The model argument has to be of class gecon_model")
    }

    if (is.null(no_eq)) {
        no_eq = (1:length(model@calibr_equations))
    }

    if (any(no_eq < 1 | no_eq > length(model@calibr_equations))) {
        stop('Invalid calibration equation index')
    }

    eq <- matrix(model@calibr_equations[no_eq], length(no_eq), 1)
    rownames(eq) <- paste('Eq. ', no_eq, ': ', sep = '')
    colnames(eq) <- ''

    print(eq)
    cat('\n')
}

# ###################################################################
# Function round2zero rounds results to zero with given tolerance
# ###################################################################
# Input
#   x - a vector or a matrix
#   tol -  the tolerance level (default 1e-10)
# Output
#   x with small (in absolute value) entries rounded to zero
# ###################################################################
round2zero <- function(x, tol = 1e-8)
{
    ind = which(abs(x) < tol)
    if (length(ind)) {
        x[ind] = 0
    }
    return (x);
}


# ###################################################################
# Function steady_state computes the steady state
# of the dynamic model or equilibrium of the static model
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   solver - the name of nonlinear equations solver.
#             Currently, only:
#               "slv1_nleqslv" - solver from nleqslv package
#             is available
#   use_jac - option to use Jacobian generated by symbolic library
#             (it can resolve plenty of numerical problems),
#             if FALSE numerical derivatives are computed
#   calibration - if TRUE, the calibrating equations are taken
#                 into account while solving for steady state
#                 or equilibrium and parameters are calibrated. 
#                 Otherwise calibrating equations are dropped.
#   options_list - solver settings.
#          Following options are allowed for 'slv1_nleqslv':
#                   to be set in the solver:
#                   -> method - (character) method of solving system
#                               can be set to "Newton" or "Broyden",
#                               default option is "Newton"
#                   -> global - (character) search strategy
#                               can be set to "dbldog", "pwldog",
#                               "qline", "gline", "none"
#                               default option is "qline"
#                   -> xscalm - (character) method of scaling x
#                               can be set to "fixed", "auto",
#                               default option is "fixed"
#                   -> max_iter - (numeric) max. number of
#                               iterations default option is 150
#                   -> tol - (numeric) tolerance of solution,
#                            default is 1e-6
#                   -> xtol - (numeric) relative tolerance of solver step,
#                            default is 1e-6
#   solver_status - if TRUE prints the status of solver
# Output
#   This function sets slot steady to computed value.
#   The solver status is also updated
# ###################################################################
steady_state <- function(model,
                         solver = 'slv1_nleqslv',
                         use_jac = TRUE,
                         calibration = TRUE,
                         options_list = NULL,
                         solver_status = FALSE)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    if (!calibration && 
        (length(model@parameters_free) != length(model@parameters_calibr)))
        model@is_calibrated <- FALSE

    if (any(is.nan(model@parameters_free_val) | 
            is.infinite(model@parameters_free_val))) {
        
        # indices of misspecified parameters
        nan_ind <- which(is.nan(model@parameters_free_val))
        inf_ind <- which(is.infinite(model@parameters_free_val))
        
        # names of misspecified parameters
        nans <- model@parameters_free[nan_ind]
        infs <- model@parameters_free[inf_ind]
                
        misspec <- c(nans, infs)
        
        # labels of misspecified parameters
        len <- (length(nans) + length(infs))
        labs <- rep("Inf", len)      
        if (length(nans))
            labs[1:length(nans)] <- "NaN"

        stop('Values of the following free parameters have been incorrectly specified: \n',
             paste(paste(misspec, paste('(= ', labs, ')', sep = '')), collapse = ', '))
    }

    
    if (any(is.na(model@parameters_free_val))) {
        if (model@is_dynamic) {
            stop('Values of the following free parameters have NOT been specified: \n',
                 paste(model@parameters_free[
                 which(is.na(model@parameters_free_val))], collapse = ' '),
                 '\n All the values of free parameters have to be specified before the steady state is computed.')
        } else {
            stop('Values of the following free parameters have NOT been specified: \n',
                 paste(model@parameters_free[
                 which(is.na(model@parameters_free_val))], collapse = ' '),
                 '\n All the values of free parameters have to be specified before the equilibrium solution is computed.')            
        }
    }

    if (!calibration && any(!model@parameters_calibr_mod_flag)) {
        stop(paste('When calibrated parameter values are NOT determined from calibrating equations',
                   'all of them have to be specified using initval_calibr_par method.',
                   'Following parameters have NOT been specified: \n',
                   paste(model@parameters_calibr[
                         which(!model@parameters_calibr_mod_flag)], 
                         collapse = ' ')))
    }

    if (!calibration && model@init_calib_pars_supplied) {
        if (model@is_dynamic) {
            warning(paste('All parameters will be treated as free parameters and',
                  'calibration equations will NOT be taken into account when determining',
                  'steady state. The initial values of calibration parameters',
                   'will be treated as their values.'))
        } else {
            warning(paste('All parameters will be treated as free parameters and',
                  'calibration equations will NOT be taken into account when determining',
                  'equilibrium. The initial values of calibration parameters',
                   'will be treated as their values.'))       
        }
    }

    # Initial values for solver
    calib_pars <- NULL
    if (model@init_vals_supplied) stead <- model@steady
    calib_pars <- model@parameters_calibr_val
    if (is.numeric(model@parameters_free_val) &&
         !is.na(model@parameters_free_val) &&
         length(model@parameters_free_val))
             free_pars  <- model@parameters_free_val
    if (!(model@init_vals_supplied)) {
        na_stead <- which(is.na(model@steady))
        stead <- model@steady
        stead[na_stead] <- rep(1, length(na_stead))
    }

    # Joining steady_state and calibration functions
    ss_calibr <- function(v, pcalibr = NULL) {
        n1 <- length(v)
        # creating output vector
        if (calibration) {
            n2 <- length(pcalibr)
            sc <- vector(length = (n1 + n2))
        } else {
            n2 <- 0
            sc <- vector(length = (n1)) 
        }
        # assigning values for output vector
        if (n2) {
            sc[1 : n1] <-
                model@ss_function(v, pcalibr, free_pars)
            sc[(n1 + 1) : (n1 + n2)] <-
                model@calibr_function(v, pcalibr, free_pars)
        } else sc[1 : n1] <-
                model@ss_function(v, calib_pars, free_pars)
        return(sc)
    }

    ss_calibr_jac <- function(v, pcalibr = NULL) {
        if (calibration) {
            full <-
                model@ss_calibr_function_jac(v, pcalibr, free_pars)
        } else {
            full <- model@ss_calibr_function_jac(v, 
                                                 calib_pars, free_pars)
        }
        # returning matrix with appropriate dimensions
        if (calibration) {
            return(full[1: (length(v) + length(pcalibr)),
                    1 : (length(v) + length(pcalibr))])
        } else {
            return(full[1: length(v), 1 : length(v)])
        }
    }
    
    # trial evaluation
    ie <- tryCatch(initial_errors <- ss_calibr(stead, calib_pars),
                    warning = function(w) w,
                    error = function(w) w)
    
    # Fast exit in case of NaN's
    if (inherits(ie, 'warning') || inherits(ie, 'error') ) {
        stop('Evaluation of steady-state function has NOT been successful',
               '(indeterminate values): \n', ie, 
                '\n Change initial values in order to find steady_state. \n') 
    }
    
    # Fast exit in case of Infs
    if (any(!is.finite(initial_errors))) {
        stop('Following equations yield non-finite residuals: \n',
                paste(which(!is.finite(initial_errors)), collapse = ' '),
                '\n Change initial values in order to find steady_state. \n') 
    }        
    
    model@init_residual_vector <- initial_errors
      
    out <- ss_par_compute(ss_calibr, 
                          ss_guess = stead,
                          param_guess = calib_pars,
                          calib_par_l = length(calib_pars),
                          calibration = calibration,
                          ss_function_jac = ss_calibr_jac,
                          use_jac = use_jac,
                          solver = solver,
                          options_list = options_list)

    out$ss_values <- round2zero(out$ss_values)
    model@steady <- out$ss_values
    model@solver_status <- out$solver_status
    
    # clearing slots

    model@loglin_var = logical(length = 0)
    model@loglin = TRUE
    model@re_solved = FALSE
    model@eig_vals = matrix(nrow=0, ncol=0)
    model@solution = list(P = NULL,
                          Q = NULL,
                          R = NULL,
                          S = NULL)
    model@state_var_indices = numeric(0)
    model@solver_exit_info = character(0)
    model@solution_resid = list(NULL)
    model@corr_mat = matrix()
    model@corr_computed = FALSE
    model@autocorr_mat = matrix()
    model@corr_variable_mat = matrix()
    model@var_position = numeric(0)
    model@var_dec = matrix()
    model@sdev = matrix() 
    
    if (calibration) {
        if (length(calib_pars)) {
            model@parameters_calibr_val <- out$par_values
            model@params[model@map_calibr_into_params] <- out$par_values
        }
        
        if (length(calib_pars))
            model@residual_vector <-
                ss_calibr(out$ss_values, out$par_values)
        else model@residual_vector <-
                ss_calibr(out$ss_values, NULL)
    } else {
        model@residual_vector <-
            ss_calibr(out$ss_values)
    }

    prec = 1e-6
    tol_index <- which(names(options_list) %in% c('tol'))
    if (length(tol_index) > 0) {
        prec = options_list[tol_index]
    }

    if (any(is.na(model@residual_vector)) || any(!is.finite(model@residual_vector))) {
        nor <- Inf
    } else {
        nor  <- max(abs(model@residual_vector))
    }
    
    if (nor < prec) {
        model@ss_solved <- TRUE
        if (model@is_dynamic) {
            cat('Steady state has been FOUND\n')
        } else {
            cat('The equilibrium has been FOUND\n')
        }
        
        if (solver_status) {
            cat('\n')
            cat('Solver exit info:\n')
            cat(model@solver_status)
        }
    } else {
        if (!(model@init_vals_supplied)) {
            warning(paste('Full vector of initial values',
                          'has NOT been supplied. In missing cases',
                          'default values of',
                          '0.9 have been used.'))
        }
        if (model@is_dynamic) {
            warning(paste('The steady state has NOT been found. 1-Norm of output is: ', nor,
                        ". It is more than requested precision. \n",
                        'Solver exit info: \n', (model@solver_status),
                        sep = ""), 
                        '\n Change initial values or check if the model has the steady state.')
        } else {
            warning(paste('The equilibrium has NOT been FOUND. 1-Norm of output is: ', nor,
                        ". It is more than requested precision. \n",
                        'Solver exit info: \n', (model@solver_status),
                        sep = ""), 
                        '\n Change initial values or check if the model allows for the existence of equilibrium.')        
        }
    }

    # returning
    return(model)
}

# ###################################################################
# The get_residuals function prints the residuals of the 
# steady state or equilibrium equations before the first and 
# after the last solver iteration
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   highest - the name(s) of parameter(s) whose values one wants to check,
#            default option is to print all values
# Output
#   returns vector of residuals
# ###################################################################
get_residuals <- function(model, highest = 5) 
{    
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    output <- list()
    ss_eq_no <- length(model@equations)
    calibr_eq_no <- length(model@calibr_equations)
    if (length(model@init_residual_vector)) {
        irv <-  model@init_residual_vector
        if (length(irv) > ss_eq_no) {
            names(irv) <- c(1:ss_eq_no, paste(1:calibr_eq_no, "calibr"))
        } else names(irv) <- 1:ss_eq_no      
        tops <- min(length(irv), highest)
        cat('\nInitial residuals:\n')
        print(round(irv, digits = 3))
        cat('\nEquations with the largest initial residuals:\n')
        irvo <- names(irv)[as.numeric(order(abs(irv), decreasing=TRUE)[1:tops])]
        cat(paste(as.character(irvo), collapse = ', '), '\n')
        output$initial <- irv
    } else {
        if (model@is_dynamic)
            stop('The steady-state function has NOT been successfully computed yet.')
        else 
            stop('The equilibrium function has NOT been successfully computed yet.')            
    }
    if (length(model@residual_vector)) {
        rv <-  model@residual_vector
        if (length(rv) > ss_eq_no) {
            names(rv) <- c(1:ss_eq_no, paste(1:calibr_eq_no, "calibr"))
        } else names(rv) <- 1:ss_eq_no  
        cat('\nFinal residuals:\n')
        print(round(rv, digits = 3))
        cat('\nEquations with the largest final residuals:\n')
        rvo <- names(rv)[as.numeric(order(abs(rv), decreasing=TRUE)[1:tops])]
        cat(paste(as.character(rvo), collapse = ', '), '\n')
        output$final <- rv
    }
    if (length(model@solution_resid) > 1) {
        cat('\nChecking functions for perturbation solution\n')
        print(model@solution_resid)
    }
    return(invisible(output))
}

# ###################################################################
# The get_par_values function retrieves parameter values
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   par_names - the name(s) of parameter(s) whose values one 
#            wants to check, default option is to print all values
#   to_tex - prints the table into tex file
# Output
#   returns vector of parameter values
# ###################################################################
get_par_values <- function (model, par_names = model@parameters,
                                to_tex = FALSE)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    if(is(par_names, 'numeric')) {
        if (any(par_names > length(model@parameters))) {
            stop(paste('Supplied index is greater than the number of parameters.',
                       'Model has:', length(model@parameters),
                       'parameters.', sep = " "))
        } else par_ind <- par_names
    } else {
        not_in_model <- par_names[which(!(par_names  %in%  model@parameters))]

        if (length(not_in_model)) {
            stop(paste('Following parameters are NOT',
                       'parameters of this model: ', not_in_model))
        } else
            par_ind <- match(par_names, model@parameters)
    }

    par_vector <- model@params[par_ind]
    par_names <- model@parameters[par_ind]

    # Printing parameter information
    par_val <- matrix(par_vector, length(par_ind), 1)
    rownames(par_val) <- par_names
    colnames(par_val) <- 'Parameters'


    if (to_tex) {
        to_tex(par_val, 
               paste(model@model_info[2], '.results.tex', sep=''),
               colnam = "Parameters",
               section = "Parameters of the model \n") 
    } else {
        cat('\nParameters of the model:\n\n')
        print(round(par_val, digits = 3))
    }
    # Returning parameter values
    names(par_vector) <- par_names
    invisible(par_vector)
}

# ###################################################################
# The get_par_names function returns names of free
# parameters stored in an object of gecon_model class.
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   free_par - logical, controls if free parameters should be
#              added to the vector of parameter names.
#   calibr_par - logical, controls if calibrated parameters should be
#                added to the vector of parameter names.
# Output
#   returns vector of parameter names chosen by the user.
# ###################################################################
get_par_names <- function(model, free_par = TRUE, calibr_par = TRUE)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    if (!free_par & !calibr_par)
       stop("At least one of the free_par and calibr_par options has to be set to TRUE.")
       
    ret_list <- character(0)
    
    if (calibr_par)
        ret_list <- model@parameters_calibr 

    if (free_par)
        ret_list <- c(ret_list, model@parameters_free)

    return(ret_list)
}

# ###################################################################
# The get_var_names function returns names of variables 
# stored in an object of gecon_model class
# ###################################################################
# Input
#   model - an object of the gecon_model class
# Output
#   returns a vector of variable names
# ###################################################################
get_var_names <- function(model)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
           
    return(model@variables)
}

# ###################################################################
# The get_shock_names function returns names of shocks 
# stored in an object of gecon_model class
# ###################################################################
# Input
#   model - an object of the gecon_model class
# Output
#   returns a vector of shock names
# ###################################################################
get_shock_names <- function(model)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")   
       
    return(model@shocks)
}

# ###################################################################
# The get_ss_values function prints and returns steady state 
# (equilibrium) values of model variables 
# in case of dynamic (static) model
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   var_names - names or indices of the variables, whose steady state
#               (or equilibrium) values are to be returned, 
#               default option is a vector of all the variables.
#   to_tex - if TRUE, prints the table into tex file.
#   silent - logical. The default value is FALSE. 
#            If set to TRUE, it suppresses console output. 
# Output
#   returns a vector of steady state (equilibrium) values.
# ###################################################################
get_ss_values <- function(model, var_names = NULL, 
                          to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    if (is.null(var_names))
        var_names = model@variables
        
    # Validation and getting indices of variables
    if (!model@ss_solved) {
        if (model@is_dynamic) 
            stop("Compute steady state first using 'steady_state' function.")
        else
            stop("Compute equilibrium first using 'steady_state' function.")            
    }
    if(is(var_names, 'numeric')) {
        if (any(var_names > length(model@variables))) {
            stop(paste('Supplied index is greater than the number of variables.',
                       'Model has:', length(model@variables),
                       'variables.', sep = " "))
        } else var_ind <- var_names
    } else {
        # Checking if all supplied variables are model variables
        not_in_model <-
            (var_names[which(!(var_names %in% model@variables))])
        if (length(not_in_model))
            stop(paste('Following variables are NOT model variables:',
                       paste(not_in_model, collapse = ' '), sep = " "))
        var_ind <- which(model@variables %in% var_names)
    }

    ss_vector <- model@steady[var_ind]
    ss_names <- model@variables[var_ind]

    # Printing SS information
    ss_val <- matrix(ss_vector, length(var_ind), 1)
    rownames(ss_val) <- ss_names
    if (model@is_dynamic) {
        colnames(ss_val) <- 'Steady-state'
    } else {
        colnames(ss_val) <- 'Equilibrium'  
    }

    if (to_tex) {
        if (model@is_dynamic) {
            to_tex(ss_val, 
                   paste(model@model_info[2], '.results.tex', sep=''),
                   colnam = "Steady-state values",
                   section = "Steady-state values \n") 
        } else {
            to_tex(ss_val, 
                   paste(model@model_info[2], '.results.tex', sep=''),
                   colnam = "Equilibrium values",
                   section = "Equilibrium values \n")         
        }
    } else {
        if (!silent) {
            if (model@is_dynamic) {
                cat('\nSteady-state values:\n\n')
            } else {
                cat('\nEquilibrium values:\n\n')    
            }
            print(round(ss_val, digits = 4)) 
        }
    }
    # Returning SS values
    names(ss_vector) <- ss_names
    invisible(ss_vector)
}


# ###################################################################
# The solve_pert function (log-)linearises model, removes zero rows,
#   solves the first order perturbation
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   loglin - a logical option. If TRUE, it allows to log-linearise 
#               perturbation. If FALSE, model is linearised only,
#   not_loglin_var - vector of variables that will NOT
#                    be log-linearised, valid only if 
#                    loglin is set to TRUE
#   norm_tol = 1e-8,
#   solver - first order perturbation solver, 
#            the default solver is Christopher Sims solver gensys.
#            Currently no other solvers available
#
# Output
#   Solves the model and sets the slots where the solution 
#   or solver information is stored to new values
# ###################################################################
solve_pert <- function(model, loglin = TRUE, not_loglin_var = NULL,
                       norm_tol = 1e-8, solver = 'sims_solver')
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    if (!model@is_dynamic) stop('Model is NOT dynamic')

    model@loglin <- loglin
    model@loglin_var <- rep(loglin, length(model@variables))
    
    if ((!is.gecon_model(model))) 
        stop("Argument has to be of class gecon_model")    
    
    if ((!model@ss_solved)) {
         stop(paste('In order to obtain perturbation',
                    'you have to compute steady-state values and parameters first'))        
    }
        
    # Linearisation or log-linearisation of a model
    if (!model@loglin) {
        abcd <- model@pert(model@steady,
                           model@parameters_calibr_val,
                           model@parameters_free_val)
        abcd <- row_reduction(abcd)
    } else if (length(not_loglin_var) == 0) {
        abcd <- model@pert(model@steady,
                           model@parameters_calibr_val,
                           model@parameters_free_val)
        abcd <- row_reduction(abcd)
        abcd <- lintolog(abcd, model@steady)
    } else {
        abcd <- model@pert(model@steady,
                           model@parameters_calibr_val,
                           model@parameters_free_val)
        abcd <- row_reduction(abcd)
        var_loglin <- model@steady
        ind <- which(model@variables  %in%  not_loglin_var)
        var_loglin[ind] <- 1
        model@loglin_var[ind] <- FALSE
        abcd <- lintolog(abcd, var_loglin)
    }
    # Solution of model
    solve_output <- switch(solver,
                           sims_solver = sims_solver(abcd,
                                is_stochastic = model@is_stochastic),
                           stop(paste('Solver ', solver,
                               'has NOT been implemented')))

    if (model@is_stochastic) {
        shock_n <- dim(solve_output$Q)[2]
        shock_indices <- c(1:shock_n)
    }

    # Names of rows and columns of matrices P, Q, R, S
    if (length(solve_output$state_var_indices)) {
        # P
        colnames(solve_output$P) <-
            paste(model@variables[solve_output$state_var_indices],
                  '[-1]', sep = "")

        rownames(solve_output$P) <-
            model@variables[solve_output$state_var_indices]

        # R
        colnames(solve_output$R) <-
            paste(model@variables[solve_output$state_var_indices],
                  '[-1]', sep = "")

        rownames(solve_output$R) <-
            model@variables[-c(solve_output$state_var_indices)]
    }

    if (model@is_stochastic) {
        if (solve_output$norm_stoch < norm_tol) {
            # Q
            if (model@is_stochastic) {
                colnames(solve_output$Q) <- model@shocks
                rownames(solve_output$Q) <-
                    model@variables[solve_output$state_var_indices]
            }
            
            # S
            if (model@is_stochastic) {
                rownames(solve_output$S) <-
                    model@variables[-c(solve_output$state_var_indices)]
                colnames(solve_output$S) <-  model@shocks
            }
            
        }
    }
    
    # assigning eigenvalues
    model@eig_vals <- matrix(solve_output$eig, ncol=3)
    if (length(model@eig_vals > 2)) {
        colnames(model@eig_vals) <- c('Mod', 'Re', 'Im')
    }

    # Check if model was solved
    if (model@is_stochastic) {
        if ((solve_output$slv_exit_code[1] == 1) &
                    (solve_output$slv_exit_code[2] == 1)) {
            # solution info                      
            cat('Model has been SOLVED \n')
            
            if ((solve_output$norm_deter > norm_tol | 
                 solve_output$norm_stoch > norm_tol)) {
                warning(paste('\n',
                              'The solution has been found by the solver',
                              'but it does NOT satisfy perturbation equations.\n',
                              'Norm of residuals in the deterministic part=',
                               solve_output$norm_deter, '\n',
                              'Norm of residuals in the stochastic part=',
                               solve_output$norm_stoch, '\n', 
                              'This discrepancy may be due to numerical roundoff errors.',
                              'Consider changing parametrisation of the model',
                              'or loglinearising variables with large steady state values.\n','\n'))
            }
            if (!length(solve_output$state_var_indices)) {
                warning('The model does NOT have any lagged variables. \n',
                        'There are no state variables. Only S matrix of the solution is not empty.')
            }
            # Passing solution to slots
            model@solution = list(P = round2zero(solve_output$P),
                                  Q = round2zero(solve_output$Q),
                                  R = round2zero(solve_output$R),
                                  S = round2zero(solve_output$S))
            model@state_var_indices = solve_output$state_var_indices
            model@solver_exit_info = solve_output$slv_exit_info
            model@solution_resid = list(norm_deter = solve_output$norm_deter,
                                        norm_stoch = solve_output$norm_stoch,
                                        norm_deter_prim = solve_output$norm_deter_prim,
                                        norm_stoch_prim = solve_output$norm_stoch_prim)

            model@re_solved <- TRUE
            
            # clearing slots
            model@corr_computed <- FALSE
            model@corr_mat = matrix()
            model@autocorr_mat = matrix()
            model@corr_variable_mat = matrix()
            model@var_position = numeric(0)
            model@var_dec = matrix()
            model@sdev = matrix() 
        } else {
            # solution info
            cat('Model has NOT been SOLVED\n')
            
            warning(paste('\n',
                          'Solver exit info:\n',
                          solve_output$slv_exit_info,
                          '\n',
                          'Norm of residuals in the deterministic part =',
                           solve_output$norm_deter, '\n',
                          'Norm of residuals in the stochastic part =',
                           solve_output$norm_stoch, '\n', '\n'))

            model@solver_exit_info = solve_output$slv_exit_info
        }
    } else {
        if ((solve_output$slv_exit_code[1] == 1) &
                    (solve_output$slv_exit_code[2] == 1)) {
            # solution info                            
            cat('Model has been SOLVED\n')
            
            if ((solve_output$norm_deter > norm_tol)) {
                warning(paste('\n',
                              'The solution has been found by the solver',
                              'but it does NOT satisfy perturbation equations.\n',
                              'Norm of residuals in the deterministic part = ',
                              solve_output$norm_deter, '\n',
                              'This discrepancy may be due to numerical roundoff errors.',
                              'Consider changing parametrisation of the model',
                              'or loglinearising variables with large steady state values.\n', '\n'))
            }
            
            if (!length(solve_output$state_var_indices)) {
                warning('The model does not have any lagged variables.\n',
                        'There are no state variables. All the solution matrices are empty.')
            }
            
            # Passing solution to slots
            model@solution = list(P = round2zero(solve_output$P),
                                  R = round2zero(solve_output$R))
            model@state_var_indices = solve_output$state_var_indices
            model@solver_exit_info = solve_output$slv_exit_info
            model@solution_resid =
                    list(norm_deter = solve_output$norm_deter,
                         norm_deter_prim = solve_output$norm_deter_prim)
            model@re_solved <- TRUE
        } else {
            # solution info
            cat('Model has NOT been SOLVED\n')
            
            if ((solve_output$norm_deter > norm_tol)) {
                warning(paste('\n',
                              'Solver exit info:',
                              '\n',
                              solve_output$slv_exit_info, '\n',
                              'Norm of residuals in the deterministic part = ',
                              solve_output$norm_deter))
            }
            model@solver_exit_info = solve_output$slv_exit_info
        }
    }

    return(model)
}


# ###################################################################
# The get_pert_solution function retrieves the solution of perturbation
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   to_tex - prints an output to tex file
#   silent - TODO
# Output
#   returns a list of P, Q, R, S matrices and prints them on the console
#       (does not return Q and S when the model is not stochastic)
# ###################################################################
get_pert_solution <- function(model, to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model))
       stop("The model argument has to be of class gecon_model")
    
    # Validation
    if (!model@is_dynamic)
        stop("Model is NOT dynamic")
        
    if (!model@re_solved)
        stop(paste("The model in its (log-)linearised form has NOT been SOLVED yet.",
                   "Solve it first using solve_pert function.", sep = " "))

    P <- model@solution$P
    Q <- model@solution$Q
    R <- model@solution$R
    S <- model@solution$S

    returns_stoch <- list(P = P, Q = Q, R = R, S = S)
    returns_deter <- list(P = P, R = R)
    
    tex_state_names = paste(string_to_tex(model@variables[model@state_var_indices], math_mode = FALSE),
                         '_{t-1}', sep='')

    tex_shock_names = paste(string_to_tex(model@shocks, math_mode = FALSE), sep='')                         
    
    if (model@is_stochastic) {      
        # tex output
        if (to_tex) {
            to_tex(returns_stoch, 
                   paste(model@model_info[2], '.results.tex', sep=''),
                   colnam = list(tex_state_names, tex_shock_names, tex_state_names, tex_shock_names),
                   section = "The solution of the perturbation",
                   to_table = FALSE) 
        } else {
            if (!silent) {
                cat('\nMatrix P:\n\n')
                print(round(P, digits = 4))
                
                cat('\n\nMatrix Q:\n\n')
                print(round(Q, digits = 4))
                
                cat('\n\nMatrix R:\n\n')
                print(round(R, digits = 4))
                
                cat('\n\nMatrix S:\n\n')
                print(round(S, digits = 4))
            }
        }
        invisible(returns_stoch)
    } else { 
        # tex output
        if (to_tex) {
            to_tex(returns_deter, 
                   paste(model@model_info[2], '.results.tex', sep=''),
                   colnam = list(tex_state_names, tex_state_names),
                   section = "The solution of the perturbation",
                   to_table = FALSE) 
        } else {
            cat('\nMatrix P:\n\n')
            print(round(P, digits = 4))
            
            cat('\n\nMatrix R:\n\n')
            print(round(R, digits = 4))
        }
        invisible(returns_deter)
    }
}


# ###################################################################
# The check_bk function checks the Blanchard Kahn conditions
# ###################################################################
# Input
#   model - an object of the gecon_model class
# Output
#   Checks the Blanchard Kahn conditions and prints info about eigenvalues
#   larger than 1 in modulus and the number of forward looking variables
# ###################################################################
check_bk <- function(model)
{
    if (is.gecon_model(model)) {
        if (!length(model@eig_vals)) {
            stop('The function solve_pert has to be called first.')
        }
        cat('\nEigenvalues of system:\n')
        print(model@eig_vals)

        # number of forward looking variables
        fl <- unique(which(as.matrix(model@var_eq_map) > 3, 
                           arr.ind=TRUE)[,2])
        cat(paste('\nThere are:', length(fl), 'forward looking variables.',
                  'There are:', length(which(model@eig_vals[, 1] > 1)),
                    'eigenvalues larger than 1 in modulus \n'))
        if (length(fl) ==  length(which(model@eig_vals[, 1] > 1)))
            cat('BK conditions have been SATISFIED \n')
        else cat('BK conditions have NOT been SATISFIED \n')
    }
    else stop("The model argument has to be of class gecon_model")
}

# ###################################################################
# The random_path function simulates series of shocks based 
# on the covariance matrix and computes the behaviour of the specified
# variables in the system
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   shock_list - a list with names of shocks that should be
#              taken into account when computing shocks path
#   var_list - a list of variables for which 
#                   the impact has to be computed
#   path_length = the length of stochastic path, default=100
# Output
#   An object of class gecon_simulation
# ###################################################################
random_path <- function(model, shock_list = NULL, 
                        var_list = NULL, path_length = 100)
{
    if (!is.gecon_model(model))
        stop("Random path can be generated only for gecon_model object")

    if (!model@shock_mat_flag)
        warning(paste("Default identity matrix has been used for stochastic simulation.",
                      "\nUser-defined variance-covariance matrix can be set using set_shock_cov_mat or set_shock_distr_par functions."))
    
    if (!model@is_stochastic)
        stop("Random path can be generated only for a stochastic model")
    
    epsns <- which(model@active_shocks)

    if (!is.null(shock_list)) {
        epsns <- match(shock_list, model@shocks)
        if (any(is.na(epsns)))
            stop("Following shock names: \n",
                 shock_list[which(is.na(epsns))],
                 "\n are NOT valid model shock names")
    }
    
    if (length(intersect(epsns, which(!model@active_shocks)))) {
        if (!length(setdiff(epsns, which(!model@active_shocks))))
            stop('To perform random path simulation, at least one active shock has to be specified')
        warning(paste('The following shocks are NOT active and will NOT be used for the random path simulation:',
                  paste(model@shocks[intersect(epsns, which(!model@active_shocks))], collapse = ', ')))
        epsns <- setdiff(epsns, which(!model@active_shocks))
    }
    
    eps_length <- length(epsns)
    sim_mat <- t(chol(model@cov_mat[epsns, epsns]))
    n_shock <- eps_length

    shock_path <- matrix(0, ncol= path_length, nrow=n_shock)
    random_vec <- rnorm(path_length * length(epsns), mean = 0, sd = 1)
    
    model_name <- deparse(substitute(model))

    for (i in 1:path_length) {
        shock_path[(1:length(epsns)), i] <- 
            sim_mat %*% random_vec[(((i-1) * length(epsns) + 1):(length(epsns)*i))]
    }
    
    simulate_model(model, shock_m = shock_path, periods= c(1:path_length), 
                   path_length = path_length,
                   shock_list = model@shocks[epsns], var_list = var_list,
                   sim_type = "Random path simulation",
                   model_name = model_name)
}

# ###################################################################
# The compute_irf function computes impulse response functions
# ###################################################################
# Input
#   model - an object of class gecon_model
#   shock_list - a list of shocks for which irfs are to be computed,
#              if the argument is missing, irfs are computed for all
#              the shocks in model
#   var_list - a list of variables for which impact of shocks has
#              to be computed
#   path_length - the number of periods for which the irfs are computed
#   cholesky - computes irfs based on uncorrelated shocks
#              when FALSE, and irfs based on correlated shocks
#              when TRUE, default value is FALSE
# Output
#   An object of class gecon_simulation
# ###################################################################
compute_irf <- function(model,  shock_list = NULL, var_list = NULL, 
                        path_length = 40, cholesky = TRUE)
{
    if (!is.gecon_model(model))
        stop("IRFs can be generated only for gecon_model object")

    if (!model@shock_mat_flag)
        warning(paste("Default identity matrix has been used for stochastic simulation.",
                      "\nUser-defined variance-covariance matrix can be set using set_shock_cov_mat or set_shock_distr_par functions."))

    if (!model@is_stochastic)
        stop("There are no shocks in the model")
    
    epsns <- which(model@active_shocks)

    if (!is.null(shock_list)) {
        epsns <- match(shock_list, model@shocks)
        if (any(is.na(epsns)))
            stop("Following shock names: \n",
                 shock_list[which(is.na(epsns))],
                 "\n are NOT valid model shock names")
    }
    
    if (length(intersect(epsns, which(!model@active_shocks)))) {
        if (!length(setdiff(epsns, which(!model@active_shocks))))
            stop('To compute the IRFs, at least one active shock has to be specified')
        warning(paste('The following shocks are NOT active and will NOT be used for the IRFs computation:',
                  paste(model@shocks[intersect(epsns, which(!model@active_shocks))], collapse = ', ')))
        epsns <- setdiff(epsns, which(!model@active_shocks))
    }
    

    shock_path <- rep(0, sum(model@active_shocks) * path_length *
                         length(epsns))
    dim(shock_path) <- c(sum(model@active_shocks), path_length,
                         length(epsns))
    epsns2 <- match(epsns, which(model@active_shocks))
    
    for (i in (1:length(epsns2))) {
        eps <- rep(0, sum(model@active_shocks))
        eps[epsns2[i]] <- 1
        
        if (is.null(model@cov_mat)) {
            corr_mat <- diag(length(eps))
        } else 
            corr_mat <- model@cov_mat[model@active_shocks, model@active_shocks]
        
        # whether to use cholesky decomposition of shocks
        if (cholesky) {
            eps <- chol(corr_mat) %*% eps 
        }
        shock_path[, 1, i] <- eps
    }    
    
    shock_path <- list(shock_path = shock_path,
                       shock_list   = model@shocks[epsns])
    model_name <- deparse(substitute(model))

    simulate_model(model, 
                   shock_m = shock_path, 
                   periods = c(1:path_length), 
                   path_length = path_length,
                   shock_list = model@shocks[model@active_shocks], 
                   var_list = var_list,
                   sim_type = "Impulse response functions",
                   model_name = model_name)
}

# ###################################################################
# The simulate_model function allows user to specify shock
# realizations and to simulate the economy based on them
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   shock_list - a list of shock names. The number of shocks
#               has to be equal to the number of 
#               shock_m matrix rows
#   var_list - a list of variables for which the impact has 
#               to be computed
#   shock_m - a matrix of shocks specified by the user
#               When random_path and simulate_model are called
#               simulate_model accepts also a list
#   periods - periods for which user specified shocks in shock_m,
#               default values are from 1 to number of columns of 
#               shock matrix
#   path_length - number of periods for which the system is simulated
#   sim_type - the type of conducted simulation, according
#               to the type of the function called
#               -> "Impulse response functions" by compute_irf
#               -> "Random path simulation" by random_path
#               -> "Simulation with user-defined shocks" when 
#                       the function simulate_model is called directly
#   model_name - the name of the variable storing model
# Output
#   An object of class gecon_simulation
# ###################################################################
simulate_model <- function(model, shock_list = NULL, var_list = NULL,
                           shock_m = NULL, periods = NULL, 
                           path_length = 40, 
                           sim_type = NULL,
                           model_name = NULL)
{
    if (!is.gecon_model(model))
        stop("simulate_model function requires gecon_model object as an argument.")
    
    if (is.null(shock_m))
        stop("Shocks have to be specified in order to simulate model.")
    
    if (!model@is_stochastic)
        stop("There are no shocks in the model")
    
    if (!length(model@state_var_indices)) {
        stop('Model does not have any state variables. \n',
             'Shocks do not have impact on variables in any periods but t.')
    }
    
    if (is.null(sim_type)) {
        sim_type <- "Simulation with user-defined shocks"
    }
        
    if (is.null(shock_list))
        shock_list <- model@shocks[1:dim(shock_m)[1]]
      
    if (is.list(shock_m)) {
        eps_names <- shock_m[[2]]
        shock_m <- shock_m[[1]]
    } else eps_names <- shock_list
    
      
    if (dim(shock_m)[1] != length(shock_list)) 
        stop('Number of shock names is different than the number of matrix rows.')
    
    if (is.null(periods))
        periods <- 1:dim(shock_m)[2]
        
    if (dim(shock_m)[2] != length(periods))
        stop('Number of periods is different than the number of matrix columns.')
       
    if(is.null(model_name))
        model_name <-  deparse(substitute(model))
        
    # Creating shock matrix
    n_shock <- length(model@shocks)        
    shock_ind <- match(shock_list, model@shocks)
    
    if (any(is.na(shock_ind)))
        stop("Shocks: \n", 
             shock_list[which(is.na(shock_ind))],
             "\n are not model shocks")

    if ((sim_type == "Simulation with user-defined shocks") & 
        length(intersect(shock_ind, which(!model@active_shocks)))) {
        eps_names <- eps_names[-which(is.na(match(shock_ind, model@active_shocks)))]
        shock_m <- shock_m[-which(is.na(match(shock_ind, model@active_shocks))), ]
        if (!length(setdiff(shock_ind, which(!model@active_shocks))))
            stop('To perform any simulations, at least one active shock has to be specified.')
        warning(paste('The following shocks are NOT active and will NOT be used for the simulations:',
                  paste(model@shocks[intersect(shock_ind, which(!model@active_shocks))], collapse = ', ')))
        shock_ind <- setdiff(shock_ind, which(!model@active_shocks))
    }
         

    
             
    total_periods <- c(1:path_length)
    period_ind <- match(periods, total_periods)    
    if (any(is.na(period_ind)))
        stop("One or more periods are not valid periods. \n",
             "Check if period numbers are positive integers smaller than path length")
             

    # Creating indices
    # selecting shock indices
    epsns <- model@shocks[c(1: length(model@shocks))]
    if (!is.null(shock_list)) {
        epsns <- model@shocks[which(model@shocks %in%
                                            shock_list)]
        miss_e <- shock_list[-which(shock_list %in% model@shocks)]
        if (!length(which(shock_list %in% model@shocks)))
            stop("None of specified shocks are model shocks")
        if (length(miss_e)) {
            warning('Following shocks are not model shocks: ',
                    paste(miss_e, collapse= ' '), 
                    '\n Simulations will not be performed for these shocks. \n') 
        }                                                        
    } 

    # selecting variables indices 
    vars <- model@variables[c(1: length(model@variables))]
    if (!is.null(var_list)) {
        vars <-  model@variables[which(model@variables %in%
                                         var_list)]
        miss <- var_list[-which(var_list %in% model@variables)]
        if (!length(which(var_list %in% model@variables)))
            stop("None of specified variables are model variables")
        if (length(miss)) {
            warning('Following variables are not model variables: ',
                    paste(miss, collapse= ' '), 
                    '\n Simulations will not be performed for these variables. \n') 
        }
    } else vars <-model@variables[model@state_var_indices]
    
    # variables that are not state variables
    jump <- which(!(vars %in%
                        model@variables[model@state_var_indices]))
    jump <-
        which(model@variables[-c(model@state_var_indices)]
              %in% vars[jump])
    jump_names <-
        model@variables[-c(model@state_var_indices)]
    jump_names <- jump_names[jump]
    
    # variables that are state variables
    stat <-  which((vars %in%
                        model@variables[model@state_var_indices]))
    stat <-
        which(model@variables[model@state_var_indices] %in%
                  vars[stat])
    stat_names <-
        model@variables[(model@state_var_indices)]
    stat_names <- stat_names[stat]

    
    variables_names <- c(stat_names, jump_names)
    # loop for computing IRFs for all passed eps
    if (length(dim(shock_m)) < 3) {
        epsilons = 1
    }    
    else 
        epsilons <- dim(shock_m)[3]
        
    mat_dim <- c(length(vars), path_length, epsilons)
    
    irfs <- rep(0, (length(vars)* path_length * epsilons))
    dim(irfs) <- mat_dim
    
        
    shock_path <- matrix(0, nrow = n_shock, ncol = path_length)
    if (length(dim(shock_m)) < 3) {
        shock_path[shock_ind, period_ind] <- shock_m
        irfs[1:length(vars), 1:path_length, 1] <-  matrix(simulate_path(P = model@solution$P,
                                                               Q = model@solution$Q,
                                                               R = model@solution$R,
                                                               S = model@solution$S,
                                                               tn = path_length, 
                                                               init_st = NULL,
                                                               shock_path = shock_path,  # number of shock - to be removed
                                                               jump = jump,              # variables to be simulated
                                                               stat = stat), 
                                                               nrow=length(vars), ncol = path_length) 

    } else {
        for (i in 1:epsilons) {
            shock_path[shock_ind, period_ind] <- as.matrix(shock_m[ , ,i])
            irfs[1:length(vars), 1:path_length, i] <- matrix(simulate_path(P = model@solution$P,
                                                                   Q = model@solution$Q,
                                                                   R = model@solution$R,
                                                                   S = model@solution$S,
                                                                   tn = path_length, 
                                                                   init_st = NULL,
                                                                   shock_path = shock_path,  # number of shock - to be removed
                                                                   jump = jump,              # variables to be simulated
                                                                   stat = stat), nrow=length(vars), ncol = path_length)            
        }
    }
    
    sim_res <- gecon_simulation(irfs, 
                                 eps_names, 
                                 variables_names, 
                                 sim_type,
                                 path_length, 
                                 model_info = model@model_info,
                                 model_variable_name = model_name) 
    
    return(sim_res)
}

# ###################################################################
# The compute_moments function computes correlation matrices based on
# spectral and simulation methods, with HP filtered or not
# HP filtered series
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   ngrid - grid for FFT computation, has to be sufficiently large
#               to compute unbiased correlations
#   filter - logical, if TRUE correlations are computed for HP-filtered variables
#   sim - logical, if TRUE simulation methods are used for correlations computations
#   nrun - the number of MC simulation runs
#   lambda - lambda for HP filter
#   ref_var - name or the index of variable in relation to which
#            correlations are computed
#   n_leadlags - the number of leads/lags for computing correlation
#                tables
# Output
#   Computes correlations, assigns them to specific slot of the gecon_model
#   class and prints them
# ###################################################################
compute_moments <- function(model, ngrid = (64 * 16),
                         filter = TRUE, sim = FALSE,
                         nrun = 100000, lambda = 1600,
                         ref_var = NULL, n_leadlags = 6)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")

    # Checking if shock matrix is set
    if (!model@shock_mat_flag)
        warning(paste("Default identity matrix has been used for stochastic simulation.",
                      "\nUser-defined variance-covariance matrix can be set using set_shock_cov_mat or set_shock_distr_par functions."))

    # Checking if shock matrix is set
    if (!model@shock_mat_flag)
        warning(paste("Default identity matrix has been used for stochastic simulation.",
                      "\nUser-defined variance-covariance matrix can be set using set_shock_cov_mat or set_shock_distr_par functions."))

    # Checking if model is solved
    if (!model@re_solved) {
        stop(paste('Model is not solved.',
             'Solve it first to obtain recursive laws of motion.'))
    }         
    
    if (!length(model@state_var_indices)) {
        stop('Model does not have any state variables. \n',
             'Shocks do not have impact on variables in any periods but t.')
    }
        
    # Creating PP , QQ matrices
    PP <- matrix(0, length(model@variables),
                 length(model@variables))
    QQ <- matrix(0, length(model@variables),
                 sum(model@active_shocks))
                 
    PP[model@state_var_indices,
       model@state_var_indices] <- model@solution$P
    PP[-c(model@state_var_indices),
       model@state_var_indices] <- model@solution$R
    QQ[model@state_var_indices, ] <- 
        model@solution$Q[, model@active_shocks]
    QQ[-c(model@state_var_indices), ]  <- 
        model@solution$S[, model@active_shocks]

    y_position = integer(0)
    if (!is.null(ref_var)) {
        if (length(ref_var) > 1) 
            stop('Only one reference variable can be specified')
        y_position = which(model@variables == ref_var)
        if (!length(y_position)) {
            stop('The reference variable is not a model variable. Choose another variable.')
        }
    }

    if (length(y_position) == 0) y_position = 1
    # saving position of variable
    model@var_position <- y_position


    # Computing covariance matrix
    if (!sim) {
        output <- covar(PP, QQ, 
                        model@cov_mat[model@active_shocks, model@active_shocks, drop = FALSE],
                        ngrid = ngrid, filter = filter,
                        lambda = lambda,
                        y_position=y_position,
                        n_leadlags=n_leadlags)
        # variance- covariance matrix
        cov_mat <- matrix(output$covmat,
                          length(model@variables),
                          length(model@variables))
        acov <- output$acov
        autcor_output <-  output$autcor_output
    } else {
        output <- covar_sim(PP, QQ, model@cov_mat[model@active_shocks, model@active_shocks, drop = FALSE],
                             filter = filter,
                             lambda = lambda, nrun,
                            y_position=y_position,
                            n_leadlags=n_leadlags)
        cov_mat <- output$covmat
        acov <- output$acov
        autcor_output <-  output$corr_variable_mat
    }   
    
    #autocorrelations
    colnames(acov) <- paste('t', seq(from = -1,
        to = -(n_leadlags - 1), by = -1), sep='')
    rownames(acov) <- model@variables
    
    #correlation with product
        
    colnames(autcor_output) <-
        paste(model@variables[y_position],'_[',
              seq(from = -(n_leadlags - 1),
                  to = (n_leadlags - 1), by = 1) ,
              ']', sep='')
    rownames(autcor_output) <- model@variables
    
    # Computing vardecomp
    PPs <- PP[, model@state_var_indices]
    err_var <- vardecomp(P = model@solution$P,
                         cov_mat = model@cov_mat[model@active_shocks, 
                                                     model@active_shocks, drop = FALSE],
                         Q = model@solution$Q[, model@active_shocks, drop = FALSE],
                         PP = PPs, QQ = QQ,
                         lambda = lambda,
                         ngrid = ngrid)
    rownames(err_var) <- model@variables
    colnames(err_var) <- model@shocks[model@active_shocks]

    
    # Passing results to gecon_model class
    # moments
    sdev <- sqrt(diag(cov_mat))
    var <- kronecker(sdev, t(sdev))
    model@sdev <- matrix(sdev, length(model@variables), 1)
    rownames(model@sdev) <- model@variables

    # correlation matrix
    model@corr_mat <- cov_mat / var
    model@corr_mat[is.nan(model@corr_mat)] <- 0
    rownames(model@corr_mat) <- model@variables
    colnames(model@corr_mat) <- model@variables

    
    
    # autocorrelation matrix
    model@autocorr_mat <- acov

    # correlation with specified variable
    model@corr_variable_mat <- autcor_output

    # variance decomposition
    model@var_dec <- err_var
    
    # checking control for true value
    model@corr_computed <- TRUE



    return(model)
}

# ###################################################################
# The get_moments function prints and returns moments of variables
# (absolute and relative to a chosen variable)
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   var_names - names of the variables, about which we want to get
#              the information 
#   relative_to - a logical, if TRUE then the function returns moments
#                 relative to one of the variables in accordance with
#                 relevant options chosen in compute_moments; if TRUE
#                 then only 'moments' and 'correalations' are active
#                 options)
#   moments - a logical, if TRUE then moments of variables: mean,
#             standard deviations and variance are returned;
#             if 'relative_to' is set to TRUE then moments ralative to
#             a chosen variable are returned
#   correlations - a logical, if TRUE a correlation matrix is returned;
#                  if relative_to is set to TRUE then corraletions of
#                  variables with lagged and leading values of a chosen
#                  variable are returned
#   autocorrelations - a logical, if TRUE then autocorrelations of variables
#                      are returned; if relative_to is set to TRUE then
#                      the option is inactive
#   var_dec - a logical, if TRUE then variance decomposition (of
#             shocks) is returned; if relative_to is set to TRUE then
#             the option is inactive
#   to_tex - a logical, enables printing output into tex file
#
# Output
#   Chosen results are printed on the console and a list of relevant
#   matrices is returned
# ###################################################################
get_moments <- function(model,
                        var_names = model@variables,
                        relative_to = FALSE,
                        moments = TRUE,
                        correlations = TRUE,
                        autocorrelations = TRUE,
                        var_dec = TRUE,
                        to_tex = FALSE)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    # Checking if correlations have been computed
    if (!model@corr_computed) {
        stop(paste('Compute correlations and moments first',
                   'using compute_moments method', sep = " "))
    }
    
    # Validation and getting indices of variables
    if(is(var_names, 'numeric')) {
        if (any(var_names > length(model@variables))) {
            stop(paste('Supplied index is greater than the number of variables.',
                       'Model has:', length(model@variables),
                       'variables.', sep = " "))
        } else var_ind <- var_names
    } else {
        # Checking if all supplied variables are model variables
        not_in_model <-
            (var_names[which(!(var_names %in% model@variables))])
        if (length(not_in_model))
            stop(paste('Following variables are not model variables:',
                       paste(not_in_model, collapse = ' '), sep = " "))
        var_ind <- which(model@variables %in% var_names)
    }
    
    if (relative_to) {
        
        # Validation and getting indices of variables
        if (length(model@var_position) == 0) {
            stop(paste('You have to compute correlation with',
                       'specified variable using compute_moments method',
                       sep = " "))
        }
        
        rel_ind <- model@var_position
        
        # Printing and returning results
        if (!any(moments, correlations)) {
            cat('\nNothing to be printed.\n')
        } else {    
            # creating data frame
            loglin_indic <- rep('Y    ', length(model@loglin_var[var_ind]))
            loglin_indic[which(!model@loglin_var[var_ind])] <- 'N    '
            mom <- round(cbind(model@steady[var_ind] / model@steady[rel_ind],
                         model@sdev[var_ind] / model@sdev[rel_ind],
                         model@sdev[var_ind] ^ 2 / model@sdev[rel_ind] ^ 2), digits = 4)
            mom_to_tex <- cbind(mom, loglin_indic)
            rownames(mom_to_tex) <- model@variables[var_ind]
            colnames(mom_to_tex) <- c(paste('Steady-state value relative to',
                                     model@variables[rel_ind]),
                               paste('Std. dev. relative to',
                                     model@variables[rel_ind]),
                               paste('Variance relative to',
                                     model@variables[rel_ind]),
                                    'Loglinear')
            mom <- as.data.frame(mom_to_tex)
            
            lead_lag <- model@corr_variable_mat[var_ind, ]
            
            # number of leads and lags
            n_ll <- (dim(lead_lag)[2] - 1) / 2
            
            if (moments & correlations) {
                results <- list(Moments_relative_to_moments_of_the_reference_variable = mom,
                                Correlations_with_the_reference_variable = lead_lag)
                if (to_tex) {
                    results_tex <- list(Moments_relative_to_moments_of_the_reference_variable = mom_to_tex,
                        Correlations_with_the_reference_variable = lead_lag)
                    nam <- names(results_tex)
                    tex_nam <- vector("list", length(nam))
                    nam_m <- which(nam == 'Moments_relative_to_moments_of_the_reference_variable') 
                    if (length(nam_m)) {
                        tex_nam[[nam_m]] <-  c(paste('Steady-state value relative to',
                                                   string_to_tex(model@variables[rel_ind])),
                                                paste('Std. dev. relative to',
                                                    string_to_tex(model@variables[rel_ind])),
                                                paste('Variance relative to',
                                                        string_to_tex(model@variables[rel_ind])),
                                                'Loglinear')
                    }  
                    nam_c <- which(nam == 'Correlations_with_the_reference_variable') 
                    if (length(nam_c)) {
                        time_ind <- paste('t',c(c(-n_ll:-1),''), sep='')
                        time_ind <- c(time_ind, paste('t+', c(1:n_ll), sep=''))
                        tex_nam[[nam_c]] <- paste('$', string_to_tex(model@variables[rel_ind], math_mode =FALSE),
                                                   '_{', time_ind, '}$', sep ='')
                    } 

                    to_tex(results_tex, 
                       paste(model@model_info[2], '.results.tex', sep=''),
                       section = "Statistics of the model",
                       colnam = tex_nam,
                       rownam = string_to_tex(model@variables[var_ind]))
                } else {
                    cat(paste('\nMoments of variables relative to',
                              model@variables[rel_ind], ': \n\n'))
                    print(mom)
                    cat(paste('\nCorrelations of variables with lead and lagged',
                              model@variables[rel_ind], ': \n\n'))
                    print(round(lead_lag, digits = 4))
                    cat('\n')
                }
                return(invisible(results))
            } else if (moments & !correlations) {
                if (to_tex) {
                    tex_nam <- vector("list", 1)
                    tex_nam <-  c(paste('Steady-state value relative to',
                                                    string_to_tex(model@variables[rel_ind])),
                                                paste('Std. dev. relative to',
                                                    string_to_tex(model@variables[rel_ind])),
                                                paste('Variance relative to',
                                                        string_to_tex(model@variables[rel_ind])),
                                                'Loglinear')
                    to_tex(mom_to_tex, 
                           paste(model@model_info[2], '.results.tex', sep=''),
                           section = paste("Moments relative to moments of the reference value", string_to_tex(model@variables[rel_ind])),
                           colnam = tex_nam,
                           rownam = string_to_tex(model@variables[var_ind]))
                } else {
                    cat(paste('\nMoments of variables relative to',
                              model@variables[rel_ind], ': \n\n'))
                    print(mom)
                    cat('\n')
                }
                return(invisible(mom))
            } else {
                if (to_tex) {
                    time_ind <- paste('t',c(c(-n_ll:-1),''), sep='')
                    time_ind <- c(time_ind, paste('t+', c(1:n_ll), sep=''))
                    tex_nam <- vector("list", 1)
                    tex_nam <- paste('$', string_to_tex(model@variables[rel_ind], math_mode =FALSE),
                                               '_{', time_ind, '}$', sep ='')
                    to_tex(lead_lag, 
                       paste(model@model_info[2], '.results.tex', sep=''),
                       section = paste("Correlations with lead and lagged", model@variables[rel_ind]),
                       colnam = tex_nam,
                       rownam = string_to_tex(model@variables[var_ind]))
                } else {
                    cat(paste('\nCorrelations of variables with lead and lagged',
                              model@variables[rel_ind], ': \n\n'))
                    print(round(lead_lag, digits = 4))
                    cat('\n')                
                }
                return(invisible(lead_lag))
            }
        }
    } else {
        
        if (!any(moments, correlations, autocorrelations, var_dec)) {
            cat('\nNothing to be printed.\n')
        } else {
            # creating data frame
            loglin_indic <- rep('Y    ', length(model@loglin_var[var_ind]))
            loglin_indic[which(!model@loglin_var[var_ind])] <- 'N    '
            mom <- round(cbind(model@steady[var_ind],  
                                   model@sdev[var_ind], 
                                   model@sdev[var_ind] ^ 2), digits = 4)

            mom_to_tex <- cbind(mom, loglin_indic)
            rownames(mom_to_tex) <- model@variables[var_ind]
            colnames(mom_to_tex) <- c('Steady-state value', 'Std. dev.', ' Variance', 'Loglinear')
            mom <- as.data.frame(mom_to_tex)
             

            correl <- model@corr_mat[, var_ind]
            autocorrel <- as.matrix(model@autocorr_mat[var_ind, ])
            var_decomp <- as.matrix(model@var_dec[var_ind, ])
            colnames(var_decomp) <- model@shocks[model@active_shocks]
            
            # Returning results
            tables <- list(moments = mom,
                           correlation_matrix = correl,
                           autocorrelations = autocorrel,
                           variance_decomposition = var_decomp)
               
            tables_tex <- list(Moments = mom_to_tex,
                               Correlation_matrix = correl,
                               Autocorrelations = autocorrel,
                               Variance_decomposition = var_decomp)
            tables_log <- c(moments, correlations,
                            autocorrelations, var_dec)
            tables_ind <- which(tables_log)
            chosen_tables <- tables[tables_ind]
            chosen_tables_tex <- tables_tex[tables_ind]
            
            # Returning tex tables
            if (to_tex) {
                nam <- names(chosen_tables_tex)
                tex_nam <- vector("list", length(nam))
                nam_m <- which(nam == 'Moments') 
                if (length(nam_m)) {
                    tex_nam[[nam_m]] <- c('Steady-state value', 'Std. dev.', 'Variance', 'Loglinear')
                }  

                nam_c <- which(nam == 'Correlation_matrix') 
                if (length(nam_c)) {
                    tex_nam[[nam_c]] <- string_to_tex(colnames(chosen_tables_tex[[nam_c]]))
                } 

                nam_a <- which(nam == 'Autocorrelations') 
                if (length(nam_a)) {
                    tex_nam[[nam_a]] <- string_to_tex(colnames(chosen_tables_tex[[nam_a]]))
                }

                nam_v <- which(nam == 'Variance_decomposition') 
                if (length(nam_v)) {
                    tex_nam[[nam_v]] <- string_to_tex(model@shocks)
                }
                
                to_tex(chosen_tables_tex, 
                       paste(model@model_info[2], '.results.tex', sep=''),
                       section = "Statistics of the model",
                       colnam = tex_nam)
            } else {
                if (moments) {
                    cat('\nMoments of variables:\n\n')
                    print(mom)
                }
                if (correlations) {
                    cat('\nCorrelations of variables:\n\n')
                    print(round(correl, digits = 4))
                }
                if (autocorrelations) {
                    cat('\nAutocorrelations of variables:\n\n')
                    print(round(autocorrel, digits = 4))
                }
                if (var_dec) {
                    cat('\nDecomposition of variance:\n\n')
                    print(round(var_decomp, digits = 4))
                }
                cat('\n')    
            }

            
            return(invisible(chosen_tables))
        }
    }
}

# ###################################################################
# The show method controls how object of class gecon_model is printed
# ###################################################################
# Input
#   object - an object of class gecon_model
# Output
#   Information about object
# ###################################################################
setMethod("show", signature(object = "gecon_model"),
function(object)
{
    model2show <- object
    if (!is.gecon_model(model2show)) {
        stop("The model argument has to be of class gecon_model")
    }
    model_name <- model2show@model_info[1]
    no_var <- length(model2show@variables)
    no_shocks <- length(model2show@shocks)
    no_par <- length(model2show@parameters)
    no_par_mod <- sum(model2show@parameters_free_mod_flag)
    no_calibr_par <- length(model2show@parameters_calibr)
    no_state <- length(model2show@state_var_indices)

    cat("\n")
    if (!model2show@is_dynamic) {
        cat('\'', model_name, '\' is a static model.\n\n', sep = '')
    } else if (model2show@is_stochastic) {
        cat('\'', model_name, '\' is a dynamic, stochastic model.\n\n', sep = '')
    } else {
        cat('\'', model_name, '\' is a dynamic, deterministic model.\n\n', sep = '')
    }

    cat('Model has been read from file: \'', to_gcn(model2show@model_info[1]), '\'',
        sep = '')
    cat('\nDate of creation:', model2show@model_info[3], '\n\n')

    if (model2show@is_dynamic)
        cat("Steady-state values of variables have ")
    else
        cat("Equilibrium has ")
        
    if (model2show@ss_solved) {
        cat("been FOUND.\n")
    } else {
        cat("NOT been FOUND.\n")
    }

    if (length(model2show@calibr_equations)) {
        if (model2show@is_calibrated && model2show@ss_solved) {
            if (model2show@is_dynamic) {
                cat("Calibrating equations have been taken \n")
                 cat("into account when solving for the steady state \n")
            } else {
                cat("Calibrating equations have been taken \n")
                 cat("into account when solving for the equilibrium\n")
            }
        } else if (!model2show@is_calibrated && model2show@ss_solved) {
            if (model2show@is_dynamic) {
                cat("Calibrating equations have not been taken \n")
                 cat("into account when solving for the steady state \n")
            } else {
                cat("Calibrating equations have not been taken \n")
                 cat("into account when solving for the equilibrium \n")        
            }
        }
    }
    
    if (model2show@re_solved) {
        if (model2show@loglin) {
            cat("The log-linearised version of the model has been SOLVED.\n")
        } else {
            cat("The linearised version of the model has been SOLVED.\n")
        }
    } else if (model2show@is_dynamic) {
        cat("The model in its (log-)linearised form has NOT been SOLVED.\n")
    } else {
        cat("\n")
    }

    if (model2show@corr_computed) {
        cat("Correlations have been COMPUTED.\n\n")
    } else if(model2show@is_dynamic) {
            cat("Correlations have NOT been COMPUTED.\n\n")
    }

    cat("Number of variables:", no_var, "\n")
    cat("Number of parameters:", no_par, "\n")
    cat("Calibrated parameters:", no_calibr_par, "\n")
    if (model2show@is_dynamic) {
        cat("Number of stochastic shocks:", no_shocks, "\n")
    }
    if (model2show@is_dynamic) {
        cat("Number of state variables: ")
        if (model2show@re_solved) {
            cat(no_state, "\n\n")
        } else {
            cat("?\n\n")
        }
    }
})


# ###################################################################
# The print method prints basic information and diagnostics of 
# gecon_model objects
# ###################################################################
# Input
#   x - an object of class gecon_model
# Output
#   Information about the model:
#       information about steady state or equilibrium and perturbation 
#       solution status, the number of variables of different sorts, 
#       parameters and shocks
# ###################################################################
setMethod("print", signature(x = "gecon_model"),
function (x, ...)
{
    model2print <- x
    if (!is.gecon_model(model2print)) {
        stop("The model argument has to be of class gecon_model")
    }

    # model_name <- deparse(substitute(x))
    model_name <- model2print@model_info[1]
    no_var <- length(model2print@variables)
    no_shocks <- length(model2print@shocks)
    no_par <- length(model2print@parameters)
    no_par_mod <- sum(model2print@parameters_free_mod_flag)
    no_calibr_par <- length(model2print@parameters_calibr)
    no_state <- length(model2print@state_var_indices)
    variables <- model2print@variables
    calibr_parameters <- model2print@parameters_calibr
    shocks <- model2print@shocks
    parameters <- model2print@parameters
    parameters_modified <-
        model2print@parameters_free[
            as.logical(model2print@parameters_free_mod_flag)]
    state_variables <-
        model2print@variables[model2print@state_var_indices]

    cat("\n")
    if (!model2print@is_dynamic) {
        cat('\'', model_name, '\' is a static model.\n\n', sep = '')
    } else if (model2print@is_stochastic) {
        cat('\'', model_name, '\' is a dynamic, stochastic model.\n\n', sep = '')
    } else {
        cat('\'', model_name, '\' is a dynamic, deterministic model.\n\n', sep = '')
    }

    cat('Model has been read from file:', to_gcn(model2print@model_info[2]))
    cat('\nDate of creation:', model2print@model_info[3], '\n\n')
    
    if (model2print@is_dynamic) {
        cat("Steady-state values of variables have ")
    } else {
        cat("Equilibrium has ")        
    }
    
    if (model2print@ss_solved) {
        cat("been FOUND.\n")
    } else {
        cat("NOT been FOUND.\n")
    }
    
    if (length(model2print@calibr_equations)) {
        if (model2print@is_calibrated && model2print@ss_solved) {
            if (model2print@is_dynamic) {
                cat("Calibrating equations have been taken \n")
                 cat("into account when solving for the steady state. \n")
            } else {
                cat("Calibrating equations have been taken \n")
                 cat("into account when solving for the equilibrium. \n")
            }
        } else if (!model2print@is_calibrated && model2print@ss_solved) {
            if (model2print@is_dynamic) {
                cat("Calibrating equations have not been taken \n")
                 cat("into account when solving for the steady state. \n")
            } else {
                cat("Calibrating equations have not been taken \n")
                 cat("into account when solving for the equilibrium. \n")        
            }
        }
    }
    
    if (model2print@re_solved) {
        if (model2print@loglin) {
            cat("The log-linearised version of the model has been SOLVED.\n")
        } else {
            cat("The linearised version of the model has been SOLVED.\n")
        }
    } else if (model2print@is_dynamic) {
            cat("The model in its (log-)linearised form has NOT been SOLVED.\n")
    } else {
        cat("\n")
    }

    if (model2print@corr_computed) {
        cat("Correlations have been COMPUTED.\n\n")
    } else {
        if (model2print@is_dynamic) {
            cat("Correlations have NOT been COMPUTED.\n\n")
        }
    }

    cat("----------------------------------------------------------", "\n\n")

    cat("Number of variables:", no_var, "\n")
    cat("\n")
    print(as.data.frame(variables))
    cat("\n----------------------------------------------------------", "\n\n")

    cat("Number of parameters:", no_par, "\n")
    cat("\n")
    if (no_par) {
        print(as.data.frame(parameters))
    }
    cat("\n----------------------------------------------------------", "\n\n")

    if(no_par_mod) {
        cat("Following parameters have been modified with respect to .gcn file \n",
            "or their values have been set in R interface:", no_par_mod, "\n")
        cat("\n")
        print(as.data.frame(parameters_modified))
        cat("\n----------------------------------------------------------", "\n\n")
    }

    cat("Calibrated parameters:", no_calibr_par, "\n")
    if (no_calibr_par != 0) {
        cat("\n")
        print(as.data.frame(calibr_parameters))
    }
    cat("\n----------------------------------------------------------", "\n\n")

    if (no_shocks != 0) {
        cat("Number of stochastic shocks:", no_shocks, "\n")
        cat("\n")
        print(as.data.frame(shocks))
    } else {
        cat("The model is deterministic. It does not have a stochastic part.\n")
    }
    cat("\n----------------------------------------------------------", "\n\n")

    if (model2print@re_solved) {
        cat("Number of state variables:", no_state, "\n")
        cat("\n")
        print(as.data.frame(state_variables))
        cat("\n----------------------------------------------------------", "\n\n")
    } else if (model2print@is_dynamic) {
            cat("In order to identify state variables you need to solve\n")
            cat("the model in its (log-)linear form first.\n\n")
    }
})


# ###################################################################
# The summary method prints and returns the summary of model solution
# ###################################################################
# Input
#   object - an object of class gecon_model
# Output
#   Summary of model solution
# ###################################################################
setMethod("summary", signature(object = "gecon_model"),
function (object)
{
    if (!is.gecon_model(object)) {
        stop("The model argument has to be of class gecon_model");
    }
    if (!object@ss_solved) {
        cat('Model has NOT been SOLVED.\n')
        return(invisible())
    }
    if (object@is_dynamic) {
        cat('\nSteady state:\n')
    } else {
        cat('\nEquilibrium:\n')    
    }
    cat_ss <- as.data.frame(object@steady)
    colnames(cat_ss) <- c("")
    rownames(cat_ss) <- object@variables
    print(round(cat_ss, digits = 6))

    if (length(object@params)) {
        cat('\n----------------------------------------------------------\n')
        cat('\nParameter values:\n')
        cat_par <- as.data.frame(object@params)
        colnames(cat_par) <- c("")
        rownames(cat_par) <- object@parameters
        print(cat_par)
    }

    if (object@is_dynamic) {
        cat('\n----------------------------------------------------------\n')
    } else {
        return(invisible())
    }

    if (!object@re_solved && object@is_dynamic) {
        cat('\nThe model in its (log-)linearised form has NOT been SOLVED.\n\n')
        return(invisible())
    } else if (!object@is_dynamic) {
        return(invisible())    
    }
    
    cat('\nLinearisation:\n')
    cat('x_t = P x_{t-1} + Q epsilon_t\n')
    cat('y_t = R x_{t-1} + S epsilon_t\n')
    if (length(object@solver_exit_info)) {
        cat('\n', object@solver_exit_info, '\n', sep = '')
    }
    

    if (!length(object@state_var_indices)) {
        cat('Model does NOT have state variables.',
            'All the solution matrices are empty.')
        return(invisible())
    }
        
    cat('\nP:\n')
    print(round(object@solution$P, digits = 6))
    
    if (object@is_stochastic) {
        cat('\nQ:\n')
        print(round(object@solution$Q, digits = 6))
    }
    
    cat('\nR:\n')
    print(round(object@solution$R, digits = 6))
    
    if (object@is_stochastic) {
        cat('\nS:\n')
        print(round(object@solution$S, digits = 6))
    }
    
    if (object@is_stochastic) {
        cat('\n----------------------------------------------------------\n')
    } else {
        return(invisible())
    }

    if (!object@corr_computed) {
        if (object@is_dynamic)
            cat('\nCorrelations have NOT been COMPUTED.\n\n')
        return(invisible())
    }

    cat('\nShock covariance matrix:\n\n')
    rownames(object@cov_mat) <- object@shocks
    colnames(object@cov_mat) <- object@shocks
    print(object@cov_mat)

    cat("\nSummary statistics:\n\n")
    loglin_indic <- rep('Y    ', length(object@loglin_var))
    loglin_indic[which(!object@loglin_var)] <- 'N    '
    cat_mom <- round(cbind(object@steady,  object@sdev, object@sdev ^ 2), digits = 4)
    cat_mom <- as.data.frame(cbind(cat_mom, loglin_indic))
    rownames(cat_mom) <- object@variables
    colnames(cat_mom) <- c('Steady-state value', 'Std. dev.', ' Variance', 'Loglinear')
    print(cat_mom)

    nzind <- which(object@sdev != 0)

    cat("\nCorrelation matrix:\n\n")
    cat_corrm <- object@corr_mat[nzind, nzind]
    cat_corrm <- as.data.frame(cat_corrm)
    N <- length(nzind)
    for (j in (1:N)) {
        cat_corrm[, j] <- as.character(round(cat_corrm[, j], digits = 6))
        if (j < N) {
            for (i in ((j+1):N)) {
                cat_corrm[i, j] <- ""
            }
        }
    }
    print(cat_corrm)

    cat("\nAutocorrelations:\n\n")
    print(round(object@autocorr_mat[nzind, ], digits = 6))

    cat("\nCross correlations:\n\n")
    print(round(object@corr_variable_mat[nzind, ], digits = 6))

    cat("\nVariance decomposition:\n\n")
    cat_vd <- object@var_dec[nzind, ]
    dim(cat_vd) <- c(length(nzind), dim(object@var_dec)[2])
    colnames(cat_vd) <- colnames(object@var_dec)
    rownames(cat_vd) <- rownames(object@var_dec)[nzind]
    print(round(cat_vd, digits = 6))
    cat("\n")
})

# ###################################################################
# Functions for removing rows and columns with zeros as only elements
# ###################################################################

# ###################################################################
# The nonempty_col_indices function finds indices of empty columns
# ###################################################################
# Input 
#   inp_mat - a matrix to be reduced
#   tol - numerical inrelevance
# Output
#   vector of empty columns indices
# ###################################################################
nonempty_col_indices <- function(inp_mat, tol = 1e-8) 
{
    a <- dim(inp_mat)
    indices <- which(colSums(abs(inp_mat)) > tol)
    return(indices)
}

# ###################################################################
# The nonempty_row_indices function finds indices of empty rows
# ###################################################################
# Input 
#   inp_mat - a matrix to be reduced
#   tol - numerical inrelevance
# Output
#   vector of empty rows indices
# ###################################################################
nonempty_row_indices <- function(inp_mat, tol = 1e-8) 
{
    return(nonempty_col_indices(t(inp_mat), tol))
}

# ###################################################################
# The to_remove function finds indices of empty rows
# ###################################################################
# Input 
#   mat - a matrix to be reduced
#   tol - numerical inrelevance
# Output
#   vector of empty rows indices
# ###################################################################
to_remove <- function(mat, tol = 1e-6)
{
    mat <- as.matrix(mat)
    size <- dim(mat)[2]
    size_c <- dim(mat)[1]
    mat <- t(mat)
    to_be_removed <- which(abs(mat[max.col(abs(t(mat))) + 
                        size * c(0: (size_c - 1))]) < tol)
    return(to_be_removed)
}

# ###################################################################
# The remo function removes specified indices from matrix
# ###################################################################
# Input 
#   mat - a matrix to be reduced
#   to_rem - indices of rows to be removed
# Output
#   matrix with removed rows
# ###################################################################
remo <- function(mat, to_rem)
{
    mat <- Matrix(mat[-c(to_rem),])
}

# ###################################################################
# The row_reduction function removes empty rows from 3 first 
#   elements of a list of matrices
# ###################################################################
# Input 
#   can - a list of matrices
# Output
#   list of matrices with removed rows
# ###################################################################
row_reduction <- function(can)
{
    can2 <- list(can[[1]], can[[2]], can[[3]])
    rem <- Reduce(intersect, lapply(can2, to_remove))
    if (length(rem) > 0) {
        can <- lapply(can, remo, rem)
    }
    return(can)
}

# ###################################################################
# The get_var_names_by_index function returns 
# the names of variables indexed by the indices specified 
# in the index_names argument
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   index_names - a character vector with the names of indices
# Output
#   A character vector with relevant variable names
# ###################################################################
get_var_names_by_index <- function(model, index_names)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    if (!is.character(index_names))
        stop("The index_names argument has to be a vector of characters")

    if (!length(model@index_sets))
        stop("The model has NOT got any index sets specified")
        
    # Check if index sets are correctly specified
    ind_sets <- unlist(model@index_sets)
    
    not_indices <- which(!(index_names %in% ind_sets))

    if (length(not_indices))
        stop(paste("The following index names have been misspelled:\n",
                   paste(index_names[not_indices], sep = ", ")))
    
    # Finding relevant variables
    index_names_alternative <- paste('(', paste(index_names,  collapse ='|'), ')', collapse = '', sep = '')
    assoc_names  <- grep(pattern = paste('__', index_names_alternative,
                                         '(__)*(?!(_?[a-zA-Z0-9]))', sep=''), 
                                x = model@variables,
                                perl = T)
    
    if (!length(assoc_names))
        stop("The model does NOT contain any variables indexed by the specified indices")
        
    return(model@variables[assoc_names])
}

# ###################################################################
# The get_par_names_by_index function returns 
# the names of parameters indexed by the indices specified 
# in the index_names argument
# ###################################################################
# Input
#   model - an object of the gecon_model class
#   index_names - a character vector with the names of indices
# Output
#   A character vector with relevant parameter names
# ###################################################################
get_par_names_by_index <- function(model, index_names)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    if (!is.character(index_names))
        stop("The index_names argument has to be a vector of characters")

    if (!length(model@index_sets))
        stop("The model has NOT got any index sets specified")
        
    # Check if index sets are correctly specified
    ind_sets <- unlist(model@index_sets)
    
    not_indices <- which(!(index_names %in% ind_sets))
    
    if (length(not_indices))
        stop(paste("The following index names have been misspelled:\n",
                   paste(index_names[not_indices], sep = ", ")))
    
    # Finding relevant parameters
    index_names_alternative <- paste('(', paste(index_names,  collapse ='|'), ')', collapse = '', sep = '')
    assoc_names  <- grep(pattern = paste('__', index_names_alternative,
                                         '(__)*(?!(_?[a-zA-Z0-9]))', sep=''), 
                                x = model@parameters,
                                perl = T)
    
    if (!length(assoc_names))
        stop("The model does NOT contain any parameters indexed by the specified indices")
        
    return(model@parameters[assoc_names])
}

# ###################################################################
# The get_shock_names_by_index function returns 
# the names of shocks indexed by the indices specified 
# in the index_names argument
# ###################################################################
# Input
#   model - an object of the gecon_model class.
#   index_names - a character vector with the names of indices.
# Output
#   A character vector with relevant shock names.
# ###################################################################
get_shock_names_by_index <- function(model, index_names)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    if (!is.character(index_names))
        stop("The index_names argument has to be a vector of characters")

    if (!length(model@index_sets))
        stop("The model has NOT got any index sets specified")
        
    # Check if index sets are correctly specified
    ind_sets <- unlist(model@index_sets)
    
    not_indices <- which(!(index_names %in% ind_sets))
    
    if (length(not_indices))
        stop(paste("The following index names have been misspelled:\n",
                   paste(index_names[not_indices], sep = ", ")))
    
    # Finding relevant parameters
    index_names_alternative <- paste('(', paste(index_names,  collapse ='|'), ')', collapse = '', sep = '')
    assoc_names  <- grep(pattern = paste('__', index_names_alternative,
                                         '(__)*(?!(_?[a-zA-Z0-9]))', sep=''), 
                                x = model@shocks,
                                perl = T)
    
    if (!length(assoc_names))
        stop("The model does NOT contain any shocks indexed by the specified indices")
    
    return(model@shocks[assoc_names])
}

# ###################################################################
# The get_index_sets function allows to retrieve a list with
# the index sets from an object of gecon_model class
# ###################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#   A list with index sets specified for the model.
# ###################################################################
get_index_sets <- function(model)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
   
    if (!length(model@index_sets))
        stop("The model has NOT got any index sets specified")
        
    return(model@index_sets)
}

# ###################################################################
# The re_solved function returns a logical value indicating
# if the perturbation has been solved.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#   A logical value. If TRUE, the model has been solved.
# ###################################################################
re_solved <- function(model)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
        
    return(model@re_solved)
}

# ###################################################################
# The ss_solved function returns a logical value indicating
# if the steady state (equilibrium) for the model has been found.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#   A logical value. If TRUE, the steady state (equilibrium)
#   for the model has been found.
# ###################################################################
ss_solved <- function(model)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
        
    return(model@ss_solved)
}

# ###################################################################
# The get_cov_mat function returns a variance-covariance 
# matrix of model shocks.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#   A a variance-covariance matrix of model shocks.
# ###################################################################
get_cov_mat <- function(model)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")
    
    cov_mat <- model@cov_mat
    rownames(cov_mat) <- model@shocks
    colnames(cov_mat) <- model@shocks

    return(cov_mat)
}

# ###################################################################
# The get_model_info function returns a character vector
# with information about the model.
# ###################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#    a character vector of length 3, containing information 
#    about the model: the input file name, the input file path, 
#    and the date of creation.
# ###################################################################
get_model_info <- function(model)
{
    if (!is.gecon_model(model))
        stop("The model argument has to be of class gecon_model")

    return(model@model_info)
}


