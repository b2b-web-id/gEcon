# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Class for storing information about model variables
# ###################################################################

# ###################################################################
# Class definition
# ###################################################################

setClass(
    Class = "gecon_var_info",
    representation = representation(model_info = "character",
                                    model_variable_name = "character",
                                    var_names = "character",
                                    is_stochastic = "logical",
                                    is_dynamic = "logical",
                                    ss_solved = "logical",
                                    re_solved = "logical",
                                    corr_computed = "logical",
                                    ss_val = "matrix",
                                    state = "logical",
                                    state_var_impact = "matrix",
                                    shock_impact = "matrix",
                                    std_dev_val = "numeric",
                                    loglin_flag = "logical",
                                    cr = "matrix",
                                    incid_mat = "Matrix")
)

# ###################################################################
# The function gecon_var_info is 
# a gecon_var_info class object constructor
# ###################################################################
# Input 
#   model_info - a character vector of length 3, containing 
#                information about the model: the input file name,
#                the input file path and the date of creation.
#   model_variable_name - a character denoting the name of variable 
#                         storing model for which the simulations 
#                         have been created.
#   var_names - a character vector of the variables names.
#   is_stochastic - logical. If TRUE, the model has stochastic shocks.
#   is_dynamic - logical. If TRUE, the model has any lead 
#                or lagged variables.
#   ss_solved - logical. If TRUE, steady state (equilibrium
#               for static models) has been found.
#   re_solved - logical. It is set to TRUE if the model 
#               has been solved. The default value is FALSE.
#   corr_computed - logical. If TRUE, indicates that 
#                   the correlations and other statistics 
#                   of variables have been computed. 
#                   The default value is FALSE.
#   ss_val - a vector of the steady state values of variables 
#            for dynamic models or equilibrium for static models. 
#            If the steady state has not been computed when the object
#            was generated, this slot contains initial values 
#            of variables.
#   state - logical. IF TRUE, indicates that the corresponding 
#           variable is a state variable.
#   state_var_impact - a matrix with the impact of state variables values
#           in the previous period on the variables values.
#   shock_impact - a matrix with the impact of exogenous variables 
#             values on the variables values.
#   std_dev_val - a numeric vector of standard deviations of variables.
#   loglin_flag - a logical vector with the length equal 
#                 to the number of the variables. 
#                 The TRUE entries denote that a corresponding variable 
#                 has been loglinearized before solving the model.
#   cr - a matrix with correlations of the variables 
#        with all the model variables.
#   incid_mat - a Matrix object with the mapping of variables 
#               into equations and calibrating equations.
# Output
#   An object of the gecon_var_info class.
# ###################################################################
gecon_var_info <- function(model_info, model_variable_name, 
                           indices, var_names, is_stochastic,
                           is_dynamic, ss_solved, re_solved,
                           corr_computed, ss_val, 
                           state, state_var_impact, shock_impact,
                           std_dev_val, loglin_flag,
                           cr, incid_mat) 
{
    var_info_object <- new('gecon_var_info')
     
    if (!is.character(model_info)) {
        stop('model_info should be of vector class')
    } else var_info_object@model_info <- model_info   

    if (!is.character(model_variable_name)) {
        stop('model_variable_name should be of vector class')
    } else var_info_object@model_variable_name <- model_variable_name   

    if (!is.character(var_names)) {
        stop('var_names should be of character class')
    } else var_info_object@var_names <- var_names  

    if (!is.logical(is_stochastic)) {
        stop('is_stochastic should be of logical class')
    } else var_info_object@is_stochastic <- is_stochastic      

    if (!is.logical(is_dynamic)) {
        stop('is_dynamic should be of logical class')
    } else var_info_object@is_dynamic <- is_dynamic
    
    if (!is.logical(ss_solved)) {
        stop('ss_solved should be of logical class')
    } else var_info_object@ss_solved <- ss_solved   
    
    if (!is.logical(re_solved)) {
        stop('re_solved should be of logical class')
    } else var_info_object@re_solved <- re_solved  
        
    if (!is.logical(corr_computed)) {
        stop('corr_computed should be of logical class')
    } else var_info_object@corr_computed <- corr_computed  

    if (!is.matrix(ss_val)) {
        stop('ss_val should be of matrix class')
    } else var_info_object@ss_val <- ss_val     
    
    if (!is.logical(state)) {
        stop('state should be of logical class')
    } else var_info_object@state <- state   
    
    if (!is.matrix(state_var_impact)) {
        stop('state_var_impact should be of matrix class')
    } else var_info_object@state_var_impact <- state_var_impact   
    
    if (!is.matrix(shock_impact)) {
        stop('shock_impact should be of matrix class')
    } else var_info_object@shock_impact <- shock_impact  

    if (!is.numeric(std_dev_val)) {
        stop('std_dev_val should be of numeric class')
    } else var_info_object@std_dev_val <- std_dev_val    

    if (!is.logical(loglin_flag)) {
        stop('loglin_flag should be of logical class')
    } else var_info_object@loglin_flag <- loglin_flag
    
    if (!is.matrix(cr)) {
        stop('cr should be of matrix class')
    } else var_info_object@cr <- cr 
    
    if (!inherits(incid_mat, 'Matrix')) {
        stop('incid_mat should be of Matrix class')
    } else var_info_object@incid_mat <- incid_mat 

    return(var_info_object)
}



# ###################################################################
# The show method controls how the gecon_var_info
# object is printed
# ###################################################################
# Input 
#   object - object of class gecon_var_info
# Output
#   Information about variables incidence and already computed
#   statistics for the variables
# ################################################################### 
setMethod(
    "show", signature(object = "gecon_var_info"),
    function(object) {
    
        # Print incidence matrix
        ch_incid_mat <- matrix(as.character(object@incid_mat), 
                               nrow = nrow(object@incid_mat),
                               ncol = ncol(object@incid_mat))
        ch_incid_mat[which(object@incid_mat == 0)] <-  '      .'
        ch_incid_mat[which(object@incid_mat == 1)] <-  '    t-1'
        ch_incid_mat[which(object@incid_mat == 2)] <-  '      t'
        ch_incid_mat[which(object@incid_mat == 3)] <-  ' t-1, t'
        ch_incid_mat[which(object@incid_mat == 4)] <-  '    t+1'
        ch_incid_mat[which(object@incid_mat == 5)] <-  't-1, t+1'
        ch_incid_mat[which(object@incid_mat == 6)] <-  ' t, t+1'
        ch_incid_mat[which(object@incid_mat == 7)] <-  't-1, t, t+1'
        ch_incid_mat[which(object@incid_mat == 8)] <-  '     ss'
        ch_incid_mat[which(object@incid_mat == 9)] <-  'ss, t-1'
        ch_incid_mat[which(object@incid_mat == 10)] <- '  ss, t'
        ch_incid_mat[which(object@incid_mat == 11)] <- 'ss, t-1, t'
        ch_incid_mat[which(object@incid_mat == 12)] <- 'ss, t+1'
        ch_incid_mat[which(object@incid_mat == 13)] <- 'ss, t-1, t+1'
        ch_incid_mat[which(object@incid_mat == 14)] <- 'ss, t, t+1'
        ch_incid_mat[which(object@incid_mat == 15)] <- 'ss, t-1, t, t+1'
        colnames(ch_incid_mat) <- colnames(object@incid_mat)
        rownames(ch_incid_mat) <- rownames(object@incid_mat)
        cat('\nIncidence info:\n\n')
        print(as.data.frame(ch_incid_mat))
        
        # Print info about value
        cat("\n----------------------------------------------------------", "\n")
        if (object@ss_solved) {
            if (object@is_dynamic) {
                cat('\nSteady-state values:\n\n')
            } else {
                cat('\nEquilibrium values:\n\n')        
            }
        } else {
            if (object@is_dynamic) {
                cat('Initial values for SS computation:\n\n')
            } else {
                cat('Initial values for equilibrium computation:\n\n')        
            }
        }
        print(round(object@ss_val, digits = 4))
        
        # Print info about loglinearized and state variables        
        if (object@re_solved) {        
            ch_state <- matrix('', length(object@var_names), 2)        
            colnames(ch_state) <- c('Is a state variable?', 'Is loglinearized?')
            rownames(ch_state) <- object@var_names
            ch_state[which(object@state), 1] <- 'Y          '
            ch_state[which(object@loglin_flag), 2] <- 'Y          '
            cat("\n----------------------------------------------------------", "\n")
            cat('\nVariable info:\n\n')
            print(as.data.frame(ch_state))
        }
        
        # Print info about the recursive laws of motion
        if (object@re_solved) {
            cat("\n----------------------------------------------------------", "\n")        
            cat('\nRecursive laws of motion for the variables\n')
                    
            cat('\nState variables impact:\n\n')
            print(round(object@state_var_impact, digits = 4))
            
            if (object@is_stochastic) {
                cat('\nShocks impact:\n\n')
                print(round(object@shock_impact, digits = 4))
            }
        }        
        
        # Print info about moments and correlations
        if (object@corr_computed) {
            loglin_indic <- rep('N    ', length(object@var_names))
            loglin_indic[object@loglin_flag] <- 'Y    '
            mom <- round(cbind(object@ss_val, 
                               object@std_dev_val,
                               object@std_dev_val ^ 2), 4)

            rownames(mom) <- object@var_names
            mom <- as.data.frame(cbind(mom, loglin_indic)) 
            colnames(mom) <- c('Steady-state value', 
                               ' Std. dev.', 
                               ' Variance', 'Loglin') 
            cat("\n----------------------------------------------------------", "\n\n")
            cat('Moments:\n\n')        
            print(mom)
        
            # Correlations
            cat('\nCorrelations:\n\n')
            print(round(object@cr, digits = 4))
            cat('\n')
        }
    }
)

# ###################################################################
# The print method prints information about gecon_var_info object.
# ###################################################################
# Input 
#   x - an object of the gecon_var_info class. 
# Output
#   Information about the gecon_var_info object.
# ###################################################################
setMethod(
    "print",
    signature(x = "gecon_var_info"), 
    function(x) 
    {       
        
        cat("\n")
        cat('The information about variables has been generated\n',
            'based on model stored in variable: ', 
            x@model_variable_name, '.\n', sep = '')
        cat("\n")
        cat('The model has been created based on:', 
            x@model_info[1], 'file. \n')
        cat("\n")
        cat('The information has been created for', 
            length(x@var_names), 
            'variable(s):\n')
        cat(cat(x@var_names, sep = ' '), ".\n", sep = '')
        
        if (x@is_dynamic) {
            if (!x@ss_solved) {
                cat('\nThe steady state values for the variables have NOT',
                      'been found.')
            } else {
                if (!x@re_solved) {
                    add <- ' for variables have been found.'
                } else if (x@corr_computed) {
                    add <- paste(' and recursive laws of motion for variables have been found.',
                                 '\n \nThe correlations have been computed.')
                } else {
                    add <- ' and recursive laws of motion for variables have been found.'
                }
                cat('\nThe steady state values',
                      add, sep = '')                
            }
        } else {
            if (!x@ss_solved) {
                cat('\nThe equilibrium values for variables have NOT',
                      'been found.')
            } else {
                cat('\nThe equilibrium values for variables have',
                      'been found.')            
            }        
        
        }
        cat("\n")
        cat("\n----------------------------------------------------------", "\n")
        show(x)
    }
)

# ###################################################################
# The summary method prints the full information about
# the specified variables
# ###################################################################
# Input 
#   object - an object of class gecon_var_info
# Output
#   Prints the full information about the specified variables
# ###################################################################
setMethod(
    "summary",
    signature(object = "gecon_var_info"), 
    function(object, ...) 
    {
        print(object)
    }
)

