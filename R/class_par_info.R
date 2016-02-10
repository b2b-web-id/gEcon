# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Class for storing information about model parameters
# ###################################################################

# ###################################################################
# Class definition
# ###################################################################

setClass(
    Class = "gecon_par_info",
    representation = representation(model_info = "character",
                                    model_variable_name = "character",
                                    par_names = "character",
                                    gcn_values = "numeric",
                                    current_values = "numeric",
                                    calibr_flag = "logical",
                                    incid_mat = "Matrix")
)

# ###################################################################
# The function gecon_par_info is 
# a gecon_par_info class object constructor.
# ###################################################################
# Input 
#   model_info - a character vector of length 3, containing
#                information about the model: the input file name,
#                the input file path and the date of creation.
#   model_variable_name - a character denoting the name of variable 
#                         storing model for which the simulations 
#                         have been created.
#   par_names - a character vector of the parameter names.
#   gcn_values - a numeric vector with the values of free 
#                parameters specified in the .gcn file.
#   current_values - a numeric vector with the current values
#                    of the parameters.
#   calibr_flag - a logical vector with the length equal 
#                 to the number of the parameters. 
#                 The TRUE entries denote that a corresponding
#                 parameter is a calibrated parameter.
#   incid_mat - a Matrix object with the mapping of parameters 
#               into equations and calibrating equations.
# Output
#   An object of the gecon_par_info class.
# ###################################################################
gecon_par_info <- function(model_info,
                                 model_variable_name,
                                 par_names,
                                 gcn_values,
                                 current_values,
                                 calibr_flag,
                                 incid_mat) 
{
    par_info_object <- new('gecon_par_info')
     
    if (!is.character(model_info)) {
        stop('model_info should be of vector class')
    } else par_info_object@model_info <- model_info   

    if (!is.character(model_variable_name)) {
        stop('model_variable_name should be of vector class')
    } else par_info_object@model_variable_name <- model_variable_name   

    if (!is.character(par_names)) {
        stop('par_names should be of character class')
    } else par_info_object@par_names <- par_names  

    if (!is.numeric(gcn_values) & !is.logical(gcn_values)) {
        stop('gcn_values should be of numeric class')
    } else par_info_object@gcn_values <- gcn_values  

    if (!is.numeric(current_values) & !is.logical(current_values)) {
        stop('current_values should be of numeric class')
    } else par_info_object@current_values <- as.numeric(current_values)  

    if (!is.logical(calibr_flag)) {
        stop('calibr_flag should be of logical class')
    } else par_info_object@calibr_flag <- calibr_flag 

    if (!inherits(incid_mat, 'Matrix')) {
        stop('incid_mat should be of Matrix class')
    } else par_info_object@incid_mat <- incid_mat 

    return(par_info_object)
}

# ###################################################################
# The show method controls how the gecon_par_info
# object is printed.
# ###################################################################
# Input 
#   object - an object of class gecon_par_info.
# Output
#   Information parameter values and incidence in equations.
# ################################################################### 
setMethod(
    "show", signature(object = "gecon_par_info"),
    function(object) {
         # Printing the incidence matrix
        ch_incid_mat <- matrix(as.character(object@incid_mat), 
                                  nrow = nrow(object@incid_mat),
                                  ncol = ncol(object@incid_mat))
                                  
        colnames(ch_incid_mat) <- colnames(object@incid_mat)
        rownames(ch_incid_mat) <- rownames(object@incid_mat)        

        ch_incid_mat[which(object@incid_mat == 0)] <- '  .'
        ch_incid_mat[which(object@incid_mat == 1)] <- '  X'
        
        cat('\nIncidence info:\n\n')
        print(as.data.frame(ch_incid_mat)) 
        
        # Printing the information about parameters
        calibr_indic <- rep('Free     ', length(object@gcn_values))
        calibr_indic[object@calibr_flag] <- 'Calibrated  '
   
        gcn_val <- as.character(round(object@gcn_values, digits = 4))
        gcn_val[which(is.na(gcn_val))] <- '.'
        
        
        current_val <- as.character(round(object@current_values, digits = 4))
        current_val[which(is.na(current_val))] <- '.'
        
        cat("\n----------------------------------------------------------", "\n")       
        cat('\nParameter info: \n\n')           
        tab   <- cbind(gcn_val,
                       current_val,
                       calibr_indic)   
        rownames(tab) <- object@par_names
        colnames(tab) <- c('gcn file value', 'Current value', 'Parameter type')
        print(as.data.frame(tab))       
    }
)


# ###################################################################
# The print method prints information the gecon_par_info
# object.
# ###################################################################
# Input 
#   x - an object of the gecon_par_info class.
# Output
#   Information about the gecon_par_info object.
# ###################################################################
setMethod(
    "print",
    signature(x = "gecon_par_info"), 
    function(x) 
    {
        cat("\n")
        cat('The information about parameters has been generated\n',
            'based on model stored in variable: ', 
            x@model_variable_name, '.\n', sep = '')
        cat("\n")
        cat('The model has been created based on:', 
            x@model_info[1], 'file.\n')
        cat("\n")  
        
        n_free <- length(x@par_names) - sum(x@calibr_flag)
        
        cat('The information has been created for:\n')
        if (sum(x@calibr_flag)) {
            cat(sum(x@calibr_flag), 'calibrated parameter(s): \n')
            cat(cat(x@par_names[which(x@calibr_flag)], sep = ' '), ',\n', sep = '')            
        } else cat(sum(x@calibr_flag), 'calibrated parameters,\n')
        
        if (n_free) {
            cat(n_free , 'free parameter(s):\n')
            cat(cat(x@par_names[which(!x@calibr_flag)], sep = ' '), '.\n', sep = '')
        } else cat( cat(n_free), 'free parameters.\n')
        cat("----------------------------------------------------------", "\n")
        show(x)
    }
)



# ###################################################################
# The summary method displays the full information about the 
# parameters.
# ###################################################################
# Input 
#   object - an object of class gecon_par_info.
# Output
#   Prints the full information about the specified parameters.
# ###################################################################
setMethod(
    "summary",
    signature(object = "gecon_par_info"), 
    function(object, ...) 
    {
        print(object)
    }
)

