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
    Class = "gecon_shock_info",
    representation = representation(model_info = "character",
                                    model_variable_name = "character",
                                    shock_names = "character",
                                    shock_matrix = "matrix",
                                    shock_matrix_flag = "logical",
                                    incid_mat = "Matrix")
)

# ###################################################################
# The function gecon_shock_info is a gecon_shock_info class
# object constructor
# ###################################################################
# Input 
#   model_info - a character vector of length 3, containing 
#                information about the model: the input file name,
#                the input file path and the date of creation.
#   model_variable_name - a character denoting the name of variable 
#                         storing model for which the simulations 
#                         have been created.
#   shock_names - a character vector of the shock names.
#   shock_matrix - a matrix object with a variance-covariance matrix
#                  of specified shocks and all the model shocks.
#   shock_matrix_flag - a logical value, set to TRUE when user defines
#                       non-default entries in variance-covariance
#                       matrix of shocks.
#   incid_mat - a Matrix object with the mapping of shocks
#               into equations.
# Output
#   An object of the gecon_shock_info class
# ###################################################################
gecon_shock_info <- function(model_info,
                             model_variable_name,
                             shock_names,
                             shock_matrix,
                             shock_matrix_flag,
                             incid_mat)
{
    shock_info_object <- new('gecon_shock_info')
     
    if (!is.character(model_info)) {
        stop('model_info should be of vector class')
    } else shock_info_object@model_info <- model_info   

    if (!is.character(model_variable_name)) {
        stop('model_variable_name should be of vector class')
    } else shock_info_object@model_variable_name <- model_variable_name   

    if (!is.character(shock_names)) {
        stop('shock_names should be of character class')
    } else shock_info_object@shock_names <- shock_names  

    if (!is.matrix(shock_matrix)) {
        stop('shock_matrix should be of matrix class')
    } else shock_info_object@shock_matrix <- shock_matrix  

    if (!is.logical(shock_matrix_flag)) {
        stop('shock_matrix_flag should be of logical class')
    } else shock_info_object@shock_matrix_flag <- shock_matrix_flag  

    if (!inherits(incid_mat, 'Matrix')) {
        stop('incid_mat should be of Matrix class')
    } else shock_info_object@incid_mat <- incid_mat 

    return(shock_info_object)
}

# ###################################################################
# The show method controls how the gecon_shock_info
# object is printed
# ###################################################################
# Input 
#   object - an object of class gecon_shock_info
# Output
#   Information about shocks incidence and variance - covariance 
#   matrix.
# ################################################################### 
setMethod(
    "show", signature(object = "gecon_shock_info"),
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
        if (object@shock_matrix_flag) {
            cat("\n----------------------------------------------------------", "\n")       
            cat('\nVariance - covariance matrix of shocks: \n\n')           
            print(as.data.frame(object@shock_matrix))  
        } else {
            cat("\n----------------------------------------------------------", "\n")       
            cat(paste('\nVariance - covariance matrix of shocks has NOT been specified.',
                      'Default identity matrix is used for the stochastic simulations of model.'))                  
        }        
    }
)

# ###################################################################
# The print method prints information about shocks
# ###################################################################
# Input 
#   x - an object of the gecon_shock_info class.
# Output
#   Information about the gecon_shock_info object.
# ###################################################################
setMethod(
    "print",
    signature(x = "gecon_shock_info"), 
    function(x) 
    {       
        cat("\n")
        cat('The information about shocks has been generated\n',
            'based on model stored in variable: ', 
            x@model_variable_name, '.\n', sep = '')
        cat("\n")
        cat('The model has been created based on:', 
            x@model_info[1], 'file. \n')
        cat("\n")
        
        cat('The information has been created for',
             length(x@shock_names), 'shock(s):\n')
        cat(cat(x@shock_names, sep = ' '), ".\n", sep = '')     
        cat("\n----------------------------------------------------------", "\n")
        show(x)      
    }
)


# ###################################################################
# The summary method displays full information about
# the specified shocks
# ###################################################################
# Input 
#   object - an object of class gecon_shock_info
# Output
#   Prints the full information about the specified shocks
# ###################################################################
setMethod(
    "summary",
    signature(object = "gecon_shock_info"), 
    function(object, ...) 
    {
        print(object)
    }
)

