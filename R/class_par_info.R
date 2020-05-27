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
# Class for storing information about model parameters
# ############################################################################

# ############################################################################
# Class definition
# ############################################################################

setClass(
    Class = "gecon_par_info",
    representation = representation(r_object_name = "character",
                                    parameters = "character",
                                    gcn_values = "numeric",
                                    current_values = "numeric",
                                    calibr_flag = "logical",
                                    calibr_init_val = "numeric",                                     
                                    incid_mat = "Matrix")
)

# ############################################################################
# The function gecon_par_info is
# a gecon_par_info class object constructor.
# ############################################################################
# Input
#   r_object_name - a character string denoting the name of a gecon_model R object
#                         storing model for which the simulations
#                         have been created.
#   parameters - a character vector of the parameters' names.
#   gcn_values - a numeric vector with the values of free
#                parameters specified in the .gcn file.
#   current_values - a numeric vector with the current values
#                    of the parameters.
#   calibr_flag - a logical vector with the length equal
#                 to the number of the parameters.
#                 The TRUE entries denote that a corresponding
#                 parameter is a calibrated parameter.
#   calibr_init_val - a numeric vector with initial values of calibrated parameters.                                                                                 
#   incid_mat - a Matrix object with the mapping of parameters
#               into equations and calibrating equations.
# Output
#   An object of the gecon_par_info class.
# ############################################################################
gecon_par_info <- function(r_object_name,
                           parameters,
                           gcn_values,
                           current_values,
                           calibr_flag,
                           calibr_init_val,                          
                           incid_mat)
{
    par_info_object <- new("gecon_par_info")

    if (!is.character(r_object_name)) {
        stop("r_object_name should be of character type")
    } else par_info_object@r_object_name <- r_object_name

    if (!is.character(parameters)) {
        stop("parameters should be of character type")
    } else par_info_object@parameters <- parameters

    if (!is.numeric(gcn_values) & !is.logical(gcn_values)) {
        stop("gcn_values should be of numeric type")
    } else par_info_object@gcn_values <- gcn_values

    if (!is.numeric(current_values) & !is.logical(current_values)) {
        stop("current_values should be of numeric type")
    } else par_info_object@current_values <- as.numeric(current_values)

    if (!is.logical(calibr_flag)) {
        stop("calibr_flag should be of logical type")
    } else par_info_object@calibr_flag <- calibr_flag

    if (!is.numeric(calibr_init_val) & !is.logical(calibr_init_val)) {
        stop("calibr_init_val should be of numeric type")
    } else par_info_object@calibr_init_val <- calibr_init_val
    
    if (!inherits(incid_mat, "Matrix")) {
        stop("incid_mat should be of Matrix class")
    } else par_info_object@incid_mat <- incid_mat

    return (par_info_object)
}

# ############################################################################
# Showing information stored in an object of gecon_par_info class.
# ############################################################################
# Input
#   object - an object of class gecon_par_info.
# Output
#   Information about parameters' values and incidence in equations.
# ############################################################################
setMethod("show", signature(object = "gecon_par_info"),
function(object) {
    # Printing the incidence matrix
    ch_incid_mat <- matrix(as.character(object@incid_mat),
                                nrow = nrow(object@incid_mat),
                                ncol = ncol(object@incid_mat))

    colnames(ch_incid_mat) <- colnames(object@incid_mat)
    rownames(ch_incid_mat) <- rownames(object@incid_mat)

    ch_incid_mat[which(object@incid_mat == 0)] <- "  ."
    ch_incid_mat[which(object@incid_mat == 1)] <- "  X"

    cat("Incidence info:\n\n")
    print(as.data.frame(ch_incid_mat))

    # Printing the information about parameters
    calibr_indic <- rep("Free     ", length(object@gcn_values))
    calibr_indic[object@calibr_flag] <- "Calibrated  "

    gcn_val <- as.character(round(object@gcn_values, digits = 4))
    gcn_val[which(is.na(gcn_val))] <- "."


    current_val <- as.character(round(object@current_values, digits = 4))
    current_val[which(is.na(current_val))] <- "."

    init_val <- as.character(round(object@calibr_init_val, digits = 4))
    init_val[which(is.na(init_val))] <- "."
    cat("\n----------------------------------------------------------", "\n")
    cat("\nParameter info: \n\n")
    tab <- cbind(gcn_val, current_val, init_val, calibr_indic)
    rownames(tab) <- object@parameters
    colnames(tab) <- c(".gcn file value", "Current value", "Calibr.initial value", 
                        "Parameter type")
    print(as.data.frame(tab))
})


# ############################################################################
# Showing information stored in an object of gecon_par_info class.
# ############################################################################
# Input
#   object - an object of class gecon_par_info.
# Output
#   Information about parameters' values and incidence in equations.
# ############################################################################
setMethod("print", signature(x = "gecon_par_info"),
function(x, ...)
{
    show(x)
})


# ############################################################################
# Showing information stored in an object of gecon_par_info class.
# ############################################################################
# Input
#   object - an object of class gecon_par_info.
# Output
#   Information about parameters' values and incidence in equations.
# ############################################################################
setMethod("summary", signature(object = "gecon_par_info"),
function(object)
{
    print(object)
})
