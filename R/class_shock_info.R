# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################
# Class for storing information about model shocks
# ############################################################################

# ############################################################################
# Class definition
# ############################################################################

setClass(
    Class = "gecon_shock_info",
    representation = representation(r_object_name = "character",
                                    shocks = "character",
                                    cov_matrix = "matrix",
                                    cov_matrix_flag = "logical",
                                    incid_mat = "Matrix")
)

# ############################################################################
# The function gecon_shock_info is a gecon_shock_info class
# object constructor
# ############################################################################
# Input
#   r_object_name - a character string denoting the name of an R object
#                   of \code{gecon_model} class storing the model for which
#                   the simulations have been performed.
#   shocks - a character vector of the shocks' names.
#   cov_matrix - a numeric matrix containing columns of shock covariance matrix
#                corresponding to selected shocks.
#   cov_matrix_flag - a logical value, set to TRUE when the user enter
#                       non-default data into a covariance
#                       matrix of shocks.
#   incid_mat - a object of Matrix class representing the mapping of shocks
#               into equations.
# Output
#   An object of the gecon_shock_info class
# ############################################################################
gecon_shock_info <- function(r_object_name,
                             shocks,
                             cov_matrix,
                             cov_matrix_flag,
                             incid_mat)
{
    shock_info_object <- new("gecon_shock_info")

    if (!is.character(r_object_name)) {
        stop("r_object_name should be of character type")
    } else shock_info_object@r_object_name <- r_object_name

    if (!is.character(shocks)) {
        stop("shocks should be of character type")
    } else shock_info_object@shocks <- shocks

    if (!is.matrix(cov_matrix)) {
        stop("cov_matrix should be of matrix type")
    } else shock_info_object@cov_matrix <- cov_matrix

    if (!is.logical(cov_matrix_flag)) {
        stop("cov_matrix_flag should be of logical type")
    } else shock_info_object@cov_matrix_flag <- cov_matrix_flag

    if (!inherits(incid_mat, "Matrix")) {
        stop("incid_mat should be of Matrix class")
    } else shock_info_object@incid_mat <- incid_mat

    return (shock_info_object)
}

# ############################################################################
# Showing information stored in an object of gecon_shock_info class.
# ############################################################################
# Input
#   object - an object of class gecon_shock_info
# Output
#   Information about shocks incidence and covariance matrix printed on the console.
# ############################################################################
setMethod("show", signature(object = "gecon_shock_info"),
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
    if (object@cov_matrix_flag) {
        cat("\n----------------------------------------------------------", "\n")
        cat("\nCovariance matrix of shocks: \n\n")
        print(as.data.frame(object@cov_matrix))
    } else {
        cat("\n----------------------------------------------------------", "\n")
        cat(paste("\nCovariance matrix of shocks has NOT been specified.",
                    "By default identity matrix is used for stochastic simulations of the model."))
    }
})

# ############################################################################
# Showing information stored in an object of gecon_shock_info class.
# ############################################################################
# Input
#   object - an object of class gecon_shock_info
# Output
#   Information about shocks incidence and covariance matrix printed on the console.
# ############################################################################
setMethod("print", signature(x = "gecon_shock_info"),
function(x, ...)
{
    show(x)
})


# ############################################################################
# Showing information stored in an object of gecon_shock_info class.
# ############################################################################
# Input
#   object - an object of class gecon_shock_info
# Output
#   Information about shocks incidence and covariance matrix printed on the console.
# ############################################################################
setMethod("summary", signature(object = "gecon_shock_info"),
function(object)
{
    show(object)
})

