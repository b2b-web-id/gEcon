# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################
# Model shocks
# ############################################################################


# ############################################################################
# The get_shock_names function returns names of shocks
# stored in an object of gecon_model class
# ############################################################################
# Input
#   model - an object of the gecon_model class
# Output
#   returns a vector of shock names
# ############################################################################
get_shock_names <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    return (model@shocks)
}

# ############################################################################
# The get_shock_names_by_index function returns
# the names of shocks indexed by the indices specified
# in the index_names argument
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   index_names - a character vector with the names of indices.
# Output
#   A character vector with relevant shock names.
# ############################################################################
get_shock_names_by_index <- function(model, index_names)
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

    # Finding relevant parameters
    index_names_alternative <- paste("(", paste(index_names,  collapse = "|"), ")", collapse = "", sep = "")
    assoc_names  <- grep(pattern = paste("__", index_names_alternative,
                                         "(__)*(?!(_?[a-zA-Z0-9]))", sep = ""),
                                x = model@shocks,
                                perl = T)

    if (!length(assoc_names))
        stop("no shocks corresponding to given indices found")

    return (model@shocks[assoc_names])
}


# ############################################################################
# The set_shock_distr_par function assigns distribution parameters
# (standard deviations, correlations, etc.) to shocks in
# an object of gecon_model class.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   distr_par - a list or vector of distribution parameters.
# Output
#   A gecon_model object with shocks' distribution parameters set
#   to values given by the user.
# ############################################################################
set_shock_distr_par <- function(model, distr_par)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    distr_par_names <- names(distr_par)
    if(is.null(distr_par_names))
        stop("distr_par argument has to be a list or a vector with named elements")

    is_proper_numeric <- tryCatch(expr = {distr_par <- as.numeric(distr_par)
                                            TRUE
                                       },
                               warning = function(w) FALSE,
                               error = function(w) FALSE)

    if (!is_proper_numeric) {
        stop(paste("the following elements of distr_par cannot be coerced into distribution parameters' values:",
                    paste(distr_par_names[union(which(as.logical(lapply(distr_par, is.character))),
                                         which(as.numeric(lapply(distr_par, length)) > 1))], collapse = ", ")))
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
        labs <- rep("Inf", length(misspec_ind))
        labs[which(is.nan(distr_par[misspec_ind]))] <- "NaN"
        labs[which(is.na(distr_par[misspec_ind]))] <- "NA"

        stop("values of the following distribution parameters have been incorrectly specified: \n",
             paste(paste0("\"", misspec_names, "\"", paste0(" (=", labs, ")")),
                   collapse = ", "))
    }

    par_type <- rep(0, length(distr_par))

    par_type <- grepl("^sd\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) +
                grepl("^var\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 2 +
                grepl("^cor\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*,\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 3 +
                grepl("^cov\\(\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*,\\s*[a-zA-Z]([_a-zA-Z0-9])*[_a-zA-Z0-9]\\s*\\)$", distr_par_names, perl = T) * 4


    if (any(par_type == 0))
        stop(paste("the following are not valid shock distribution parameters:\n",
                   list2str(distr_par_names[which(par_type == 0)], "\"")))

    single_var <- which(par_type %in% c(1, 2))
    mult_var <- which(par_type %in% c(3, 4))

    # validating values
    if (any(distr_par[single_var] < 0))
        stop("variances and standard deviations cannot be negative")

    if (any(distr_par[which(par_type == 3)] > 1 | distr_par[which(par_type == 3)] < -1))
        stop("correlations have to lie in the interval [-1, 1]")

    # retrieving shock names
    shock_names_matrix <- matrix(data = character(), length(distr_par_names), 2)

    shock_names <- regexpr("(?<=\\()(_?[\\s,_a-zA-Z0-9])*(?=\\))", distr_par_names, perl = T)
    shock_names <- regmatches(distr_par_names, shock_names)
    shock_names <- gsub("\\s", "", shock_names)

    shock_names_matrix[single_var, ] <- shock_names[single_var]

    name_vec <- regexpr("(?<=,)(_?[\\s,_a-zA-Z0-9])*", shock_names[mult_var], perl = T)
    shock_names_matrix[mult_var, 1] <- regmatches(shock_names[mult_var], name_vec)

    name_vec <- regexpr("(_?[\\s,_a-zA-Z0-9])*(?=,)", shock_names[mult_var], perl = T)
    shock_names_matrix[mult_var, 2] <- regmatches(shock_names[mult_var], name_vec)

    shock_names_matrix <- matrix(match(shock_names_matrix, model@shocks),
                                 length(distr_par_names), 2)

    # validation of inputs
    wrong_names <- unique(which(is.na(shock_names_matrix), arr.ind = 1)[, 1])

    if(length(wrong_names))
        stop(paste("the following entries do not correspond to any shock distribution parameters:",
                   list2str(distr_par_names[wrong_names], "\"")))

    complying_ind <- shock_names_matrix[, 1] == shock_names_matrix[, 2]

    corr_err <- which(par_type == 3 & complying_ind)
    var_spec <- which(par_type == 4 & complying_ind)

    if (length(corr_err))
        stop(paste("correlations have to be specified for pairs of different shocks; the following are not valid specifications: ",
                   list2str(distr_par_names[corr_err],  "\"")))

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
        stop(paste("the following covariances/correlations have been specified more than once:",
                    paste("\"", model@shocks[double_entries[, 1]], "\" and \"", model@shocks[double_entries[, 2]], "\"", collapse = ", ")))
    duplicates <- duplicated(shock_names_matrix[which(par_type < 3), 1])
    duplicates <- model@shocks[unique(shock_names_matrix[duplicates, 1])]
    if (length(duplicates))
        stop(paste("variance/standard deviation has been specified more than once for the following shocks:",
                   list2str(duplicates, "\"")))

    # updating shock matrix
    cov_matrix <- model@shock_cov_mat
    sd_init <- sqrt(diag(cov_matrix))
    correlation_matrix <- matrix(0, dim(model@shock_cov_mat)[1], dim(model@shock_cov_mat)[1])
    correlation_matrix[model@active_shocks, model@active_shocks] <- cov_matrix[model@active_shocks, model@active_shocks] /
                            kronecker(sd_init, t(sd_init))[model@active_shocks, model@active_shocks]
    diag(correlation_matrix) <- 1

    # updating diagonal of the covariance matrix
    distr_par[which(par_type == 1)] <- distr_par[which(par_type == 1)] ^ 2
    cov_matrix[cbind(shock_names_matrix[which(par_type < 3)], shock_names_matrix[par_type < 3])] <-
                                                distr_par[which(par_type < 3)]
    sd_fin <- sqrt(diag(cov_matrix))

    # updating correlation matrix
    col_ind <-  c(shock_names_matrix[which(par_type == 3), 1],
                      shock_names_matrix[which(par_type == 3), 2])
    row_ind <-  c(shock_names_matrix[which(par_type == 3), 2],
                      shock_names_matrix[which(par_type == 3), 1])
    correlation_matrix[(col_ind - 1) * dim(model@shock_cov_mat)[1] + row_ind] <-
                                               c(distr_par[which(par_type == 3)],
                                                 distr_par[which(par_type == 3)])
    cov_matrix <- kronecker(sd_fin, t(sd_fin)) *  correlation_matrix


    # updating covariances
    pure_covariances <- which(par_type == 4)

    col_ind <-  shock_names_matrix[pure_covariances, 1]
    row_ind <-  shock_names_matrix[pure_covariances, 2]

    cov_matrix[(col_ind - 1) * dim(model@shock_cov_mat)[1] + row_ind] <- distr_par[pure_covariances]
    cov_matrix[(row_ind - 1) * dim(model@shock_cov_mat)[1] + col_ind] <- distr_par[pure_covariances]

    active_shocks <- rep(TRUE, length(model@shocks))
    zero_shocks <- which(diag(cov_matrix) == 0)
    if (length(zero_shocks) == length(model@shocks))
        stop("at least one shock has to have variance greater than zero")

    if (length(zero_shocks)) {
        active_shocks[zero_shocks] <- FALSE
        cov_matrix[zero_shocks, zero_shocks] <- 0
        warning(paste("the following shocks will not be taken into account in model simulations: ",
                      list2str(model@shocks[zero_shocks], "\"")))
    }
    # check if matrix is positive definite
    is_pos_def <- tryCatch(expr = {chol(cov_matrix[active_shocks, active_shocks])
                                    TRUE
                                   },
                           warning = function(w) FALSE,
                           error = function(w) FALSE)

    if (!is_pos_def) {
        stop("covariance matrix is not positive definite")
    }

    model@active_shocks <- active_shocks
    model@shock_cov_mat <- cov_matrix
    model@shock_cov_mat_flag <- TRUE

    # clearing slots
    model@corr_computed <- FALSE
    model@corr_mat <- matrix(nrow = 0, ncol = 0)
    model@autocorr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_corr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_idx <- 0L
    model@var_dec <- matrix(nrow = 0, ncol = 0)
    model@sdev <- matrix(nrow = 0, ncol = 0)

    return (model)
}


# ############################################################################
# The set_shock_cov_mat function sets shock covariance matrix to one specified
# by the user.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   cov_matrix -  a positive definite matrix with dimensions
#                   (n x n) where n is the number of shocks
#                   in the model.
#   shock_order - a character vector declaring the ordering
#                 of shocks in the cov_matrix. If not specified,
#                 it is assumed that the ordering is in accordance
#                 with the internal ordering of model.
#                 The default ordering can be displayed
#                 with the shock_info function with
#                 the all_shocks argument set to TRUE.
# Output
#   A gecon_model object with the covariance matrix set to
#   the values given by the user.
# ############################################################################
set_shock_cov_mat <- function(model, cov_matrix, shock_order = NULL)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    if ((!is.matrix(cov_matrix) | !is.numeric(cov_matrix))) {
        stop("cov_matrix should be a numeric matrix")
    }

    if (!isSymmetric(cov_matrix)) {
        stop("covariance matrix has to be symmetric")
    }
    smd <- dim(cov_matrix)
    if (smd[1] != length(model@shocks)) {
        stop(paste("covariance matrix dimensions do not",
                   "agree with the number of shocks in the model"))
    }

    if (length(as.numeric(cov_matrix)) == 1 & (length(model@shocks) == 1)) {
        shock_order <- model@shocks
    }

    if (is.null(shock_order)) {
        shock_order <- model@shocks
        if (length(model@shocks) > 1) {
            warning(paste("shock_order argument has not been specified,",
                        "shocks will be sorted according to their internal ordering",
                        "(as reported by the \'get_shock_names\' function)"))
        }
    }

    if (!is.character(shock_order) | length(shock_order) != length(model@shocks)) {
            stop(paste("shock_order argument has to be a character vector",
                       "specifying the ordering of all shocks in the model"))
    }

    shock_perm <- match(model@shocks, shock_order)

    if (any(is.na(shock_perm)))
        stop("invalid shock name(s) in shock_order argument")

    cov_matrix <- cov_matrix[shock_perm, shock_perm, drop = FALSE]

    # shocks with zero variance
    active_shocks <- rep(TRUE, length(model@shocks))
    zero_shocks <- which(diag(cov_matrix) == 0)
    if (any(abs(cov_matrix[zero_shocks,]) != 0))
        stop("all covariances of shocks with zero variance have to be zero")

    if (length(zero_shocks) == length(model@shocks))
        stop("at least one shock has to have variance greater than zero")

    if (length(zero_shocks)) {
        active_shocks[zero_shocks] <- FALSE
        cov_matrix[zero_shocks, zero_shocks] <- 0
        warning(paste("the following shocks will not be taken into account in model simulations:",
                      list2str(model@shocks[zero_shocks], "\"")))
    }

    is_pos_def <- tryCatch(expr = {chol(cov_matrix[active_shocks, active_shocks])
                                    TRUE
                                  },
                            warning = function(w) FALSE,
                            error = function(w) FALSE)

    if (!is_pos_def) {
        stop("covariance matrix is not positive definite")
    }

    rownames(cov_matrix) <- model@shocks
    colnames(cov_matrix) <- model@shocks
    model@shock_cov_mat <- cov_matrix
    model@shock_cov_mat_flag <- TRUE
    model@active_shocks <- active_shocks

    # clearing slots
    model@corr_computed <- FALSE
    model@corr_mat <- matrix(nrow = 0, ncol = 0)
    model@autocorr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_corr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_idx <- 0L
    model@var_dec <- matrix(nrow = 0, ncol = 0)
    model@sdev <- matrix(nrow = 0, ncol = 0)

    return (model)
}




# ############################################################################
# The get_shock_cov_mat function returns covariance matrix of model shocks.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
# Output
#   The covariance matrix of model shocks.
# ############################################################################
get_shock_cov_mat <- function(model)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    shock_cov_mat <- model@shock_cov_mat
    rownames(shock_cov_mat) <- model@shocks
    colnames(shock_cov_mat) <- model@shocks

    return (shock_cov_mat)
}


# ############################################################################
# The shock_info function prints information about model shocks and allows to create
# an object of gecon_shock_info class.
# ############################################################################
# Input
#   model - an object of the gecon_model class.
#   shocks - the names or indices of the shocks of interest.
#   all - logical value. If set to TRUE, the shocks
#                    argument is overwritten with all the shocks
#                    appearing in the model.
#                    The default value is FALSE.
# Output
#   An object of gecon_shock_info class.
# ############################################################################
shock_info <- function(model, shocks = NULL, all = FALSE)
{
    if (!is.gecon_model(model))
        stop("model argument should be of gecon_model class")

    if (!model@is_stochastic)
        stop("model is deterministic, it does not have any shocks")

    if (is.null(shocks) & !all)
        stop("no shocks have been specified")

    if (all) {
        if (!is.null(shocks)) {
            warning("ignoring shocks argument when all option is set to TRUE")
        }
        shocks <- model@shocks
        shock_ind <- 1:length(shocks)
    } else {
        shock_ind <- list2ind(shocks, model@shocks, "shock")
        shocks <- model@shocks[shock_ind]
    }

    # Information about variable incidence
    shock_eq_map <- model@shock_eq_map[, shock_ind, drop = FALSE]
    incid_eq <- sort(unique(shock_eq_map@i)) + 1
    shock_eq_map <- shock_eq_map[incid_eq, , drop = FALSE]

    rownames(shock_eq_map) <- paste("Eq. ", incid_eq, sep = "")
    colnames(shock_eq_map) <- model@shocks[shock_ind]

    # Information about shocks correlation matrix
    shock_cov_mat <- model@shock_cov_mat[, shock_ind, drop = FALSE]
    rownames(shock_cov_mat) <- model@shocks
    colnames(shock_cov_mat) <- model@shocks[shock_ind]

    info <- gecon_shock_info(r_object_name = deparse(substitute(model)),
                             shocks = model@shocks[shock_ind],
                             cov_matrix = shock_cov_mat,
                             cov_matrix_flag = model@shock_cov_mat_flag,
                             incid_mat = shock_eq_map)

    return (info)
}
