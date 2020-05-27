# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Model statistics and simulation
# ############################################################################


# ############################################################################
# The function computes statistics of the model using spectral (FFT) or
# simulation methods.
# ############################################################################
# Input
#   model - an object of gecon_model class.
#   n_leadlags - the number of leads and lags of the model's variables for which
#                correlations are to be computed.
#   ref_var - the name or the index of the reference variable with respect to which
#             correlations are to be computed.
#   lambda - HP filter lambda, if it is set to 0 no filtering is performed,
#            1600 is default value (quarterly data).
#   ngrid - a numeric value, the density of grid used by the Fast Fourier transform
#           (used only if the \code{sim} option is set to FALSE). It has to be a multiple
#           of 8 and has to be large enough to guarantee unbiased results.
#   sim - a logical value, if TRUE simulation is used for computing correlations, if FALSE,
#         the Fast Fourier transform is used.
#   sim_length - the length of simulation path (used only if the \code{sim} option
#                is set to TRUE).
# Output
#   An object of gecon_model class representing the model.
# ############################################################################
compute_model_stats <- function(model, n_leadlags = 5, ref_var = NULL, lambda = 1600,
                                ngrid = 64 * 16,
                                sim = FALSE, sim_length = 1e5)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }
    if (!model@is_stochastic) stop("the model is not stochastic")
    if (!model@re_solved) stop("solve the 1st order perturbation first using 'solve_pert' function")
    if (!length(model@state_var_indices)) stop("the model does not have any state variables")
    if (!model@shock_cov_mat_flag)
        warning(paste0("by default, identity matrix has been assumed as the shock covariance matrix, ",
                       "covariance matrix can be set through 'set_shock_cov_mat' ",
                       "or 'set_shock_distr_par' functions."))

    # PP and QQ matrices
    nvar <- length(model@variables)
    PP <- matrix(0, nvar, nvar)
    PP[model@state_var_indices, model@state_var_indices] <- model@solution$P
    PP[-c(model@state_var_indices), model@state_var_indices] <- model@solution$R
    QQ <- matrix(0, nvar, sum(model@active_shocks))
    QQ[model@state_var_indices, ] <- model@solution$Q[, model@active_shocks]
    QQ[-c(model@state_var_indices), ] <- model@solution$S[, model@active_shocks]

    if (!is.null(ref_var)) {
        if ((length(ref_var) > 1) || (!is.numeric(ref_var) && !is.character(ref_var))) {
            stop("invalid 'ref_var' argument", call. = FALSE)
        } else if (is.numeric(ref_var)) {
            if ((ref_var < 1) || (ref_var > nvar)) {
                stop("invalid variable index", call. = FALSE)
            }
            ref_var_idx <- as.integer(ref_var)
        } else {
            ref_var_idx <- which(model@variables == ref_var)
            if (!length(ref_var_idx)) {
                stop(paste0("\"", ref_var, "\" is not a model variable"), call. = FALSE)
            }
        }
    } else {
        ref_var_idx <- 0L
    }
    model@ref_var_idx <- ref_var_idx

    # Computing covariance matrix
    if (!sim) {
        output <- covar_fft(PP, QQ,
                            model@shock_cov_mat[model@active_shocks,
                                                model@active_shocks,
                                                drop = FALSE],
                            n_leadlags, ref_var_idx, lambda, ngrid)
    } else {
        output <- covar_sim(PP, QQ,
                            model@shock_cov_mat[model@active_shocks,
                                                model@active_shocks,
                                                drop = FALSE],
                            n_leadlags, ref_var_idx, lambda, sim_length)
    }
    cov_mat<- output$covmat
    acor <- output$acor
    ref_var_cor <-  output$ref_var_cor

    # standard deviations
    sdev <- sqrt(diag(cov_mat))
    var <- kronecker(sdev, t(sdev))
    model@sdev <- matrix(sdev, nvar, 1)
    rownames(model@sdev) <- model@variables

    # correlation matrix
    model@corr_mat <- cov_mat / var
    rownames(model@corr_mat) <- model@variables
    colnames(model@corr_mat) <- model@variables

    #autocorrelations
    colnames(acor) <- paste0("Lag ", 1:n_leadlags)
    rownames(acor) <- model@variables
    model@autocorr_mat <- acor

    # correlation with specified variable
    if (ref_var_idx) {
        if (n_leadlags) {
            cnames <- c((-n_leadlags):(-1), "", 1:n_leadlags)
        } else {
            cnames <- ""
        }
        cnames <- paste0(model@variables[ref_var_idx], "[", cnames, "]")
        colnames(ref_var_cor) <- cnames
        rownames(ref_var_cor) <- paste0(model@variables, "[]")
        model@ref_var_corr_mat <- ref_var_cor
    } else {
        model@ref_var_corr_mat <- matrix(nrow = 0, ncol = 0)
    }

    if (sim) {
        warning(paste0("when the 'sim' option is set to TRUE (statistics are determined ",
                       "through simulation method) variance decomposition is not computed"))
    } else {
        # Computing vardecomp
        PPs <- PP[, model@state_var_indices]
        vardec <- vardecomp_fft(model@solution$P,
                                model@shock_cov_mat[model@active_shocks,
                                                    model@active_shocks,
                                                    drop = FALSE],
                                model@solution$Q[, model@active_shocks, drop = FALSE],
                                PPs, QQ,
                                lambda, ngrid)
        rownames(vardec) <- model@variables
        colnames(vardec) <- model@shocks[model@active_shocks]
        model@var_dec <- vardec
    }

    model@corr_computed <- TRUE

    return (model)
}


# ############################################################################
# The get_model_stats function prints and returns the statistics of the model
# (absolute and relative to the reference variable)
# ############################################################################
# Input
#   model - an object of gecon_model class.
#   variables - the names or indices of the variables of interest.
#   basic_stats - a logical value, if TRUE, the following information is returned
#                                  for selected variables: steady-state value,
#                                  standard deviation, variance, and information
#                                  whether a variable has been log-linearised.
#   corr - a logical value, if TRUE, a correlation matrix is returned. If a reference
#                           variable was not NULL while invoking the 'compute_model_stats'
#                           function, then correlations of selected variables
#                           with leads and lags of the reference variable are also returned.
#   autocorr - a logical value, if TRUE, autocorrelations of selected variables are returned.
#   var_dec - a logical value, if TRUE, variance decomposition (shocks' contributions
#                              to the variables' variances) is returned.
#   to_tex - logical, if TRUE results are written to a .tex file
#   silent - logical, if TRUE, console output is suppressed (FALSE by default).
#
# Output
#   a list of model statistics
# ############################################################################
get_model_stats <- function(model,
                            variables = NULL,
                            basic_stats = TRUE,
                            corr = TRUE,
                            autocorr = TRUE,
                            var_dec = TRUE,
                            to_tex = FALSE,
                            silent = FALSE)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }
    if (!model@is_stochastic) stop("the model is not stochastic")
    if (!model@corr_computed) stop("compute model statistics first using the 'compute_model_stats' function")

    if (is.null(variables)) {
        variables <- model@variables
        var_ind <- 1:length(variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
        variables <- model@variables[var_ind]
    }

    if (var_dec && !(dim(model@var_dec)[2])) {
        stop("variance decomposition has not been determined")
    }
    nzind <- which(model@sdev != 0)
    if (!length(intersect(nzind, var_ind))) {
        if (any(corr, autocorr, var_dec)) {
            warning("all selected variables have zero variance", call. = FALSE)
        }
        corr <- FALSE
        autocorr <- FALSE
        var_dec <- FALSE
    }
    if (!any(basic_stats, corr, autocorr, var_dec)) {
        stop("no information selected")
    }

    refv_idx <- model@ref_var_idx
    res <- list()
    if (to_tex) {
        texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
        tex <- "\\section{Model statistics}\n"
    }

    if (basic_stats) {
        loglin_ind <- rep("Y   ", length(model@loglin_var[var_ind]))
        loglin_ind[which(!model@loglin_var[var_ind])] <- "N   "
        loglin_ind <- as.data.frame(loglin_ind, stringsAsFactors = FALSE)
        basic_stats <- round(cbind(model@variables_ss_val[var_ind],
                                   model@sdev[var_ind],
                                   model@sdev[var_ind] ^ 2),
                             digits = 4)
        basic_stats <- as.data.frame(basic_stats)
        basic_stats <- cbind(basic_stats, loglin_ind)
        
        if (model@ss_solved) {
            colnames(basic_stats) <- c("Steady-state value", "Std. dev.",
                                       "Variance", "Loglin")
        } else if (!all(is.na(model@ss_val))) {
            colnames(basic_stats) <- c("Last solver iteration value", "Std. dev.",
                                       "Variance", "Loglin")
        }
        rownames(basic_stats) <- model@variables[var_ind]

        res <- c(res, list(basic_stats = basic_stats))

        if (!silent) {
            cat("Basic statistics:\n\n")
            print(basic_stats)
        }
        if (to_tex) {
            tex <- paste0(tex, "\n\\subsection{Basic statistics}\n")
            rownames(basic_stats) <- paste0("$", model@variables_tex[var_ind], "$")
            tex <- paste0(tex, df2textable(basic_stats))
        }
    }

    var_ind <- intersect(nzind, var_ind)
    N <- length(var_ind)

    if (corr) {
        corr <- model@corr_mat[var_ind, var_ind, drop = FALSE]
        res <- c(res, list(corr = corr))
        if (refv_idx) {
            corr_refvar <- model@ref_var_corr_mat[var_ind, , drop = FALSE]
            corr_refvar <- cbind(model@sdev[var_ind] / model@sdev[refv_idx], corr_refvar)
            colnames(corr_refvar)[1] <- paste0("Std. dev. rel. to ", model@variables[refv_idx])
            res <- c(res, list(corr_refvar = corr_refvar))
        }

        if (!silent || to_tex) {
            corr <- round(corr, digits = 3)
            corr <- as.data.frame(corr)
            for (j in (1:N)) {
                if (j < N) {
                    corr[(j + 1):N, j] <- ""
                }
            }
        }

        if (!silent) {
            cat("\nCorrelation matrix:\n\n")
            print(corr)
            if (refv_idx) {
                cat(paste0("\nCross correlations with the reference variable (",
                           model@variables[refv_idx], "):\n\n"))
                print(round(corr_refvar, digits = 3))
            }
        }
        if (to_tex) {
            tex <- paste0(tex, "\n\\subsection{Correlation matrix}\n")
            rownames(corr) <- paste0("$", model@variables_tex[var_ind], "$")
            colnames(corr) <- paste0("$", model@variables_tex[var_ind], "$")
            tex <- paste0(tex, df2textable(corr, FALSE))
            if (refv_idx) {
                tex <- paste0(tex, "\n\\subsection{Cross correlations with the reference variable ($",
                              model@variables_tex[refv_idx], "$)}\n")
                corr_refvar <- round(corr_refvar, digits = 3)
                rownames(corr_refvar) <- paste0("$", model@variables_tex[var_ind], "_{t}$")
                nl <- dim(model@ref_var_corr_mat)[2]
                nl <- (nl - 1) / 2
                cn <- c((-nl):(-1), "", paste0("+", 1:nl))
                cn <- paste0("$", model@variables_tex[refv_idx], "_{t", cn, "}$")
                cn <- c(paste0("$\\sigma[\\cdot]$ rel. to $\\sigma[",
                               model@variables_tex[refv_idx], "]$"),
                        cn)
                colnames(corr_refvar) <- cn
                tex <- paste0(tex, df2textable(corr_refvar))
            }
        }
    }

    if (autocorr) {
        autocorr <- model@autocorr_mat[var_ind, , drop = FALSE]
        res <- c(res, list(autocorr = autocorr))
        if (!silent) {
            cat("\nAutocorrelations:\n\n")
            print(round(autocorr, digits = 3))
        }
        if (to_tex) {
            tex <- paste0(tex, "\n\\subsection{Autocorrelations}\n")
            autocorr <- round(autocorr, digits = 3)
            rownames(autocorr) <- paste0("$", model@variables_tex[var_ind], "$")
            tex <- paste0(tex, df2textable(autocorr, FALSE))
        }
    }

    if(var_dec) {
        var_dec <- model@var_dec[var_ind, , drop = FALSE]
        res <- c(res, list(var_dec = var_dec))
        if (!silent) {
            cat("\nVariance decomposition:\n\n")
            print(round(var_dec, digits = 3))
        }
        if (to_tex) {
            tex <- paste0(tex, "\n\\subsection{Variance decomposition}\n")
            var_dec <- round(var_dec, digits = 3)
            rownames(var_dec) <- paste0("$", model@variables_tex[var_ind], "$")
            colnames(var_dec) <- paste0("$", model@shocks_tex[model@active_shocks], "$")
            tex <- paste0(tex, df2textable(var_dec, FALSE))
        }
    }

    if (to_tex) {
        write_tex(texf, tex, !silent)
    }

    return (invisible(res))
}



# ############################################################################
# The function generates random shock paths based on the shock covariance matrix
# specified by the user and simulates the behaviour of the system.
# ############################################################################
# Input
#   model - an object of gecon_model class
#   variables - the names or indices of variables whose paths are to be simulated
#   sim_length - the length of simulation path, the default value is 40
# Output
#   an object of gecon_simulation class
# ############################################################################
random_path <- function(model, variables = NULL, sim_length = 40)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }
    if (!model@is_stochastic) stop("the model is not stochastic")
    if (!model@re_solved) stop("solve the 1st order perturbation first using the 'solve_pert' function")
    if (!model@shock_cov_mat_flag)
        warning(paste0("by default, identity matrix has been assumed as the shock covariance matrix, ",
                       "covariance matrix can be set through 'set_shock_cov_mat' ",
                       "or 'set_shock_distr_par' functions."))

    if (is.null(variables)) {
        variables <- model@variables
        var_ind <- 1:length(variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
        variables <- model@variables[var_ind]
    }

    shock_ind <- which(model@active_shocks)
    n_shock <- length(shock_ind)
    A <- t(chol(model@shock_cov_mat[shock_ind, shock_ind]))
    eps <- matrix(rnorm(sim_length * n_shock), nrow = n_shock, ncol = sim_length)
    eps <- A %*% eps
    shock_path <- matrix(0, nrow = length(model@shocks), ncol = sim_length)
    shock_path[shock_ind, ] <- eps

    res <- simulate_path(P = model@solution$P,
                         Q = model@solution$Q,
                         R = model@solution$R,
                         S = model@solution$S,
                         sim_length = sim_length,
                         init_st = NULL,
                         shock_path = shock_path,
                         var_ind = var_ind,
                         state_var_ind = model@state_var_indices)
    rownames(res) <- variables
    colnames(res) <- 1:sim_length

    sim_res <- gecon_simulation(sim = res,
                                shocks = model@shocks[shock_ind],
                                shocks_tex = model@shocks_tex[shock_ind],
                                variables = variables,
                                variables_tex = model@variables_tex[var_ind],
                                sim_name = "Random path simulation",
                                model_info = model@model_info,
                                r_object_name = deparse(substitute(model)))

    return (sim_res)
}


# ############################################################################
# The simulate_model function simulates the impact of shock paths specified
# by the user on the model's variables.
# ############################################################################
# Input
#   model - an object of gecon_model class
#   variables - the names or indices of variables whose paths are to be simulated
#   shocks - the names or indices of shocks corresponding to consecutive rows
#            of the shock_path matrix. If missing, the \code{rownames} of
#            the shock_path matrix are used.
#   shock_path - a matrix simulated paths of \code{shocks} in rows
#   sim_length - the length of simulation path, the default value is 40
# Output
#   an object of gecon_simulation class
# ############################################################################
simulate_model <- function(model, variables = NULL, shocks = NULL,
                           shock_path, sim_length = 40)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }
    if (!model@is_stochastic) stop("the model is not stochastic")
    if (!model@re_solved) stop("solve the 1st order perturbation first using the 'solve_pert' function")

    if (!length(model@state_var_indices)) {
        stop(paste0("the model does not have any state variables; shocks do not ",
                    "have impact on variables in any periods but t"))
    }

    if (length(dim(shock_path)) != 2) {
        stop("invalid shock_path argument")
    }

    if (is.null(variables)) {
        variables <- model@variables
        var_ind <- 1:length(variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
        variables <- model@variables[var_ind]
    }

    if (is.null(shocks)) {
        shocks <- rownames(shock_path)
    } else if (dim(shock_path)[1] != length(shocks)) {
        stop("the number of shocks is different than the number of rows in shock_path argument")
    }

    shock_ind <- list2ind(shocks, model@shocks, "shock")
    shock_path0 <- matrix(0, nrow = length(model@shocks), ncol = sim_length)
    shock_path0[shock_ind, 1:min(sim_length, dim(shock_path)[2])] <- shock_path
    shock_path <- shock_path0

    res <- simulate_path(P = model@solution$P,
                         Q = model@solution$Q,
                         R = model@solution$R,
                         S = model@solution$S,
                         sim_length = sim_length,
                         init_st = NULL,
                         shock_path = shock_path,
                         var_ind = var_ind,
                         state_var_ind = model@state_var_indices)
    rownames(res) <- variables
    colnames(res) <- 1:sim_length

    sim_res <- gecon_simulation(sim = res,
                                shocks = model@shocks[shock_ind],
                                shocks_tex = model@shocks_tex[shock_ind],
                                variables = variables,
                                variables_tex = model@variables_tex[var_ind],
                                sim_name = "Simulation with user-supplied shock path",
                                model_info = model@model_info,
                                r_object_name = deparse(substitute(model)))

    return (sim_res)
}


# ############################################################################
# The compute_irf function computes impulse response functions
# ############################################################################
# Input
#   model - an object of gecon_model class
#   variables - the names or indices of variables whose responses are to be simulated
#   shocks - the names or indices of shocks for which IRFs are to be computed,
#            if missing, IRFs are computed for all shocks with non-zero variance
#   sim_length - the length of simulation path, the default value is 40
#   cholesky - if FALSE, IRFs are computed based on uncorrelated shocks,
#              otherwise, irfs are based on correlated shocks,
#              the default value is FALSE
# Output
#   an object of gecon_simulation class
# ############################################################################
compute_irf <- function(model, variables = NULL, shocks = NULL,
                        sim_length = 40, cholesky = TRUE)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }
    if (!model@is_stochastic) stop("the model is not stochastic")
    if (!model@re_solved) stop("solve the 1st order perturbation first using the 'solve_pert' function")

    if (cholesky) {
        if (!model@shock_cov_mat_flag)
            warning(paste0("by default, identity matrix has been assumed as the shock covariance matrix, ",
                           "setting cholesky option to TRUE will have no effect on the results in this case, ",
                           "covariance matrix can be set through 'set_shock_cov_mat' ",
                           "or 'set_shock_distr_par' functions."))
    }

    if (is.null(variables)) {
        variables <- model@variables
        var_ind <- 1:length(variables)
    } else {
        var_ind <- list2ind(variables, model@variables, "variable")
        variables <- model@variables[var_ind]
    }

    if (is.null(shocks)) {
        shock_ind <- which(model@active_shocks)
        shocks <- model@shocks[shock_ind]
    } else {
        shock_ind <- list2ind(shocks, model@shocks, "shock")
        shock_ind_ia <- intersect(shock_ind, which(!model@active_shocks))
        if (length(shock_ind_ia)) {
            warning(paste0("the following shocks are not active (have zero variance) ",
                           "and will not be used for the IRFs computation: ",
                           list2str(model@shocks[shock_ind_ia])))
        }
        shock_ind <- setdiff(shock_ind, shock_ind_ia)
        shocks <- model@shocks[shock_ind]
    }
    n_shock <- length(shock_ind)
    if (!n_shock) {
        stop("no (active) shocks selected")
    }

    if (cholesky) {
        shock_ind_a <- which(model@active_shocks)
        A <- t(chol(model@shock_cov_mat[shock_ind_a, shock_ind_a]))
        shock_ind_map_a <- match(shock_ind, shock_ind_a)
    }

    res <- array(0, dim = c(length(variables), sim_length, n_shock),
                 dimnames = list(variables, 1:sim_length, shocks))

    for (s in 1:n_shock) {
        shock_path <- matrix(0, nrow = length(model@shocks), ncol = sim_length)
        if (cholesky) {
            one <- rep(0, length(shock_ind_a))
            one[shock_ind_map_a[s]] <- 1
            shock_path[shock_ind_a, 1] <- A %*% one
        } else {
            shock_path[shock_ind[s], 1] <- 1
        }

        res[, , s] <- simulate_path(P = model@solution$P,
                                    Q = model@solution$Q,
                                    R = model@solution$R,
                                    S = model@solution$S,
                                    sim_length = sim_length,
                                    init_st = NULL,
                                    shock_path = shock_path,
                                    var_ind = var_ind,
                                    state_var_ind = model@state_var_indices)
    }

    sim_res <- gecon_simulation(sim = res,
                                shocks = model@shocks[shock_ind],
                                shocks_tex = model@shocks_tex[shock_ind],
                                variables = variables,
                                variables_tex = model@variables_tex[var_ind],
                                sim_name = "Impulse response functions",
                                model_info = model@model_info,
                                r_object_name = deparse(substitute(model)))

    return (sim_res)
}
