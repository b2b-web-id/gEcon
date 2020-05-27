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
# Solving the first order perturbation
# ############################################################################

# ############################################################################
# The solve_pert function (log-)linearises model, removes zero rows, and
# solves the first order perturbation
# ############################################################################
# Input
#   model - an object of gecon_model class
#   loglin - a logical value. If TRUE, all variables are selected for log-linearisation
#   loglin_var - vector of variables that are to be log-linearised
#   not_loglin_var - vector of variables that are not to be
#                    log-linearised (overrides previous settings)
#   tol - tolerance level of a solution (1e-6 by default)
#   solver - first order perturbation solver,
#            the default solver is Christopher Sims solver gensys
#
# Output
#   an object of \code{gecon_model} class representing the model
# ############################################################################
solve_pert <- function(model,
                       loglin = TRUE, loglin_var = NULL, not_loglin_var = NULL,
                       tol = 1e-6, solver = "gensys")
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }

    if (!is.logical(loglin) || (length(loglin) != 1)) {
        stop("loglin argument should be a logical value")
    }

    if (!model@is_dynamic) stop("the model is not dynamic")
    if (!model@ss_solved) stop("find the steady state first using 'steady_state' function")

    abcd <- model@pert(model@variables_ss_val,
                       model@parameters_calibr_val,
                       model@parameters_free_val)
    abcd <- row_reduction(abcd)
    # log-linearisation of a model
    model@loglin_var <- rep(loglin, length(model@variables))
    if (!is.null(loglin_var)) {
        ind <- list2ind(loglin_var, model@variables, "variable", "loglin_var")
        model@loglin_var[ind] <- TRUE
    }
    if (!is.null(not_loglin_var)) {
        ind <- list2ind(not_loglin_var, model@variables, "variable", "not_loglin_var")
        model@loglin_var[ind] <- FALSE
    }
    loglin_var <- abs(model@variables_ss_val)
    ind <- which(loglin_var < 1e-4)
    indl <- which(model@loglin_var)
    ind0l <- intersect(ind, indl)
    if (length(ind0l)) {
        mes <- paste0("the following variables will not be log-linearised as ",
                      "their steady-state values are equal or close to zero: ",
                      list2str(model@variables[ind0l], "\""))
        warning(mes, call. = FALSE)
        model@loglin_var[ind] <- FALSE
    }
    loglin_var[!model@loglin_var] <- 1
    abcd <- lintolog(abcd, loglin_var)

    # solution
    solve_stat <- tryCatch(solver_output <- pert_solver(abcd, model@is_stochastic,
                                                        tol, solver),
                           error = function(e) e)
    if (inherits(solve_stat, "error")) {
        stop(solve_stat)
    }

    # clear slots
    model@re_solved <- FALSE
    model@eig_vals <- matrix(nrow = 0, ncol = 0)
    model@solution <- list(P = NULL, Q = NULL, R = NULL, S = NULL)
    model@state_var_indices <- numeric(0)
    model@solver_exit_info <- character(0)
    model@solution_resid <- list(NULL)
    model@corr_mat <- matrix(nrow = 0, ncol = 0)
    model@corr_computed <- FALSE
    model@autocorr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_corr_mat <- matrix(nrow = 0, ncol = 0)
    model@ref_var_idx <- 0L
    model@var_dec <- matrix(nrow = 0, ncol = 0)
    model@sdev <- matrix(nrow = 0, ncol = 0)

    # eigenvalues
    model@eig_vals <- solver_output$eig
    # solver exit info
    model@solver_exit_info <- solver_output$exit_info

    # not solved
    if ((solver_output$exit_code != 0)) {
        mes <- paste0("Model has NOT been SOLVED.\nSolver exit info:\n",
                      solver_output$exit_info)
        warning(mes, call. = FALSE)
        return (model)
    }

    # solved
    model@re_solved <- TRUE
    cat("Model has been SOLVED\n")

    # residuals
    model@solution_resid <- list(norm_deter = solver_output$norm_deter,
                                 norm_stoch = solver_output$norm_stoch)
    if (model@is_stochastic) {
        if ((solver_output$norm_deter > tol ||
             solver_output$norm_stoch > tol)) {
            mes <- paste0("The solution has been found ",
                          "but it does NOT satisfy perturbation equations with a given tolerance.\n",
                          "Norm of residuals in the deterministic part: ",
                          solver_output$norm_deter, "\n",
                          "Norm of residuals in the stochastic part: ",
                          solver_output$norm_stoch, "\n",
                          "This may be caused by numerical roundoff errors. ",
                          "Consider changing the parametrisation of the model ",
                          "or log-linearising the variables with large steady-state values.")
                warning(mes, call. = FALSE)
        }
    } else if ((solver_output$norm_deter > tol)) {
        mes <- paste0("The solution has been found ",
                      "but it does NOT satisfy perturbation equations with a given tolerance.\n",
                      "Norm of residuals: ",
                      solver_output$norm_deter, "\n",
                      "This may be caused by numerical roundoff errors. ",
                      "Consider changing the parametrisation of the model ",
                      "or log-linearising the variables with large steady-state values.")
        warning(mes, call. = FALSE)
    }

    # state variables, naming rows and columns of solution matrices
    model@state_var_indices <- solver_output$state_var_indices
    if (!length(solver_output$state_var_indices)) {
        if (model@is_stochastic) {
            mes <- "There are no state variables. Only S matrix of the solution is not empty."
        } else {
            mes <- "There are no state variables. Solution matrices are empty."
        }
        warning(mes, call. = FALSE)
    } else {
        # P
        colnames(solver_output$P) <-
            paste0(model@variables[solver_output$state_var_indices], "[-1]")
        rownames(solver_output$P) <-
            paste0(model@variables[solver_output$state_var_indices], "[]")

        # R
        colnames(solver_output$R) <-
            paste0(model@variables[solver_output$state_var_indices], "[-1]")
        if (length(model@variables[-c(solver_output$state_var_indices)])) {
            rownames(solver_output$R) <-
                paste0(model@variables[-c(solver_output$state_var_indices)], "[]")
        }
        # Q
        if (model@is_stochastic) {
            colnames(solver_output$Q) <- model@shocks
            rownames(solver_output$Q) <-
                model@variables[solver_output$state_var_indices]
        }
    }
    # S
    if (model@is_stochastic) {
        if (length(model@variables[-c(solver_output$state_var_indices)])) {
            rownames(solver_output$S) <-
                model@variables[-c(solver_output$state_var_indices)]
        }
        colnames(solver_output$S) <-  model@shocks
    }

    # save solution
    model@solution <- list(P = solver_output$P,
                           Q = solver_output$Q,
                           R = solver_output$R,
                           S = solver_output$S)

    return (model)
}


# ############################################################################
# The get_pert_solution function prints and returns the recursive laws of motion
# of the model's variables
# ############################################################################
# Input
#   model - an object of gecon_model class
#   to_tex - logical, if TRUE results are written to a .tex file (FALSE by default)
#   silent - logical, if TRUE, console output is suppressed (FALSE by default)
# Output
#   returns a list of P, Q, R, S matrices (P and R matrices
#           in case of deterministic model), a vector of variables' steady-state
#           values, a vector of flags indicating which variables
#           were log-linearised, and a vector of indices of state variables
# ############################################################################
get_pert_solution <- function(model, to_tex = FALSE, silent = FALSE)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }

    if (!model@is_dynamic) stop("model is not dynamic")
    if (!model@re_solved)
        stop(paste0("the model in its (log-)linearised form has not been solved ",
                    "(solve it first using \'solve_pert\' function)"))

    P <- model@solution$P
    Q <- model@solution$Q
    R <- model@solution$R
    S <- model@solution$S
    ss_val <- model@variables_ss_val
    names(ss_val) <- model@variables
    loglin <- model@loglin_var
    names(loglin) <- model@variables
    state_ind <- model@state_var_indices

    if (model@is_stochastic) {
        returns_stoch <- list(P = P, Q = Q, R = R, S = S,
                              ss_val = ss_val, loglin = loglin, state_ind = state_ind)

        if (!silent) {
            cat("\nMatrix P:\n\n")
            print(round(P, digits = 4))
            cat("\n\nMatrix Q:\n\n")
            print(round(Q, digits = 4))
            cat("\n\nMatrix R:\n\n")
            print(round(R, digits = 4))
            cat("\n\nMatrix S:\n\n")
            print(round(S, digits = 4))
        }

        if (to_tex) {
            texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
            tex <- "\\section{The solution of the 1st order perturbation}\n\n"
            # P
            colnames(P) <- paste0(model@variables_tex[model@state_var_indices], "_{t-1}")
            rownames(P) <- paste0(model@variables_tex[model@state_var_indices], "_{t}")
            tex <- paste0(tex, "\\subsection*{Matrix $P$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(P, digits = 4)))
            # Q
            colnames(Q) <- model@shocks_tex
            rownames(Q) <- model@variables_tex[model@state_var_indices]
            tex <- paste0(tex, "\\subsection*{Matrix $Q$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(Q, digits = 4)))
            # R
            colnames(R) <- paste0(model@variables_tex[model@state_var_indices], "_{t-1}")
            rownames(R) <- paste0(model@variables_tex[-model@state_var_indices], "_{t}")
            tex <- paste0(tex, "\\subsection*{Matrix $R$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(R, digits = 4)))
            # S
            colnames(S) <- model@shocks_tex
            rownames(S) <- model@variables_tex[-model@state_var_indices]
            tex <- paste0(tex, "\\subsection*{Matrix $S$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(S, digits = 4)))

            write_tex(texf, tex, !silent)
        }

        return(invisible(returns_stoch))
    } else {
        returns_deter <- list(P = P, R = R,
                              ss_val = ss_val, loglin = loglin, state_ind = state_ind)

        if (!silent) {
            cat("\nMatrix P:\n\n")
            print(round(P, digits = 4))
            cat("\n\nMatrix R:\n\n")
            print(round(R, digits = 4))
        }

        if (to_tex) {
            texf <- paste0(rm_gcn_ext(model@model_info[2]), ".results.tex")
            tex <- "\\section{The solution of the 1st order perturbation}\n\n"
            # P
            colnames(P) <- paste0(model@variables_tex[model@state_var_indices], "_{t-1}")
            rownames(P) <- paste0(model@variables_tex[model@state_var_indices], "_{t}")
            tex <- paste0(tex, "\\subsection*{Matrix $P$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(P, digits = 4)))
            # R
            colnames(R) <- paste0(model@variables_tex[model@state_var_indices], "_{t-1}")
            rownames(R) <- paste0(model@variables_tex[-model@state_var_indices], "_{t}")
            tex <- paste0(tex, "\\subsection*{Matrix $R$}\n\n")
            tex <- paste0(tex, mat2texbmat(round(R, digits = 4)))

            write_tex(texf, tex, !silent)
        }

        return(invisible(returns_deter))
    }
}


# ############################################################################
# The check_bk function checks the Blanchard-Kahn conditions
# ############################################################################
# Input
#   model - an object of gecon_model class
# Output
#   Checks the Blanchard-Kahn conditions and prints info about eigenvalues
#   larger than 1 in modulus and the number of forward looking variables
# ############################################################################
check_bk <- function(model)
{
    if (!is.gecon_model(model)) {
        stop("model argument should be of gecon_model class")
    }

    if (!length(model@eig_vals)) {
        stop("function \'solve_pert\' has to be called first")
    }

    # number of forward looking variables
    matrices <- model@pert(model@variables_ss_val,
                           model@parameters_calibr_val,
                           model@parameters_free_val)

    tol <- 1e-7
    matrices <- row_reduction(matrices)
    C <- as.matrix(matrices[[3]])
    exp_var_ind <- nonempty_col_indices(C, tol)
    fl <- unique(exp_var_ind)

    cat("\nEigenvalues of system:\n")
    print(model@eig_vals)

    cat(paste0("\nThe model has ",
               numnoun(length(fl), "forward looking variable"),
               " and ",
               numnoun(length(which(model@eig_vals[, 1] > 1)), "eigenvalue"),
               " larger than 1 in modulus.\n"))
    if (length(fl) ==  length(which(model@eig_vals[, 1] > 1)))
        cat("BK conditions have been SATISFIED.\n")
    else cat("BK conditions have NOT been SATISFIED.\n")
}


# ############################################################################
# The lintolog function converts 4 sparse matrices (class Matrix) of
# a linerized system of equations - corresponding to variable vectors:
# y(t-1), y(t), E(t)y(t+1) and eps(t+1) - into 4 sparse
# matrices (class Matrix) of a log-linearized system.
# When the steady-state value is equal or nearly equal to 0, equation
# is not log-linearized
# ############################################################################
# Input
#   matrices - a list of 4 matrices,
#   y_ss - a vector of steady-state values of variables (or 1's)
# Output
#   a list of 4 matrices
# ############################################################################
lintolog <- function(matrices, y_ss)
{
    no_var <- dim(matrices[[1]])[2]
    D_ss <- Diagonal(n = no_var, x = y_ss)
    A <- matrices[[1]] %*% D_ss
    B <- matrices[[2]] %*% D_ss
    C <- matrices[[3]] %*% D_ss
    D <- matrices[[4]]

    return (list(A, B, C, D))
}


# ############################################################################
# The to_remove function finds indices of empty rows
# ############################################################################
# Input
#   mat - a matrix to be reduced
#   tol - numerical inrelevance
# Output
#   vector of empty rows indices
# ############################################################################
to_remove <- function(mat, tol = 1e-6)
{
    mat <- as.matrix(mat)
    size <- dim(mat)[2]
    size_c <- dim(mat)[1]
    mat <- t(mat)
    to_be_removed <- which(abs(mat[max.col(abs(t(mat))) +
                        size * c(0: (size_c - 1))]) < tol)
    return (to_be_removed)
}



# ############################################################################
# The remo function removes specified indices from matrix
# ############################################################################
# Input
#   mat - a matrix to be reduced
#   to_rem - indices of rows to be removed
# Output
#   matrix with removed rows
# ############################################################################
remo <- function(mat, to_rem)
{
    mat <- Matrix(mat[-c(to_rem),])
}

# ############################################################################
# The row_reduction function removes empty rows from 3 first
#   elements of a list of matrices
# ############################################################################
# Input
#   can - a list of matrices
# Output
#   list of matrices with removed rows
# ############################################################################
row_reduction <- function(can)
{
    can2 <- list(can[[1]], can[[2]], can[[3]])
    rem <- Reduce(intersect, lapply(can2, to_remove))
    if (length(rem) > 0) {
        can <- lapply(can, remo, rem)
    }
    return (can)
}





