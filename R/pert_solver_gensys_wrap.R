# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################

# ############################################################################
# The gensys_wrap function solves perturbation of the first
# order represented by four sparse matrices using gensys solver developed
# by Christopher A. Sims
# ############################################################################
# Input
#   matrices [class Matrix] - model matrices (corresponding y[t-1], y[t],
#                             E[t](y[t+1]) and eps[t])
#   is_stochastic - a logical value, TRUE if model is stochastic and FALSE,
#                                    if model is deterministic
#   tol - a numeric value, the tolerance level for treating column elements
#                          as zeros
# Output
#   A list containing the following fields:
#       eig - eigenvalues of the system
#       P - matrix of impact of state variables on state variables,
#       Q - matrix of impact of exog shocks on state variables
#               (only when the model is stochastic),
#       R - matrix of impact of state variables on non-state variables,
#       S - matrix of impact of exog shocks on non-state
#                variables (only when the model is stochastic),
#       state_var_indices - indices of state variables
#       exit_code - 0 - solution exists and is unique,
#                   1 - solution exists but is not unique (sunspots),
#                   2 - solution not found,
#       exit_info - information returned by the solver
#       norm_deter - the norm of residuals of the model's deterministic part
#                        (solution check),
#       norm_stoch - the norm of residuals of the model's stochastic part
#                        (solution check)
# ############################################################################
gensys_wrap <- function(matrices, is_stochastic, tol)
{
    # conversion to dense
    A <- as.matrix(matrices[[1]])
    B <- as.matrix(matrices[[2]])
    C <- as.matrix(matrices[[3]])
    D <- as.matrix(matrices[[4]])

    # variables and equations
    no_var <- dim(A)[2]
    no_eq <- dim(A)[1]
    exp_var_ind <- nonempty_col_indices(C, tol)
    no_exp_var <- length(exp_var_ind)

    if (is_stochastic) {
        no_eps <- dim(D)[2]
        shock_indices <- c(1: no_eps)
    }

    if (no_exp_var) {
        g0 <- matrix(nrow = no_eq + no_exp_var, ncol = no_var + no_exp_var)
        g0[1:no_eq, 1:no_var] <- B
        g0[1:no_eq, (no_var + 1):(no_var + no_exp_var)] <- C[, exp_var_ind]
        exp_1 <- matrix(0, nrow = no_exp_var, ncol = no_var)
        constr <- array(c(1:no_exp_var, exp_var_ind) , dim = c(no_exp_var, 2))
        exp_1[constr] <- -1
        g0[(no_eq + 1):(no_eq + no_exp_var), 1:no_var] <- exp_1
        g0[(no_eq + 1):(no_eq + no_exp_var), (no_var + 1):
            (no_var + no_exp_var)] <-
                        matrix(0, nrow = no_exp_var, ncol = no_exp_var)
    } else {
        g0 <- matrix(nrow = no_eq, ncol = no_var)
        g0[1:no_eq, 1:no_var] <- B
    }

    # g1
    if (no_exp_var) {
        g1 <- matrix(nrow = no_eq + no_exp_var, ncol = no_var + no_exp_var)
        g1[1:no_eq, 1:no_var] <- A
        g1[1:no_eq, (no_var + 1):(no_var + no_exp_var)] <-
                                    matrix(0, nrow = no_eq, ncol = no_exp_var)
        g1[(no_eq + 1):(no_eq + no_exp_var), 1:no_var] <-
                                    matrix(0, nrow = no_exp_var, ncol = no_var)
        g1[(no_eq + 1):(no_eq + no_exp_var),
                        (no_var + 1):(no_var + no_exp_var)] <- diag(no_exp_var)
    } else {
        g1 <- A
    }

    # Pi matrix
    if (no_exp_var) {
        PI <- matrix(nrow = no_eq + no_exp_var, ncol = no_exp_var)
        PI[1:no_eq,] <- matrix(0, nrow = no_eq, ncol = no_exp_var)
        PI[(no_eq + 1):(no_eq + no_exp_var),] <- diag(no_exp_var)
    } else PI <- matrix(0, nrow = no_eq, ncol = 1)
    # Psi matrix
    if (is_stochastic) {
        if (no_exp_var) {
            Psi <- matrix(0, nrow = no_eq + no_exp_var, ncol = no_eps)
            Psi[1:no_eq,] <- D
            Psi[(no_eq + 1):(no_eq + no_exp_var),] <-
                    matrix(0, nrow = no_exp_var, ncol = no_eps)
        } else Psi <- D
    } else {
        if (no_exp_var) {
            Psi <- matrix(0, nrow = no_eq + no_exp_var, ncol = 1)
        } else {
            Psi <- matrix(0, nrow = no_eq, ncol = 1)
        }
    }

    # ##########
    # SOLUTION #
    # ##########
    suppressWarnings(sol <- gensys(-g0, g1, c0 = matrix(0, dim(g0)[1], 1), Psi, PI))

    # Eigenvalues
    eigenval <- (sol$gev[, 2] / sol$gev[, 1])
    eig <- matrix(0, length(which(Mod(eigenval) > 0)), 3)
    colnames(eig) <- c("Mod", "Re", "Im")
    eig[, 1] <- Mod(eigenval[which(Mod(eigenval) > 0)])
    eig[, 2] <- Re(eigenval[which(Mod(eigenval) > 0)])
    eig[, 3] <- Im(eigenval[which(Mod(eigenval) > 0)])
    eig <- eig[order(eig[,1]), , drop = FALSE]

    # gensys exit codes
    if ((sol$eu[1] == 1) & (sol$eu[2] == 1)) {
        exit_code <- 0
        exit_info <- "Model solved. Solution exists and is unique."
    } else if ((sol$eu[1] == 1) & (sol$eu[2] == 0)) {
        exit_code <- 1
        exit_info <- paste("Solution exists but is not unique (sunspots).",
                           "Model contains loose endogenous variables.")
    } else if ((sol$eu[1] == 0) & (sol$eu[2] == 1)) {
        exit_code <- 2
        exit_info <- "Solution does not exist."
    } else {
        exit_code <- 2
        exit_info <- paste("Coincident zeros.",
                           "Indeterminacy and/or nonexistence.\n")
    }

    if (exit_code != 0) {
        return (list(eig = eig,
                     P = NULL,
                     Q = NULL,
                     R = NULL,
                     S = NULL,
                     state_var_indices = numeric(0),
                     exit_code = exit_code,
                     exit_info = exit_info,
                     norm_deter = NA,
                     norm_stoch = NA))
    }

    # submatrices containing solution
    G1 <- matrix(sol$G1[1:no_var, 1:no_var],
                            nrow = no_var, ncol = no_var)
    if (is_stochastic) {
        impact <- matrix(sol$impact[1:no_var, 1:no_eps],
                            nrow = no_var, ncol = no_eps)
    }

    # finding state variables and shocks
    state_var_indices <- which(abs(G1[max.col(t(abs(G1))) +
                                no_var * (0:(no_var - 1))]) >= tol)
    state_n <- length(state_var_indices)
    if (is_stochastic) {
        shock_n <- length(shock_indices)
    }

    # finding matrices PP & QQ
    PP <- round2zero(matrix(G1, nrow = no_var, ncol = no_var), tol)
    if (is_stochastic) {
        QQ <- round2zero(matrix(impact, nrow = no_var, ncol = shock_n), tol)
    }

    # P, Q, R, S
    P <- matrix(PP[state_var_indices, state_var_indices],
                                    nrow = state_n, ncol = state_n)
    if (is_stochastic) {
        Q <- matrix(QQ[state_var_indices, shock_indices],
                    nrow = state_n, ncol = shock_n)
    } else {
        Q <- NULL
    }
    R <- matrix(PP[-c(state_var_indices), state_var_indices],
                                nrow = no_var - state_n, ncol = state_n)
    if (is_stochastic) {
        if (length(state_var_indices)) {
            S <- matrix(QQ[-c(state_var_indices), shock_indices],
                                nrow = no_var - state_n, ncol = shock_n)
        } else {
            S <- matrix(QQ[, shock_indices],
                                nrow = no_var - state_n, ncol = shock_n)
        }
    } else {
        S <- NULL
    }

    A_prim <- matrix(A[, state_var_indices],
                                nrow = no_eq, ncol = state_n)
    R_prim <- matrix(PP[, state_var_indices],
                                nrow = no_var, ncol = state_n)

    if (is_stochastic) {
        S_prim <- matrix(QQ[, shock_indices],
                                nrow = no_var, ncol = shock_n)
    }

    # Check deterministic part
#     norm_deter <- norm(A + B %*% PP + C %*% PP %*% PP)
    norm_deter <- norm(A_prim + B %*% R_prim + C %*% R_prim %*% P)
    # Check stochastic part
    if (is_stochastic) {
#         norm_stoch <- norm(B %*% QQ + C %*% PP %*% QQ + D)
        norm_stoch <- norm(B %*% S_prim + C %*% R_prim %*% Q + D)
    } else {
        norm_stoch <- NA
    }

    # Output
    return (list(eig = eig,
                 P = P,
                 Q = Q,
                 R = R,
                 S = S,
                 state_var_indices = state_var_indices,
                 exit_code = exit_code,
                 exit_info = exit_info,
                 norm_deter = norm_deter,
                 norm_stoch = norm_stoch))
}



# ############################################################################
# The nonempty_col_indices function finds indices of empty columns
# ############################################################################
# Input
#   inp_mat - a matrix to be reduced
#   tol - numerical inrelevance
# Output
#   a vector of empty columns' indices
# ############################################################################
nonempty_col_indices <- function(inp_mat, tol)
{
    a <- dim(inp_mat)
    indices <- which(colSums(abs(inp_mat)) > tol)
    return (indices)
}

