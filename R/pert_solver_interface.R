# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# The sims_solver() function solves perturbation of the first
#     order represented by four sparse matrices 
# ###################################################################
# Input
#   matrices [class Matrix] - corresponding to variable vectors 
#                           y(t-1), y(t), E(t)y(t+1) and eps(t) -
#                           by means of a solver developed by 
#                           Christopher A. Sims
#   is_stochastic - a logical value denoting true is 
#                   model is stochastic and false if model is 
#                   deterministic
#   tol - a numeric, the tolerance for treating column elements as zeros
# Output
#   List with following fields:
#       norm_deter - a norm of deterministic part of model 
#                           (solution check),
#       norm_stoch - a norm of stochastic part of model 
#                               (solution check),
#       norm_deter_prim - a norm of deterministic part of model 
#                  computed in a bit different way 
#                               (solution check),
#       norm_stoch_prim = a norm of stochastic part of model 
#                  computed in a bit different way (solution check),
#       P - a matrix of impact of state variables on state variables,
#       Q - a matrix of impact of exog shocks variables on state variables
#               (only when the model is stochastic),
#       R - a matrix of impact of state variables on non-state variables,
#       S - a matrix of impact of exog shocks variables on non-state 
#                variables (only when the model is stochastic),
#       slv_exit_code = an exit code of Sims solver
#       state_var_indices = indices of state variables
#       eig - eigenvalues of the system
# ###################################################################
sims_solver <- function(matrices, is_stochastic, tol = 1e-8)
{
    # conversion to dense
    A <- as.matrix(matrices[[1]])
    B <- as.matrix(matrices[[2]])
    C <- as.matrix(matrices[[3]])
    D <- as.matrix(matrices[[4]])
    
    # variables and equations
    no_var <- dim(A)[2]
    no_eq <- dim(A)[1]
    exp_var_ind <- nonempty_col_indices(C)
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
            Psi = matrix(0, nrow = no_eq + no_exp_var, ncol = 1)
        } else {
            Psi = matrix(0, nrow = no_eq, ncol = 1)
        }
    }

    # ##########
    # SOLUTION #
    # ##########
    Sims <- gensys(-g0, g1, c0 = matrix(0, dim(g0)[1], 1), Psi, PI)
   
    # Eigenvalues
    eigenval <- (Sims$gev[, 2] / Sims$gev[, 1])
    eig <- matrix(0, length(which(Mod(eigenval) > 0)), 3)
    colnames(eig) <- c('Mod', 'Re', 'Im')
    eig[, 1] <- Mod(eigenval[which(Mod(eigenval) > 0)])
    eig[, 2] <- Re(eigenval[which(Mod(eigenval) > 0)])
    eig[, 3] <- Im(eigenval[which(Mod(eigenval) > 0)])
    eig <- eig[order(eig[,1]), ]
        
    if (Sims$eu[1] == -2 & Sims$eu[2] == -2) {
        #stop(paste('Model was not solved. Exit codes are \n Existence:',
        #           Sims$eu[1], '\n Unique:', Sims$eu[2]))
        if (is_stochastic) {
            output <- list(norm_deter = Inf,
                       norm_stoch = Inf,
                       norm_deter_prim = Inf,
                       norm_stoch_prim = Inf,
                       P = matrix(0),
                       Q = matrix(0),
                       R = matrix(0),
                       S = matrix(0),
                       slv_exit_info = paste("Coincident zeros.",
                                        "Indeterminacy and/or nonexistence.\n"),
                       slv_exit_code = Sims$eu,
                       state_var_indices = numeric(0),
                       eig = eig)
        } else {
            output <- list(norm_deter = Inf,
                           norm_deter_prim = Inf,
                           P = matrix(0),
                           R = matrix(0),
                           slv_exit_info = paste("Coincident zeros.",
                                                 "Indeterminacy and/or nonexistence.\n"),
                           slv_exit_code = Sims$eu,
                           state_var_indices = numeric(0),
                           eig = eig)            
        }
        return(output)
    } else {        
        # submatrices containing solution
        G1 <- matrix(Sims$G1[1:no_var, 1:no_var], 
                                nrow = no_var, ncol = no_var)
        if (is_stochastic) {
            impact <- matrix(Sims$impact[1:no_var, 1:no_eps], 
                                nrow = no_var, ncol = no_eps) 
        }
 
        # finding state variables and shocks
        state_var_indices <- which(abs(G1[max.col(t(abs(G1))) +
                                    no_var * (0:(no_var - 1))]) > tol)        
        state_n <- length(state_var_indices)
        if (is_stochastic) {
            shock_n <- length(shock_indices)
        }
        
        # finding matrices PP & QQ
        PP <- matrix(G1, nrow = no_var, ncol = no_var)
        
        if (is_stochastic) {
            QQ <- matrix(impact, nrow = no_var, ncol = shock_n)
        }
        
        # Check deterministic part
        n1 <- norm(A + B %*% PP + C %*% PP %*% PP)
        
        # Check stochastic part
        if (is_stochastic) {
            n2 <- norm(B %*% QQ + C %*% PP %*% QQ + D)
        }

        # P, Q, R, S
        P = matrix(PP[state_var_indices, state_var_indices], 
                                        nrow = state_n, ncol = state_n)
        
        if (is_stochastic) {
            Q = matrix(QQ[state_var_indices, shock_indices], 
                        nrow = state_n, ncol = shock_n)
        }
        
        R = matrix(PP[-c(state_var_indices), state_var_indices], 
                                 nrow = no_var - state_n, ncol = state_n)
        if (is_stochastic) {
            if (length(state_var_indices)) {
                S = matrix(QQ[-c(state_var_indices), shock_indices], 
                                    nrow = no_var - state_n, ncol = shock_n)
            } else {
                S = matrix(QQ[, shock_indices], 
                                    nrow = no_var - state_n, ncol = shock_n)            
            }
        }
        

        A_prim = matrix(A[, state_var_indices], 
                                    nrow = no_eq, ncol = state_n)
        R_prim = matrix(PP[, state_var_indices], 
                                    nrow = no_var, ncol = state_n)
        
        if (is_stochastic) {
            S_prim = matrix(QQ[, shock_indices], 
                                    nrow = no_var, ncol = shock_n)
        }
        # Check deterministic part
        n3 <- norm(A_prim + B %*% R_prim + C %*% R_prim %*% P)
        
        # Check stochastic part
        if (is_stochastic) {        
            n4 <- norm(B %*% S_prim + C %*% R_prim %*% Q + D)
        }        

        
        # Decoding Sims solver error messages
        if ((Sims$eu[1] == 1) & (Sims$eu[2] == 1)) {
            slv_exit_info = "Model solved. Solution exists and is unique."
        } else if ((Sims$eu[1] == 1) & (Sims$eu[2] == 0)) {
            slv_exit_info = paste("Solution exists but is not unique (sunspot).",
                                   "Model contains loose endogenous variables.")
        } else if ((Sims$eu[1] == 0) & (Sims$eu[2] == 1)) {
            slv_exit_info = "Solution doesn't exist."            
        } else {
            slv_exit_info = paste("Coincident zeros.",
                                "Indeterminacy and/or nonexistence.\n")
        }
        
        # Output   
        if (is_stochastic) {
        output <- list(norm_deter = n1,
                       norm_stoch = n2,
                       norm_deter_prim = n3,
                       norm_stoch_prim = n4,
                       P = P,
                       Q = Q,
                       R = R,
                       S = S,
                       slv_exit_info = slv_exit_info,
                       slv_exit_code = Sims$eu,
                       state_var_indices = state_var_indices,
                       eig = eig)
        } else {
            output <- list(norm_deter = n1,
                           norm_deter_prim = n3,
                           P = P,
                           R = R,
                           slv_exit_info = slv_exit_info,
                           slv_exit_code = Sims$eu,
                           state_var_indices = state_var_indices,
                           eig = eig)            
        }
        
        return(output)   
    }  
}

# ###################################################################
# The lintolog function converts 4 sparse matrices (class Matrix) of 
# a linerized system of equations - corresponding to variable vectors: 
# y(t-1), y(t), E(t)y(t+1) and eps(t+1) - into 4 sparse 
# matrices (class Matrix) of a log-linearized system.  
# When the steady state value is equal or nearly equal 0, equation 
# is not loglinearized
# ###################################################################
# Input
#   matrices - a list of 4 matrices, 
#   y_ss - a vector of steady-state values 
#               of variables (variables' order in 
#               a variable vector and a vector containing 
#               ss values should be the same)
#   tol - tolerance to treat given steady state value as zero
# Output
#   a list of 4 matrices
# ###################################################################
lintolog <- function(matrices, y_ss, tol = 1e-6)
{
    no_var <- dim(matrices[[1]])[2]
    
    ind <- which(abs(y_ss) < tol)
    
    if (length(ind) != 0) y_ss[ind] <- 1
 
    D_ss <- Diagonal(n = no_var, x = y_ss)

    A <- matrices[[1]] %*% D_ss
    B <- matrices[[2]] %*% D_ss
    C <- matrices[[3]] %*% D_ss
    D <- matrices[[4]]
    
    return(list(A, B, C, D))
}
