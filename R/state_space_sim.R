# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################
# Function simulate_path simulates path given shock realisations
# ############################################################################
# Inputs
#   P, Q, R, S - matrices describing the recursive laws of motion
#   sim_length - simulation horizon
#   init_st - initial deviations from the steady state
#   shock_path - shocks path
#   var_ind - indices of variables to be simulated
#   state_var_ind - indices of state variables
# Output
#   matrix with simulation results (variables in rows,
#   consecutive periods in columns)
# ############################################################################
simulate_path <- function(P, Q, R, S, sim_length = 40, init_st = NULL,
                          shock_path, var_ind, state_var_ind)
{
    if (any(nrow(P) != ncol(P), nrow(P) != nrow(Q),
            nrow(P) != ncol(R), nrow(R) != nrow(S),
            ncol(Q) != ncol(S))) {
        stop("inconsistent sizes of matrices")
    }

    # computing indices, ugly, but works...
    all <- 1:sum(dim(R))
    jump <- which(!(var_ind %in% state_var_ind))
    jump <- which(all[-state_var_ind] %in% var_ind[jump])
    jumpind <- (all[-state_var_ind])[jump]
    state <- which(var_ind %in% state_var_ind)
    state <- which(state_var_ind %in% var_ind[state])
    stateind <- state_var_ind[state]

    x <- matrix(data = NA, nrow = nrow(P), ncol = sim_length + 1)

    # if any of non-state variables is picked
    if (length(jump)) {
        y <- matrix(0, nrow = length(jump), ncol = sim_length + 1)

        # setting initial values
        if (!is.null(init_st)) {
            x[ , 1] <- init_st
        } else {
            x[ , 1] <- 0
        }

        RR <- matrix(R[jump,], nrow = length(jump), ncol = dim(R)[2])
        SS <- matrix(S[jump,], nrow = length(jump), ncol = dim(S)[2])

        # first step of simulation
        x[, 2] <- P %*% x[, 1] + Q %*% shock_path[, 1, drop = FALSE]
        y[, 2] <- RR %*% x[, 1] + SS %*% shock_path[, 1, drop = FALSE]

        # following steps
        if (sim_length > 1) {
            for (t_ind in 2:(sim_length)) {
                x[, t_ind + 1] <- P %*% x[, t_ind] + Q %*% shock_path[, t_ind, drop = FALSE]
                y[, t_ind + 1] <- RR %*% x[, t_ind] + SS %*% shock_path[, t_ind, drop = FALSE]
            }
        }

        if (length(state)) {
            res <- rbind(matrix(x[state, 2:(sim_length + 1)], ncol = sim_length),
                         matrix(y[, 2:(sim_length + 1)], ncol = sim_length))
            rownames(res) <- c(stateind, jumpind)
            roword <- match(var_ind, c(stateind, jumpind))
            res <- res[roword, , drop = FALSE]
            return (res)
        } else {
            res <- matrix(y[, 2:(sim_length + 1)], ncol = sim_length)
            rownames(res) <- jumpind
            roword <- match(var_ind, jumpind)
            res <- res[roword, , drop = FALSE]
            return (res)
        }

    } else { # only state variables are simulated
        if (!is.null(init_st) && length(init_st) == dim(P)[2]) {
            x[ , 1] <- init_st
        } else x[ , 1] <- 0

        x[, 2] <- P %*% x[, 1] + Q  %*% shock_path[, 1]
        for (t_ind in 2:sim_length) {
            x[, t_ind + 1] <- P %*% x[, t_ind] + Q %*% shock_path[, t_ind, drop = FALSE]
        }
        res <- matrix(x[state, 2:(t_ind + 1)], ncol = sim_length)
        rownames(res) <- stateind
        roword <- match(var_ind, stateind)
        res <- res[roword, , drop = FALSE]
        return (res)
    }
}

