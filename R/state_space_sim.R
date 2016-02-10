# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Function simulate_path simulates path given shock realisations
# ###################################################################
# Inputs
#   P, Q, R, S - matrices determining recursive laws of motion
#   tn - time horizon in which y_i is to be observed starting t=0
#   init_st - initial deviations from steady state
#   shock_path - shocks path
#   jump - indices non state variables to be simulated
#   stat - indices of state variables to be simulated
# Output
#   matrix with simulation results (variables in rows,
#   time in columns)
# ##################################################################
simulate_path <- function(P, Q, R, S, tn = 40, init_st = NULL,
                          shock_path, jump, stat) 
{
    # verifying consistency of the input data
    if (any(nrow(P) != ncol(P), nrow(P) != nrow(Q),
            nrow(P) != ncol(R), nrow(R) != nrow(S),
            ncol(Q) != ncol(S))) {
        stop("Inconsistent sizes of matrices")
    } else {
    
        # defining modes and sizes of the objects
        tn <- tn
        x <- matrix(data = NA, nrow= nrow(P), ncol = tn + 1)
        
        # if any of non state variables is picked to compute IRFs
        if (length(jump)) {
            y <-  matrix(0, nrow = length(jump), ncol = tn + 1)
            
            # setting initial values for computing IRFs
            if ( !is.null(init_st)  && length(init_st) == dim(P)[2]) {
                x[ , 1] <- init_st
            } else {
                x[ , 1] <- 0
            }
            
            RR <- matrix(R[jump,], nrow = length(jump), ncol = dim(R)[2])
            SS <- matrix(S[jump,], nrow = length(jump), ncol = dim(S)[2])
            
            # first step of simulation
            x[, 2] <- P %*% x[, 1] + Q %*% as.matrix(shock_path[, 1])
            y[, 2] <- RR %*% x[, 1] + SS %*% as.matrix(shock_path[, 1])
            
            # following steps
            if (tn > 1) {
                for (t_ind in 2:(tn)) { 
                    x[, t_ind + 1] <- P %*% x[, t_ind] + Q %*% as.matrix(shock_path[, t_ind])
                    y[, t_ind + 1] <- RR %*% x[, t_ind] + SS %*% as.matrix(shock_path[, t_ind])
                }
            }
            
            if (length(stat))
                return(rbind(matrix(x[stat, (2: (tn + 1))], ncol = tn),
                             matrix(y[ , (2: (tn + 1))], ncol = tn)))            
            else
                return(matrix(y[ , (2: (tn + 1))], ncol = (2: (tn + 1))))

        # simulates only state variables if such irfs are needed
        } else {
            # setting initial values for computing IRFs
            if ( !is.null(init_st)  && length(init_st) == dim(P)[2]) {
                x[ , 1] <- init_st
            } else x[ , 1] <- 0
            
            x[, 2] <- P %*% x[, 1] + Q  %*% shock_path[, 1]
            for (t_ind in 2:(tn)) { 
                x[, t_ind + 1] <- P %*% x[, t_ind] + Q %*% as.matrix(shock_path[, t_ind])
            }            
            return(x[stat, (2: (t_ind + 1))])
        }
    }
}
