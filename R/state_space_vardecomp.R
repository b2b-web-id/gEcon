# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# The vardecomp function for theoretical variance decomposition
# ###################################################################
# Inputs
#   P - a matrix[n_var, n_state_var] impact of lagged state variables
#       on all the variables in model
#   cov_mat - a matrix[n_shocks, n_shocks] shock variance covariance
#                   matrix
#   Q - a matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   PP - a matrix[n_state_var, n_state_var] an impact of lagged 
#       state variables on state variables
#   QQ - a matrix[n_var, n_shocks] instantaneous an impact of shocks
#       on state variables in model
#   lambda - a numeric, lambda to HP filter time series
#   ngrid - a numeric, the grid for FFT computation
# Outputs
#   a matrix with variance decomposition
# ##################################################################
vardecomp <- function(P, cov_mat, Q, PP, QQ, 
                lambda = 1600, ngrid = 64 * 16)
{
    nexo <- ncol(Q);
    SS = cov_mat +  1e-14 * diag(nexo) ;
    cs = t(chol(SS));
    SS = cs %*% t(cs);

    freqs = seq(from = 0, by = ((2*pi)/ngrid),
                to = (2*pi*(1 - .5/ngrid)));
    tpos  = exp(1i * freqs);
    tneg  =  exp(-1i * freqs);
    
    IA = diag(dim(P)[2])
    IE = diag(nexo)
    math_col <- matrix(0, (ngrid), dim(PP)[1] ^ 2)
    hp1 = 4*lambda*(1 - cos(freqs))^2 / 
            (1 + 4 * lambda*(1 - cos(freqs))^2)
    
    # computing total variance
    for (ig in 1:ngrid) {
        # state variables
        f_omega  =(1/(2*pi)) * 
            ( rbind(solve(IA - P * tneg[ig], Q), IE) %*% SS %*% 
            cbind(t(Q) %*% solve(IA-t(P) * tpos[ig]), IE)); 
        # selected variables
        g_omega = cbind(PP * tneg[ig], QQ) %*%
                        f_omega %*% rbind(t(PP) *tpos[ig], t(QQ)); 
        f_hp = hp1[ig] ^ 2 * g_omega;
        math_col[ig,] <- t(matrix(f_hp))
    }
    
    imathp_col <- Re(mvfft(math_col, inverse = TRUE) / 
                    nrow(math_col)) * 2 * pi;
    vv <- diag(matrix(imathp_col[1, ], nrow(PP) , nrow(PP) ));
    var_decom <- matrix(0, nrow(PP) ,dim(Q)[2])
    
    
    # computing amount of variance for which each shock accounts
    for (i in 1:dim(Q)[2]) {       
        SSi = cs[ ,i] %*% t(cs[,i]);
        for (ig in 1:ngrid) {
            # state variables
            f_omega  =(1/(2*pi)) * 
                ( rbind(solve(IA - P * tneg[ig], Q), IE) %*% 
                 SSi %*% cbind(t(Q) %*% 
                 solve(IA-t(P) * tpos[ig]), IE));
            # selected variables
            g_omega = cbind(PP * tneg[ig], QQ) %*% f_omega %*%
                            rbind(t(PP) *tpos[ig], t(QQ)); 
            f_hp = hp1[ig] ^ 2 * g_omega;
            math_col[ig,] <- t(matrix(f_hp));
        }
        imathp_col <- Re(mvfft(math_col, inverse = TRUE) /
                            nrow(math_col)) * 2 * pi;
        # computing ratio of variance for each shock and all shocks
        var_decom[, i] = abs(diag(matrix(imathp_col[1, ],
                                nrow(PP) , nrow(PP) ))) / vv;
    }   
    return(var_decom)
}


