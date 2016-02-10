# ###################################################################
# (c) Kancelaria Prezesa Rady Ministrów 2012-2015                   #
# Treść licencji w pliku 'LICENCE'                                  #
#                                                                   #
# (c) Chancellery of the Prime Minister 2012-2015                   #
# Licence terms can be found in the file 'LICENCE'                  #
#                                                                   #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                #
# ###################################################################
# Functions for computing theoretical and simulated covariance 
# matrix, correlations with chosen variable and autocovariances
# on either filtered or nonfiltered series
# ###################################################################

# ###################################################################
# The covar function computes covariance matrix 
# using spectral methods, correlations with chosen variable 
# and autocovariances both on filtered and nonfiltered series
# ###################################################################
# Inputs
#   P - a matrix[n_var, n_state_var] impact of lagged state variables
#       on all the variables in model
#   Q - a matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   cov_mat - a matrix[n_shocks, n_shocks] shock variance covariance
#               matrix
#   ngrid - grid for FFT computation
#   lambda - if series has to be filtered, sets lambda for HP filter
#   filter - a logical, TRUE means that series has to be HP-filtered
#   y_position - position of reference variable. For each model
#                variable the correlations with lead and lagged
#                values of the reference variable are computed.
#   n_leadlags - number of leads and lags when analysing correlation,
#                default = 6
# Output:
#   list with the following elements:
#       covmat  - a matrix of covariances 
#       autcor_output - a matrix of correlations with given variable
#       acov -  a matrix of autocovariances
# ###################################################################
covar <- function(P, Q, cov_mat, ngrid = 64 * 16, lambda = 1600,
                  filter = TRUE, y_position = 1, n_leadlags = 6)
{
    HP_LAMBDA <- lambda
    k <- ncol(P)
    N_GRID <- ngrid  # gridpoints for inverse Fourier transform. 
    MOM_TOL <-  0.0001

    # computing correlation matrix
    freq =seq(0, (2 * pi * (1 - 0.5 / N_GRID)), ((2 * pi) / N_GRID))
    hp = 4 * HP_LAMBDA * (1 - cos(freq))^2 /
                    (1 + 4 * HP_LAMBDA * (1 - cos(freq))^2)
    svv <- matrix(0, N_GRID, (ncol(P) *ncol(P)))
    P_may_be_singular <- (min(abs(eigen(P)$values)) < MOM_TOL)

    QSigtQ <- Q %*% cov_mat %*% t(Q)        

    for (pts in 1:N_GRID) {
        z <- exp(1i * freq[pts])
        zi <- exp(-1i * freq[pts])
        X = solve(diag(k) - P * zi);
        XH = solve(diag(k) - t(P) * z);
        s <-   1/ (2 * pi) * X %*% QSigtQ %*% XH       
        ssd <- s
        
        if (filter) {
            ssd <- hp[pts]^2 * ssd;
        }
        svv[pts, ] <- t(matrix(ssd))
    }
    
    autcov <- Re(mvfft(svv, inverse = T) / nrow(svv)) * (2 * pi)
    covmat <- autcov[1,]
    covmat <- matrix(covmat)
    n_select <- dim(P)[1]
    GNP_INDEX <- y_position; 
    N_LEADS_LAGS <- n_leadlags;
    
    # #######################################################
    # correlations with y or any other chosen variable
    # #######################################################
    
    # selecting indices of gdp variable from vectorization
    gnp_select = n_select*(c(0:(n_select-1))) + GNP_INDEX;
    # diagonal of vectorization
    diag_select = (n_select + 1)*c(0:(n_select-1)) + 1;

    # variance of gdp
    cov00 = autcov[1, ((GNP_INDEX-1)*n_select+GNP_INDEX)];  
    cov0 = (sqrt(cov00) * matrix(1, N_LEADS_LAGS, 1) ) %*% 
                sqrt(autcov[1, diag_select]);
                
    # standarizing
    autcor1 = t(autcov[1:N_LEADS_LAGS, gnp_select] / cov0)
    gnp_select <- seq(from = ((GNP_INDEX - 1) * n_select + 1) ,
                        to = (GNP_INDEX*n_select), by = 1)
    diag_select <- (n_select+1)*(0:(n_select - 1)) + 1; 
    autcor2 = t(autcov[1:N_LEADS_LAGS, gnp_select] / cov0);
    leadlags = seq(from = (1 - N_LEADS_LAGS),
                                to = (N_LEADS_LAGS - 1), by = 1)
    autcor_output = matrix(0, n_select, (N_LEADS_LAGS * 2 - 1))
    
    autcor_output[, (1: N_LEADS_LAGS)] <- autcor1[, ncol(autcor1):1] 
    autcor_output[, ((N_LEADS_LAGS + 1):(N_LEADS_LAGS * 2 - 1))] <-
                                autcor2[, 2:ncol(autcor2)]
    
    # #######################################################
    # Computing autocovariance matrix
    # #######################################################
      
    # autocovariance matrices 
    variances <- sqrt(diag(matrix(covmat, n_select, n_select)))
    acov <- matrix(0, n_select, (N_LEADS_LAGS - 1))
    for (i  in 2: N_LEADS_LAGS) {
        acov[, (i - 1) ] <- diag((matrix(autcov[i, ], n_select, n_select)
                                    / variances %*% t(variances)))
    }
  
    # return list
    cov_ret <- list(covmat = covmat, autcor_output = autcor_output, 
                        acov = acov)
    return(cov_ret)
}
               

# ###################################################################
# The covar_sim function computes simulated covariance matrix
# ###################################################################
# Inputs
#   P - matrix[n_var, n_state_var] impact of lagged state variables
#       on all the variables in model
#   Q - matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   cov_mat - matrix[n_shocks, n_shocks] shock variance covariance
#                   matrix
#   filter - logical, TRUE means that series has to be HP-filtered
#   lambda - if series has to be filtered, sets lambda for HP filter
#   nrun - number of MC simulations of series
#   y_position - position of the reference variable
#   n_leadlags - number of leads and lags when analyzing correlation,
#                default = 6
# Output:
#   covariance matrix of analyzed variables
# ###################################################################
covar_sim <- function(P, Q, cov_mat, filter = FALSE, 
                      lambda = 1600, nrun = 1000,
                      y_position = 1, n_leadlags = 6)
{
    # Simulation of x
    x <- matrix(0, dim(P)[1], nrun)
    for (i in 2: nrun) {
        eps <- chol(cov_mat) %*% matrix(rnorm(n = ncol(cov_mat),
                mean = 0, sd = 1), nrow = ncol(cov_mat), ncol = 1)
        x[, i] <- P %*% x[, i - 1]+ Q  %*% eps
    }
    
    # Filtering
    if (!filter) {
        a <- t(x)
    } else {
        a <- t(x)
        
        sparse_hp_wrap <- function(y, lambda)
        {
            return(sparse_hp(y, lambda)$cycle)
        }
        # Filter variables
        a <- apply(a, MARGIN = 2, FUN = sparse_hp_wrap, lambda)
    }    
    
    # compute autocovariances
    acf_wrapper <- function(t1) {
        return(acf(t1, plot = FALSE)$acf)
    }
    
    acov <- apply(a, MARGIN = 2, FUN = acf_wrapper)
    acov <- t(acov[2:(n_leadlags),])
    
    # compute covariance matrix
    covmat <- var(a)
    
    # compute covariances with specified variable
    a_small <- a[n_leadlags:(nrow(a) - n_leadlags), ]
    obs_numb <- nrow(a_small)
    y_a <- matrix(0, nrow = obs_numb, ncol = 2*(n_leadlags - 1) + 1)
    for (i in (2 * (n_leadlags - 1) + 1):1) {
        y_a[, i] <- a[i:(i + obs_numb - 1), y_position]
    }
    scale <- apply(a_small, MARGIN = 2, sd) %*% t(sd(a[, y_position])) %*% 
             matrix(rep(1, 2 * (n_leadlags - 1) + 1), 
                    nrow = 1, ncol = 2 * (n_leadlags - 1) + 1)
    corr_variable_mat <- (var(a_small, y_a) / scale)
    corr_variable_mat <- corr_variable_mat[, ncol(corr_variable_mat):1] 
           
    return(list(covmat = covmat,
                acov = acov,
                corr_variable_mat = corr_variable_mat))
}


# ###################################################################
# Function sparse_hp HP-filters input time series (uses sparse 
# matrices to create problem matrix and solve it).
# ###################################################################
# Inputs:
#   y - series to be filtered (numeric vector)
#   lambda - smoothing parameter lambda, default value  1600
# Outputs:
#   list with two series: trend and cyclical part
# ###################################################################
sparse_hp <- function(y, lambda = 1600) {
    
    # converting into numeric
    y <- as.numeric(y)
    Tspan <- length(y)
    
    # creating sparse matrix A
    diag1 <- lambda * c(1, 5, array(6, Tspan - 4), 5, 1) + array(1, Tspan)
    diag2 <- lambda * c(-2, array(-4, Tspan - 3), -2)
    diag3 <- lambda * array(1, Tspan - 2)
    A <- bandSparse(Tspan, k = -2:2, diagonals = list(diag3, diag2, diag1, diag2, diag3))
    
    # extracting trend and cycle
    trend <- solve(A, y)
    trend <- as.numeric(trend)
    cycle <- y - trend
    
    return(list(trend = trend,
                cycle = cycle))
}
