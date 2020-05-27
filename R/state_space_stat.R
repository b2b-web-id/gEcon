# ############################################################################
# This file is a part of gEcon.                                              #
#                                                                            #
# (c) Chancellery of the Prime Minister of the Republic of Poland 2012-2015  #
# (c) Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak 2015-2018    #
# License terms can be found in the file 'LICENCE'                           #
#                                                                            #
# Authors: Karol Podemski, Kaja Retkiewicz-Wijtiwiak                         #
# ############################################################################
# Functions for computing variables' covariances, autocorrelations,
# correlations with reference variable and variance decomposition
# on either filtered or unfiltered series.
# ############################################################################

# ############################################################################
# The covar_fft function computes model statistics using spectral methods.
# ############################################################################
# Inputs
#   P - matrix[n_var, n_var] impact of lagged variables on current values
#   Q - matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   Sigma - shock covariance matrix
#   n_leadlags - number of leads and lags when analysing correlation,
#                default 6
#   n_leadlags - the number of leads/lags for computing correlations.
#   ref_var_idx - the index of the variable with which correlations are computed.
#   lambda - HP filter lambda, if 0 no filtering is performed, default is 1600 (quarterly data).
#   ngrid - numeric value, density of grid used by the Fast Fourier transform (used only if \code{sim}
#                          option is set to FALSE). It has to be a multiple of 8 and has
#                          to be sufficiently large to obtain unbiased results.
# Output:
#   list with the following elements:
#       covmat  - a matrix of covariances
#       acov -  a matrix of autocovariances
#       ref_var_cor - a matrix of correlations with given variable
# ############################################################################
covar_fft <- function(P, Q, Sigma, n_leadlags, ref_var_idx, lambda, ngrid)
{
    # MOM_TOL <- 0.0001
    # P_may_be_singular <- (min(abs(eigen(P)$values)) < MOM_TOL)
    freq <- seq(0, (2 * pi * (1 - 0.5 / ngrid)), ((2 * pi) / ngrid))
    if (lambda) {
        hp <- 4 * lambda * (1 - cos(freq))^2 / (1 + 4 * lambda * (1 - cos(freq))^2)
    } else {
        hp <- rep(1, length(freq))
    }
    svv <- matrix(0, ngrid, (ncol(P) * ncol(P)))

    QSigtQ <- Q %*% Sigma %*% t(Q)

    k <- ncol(P)
    for (pts in 1:ngrid) {
        z <- exp(1i * freq[pts])
        zi <- exp(-1i * freq[pts])
        X <- solve(diag(k) - P * zi)
        XH <- solve(diag(k) - t(P) * z)
        s <- 1 / (2 * pi) * X %*% QSigtQ %*% XH
        ssd <- s
        ssd <- hp[pts]^2 * ssd
        svv[pts, ] <- t(matrix(ssd))
    }

    autcov <- Re(mvfft(svv, inverse = T) / nrow(svv)) * (2 * pi)
    autcov <- round2zero(autcov, 1e-6)
    covmat <- autcov[1,]
    covmat <- matrix(covmat, k, k)
    n_select <- dim(P)[1]

    # correlations with the reference variable
    if (ref_var_idx) {
        ref_var_select <- n_select*(c(0:(n_select-1))) + ref_var_idx
        diag_select <- (n_select + 1) * c(0:(n_select-1)) + 1

        cov00 <- autcov[1, ((ref_var_idx - 1) * n_select + ref_var_idx)]
        cov0 <- (sqrt(cov00) * matrix(1, n_leadlags + 1, 1) ) %*%
                    sqrt(autcov[1, diag_select])

        autcor1 <- t(autcov[1:(n_leadlags + 1), ref_var_select] / cov0)
        ref_var_select <- seq(from = ((ref_var_idx - 1) * n_select + 1),
                              to = (ref_var_idx * n_select), by = 1)
        autcor2 <- t(autcov[1:(n_leadlags + 1), ref_var_select] / cov0)
        leadlags <- seq(from = -n_leadlags, to = n_leadlags, by = 1)
        ref_var_cor <- matrix(0, n_select, (n_leadlags * 2 + 1))
        ref_var_cor[, 1:(n_leadlags + 1)] <- autcor1[, ncol(autcor1):1]
        ref_var_cor[, ((n_leadlags + 2):(n_leadlags * 2 + 1))] <-
                                    autcor2[, 2:ncol(autcor2)]
    } else {
        ref_var_cor <- NULL
    }

    # autocovariance matrices
    variances <- sqrt(diag(matrix(covmat, n_select, n_select)))
    acor <- matrix(0, n_select, n_leadlags)
    for (i in 2:(n_leadlags + 1)) {
        acor[, (i - 1) ] <- diag((matrix(autcov[i, ], n_select, n_select)
                                    / variances %*% t(variances)))
    }

    # return list
    cov_ret <- list(covmat = covmat,
                    acor = acor,
                    ref_var_cor = ref_var_cor)
    return (cov_ret)
}

# ############################################################################
# The vardecomp_fft determines variance decomposition using spectral methods.
# ############################################################################
# Inputs
#   P - matrix[n_var, n_state_var] impact of lagged state variables
#       on all the variables in model
#   Sigma - a matrix[n_shocks, n_shocks] shock covariance
#                   matrix
#   Q - matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   PP - matrix[n_state_var, n_state_var] an impact of lagged
#       state variables on state variables
#   QQ - matrix[n_var, n_shocks] instantaneous an impact of shocks
#       on state variables in model
#   lambda - HP filter lambda, if 0 no filtering is performed, default is 1600 (quarterly data).
#   ngrid - numeric value, density of grid used by the Fast Fourier transform (used only if \code{sim}
#                          option is set to FALSE). It has to be a multiple of 8 and has
#                          to be sufficiently large to obtain unbiased results.
# Outputs
#   a matrix with variance decomposition
# ############################################################################
vardecomp_fft <- function(P, Sigma, Q, PP, QQ, lambda, ngrid)
{
    nexo <- ncol(Q)
    SS <- Sigma + 1e-14 * diag(nexo)
    cs <- t(chol(SS))
    SS <- cs %*% t(cs)

    freqs <- seq(from = 0, by = ((2 * pi) / ngrid), to = (2 * pi * (1 - .5 / ngrid)))
    tpos <- exp(1i * freqs)
    tneg <- exp(-1i * freqs)

    IA <- diag(dim(P)[2])
    IE <- diag(nexo)
    math_col <- matrix(0, (ngrid), dim(PP)[1] ^ 2)
    if (lambda) {
        hp1 <- 4 * lambda * (1 - cos(freqs))^2 / (1 + 4 * lambda * (1 - cos(freqs))^2)
    } else {
        hp1 <- rep(1, length(freqs))
    }

    # computing total variance
    for (ig in 1:ngrid) {
        # state variables
        f_omega <- (1/(2*pi)) *
                   (rbind(solve(IA - P * tneg[ig], Q), IE) %*% SS %*%
                          cbind(t(Q) %*% solve(IA - t(P) * tpos[ig]), IE))
        # selected variables
        g_omega <- cbind(PP * tneg[ig], QQ) %*%
                         f_omega %*% rbind(t(PP) * tpos[ig], t(QQ))
        f_hp <- hp1[ig]^2 * g_omega
        math_col[ig,] <- t(matrix(f_hp))
    }

    imathp_col <- Re(mvfft(math_col, inverse = TRUE) /
                     nrow(math_col)) * 2 * pi
    vv <- diag(matrix(imathp_col[1,], nrow(PP), nrow(PP)))
    var_decom <- matrix(0, nrow(PP) ,dim(Q)[2])

    # computing amount of variance for which each shock accounts
    for (i in 1:dim(Q)[2]) {
        SSi <- cs[,i] %*% t(cs[,i])
        for (ig in 1:ngrid) {
            # state variables
            f_omega <- (1 / (2 * pi)) *
                       (rbind(solve(IA - P * tneg[ig], Q), IE)
                        %*% SSi %*% cbind(t(Q) %*% solve(IA-t(P) * tpos[ig]), IE))
            # selected variables
            g_omega <- cbind(PP * tneg[ig], QQ) %*% f_omega %*%
                             rbind(t(PP) * tpos[ig], t(QQ))
            f_hp <- hp1[ig] ^ 2 * g_omega
            math_col[ig,] <- t(matrix(f_hp))
        }
        imathp_col <- Re(mvfft(math_col, inverse = TRUE) /
                         nrow(math_col)) * 2 * pi
        # computing ratio of variance for each shock and all shocks
        var_decom[, i] <- abs(diag(matrix(imathp_col[1,],
                              nrow(PP), nrow(PP)))) / vv
    }
    return (var_decom)
}


# ############################################################################
# The covar_sim function computes model statistics through model simulation.
# ############################################################################
# Inputs
#   P - matrix[n_var, n_var] impact of lagged variables on current values
#   Q - matrix[n_var, n_shocks] instantaneous impact of shocks
#       on all the variables in model
#   Sigma - shock covariance matrix
#   n_leadlags - the number of leads/lags for computing correlations.
#   ref_var_idx - the index of the variable with which correlations are computed.
#   lambda - HP filter lambda, if 0 no filtering is performed, default is 1600 (quarterly data).
#   sim_length - length of simulation
# Output:
#   list with the following elements:
#       covmat  - a matrix of covariances
#       acov -  a matrix of autocovariances
#       ref_var_cor - a matrix of correlations with given variable
# ############################################################################
covar_sim <- function(P, Q, Sigma, n_leadlags, ref_var_idx, lambda, sim_length)
{
    # simulation
    x <- matrix(0, dim(P)[1], sim_length)
    for (i in 2: sim_length) {
        eps <- chol(Sigma) %*% matrix(rnorm(n = ncol(Sigma),
                mean = 0, sd = 1), nrow = ncol(Sigma), ncol = 1)
        x[, i] <- P %*% x[, i - 1] + Q %*% eps
    }

    a <- t(x)
    sparse_hp_wrap <- function(y, lambda) return (sparse_hp(y, lambda)$cycle)
    a <- apply(a, MARGIN = 2, FUN = sparse_hp_wrap, lambda)

    # compute autocovariances
    acf_wrapper <- function(t1) {
        return (acf(t1, plot = FALSE)$acf)
    }

    acor <- apply(a, MARGIN = 2, FUN = acf_wrapper)
    acor <- t(acor[2:(n_leadlags + 1),])

    # compute covariance matrix
    covmat <- round2zero(var(a), 1e-6)

    # compute covariances with specified variable
    if (ref_var_idx) {
        a_small <- a[(n_leadlags + 1):(nrow(a) - n_leadlags - 1), ]
        obs_numb <- nrow(a_small)
        y_a <- matrix(0, nrow = obs_numb, ncol = 2 * n_leadlags + 1)
        for (i in (2 * n_leadlags + 1):1) {
            y_a[, i] <- a[i:(i + obs_numb - 1), ref_var_idx]
        }
        scale <- apply(a_small, MARGIN = 2, sd) %*% t(sd(a[, ref_var_idx])) %*%
                matrix(rep(1, 2 * n_leadlags + 1),
                        nrow = 1, ncol = 2 * n_leadlags + 1)
        ref_var_cor <- (var(a_small, y_a) / scale)
        ref_var_cor <- ref_var_cor[, ncol(ref_var_cor):1]
    } else {
        ref_var_cor <- NULL
    }

    return (list(covmat = covmat,
                 acor = acor,
                 ref_var_cor = ref_var_cor))
}


# ############################################################################
# Function sparse_hp HP-filters input time series (uses sparse
# matrices to create problem matrix and solve it).
# ############################################################################
# Inputs:
#   y - series to be filtered (numeric vector)
#   lambda - smoothing parameter lambda, default value 1600
# Outputs:
#   list with two series: trend and cyclical part
# ############################################################################
sparse_hp <- function(y, lambda = 1600)
{
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

    return (list(trend = trend, cycle = cycle))
}
