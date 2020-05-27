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
# QZ decomposition of complex matrices A, B using the Lapack routine zgges
# ############################################################################
# INPUTS:
#        A, B - complex matrices (size, size)
# OUTPUTS
#        a, b, q, z - complex matrices of dimensions (size, size)
#                         such that:
#                             A = q %*% a %*% t(ZZ)
#                             B = q %*% b %*% t(ZZ)
#                         holds for real QQ, ZZ
#
#                         For complex matrices q and z t,
#                         t(z) should be interpreted as
#                         hermitian transpose:
#                             A = q %*% a %*% Conj(t(z))
#                             B = q %*% b %*% Conj(t(z))
# ############################################################################
qz <- function(A, B)
{
    size <- dim(A)[1]

    ALPHA <- vector("complex", size)
    BETA <- vector("complex", size)
    VSL <- rep( 0 + 0i, size^2)
    VSR <- rep( 0 + 0i, size^2)
    info <- 0

    # Calling ZGGES using zgges_interface wrapper
    out <- .C("zgges_interface",
              as.integer(size),
              as.integer(info),
              as.complex(A),
              as.complex(B),
              as.complex(ALPHA),
              as.complex(BETA),
              as.complex(VSL),
              as.complex(VSR),
              PACKAGE = "gEcon")

    # Extracting relevant outputs
    output <- list(a = matrix(out[[3]], size, size),
                   b = matrix(out[[4]], size, size),
                   q = matrix(out[[7]], size, size),
                   z = matrix(out[[8]], size, size),
                   gev = matrix(c(out[[5]] ,out[[6]]),
                                nrow = size, ncol = 2),
                   rc = out[[2]])

    return (output)
}
