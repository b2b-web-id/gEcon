qzdivct <- function(stake,qzlist) {
  ## Takes qz decomp output, rearranges a and b
  ## so that all cases of real(b(i,i)/a(i,i))>stake are in lower right
  ## corner, while preserving U.T. and orthonormal properties and q'az' and
  ## q'bz'.  abs(a(i,i)) < realsmall is interpreted as a zero and as generating
  ## an infinitely positive real part of the ratio.  All i's for which this
  ## criterion are satisfied are grouped together in the lower right corner
  ## of the lower right corner, with the non-zero roots above them.  This
  ## version differs from
  ## qzdiv in that it works on the real part's value, as is appropriate for
  ## continuous time models, instead of on the absolute value, as is
  ## appropriate for discrete time models.
  ##
  realsmall <- sqrt(.Machine$double.eps)*10
  ##realsmall <- 1e-3
  n   <-  dim(qzlist$a)[1]
  root  <-  cbind(diag(qzlist$a), diag(qzlist$b))
  ## first sort on the non-zero root criterion
  xdown0  <-  abs(root[, 1]) < realsmall
  xdown  <-  xdown0 | (Re(root[ , 2 ] / (xdown0+root[ , 1])) > stake)
  for (i in n:1) {
    m <- 0
    for ( j in seq(i, 1, by=-1) ) {
      if (xdown0[j]) { 
        m <- j
        break
      }
    }
    if ( m == 0) break               #This means we're done with zeros
    if (i > m) {
      for ( k in m:(i-1) ) {
        qzlist  <-  qzswitch(k,qzlist)
        root <- cbind(diag(qzlist$a), diag(qzlist$b))
        xdown0[k:(k+1)] <- xdown0[(k+1):k]
        xdown[k:(k+1)] <- xdown[(k+1):k]
      }
    }
  }
  ## now repeat, using the stake criterion
  for (i in n:1) {
    m <- 0
    for (j in i:1) {
      if (xdown[j]) {
        m <- j
        break
      }
    }
    if (m == 0) return(qzlist)
    if (i > m) {
      for (k in m:(i-1) ) {
        gevOld <- root[k:(k+1), ]
        qzlist  <-  qzswitch(k,qzlist)
        root <- cbind(diag(qzlist$a), diag(qzlist$b))
        xdown0[k:(k+1)] <- xdown0[(k+1):k]
        xdown[k:(k+1)] <- xdown[(k+1):k]
        gev <- root[k:(k+1), ]
        print(rbind(cbind(gevOld, gevOld[ , 2]/gevOld[, 1]), cbind(gev, gev[, 2]/gev[, 1])))
      }
    }
  }
  xdown0f <- abs(root[,1]) < realsmall
  xdownf <- xdown0f | (Re(root[ , 2] / (xdown0f + root[ , 1])) > stake)
  d0 <- which(xor(xdown0f, xdown0))
  d <- which(xor(xdownf, xdown))
  #modified
  if (length(d0) + length(d) > 0 ) warning(paste("root classification shifted with sort: d0 = ", d0, "d1 = ",d))
  return(qzlist)
}
