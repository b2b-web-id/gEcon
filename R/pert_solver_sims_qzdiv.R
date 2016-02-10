qzdiv <- function (stake,qzlist)  {
  ##
  ## Takes U.T. matrices a, b, orthonormal matrices q,z, rearranges them
  ## so that all cases of abs(b(i,i)/a(i,i))>stake are in lower right 
  ## corner, while preserving U.T. and orthonormal properties and qaz' and
  ## qbz'.  
  ##
  ## by Christopher A. Sims, 2/22/2004, based on earlier matlab code finished
  ## 7/27/00
  a <- qzlist$a
  b <- qzlist$b
  q <- qzlist$q
  z <- qzlist$z
  n  <- dim(a)[1]
  root <- abs(cbind(diag(a), diag(b)))
  root[,1] <- root[,1]-(root[,1]<1.e-13)*(root[,1]+root[,2])
  root[,2] <- root[,2]/root[,1]
  for (i in  n:1) {
    m <- 0
    for (j in i:1) {
      if (root[j,2] > stake || root[j,2] < -.1) {
        m <- j
        break                           #found an unstable root.  Now check for stable ones below.
      }
    }
    if (m == 0) {
      break                             #quit sort because only stable roots left
    } else {
      if (m<i) {
        for (k in m:(i-1)) {
          qzlist <- qzswitch(k,qzlist)
          tmp <- root[k,2]
          root[k,2] <- root[k+1,2]
          root[k+1,2] <- tmp
        }
      }         
    }
  }
  return(qzlist)
}
