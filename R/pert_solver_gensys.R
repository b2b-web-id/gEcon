gensys <- function(g0, g1, c0=matrix(0,dim(g0)[1]), psi, pi, div=-1)
  {
    ##System given as
    ##        g0*y(t)=g1*y(t-1)+c0+psi*z(t)+pi*eta(t),
    ##with z an exogenous variable process and eta being endogenously determined
    ##one-step-ahead expectational errors.  Returned system is
    ##       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) + loose*eta .
    ## If z(t) is i.i.d., the term involving fmat and fwt drops out.
    ## If the solution is unique (eu[2]==1) there is no "loose" term.  Otherwise
    ## loose characterizes the dimensions along which there is non-uniqueness.
    ## If div is omitted from argument list, a div>1 is calculated.
    ## eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu=[-2,-2] for coincident zeros.
    ## By Christopher A. Sims 2/24/2004, from earlier matlab code of same author.
    eu <- c(0,0)
    realsmall <- 1e-7
    fixdiv <- (div>0)
    n <- dim(g0)[1]
    nshock <- if (is.matrix(psi)) dim(psi)[2] else if (is.null(psi)) 0 else 1
    qzl <- qz(g0,g1)

    # Checks if any of a, b matrices diagonals from QZ decomposition are zeros
    zxz <- any((abs(diag(qzl$a))<realsmall) & (abs(diag(qzl$b))<realsmall))
    if (zxz) {
      "Coincident zeros.  Indeterminacy and/or nonexistence.\n"
      eu <- c(-2,-2)
      gev <- qzl$gev
      # any of a, b are diagonal
      return (list(eu=eu,gev=gev))
    }

    zeroax <- abs(diag(qzl$a)) < realsmall
    unstabx <- abs(diag(qzl$a)) < (1-realsmall)*abs(diag(qzl$b)) # near unit roots don't count
    unstabx <- (! zeroax) & unstabx
    if (! fixdiv) {
      if (! any(unstabx)){
        div <- 1.01
      } else
      {
        div <- .5*(min(abs(diag(qzl$b)[unstabx]/diag(qzl$a)[unstabx]))+1)
      }
    }
    unstabx <- div*abs(diag(qzl$a))<= abs(diag(qzl$b))
    nunstab <- sum(unstabx)
    qzl <- qzdiv(div,qzl)
    qq <- t(Conj(qzl$q))                # to match matlab convention
    gev <- qzl$gev

    ## note that this means that gev is not simply the diagonals of a nd b.  qzdiv
    ## changes the numbers on the diagonals (though not their ratios), but merely reorders
    ## the original gev.
    if (nunstab==n){
      six <- NULL
      uix <- 1:n
    } else
    {
      if (nunstab==0){
        uix <- NULL
        six <- 1:n
      } else
      {
        uix <- (n-nunstab+1):n
        six <- 1:(n-nunstab)
      }
    }
    q1 <- qq[six,,drop=FALSE]
    q2 <- qq[uix,,drop=FALSE]
    z1 <- t(Conj(qzl$z[,six,drop=FALSE]))
    z2 <- t(Conj(qzl$z[,uix,drop=FALSE]))
    a2 <- qzl$a[uix,uix,drop=FALSE]
    b2 <- qzl$b[uix,uix,drop=FALSE]
    ## debug
    ## browser()
    etawt <- q2 %*% pi
    neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
    ndeta <- min(nunstab,neta)
    if(ndeta==0){
      ueta <- matrix(0,nunstab,0)
      deta <- vector("numeric",0)
      veta <- matrix(0,neta,0)
      bigev <- vector("logical",0)
    } else {
      sd <- svd(etawt)
      ueta <- sd$u; deta <- sd$d; veta <- sd$v
      bigev <- deta>realsmall
      ueta<-ueta[,bigev,drop=FALSE]
      veta<-veta[,bigev,drop=FALSE]
      deta<-deta[bigev]
    }
    eu[1] <- sum(bigev) >= nunstab
    ##----------------------------------------------------
    ## Note that existence and uniqueness are not just matters of comparing
    ## numbers of roots and numbers of endogenous errors.  These counts are
    ## reported below because usually they point to the source of the problem.
    ##------------------------------------------------------
    etawt1 <- q1 %*% pi
    ndeta1 <- min(n-nunstab,neta)
    if(ndeta1==0){
      ueta1 <- matrix(0,n-nunstab,0)
      deta1 <- vector("numeric",0)
      veta1 <- matrix(0,neta,0)
      bigev1 <- vector("logical",0)
    } else {
      sd <- svd(etawt1)
      ueta1<-sd$u
      deta1 <- sd$d
      veta1 <- sd$v
      bigev1 <- deta1 > realsmall
    }
    if (any(bigev1)) { #needed because empty dimensions are dropped after select
      ueta1 <- ueta1[,bigev1,drop=FALSE]
      veta1 <- veta1[,bigev1,drop=FALSE]
      deta1 <- deta1[bigev1]
      loose <- veta1-veta %*% t(Conj(veta)) %*% veta1
      svdl <- svd(loose)
      loose <- sum(abs(svdl$d)>realsmall*n)
      unq <- (loose==0)
    } else {
      ueta1 <- matrix(1,n-nunstab,0)
      veta1 <- matrix(1,neta,0)
      deta1 <- vector("complex",0)
      unq <- TRUE
    }
    if (unq) {
      eu[2] <- 1
    } else
    {
      cat("Indeterminacy.", loose, "loose endog errors.\n")
    }
    ## Note: if v is a vector of length n and m is an nxp matrix,
    ## v*m==diag(v)%*%m, m/v==solve(diag(v),m)==diag(v)\m (matlab notation)
    ##
    tmat <- cbind(diag(n-nunstab),
                  -t(Conj((ueta %*% (t(Conj(veta))/deta)) %*% veta1 %*% (deta1 * t(Conj(ueta1)))))  )
    G0<- rbind( tmat %*% qzl$a, cbind(matrix(0,nunstab,n-nunstab), diag(nunstab)))
    G1<- rbind(tmat %*% qzl$b, matrix(0,nunstab,n))
    ##----------------------
    ## G0 is always non-singular because by construction there are no zeros on
    ## the diagonal of a[1:(n-nunstab),1:(n-nunstab)], which forms G0's ul corner.
    ##-----------------------
    G0I <- solve(G0)
    G1 <- G0I%*%G1
    ##----------- uix can be empty, e.g. in indeterminate systems with no unstable roots ------------
    if(is.null(uix)){
      C <- G0I %*% tmat %*% qq %*% c0
      fmat <- matrix(0,0,0)
      fwt <- matrix(0, 0, nshock)
      impact <- G0I %*% tmat %*% qq %*% psi
    }else{
      C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(qzl$a[uix,uix,drop=FALSE]-qzl$b[uix,uix,drop=FALSE],q2%*%c0) )
      impact <- G0I %*% rbind(tmat %*% qq %*% psi, matrix(0,nunstab, nshock))
      fmat <- solve(qzl$b[uix,uix,drop=FALSE],qzl$a[uix,uix,drop=FALSE])
      fwt <- -solve(qzl$b[uix,uix,drop=FALSE],q2 %*% psi)
    }
    ywt <- G0I[,uix,drop=FALSE]
    ##loose <- G0I %*% etawt1 %*% (diag(neta) - veta %*% t(Conj(veta)))
    loose <- G0I %*% qq %*% pi %*% (diag(neta) - veta %*% t(Conj(veta)))
    ## loose <- G0I %*% rbind(loose,matrix(0,nunstab,neta))  #(think the above is a mistaken remnant)
    ##-------------------- above are output for system in terms of z'y -------
    G1 <- Re(qzl$z %*% G1 %*% t(Conj(qzl$z)))
    C <- Re(qzl$z%*%C)
    impact <- Re(qzl$z%*%impact)
    ywt <- qzl$z%*%ywt
    loose <- Re(qzl$z %*% loose)
    vn <- dimnames(g0)[[2]]
    dimnames(G1) <- list(vn,vn)
    dimnames(C) <- list(vn,NULL)
    dimnames(impact)[[1]] <- vn
    dimnames(ywt)[[1]] <- vn
    dimnames(loose)[[1]] <- vn
    ## qzl added
    return (list(G1=G1,C=C,impact=impact,fmat=fmat,fwt=fwt,
                ywt=ywt,gev=gev,eu=eu,loose=loose, qzl=qzl, nunstab = nunstab))
  }

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
  return (qzlist)
}

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
    if (m == 0) return (qzlist)
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
  return (qzlist)
}

qzswitch <- function(i=1,qzlist)
  {
    ## Takes U.T. matrices a, b, orthonormal matrices q,z, interchanges
    ## diagonal elements i and i+1 of both a and b, while maintaining
    ## qaz' and qbz' unchanged.  If diagonal elements of a and b
    ## are zero at matching positions, the returned a will have zeros at both
    ## positions on the diagonal.  This is natural behavior if this routine is used
    ## to drive all zeros on the diagonal of a to the lower right, but in this case
    ## the qz transformation is not unique and it is not possible simply to switch
    ## the positions of the diagonal elements of both a and b.
    realsmall <- 1e-7;
    ##realsmall<-1e-3;
    a <- qzlist$a
    b <- qzlist$b
    q <- qzlist$q
    z <- qzlist$z
    gev <- qzlist$gev
    q <- t(Conj(q))                     #This is needed because the code was originally for matlab, where it is q'az' that
                                        # is preserved.
    A <- a[i,i]; d <- b[i,i]; B <- a[i,(i+1)]; e <- b[i,(i+1)]
    g <- a[i+1,i+1]; f <- b[i+1,i+1]
    w2 <- rep(0,2)
    ## a[i:(i+1),i:(i+1)]<-[A B; 0 g]
    ## b[i:(i+1),i:(i+1)]<-[d e; 0 f]
    if (abs(g)<realsmall & abs(f)<realsmall)
      {
        if (abs(A)<realsmall)
          {
            ## l.r. coincident 0's with u.l. of a<-0; do nothing
            return (list(a=a,b=b,q=q,z=z))
          } else
        {
          ## l.r. coincident zeros; put 0 in u.l. of a
          wz <- c(B, -A)
          wz <- wz/sqrt(sum(Conj(wz)*wz))
          wz <- array(c(wz,Conj(wz[2]),-Conj(wz[1])),dim=c(2,2))
          xy <- diag(2)
        }
      } else
    {
      if (abs(A)<realsmall && abs(d)<realsmall)
        {
          if (abs(g)<realsmall)
            {
              ## u.l. coincident zeros with l.r. of a<-0; do nothing
              return (list(a=a,b=b,q=q,z=z))
            } else
          {
            ## u.l. coincident zeros; put 0 in l.r. of a
            wz <- diag(2)
            xy <- c(g,-B)
            xy <- xy/sqrt(sum(xy*Conj(xy)))
            xy <- t(matrix(c(Conj(xy[2]), -Conj(xy[1]),xy),nrow=2,ncol=2))
          }
        } else
      {
        ## usual case
        wz <- c(g*e-f*B, Conj(g*d-f*A))
        xy <- c(Conj(B*d-e*A), Conj(g*d-f*A))
        n <- sqrt(wz %*% Conj(wz))
        m <- sqrt(xy %*% Conj(xy));
        if (Re(m)<1e-12*100)
          {
            ## all elements of a and b proportional
            return (list(a=a,b=b,q=q,z=z))
          }
        wz <- wz/n
        xy <- xy/m
        wz <- matrix(c(wz, -Conj(wz[2]),Conj(wz[1])),byrow=TRUE,ncol=2,nrow=2)
        xy <- matrix(c(xy,-Conj(xy[2]), Conj(xy[1])),byrow=TRUE,ncol=2,nrow=2)
      }
    }
    a[i:(i+1),] <- xy %*% a[i:(i+1),]
    b[i:(i+1),] <- xy %*% b[i:(i+1),]
    a[,i:(i+1)] <- a[,i:(i+1)] %*% wz
    b[,i:(i+1)] <- b[,i:(i+1)] %*% wz
    z[,i:(i+1)] <- z[,i:(i+1)] %*% wz
    q[i:(i+1),] <- xy %*% q[i:(i+1),]
    q <- t(Conj(q))
    w2 <- gev[i,]
    gev[i,] <- gev[i+1,]
    gev[i+1,] <- w2
    qzlist$a <- a
    qzlist$b <- b
    qzlist$q <- q
    qzlist$z <- z
    qzlist$gev <- gev
    return (qzlist)
  }
