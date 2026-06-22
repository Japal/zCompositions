lrEM <- function(X,label=NULL,dl=NULL,rob=FALSE,ini.cov=c("complete.obs","multRepl"),frac=0.65,tolerance=0.0001,
                 max.iter=50,rlm.maxit=150,imp.missing=FALSE,suppress.print=FALSE,
                 closure=NULL,z.warning=0.8,z.delete=TRUE,delta=NULL){
  
  if (any(X<0, na.rm=TRUE)) stop("X contains negative values")
  
  if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
    stop("X must be a data matrix with at least two rows and two columns")
  }
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (!any(X==label,na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=TRUE)) stop("Zero values not labelled as censored or missing values were found in the data set")
    if (any(is.na(X))) stop(paste("NA values not labelled as censored or missing values were found in the data set"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=TRUE)) stop("Zero values not labelled as censored or missing values were found in the data set")
    if (!any(is.na(X),na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
  }
  
  if (imp.missing==FALSE){
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.null(dl)){ 
      if (is.na(label)) {
        dl <- apply(X, 2, function(x) min(x[!is.na(x)], na.rm = TRUE))
      } else {
        dl <- apply(X, 2, function(x) min(x[x != label], na.rm = TRUE))
      }
      warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
    }
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    dl <- as.matrix(dl)
    if (!is.numeric(dl)) {
      stop("dl must be numeric")
    }
  }
  
  if (imp.missing==FALSE){
    if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
    if ((nrow(dl)>1) && (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  }
  
  if (!missing(delta)){
    warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
    frac <- delta
  }
  
  if (!is.numeric(frac) || length(frac) != 1 || is.na(frac) || frac <= 0 || frac > 1) {
    stop("frac must be a single numeric value in the interval (0, 1]")
  }
  
  if (!is.numeric(z.warning) || length(z.warning) != 1 || is.na(z.warning) ||
      z.warning < 0 || z.warning > 1) {
    stop("z.warning must be a single numeric value between 0 and 1")
  }
  
  if (!is.numeric(tolerance) || length(tolerance) != 1 || is.na(tolerance) || tolerance <= 0) {
    stop("tolerance must be a positive numeric value")
  }
  
  if (!is.numeric(max.iter) || length(max.iter) != 1 || is.na(max.iter) ||
      max.iter < 1 || max.iter != as.integer(max.iter)) {
    stop("max.iter must be a positive integer")
  }
  
  if (!is.numeric(rlm.maxit) || length(rlm.maxit) != 1 || is.na(rlm.maxit) ||
      rlm.maxit < 1 || rlm.maxit != as.integer(rlm.maxit)) {
    stop("rlm.maxit must be a positive integer")
  }
  
  ini.cov <- match.arg(ini.cov)
  
  lm.sweep <- function(M,C,varobs){
    
    sweep.matrix <- function(A,ind){ 
      
      D <- ncol(A)
      S <- A
      
      for (j in ind){                       
        S[j,j] <- -1/A[j,j]     
        for (i in 1:D) {           
          if (i != j){
            S[i,j] <- -A[i,j]*S[j,j]
            S[j,i] <- S[i,j]
          }
        }
        for (i in 1:D){           
          if (i != j){
            for (k  in 1:D){
              if (k != j){
                S[i,k] <- A[i,k] - S[i,j]*A[j,k]
                S[k,i] <- S[i,k]
              }
            }
          }
        }
        A <- S 
      }
      return(A)
    }  
    
    D <- length(M)                      
    q <- length(varobs)                     
    i <- rep(1,D)                  
    i[varobs] <- i[varobs]-1
    dep <- which(i!=0)
    
    A <- matrix(0,D+1,D+1)             
    A[1,1] <- -1                    
    A[1,2:(D+1)] <- M                
    A[2:(D+1),1] <- matrix(M,ncol=1)
    A[2:(D+1),2:(D+1)] <- C             
    
    reor <- c(1,varobs+1,dep+1)     
    A <- A[reor,reor]              
    A <- sweep.matrix(A,1:(q+1))             
    
    B <- A[1:(q+1),(q+2):(D+1)]       
    CR <- A[(q+2):(D+1),(q+2):(D+1)]       
    
    return(list(betas=B,resid=CR))
  }
  
  inv.alr <- function(x,pos){
    
    ad<-1/(rowSums(exp(x))+1)
    ax<-exp(x)*ad
    if(pos==1) {
      a<-cbind(ad,ax)
    }
    else { 
      if (dim(x)[2] < pos){
        a<-cbind(ax,ad)
      }   
      else {
        a<-cbind(ax[,1:(pos-1)],ad,ax[,pos:(dim(x)[2])])
      }
    }
    return(a)
  }
  
  ilr <- function (x){
    
    D <- length(x)
    z <- vector(mode="double",length=D-1)
    
    for (i in 1:(D-1)){
      z[i] <- sqrt((D-i)/(D-i+1))*log(x[i]/((prod(x[(i+1):D]))^(1/(D-i))))
    }
    return(as.numeric(z))
  }
  
  inv.ilr <- function(z){
    
    D <- length(z) + 1
    x <- vector(mode="double",length=D)
    
    x[1] <- exp(sqrt((D-1)/D)*z[1])
    for (i in 2:(D-1)){
      x[i] <- exp(-sum((1/sqrt((D-1:(i-1)+1)*(D-1:(i-1))))*z[1:(i-1)]) + sqrt(D-i)/sqrt(D-i+1)*z[i]) 
    }
    x[D] <- exp(-sum((1/sqrt((D-1:(D-1)+1)*(D-1:(D-1))))*z[1:(D-1)]))
    
    x <- x/sum(x)
    return(x)
  }
  
  X <- as.data.frame(X)
  nn <- nrow(X); D <- ncol(X)
  if (nn <= D) stop("The lrEM algorithm works on regular data sets (no. rows > no. columns). You can consider lrSVD for wide data sets.")
  
  if (!is.na(label)) {
    X[X == label] <- NA
  }
  X <- as.data.frame(apply(X,2,as.numeric))
  c <- rowSums(X, na.rm=TRUE)
  
  checkNumZerosCol <- colSums(is.na(X))
  
  if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
    cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (ncol(X)-2)) {
        stop(paste("Almost all columns contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }      
      X <- X[,-cases, drop = FALSE]
      if (imp.missing == FALSE) {
        dl <- dl[, -cases, drop = FALSE]
      }
      action <- "deleted"
      
      warning(paste("Column no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                    sep=""))
    } else {      
      action <- "found"      
      warning(paste("Column no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                    sep=""))      
    }
  }
  
  checkNumZerosRow <- rowSums(is.na(X))  
  if (any(checkNumZerosRow/ncol(X) > z.warning)) {    
    cases <- which(checkNumZerosRow/ncol(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (nrow(X)-2)) {
        stop(paste("Almost all rows contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }
      X <- X[-cases, , drop = FALSE]   
      if (imp.missing == FALSE && nrow(dl) > 1) {
        dl <- dl[-cases, , drop = FALSE]
      }
      action <- "deleted"
      
      warning(paste("Row no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                    sep=""))
    } else {      
      action <- "found"      
      warning(paste("Row no. ", cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action,
                    " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                    sep=""))      
    }
  }
  
  # Update dimensions after potential deletion of rows/columns
  nn <- nrow(X)
  D <- ncol(X)
  c <- rowSums(X, na.rm = TRUE)
  
  if (nn <= D) {
    stop("The lrEM algorithm works on regular data sets (no. rows > no. columns). You can consider lrSVD for wide data sets.")
  }
  
  
  if (imp.missing==FALSE) {if (nrow(dl)==1) dl <- matrix(dl, nrow=nn, ncol=D, byrow=TRUE)}
  
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
  
  misspat <- as.data.frame(is.na(X)*1)
  misspat <- as.factor(do.call(paste,c(misspat,sep="")))
  levels(misspat) <- 1:(length(levels(misspat)))
  
  if (any(X <= 0, na.rm = TRUE)) {
    stop("Observed values must be strictly positive before log-ratio transformation")
  }
  
  if (imp.missing == FALSE && any(dl < 0, na.rm = TRUE)) {
    stop("Detection limits in dl must be strictly positive")
  }
  
  if (rob==FALSE){
    
    pos <- which(!is.na(colSums(X)))[1]
    if (is.na(pos)) stop("lrEM based on alr requires at least one complete column")
    
    if (imp.missing==FALSE){
      cpoints <- log(dl)-log(X[,pos])-.Machine$double.eps
      cpoints <- cpoints[,-pos, drop = FALSE]
    }
    
    X_alr <- log(X)-log(X[,pos]); X_alr <- as.matrix(X_alr[,-pos, drop = FALSE])
    nn <- nrow(X_alr); D <- ncol(X_alr)
    
    if (ini.cov != "multRepl"){
      if (inherits(try(solve(cov(X_alr,use=ini.cov)),silent=TRUE),"try-error"))
        stop("ini.cov: singular initial covariance matrix. Probably too few complete rows in data set for using 'complete.obs'")
      M <- matrix(colMeans(X_alr,na.rm=TRUE),ncol=1)
      C <- cov(X_alr,use=ini.cov)}
    else {
      X.mr <- multRepl(X,label=NA,dl=dl,frac=frac,imp.missing=imp.missing,closure=closure,z.warning=z.warning,z.delete=z.delete)
      if (any(X.mr < 0)) {stop("ini.cov: negative values produced using multRepl (please check out closure argument and multRepl help for advice)")}
      X.mr_alr <- (log(X.mr) - log(X.mr[, pos]))[, -pos, drop = FALSE]
      M <- matrix(colMeans(X.mr_alr,na.rm=TRUE),ncol=1)
      C <- cov(X.mr_alr)
    }  
    
    iter_again <- 1
    niters <- 0
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    while (iter_again == 1){
      
      niters <- niters + 1
      Mnew <- M                       
      Cnew <- C
      Y <- X_alr                              
      v <- matrix(0,D,D)
      
      for (npat in 1:length(levels(misspat))){
        i <- which(misspat==npat) 
        varmiss <- which(is.na(X_alr[i[1],]))
        if (length(varmiss) == 0) {next} 
        varobs <- which(!is.na(X_alr[i[1],]))
        if (length(varobs) == 0){
          alt.in <- TRUE
          if (imp.missing == FALSE) {
            temp <- multRepl(
              X[i, , drop = FALSE],
              label = NA,
              dl = dl[i, , drop = FALSE],
              frac = frac,
              imp.missing = imp.missing,
              closure = closure,
              z.warning = z.warning,
              z.delete = z.delete
            )
          } else {
            temp <- multRepl(
              X[i, , drop = FALSE],
              label = NA,
              frac = frac,
              imp.missing = imp.missing,
              closure = closure,
              z.warning = z.warning,
              z.delete = z.delete
            )
          }
          Y[i,] <- (log(temp) - log(temp[, pos]))[, -pos, drop = FALSE]
          if (niters == 1){
            alt.pat <- c(alt.pat,npat)
            alt.mr <- list(alt.mr,i)
          }
          break
        }
        sigmas <- matrix(0,ncol=D)
        B <- matrix(lm.sweep(M,C,varobs)[[1]],ncol=length(varmiss))
        CR <- lm.sweep(M,C,varobs)[[2]]
        Y[i,varmiss] <- matrix(1,nrow=length(i))%*%B[1,] + X_alr[i, varobs, drop=FALSE]%*%B[2:(length(varobs)+1),]
        sigmas[varmiss] <- sqrt(diag(as.matrix(CR)))
        if (imp.missing==FALSE){
          for (j in 1:length(varmiss)){
            sigma <- sigmas[varmiss[j]]
            fdN01 <- dnorm((cpoints[i,varmiss[j]]-Y[i,varmiss[j]])/sigma)
            fdistN01 <- pnorm((cpoints[i,varmiss[j]]-Y[i,varmiss[j]])/sigma)
            ratio <- fdN01 / pmax(fdistN01, .Machine$double.eps)
            Y[i, varmiss[j]] <- Y[i, varmiss[j]] - sigma * ratio
          }
        }
        v[varmiss,varmiss] <- v[varmiss,varmiss] + CR*length(i)
      }
      
      M <- matrix(colMeans(Y),ncol=1)
      dif <- Y - matrix(1,nrow=nn)%*%t(M)              
      PC <- t(dif)%*%dif                       
      C <- (PC+v)/(nn-1)   
      
      Mdif <- max(abs(M-Mnew))    
      Cdif <- max(max(abs(C-Cnew)))  
      if ((max(c(Mdif,Cdif)) < tolerance) || (niters == max.iter)) iter_again <- 0
    }
    
    Y <- inv.alr(Y,pos)
    
    for (i in 1:nn){
      if (any(is.na(X[i,]))){
        vbdl <- which(is.na(X[i,]))
        X[i,vbdl] <- (X[i,pos]/Y[i,pos])*Y[i,vbdl]
      }
    }     
  } 
  
  if (rob==TRUE){
    
    if (ini.cov == "multRepl"){
      if (imp.missing == TRUE){
        X.mr <- multRepl(X,label=NA,imp.missing=TRUE,closure=closure,z.warning=z.warning,z.delete=z.delete)
        if (any(X.mr < 0)) {stop("ini.cov: negative values produced using multRepl (please check out closure argument and multRepl help for advice)")}
      }
      else {X.mr <- multRepl(X,label=NA,dl=dl,frac=frac,closure=closure,z.warning=z.warning,z.delete=z.delete)
      if (any(X.mr < 0)) {stop("ini.cov: negative values produced using multRepl (please check out closure argument and multRepl help for advice)")}
      }
    }
    
    miss <- by(X,misspat,function(x) which(is.na(x[1,])))
    obs <- by(X,misspat,function(x) which(!is.na(x[1,])))
    
    iter_again <- 1
    niters <- 0
    X.old <- X
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    nnn <- 0
    
    while (iter_again == 1){
      
      niters <- niters+1
      if (niters > 1) {X.old <- X; C.old <- C}
      
      for (npat in 1:length(levels(misspat))){
        if (length(miss[[npat]]) == 0) {next} 
        if ((length(obs[[npat]]) == 1) && (!any(npat==alt.pat))){
          alt.in <- TRUE
          if (imp.missing==FALSE){
            X[misspat==npat,] <- multRepl(X.old[misspat==npat,,drop=FALSE],
                                          label=NA,dl=dl[misspat==npat,,drop=FALSE],
                                          frac=frac,closure=closure,z.warning=z.warning,z.delete=z.delete)    
          }
          if (imp.missing==TRUE){
            stop("Please remove samples with only one observed component (check it out using zPatterns).")
          }
          alt.pat <- c(alt.pat,npat)
          alt.mr <- list(alt.mr,which(misspat==npat))
        }
        if (length(obs[[npat]]) > 1) {
          feeder <- X.old[,obs[[npat]]]
          for (m in 1:length(miss[[npat]])){
            p <- miss[[npat]][m]
            target <- X.old[,p]
            if (imp.missing==FALSE){
              phi <- t(apply(cbind(dl=dl[misspat==npat,p],feeder[misspat==npat,]),1,ilr))
            }
            regbasis <- as.data.frame(t(apply(cbind(target,feeder),1,ilr)))
            
            if (niters == 1){
              if (ini.cov == "complete.obs"){
                if (nrow(regbasis[misspat==1,]) > ncol(regbasis[misspat==1,]))
                  robreg <- rlm(V1 ~ .,data=regbasis[misspat==1,],method="MM",maxit = rlm.maxit)
                else
                  stop("ini.cov: singular initial covariance matrix. Probably too few complete rows in data set. Use ini.cov = 'multRepl' instead")
              }
              if (ini.cov == "multRepl"){
                target <- X.mr[,p]
                feeder <- X.mr[,obs[[npat]]]
                regbasis.mr <- as.data.frame(t(apply(cbind(target,feeder),1,ilr)))
                robreg <- rlm(V1 ~ .,data=regbasis.mr,method="MM",maxit = rlm.maxit)
              }  
            }
            else
              robreg <- rlm(V1 ~ .,data=regbasis,method="MM",maxit = rlm.maxit)
            
            B <- matrix(robreg$coefficients,ncol=1)
            sigma <- robreg$s
            est <- cbind(V1=B[1,] + as.matrix(regbasis[misspat==npat,-1])%*%B[-1,],regbasis[misspat==npat,-1])
            if (imp.missing==FALSE){
              
              if (!is.finite(sigma) || sigma <= 0) {
                stop("Non-positive or non-finite sigma encountered in robust censored expectation step")
              }
              
              z <- (phi[, 1] - est[, 1]) / sigma
              pz <- pnorm(z)
              dz <- dnorm(z)
              ratio <- dz / pmax(pz, 1e-12) # Instead of using .Machine$double.eps for increased stability
              est[, 1] <- est[, 1] - sigma * ratio
            }
            est <- t(apply(est,1,inv.ilr))
            est <- est[,1]*(feeder[misspat==npat,1]/est[,2])
            X[misspat==npat,p] <- est
          }
        }    
      }
      
      C <- cov(t(apply(X,1,ilr)))
      
      if (niters >= max.iter) {
        iter_again <- 0
      } else if (niters > 1 && norm(C - C.old, type = "F") < tolerance) {
        iter_again <- 0
      }
      
    } 
  } 
  
  if (closed==1){
    X <- (X / rowSums(X)) * c[1]
  }
  
  if (suppress.print==FALSE){
    if (alt.in) {
      if (imp.missing==FALSE){
        warning("Censoring patterns with only one observed component in the data set.")
        cat("Censored samples with only one observed component imputed by simple multiplicative replacement. \n")
        for (i in 2:length(alt.pat)){
          cat("Row numbers: "); cat(alt.mr[[i]]); cat("\n\n")
        }
      }
    }
    cat(paste("No. iterations to converge: ",niters,"\n\n"))
  }
  
  return(as.data.frame(X))  
  
}