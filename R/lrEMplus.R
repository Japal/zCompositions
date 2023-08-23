lrEMplus <- function(X, dl = NULL, rob = FALSE, ini.cov = c("complete.obs", "multRepl"), frac = 0.65,
                     tolerance = 0.0001, max.iter = 50,
                     rlm.maxit=150, suppress.print = FALSE, closure=NULL, z.warning=0.8, delta=NULL){
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.character(dl) || is.null(dl)) stop("dl must be a numeric vector or matrix")
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)
  dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")

  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  
  if (any(is.na(X))==FALSE) stop("No missing data were found in the data set")
  if (any(X==0, na.rm=T)==FALSE) stop("No zeros were found in the data set")

  if (!missing("delta")){
    warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
    frac <- delta
  }
  
  ini.cov <- match.arg(ini.cov)
  
  gm <- function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  
  ## Preliminaries ----

  X <- as.data.frame(X,stringsAsFactors=TRUE)
  nn <- nrow(X); D <- ncol(X)
  X <- as.data.frame(apply(X,2,as.numeric),stringsAsFactors=TRUE)
  c <- apply(X,1,sum,na.rm=TRUE)
  
  # Number of zeros or missing per column for warning
  checkNumZerosCol <- apply(X,2,function(x) sum(is.na(x) | (x==0)))
  if (any(checkNumZerosCol/nrow(X) >= z.warning)) {
    cases <- which(checkNumZerosCol/nrow(X) >= z.warning)
    X <- X[,-cases]
    warning(paste("Column ",cases," containing more than ",z.warning*100,"% zeros/unobserved values was deleted (pre-check out using function zPatterns/modify threshold using argument z.warning).\n",sep=""))
  }
  
  checkNumZerosRow <- apply(X,1,function(x) sum(is.na(x) | (x==0)))
  if (any(checkNumZerosRow/ncol(X) >= z.warning)) {
    cases <- which(checkNumZerosRow/ncol(X) >= z.warning)
    X <- X[-cases,]
    warning(paste("Row ",cases," containing more than ",z.warning*100,"% zeros/unobserved values was deleted (pre-check out using function zPatterns/modify threshold using argument z.warning).\n",sep=""))
  }

  if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl

  # Check for closure
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1

  if (sum(is.na(X)) > sum(X==0,na.rm=T)){
    X.old <- X
    # Initial simple imputation of zero
    for (i in 1:nn){
      if (any(X.old[i, ]==0,na.rm=T)){
        z <- which(X.old[i, ]==0)
        X.old[i,z] <- frac*dl[i,z]
      }
    }
    # Initial lrEM imputation of missing data
    X.old <- lrEM(X.old, label = NA, imp.missing = TRUE, ini.cov = ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit,
                  suppress.print = TRUE, closure = closure)
  }
  
  if (sum(is.na(X)) <= sum(X==0,na.rm=T)){
    X.old <- X
    # Initial ordinary geo mean imputation of missing (ignores 0s in the column if any)
    gmeans <- apply(X.old,2,function(x) gm(x[x!=0]))
    for (i in 1:nn){
      if (any(is.na(X.old[i, ]))){
        z <- which(is.na(X.old[i, ]))
        X.old[i,z] <- gmeans[z]
      }
    }
    # Initial lrEM imputation of zeros
    X.old <- lrEM(X.old, label = 0, dl = dl, ini.cov = ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit,
                  suppress.print = TRUE, closure = closure)
  }

  # Initial parameter estimates
  X.old_alr <- log(X.old)-log(X.old[,D])
  X.old_alr <- as.matrix(X.old_alr[,-D])
  M.old <- matrix(colMeans(X.old_alr),ncol=1)
  C.old <- cov(X.old_alr)

  iter_again <- 1
  niters <- 0

  while (iter_again == 1){

    niters <- niters+1
    if (niters > 1) {X.old <- X.new; M.old <- M.new; C.old <- C.new}

    X.old <- as.matrix(X.old)
    X.old[which(X==0)] <- 0
    X.new <- lrEM(X.old, label = 0, dl = dl, ini.cov =  ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit, suppress.print = TRUE,
                  closure = closure)
    X.new[is.na(X)] <- NA
    X.new <- lrEM(X.new, label = NA, imp.missing = TRUE, ini.cov =  ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit, suppress.print = TRUE,
                  closure = closure)

    X.new_alr <- log(X.new)-log(X.new[,D])
    X.new_alr <- as.matrix(X.new_alr[,-D])
    M.new <- matrix(colMeans(X.new_alr),ncol=1)
    C.new <- cov(X.new_alr)
    
    # Convergence check
    Mdif <- max(abs(M.new-M.old))    
    Cdif <- max(max(abs(C.new-C.old)))  
    if ((max(c(Mdif,Cdif)) < tolerance) | (niters == max.iter)) iter_again <- 0
    
  }

  ## Final section ----
  
  if (closed==1) X.new <- t(apply(X.new,1,function(x) x/sum(x)*c[1])) # If not closed lrEM above takes care of it

  if (suppress.print==FALSE) cat(paste("No. iterations to converge: ",niters,"\n\n"))
  
  return(as.data.frame(X.new,stringsAsFactors=TRUE))
  
}