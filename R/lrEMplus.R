lrEMplus <- function(X, dl = NULL, rob = FALSE, ini.cov = c("complete.obs", "multRepl"), frac = 0.65,
                     tolerance = 0.0001, max.iter = 50,
                     rlm.maxit=150, suppress.print = FALSE, closure=NULL, z.warning=0.8, z.delete=TRUE, delta=NULL){
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.null(dl)){ 
    dl <- apply(X, 2, function(x) min(x[!(x %in% c(0, NA))], na.rm = TRUE))
    warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
  }
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)
  dl <- as.matrix(dl) 
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
    n_valid <- if(na.rm) sum(!is.na(x)) else length(x)
    exp(sum(log(x), na.rm=na.rm) / n_valid)
  }
  
  X <- as.data.frame(X,stringsAsFactors=TRUE)
  nn <- nrow(X); D <- ncol(X)
  X <- as.data.frame(apply(X,2,as.numeric),stringsAsFactors=TRUE)
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
      X <- X[,-cases]      
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
      X <- X[-cases,]      
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
  
  if (nrow(dl)==1) dl <- matrix(dl, nrow=nn, ncol=D, byrow=TRUE)
  
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
  
  num_na <- sum(is.na(X))
  num_zeros <- sum(X==0, na.rm=TRUE)
  
  if (num_na > num_zeros){
    X.old <- X
    zeros_mask <- X.old == 0 & !is.na(X.old)
    X.old[zeros_mask] <- frac * dl[zeros_mask]
    
    X.old <- lrEM(X.old, label = NA, imp.missing = TRUE, ini.cov = ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit,
                  suppress.print = TRUE, closure = closure, z.warning = z.warning,
                  z.delete = z.delete)
  }
  
  if (num_na <= num_zeros){
    X.old <- X
    gmeans <- apply(X.old, 2, function(x) gm(x[!(x %in% c(0, NA))], na.rm=TRUE))
    na_mask <- is.na(X.old)
    gmeans_mat <- matrix(gmeans, nrow=nn, ncol=D, byrow=TRUE)
    X.old[na_mask] <- gmeans_mat[na_mask]
    
    X.old <- lrEM(X.old, label = 0, dl = dl, ini.cov = ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit,
                  suppress.print = TRUE, closure = closure, z.warning = z.warning,
                  z.delete = z.delete)
  }
  
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
                  closure = closure, z.warning = z.warning, z.delete = z.delete)
    X.new[is.na(X)] <- NA
    X.new <- lrEM(X.new, label = NA, imp.missing = TRUE, ini.cov =  ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit, suppress.print = TRUE,
                  closure = closure, z.warning = z.warning, z.delete = z.delete)
    
    X.new_alr <- log(X.new)-log(X.new[,D])
    X.new_alr <- as.matrix(X.new_alr[,-D])
    M.new <- matrix(colMeans(X.new_alr),ncol=1)
    C.new <- cov(X.new_alr)
    
    Mdif <- max(abs(M.new-M.old))    
    Cdif <- max(max(abs(C.new-C.old)))  
    if ((max(c(Mdif,Cdif)) < tolerance) | (niters == max.iter)) iter_again <- 0
    
  }
  
  if (closed==1) X.new <- (X.new / rowSums(X.new)) * c[1]
  
  if (suppress.print==FALSE) cat(paste("No. iterations to converge: ",niters,"\n\n"))
  
  return(as.data.frame(X.new,stringsAsFactors=TRUE))
  
}