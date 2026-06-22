lrEMplus <- function(X, dl = NULL, rob = FALSE, ini.cov = c("complete.obs", "multRepl"), frac = 0.65,
                     tolerance = 0.0001, max.iter = 50,
                     rlm.maxit=150, suppress.print = FALSE, closure=NULL, z.warning=0.8, z.delete=TRUE, delta=NULL){
  
  if (any(X<0, na.rm=TRUE)) stop("X contains negative values")
  if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
    stop("X must be a data matrix with at least two rows and two columns")
  }
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.null(dl)){ 
    dl <- apply(X, 2, function(x) min(x[!is.na(x) & x != 0], na.rm = TRUE))
    warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
  }
  
  if (is.vector(dl)) dl <- matrix(dl, nrow = 1)
  dl <- as.matrix(dl)
  
  if (!is.numeric(dl)) {
    stop("dl must be numeric")
  }
  
  
  if (any(dl < 0, na.rm = TRUE) || any(!is.finite(dl), na.rm = TRUE)) {
    stop("Detection limits in dl must be strictly positive and finite")
  }
  
  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if (nrow(dl) > 1 && nrow(dl) != nrow(X)) stop("The number of rows in X and dl do not agree")
  
  if (any(is.na(X))==FALSE) stop("No missing data were found in the data set")
  if (any(X==0, na.rm=TRUE)==FALSE) stop("No zeros were found in the data set")
  
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
  
  gm <- function(x, na.rm=TRUE){
    n_valid <- if(na.rm) sum(!is.na(x)) else length(x)
    exp(sum(log(x), na.rm=na.rm) / n_valid)
  }
  
  X <- as.data.frame(X)
  nn <- nrow(X); D <- ncol(X)
  X <- as.data.frame(apply(X,2,as.numeric))
  c <- rowSums(X, na.rm=TRUE)
  
  checkNumZerosCol <- colSums(is.na(X) | X == 0, na.rm = TRUE)
  
  if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
    cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (ncol(X)-2)) {
        stop(paste("Almost all columns contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }      
      
      X <- X[, -cases, drop = FALSE]
      dl <- dl[, -cases, drop = FALSE]
      
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
  
  checkNumZerosRow <- rowSums(is.na(X) | X == 0, na.rm = TRUE)
  if (any(checkNumZerosRow/ncol(X) > z.warning)) {    
    cases <- which(checkNumZerosRow/ncol(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (nrow(X)-2)) {
        stop(paste("Almost all rows contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }
      
      X <- X[-cases, , drop = FALSE]
      
      if (nrow(dl) > 1) {
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
  
  # Update nn, D, and c after potential row/column deletions
  nn <- nrow(X)
  D <- ncol(X)
  c <- rowSums(X, na.rm = TRUE)
  
  if (nn <= D) {
    stop("The lrEMplus algorithm works on regular data sets (no. rows > no. columns). You can consider lrSVD for wide data sets.")
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
                  z.delete = FALSE) # Already handled by lrEMplus
  }
  
  if (num_na <= num_zeros){
    X.old <- X
    gmeans <- apply(X.old, 2, function(x) gm(x[!is.na(x) & x != 0], na.rm = TRUE))
    
    if (any(!is.finite(gmeans))) {
      stop("Cannot compute finite geometric means; at least one column has no positive observed values")
    }
    
    na_mask <- is.na(X.old)
    gmeans_mat <- matrix(gmeans, nrow=nn, ncol=D, byrow=TRUE)
    X.old[na_mask] <- gmeans_mat[na_mask]
    
    X.old <- lrEM(X.old, label = 0, dl = dl, ini.cov = ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit,
                  suppress.print = TRUE, closure = closure, z.warning = z.warning,
                  z.delete = FALSE)
  }
  
  if (any(X.old <= 0, na.rm = TRUE)) {
    stop("X.old contains non-positive values before log-ratio transformation")
  }
  
  X.old_alr <- log(X.old)-log(X.old[,D])
  X.old_alr <- as.matrix(X.old_alr[,-D, drop = FALSE])
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
                  closure = closure, z.warning = z.warning, z.delete = FALSE)
    X.new[is.na(X)] <- NA
    X.new <- lrEM(X.new, label = NA, imp.missing = TRUE, ini.cov =  ini.cov, rob = rob,
                  tolerance = tolerance, max.iter = max.iter, rlm.maxit = rlm.maxit, suppress.print = TRUE,
                  closure = closure, z.warning = z.warning, z.delete = FALSE)
    
    if (any(X.new <= 0, na.rm = TRUE)) {
      stop("X.new contains non-positive values before log-ratio transformation")
    }
    
    X.new_alr <- log(X.new)-log(X.new[,D])
    X.new_alr <- as.matrix(X.new_alr[,-D, drop = FALSE])
    M.new <- matrix(colMeans(X.new_alr),ncol=1)
    C.new <- cov(X.new_alr)
    
    Mdif <- max(abs(M.new-M.old))    
    Cdif <- max(max(abs(C.new-C.old)))  
    if (max(c(Mdif, Cdif)) < tolerance || niters == max.iter) iter_again <- 0
    
  }
  
  if (closed==1) X.new <- (X.new / rowSums(X.new)) * c[1]
  
  if (suppress.print==FALSE) cat(paste("No. iterations to converge: ",niters,"\n\n"))
  
  return(as.data.frame(X.new))
  
}