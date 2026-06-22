multRepl <-
  function(X,label=NULL,dl=NULL,frac=0.65,imp.missing=FALSE,closure=NULL,z.warning=0.8,z.delete=TRUE,delta=NULL){
    
    if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")
    
    if (is.character(X)) stop("X is not a valid data matrix or vector.")
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X == label, na.rm = TRUE)) stop(paste("Label",label,"was not found in the data set"))
      if (label != 0 && any(X == 0, na.rm = TRUE)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored or missing values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X == 0, na.rm = TRUE)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (!any(is.na(X),na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
    }
    
    if (!is.logical(imp.missing) || length(imp.missing) != 1 || is.na(imp.missing)) {
      stop("imp.missing must be TRUE or FALSE")
    }
    
    if (imp.missing==FALSE){
      if (is.character(dl)) stop("dl must be a numeric vector or matrix")
      if (is.null(dl)){
        
        X_for_dl <- if (is.vector(X)) {
          matrix(X, nrow = 1)
        } else {
          X
        }
        
        if (is.na(label)) {
          dl <- apply(X_for_dl, 2, function(x) min(x[!is.na(x)], na.rm = TRUE))
        } else {
          dl <- apply(X_for_dl, 2, function(x) min(x[x != label], na.rm = TRUE))
        }
        
        warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
      }
      if (is.vector(dl)) dl <- matrix(dl,nrow=1)
      dl <- as.matrix(dl) 
      
      if (!is.numeric(dl)) {
        stop("dl must be numeric")
      }
      
      if (any(dl < 0, na.rm = TRUE) || any(!is.finite(dl), na.rm = TRUE)) {
        stop("Detection limits in dl must be non-negative and finite")
      }
      
    }
    if (is.vector(X)){
      if (imp.missing==TRUE) stop("Data matrix required: missing values cannot be imputed in single vectors")
      if (ncol(dl)!=ncol(as.data.frame(matrix(X,ncol=length(X))))) stop("The number of columns in X and dl do not agree")
    }
    if (!is.vector(X)){
      if (imp.missing==FALSE){
        if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
        
        if (nrow(dl) > 1 && nrow(dl) != nrow(X)) {
          stop("The number of rows in X and dl do not agree")
        }
        
      }
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
    
    if (!is.logical(z.delete) || length(z.delete) != 1 || is.na(z.delete)) {
      stop("z.delete must be TRUE or FALSE")
    }
    
    
    gm <- function(x, na.rm=TRUE){
      n_valid <- if(na.rm) sum(!is.na(x)) else length(x)
      exp(sum(log(x), na.rm=na.rm) / n_valid)
    }
    
    nam <- NULL
    
    if (is.vector(X)) {
      if (!is.null(names(X))) {
        nam <- names(X)
      }
      X <- as.data.frame(matrix(X, ncol = length(X)))
    } else {
      if (!is.null(colnames(X))) {
        nam <- colnames(X)
      }
    }
    
    rnames <- rownames(X)
    
    if (!is.na(label)) {
      X[X == label] <- NA
    }
    
    X <- apply(X,2,as.numeric)
    if (!is.null(rownames(X))) {rownames(X) <- rnames} 
    if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)))
    
    if (nrow(X) > 1){
      checkNumZerosCol <- colSums(is.na(X))
      
      if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
        cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
        if (z.delete == TRUE) {
          if (length(cases) > (ncol(X)-2)) {
            stop(paste("Almost all columns contain >", z.warning*100,
                       "% zeros/unobserved values (see arguments z.warning and z.delete).",
                       sep=""))
          }      
          
          X <- X[, -cases, drop = FALSE]
          
          if (imp.missing == FALSE) {
            dl <- dl[, -cases, drop = FALSE]
          }
          
          if (!is.null(nam)) {
            nam <- nam[-cases]
          }
          
          action <- "deleted"
          
          warning(paste("Column no. ", paste(cases, collapse = ", "),
                        " containing >", z.warning * 100,
                        "% zeros/unobserved values ", action,
                        " (see arguments z.warning and z.delete).\n",
                        sep = ""))
        } else {      
          action <- "found"      
          warning(paste("Column no. ", paste(cases, collapse = ", "),
                        " containing >", z.warning * 100,
                        "% zeros/unobserved values ", action,
                        " (see arguments z.warning and z.delete).\n",
                        sep = ""))     
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
          
          warning(paste("Row no. ", paste(cases, collapse = ", "),
                        " containing >", z.warning * 100,
                        "% zeros/unobserved values ", action,
                        " (see arguments z.warning and z.delete).\n",
                        sep = ""))
        } else {      
          action <- "found"      
          warning(paste("Row no. ", paste(cases, collapse = ", "),
                        " containing >", z.warning * 100,
                        "% zeros/unobserved values ", action,
                        " (see arguments z.warning and z.delete).\n",
                        sep = ""))     
        }
      }
    }
    
    nn <- nrow(X); D <- ncol(X)
    
    if (nn < 1 || D < 2) {
      stop("After deleting rows/columns, X must still have at least one row and two columns")
    }
    
    c_base <- rowSums(X, na.rm=TRUE)
    
    closed <- 0
    if (all(abs(c_base - mean(c_base)) < .Machine$double.eps^0.3)) closed <- 1
    
    if (imp.missing==FALSE){
      if (nrow(dl)==1) dl <- matrix(dl, nrow=nn, ncol=D, byrow=TRUE)
    }
    
    if (!is.null(closure)){
      if (!is.numeric(closure) || length(closure) != 1 || is.na(closure) ||
          !is.finite(closure) || closure <= 0) {
        stop("closure must be a single positive finite numeric value")
      }
      
      if (closed == 1) {
        stop("closure: The data are already closed to ", c_base[1])
      }
      
      c_adj <- rep(closure, nn)
    } else {
      c_adj <- c_base
    }
    
    X_na <- is.na(X)
    
    
    if (imp.missing == FALSE && any(dl[X_na] <= 0, na.rm = TRUE)) {
      stop("Detection limits for censored values must be strictly positive")
    }
    
    
    if (imp.missing==FALSE){
      est <- frac * dl
    } else {
      gms <- apply(X, 2, function(x) gm(x[!is.na(x) & x > 0], na.rm = TRUE))
      
      if (any(!is.finite(gms)) || any(gms <= 0)) {
        stop("Cannot compute finite positive geometric means for missing-value imputation")
      }
      
      est <- matrix(gms, nrow = nn, ncol = D, byrow = TRUE)
    }
    
    imputed_only <- est
    imputed_only[!X_na] <- 0
    sum_est <- rowSums(imputed_only)
    
    if (any(!is.finite(c_adj)) || any(c_adj <= 0)) {
      stop("Closure/row sums must be positive and finite")
    }
    
    problem_rows <- which(sum_est >= c_adj)
    
    if (length(problem_rows) > 0) {
      if (imp.missing == FALSE) {
        stop("Estimated replacement values exceed or equal the row total; multiplicative replacement cannot proceed")
      } else {
        warning("Some missing-value estimates exceed or equal the row total; unadjusted geometric-mean estimates used for those rows.")
      }
    }
    
    adjustment <- 1 - (sum_est / c_adj)
    
    if (length(problem_rows) > 0 && imp.missing == TRUE) {
      adjustment[problem_rows] <- 1
    }
    
    if (any(!is.finite(adjustment)) || any(adjustment <= 0)) {
      stop("Non-positive or non-finite adjustment factor encountered")
    }
    
    
    X[X_na] <- (est / adjustment)[X_na]
    
    if (!is.null(nam) && length(nam) == ncol(X)) {
      colnames(X) <- nam
    }
    
    if (closed==1){
      X <- (X / rowSums(X)) * c_base[1]
    } 
    
    if (any(X < 0, na.rm = TRUE)) warning("multRepl: negative imputed values were generated (please check out help for advice)")
    
    return(as.data.frame(X))
  }
