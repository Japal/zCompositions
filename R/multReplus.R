multReplus <- function(X, dl = NULL, frac = 0.65, closure = NULL, z.warning = 0.8, z.delete = TRUE, delta = NULL){
  
  
  if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")
  
  if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
    stop("X must be a data matrix with at least two rows and two columns")
  }
  
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  
  if (is.null(dl)){ 
    dl <- apply(X, 2, function(x) min(x[!is.na(x) & x != 0], na.rm = TRUE))
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
  
  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  
  if (nrow(dl) > 1 && nrow(dl) != nrow(X)) {
    stop("The number of rows in X and dl do not agree")
  }
  
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
  
  if (!is.logical(z.delete) || length(z.delete) != 1 || is.na(z.delete)) {
    stop("z.delete must be TRUE or FALSE")
  }
  
  if (!is.null(closure)) {
    if (!is.numeric(closure) || length(closure) != 1 || is.na(closure) ||
        !is.finite(closure) || closure <= 0) {
      stop("closure must be a single positive finite numeric value")
    }
  }
  
  gm <- function(x, na.rm=TRUE){
    n_valid <- if(na.rm) sum(!is.na(x)) else length(x)
    exp(sum(log(x), na.rm=na.rm) / n_valid)
  }
  
  
  nam <- NULL
  if (!is.null(colnames(X))) {
    nam <- colnames(X)
  }
  
  if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)))
  
  X <- as.data.frame(X)
  X <- as.data.frame(apply(X, 2, as.numeric))
  
  zero_or_na <- is.na(X) | X == 0
  checkNumZerosCol <- colSums(zero_or_na)
                              
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
  
  
  zero_or_na <- is.na(X) | X == 0
  checkNumZerosRow <- rowSums(zero_or_na)
  
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
  
  nn <- nrow(X)
  D <- ncol(X)
  
  if (nn < 2 || D < 2) {
    stop("After deleting rows/columns, X must still have at least two rows and two columns")
  }
  
  c_base <- rowSums(X, na.rm = TRUE)

  if (!any(is.na(X))) {
    warning("No missing values remain after deleting rows/columns.")
  }
  
  if (!any(X == 0, na.rm = TRUE)) {
    warning("No zeros remain after deleting rows/columns; skipping zero replacement.")
  }
  
  if (nrow(dl) == 1) {
    dl <- matrix(dl, nrow = nn, ncol = D, byrow = TRUE)
  }
  
  closed <- 0
  if (all(abs(c_base - mean(c_base)) < .Machine$double.eps^0.3)) closed <- 1
  
  if (!is.null(closure)){
    if (closed == 1) {stop("closure: The data are already closed to ", c_base[1])}
    c_adj <- rep(closure, nn)
  } else {
    c_adj <- c_base
  }
  
  X_na <- is.na(X)
  
  
  est <- matrix(0, nrow = nn, ncol = D)
  
  miss_cols <- which(colSums(X_na) > 0)
  
  if (length(miss_cols) > 0) {
    gms <- apply(
      X[, miss_cols, drop = FALSE],
      2,
      function(x) gm(x[!is.na(x) & x > 0], na.rm = TRUE)
    )
    
    if (any(!is.finite(gms)) || any(gms <= 0)) {
      stop("Cannot compute finite positive geometric means for missing-value imputation")
    }
    
    est[, miss_cols] <- matrix(
      gms,
      nrow = nn,
      ncol = length(miss_cols),
      byrow = TRUE
    )
  }
  
  imputed_only <- est
  imputed_only[!X_na] <- 0
  sum_est <- rowSums(imputed_only)
  
  if (any(!is.finite(c_adj)) || any(c_adj <= 0)) {
    stop("Closure/row sums must be positive and finite")
  }
  
  problem_rows <- which(sum_est >= c_adj)
  
  if (length(problem_rows) > 0) {
    warning("Some missing-value estimates exceed or equal the row total; unadjusted geometric-mean estimates used for those rows.")
  }
  
  adjustment <- 1 - (sum_est / c_adj)
  
  if (length(problem_rows) > 0) {
    adjustment[problem_rows] <- 1
  }
  
  if (any(!is.finite(adjustment)) || any(adjustment <= 0)) {
    stop("Non-positive or non-finite adjustment factor encountered")
  }
  
  X[X_na] <- (est / adjustment)[X_na]
  
  if (closed==1){
    X <- (X / rowSums(X)) * c_base[1]
  } 
  
  
  if (any(X < 0, na.rm = TRUE)) {
    stop("multRepl: negative imputed values were generated (please check out help for advice)")
  }
  
  
  
  if (any(X == 0, na.rm = TRUE)) {
    X <- multRepl(
      X,
      label = 0,
      dl = dl,
      frac = frac,
      closure = closure,
      z.warning = z.warning,
      z.delete = FALSE
    )
  }
  
  
  
  if (!is.null(nam) && length(nam) == ncol(X)) {
    colnames(X) <- nam
  }
  
  
  return(as.data.frame(X))
  
}