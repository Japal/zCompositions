multKM <-
  function (X,label=NULL,dl=NULL,n.draws=1000,n.knots=NULL,z.warning=0.8,z.delete=TRUE)
  {
    
    if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")
    
    if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
      stop("X must be a data matrix with at least two rows and two columns")
    }
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
      if (label!=0 & any(X==0,na.rm=TRUE)) stop("Zero values not labelled as censored values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X==0,na.rm=TRUE)) stop("Zero values not labelled as censored values were found in the data set")
      if (!any(is.na(X),na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
    }
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
    
    if (any(dl < 0, na.rm = TRUE) || any(!is.finite(dl), na.rm = TRUE)) {
      stop("Detection limits in dl must be strictly positive and finite")
    }
    
    if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
    
    if (nrow(dl) > 1 && nrow(dl) != nrow(X)) {
      stop("The number of rows in X and dl do not agree")
    }
    
    if (is.null(n.knots)) {
      n.knots <- rep(list(NULL), ncol(X))
    } else if (length(n.knots) == 1) {
      n.knots <- rep(list(n.knots), ncol(X))
    } else if (length(n.knots) == ncol(X)) {
      n.knots <- as.list(n.knots)
    } else {
      stop("The dimensions of n.knots and X do not agree")
    }
    
    if (!is.numeric(n.draws) || length(n.draws) != 1 || is.na(n.draws) ||
        n.draws < 1 || n.draws != as.integer(n.draws)) {
      stop("n.draws must be a positive integer")
    }
    
    if (!is.numeric(z.warning) || length(z.warning) != 1 || is.na(z.warning) ||
        z.warning < 0 || z.warning > 1) {
      stop("z.warning must be a single numeric value between 0 and 1")
    }
    
    if (!is.logical(z.delete) || length(z.delete) != 1 || is.na(z.delete)) {
      stop("z.delete must be TRUE or FALSE")
    }
    
    km.imp <- function(x,dl,n.knots.part,n.draws,...){
      
      who <- is.na(x); w <- which(who)
      
      xcen <- ifelse(who,TRUE,FALSE)
      x[who] <- dl[who]
      
      if (sum(!who) < 1) {
        x[who] <- dl[who]
        return(as.numeric(x))
      }
      
      km.ecdf <- cenfit(x,xcen)
      x.km <- rev(km.ecdf@survfit$time) 
      y.km <- rev(km.ecdf@survfit$surv)
      if (is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km)}
      if (!is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km,nknots=n.knots.part)}
      scdf.fun <- approxfun(scdf$x,scdf$y)
      inv.scdf <- approxfun(scdf$y,scdf$x)
      
      for (i in 1:length(w)){
        if (dl[w[i]] > min(x[!who])){
    
          temp <- inv.scdf(runif(n.draws, 0, scdf.fun(dl[w[i]])))
          temp <- temp[is.finite(temp) & temp > 0]
          
          if (length(temp) == 0) {
            temp <- inv.scdf(runif(max(n.draws, 100), 0, scdf.fun(dl[w[i]])))
            temp <- temp[is.finite(temp) & temp > 0]
          }
          
          if (length(temp) == 0) {
            x[w[i]] <- dl[w[i]]
          } else {
            x[w[i]] <- exp(mean(log(temp), na.rm = TRUE))
          }
        }
      }
      return(as.numeric(x))
    }
    
    rnames <- rownames(X)
    if (!is.na(label)) {X[X==label] <- NA}
    X <- apply(X,2,as.numeric)
    rownames(X) <- rnames
    
    checkNumZerosCol <- colSums(is.na(X))
    
    if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
      cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
      if (z.delete == TRUE) {
        if (length(cases) > (ncol(X)-2)) {
          stop(paste("Almost all columns contain >", z.warning*100,
                     "% zeros/unobserved values (see arguments z.warning and z.delete).",
                     sep="", collapse = ", "))
        }      
        
        X <- X[, -cases, drop = FALSE]
        dl <- dl[, -cases, drop = FALSE]
        n.knots <- n.knots[-cases]
        action <- "deleted"
        
        warning(paste("Column no. ", paste(cases, collapse = ", "),
                      " containing >", z.warning*100,
                      "% zeros/unobserved values ", action,
                      " (see arguments z.warning and z.delete).\n",
                      sep = ""))
        
      } else {      
        action <- "found"      
        warning(paste("Column no. ", paste(cases, collapse = ", "),
                      " containing >", z.warning*100,
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
                     sep="", collapse = ", "))
        }
        
        X <- X[-cases, , drop = FALSE]
        
        if (nrow(dl) > 1) {
          dl <- dl[-cases, , drop = FALSE]
        }
        
        action <- "deleted"
        
        warning(paste("Row no. ", paste(cases, collapse = ", "),
                      " containing >", z.warning*100,
                      "% zeros/unobserved values ", action,
                      " (see arguments z.warning and z.delete).\n",
                      sep = ""))
      } else {      
        action <- "found"      
        warning(paste("Row no. ", paste(cases, collapse = ", "),
                      " containing >", z.warning*100,
                      "% zeros/unobserved values ", action,
                      " (see arguments z.warning and z.delete).\n",
                      sep = ""))    
      }
    }
    
    nn <- nrow(X); p <- ncol(X)
    c <- rowSums(X, na.rm = TRUE)
    
    if (nn < 2 || p < 2) {
      stop("After deleting rows/columns, X must still have at least two rows and two columns")
    }
    
    if (length(n.knots) != p) {
      stop("The dimensions of n.knots and X do not agree after deleting columns")
    }
  
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
    
    if (nrow(dl)==1){
      dl <- dl[rep(1, nn), , drop = FALSE]
      est <- dl
    }
    else est <- dl
    
    for (part in 1:p)
    {
      if (any(is.na(X[,part]))) 
      {
        n.knots.part <- n.knots[[part]]
        est[,part] <- km.imp(X[,part],dl[,part],n.knots.part,n.draws)
      }
      else {est[,part] <- 0}
    }
    
    Y <- X
    
    for (i in 1:nn){
      if (any(is.na(X[i,]))){
        z <- which(is.na(X[i,]))
        observed <- setdiff(seq_len(p), z)
        
        if (length(observed) == 0) {
          stop("At least one row has no observed components")
        }
        
        if (!is.finite(c[i]) || c[i] <= 0) {
          stop("At least one row has zero or non-finite observed sum")
        }
        
        if (sum(est[i, z], na.rm = TRUE) >= c[i]) {
          stop("Estimated censored values exceed or equal the row total; multiplicative adjustment cannot proceed")
        }
        
        Y[i, z] <- est[i, z]
        Y[i, observed] <- (1 - sum(Y[i, z]) / c[i]) * X[i, observed]
        
        if (!is.finite(Y[i, observed][1]) || Y[i, observed][1] <= 0) {
          stop("Non-positive scaling denominator encountered in multiplicative adjustment")
        }
        
        X[i, z] <- as.numeric(X[i, observed][1] / Y[i, observed][1]) * Y[i, z]
      }
    } 
    
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    }
    
    return(as.data.frame(X))
  }  
