multLN <-
  function (X,label=NULL,dl=NULL,rob=FALSE,random=FALSE,z.warning=0.8,z.delete=TRUE)
  {
    
    
    if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")
    
    if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
      stop("X must be a data matrix with at least two rows and two columns")
    }
    
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=TRUE)) stop(paste("Label",label,"was not found in the data set"))
      if (label != 0 && any(X == 0, na.rm = TRUE)) {
        stop("Zero values not labelled as censored values were found in the data set")
      }
      
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
    
    
    if (!is.numeric(z.warning) || length(z.warning) != 1 || is.na(z.warning) ||
        z.warning < 0 || z.warning > 1) {
      stop("z.warning must be a single numeric value between 0 and 1")
    }
    
    if (!is.logical(z.delete) || length(z.delete) != 1 || is.na(z.delete)) {
      stop("z.delete must be TRUE or FALSE")
    }
    
    
    rnames <- rownames(X)
    if (!is.na(label)) X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    rownames(X) <- rnames
    
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
        dl <- dl[, -cases, drop = FALSE]
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
    
    nn <- nrow(X); p <- ncol(X)

    if (nn < 2 || p < 2) {
      stop("After deleting rows/columns, X must still have at least two rows and two columns")
    }
    
    c <- rowSums(X, na.rm=TRUE)
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
    
    if (nrow(dl)==1){
      dl <- dl[rep(1, nn), , drop = FALSE]
      est <- dl
    }
    else est <- dl
    
    if (random==FALSE){
      
      cenGeoMean <- function(x, dl, ...) {
        
        xcen <- is.na(x)
        
        if (any(dl[xcen] <= 0, na.rm = TRUE) || any(!is.finite(dl[xcen]), na.rm = TRUE)) {
          stop("Detection limits for censored values must be strictly positive and finite")
        }
        
        x[xcen] <- dl[xcen]
        
        if (any(x <= 0, na.rm = TRUE)) {
          stop("Positive values are required for lognormal censored imputation")
        }
        
        if (sum(!xcen) < 2) {
          stop("Lognormal censored imputation requires at least two observed values per part")
        }
        
        if (rob) {
          fit <- summary(cenros(x, xcen))$coefficients
          ymean <- fit[1]
          ysd <- fit[2]
        } else {
          fit <- suppressWarnings(cenmle(log(x), xcen, dist = "gaussian"))
          ymean <- mean(fit)[1]
          ysd <- sd(fit)[1]
        }
        
        if (!is.finite(ymean) || !is.finite(ysd) || ysd <= 0) {
          stop("Non-finite or non-positive lognormal standard deviation encountered")
        }
        
        gmeancen <- numeric(length(x))
        
        z <- (log(dl[xcen]) - ymean) / ysd
        fdl <- dnorm(z)
        Pdl <- pnorm(z)
        
        ratio <- fdl / pmax(Pdl, 1e-12)
        gmeancen[xcen] <- exp(ymean - ysd * ratio)
        
        return(as.numeric(gmeancen))
      }
      
      for (part in 1:p)
      {
        if (any(is.na(X[,part]))) 
        {
          est[,part] <- cenGeoMean(X[,part],dl[,part],rob)
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
          
          Y[i, z] <- est[i, z]
          
          if (sum(Y[i, z], na.rm = TRUE) >= c[i]) {
            stop("Estimated censored values exceed or equal the row total; multiplicative adjustment cannot proceed")
          }
          
          Y[i, observed] <- (1 - sum(Y[i, z]) / c[i]) * X[i, observed]
          
          if (!is.finite(Y[i, observed][1]) || Y[i, observed][1] <= 0) {
            stop("Non-positive scaling denominator encountered in multiplicative adjustment")
          }
          
          X[i, z] <- as.numeric(X[i, observed][1] / Y[i, observed][1]) * Y[i, z]
        }
      }  
    } # End if not random
    
    else{ # If random
      
      meanln <- rep(0,p); sdln <- rep(0,p)
      
      for (j in 1:p){
        x <- X[,j]
        xcen <- ifelse(is.na(X[,j]),TRUE,FALSE)
        x[is.na(X[,j])] <- dl[is.na(X[,j]),j]
        
        if (rob) {ymean <- summary(cenros(x,xcen))$coefficients[1];
        ysd <- summary(cenros(x,xcen))$coefficients[2]} 
        else 
        {ymean <- mean(suppressWarnings(cenmle(log(x),xcen,dist="gaussian")))[1];
        ysd <- sd(suppressWarnings(cenmle(log(x),xcen,dist="gaussian")))[1]}
        
        meanln[j] <- ymean
        sdln[j] <- ysd
      }
      
      Y <- X
      
      for (i in 1:nn){
        if (any(is.na(X[i,]))){
          z <- which(is.na(X[i,]))
          for (j in 1:length(z)){
            
            if (!is.finite(dl[i, z[j]]) || dl[i, z[j]] <= 0) {
              stop("Detection limits must be strictly positive for random lognormal imputation")
            }
            
            draw <- rtruncnorm(
              1,
              -Inf,
              log(dl[i, z[j]]),
              meanln[z[j]],
              sdln[z[j]]
            )
            
            if (!is.finite(draw)) {
              stop("Random lognormal imputation produced a non-finite draw")
            }
            
            Y[i, z[j]] <- exp(draw)
            
          }
          observed <- setdiff(seq_len(p), z)
          
          if (length(observed) == 0) {
            stop("At least one row has no observed components")
          }
          
          if (!is.finite(c[i]) || c[i] <= 0) {
            stop("At least one row has zero or non-finite observed sum")
          }
          
          if (sum(Y[i, z], na.rm = TRUE) >= c[i]) {
            stop("Estimated censored values exceed or equal the row total; multiplicative adjustment cannot proceed")
          }
          
          Y[i, observed] <- (1 - sum(Y[i, z]) / c[i]) * X[i, observed]
          
          if (!is.finite(Y[i, observed][1]) || Y[i, observed][1] <= 0) {
            stop("Non-positive scaling denominator encountered in multiplicative adjustment")
          }
          
          X[i, z] <- as.numeric(X[i, observed][1] / Y[i, observed][1]) * Y[i, z]
        }
      }  
    } # End if random
    
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    }
    
    return(as.data.frame(X))
  }