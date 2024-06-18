multRepl <-
  function(X,label=NULL,dl=NULL,frac=0.65,imp.missing=FALSE,closure=NULL,z.warning=0.8,z.delete=TRUE,delta=NULL){
    
    if (any(X<0, na.rm=T)) stop("X contains negative values")

    if (is.character(X)) stop("X is not a valid data matrix or vector.")
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
      if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored or missing values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    }
    if (imp.missing==FALSE){
      if (is.character(dl)) stop("dl must be a numeric vector or matrix")
      if (is.null(dl)){ # If dl not given use min per column
        dl <- apply(X,2, function(x) min(x[x!=label]))
        warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
      }
      if (is.vector(dl)) dl <- matrix(dl,nrow=1)
      dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
    }
    if (is.vector(X)){
      if (imp.missing==TRUE) stop("Data matrix required: missing values cannot be imputed in single vectors")
      if (ncol(dl)!=ncol(as.data.frame(matrix(X,ncol=length(X)),stringsAsFactors=TRUE))) stop("The number of columns in X and dl do not agree")
    }
    if (!is.vector(X)){
      if (imp.missing==FALSE){
        if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
        if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
      }
    }
    
    if (!missing("delta")){
      warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
      frac <- delta
    }
    
    gm <- function(x, na.rm=TRUE){
      exp(sum(log(x), na.rm=na.rm) / length(x[!is.na(x)]))
    }
    
    nam <- NULL
    if (!is.null(names(X))) nam <- names(X)
    if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)),stringsAsFactors=TRUE)
    
    rnames <- rownames(X)
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    if (!is.null(rownames(X))) {rownames(X) <- rnames} # Needed for the vector case
    if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)),stringsAsFactors=TRUE)
    
    if (nrow(X) > 1){
      checkNumZerosCol <- apply(X, 2, function(x) sum(is.na(x)))
      
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
      
      checkNumZerosRow <- apply(X, 1, function(x) sum(is.na(x)))  
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
    }
    
    nn <- nrow(X); D <- ncol(X)
    c <- apply(X,1,sum,na.rm=TRUE)
      
    # Check for closure
    closed <- 0
    if (all(abs(c - mean(c)) < .Machine$double.eps^0.3)) closed <- 1
    
    if (imp.missing==FALSE){
      if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl
    }
    
    Y <- X
    
    if (!is.null(closure)){
      if (closed == 1) {stop("closure: The data are already closed to ",c[1])}
      resid <- apply(X,1, function(x) closure-sum(x, na.rm = TRUE))
      Xresid <- cbind(X,resid,stringsAsFactors=TRUE)
      c <- rep(closure,nn)
      Y <- Xresid
    }
    
    if (imp.missing==FALSE){
      for (i in 1:nn){
        if (any(is.na(X[i,]))){
          z <- which(is.na(X[i,]))
          Y[i,z] <- frac*dl[i,z]
          if (!is.null(closure)){
            Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*Xresid[i,-z]
            tmp <- Y[i,-(D+1)]
            X[i,z] <- as.numeric((X[i,-z][1]/tmp[-z][1]))*Y[i,z]
          }
          else{
          Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*X[i,-z]
          X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
          }
        }
      }
    }
    
    if (imp.missing==TRUE){
      gms <- apply(X,2,gm)
      for (i in 1:nn){
        if (any(is.na(X[i,]))){
          z <- which(is.na(X[i,]))
          Y[i,z] <- gms[z]
          if (!is.null(closure)){
            Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(Xresid[i,-z]))*Xresid[i,-z]
            tmp <- Y[i,-(D+1)]
            X[i,z] <- as.numeric((X[i,-z][1]/tmp[-z][1]))*Y[i,z]
          }
          else{
          Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(X[i,-z]))*X[i,-z]
          X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
          }
        }
      }      
    }

    if (!is.null(nam)) names(X) <- nam
       
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    } 
    
    if (any(X < 0)) warning("multRepl: negative imputed values were generated (please check out help for advice)")
    
    return(as.data.frame(X,stringsAsFactors=TRUE))
  }
