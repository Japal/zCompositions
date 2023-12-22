multLN <-
  function (X,label=NULL,dl=NULL,rob=FALSE,random=FALSE,z.warning=0.8,z.delete=TRUE)
  {
    
    if (any(X<0, na.rm=T)) stop("X contains negative values")
    if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
      if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
      if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    }
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.null(dl)){ # If dl not given use min per column
      dl <- apply(X,2, function(x) min(x[x!=label]))
      warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
    }
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
    if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
    if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
    
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    
    checkNumZerosCol <- apply(X,2,function(x) sum(is.na(x)))
    if (any(checkNumZerosCol/nrow(X) >= z.warning)) {
      cases <- which(checkNumZerosCol/nrow(X) >= z.warning)
      if (z.delete == TRUE){
        if (length(cases) > (ncol(X)-2)) {stop(paste("More than 2 columns contain >",z.warning*100,"% zeros/unobserved values (see arguments z.warning and z.delete).",sep=""))}
        X <- X[,-cases]
        action <- "deleted"
        warning(paste("Column no. ",cases," containing >",z.warning*100,"% zeros/unobserved values ",action," (see arguments z.warning and z.delete).\n",sep=""))
      }
      else{
        action <- "found"
        stop(paste("Column no. ",cases," containing >",z.warning*100,"% zeros/unobserved values ",action," (see arguments z.warning and z.delete. Check out with zPatterns()).\n",sep=""))
      }
    }
    
    checkNumZerosRow <- apply(X,1,function(x) sum(is.na(x)))
    if (any(checkNumZerosRow/ncol(X) >= z.warning)) {
      cases <- which(checkNumZerosRow/ncol(X) >= z.warning)
      if (z.delete == TRUE){
        if (length(cases) > (nrow(X)-2)) {stop(paste("More than 2 rows contain >",z.warning*100,"% zeros/unobserved values (see arguments z.warning and z.delete).",sep=""))}
        X <- X[,-cases]
        action <- "deleted"
        warning(paste("Column no. ",cases," containing >",z.warning*100,"% zeros/unobserved values ",action," (see arguments z.warning and z.delete).\n",sep=""))
      }
      else{
        action <- "found"
        stop(paste("Column no. ",cases," containing >",z.warning*100,"% zeros/unobserved values ",action," (see arguments z.warning and z.delete. Check out with zPatterns()).\n",sep=""))
      }
    }
    
    nn <- nrow(X); p <- ncol(X)
    c <- apply(X,1,sum,na.rm=TRUE)
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
    
    if (nrow(dl)==1){
      dl <- matrix(rep(1,nn),ncol=1)%*%dl
      est <- dl
    }
    else est <- dl
    
    if (random==FALSE){
      
      cenGeoMean <- function(x,dl,...){ 
        
        xcen <- ifelse(is.na(x),TRUE,FALSE)
        x[is.na(x)] <- dl[is.na(x)]
        
        if (rob) {ymean <- summary(cenros(x,xcen))$coefficients[1];
                  ysd <- summary(cenros(x,xcen))$coefficients[2]} 
        else 
        {ymean <- mean(suppressWarnings(cenmle(log(x),xcen,dist="gaussian")))[1];
         ysd <- sd(suppressWarnings(cenmle(log(x),xcen,dist="gaussian")))[1]}
        
        fdl <- dnorm((log(dl)-ymean)/ysd, mean = 0, sd = 1, log = FALSE)
        Pdl <- pnorm((log(dl)-ymean)/ysd, mean = 0, sd = 1, log.p = FALSE)
        gmeancen <- exp(ymean-ysd*(fdl/Pdl))
        
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
          Y[i,z] <- est[i,z]
          Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*X[i,-z]
          X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
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
            Y[i,z[j]] <- exp(rtruncnorm(1,-Inf,log(dl[i,z[j]]),meanln[z[j]],sdln[z[j]]))
          }
          Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*X[i,-z]
          X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
        }
      }  
    } # End if random
    
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    }
    
    return(as.data.frame(X,stringsAsFactors=TRUE))
  }
