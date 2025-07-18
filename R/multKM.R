multKM <-
  function (X,label=NULL,dl=NULL,n.draws=1000,n.knots=NULL,z.warning=0.8,z.delete=TRUE)
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
    
    if ((!is.null(n.knots)) & (length(n.knots)!=1) & (length(n.knots)!=ncol(X))) stop("The dimensions of n.knots and X do not agree")
    if ((!is.null(n.knots)) & (length(n.knots)==1)) {n.knots <- rep(list(n.knots),ncol(X))}
    
    # Standalone Replication of NADA::cenfit
    #
    # Computes an estimate of an empirical cumulative distribution function (ECDF)
    # for left-censored data using the Kaplan-Meier method, by "flipping" the
    # data to be compatible with the survival::survfit function.
    
    cenfit_standalone <- function(obs, censored,...) {
      
      flip_factor <- max(obs, na.rm = TRUE) + (diff(range(obs, na.rm = TRUE)) / 2)
      flipped_obs <- flip_factor - obs
      event_status <-!censored
      surv_obj <- survival::Surv(time = flipped_obs, event = event_status, type = "right")
      fit <- survival::survfit(surv_obj ~ 1,...)
      fit$time <- flip_factor - fit$time
      ord <- order(fit$time)
      fit$time <- fit$time[ord]; fit$surv <- fit$surv[ord]; fit$n.risk <- fit$n.risk[ord]
      return(fit)
      
    }
    
    km.imp <- function(x,dl,...){
      
      who <- is.na(x); w <- which(who)
      
      xcen <- ifelse(who,TRUE,FALSE)
      x[who] <- dl[who]
      
      km.ecdf <- cenfit_standalone(x,xcen)
      x.km <- rev(km.ecdf$time) 
      y.km <- rev(km.ecdf$surv)
      if (is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km)}
      if (!is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km,nknots=n.knots.part)}
      scdf.fun <- approxfun(scdf$x,scdf$y)
      inv.scdf <- approxfun(scdf$y,scdf$x)
      
      for (i in 1:length(w)){
        if (dl[w[i]] > min(x[!who])){
          temp <- inv.scdf(runif(n.draws,0,scdf.fun(dl[w[i]])))
          x[w[i]] <- exp(mean(log(temp),na.rm=T))
        }
      }
      return(as.numeric(x))
    }
    
    rnames <- rownames(X)
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    rownames(X) <- rnames
    
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
    
    for (part in 1:p)
    {
      if (any(is.na(X[,part]))) 
      {
        n.knots.part <- n.knots[[part]]
        est[,part] <- km.imp(X[,part],dl[,part],n.draws,n.knots.part)
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
  
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    }
    
    return(as.data.frame(X,stringsAsFactors=TRUE))
  }  
    