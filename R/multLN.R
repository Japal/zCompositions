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
    
    # Standalone version of cenmle function
    # This function is a recreation of `cenmle` from the R package `NADA`
    # by D.R. Helsel and L.H. Lee. All credit for the methodology goes to the original authors
    
    cenmle_standalone <- function(obs, censored, dist = "normal") {
    
      # Internal Log-Likelihood Function (for Normal Distribution) ---

      llike.norm <- function(params, x, cen) {
        mean.val <- params[1]
        sd.val <- params[2]
        
        # Return a very small number if sd is non-positive to avoid errors
        if (sd.val <= 0) {
          return(-1e9)
        }
        
        uncensored_obs <- x[!cen]
        censored_obs <- x[cen]
        
        # Log-likelihood for the uncensored values
        ll_uncensored <- sum(dnorm(uncensored_obs, mean = mean.val, sd = sd.val, log = TRUE))
        
        # Log-likelihood for the censored values
        # This is the log of the cumulative distribution function (pnorm)
        ll_censored <- sum(pnorm(censored_obs, mean = mean.val, sd = sd.val, log.p = TRUE))
        
        # Total log-likelihood is the sum of the two parts
        total_ll <- ll_uncensored + ll_censored
        return(total_ll)
      }

      if (dist == "lognormal") {
        data_for_optim <- log(obs)
      } else {
        data_for_optim <- obs
      }
      
      # Parameter Estimation 
      start_mean <- mean(data_for_optim)
      start_sd <- sd(data_for_optim)
      
      if(start_sd == 0){
        start_sd <- 1
      }
      
      # Maximize the log-likelihood.
      # `control = list(fnscale = -1)` tells optim to maximize instead of minimize.
      mle_results <- optim(
        par = c(start_mean, start_sd),
        fn = llike.norm,
        x = data_for_optim,
        cen = censored,
        control = list(fnscale = -1),
        hessian = FALSE # Hessian matrix not needed for this implementation
      )
      
      mu_hat <- mle_results$par[1]
      sigma_hat <- mle_results$par[2]
      
      if (dist == "lognormal") {
        # For lognormal, back-transform the estimated mean (mu) and sd (sigma)
        # from the log-scale to the original data scale.
        final_mean <- exp(mu_hat + 0.5 * sigma_hat^2)
        final_sd <- sqrt(exp(2 * mu_hat + sigma_hat^2) * (exp(sigma_hat^2) - 1))
        param_names <- c("Mean (lognormal)", "SD (lognormal)")
      } else {
        # For normal, the results are already on the correct scale.
        final_mean <- mu_hat
        final_sd <- sigma_hat
        param_names <- c("Mean (normal)", "SD (normal)")
      }
      
      output <- list(
        parameters = c(final_mean, final_sd)
      )
      
      return(output)
    }
    
    # Standalone version of cenros function
    # This function is a recreation of `cenros` from the R package `NADA`
    # by D.R. Helsel and L.H. Lee. All credit for the methodology goes to the original authors
    
    # Helper functions
    hc.ppoints <- function(obs, censored) {
      if (length(obs) != length(censored))
        stop("obs and censored must be the same length.")
      
      n <- length(obs)
      n.cen <- sum(censored)
      n.uncen <- n - n.cen
      
      if (n.uncen == 0)
        stop("All data are censored. Cannot compute plotting positions.")
      
      # Order the data
      sort_order <- order(obs)
      obs_sorted <- obs[sort_order]
      censored_sorted <- censored[sort_order]
    
      pp <- numeric(n)
      last_pp <- 0
      last_val <- -Inf
      
      for (i in 1:n) {
        if (!censored_sorted[i]) {
          # For uncensored values
          val <- obs_sorted[i]
          if (val > last_val) {
            # New uncensored value
            A <- sum(obs_sorted < val)
            N <- sum(obs_sorted == val & !censored_sorted)
            C <- (A + 1 - 0.5) / n
            D <- (A + N - 0.5) / n
            pp[i] <- (C + D) / 2
          } else {
            # Tied with previous uncensored value
            pp[i] <- last_pp
          }
          last_pp <- pp[i]
          last_val <- val
        }
      }
      
      # For censored values, we estimate their position
      for (i in 1:n) {
        if (censored_sorted[i]) {
          val <- obs_sorted[i]
          A <- sum(obs_sorted < val)
          N <- sum(obs_sorted == val & !censored_sorted)
          if (N > 0) {
            # If censored value is tied with uncensored values
            pp[i] <- (A + 1 - 0.5) / n
          } else {
            # If censored value is not tied
            pp[i] <- (A - 0.5) / n
          }
        }
      }
      
      # Un-sort the plotting positions to match original data order
      pp[sort_order] <- pp
      
      # Handle potential issues with plotting positions for values > 1 or < 0
      # This can happen with tied values at the upper/lower boundaries.
      pp[pp >= 1] <- 1 - .Machine$double.eps
      pp[pp <= 0] <- .Machine$double.eps
      
      return(pp)
    }
    trueT <- function(x) {
      return(x)
    }
    
    # Core cenros function
    cenros_standalone <- function(obs, censored, forwardT = "log", reverseT = "exp") {
    
      # Get transformation functions from their string names ---
      # If NULL, use the identity function `trueT`
      if (is.null(forwardT) || is.null(reverseT)) {
        forwardT <- "trueT"
        reverseT <- "trueT"
      }
      forwardT_func <- get(forwardT)
      reverseT_func <- get(reverseT)
      
      # ROS Algorithm
      
      pp <- hc.ppoints(obs, censored)

      uncensored_obs <- obs[!censored]
      uncensored_pp <- pp[!censored]

      obs_transformed <- forwardT_func(uncensored_obs)
      pp_quantiles <- qnorm(uncensored_pp)
      
      if (any(is.infinite(obs_transformed)) || any(is.infinite(pp_quantiles))) {
        stop("Infinite values produced during transformation. Check data or transformation functions.")
      }
      
      ros_lm <- lm(obs_transformed ~ pp_quantiles)
      
      censored_pp <- pp[censored]
      censored_quantiles <- qnorm(censored_pp)
      
      new_data <- data.frame(pp_quantiles = censored_quantiles)
      predicted_transformed_values <- predict(ros_lm, newdata = new_data)
      
      modeled_censored_values <- reverseT_func(predicted_transformed_values)
      
      final_modeled_data <- numeric(length(obs))
      final_modeled_data[!censored] <- obs[!censored] # Use original uncensored values
      final_modeled_data[censored] <- modeled_censored_values # Use modeled censored values
      
      result <- list(
        modeled = final_modeled_data
      )

      return(result)
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
    
    if (random==FALSE){
      
      cenGeoMean <- function(x,dl,...){ 
        
        xcen <- ifelse(is.na(x),TRUE,FALSE)
        x[is.na(x)] <- dl[is.na(x)]
        
        if (rob) {ymean <- mean(cenros_standalone(log(x),xcen)$modeled);
                    ysd <- sd(cenros_standalone(log(x),xcen)$modeled)} 
        else 
        {ymean <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[1];
           ysd <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[2]}
        
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
        
        if (rob) {ymean <- mean(cenros_standalone(log(x),xcen)$modeled);
                  ysd <- sd(cenros_standalone(log(x),xcen)$modeled)} 
        else 
        {ymean <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[1];
         ysd <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[2]}
        
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
