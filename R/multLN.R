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
    if (is.null(dl)){ 
      dl <- apply(X, 2, function(x) min(x[!(x %in% label)], na.rm = TRUE))
      warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
    }
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    dl <- as.matrix(dl) 
    if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
    if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
    
    # Standalone version of cenmle function
    cenmle_standalone <- function(obs, censored, dist = "normal") {
      
      llike.norm <- function(params, x, cen) {
        mean.val <- params[1]
        sd.val <- params[2]
        
        if (sd.val <= 0) {
          return(-1e9)
        }
        
        uncensored_obs <- x[!cen]
        censored_obs <- x[cen]
        
        ll_uncensored <- sum(dnorm(uncensored_obs, mean = mean.val, sd = sd.val, log = TRUE))
        ll_censored <- sum(pnorm(censored_obs, mean = mean.val, sd = sd.val, log.p = TRUE))
        
        total_ll <- ll_uncensored + ll_censored
        return(total_ll)
      }
      
      if (dist == "lognormal") {
        data_for_optim <- log(obs)
      } else {
        data_for_optim <- obs
      }
      
      start_mean <- mean(data_for_optim)
      start_sd <- sd(data_for_optim)
      
      if(start_sd == 0){
        start_sd <- 1
      }
      
      mle_results <- optim(
        par = c(start_mean, start_sd),
        fn = llike.norm,
        x = data_for_optim,
        cen = censored,
        control = list(fnscale = -1),
        hessian = FALSE 
      )
      
      mu_hat <- mle_results$par[1]
      sigma_hat <- mle_results$par[2]
      
      if (dist == "lognormal") {
        final_mean <- exp(mu_hat + 0.5 * sigma_hat^2)
        final_sd <- sqrt(exp(2 * mu_hat + sigma_hat^2) * (exp(sigma_hat^2) - 1))
        param_names <- c("Mean (lognormal)", "SD (lognormal)")
      } else {
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
    hc.ppoints <- function(obs, censored) {
      if (length(obs) != length(censored))
        stop("obs and censored must be the same length.")
      
      n <- length(obs)
      n.cen <- sum(censored)
      n.uncen <- n - n.cen
      
      if (n.uncen == 0)
        stop("All data are censored. Cannot compute plotting positions.")
      
      sort_order <- order(obs)
      obs_sorted <- obs[sort_order]
      censored_sorted <- censored[sort_order]
      
      pp <- numeric(n)
      last_pp <- 0
      last_val <- -Inf
      
      for (i in 1:n) {
        if (!censored_sorted[i]) {
          val <- obs_sorted[i]
          if (val > last_val) {
            A <- sum(obs_sorted < val)
            N <- sum(obs_sorted == val & !censored_sorted)
            C <- (A + 1 - 0.5) / n
            D <- (A + N - 0.5) / n
            pp[i] <- (C + D) / 2
          } else {
            pp[i] <- last_pp
          }
          last_pp <- pp[i]
          last_val <- val
        }
      }
      
      for (i in 1:n) {
        if (censored_sorted[i]) {
          val <- obs_sorted[i]
          A <- sum(obs_sorted < val)
          N <- sum(obs_sorted == val & !censored_sorted)
          if (N > 0) {
            pp[i] <- (A + 1 - 0.5) / n
          } else {
            pp[i] <- (A - 0.5) / n
          }
        }
      }
      
      pp[sort_order] <- pp
      
      pp[pp >= 1] <- 1 - .Machine$double.eps
      pp[pp <= 0] <- .Machine$double.eps
      
      return(pp)
    }
    
    trueT <- function(x) {
      return(x)
    }
    
    cenros_standalone <- function(obs, censored, forwardT = "log", reverseT = "exp") {
      
      if (is.null(forwardT) || is.null(reverseT)) {
        forwardT <- "trueT"
        reverseT <- "trueT"
      }
      forwardT_func <- get(forwardT)
      reverseT_func <- get(reverseT)
      
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
      final_modeled_data[!censored] <- obs[!censored] 
      final_modeled_data[censored] <- modeled_censored_values 
      
      result <- list(
        modeled = final_modeled_data
      )
      
      return(result)
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
    
    nn <- nrow(X); p <- ncol(X)
    c <- rowSums(X, na.rm=TRUE)
    
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
    
    if (nrow(dl)==1){
      dl <- matrix(dl, nrow=nn, ncol=p, byrow=TRUE)
      est <- dl
    }
    else est <- dl
    
    X_na <- is.na(X)
    
    if (random==FALSE){
      
      cenGeoMean <- function(x,dl,...){ 
        
        xcen <- is.na(x)
        x[xcen] <- dl[xcen]
        
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
        if (any(X_na[,part])) 
        {
          est[,part] <- cenGeoMean(X[,part],dl[,part],rob)
        }
        else {est[,part] <- 0}
      }
      
    } else { 
      
      meanln <- rep(0,p); sdln <- rep(0,p)
      
      for (j in 1:p){
        x <- X[,j]
        xcen <- X_na[,j]
        
        if (any(xcen)) {
          x[xcen] <- dl[xcen,j]
          
          if (rob) {ymean <- mean(cenros_standalone(log(x),xcen)$modeled);
          ysd <- sd(cenros_standalone(log(x),xcen)$modeled)} 
          else 
          {ymean <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[1];
          ysd <- cenmle_standalone(log(x),xcen,dist="normal")$parameters[2]}
          
          meanln[j] <- ymean
          sdln[j] <- ysd
          
          # Vectorized generation of random estimates for missing values in this column
          est[xcen, j] <- exp(rtruncnorm(sum(xcen), -Inf, log(dl[xcen, j]), meanln[j], sdln[j]))
        } else {
          est[, j] <- 0
        }
      }
    } 
    
    # Unified Vectorized Closure Adjustment 
    imputed_only <- est
    imputed_only[!X_na] <- 0
    sum_est <- rowSums(imputed_only)
    
    adjustment <- 1 / (1 - (sum_est / c))
    
    X[X_na] <- (est * adjustment)[X_na]
    
    if (closed==1){
      X <- (X / rowSums(X)) * c[1]
    }
    
    return(as.data.frame(X,stringsAsFactors=TRUE))
  }
