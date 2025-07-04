splineKM <- function(x,label=NULL,dl=NULL,n.knots=NULL,
                     legend.pos="bottomright",
                     ylab="ECDF",
                     xlab="Value",
                     col.km="black",lty.km=1,lwd.km=1,
                     col.sm="red",lty.sm=2,lwd.sm=2,...){
  
  if (is.character(dl) || is.null(dl)) stop("dl must be a numeric vector or matrix")
  if (length(dl)!=length(x)) stop("x and dl must be two vectors of the same length")
  
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (label!=0 & any(x==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data")
    if (any(is.na(x))) stop(paste("NA values not labelled as censored values were found in the data"))
  }
  if (is.na(label)){
    if (any(x==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data")
    if (!any(is.na(x),na.rm=T)) stop(paste("Label",label,"was not found in the data"))
  }
  
  if ((!is.null(n.knots)) & (length(n.knots)!=1)) stop("n.knots must contain a single value")
  
  # Standalone Replication of NADA::cenfit (same as cenfit_standalone in multKM but issue if same name)
  #
  # Computes an estimate of an empirical cumulative distribution function (ECDF)
  # for left-censored data using the Kaplan-Meier method, by "flipping" the
  # data to be compatible with the survival::survfit function.
  
  cenfit_standalone2 <- function(obs, censored,...) {
    
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
  
  x[x==label] <- NA
  
  who <- is.na(x); w <- which(who)
  
  xcen <- ifelse(who,TRUE,FALSE)
  x[who] <- dl[who]
  
  dat <- data.frame(x,xcen,stringsAsFactors=TRUE)
  km.ecdf <- cenfit_standalone2(dat$x,dat$xcen)

  x <- rev(km.ecdf$time) 
  y <- rev(km.ecdf$surv)
  
  if (is.null(n.knots)) {scdf <- smooth.spline(x,y)}
  if (!is.null(n.knots)) {scdf <- smooth.spline(x,y,nknots=n.knots)}
  scdf <- approxfun(scdf$x,scdf$y)
  
  plot(km.ecdf,conf.int=FALSE,ylab=ylab,xlab=xlab,
       col=col.km,lty=lty.km,lwd=lwd.km, ...)
  lines(x,scdf(x),type="l",
       col=col.sm,lty=lty.sm,lwd=lwd.sm)
  abline(h=1,col="white",lwd=4)
  legend(legend.pos,bty="n",
         legend=c("KM estimate","KMSS estimate"),
         lty=c(lty.km,lty.sm),col=c(col.km,col.sm),lwd=c(lwd.km,lwd.sm)) 

}