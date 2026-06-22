splineKM <- function(x, label = NULL, dl = NULL, n.knots = NULL,
                     legend.pos = "bottomright", ylab = "ECDF",
                     xlab = "Value", col.km = "black", lty.km = 1, lwd.km = 1,
                     col.sm = "red", lty.sm = 2, lwd.sm = 2, ...) {
  
  if (!is.vector(x) || is.list(x)) {
    stop("x must be a numeric vector")
  }
  
  if (is.character(x)) {
    stop("x must be a numeric vector")
  }
  
  x <- as.numeric(x)
  
  if (any(x < 0, na.rm = TRUE)) {
    stop("x contains negative values")
  }
  
  if (is.null(label)) {
    stop("A value for label must be given")
  }
  
  if (is.character(dl) || is.null(dl)) {
    stop("dl must be a numeric vector")
  }
  
  dl <- as.numeric(dl)
  
  if (length(dl) != length(x)) {
    stop("x and dl must be two vectors of the same length")
  }
  
  if (any(dl < 0, na.rm = TRUE) || any(!is.finite(dl), na.rm = TRUE)) {
    stop("Detection limits in dl must be non-negative and finite")
  }
  
  if (!is.na(label)) {
    if (!any(x == label, na.rm = TRUE)) {
      stop(paste("Label", label, "was not found in the data"))
    }
    
    if (label != 0 && any(x == 0, na.rm = TRUE)) {
      stop("Zero values not labelled as censored values were found in the data")
    }
    
    if (any(is.na(x))) {
      stop("NA values not labelled as censored values were found in the data")
    }
    
    x[x == label] <- NA
  } else {
    if (any(x == 0, na.rm = TRUE)) {
      stop("Zero values not labelled as censored values were found in the data")
    }
    
    if (!any(is.na(x))) {
      stop(paste("Label", label, "was not found in the data"))
    }
  }
  
  if (!is.null(n.knots)) {
    if (!is.numeric(n.knots) || length(n.knots) != 1 || is.na(n.knots) || n.knots < 1) {
      stop("n.knots must be a single positive numeric value")
    }
  }
  
  who <- is.na(x)
  w <- which(who)
  
  if (length(w) == 0) {
    stop("No censored values found")
  }
  
  if (any(dl[who] <= 0, na.rm = TRUE)) {
    stop("Detection limits for censored values must be strictly positive")
  }
  
  xcen <- who
  x[who] <- dl[who]
  
  if (sum(!xcen) < 2) {
    stop("Kaplan-Meier smoothing requires at least two observed values")
  }
  
  km.ecdf <- cenfit(x, xcen)
  
  x.km <- rev(km.ecdf@survfit$time)
  y.km <- rev(km.ecdf@survfit$surv)
  
  if (length(unique(x.km)) < 2 || length(unique(y.km)) < 2) {
    stop("Kaplan-Meier curve has too few unique points for spline smoothing")
  }
  
  if (is.null(n.knots)) {
    scdf <- smooth.spline(x.km, y.km)
  } else {
    scdf <- smooth.spline(x.km, y.km, nknots = n.knots)
  }
  
  scdf.fun <- approxfun(scdf$x, scdf$y)
  
  plot(km.ecdf@survfit,
       conf.int = FALSE,
       ylab = ylab,
       xlab = xlab,
       col = col.km,
       lty = lty.km,
       lwd = lwd.km,
       ...)
  
  lines(x.km, scdf.fun(x.km),
        type = "l",
        col = col.sm,
        lty = lty.sm,
        lwd = lwd.sm)
  
  abline(h = 1, col = "white", lwd = 4)
  
  legend(legend.pos,
         bty = "n",
         legend = c("KM estimate", "KMSS estimate"),
         lty = c(lty.km, lty.sm),
         col = c(col.km, col.sm),
         lwd = c(lwd.km, lwd.sm))
}