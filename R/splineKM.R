splineKM <- function(x, label = NULL, dl = NULL, n.knots = NULL,
                     legend.pos = "bottomright", ylab = "ECDF",
                     xlab = "Value", col.km = "black", lty.km = 1, lwd.km = 1,
                     col.sm = "red", lty.sm = 2, lwd.sm = 2, ...) {
  
  if (is.character(dl) || is.null(dl)) stop("dl must be a numeric vector or matrix")
  if (length(dl) != length(x)) stop("x and dl must be two vectors of the same length")
  if (is.null(label)) stop("A value for label must be given")
  
  if (!is.na(label)) {
    if (label != 0 && any(x == 0, na.rm = TRUE)) stop("Zero values not labelled as censored values were found in the data")
    if (any(is.na(x))) stop("NA values not labelled as censored values were found in the data")
    who <- !is.na(x) & x == label
  } else {
    if (any(x == 0, na.rm = TRUE)) stop("Zero values not labelled as censored values were found in the data")
    if (!any(is.na(x))) stop(paste("Label", label, "was not found in the data"))
    who <- is.na(x)
  }
  
  if (!is.null(n.knots) && length(n.knots) != 1) stop("n.knots must contain a single value")
  
  x[who] <- dl[who]
  event_status <- !who
  
  valid_obs <- x[!is.na(x)]
  flip_factor <- max(valid_obs) + (diff(range(valid_obs)) / 2)
  flipped_obs <- flip_factor - x
  
  surv_obj <- survival::Surv(time = flipped_obs, event = event_status, type = "right")
  
  # CRITICAL FIX: No `...` here. survfit cannot handle graphical parameters like xlim.
  fit <- survival::survfit(surv_obj ~ 1) 
  
  fit$time <- flip_factor - fit$time
  ord <- order(fit$time)
  fit$time <- fit$time[ord]
  fit$surv <- fit$surv[ord]
  fit$n.risk <- fit$n.risk[ord]
  
  t_vals <- rev(fit$time)
  s_vals <- rev(fit$surv)
  
  if (is.null(n.knots)) {
    scdf_fit <- smooth.spline(t_vals, s_vals)
  } else {
    scdf_fit <- smooth.spline(t_vals, s_vals, nknots = n.knots)
  }
  
  scdf <- approxfun(scdf_fit$x, scdf_fit$y)
  
  # `...` is safely passed to plot(), where xlim belongs!
  plot(fit, conf.int = FALSE, ylab = ylab, xlab = xlab,
       col = col.km, lty = lty.km, lwd = lwd.km, ...)
  
  lines(t_vals, scdf(t_vals), type = "l",
        col = col.sm, lty = lty.sm, lwd = lwd.sm)
  abline(h = 1, col = "white", lwd = 4)
  legend(legend.pos, bty = "n",
         legend = c("KM estimate", "KMSS estimate"),
         lty = c(lty.km, lty.sm), col = c(col.km, col.sm), lwd = c(lwd.km, lwd.sm))
}