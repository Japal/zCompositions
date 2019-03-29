\name{splineKM}
\alias{splineKM}
\title{
Display Kaplan-Meier empirical cumulative distribution function and smoothing spline curve fit
}
\description{
This function shows the empirical cumulative distribution function (ECDF) for left-censored data as estimated by the Kaplan-Meier (KM) method and a cubic smoothing spline fitted to it (KMSS method, see \code{\link{multKM}}).
}
\usage{
splineKM(x, label = NULL, dl = NULL, n.knots = NULL,
            legend.pos = "bottomright",
            ylab = "ECDF", xlab = "Value",
            col.km = "black", lty.km = 1, lwd.km = 1,
            col.sm = "red", lty.sm = 2, lwd.sm = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numerical data vector (\code{\link{vector}} class).
}
  \item{label}{
Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote left-censored values in \code{x}.
}
  \item{dl}{
Numeric vector of detection limits/thresholds for each element of \code{x} (same length as \code{x}). These must be given on the same scale as \code{x} (use e.g. 0 for detected data). 
}
  \item{n.knots}{Integer or function giving the number of knots used for fitting a cubic smoothing spline to the KM ECDF (see \code{\link{smooth.spline}} for default value).
}
  \item{legend.pos}{
Location of the graph legend. Choose one amongst \code{"bottomleft"}, \code{"bottomright"} (default), \code{"topleft"} or \code{"topright"}.
}
  \item{ylab}{
Title for y-axis.
}
  \item{xlab}{
Title for x-axis.
}
  \item{col.km}{
Plotting color for KM ECDF (see base graphical parameters \code{\link{par}}).
}
  \item{lty.km}{
Line type for KM ECDF (see base graphical parameters \code{\link{par}}).
}
  \item{lwd.km}{
Line width for KM ECDF (see base graphical parameters \code{\link{par}}).
}
  \item{col.sm}{
Plotting color for smoothing spline curve.
}
  \item{lty.sm}{
Line style for smoothing spline curve.
}
  \item{lwd.sm}{
Line width for smoothing spline curve.
}
  \item{\dots}{
Other graphical parameters.
}
}

\value{
Graphical output.
}

\examples{
data(Water)
data(mdl)

# Examine default spline smoothed KM ECDF fit for Potassium and Sulphate
splineKM(Water[,1],label=0,mdl[,1])
splineKM(Water[,4],label=0,mdl[,4],xlim=c(28,41))

# Reduce to 5 knots for Potassium
splineKM(Water[,1],label=0,mdl[,1],n.knots=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
