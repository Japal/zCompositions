\name{multKM}
\alias{multKM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative Kaplan-Meier smoothing spline (KMSS) replacement
}
\description{
This function implements non-parametric multiplicative KMSS imputation of left-censored values (e.g. values below detection limit, rounded zeros) in compositional data sets. It is based on simulation from a smoothing spline fitted to the Kaplan-Meier (KM) estimate of the empirical cumulative distribution function (ECDF) of the data.
}
\usage{
multKM(X, label = NULL, dl = NULL, n.draws = 1000, n.knots = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
}
  \item{n.draws}{
Number of random draws from the inverse KM ECDF generated to produce an averaged imputed value (\code{n.draws=1000}, default).
}
  \item{n.knots}{Integer or function giving the number of knots used for fitting a cubic smoothing spline to the KM ECDF (see \code{\link{smooth.spline}} for default value). It allows for a vector or list of settings per column of \code{X}.
}
}
\details{
This function imputes left-censored compositional values by averaging (geometric mean) n random draws (\code{n.draws} argument) from a cubic smoothing spline curve fitting the inverse KM ECDF below the corresponding limit of detection or censoring threshold. It then applies a multiplicative adjustment to preserve the multivariate compositional properties of the samples. It allows for either single (\code{vector} form) or multiple (\code{matrix} form, same size as \code{X}) limits of detection by component. Although note that it is equivalent to simple substitution by the limit of detection for singly censored components. Any threshold value can be set for non-censored elements (e.g. use 0 if no threshold for a particular column or element of the data matrix).

It produces an imputed data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.

The level of smoothing of the estimated spline can be controlled by the \code{n.knots} argument. The function \code{\link{splineKM}} can assist in choosing a finer value, although the default setting works generally well.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set in the same scale as the original.
}

\references{
Palarea-Albaladejo J. and Martin-Fernandez J.A. zCompositions -- R package
for multivariate imputation of left-censored data under a compositional approach.
Chemometrics and Intelligent Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{splineKM}}, \code{\link{lrEM}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{cmultRepl}}
}
\examples{
data(Water)
data(mdl) # matrix of limits of detection for Water

Water_multKM <- multKM(Water,label=0,dl=mdl)

# Different smoothing degree by component
Water_multKM2 <- multKM(Water,label=0,dl=mdl,n.knots=c(25,50,30,75))

# Easy to use for KM multiple imputation (m = 10)
Water.mi <- vector("list",length=10)
for (m in 1:10){
  Water.mi[[m]] <- multKM(Water,label=0,dl=mdl,n.draws=1)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
