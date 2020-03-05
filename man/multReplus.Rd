\name{multReplus}
\alias{multReplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative simple replacement (plus) 
}
\description{
This function implements an extended version of multiplicative simple imputation (\code{multRepl} function) to simultaneously deal with both zeros (i.e. data below detection limit, rounded zeros) and missing data in compositional data sets.

Note: zeros and missing data must be labelled using 0 and \code{NA} respectively to use this function.
}
\usage{
multReplus(X, dl = NULL, frac = 0.65,
              suppress.print = FALSE, closure = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
}
  \item{frac}{
Fraction of the detection limit/threshold used for imputation (default = 0.65, expressed as a proportion).
}
  \item{suppress.print}{
Suppress printed feedback (\code{suppress.print = FALSE}, default).
}
  \item{closure}{Closure value used to add a residual part for imputation (see below).}
}

\details{
The procedure firstly replaces missing data using the estimated geometric mean based on the observed values and then zeros using \code{frac*dl}. The observed components are applied a multiplicative adjustment to preserve the multivariate compositional properties of the samples.

Note: negative values can be generated when unobserved components are a large portion of the composition, which is more likely for missing data (e.g in major chemical elements) and non-closed compositions. A workaround is to add a residual filling the gap up to the closure/total when possible. This is done internally when a value for \code{closure} is specified (e.g. \code{closure=10^6} if ppm or \code{closure=100} if percentages). The residual is discarded after imputation.

See \code{?multRepl} for more details.

}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set in the same scale as the original.
}
\references{
Martin-Fernandez JA, Barcelo-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. Mathematical Geology 2003; 35: 253-78.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43. DOI: http://dx.doi.org/10.1016/j.aca.2012.12.029.

Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multRepl}}, \code{\link{lrEMplus}}
}
\examples{
# Data set closed to 100 (percentages, common dl = 1%)
# (Note that zeros and missing in the same row are allowed)
X <- matrix(c(26.91,8.08,12.59,31.58,6.45,14.39,
              39.73,41.42,0.00,NA,6.80,12.05,
              NA,35.13,7.96,14.28,35.12,7.51,
              10.85,46.40,31.89,10.86,0.00,0.00,
              10.85,16.27,NA,9.16,19.57,44.15,
              38.09,7.62,23.68,9.70,20.91,0.00,
              NA,9.89,18.04,44.30,9.04,18.73,
              44.41,15.04,7.95,0.00,10.82,21.78,
              11.50,30.33,6.85,13.92,30.82,6.58,
              19.04,42.59,0.00,38.37,0.00,0.00),byrow=TRUE,ncol=6)
              
X_multReplus <- multReplus(X,dl=rep(1,6))

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_multReplus2 <- multReplus(X,dl=mdl)

# Non-closed compositional data set
data(LPdataZM) # (in ppm; 0 is nondetect and NA is missing data)

dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdataZM2 <- subset(LPdataZM,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

\dontrun{
LPdataZM2_multReplus <- multReplus(LPdataZM2,dl=dl2)
# Negative values generated (see e.g. K in sample #64)
}

# Workaround: use residual part to fill up the gap to 10^6 for imputation
LPdataZM2_multReplus <- multReplus(LPdataZM2,dl=dl2,closure=10^6)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
