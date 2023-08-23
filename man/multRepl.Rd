\name{multRepl}
\alias{multRepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiplicative simple replacement
}
\description{
This function implements non-parametric multiplicative simple imputation of left-censored (e.g. values below detection limit, rounded zeros) and missing data in compositional data sets.
}
\usage{
multRepl(X, label = NULL, dl = NULL, frac = 0.65, imp.missing = FALSE, closure = NULL,
            z.warning=0.8, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional vector (\code{numeric} class) or data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zeros/unobserved values in \code{X}.}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
}
  \item{frac}{
Fraction of the detection limit/threshold used for imputation (default = 0.65, expressed as a proportion).
}
 \item{imp.missing}{If \code{TRUE} then unobserved values identified by \code{label} are treated as missing data (default = \code{FALSE}).}
  \item{closure}{Closure value used to add a residual part for imputation (see below).}
 \item{z.warning}{Threshold used to delete individual rows or columns including an excess of zeros/unobserved values (to be specify in proportions, default \code{z.warning=0.8}).}
  \item{delta}{This argument has been deprecated and replaced by \code{frac} (see package's NEWS for details).}
}
\details{
This function imputes left-censored compositional values by a given fraction \code{frac} of the corresponding limit of detection and applies a multiplicative adjustment to preserve the multivariate compositional properties of the samples. It allows for either single (\code{vector} form) or multiple (\code{matrix} form, same size as \code{X}) limits of detection by component. Any threshold value can be set for non-censored elements (e.g. use 0 if no threshold for a particular column or element of the data matrix).

\emph{Missing data imputation}: missing data can be imputed by setting \code{imp.missing = TRUE}. They are replaced by the estimated column geometric mean from observed values. The non-missing parts in the composition are applied multiplicative adjustment. The argument \code{dl} and \code{frac} are ignored and \code{X} is require to be a data matrix in this case.

Note: negative values can be generated when unobserved components are a large portion of the composition, which is more likely for missing data (e.g in major chemical elements) and non-closed compositions. A workaround is to add a residual filling the gap up to the closure/total when possible. This is done internally when a value for \code{closure} is specified (e.g. \code{closure=10^6} if ppm or \code{closure=100} if percentages). The residual is discarded after imputation.

This function produces an imputed data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered. Note that this adjustment only applies to data sets and not when a single composition is entered. In this latter case, the composition is treated as a closed vector.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional vector or data set expressed in the original scale.
}
\references{
Martin-Fernandez JA, Barcelo-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. Mathematical Geology 2003; 35: 253-78.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43.

Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{lrEM}}, \code{\link{lrSVD}}, \code{\link{lrDA}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
}
\examples{
# A compositional vector (NA indicates nondetect)
y <- c(0.6,NA,0.25,0.03,0.12,NA)
dl <- c(0,0.01,0,0,0,0.005)
# Using the default frac = 0.65
yr <- multRepl(y,label=NA,dl=dl)
round(yr,4)

# Data set closed to 100 (percentages, common dl = 1%)
X <- matrix(c(26.91,8.08,12.59,31.58,6.45,14.39,
              39.73,26.20,0.00,15.22,6.80,12.05,
              10.76,31.36,7.10,12.74,31.34,6.70,
              10.85,46.40,31.89,10.86,0.00,0.00,
              7.57,11.35,30.24,6.39,13.65,30.80,
              38.09,7.62,23.68,9.70,20.91,0.00,
              27.67,7.15,13.05,32.04,6.54,13.55,
              44.41,15.04,7.95,0.00,10.82,21.78,
              11.50,30.33,6.85,13.92,30.82,6.58,
              19.04,42.59,0.00,38.37,0.00,0.00),byrow=TRUE,ncol=6)
              
X_multRepl <- multRepl(X,label=0,dl=rep(1,6))

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_multRepl2 <- multRepl(X,label=0,dl=mdl)

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdata_multRepl <- multRepl(LPdata,label=0,dl=dl)

# Two subsets of limits of detection
data(LPdata)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
 # DLs for first 50 samples of LPdata
dl1 <- matrix(rep(1,50),ncol=1)\%*\%dl
 # DLs for last 46 samples of LPdata
dl2 <- matrix(rep(1,46),ncol=1)\%*\%c(1,0.5,0,0,2.5,0,5.5,0.75,0.3,1.5,1,0,0,600,8) 

mdl <- rbind(dl1,dl2)
LPdata_multRepl2 <- multRepl(LPdata,label=0,dl=mdl)

# Data set with missing values closed to 100 (percentages)
X <- matrix(c(10.47,8.58,59.72,19.30,1.93,
              12.13,7.44,62.87,16.37,1.19,
                 NA,7.30,75.91,16.79,NA,
               9.77,7.80,65.68,14.78,1.97,
              10.79,9.55,65.87,12.41,1.38,
              14.54,8.18,64.55,12.73,NA,
              12.28,7.58,66.01,12.93,1.20,
              28.09,22.92,NA,40.11,8.88,
               7.02,6.30,75.65,11.03,NA),byrow=TRUE,ncol=5)
              
X_multReplMiss <- multRepl(X,label=NA,imp.missing=TRUE)

# Non-closed compositional data set
data(LPdata) # (in ppm units)
# Treating zeros as missing data for illustration purposes only
LPdata_multReplMiss <- multRepl(LPdata,label=0,imp.missing=TRUE)
# Negative values generated (see e.g. K and Rb in sample #60)

# Workaround: use residual part to fill up the gap to 10^6 for imputation
LPdata_multReplMiss2 <- multRepl(LPdata,label=0,imp.missing=TRUE,closure=10^6)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
