\name{multLN}
\alias{multLN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative lognormal replacement
}
\description{
This function implements model-based multiplicative lognormal imputation of left-censored values (e.g. values below detection limit, rounded zeros) in compositional data sets.
}
\usage{
multLN(X, label = NULL, dl = NULL, rob = FALSE, random = FALSE,
          z.warning = 0.8, z.delete = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zeros/unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}. If \code{NULL} the column minima are used as thresholds.
}
  \item{rob}{
Logical value. \code{FALSE} provides maximum-likelihood estimates of model parameters (default), \code{TRUE} provides robust estimates (see \code{NADA} package for details).
}
  \item{random}{
Logical value. Values imputed using either estimated geometric mean (\code{FALSE}, default) or random values (\code{TRUE}) below the limit of detection.
}
 \item{z.warning}{Threshold used to identify individual rows or columns including an excess of zeros/unobserved values (to be specify in proportions, default \code{z.warning=0.8}).}
  \item{z.delete}{Logical value. If set to \code{TRUE}, rows/columns identified by \code{z.warning} are omitted in the imputed data set. Otherwise, the function stops in error when rows/columns are identified by \code{z.warning} (default \code{z.delete=TRUE}).}
}
\details{
By default, this function imputes left-censored compositional values by the estimated geometric mean of the values below the corresponding limit of detection or censoring threshold and applies a multiplicative adjustment to preserve the multivariate compositional properties of the samples. Alternatively, imputation can be carried out by random values below the limit of detection (\code{random = TRUE}) based on a normal distribution on the positive real line (see below).

It depends on package \code{NADA} to produce the required model parameter estimates (either maximum likelihood or robust regression on order statistics). It allows for either single (\code{vector} form) or multiple (\code{matrix} form, same size as \code{X}) limits of detection by component. Any threshold value can be set for non-censored elements (e.g. use 0 if no threshold for a particular column or element of the data matrix).

It produces an imputed data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered. Note that a normal distribution on the positive real line is considered. That is, it is defined with respect to a measure according to own geometry of the positive real line, instead of the standard lognormal based on the Lebesgue measure in real space.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set expressed in the original scale.
}
\references{
Mateu-Figueras G, Pawlowsky-Glahn V, Egozcue JJ. The normal distribution in some constrained sample spaces. SORT 2013; 37: 29-56.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43.

Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{lrEM}}, \code{\link{lrSVD}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
}
\examples{

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
              
X_multLN <- multLN(X,label=0,dl=rep(1,6))

# Using ROS for parameter estimation
X_multLNrob <- multLN(X,label=0,dl=rep(1,6),rob=TRUE)

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_multLN2 <- multLN(X,label=0,dl=mdl)

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)

# Using ML for parameter estimation
LPdata_multLN <- multLN(LPdata,label=0,dl=dl) 
# For comparison
LPdata[30:35,1:10]
round(LPdata_multLN[30:35,1:10],1)

# Using random values < dl
LPdata_multRLN <- multLN(LPdata,label=0,dl=dl,random=TRUE)
round(LPdata_multRLN[30:35,1:10],1)

# Two subsets of limits of detection (using e.g. ML parameter estimation)
data(LPdata)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
 # DLs for first 50 samples of LPdata
dl1 <- matrix(rep(1,50),ncol=1)\%*\%dl
 # DLs for last 46 samples of LPdata
dl2 <- matrix(rep(1,46),ncol=1)\%*\%c(1,0.5,0,0,2.5,0,5.5,0.75,0.3,1.5,1,0,0,600,8) 

mdl <- rbind(dl1,dl2)
LPdata_multLN2 <- multLN(LPdata,label=0,dl=mdl)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
