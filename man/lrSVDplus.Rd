\name{lrSVDplus}
\alias{lrSVDplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-ratio SVD algorithm (plus) 
}
\description{
This function implements an extended version of the log-ratio SVD algorithm (\code{lrSVD} function) to simultaneously deal with both zeros (i.e. data below detection limit, rounded zeros) and missing data in compositional data sets.

Note: zeros and missing data must be labelled using 0 and \code{NA} respectively to use this function.
}
\usage{
lrSVDplus(X, dl = NULL, frac = 0.65,
             ncp = 2, beta = 0.5, method = c("ridge", "EM"), row.w = NULL,
             coeff.ridge = 1, threshold = 1e-04, seed = NULL, nb.init = 1,
             max.iter = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).}

\item{dl}{Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.}

\item{frac}{Parameter for initial multiplicative simple replacement of left-censored data (see \code{\link{multRepl}}) (default = 0.65).}

\item{ncp}{Number of components in low-rank matrix approximation (default = 2).}

\item{beta}{Weighting parameter, balance between the two conditions in objective function (default = 0.5).}

\item{method}{Parameter estimation method for the iterative algorithm (\code{method = "ridge"}, default).}

\item{row.w}{row weights (default = NULL, a vector of 1 for uniform row weights).}

\item{coeff.ridge}{Used when \code{method = "ridge"} (default = 1).}

\item{threshold}{Threshold for assessing convergence (default = 1e-04).}

\item{seed}{Seed for random initialisation of the algorithm (default \code{seed = NULL}, unobserved values initially imputed by the column mean).}

\item{nb.init}{Number of random initialisations (default = 1).}

\item{max.iter}{Maximum number of iterations for the algorithm (default = 1000).}

\item{...}{Further arguments.}
}

\details{
The procedure starts with an initial imputation of zeros (using simple replacement with \code{frac*dl}) and missing values (using geometric mean imputation from observed data). Subsequently, the iterative algorithm is run until convergence (see \code{?lrSVD} for more details).
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set in the same scale as the original.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lrSVD}}
}
\examples{
# Data set closed to 100 (percentages, common dl = 1%)
# (Note that zeros and missing in the same row or column are allowed)
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
              
X_lrSVDplus <- lrSVDplus(X,dl=rep(1,6))

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_lrSVDplus2 <- lrSVDplus(X,dl=mdl)

# Non-closed compositional data set
data(LPdataZM) # (in ppm; 0 is nondetect and NA is missing data)

dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdataZM2 <- subset(LPdataZM,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

LPdataZM2_lrSVDplus <- lrSVDplus(LPdataZM2,dl=dl2)
}
