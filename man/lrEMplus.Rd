\name{lrEMplus}
\alias{lrEMplus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-ratio EM algorithm (plus) 
}
\description{
This function implements an extended version of the log-ratio EM algorithm (\code{lrEM} function) to simultaneously deal with both zeros (i.e. data below detection limit, rounded zeros) and missing data in compositional data sets.

Note: zeros and missing data must be labelled using 0 and \code{NA} respectively to use this function.
}
\usage{
lrEMplus(X, dl = NULL, rob = FALSE,
            ini.cov = c("complete.obs", "multRepl"), frac = 0.65,
            tolerance = 0.0001, max.iter = 50, rlm.maxit = 150,
            suppress.print = FALSE, closure = NULL,
            z.warning=0.8, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
}
  \item{rob}{
Logical value. \code{FALSE} provides maximum-likelihood estimates of model parameters (default), \code{TRUE} provides robust parameter estimates.
}
  \item{ini.cov}{Initial estimation of either the log-ratio covariance matrix (ML estimation) or unobserved data (robust estimation). It can be based on either complete observations (\code{"complete.obs"}, default) or multiplicative simple replacement (\code{"multRepl"}).}
  \item{frac}{
If \code{ini.cov="multRepl"}, parameter for initial multiplicative simple replacement of left-censored data (see \code{\link{multRepl}}) (default = 0.65).
}
  \item{tolerance}{Convergence criterion (default = 0.0001).}
  \item{max.iter}{
Maximum number of iterations (default = 50).
}
  \item{rlm.maxit}{
If \code{rob=TRUE}, maximum number of iterations for the embedded robust regression estimation (default = 150; see \code{\link{rlm}} in \code{MASS} package for details).
}
  \item{suppress.print}{
Suppress printed feedback (\code{suppress.print = FALSE}, default).
}
  \item{closure}{Closure value used to add a residual part if needed when \code{ini.cov="multRepl"} is used (see \code{?multRepl}).}
  \item{z.warning}{
Print a warning message if the proportion of zeros/unobserved values in any column or row exceeds a given threshold (default \code{z.warning=0.8}).
}
  \item{delta}{This argument has been deprecated and replaced by \code{frac} (see package's NEWS for details).}
}


\details{
The procedure starts with an initial imputation of either zeros (using simple replacement with \code{frac*dl}) or missing values (using geometric mean imputation from observed data) depending of which problem is the least frequent in the data set. Subsequently, iterative calls to \code{lrEM} replace zeros and missing data alternately until convergence to a stable solution or the maximum number of iterations is reached.

See \code{?lrEM} for more details.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set in the same scale as the original. The number of iterations required for convergence is also printed (this can be suppressed by setting \code{suppress.print=TRUE}).
}
\references{
Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. Model-based replacement of rounded zeros in compositional data: classical and robust approaches. Computational Statistics & Data Analysis 2012; 56: 2688-2704.

Palarea-Albaladejo J, Martin-Fernandez JA, Gomez-Garcia J. A parametric approach for dealing with compositional rounded zeros. Mathematical Geology 2007; 39: 625-45.

Palarea-Albaladejo J, Martin-Fernandez JA. A modified EM alr-algorithm for replacing rounded zeros in compositional data sets. Computers & Geosciences 2008; 34(8): 902-917.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43. DOI: \url{https://doi.org/10.1016/j.aca.2012.12.029}.

Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lrEM}}
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
              
X_lrEMplus <- lrEMplus(X,dl=rep(1,6),ini.cov="multRepl")
X_roblrEMplus <- lrEMplus(X,dl=rep(1,6),ini.cov="multRepl",rob=TRUE,max.iter=4)

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_lrEMplus2 <- lrEMplus(X,dl=mdl,ini.cov="multRepl")

# Non-closed compositional data set
data(LPdataZM) # (in ppm; 0 is nondetect and NA is missing data)

dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdataZM2 <- subset(LPdataZM,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

LPdataZM2_lrEMplus <- lrEMplus(LPdataZM2,dl=dl2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
