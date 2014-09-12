\name{lrDA}
\alias{lrDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-ratio DA algorithm
}
\description{
This function implements a simulation-based Data Augmentation (DA) algorithm to impute left-censored values (e.g. values below detection limit, rounded zeros) via coordinates representation of compositional data sets which incorporate the information of the relative covariance structure. Multiple imputation estimates can be also obtained from the output.
}
\usage{
lrDA(X, label = NULL, dl = NULL,
        ini.cov=c("lrEM", "complete.obs", "multRepl"),
        delta = 0.65, n.iters = 1000, m = 1, store.mi = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
}
  \item{ini.cov}{Initial estimation of the log-ratio covariance matrix. It can be based on lrEM estimation (\code{"lrEM"}, default), complete observations (\code{"complete.obs"}) or multiplicative simple replacement (\code{"multRepl"}).}
  \item{delta}{
If \code{ini.cov="multRepl"}, delta parameter for initial multiplicative simple replacement (\code{\link{multRepl}}) in proportions (default = 0.65).
}
  \item{n.iters}{
Number of iterations for the DA algorithm (default = 1000).
}
  \item{m}{
Number of multiple imputations (default = 1).
}
  \item{store.mi}{
Logical value. If \code{m>1} creates a list with \code{m} imputed data matrices. (\code{store.mi=FALSE}, default).
}
}


\details{
After convergence of the Markov chain Monte Carlo (MCMC) iterative process to its steady state, this function imputes left-censored compositional parts by simulated values from their posterior predictive distributions through coordinates representation, given the information from the observed data and the censoring thresholds. It allows for either single (\code{vector} form) or multiple (\code{matrix} form, same size as \code{X}) limits of detection by component. Any threshold value can be set for non-censored elements (e.g. use 0 if no threshold for a particular column or element of the data matrix).

It produces imputed data sets on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.

The common conjugate normal inverted-Wishart distribution with non-informative prior has been assumed for the model parameters in the coordinates space. Under this setting, convergence is expected to be fast (\code{n.iters} set to 1000 by default). Besides, considering EM parameter estimates as initial point for the DA algorithm (\code{ini.cov="lrEM"}) assures faster convergence by starting near the centre of the posterior distribution. Note that the procedure is based on the oblique additive log-ratio (alr) transformation to simplify calculations and alleviates computational burden.

By setting \code{m} greater than 1, the procedure also allows for multiple imputations of the censored values drawn at regular intervals after convergence. In this case, in addition to the burn-in period for convergence, \code{n.iters} determines the gap, large enough to prevent from correlated values, between successive imputations. The total number of iterations is then \code{n.iters*m}. By default, a single imputed data set results from averaging the \code{m} imputations in the space of coordinates. If \code{store.mi=TRUE}, a list with \code{m} imputed data sets is generated instead.

In the case of censoring patterns involving samples containing only one observed component, these are imputed by multiplicative simple replacement (\code{\link{multRepl}}) and a warning message identifying them is printed.

%Note that when \code{ini.cov="lrEM"} is chosen, the embedded lrEM procedure (\code{\link{lrEM}}) will consider multiplicative replacement (\code{\link{multRepl}}), with delta parameter equal to \code{delta}, for initial covariance estimate. This is simply done in order to prevent it from potentially failing due to matrix singularity issues if working with few complete samples in the data set.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set or a \code{\link{list}} of imputed data sets if multiple imputation is carried out (\code{m>1}) and \code{store.mi=TRUE}.
}
\references{
Palarea-Albaladejo J, Martin-Fernandez JA, Olea, RA. A bootstrap estimation scheme for chemical compositional data with nondetects. Journal of Chemometrics 2014; 28: 585-599.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{lrEM}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
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

# Imputation by single simulated values
X_lrDA <- lrDA(X,label=0,dl=rep(1,6),ini.cov="multRepl",n.iters=150)

# Imputation by multiple imputation (m = 5, one imputation every 150 iterations)
X_milrDA <- lrDA(X,label=0,dl=rep(1,6),ini.cov="multRepl",m=5,n.iters=150)

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_lrDA2 <- lrDA(X,label=0,dl=mdl,ini.cov="multRepl",n.iters=150)

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdata2 <- subset(LPdata,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

\dontrun{ # May take a little while
LPdata_lrDA <- lrDA(LPdata2,label=0,dl=dl2)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
