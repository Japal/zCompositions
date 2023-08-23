\name{lrSVD}
\alias{lrSVD}
\title{Log-ratio SVD algorithm}
\description{
This function implements an iterative algorithm to impute left-censored data (e.g. values below detection limit,
rounded zeros) based on the singular value decomposition (SVD) of a compositional data set.

This function can be also used to impute missing data instead by setting \code{imp.missing = TRUE} (see \code{\link{lrSVDplus}} to treat censored and missing data simultaneously).
}
\usage{
lrSVD(X, label = NULL, dl = NULL, frac = 0.65, ncp = 2, 
         imp.missing=FALSE, beta = 0.5, method = c("ridge", "EM"),
         row.w = NULL, coeff.ridge = 1, threshold = 1e-04, seed = NULL,
         nb.init = 1, max.iter = 1000, z.warning = 0.8, ...)
}
\arguments{
\item{X}{Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).}

\item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zeros/unobserved values in \code{X}.}

\item{dl}{Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.}

\item{frac}{Parameter for initial multiplicative simple replacement of left-censored data (see \code{\link{multRepl}}) (default = 0.65).}

\item{ncp}{Number of components for low-rank matrix approximation (default = 2).}

\item{imp.missing}{If \code{TRUE} then unobserved data identified by \code{label} are treated as missing data (default = \code{FALSE}).}

\item{beta}{Weighting parameter, balance between the two conditions in objective function (default = 0.5).}

\item{method}{Parameter estimation method for the iterative algorithm (\code{method = "ridge"}, default).}

\item{row.w}{row weights (default = NULL, a vector of 1 for uniform row weights).}

\item{coeff.ridge}{Used when \code{method = "ridge"} (default = 1).}

\item{threshold}{Threshold for assessing convergence (default = 1e-04).}

\item{seed}{Seed for random initialisation of the algorithm (default \code{seed = NULL}, unobserved values initially imputed by the column mean).}

\item{nb.init}{Number of random initialisations (default = 1).}

\item{max.iter}{Maximum number of iterations for the algorithm (default = 1000).}

\item{z.warning}{Threshold used to delete individual rows or columns including an excess of zeros/unobserved values (to be specify in proportions, default \code{z.warning=0.8}).}

\item{...}{Further arguments.}

}
\details{
This function implements an efficient imputation algorithm particularly suitable for the case of continuous high-dimensional and wide compositional data sets (more columns than rows), although it is equally applicable to regular data sets. It is based on a low-rank representation of the data set by a principal components (PC) model as derived by singular value decomposition (SVD) of the data matrix, extending recent work on principal component imputation and matrix completion methods to the case of censored compositional data (the code builds on the function \code{imputePCA}; see \code{missMDA} package for more details). A preliminary imputation by multiplicative replacement (see \code{\link{multRepl}}) is conducted to initiate the iterative algorithm in log-ratio coordinates. Two steps, estimation of latent PC model loadings and imputation of empty data matrix cells using the model, are iteratively repeated until convergence. Parameter fitting in this context is performed by a regularisation method (ridge regression in this case) or by the expectation-maximisation (EM) algorithm. Regularization has been shown generally preferable and it is set as default method (note the regularisation parameter \code{coeff.ridge} set to 1 by default. If it is < 1 the result is closer to EM estimation, whereas for values > 1 it is closer to mean estimation).

An imputed data set is produced on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.

\emph{Missing data imputation}

When \code{imp.missing = TRUE}, unobserved values are treated as general missing data. For this case, the argument \code{label} indicates the unique label for missing values and the argument \code{dl} is ignored.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set expressed in the original scale.
}
\references{
Palarea-Albaladejo, J, Antoni Martín-Fernández, J, Ruiz-Gazen, A, Thomas-Agnan, C. lrSVD: An efficient imputation algorithm for incomplete high-throughput compositional data. Journal of Chemometrics 2022; 36: e3459.
}
\seealso{
\code{\link{zPatterns}}, \code{\link{lrSVD}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}, \code{\link{lrSVDplus}}
}
\examples{
 # Data set closed to 100 (percentages, common dl = 1\%)
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
 
 X_lrSVD<- lrSVD(X,label=0,dl=rep(1,6))
 
 # Multiple limits of detection by component
 mdl <- matrix(0,ncol=6,nrow=10)
 mdl[2,] <- rep(1,6)
 mdl[4,] <- rep(0.75,6)
 mdl[6,] <- rep(0.5,6)
 mdl[8,] <- rep(0.5,6)
 mdl[10,] <- c(0,0,1,0,0.8,0.7)
 
 X_lrSVD2 <- lrSVD(X,label=0,dl=mdl)
 
 # Non-closed compositional data set
 data(LPdata) # data (ppm/micrograms per gram)
 dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
 LPdata2 <- subset(LPdata,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
 dl2 <- dl[-c(5,7,10)]
 
 LPdata2_lrSVD <- lrSVD(LPdata2,label=0,dl=dl2)
 
 # Treating zeros as general missing data for illustration purposes only
 LPdata2_miss <- lrSVD(LPdata2,label=0,imp.missing=TRUE)
}
