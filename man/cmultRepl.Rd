\name{cmultRepl}
\alias{cmultRepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian-Multiplicative replacement of count zeros
}
\description{
This function implements methods for imputing zeros in compositional count data sets based on a Bayesian-multiplicative replacement.
}
\usage{
cmultRepl(X, label = 0,
          method = c("GBM","SQ","BL","CZM","user"), output = c("prop","counts"),
          delta = 0.65, threshold = 0.5, correct = TRUE, t = NULL, s = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Count data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote count zeros in \code{X} (default \code{label=0}).}
  \item{method}{
Geometric Bayesian multiplicative (\code{GBM}, default); square root BM (\code{SQ}); Bayes-Laplace BM (\code{BL}); count zero multiplicative (\code{CZM}); user-specified hyper-parameters (\code{user}).
}
  \item{output}{Output format: imputed proportions (\code{prop}, default) or \emph{pseudo}-counts (\code{counts}). 
}
  \item{delta}{
If \code{method="CZM"}, fraction of the upper threshold used to impute zeros (default \code{delta=0.65}). Also, fraction of the lowest estimated probability used to correct imputed proportions falling above it (when \code{correct=TRUE}).
}
  \item{threshold}{
For a vector of counts, factor applied to the quotient 1 over the number of trials (sum of the counts) used to produce an upper limit for replacing zero counts by the \code{CZM} method (default \code{threshold=0.5}).
}
  \item{correct}{
Logical vector setting whether imputed proportions falling above the lowest estimated probability for a multinomial part must be corrected or not (default \code{correct=TRUE}).
}
  \item{t}{
If \code{method="user"}, user-specified \emph{t} hyper-parameter of the Dirichlet prior distribution for each count vector (row) in \code{X}. It must be a matrix of the same dimensions as \code{X}.
}
  \item{s}{If \code{method="user"}, user-specified \emph{s} hyper-parameter of the Dirichlet prior distribution for each count vector (row) in \code{X}. It must be a vector of length equal to the number of rows of \code{X}.
}
}
\details{
Zero counts, assumed to be a consequence of the sampling process, are imputed under a Bayesian paradigm (\code{GBM}, \code{SQ} or \code{BL} method) by posterior estimates of the multinomial probabilities generating the counts, assuming a Dirichlet prior distribution. The argument \code{method} sets the Dirichlet hyper-parameters \code{t} (priori estimates of multinomial probabilities) and \code{s} (\emph{strength}). The user can specify their own by setting \code{method="user"} and entering them as \code{t} and \code{s} arguments. Note that, under certain circumstances (see references for details), these methods can generate imputed proportions falling above the lowest estimated probability  of a multinomial part (\emph{c}/\emph{n}, where \emph{c} is the count and \emph{n} is the number of trials). In such cases, the imputation is corrected by using a fraction (\code{delta}) of the minimum \emph{c}/\emph{n} for that part. Lastly, the non-zero parts are multiplicatively adjusted according to their compositional nature.

On the other hand, \code{method="CZM"} uses multiplicative simple replacement (\code{\link{multRepl}}) on the matrix of estimated probabilities. The upper limit and the fraction delta used are specified by, respectively, the arguments \code{threshold} and \code{delta}. Suggested values are \code{threshold=0.5} (so the upper limit for a multinomial probability turns out to be 0.5/\emph{n}), and \code{delta=0.65} (so the imputed proportion is 65\% of the upper limit).
}
\value{
By default (\code{output="prop"}) the function returns an imputed data set (\code{\link{data.frame}} class) in proportions (estimated probabilities). Alternatively (\code{output="counts"}), these proportions can be re-scaled to produce a  compositionally-equivalent matrix of \emph{pseudo}-counts (\code{output="counts"}).

When \code{correct=TRUE}, the number of times, if any, an imputed proportion was corrected to fall below the minimum estimated multinomial probability is printed.
}
\references{
Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. Bayesian-multiplicative treatment of count zeros in compositional data sets. Statistical Modelling 2015; 15 (2): 134-158.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}}

\examples{
data(Pigs)

# GBM method and matrix of estimated probabilities
Pigs.GBM <- cmultRepl(Pigs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
