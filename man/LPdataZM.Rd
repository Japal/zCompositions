\name{LPdataZM}
\alias{LPdataZM}
\docType{data}
\title{
La Paloma data set (incl. zeros and missing data)
}
\description{
96 samples of a 15-part geochemical composition in micrograms/gram from La Paloma stream (Venezuela). For more details see Montero-Serrano et al. (2010).

Duplicate of the \code{LPdata} data set including 2.43\% missing at random cells (36.45\% samples with missing data; coded as \code{NA}) along with 6.11\% values below the limit of detection (coded as 0).

}
\usage{data(LPdataZM)}
\format{
  A \code{\link{data.frame}} with 96 observations on the following 15 variables.
  \describe{
    \item{\code{Cr}}{a numeric vector}
    \item{\code{B}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{V}}{a numeric vector}
    \item{\code{Cu}}{a numeric vector}
    \item{\code{Ti}}{a numeric vector}
    \item{\code{Ni}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{Sr}}{a numeric vector}
    \item{\code{La}}{a numeric vector}
    \item{\code{Ce}}{a numeric vector}
    \item{\code{Ba}}{a numeric vector}
    \item{\code{Li}}{a numeric vector}
    \item{\code{K}}{a numeric vector}
    \item{\code{Rb}}{a numeric vector}
  }
}
\references{
Montero-Serrano JC, Palarea-Albaladejo J, Martin-Fernandez JA, and Martinez-Santana M and Gutierrez-Martin JV. Multivariate analysis applied to chemostratigraphic data: identification of chemofacies and stratigraphic correlation, Sedimentary Geology 2010; 228(3-4): 218-228 .
}

\seealso{
\code{\link{LPdata}}
}

\examples{
data(LPdataZM)

zPatterns(LPdataZM,label=0) # Includes zero + missing data patterns
}
\keyword{datasets}
