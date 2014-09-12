\name{Water}
\alias{Water}
\docType{data}
\title{
Water data set
}
\description{
100 simulated samples of a 4-part groundwater composition in percentage subject to multiple limits of detection by component. The associated matrix of limits of detection is stored in \code{mdl}.
}
\usage{data(Water)}
\format{
  A \code{\link{data.frame}} with 100 observations on the following 4 variables.
  \describe{
    \item{\code{Potassium}}{a numeric vector}
    \item{\code{Arsenic}}{a numeric vector}
    \item{\code{Calcium}}{a numeric vector}
    \item{\code{Sulphate}}{a numeric vector}

  }
}

\details{Three limits of detection (0.75, 1 and 1.25) were considered for Potassium, four for Arsenic (1.5, 3, 4 and 5), two for Sulphate (29 and 35) and no one for Calcium. In the case of Sulphate, the detection limit equal to 29 is the minimum value registered for that component. All nondetects coded as 0.}

\examples{
data(Water)
zPatterns(Water,label=0)
}
\keyword{datasets}
