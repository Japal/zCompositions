\name{zPatterns}
\alias{zPatterns}
\title{
Find and display patterns of zeros or unobserved values in a data set
}
\description{
This function summarises the patterns of zero or unobserved values in a data set and returns a vector of pattern numbers.
}
\usage{
zPatterns(X, label = NULL, plot = TRUE,
             axis.labels = c("Component", "Pattern number"),
             bar.colors = c("red3", "red3"), bar.labels = FALSE,
             show.means = FALSE, round.means = 2, cex.means = 1,
             type.means = c("cgm","am"),
             cell.colors = c("dodgerblue", "white"),
             cell.labels = c(label, paste("No", label)), cex.axis = 1.1,
             grid.color = "black", grid.lty = "dotted",
             legend = TRUE, suppress.print = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data set (\code{\link{matrix}} or \code{\link{data.frame}} class).}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zero or unobserved values in \code{X}.}
  \item{plot}{Logical value indicating whether a graphical summary of the patterns is produced or not (default \code{plot=TRUE}).}
  \item{axis.labels}{Vector of axis labels for the table of patterns (format \code{c("x-axis","y-axis")}).}
  \item{bar.colors}{Colors for the margin barplots (format \code{c("col.top","col.right")}).}
  \item{bar.labels}{Logical value indicating if labels showing percentages must be added to the margin barplots (default \code{bar.labels=FALSE}).}
  \item{show.means}{Logical value indicating if mean values by pattern are shown on the graphical summary table (default \code{show.means=FALSE}).}
  \item{round.means}{When \code{show.means=TRUE}, number of decimal places for the mean values shown (2=default).}
  \item{cex.means}{When \code{show.means=TRUE}, numeric character expansion factor; character size for the mean values shown (1=default).}
  \item{type.means}{When \code{show.means=TRUE}, statistic used for computing the means. Either compositional geometric mean (\code{type.means=cgm}, in percentage units, default) or standard arithmetic mean (\code{type.means=am}).}
  \item{cell.colors}{Vector of colors for the table cells (format \code{c("col.unobserved","col.observed")}).}
  \item{cell.labels}{Labels for the cells (format \code{c("Unobserved","Observed")}, default \code{c(label,paste("No",label))}).}
  \item{cex.axis}{Axis labels scaling factor relative to default.}
  \item{grid.color}{Color of the grid lines (default \code{"black"}).}
  \item{grid.lty}{Style of the grid lines (default \code{"dotted"}, see \code{lty} in \code{\link{par}}).}
  \item{legend}{Logical value indicating if a legend must be included (default \code{legend=TRUE}).}
  \item{suppress.print}{Suppress printed feedback (default \code{suppress.print=FALSE}).}
  \item{...}{Other graphical parameters.}
}
\value{
Vector of pattern IDs corresponding to each row of \code{X}.

It also prints a summary table (if \code{suppress.print=FALSE}) showing zero or non-observation patterns and some summary statistics: number of zero/unobserved components by pattern (No.Unobs), pattern frequency (absolute and percentage), percentage zero/unobserved values by component (column) and overall percentage of zero/unobserved values in the data set. The symbols + and - indicate, respectively, zero/unobserved and observed components within each pattern. A graphical version of the summary table is returned including barplots on the margins displaying percentage zero/unobserved and percentage compositional geometric means by pattern (if \code{show.means=TRUE}). Note that the posibility of computing ordinary arithmetic means (\code{type.means="am"}) is also included to give the function a broader use, however this measure should not be used if the data are compositional.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lrEM}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
}
\examples{
data(LPdata)

pattern.ID <- zPatterns(LPdata,label=0)

LPdata[pattern.ID==5,]
LPdata[pattern.ID==7,]
LPdata[pattern.ID==10,]

# Modify cell labels and show percentages along with barplots
pattern.ID <- zPatterns(LPdata,label=0,
              cell.labels=c("Zero","Non-zero"),bar.labels=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
