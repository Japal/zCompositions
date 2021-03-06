\name{zPatterns}
\alias{zPatterns}
\title{
Find and display patterns of zeros/missing values in a data set
}
\description{
This function summarises the patterns of zero and/or missing values in a data set and returns a vector of pattern numbers.
}
\usage{
zPatterns(X, label = NULL, plot = TRUE,
             axis.labels = c("Component", "Pattern ID"),
             bar.ordered = as.character(c(FALSE,FALSE)),
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
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to identify zeros/unobserved values in \code{X}.}
  \item{plot}{Logical value indicating whether a graphical summary of the patterns is produced or not (default \code{plot=TRUE}).}
  \item{axis.labels}{Vector of axis labels for the table of patterns (format \code{c("x-axis","y-axis")}).}
  \item{bar.ordered}{Vector of logical values to order table of patterns according to frequencies by patterns, component or both; with the first element referring to the patterns and the second to the components (default \code{c(FALSE,FALSE)}).}
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
Vector (factor type) of pattern IDs corresponding to each row of \code{X}.

By default, a summary table is printed showing patterns in the data according to \code{label} and some summary statistics: number of zero/missing components by pattern (No.Unobs), pattern frequency in percentage, percentage zero/missing values by component (column) and overall percentage of zero/missing values in the data set. The symbols + and - indicate, respectively, zero/missing and observed components within each pattern. A graphical version of the summary table is returned including barplots on the margins displaying percentage zero/missing and compositional geometric means by pattern (if \code{show.means=TRUE}; expressed in percentage scale). Common arithmetic means can be also shown for the case of ordinary data (\code{type.means="am"}), however this is not recommended for compositional data.

The patterns are assigned ID number and by default arranged in the table in the same order as they are found in the data set. The argument \code{bar.ordered} can be used to re-arrange the display according to frequencies of patterns, of unobserved values by component or both.

A warning message is shown if zeros or NA values not identified by \code{label} are present in the data set. These will be ignored for the graphical display and numerical summaries of patterns, which will be only based on \code{label}.

Check out 'plus' functions to deal with zeros and missing data simultaneously.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

\seealso{
\code{\link{lrEM}}, \code{\link{lrEMplus}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multReplus}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
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

# Show compositional geometric means (in \%) per zero pattern
zPatterns(LPdata,label=0,show.means=TRUE)

# Same but orderer by pattern frequency and incidence of zeros by component
zPatterns(LPdata,label=0,bar.ordered=c(TRUE,TRUE),,bar.labels=TRUE,show.means=TRUE)

# Data set with zeros and missing data (0 = zero; NA = missing) (see lrEMplus function).
data(LPdataZM)

# Show missingness patterns only
zPatterns(LPdataZM,label=NA)

# Show zero patterns only and means by pattern based on available data
# (blanks indicate not enough data available for computation)
zPatterns(LPdataZM,label=0,show.means=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
