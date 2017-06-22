#' @title Variation array for grouped data
#'
#' @description This function returns overall and separate variation arrays for groups
#' in a compositional data set. Groups can be defined by either zero/unobserved data patterns or
#' by a grouping factor in fully observed zero-free data sets.
#'
#' @details This function is mainly aimed to investigate heterogeneous relative variation
#' structures in compositional data sets containing zeros or unobserved values. For each pattern of zero or unobserved values,
#' log-ratio variances (upper triangle of variation matrix) and means (lower triangle of variation matrix) are computed from the
#' available data. Note that (1) NAs are produced for log-ratio variances and means in groups containing less than two observations,
#' and (2) at least two components must be available in each group to compute log-ratios.
#'
#' The overall estimate is obtained across groups by pairwise deletion. Note that, unlike the ordinary \code{\link{var}}
#' function, maximum likelihood estimates of the variances are computed. That is,
#' the observed sum of squares is divided by the corresponding number of observations n and not by n-1.
#'
#' Group-wise variation arrays can be obtained from fully observed zero-free data by setting a grouping factor
#' using the argument \code{groups}.
#'
#' @seealso \code{\link{zPatterns}}.
#'
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param label Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zero or unobserved data in \code{X} (\code{label = 0}, default).
#' @param groups Grouping factor in fully observed zero-free data sets (\code{groups = NULL}, default).
#' @param suppress.print Suppress printed feedback (\code{suppress.print = FALSE}, default).
#' @return List of variation arrays by pattern/group and overall.
#'
#' @examples
#' data(Water)
#' zPatterns(Water, label = 0)
#' zVarArray(Water)
#'
#' # From a completed data set
#'
#' data(mdl) # matrix of limits of detection for Water
#' Water_multKM <- multKM(Water,label=0,dl=mdl) # nondetects imputation
#'
#' # Results split by two ficticious groups A and B
#' zVarArray(Water_multKM,groups=rep(c("A","B"),each=50))

zVarArray <- function(X, label = 0, groups = NULL, suppress.print = FALSE)
{

  if (is.vector(X))
    stop("X must be a matrix or data.frame class object")
  if (is.null(label))
    stop("A value for label must be given")
  if (!is.null(groups)){
    if (any(X == label, na.rm = T))
      stop(paste("Label", label, "was found in the data set. No zeros or unobserved values are allowed when a grouping factor is specified"))
  }
  if (!is.na(label)) {
    if (!any(X == label, na.rm = T) & (is.null(groups)))
      stop(paste("Label", label, "was not found in the data set"))
    if (label != 0 & any(X == 0, na.rm = T))
      stop("Zero values not labelled as such were found in the data set")
    if (any(is.na(X)))
      stop(paste(
        "NA values not labelled as zero values were found in the data set"
      ))
  }
  if (is.na(label)) {
    if (any(X == 0, na.rm = T))
      stop("Zero values not labelled as such were found in the data set")
    if (!any(is.na(X), na.rm = T) & (is.null(groups)))
      stop(paste("Label", label, "was not found in the data set"))
  }

  X <- as.data.frame(X)
  if (is.null(groups)){
    g <- zPatterns(X,label = label, plot = FALSE, suppress.print = TRUE)
    ifelse(is.na(label), X[is.na(X)] <- 0, X[X == label] <- 0)
  }
  else{
    g <- as.factor(groups)
    levNames <- levels(g)
    levels(g) <- 1:length(levels(g))
  }

  unobs <- sapply(split(X,g),function(x) sum(1*(x[1,]==0)),simplify = TRUE)
  if (any(unobs > (ncol(X)-2))) {
    warning("Some groups have less than two components available and NAs were produced (use zPatterns to check out)")
  }

  if (any(table(g) < 2)) warning("Some groups contain less than two observations and NAs were produced (use zPatterns to check out)")

  ni <- table(g); col <- ncol(X); nind <- nrow(X)
  numPat <- length(levels(g))
  p <- as.numeric(levels(g))
  pi <- ni / nind # % obs in each pattern

  VarArrByP <- array(NA, c(col, col, numPat)) # variation by pattern
  VarArrTot <- nvartot <- matrix(0, col, col) # overall variation matrix
  colnames(VarArrTot) <- rownames(VarArrTot) <- colnames(X)

  for (pat in 1:numPat)  {diag(VarArrByP[,,pat])<-rep(0,col)}

  # Variation array by pair of logratios
  for (di in 1:(col - 1)) {
    # by rows denominator di VariatMat

    Xdi <- X[X[, di] > 0, ] # subset of observed Xdi
    gdi <- g[X[, di] > 0] # subset pattern number of obs Xdi

    for (nj in (di + 1):col) {
      # by columns numerator nj VariatMat

      if (any(Xdi[, nj] > 0)) {
        Xnj <- Xdi[Xdi[, nj] > 0, c(di, nj)] # subset common di and nj
        gnj <- factor(gdi[Xdi[, nj] > 0]) # subset common patterns
        ngnj <- table(gnj)
        lxdinj <- lxdinjV <- log(Xnj[, 2] / Xnj[, 1])
        p <- as.numeric(levels(gnj))

        # by pattern
        # exp and var
        EbP <- -tapply(lxdinj, gnj, mean) #expectation by Pattern
        VbP <- tapply(lxdinj, gnj, var) #variance by Pattern
        EbP[is.na(VbP)] <- NA # NA for both if based on only one value

        VarArrByP[di, nj, p[!is.na(VbP)]] <-
          ((ngnj[!is.na(VbP)] - 1) * VbP[!is.na(VbP)]) / (ngnj[!is.na(VbP)])
        VarArrByP[nj, di, p[!is.na(EbP)]] <- EbP[!is.na(EbP)]

        # whole

        # center to zero each group for variance test
        numpatj <- nlevels(gnj)
        for (k in 1:numpatj) {
          lxdinjV[gnj == p[k]] <- scale(lxdinj[gnj == p[k]], TRUE, FALSE)
        }
        #
        VarArrTot[di, nj] <- as.numeric((length(lxdinjV) - 1) * var(lxdinjV) / length(lxdinjV))
        VarArrTot[nj, di] <- -as.numeric(mean(lxdinj))

      } # end if: control if common no-zero data data exists
    }# end for 2 by columns numerator nj VariatMat
  }# end for 1: by rows denominato di VariatMat

  # Formatting
  VarArrByP <- lapply(seq(dim(VarArrByP)[3]), function(x) VarArrByP[ , , x])
  for (i in 1:length(VarArrByP)){
    colnames(VarArrByP[[i]]) <- colnames(X)
    rownames(VarArrByP[[i]]) <- colnames(X)
  }
  if (is.null(groups)){
    names(VarArrByP) <- paste("Pattern",seq(length(VarArrByP)),sep="")
  }
  else{
    names(VarArrByP) <- levNames
  }
  VarArrByP[["Overall"]] <- VarArrTot

  result <- VarArrByP

  if (suppress.print == FALSE) {
    print(lapply(result,round,4))
    }

  invisible(result)

}
