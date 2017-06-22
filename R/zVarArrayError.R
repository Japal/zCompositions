#' @title Variation array relative error
#'
#' @description This function computes squared relative errors of variation arrays per group with respect to the overall variation array
#' based on observed data in a compositional data set. Groups can be defined by either zero/unobserved data patterns or
#' by a grouping factor in fully observed zero-free data sets.
#'
#' @details Squared relative errors (SRE) are calculated by confronting variation arrays (log-ratio variances and means) obtained per group and
#' the overall variation array based on observed data. Raw SREs are computed for each available pair-wise log-ratio. The weighted version uses
#' the corresponding group sizes to weight raw SREs. Total SRE is obtained as the sum of weighted SREs for each log-ratio. Further details by group are
#' provided by setting \code{breakdown = TRUE}.
#' 
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param label Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zero or unobserved data in \code{X} (\code{label = 0}, default).
#' @param groups Grouping factor in fully observed zero-free data sets (\code{groups = NULL}, default).
#' @param breakdown Logical value. Show results broken down by group (\code{breakdown = FALSE}, default).
#' @param suppress.print Suppress printed feedback (\code{suppress.print = FALSE}, default).
#' @return 1. SRE for each log-ratio variance and mean.
#' 2. Weighted SRE for each log-ratio variance and mean.
#' 3. Total SRE across log-ratio variances and means.
#' 4. Percentage contribution of each log-ratio to SRE in log-ratio variances and means.
#' If \code{breakdown = TRUE}:
#' 4. SREs per group.
#' 5. Weighted SREs per group.
#' 6. Percentage contribution of each group to total SRE.
#' 
#' @seealso \code{\link{zPatterns}}, \code{\link{zVarArray}}.
#' 
#' @examples
#' data(Water)
#' zPatterns(Water, label = 0)
#' zVarArrayError(Water)
#' zVarArrayError(Water, breakdown = TRUE)
#'
#' # From a completed data set
#'
#' data(mdl) # matrix of limits of detection for Water
#' Water_multKM <- multKM(Water,label=0,dl=mdl) # nondetects imputation
#'
#' # Results split by two ficticious groups A and B
#' zVarArrayError(Water_multKM,groups=rep(c("A","B"),each=50))

zVarArrayError <- function(X, label = 0, groups = NULL, breakdown = FALSE, suppress.print = FALSE)
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


  ni <- table(g); col <- ncol(X); nind <- nrow(X)
  numPat <- length(levels(g))
  p <- as.numeric(levels(g))
  pi <- ni / nind # % obs in each group

# Error matrices
ErrMat <- matrix(NA,col,col) # raw relative error
ErrPropMat<-matrix(NA,col,col) # relative error proportional to number of obs in group (pi)

# Raw relative error by group
ErrArr <- array(NA,c(col, col, numPat))
# Relative error proportional to pi by group
ErrPropArr <- array(NA,c(col, col, numPat))

# Names
rownames(ErrPropMat) <- colnames(ErrPropMat) <- colnames(X)
rownames(ErrMat) <- colnames(ErrMat) <-colnames(X)

## Call zVarArray
  VMZ <- zVarArray(X = X, label = label, groups = groups, suppress.print = TRUE)

  vErrExp <- vErrVar <- 0

## Squared relative errors for each pair-wise log-ratio

# Compute VARIATION ARRAY BY pair of components

  for (di in 1:(col-1)){ # by rows denominator di
    for (nj in (di+1):col){ # by columns numerator nj

  # Overall
  VarArrTotV <- VMZ$Overall[di,nj] # Var
  VarArrTotE <- VMZ$Overall[nj,di] # Exp

  # By group (values for the given log-ratio across groups)
  VbP <- sapply(VMZ[1:(length(VMZ)-1)],function(x) x[di,nj],simplify=TRUE) # Var
  EbP <- sapply(VMZ[1:(length(VMZ)-1)],function(x) x[nj,di],simplify=TRUE) # Exp

  # Raw relative error by group
  ErrArr[di,nj,p[!is.na(VbP)]] <- (1-(VbP[!is.na(VbP)]/VarArrTotV))^2 # Var
  ErrArr[nj,di,p[!is.na(EbP)]] <- (1-(EbP[!is.na(EbP)]/VarArrTotE))^2 # Exp
  
  # Raw relative error (the previous added up across groups)
  ErrMat[di,nj] <- sum(ErrArr[di,nj,p[!is.na(VbP)]]) # Var
  ErrMat[nj,di] <- sum(ErrArr[nj,di,p[!is.na(EbP)]]) # Exp

  # Relative error proportional to number of obs per group (pi)
    # By group
  ErrPropArr[di,nj,p[!is.na(VbP)]] <- pi[p[!is.na(VbP)]]*ErrArr[di,nj,p[!is.na(VbP)]]
  ErrPropArr[nj,di,p[!is.na(EbP)]] <- pi[p[!is.na(EbP)]]*ErrArr[nj,di,p[!is.na(EbP)]]
    # Summed up across groups
  ErrPropMat[di,nj] <- sum(ErrPropArr[di,nj,p[!is.na(VbP)]])
  ErrPropMat[nj,di] <- sum(ErrPropArr[nj,di,p[!is.na(EbP)]])

  # Total squared relative error for log-ratio vars (Var) and means (Exp) (based on weighted SRE)
  vErrVar <- vErrVar + ErrPropMat[di,nj]
  vErrExp <- vErrExp + ErrPropMat[nj,di]

    }# end for 2 by columns numerator nj VariatMat
  }# end for 1 by rows denominator di VariatMat
  
# Percentage contribution of each log-ratio var and mean
A <- ErrPropMat
A[upper.tri(A)] <- A[upper.tri(A)]/sum(A[upper.tri(A)])*100
A[lower.tri(A)] <- A[lower.tri(A)]/sum(A[lower.tri(A)])*100

# Percentage error of each group
sExp <- sVar <- rep(NA,numPat)
for (pat in 1:numPat){
  RL <- RU <- ErrPropArr[,,pat]
  # Var log-ratio
  RU[!upper.tri(ErrPropArr[,,pat])] <- NA
  # Exp log-ratio
  RL[!lower.tri(ErrPropArr[,,pat])] <- NA
  # Total
  if (sum(!is.na(RU))>0) {sVar[pat] <- sum(RU[!is.na(RU)])}
  if (sum(!is.na(RL))>0) {sExp[pat] <- sum(RL[!is.na(RL)])}
}# end for Pat

# In percentage
sVarPer <- sVar/sum(sVar,na.rm=T)*100; names(sVarPer) <- names(VMZ)[1:(length(VMZ)-1)]
sExpPer <- sExp/sum(sExp,na.rm=T)*100; names(sExpPer) <- names(VMZ)[1:(length(VMZ)-1)]

# Formatting ErrArr
ErrArr <- lapply(seq(dim(ErrArr)[3]), function(x) ErrArr[ , , x])
for (i in 1:length(ErrArr)){
  colnames(ErrArr[[i]]) <- colnames(X)
  rownames(ErrArr[[i]]) <- colnames(X)
}

names(ErrArr) <- names(VMZ)[1:(length(VMZ)-1)]

# Formatting ErrPropArr
ErrPropArr <- lapply(seq(dim(ErrPropArr)[3]), function(x) ErrPropArr[ , , x])
for (i in 1:length(ErrPropArr)){
  colnames(ErrPropArr[[i]]) <- colnames(X)
  rownames(ErrPropArr[[i]]) <- colnames(X)
}

names(ErrPropArr) <- names(VMZ)[1:(length(VMZ)-1)]


 if (breakdown==FALSE){
   result <-  list(SRE=ErrMat,
              WeightedSRE=ErrPropMat,
              TotalSREvars=vErrVar,
              TotalSREmeans=vErrExp,
              PercContribSRE=A)
 }
  else{
    result <- list(SRE=ErrMat,
                WeightedSRE=ErrPropMat,
                TotalSREvars=vErrVar,
                TotalSREmeans=vErrExp,
                PercContribSRE=A,
                SREbyGroup=ErrArr,
                WeightedSREbyGroup=ErrPropArr,
                PercContribTotalSREbyGroupVars=sVarPer,
                PercContribTotalSREbyGroupMeans=sExpPer)
  }
    
if (suppress.print == FALSE) {
  print(lapply(result,function(x){
   if (class(x)=="list") {lapply(x,round,4)}
    else {round(x,4)}
  }))
}

invisible(result)
         
}