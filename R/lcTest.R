#' @title Log-contrast homogeneity test
#'
#' @description This function tests for homogeneity across groups of means and variances of
#' user-defined log-contrasts. Groups can be defined by either zero/unobserved data patterns or by a grouping
#' factor in fully observed zero-free data sets.
#'
#' @details Homogeneity of log-contrast means and variances across groups is tested using either parametric or non-parametric tests. When
#' \code{method = "parametric"}, ordinary analysis of variance and Bartlett's tests are used. Alternatively,
#' Kruskal-Wallis and Fligner-Killen tests are used instead when \code{method = "nonparametric"}. The results of a permutation test of homogeneity of variation
#' arrays based on total weighted squared relative errors are also provided (see \code{\link{zVarArrayTest}} for more details).
#' The log-contrast is specified by the \code{lc} argument using a vector of codes 1, -1 and 0 for components
#' in the numerator, denominator and omitted respectively.
#'
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param label Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zero or unobserved data in \code{X} (\code{label = 0}, default).
#' @param groups Grouping factor in fully observed zero-free data sets (\code{groups = NULL}, default).
#' @param lc User-defined log-contrast (see details below).
#' @param method Approach used for mean and variance homogeneity testing (\code{method = "parametric"}, default).
#' @param b Number of bootstrap resamples used by permutation test (\code{b = 1000}, default).
#'
#' @return Test p-values for log-contrast means and variances.
#'
#' @seealso \code{\link{zPatterns}}, \code{\link{zVarArray}}, \code{\link{zVarArrayError}}.
#'
#' @examples
#' data(Water)
#' zPatterns(Water, label = 0)
#'
#' # Test of homogeneity in log-contrast Potassium/Arsenic*Calcium
#' lcTest(Water, label = 0, lc = c(1,-1,-1,0))

lcTest <- function(X, label = 0, groups = NULL, lc = NULL, method = c("parametric", "nonparametric"), b = 1000){

  if (any(X<0, na.rm=T)) stop("X contains negative values")
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

  if(is.null(lc)){stop("A sensible log-contrast must be specified to use this function")}
  if(length(lc) != ncol(X)) {
    stop("The number of columns in X and lc do not agree")}
  if((all(lc >= 0)) | (all(lc <= 0))) {
    stop(paste("The log-contrast is not correctly defined"))}

  method <- match.arg(method)
  
  X <- as.data.frame(X,stringsAsFactors=TRUE)
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

  ni <- table(g); col <- ncol(X); nind <- nrow(X)
  numPat <- length(levels(g))
  p <- as.numeric(levels(g))
  pi <- ni / nind # % obs in each group

  usepart <- (lc!=0)
  usePat <- vector("numeric")
  Xfeas <- vector("numeric")
  nifeas <- pifeas <- gfeas <- vector("numeric")
  for (pat in 1:numPat){
    x <- X[g == pat,usepart]

    if (all(x[1,] > 0)){
      usePat <- cbind(usePat,pat)
      Xfeas <- rbind(Xfeas,x)
      gfeas <- c(gfeas,g[g==pat])
      nifeas <- c(nifeas,ni[pat])
      pifeas <- c(pifeas,pi[pat])
      }#end if
  }# end for
  gfeas <- as.factor(gfeas)
  pfeas <- as.numeric(levels(gfeas))

  Apval <- Bpval <- pvalExp <- pvalVar <- NA

  # check if there is a pattern
  if (length(usePat) > 0){
    # log-contrasts of data
    bal <- rep(0,ncol(Xfeas))
    den <- sum(lc==-1)
    num <- sum(lc==1)
    bal[lc[usepart]==1] <- sqrt(den/((den+num)*num))
    bal[lc[usepart]==-1] <- -sqrt(num/((den+num)*den))
    Yfeas <- YfeasV <- as.matrix(log(Xfeas))%*%bal

    # center to zero each group for variance test
    numpatj <- nlevels(gfeas)
    for (k in 1:numpatj){
      YfeasV[gfeas==pfeas[k]] <- scale(Yfeas[gfeas==pfeas[k]],TRUE,FALSE)
    }

    if (length(usePat) > 1){
      if (method == "parametric"){
        # anova test of location
        if (min(nifeas)>1) {Apval<-anova(lm(Yfeas~gfeas))$"Pr(>F)"[1]}
        # Bartlett test of variances
        if (min(nifeas)>1)  {Bpval<-bartlett.test(YfeasV~gfeas)$p.value}
        }
      if (method == "nonparametric"){
        # Kruskal-Wallis test of location
        if (min(nifeas)>1) {Apval <- kruskal.test(Yfeas ~ gfeas)$p.value}
        # Fligner-Killeen test of variances
        if (min(nifeas)>1)  {Bpval <- fligner.test(YfeasV ~ gfeas)$p.value}
        }

    if (is.na(Apval) | is.na(Bpval)) warning("Log-contrast present in groups of < 2 observations and NAs produced (use zPatterns to check out)")

      # Permutation test
      ErrVar <- ErrExp<-rep(0,b)
      for (rept in 1:b){

        lx <- sample(Yfeas)
        lxV <- sample(YfeasV)
        # exp and var
        EbP <- tapply(lx,gfeas,mean) #expectation by Pattern
        VbP <- tapply(lxV,gfeas,var) #variance by Pattern
        EbP[is.na(VbP)] <- NA # NA for both if based on only one value

        # By pattern
        VarByP <- ((nifeas[!is.na(VbP)]-1)*VbP[!is.na(VbP)])/(nifeas[!is.na(VbP)])
        ExpByP <- EbP[!is.na(EbP)]

        # Overall
        VarTot <- as.numeric((length(YfeasV)-1)*var(YfeasV)/length(YfeasV)) # Var
        ExpTot <- as.numeric(mean(Yfeas)) # Exp
        # squared relative error
        ErrVar[rept] <- sum((pi[pfeas[!is.na(VbP)]])*((1-(VarByP[!is.na(VbP)]/VarTot))^2))
        ErrExp[rept] <- sum((pi[pfeas[!is.na(EbP)]])*((1-(ExpByP[!is.na(EbP)]/ExpTot))^2))
      }# end for Permutation test



# Original errors
EbP <- tapply(Yfeas,gfeas,mean) #expectation by Pattern
VbP <- tapply(YfeasV,gfeas,var) #variance by Pattern
EbP[is.na(VbP)] <- NA # NA for both if based on only one value

# By pattern
VarByP <- ((nifeas[!is.na(VbP)]-1)*VbP[!is.na(VbP)])/(nifeas[!is.na(VbP)])
ExpByP <- EbP[!is.na(EbP)]

# Overall
VarTot <- as.numeric((length(YfeasV)-1)*var(YfeasV)/length(YfeasV)) # Var
ExpTot <- as.numeric(mean(Yfeas)) # Exp

# Squared relative error
ErrVarOr <- sum((pi[pfeas[!is.na(VbP)]])*((1-(VarByP/VarTot))^2))
ErrExpOr <- sum((pi[pfeas[!is.na(EbP)]])*((1-(ExpByP/ExpTot))^2))

# p-value (add 1 for the original sample, it is included)
pvalExp <- (sum(ErrExp>ErrExpOr)+1)/(b+1)
pvalVar <- (sum(ErrVar>ErrVarOr)+1)/(b+1)

  } else{stop("Log-contrast only available for one pattern or group")}
} else{stop("Log-contrast not available for any pattern or group")}

  lcstringn <- vector("character")
  lcstringd <- vector("character")
  lcnamesn <- names(X)[which(lc==1)]
  for (i in 1:length(lcnamesn)) lcstringn <- paste(lcstringn,lcnamesn[i],sep="*")
  lcstringn <- substring(lcstringn,2,nchar(lcstringn))
  lcnamesd <- names(X)[which(lc==-1)]
  for (i in 1:length(lcnamesd)) lcstringd <- paste(lcstringd,lcnamesd[i],sep="*")
  lcstringd <- substring(lcstringd,2,nchar(lcstringd))

  cat("\n")
  cat("Log-contrast homogeneity tests \n")
  cat("------------------------------ \n")
  cat(paste("Number of groups:",nlevels(g),"\n"))
  cat(paste("Log-contrast: ","(",lcstringn,")"," / ","(",lcstringd,")","\n",sep=""))
  if (method == "nonparametric"){
    cat(paste("Kruskal-Wallis test of log-contrast means:",round(Apval,4),"\n"))
    cat(paste("Fligner-Killeen test of log-contrast variances:",round(Bpval,4),"\n"))
  }
  if (method == "parametric"){
    cat(paste("ANOVA test of log-contrast means:",round(Apval,4),"\n"))
    cat(paste("Bartlett test of log-contrast variances:",round(Bpval,4),"\n"))
  }
  cat(paste("Permutation test of total weighted SRE in log-contrast means:",round(pvalExp,4),"\n"))
  cat(paste("Permutation test of total weighted SRE in log-contrast variances:",round(pvalVar,4),"\n"))

}
