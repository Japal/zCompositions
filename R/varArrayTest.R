#' @title Variation array homogeneity test
#'
#' @description This function performs a permutation test of the homogeneity of group-wise and overall variation arrays from all
#' pair-wise log-ratios in a compositional data set. Groups can be defined by either zero/unobserved data patterns or by a grouping
#' factor in fully observed zero-free data sets.
#'
#' @details The permutation test of homogeneity is based on total weighted squared relative errors (SRE) reflecting on divergence
#' between group-wise variation arrays and overall (see \code{\link{varArrayError}} and
#' \code{\link{varArray}} for more details). Note that for groups including less than two observations SRE is set to NA.
#' 
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param label Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zero or unobserved data in \code{X} (\code{label = 0}, default).
#' @param groups Grouping factor in fully observed zero-free data sets (\code{groups = NULL}, default).
#' @param b Number of bootstrap resamples used (\code{b = 1000}, default).
#'
#' @return Test p-values for log-ratio variances and means.
#' 
#' @seealso \code{\link{zPatterns}}, \code{\link{varArray}}, \code{\link{varArrayError}}.
#'
#' @examples
#' data(Water)
#' zPatterns(Water, label = 0)
#' varArrayTest(Water)

varArrayTest <- function(X, label = 0, groups = NULL, b = 1000){
  
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
  
  # Weighted SRE for the original data set
  resE <- varArrayError(X, label = label, groups = groups, suppress.print = TRUE)
  
  vErrExp <- vErrVar <- rep(0,b)
  
  ## Weighted SRE for EACH LOG-RATIO  
    # VARIATION ARRAY BY pair of logratios
  for (di in 1:(col-1)){ # by rows denominator di VariatMat

    Xdi <- X[X[,di]>0,] # subset of observed Xdi
    gdi <- g[X[,di]>0] # subset pattern number of obs Xdi
    
    for (nj in (di+1):col){ # by columns numerator nj VariatMat
  
      vErrExpij <- vErrVarij <- rep(0,b)
      
      if (any(Xdi[,nj]>0)){
        Xnj <- Xdi[Xdi[,nj]>0,c(di,nj)] # subset common di and nj
        gnj <- factor(gdi[Xdi[,nj]>0])
        ngnj <- table(gnj)
        p <- as.numeric(levels(gnj))
        lxdinjOr <- lxdinjOrV <- log(Xnj[,2]/Xnj[,1])
       
        # center to zero each group for variance test
        numpatj <- nlevels(gnj)
        for (k in 1:numpatj){
          lxdinjOrV[gnj==p[k]] <- scale(lxdinjOr[gnj==p[k]],TRUE,FALSE)
        }
      # PERMUTATION TEST
      for (rept in 1:b){
        
        lxdinj <- sample(lxdinjOr)
        lxdinjV <- sample(lxdinjOrV)
        # var and exp
        VbP <- tapply(lxdinjV,gnj,var) #variance by group
        EbP <- -tapply(lxdinj,gnj,mean) #expectation by group
        EbP[is.na(VbP)] <- NA # NA for both if based on only one value
        
       # By group
        VarArrByPV <- ((ngnj[!is.na(VbP)]-1)*VbP[!is.na(VbP)])/(ngnj[!is.na(VbP)])
        VarArrByPE <- EbP[!is.na(EbP)]
       
       # Overall
       VarArrTotV <- as.numeric((length(lxdinjOrV)-1)*var(lxdinjOrV)/length(lxdinjOrV)) # Var
       VarArrTotE <- -as.numeric(mean(lxdinjOr)) # Exp
       
        # Weighted squared relative error
        vErrVarij[rept] <- sum((pi[p[!is.na(VbP)]])*((1-(VarArrByPV/VarArrTotV))^2)) 
        vErrExpij[rept] <- sum((pi[p[!is.na(EbP)]])*((1-(VarArrByPE/VarArrTotE))^2)) 
       
       }# end of Permutation test
     
      vErrExp <- vErrExp + vErrExpij
      vErrVar <- vErrVar + vErrVarij
      } # end if: control if common no-zero data data exists  
      
    }# end for 2 by columns numerator nj VariatMat
  } # end for 1 by rows denominator di VariatMat
  
  # p-value (add 1 for the original sample, it is included)
  pvalExp <- (sum(vErrExp >= resE$TotalSREmeans)+1)/(b+1)
  pvalVar <- (sum(vErrVar >= resE$TotalSREvars)+1)/(b+1)
  
  cat("\n")
  cat("Variance array homogeneity test \n")
  cat("------------------------------- \n")
  cat(paste("Number of groups:",nlevels(g),"\n"))
  cat(paste("P-value for homogeneity of log-ratio variances:",round(pvalVar,4),"\n"))
  cat(paste("P-value for homogeneity of log-ratio means:",round(pvalExp,4),"\n"))
}