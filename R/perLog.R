#' Test of differences in group means for compositional data
#'
#' @description A nonparametric permutation test based on pairwise logratios to assess the hypothesis of equality of means/locations between subsets of observations according to an externally or internally defined factor. If any, zero patterns are considered as default internal grouping factor.
#'
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param groups Factor variable indicating the grouping structure. If NULL (default), any zero patterns in the data will be used as internal grouping factor.
#' 
#' Note that if a grouping factor is set by the user, then any zeros in the data must be previously dealt with, e.g. by imputation.
#' @param p Power parameter used in overall dissimilarity test statistic. (default = 10).
#' @param alpha Significance level parameter (default = 0.05).
#' @param R Number of permutation resamples (default = 1000).
#' @param posthoc.g Logical. If TRUE, performs post-hoc analysis for pairs of groups (default = FALSE).
#' @param posthoc.lr Logical. If TRUE, performs post-hoc analysis for logratios (default = FALSE).
#' @param mAdj Adjustment of p-values for multiple post-hoc testing (see \code{\link{p.adjust}}). Default is Benjamini and Hochberg's FDR method (default = "BH").
#'
#' @return A list object of class "printres" containing summaries of results:
#' \item{disOv}{Overall dissimilarity test statistic.}
#' \item{pvalOv}{Overall permutation p-value.}
#' \item{posthoc.groups}{If `posthoc.g = TRUE` and the main test is significant, results of the post-hoc analysis for pairs of groups.}
#' \item{posthoc.logratios}{If `posthoc.lr = TRUE` and the main test is significant, results of the post-hoc analysis for pairwise logratios.}
#' \item{wts}{Welch's t-statistic.}
#' \item{disEl}{Pairwise logratio elemental dissimilarity.}
#' \item{rcElBg}{Relative contribution of elemental dissimilarity to between-group dissimilarity.}
#' \item{rcElOv}{Relative contribution of elemental dissimilarity to overall dissimilarity.}
#' \item{pvalElAdj}{Adjusted p-value in post-hoc comparison.}
#' \item{parts0}{If `groups = NULL`, list containing original names of zero parts in the respective zero patterns.}
#' @details The test relies on the unique pairwise logratios between parts of the given composition. It assesses whether the observed overall dissimilarity is significantly different from that expected under the null hypothesis of equal group means. If so, it can perform post-hoc analyses by pairs of groups and pairwise logratios, evaluating their relative contributions to dissimilarity at group and overall levels. The p-values in post-hoc testing are adjusted for multiple comparisons using the specified method.
#' In the case of internal grouping defined by zero patterns, strings of binary codes are used to label each pattern in the output, with 0 indicating no zero part and 1 indicating zero part.
#' @export
#' @seealso \code{\link{zPatterns}}
#' @examples
#' # Load the Water data set 
#' data(Water)
#' # Visualise zero patterns and select the first three for illustration
#' tmp <- zPatterns(Water, label = 0)
#' Water2 <- Water[tmp %in% c(1,2,3),]
#' # Test overall differences by zero pattern on the selected data set
#' zPatterns(Water2, label = 0)
#' perLog(Water2)

perLog <- function(X, groups = NULL, p = 10, alpha = 0.05, R = 1000, 
                      posthoc.g = FALSE, posthoc.lr = FALSE, mAdj = "BH") {
  
  ## Auxiliary functions
  
  # Welch's t-statistic for 2 groups of LR
  wtstat <- function(lr1, lr2){
    lr1 <- lr1[!is.na(lr1)]
    lr2 <- lr2[!is.na(lr2)]
    n1 <- length(lr1)
    n2 <- length(lr2)
    m1 <- mean(lr1)
    m2 <- mean(lr2)
    nmr <- m1-m2
    v1 <- var(lr1)
    v2 <- var(lr2)
    s1 <- v1/n1
    s2 <- v2/n2
    dnm <- sqrt(s1+s2)
    wts <- nmr/dnm
    df1i <- 1/(n1-1)
    df2i <- 1/(n2-1)
    df <- dnm^4/(df1i*s1^2+df2i*s2^2)
    return(list(wts = wts, df = df))
  }
  
  # Elemental dissimilarities for a pair of groups
  disElPair <- function(LR1, LR2, alpha = 0.05, rwts = FALSE, pgNAlri = NULL){
    ll <- length(pgNAlri)
    LR <- rbind(LR1, LR2)
    if (ll!=0) { # case with an internal factor
      Dd <- ncol(LR)
      LR <- LR[, -pgNAlri]
    }
    if (length(LR) == 0) { # no of the D*(D-1)/2 LR available
      disEP <- wts <- rep(NA, Dd)
    } else if (is.null(ncol(LR))) { # 1 of the D*(D-1)/2 LR available
      n1 <- nrow(LR1)
      n <- length(LR)
      Wdf <- wtstat(LR[1:n1], LR[(n1+1):n])
      wts <- Wdf[[1]]
      awts <- abs(wts)
      df <- Wdf[[2]]
      qN <- qnorm(pt(awts, df))
      if (is.infinite(qN)) qN <- qnorm(0.9999999999999999) # replacing extreme quantiles by the maximal non-Inf obtainable value
      qNs <- qnorm(1-alpha/2)
      disEP <- qN/qNs
    } else { # more than 1 of the D*(D-1)/2 LR available
      n1 <- nrow(LR1)
      n <- nrow(LR)
      ind1 <- 1:n1
      ind2 <- (n1+1):n
      Wdf <- apply(LR, 2, function(x) wtstat(x[ind1], x[ind2]))
      wts <- unlist(sapply(Wdf, function(x) x[1]))
      awts <- abs(wts)
      df <- unlist(sapply(Wdf, function(x) x[2]))
      qN <- qnorm(pt(awts, df))
      if(is.infinite(sum(qN))) qN[is.infinite(qN)] <- qnorm(0.9999999999999999) # replacing extreme quantiles by the maximal non-Inf obtainable value
      qNs <- qnorm(1-alpha/2)
      disEP <- qN/qNs
    }
    if (ll!=0) { # case with an internal factor
      disEP0 <- disEP
      disEP <- rep(NA, Dd)
      disEP[-pgNAlri] <- disEP0
      if (rwts == TRUE) { # for the original data
        wts0 <- wts
        wts <- rep(NA, Dd)
        wts[-pgNAlri] <- wts0
      }
    }
    if (rwts == TRUE) { # for the original data
      return(list(disEP = disEP, wts = wts))
    } else{ # for the permuted data
      return(disEP = disEP)
    }
  }
  
  # Elemental dissimilarities for all pairs of groups
  disElAll <- function(LR, groups, alpha = 0.05, rwts = FALSE, pgNAlri = NULL){
    Dd <- ncol(LR)
    cnlr <- colnames(LR)
    cmbg <- combn(levels(groups), 2, simplify = FALSE)
    Gg <- length(cmbg)
    grcmb <- sapply(cmbg, function(x) paste(x, collapse = " vs "))
    ll <- length(pgNAlri)
    if (ll!=0) { # case with an internal factor
      if (rwts == TRUE) { # for the original data
        disEAw <- lapply(1:Gg, function(x) disElPair(LR[groups==cmbg[[x]][1], ], LR[groups==cmbg[[x]][2], ],
                                                     alpha = alpha, rwts = TRUE, pgNAlri = pgNAlri[[x]]))
        disEA <- matrix(unlist(sapply(disEAw, function(x) x[1])), length(cmbg), Dd, byrow = TRUE)
        wts <- matrix(unlist(sapply(disEAw, function(x) x[2])), length(cmbg), Dd, byrow = TRUE)
        rownames(disEA) <- rownames(wts) <- grcmb
        colnames(disEA) <- colnames(wts) <- cnlr
        return(list(disEA = disEA, wts = wts))
      } else { # for the permuted data
        disEA <- t(sapply(1:Gg, function(x) disElPair(LR[groups==cmbg[[x]][1], ], LR[groups==cmbg[[x]][2], ],
                                                      alpha = alpha, pgNAlri = pgNAlri[[x]])))
        rownames(disEA) <- grcmb
        colnames(disEA) <- cnlr
        return(disEA = disEA)
      }
    } else { # case with an external factor
      if (rwts == TRUE) { # for original data
        disEAw <- lapply(cmbg, function(x) disElPair(LR[groups==x[1], ], LR[groups==x[2], ],
                                                     alpha = alpha, rwts = TRUE))
        disEA <- matrix(unlist(sapply(disEAw, function(x) x[1])), length(cmbg), Dd, byrow = TRUE)
        wts <- matrix(unlist(sapply(disEAw, function(x) x[2])), length(cmbg), Dd, byrow = TRUE)
        rownames(disEA) <- rownames(wts) <- grcmb
        colnames(disEA) <- colnames(wts) <- cnlr
        return(list(disEA = disEA, wts = wts))
      } else { # for the permuted data
        disEA <- t(sapply(cmbg, function(x) disElPair(LR[groups==x[1], ], LR[groups==x[2], ],
                                                      alpha = alpha)))
        rownames(disEA) <- grcmb
        colnames(disEA) <- cnlr
        return(disEA = disEA)
      }
    }
  }
  
  # Elemental dissimilarities for all pairs of groups in a permuted dataset
  disElAllPerm <- function(LR, groups, alpha = 0.05, pgNAlri = NULL, gRlri = NULL, minAR = 5){
    lev <- levels(groups)
    n <- nrow(LR)
    indR <- sample(1:n, n, replace = FALSE)
    LRR <- LR[indR, ]
    ll <- length(gRlri)
    if (ll!=0) { # case with an internal factor
      nmbA <- lapply(1:ll, function(x) apply(data.frame(LRR[groups==lev[x], gRlri[[x]]]), 2,
                                             function(y) length(y[!is.na(y)]))) # i-th element ... numbers of available required LR in the i-th group
      minA <- min(unlist(nmbA))
      while (minA < minAR) {
        indR <- sample(1:n, n, replace = FALSE)
        LRR <- LR[indR, ]
        nmbA <- lapply(1:ll, function(x) apply(data.frame(LRR[groups==lev[x], gRlri[[x]]]), 2,
                                               function(y) length(y[!is.na(y)])))
        minA <- min(unlist(nmbA))
      }
      disEAP <- disElAll(LRR, groups, alpha = alpha, pgNAlri = pgNAlri)
    } else { # case with an external factor
      disEAP <- disElAll(LRR, groups, alpha = alpha)
    }
    return(disEAP)
  }
  
  # Selection of results to print, only overall dissimilarity and p-value by default
  print.printres <- function(x, ...) {
    cat(paste("\n","Overall dissimilarity:","\n\n","Test statistic E = ", round(x$disOv,4), 
              ", p-value = ", round(x$pvalOv,4),sep=""))
    #         wts ... a matrix (of size G*(G-1)/2 x D*(D-1)/2, where rows correspond to pairs of groups and columns to LR) 
    #                 containing the Welch's t-statistics (t_{ij}^{(k, l)}) in the original data set
    #         disEl ... a matrix (of size G*(G-1)/2 x D*(D-1)/2, where rows correspond to pairs of groups and columns to LR) 
    #                   containing the elemental dissimilarities (e_{ij}^{(k, l)}) in the original data set
    #         disEl_Prm ... a list (of length R, where each item corresponds to different permuted data set)
    #                       containing matrices (of size G*(G-1)/2 x D*(D-1)/2, where rows correspond to pairs of groups and columns to LR) 
    #                       containing the elemental dissimilarities in the respective permuted data set
    #         parts0 ... if groups = NULL, a list (of length G, where each item corresponds to different zero pattern)
    #                    containing names of zero parts in the respective zero patterns
  cat("\n")
    if (!is.null(x$posthoc.groups)) {
      cat("\n")
      cat("Post-hoc analysis by pairs of groups:\n")
      cat("\n")
      print(x$posthoc.groups$summary)
      if (!is.null(x$posthoc.groups$parts0)) {
        cat(paste("\n","Zero patterns and zero parts therein:","\n\n",sep=""))
        print(x$posthoc.groups$parts0)
      }
    }
    if (!is.null(x$posthoc.logratios)) {
      cat("\n")
      cat("Post-hoc analysis for pairwise logratios:\n")
      for (i in seq_along(x$posthoc.logratios$summary)) {
        cat(paste("\n","Pair:", names(x$posthoc.logratios$summary)[i], "\n"))
        print(round(x$posthoc.logratios$summary[[i]], 4))
      }
    }
    
  }
  
  # Post-hoc analysis by pairs of groups
  perlogPHg <- function(perlogR, p = 10, alpha = 0.05, mAdj = "BH"){
    
    # INPUT: perlogR ... output of the perlog function
    #        p ... a power p for computing  between-groups dissimilarities from elemental dissimilarities; default is 10
    #        alpha ... significance level; default is 0.05
    #        mAdj ... a method for adjusting p-values - options as in p.adjust function; default is "BH" (Benjamini & Hochberg)
    # OUTPUT: summary ... a matrix (of size G*(G-1)/2 x 3) with columns corresponding to different results:
    #                     disBg ... the between-groups dissimilarities (E^{(kl)}) in the original data set 
    #                     rcBgOv ... the relative contributions (in %) of the between-groups dissimilarities to the overall dissimilarity (100*E^{(kl)}/E)
    #                     pvalBgAdj ... the adjusted between-groups p-values 
    #         parts0 ... if zero parts ar considered for grouping, 
    #                    a list (of length G, where each item corresponds to different zero pattern)
    #                    containing names of zero parts in the respective zero patterns
    
    disEl <- perlogR$disEl
    Gg <- nrow(disEl)
    if (Gg == 1) { # case with 2 groups
      summary <- c("Bet-Grp diss" = perlogR$disOv, 
                  "%Ctb overall diss" = 100, 
                  "Adj p-value" =  perlogR$pvalOv)
      if (is.null(perlogR$parts0)) {
        res <- list(summary = summary)
      } else {
        res <- list(summary = summary,
                   parts0 = perlogR$parts0)
      }
    } else { # case with more than 2 groups
      disEl_Prm <- perlogR$disEl_Prm
      disBg <- rowSums(disEl^p, na.rm = T)
      rcBgOv <- 100*disBg/sum(disBg, na.rm = T)
      disBg_Prm <- t(sapply(disEl_Prm , function(x) rowSums(x^p, na.rm = T)))
      disBg_All <- rbind(disBg, disBg_Prm)
      pvalBg <- apply(disBg_All, 2, function(x) mean(x[-1] >= x[1]))
      pvalBgAdj <- p.adjust(pvalBg, method = mAdj)
      summary <- cbind("Bet-Grp diss" = round(disBg, 4),
                      "%Ctb overall diss" = round(rcBgOv, 4),
                      "Adj p-value" = round(pvalBgAdj, 4))
      if (is.null(perlogR$parts0)) {
        res <- list(summary = summary)
      } else {
        res <- list(summary = summary,
                   parts0 = perlogR$parts0)
      }
    }
    return(res)
  }
  
  # Post-hoc analysis for the pairwise logratios for significantly different pairs of groups
  perlogPHlr <- function(perlogR, p = 10, alpha = 0.05, mAdj = "BH"){
    
    # INPUT: perlogR ... output of the perlog function
    #        p ... a power p for computing  between-groups and overall dissimilarities from elemental dissimilarities; default is 10
    #        alpha ... significance level; default is 0.05
    #        mAdj ... a method for adjusting p-values - options as in p.adjust function; default is "BH" (Benjamini & Hochberg)
    # OUTPUT: summary ... a list (of size equal to the number of the significantly different pairs of groups) 
    #                     with matrices (of size D*(D-1)/2 x 5) with columns corresponding to different results:
    #                     wts ... the Welch's t-statistics (t_{ij}^{(k, l)}) in the original data set
    #                     disEl ... the elemental dissimilarities (e_{ij}^{(k, l)}) in the original data set
    #                     rcElBg ... the relative contributions (in %) of the elemental dissimilarities to the between-group dissimilarities (100*e_{ij}^{(kl)}/E^{(kl)})
    #                     rcElOv ... the relative contributions (in %) of the elemental dissimilarities to the overall dissimilarity (100*e_{ij}^{(kl)}/E)
    #                     pvalElAdj ... the adjusted elemental p-values
    #         significance ... a list (of size equal to the number of the significantly different pairs of groups)
    #                          with matrices (of size D x D, with rows corresponding to numerator and columns to denominator of logratio) 
    #                          indicating (non)significance of the logratios together with the direction of the significance:
    #                          0 ... nonsignificant logratio
    #                          1 ... significant logratio in positive direction 
    #                          -1 ... significant logratio in negative direction 
    #         parts0 ... if zero parts ar considered for grouping, 
    #                    a list (of length G, where each item corresponds to different zero pattern)
    #                    containing names of zero parts in the respective zero patterns
    
    disEl <- perlogR$disEl
    pr <- rownames(disEl)
    npr <- nrow(disEl)
    if (npr==1) { # case with 2 groups
      wts <- perlogR$wts[1,]
      disEl <- disEl[1,]
      disEl.p <- disEl^p
      rcElBg <- rcElOv <- c(100*disEl.p/perlogR$disOv)
      disEl_Prm <- perlogR$disEl_Prm
      disEl.p_Prm <- t(sapply(disEl_Prm, function(x) x^p))
      disEl.p_All <- rbind(disEl.p, disEl.p_Prm)
      pvalEl <- apply(disEl.p_All, 2, function(x) mean(x[-1] >= x[1]))
      pvalElAdj <- p.adjust(pvalEl, method = mAdj)
      summary <- list(cbind("wts" = wts,
                           "Elem diss" = disEl,
                           "%Ctb bet-grp diss" = rcElBg,
                           "%Ctb overall diss" = rcElOv,
                           "Adj p-value" = pvalElAdj))
      names(summary) <- pr
    } else { # case with more than 2 groups
      perlogPHgR <- perlogPHg(perlogR)
      iS <- which(perlogPHgR$summary[, 3]<alpha)
      npr <- length(iS)
      if (npr==0) stop("No significantly different pair.")
      pr <- pr[iS]
      wts <- perlogR$wts[iS, ]
      disEl <- disEl[iS, ]
      disEl.p <- disEl^p
      rcElBg <- 100*disEl.p/perlogPHgR$summary[iS, 1]
      rcElOv <- 100*disEl.p/perlogR$disOv
      disEl_Prm <- lapply(perlogR$disEl_Prm, function(x) x[iS, ])
      if (npr==1) { # one significantly different pair
        disEl.p_Prm <- t(sapply(disEl_Prm, function(x) x^p))
        disEl.p_All <- rbind(disEl.p, disEl.p_Prm)
        pvalEl <- apply(disEl.p_All, 2, function(x) mean(x[-1] >= x[1]))
        pvalElAdj <- p.adjust(pvalEl, method = mAdj)
        summary <- list(cbind("wts" = wts,
                             "Elem diss" = disEl,
                             "%Ctb bet-grp diss" = rcElBg,
                             "%Ctb overall diss" = rcElOv,
                             "Adj p-value" = pvalElAdj))
        names(summary) <- pr
      } else { # more than one significantly different pairs
        disEl.p_Prm <- lapply(disEl_Prm, function(x) x^p)
        disEl.p_All <- c(list(disEl.p), disEl.p_Prm)
        disEl.p_All <- lapply(1:npr, function(i) t(sapply(disEl.p_All, function(x) x[i, ])))
        pvalEl <- t(sapply(disEl.p_All, function(y) apply(y, 2, function(x) mean(x[-1] >= x[1]))))
        rownames(pvalEl) <- pr
        pvalElAdj <- pvalEl
        pvalElAdj[] <- p.adjust(c(pvalEl), method = mAdj)
        summary <- lapply(1:npr, function(i) cbind("wts" = wts[i, ], 
                                                  "Elem diss" = disEl[i, ], 
                                                  "%Ctb bet-grp diss" = rcElBg[i, ], 
                                                  "%Ctb overall diss" = rcElOv[i, ], 
                                                  "Adj p-value" = pvalElAdj[i, ]))
        }
    } 
    names(summary) <- pr
    significance <- lapply(summary, function(x) {
      cnlr <- rownames(x)
      Dd <- length(cnlr)
      D <- (1+sqrt(1+8*Dd))/2
      cn <- unique(unlist(strsplit(cnlr, "/")))
      mt <- matrix(NA, D, D)
      sgn <- rep(0, Dd)
      sgn[which(x[,5]<alpha)] <- 1
      sgn <- sign(x[,1])*sgn
      mt[lower.tri(mt)] <- -sgn
      mtt <- -t(mt)
      mt[upper.tri(mt)] <- mtt[upper.tri(mtt)]
      rownames(mt) <- colnames(mt) <- cn
      return(mt)
    })
    if (is.null(perlogR$parts0)) {
      res <- list(summary = summary,
                 significance = significance)
    } else {
      res <- list(summary = summary,
                 significance = significance,
                 parts0 = perlogR$parts0)
    }
    return(res)
  }
  
  ################################################################################################
  
  # Main function body
  
  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a matrix or data.frame class object")
  if (any(X==0, na.rm = TRUE) & (!is.null(groups))) stop("User-defined factor set but zero values found in the data set")
  if (any(X < 0, na.rm = T)) stop("X contains negative values")
  if (any(is.na(X))) stop("X contains missing values")
  groups0 <- groups
  if (is.null(groups0)) { # case with an internal factor
    zr <- 1*(X==0)
    groups <- apply(zr, 1, function(x) paste(x, collapse = "")) # 1 signifies zero value, 0 nonzero value
    zrg <- apply(do.call(rbind, strsplit(sort(unique(groups)), split = "")), 2, as.integer)
  }
  groups <- as.factor(groups)
  nk <- table(groups)
  minnk <- min(nk)
  minAR <- min(minnk, 5) # specifying the minimal number of logratios that must be present for the required combinations (k, l, i, j) in a permuted data set
  if (minnk<5) warning("At least 1 group has less than 5 observations")
  lev <- levels(groups)
  G <- length(lev)
  cmbg <- combn(lev, 2, simplify = FALSE) # G*(G-1) combinations of groups
  Gg <- length(cmbg)
  cn <- colnames(X)
  D <- ncol(X)
  if (is.null(cn)) { # naming compositional parts if no names are provided
    cn <- paste("P", 1:D, sep = "")
    colnames(X) <- cn
  }
  cmbp <- combn(cn, 2, simplify = FALSE) # D*(D-1) combinations of compositional parts
  cnlr <- sapply(cmbp, function(j) paste(j, collapse = "/"))
  LR <- sapply(cmbp, function(j, x) log(x[, j[1]]/x[, j[2]]), x = X)
  LR[is.infinite(LR)] <- NA
  LR[is.nan(LR)] <- NA
  colnames(LR) <- cnlr
  if (is.null(groups0)) { # case with an internal factor
    Dd <- length(cnlr)
    lri <- matrix(unlist(sapply(cn, function(j, x) which(apply(x == j, 1, any)),
                                x = do.call(rbind, cmbp), simplify = FALSE)), D, D-1, byrow = TRUE
    ) # i-th row ... indices of LR where the i-th part appears
    pgNAp <- lapply(cmbg, function(x) as.integer(unlist(strsplit(x[1], split = ""))) + as.integer(unlist(strsplit(x[2], split = "")))
    ) # i-th item ... nonzero number indicates unavailability of the respective part when considering the i-th pair of groups
    pgNAlri <- lapply(pgNAp, function(x) sort(unique(c(lri[which(x>0),])))
    ) # i-th item ... indices of the unavailable LR when considering the i-th pair of groups
    gNAlri <- lapply(1:G, function(x) sort(unique(c(lri[which(zrg[x,]>0),])))
    ) # i-th item ... indices of the unavailable LR in th i-th group
    gNRlri <- lapply(seq_along(gNAlri), function(i) {
      oth <- gNAlri[-i]
      com <- Reduce(intersect, oth)
      mis <- setdiff(com, gNAlri [[i]])
      return(c(gNAlri[[i]], mis))
    }) # i-th item ... indices of the non-required LR in th i-th group
    empt <- which(sapply(gNRlri, length) == 0)
    if (length(empt)!=0) gNRlri[[empt]] <- Dd+1
    gRlri <- lapply(gNRlri, function(x) c(1:Dd)[-x]
    ) # i-th item ... indices of the required LR in the i-th group
    disElw <- disElAll(LR, groups, alpha = alpha, rwts = TRUE, pgNAlri = pgNAlri)
    disEl_Prm <- lapply(1:R, function(x) disElAllPerm(LR, groups, alpha = alpha,
                                                      pgNAlri = pgNAlri, gRlri = gRlri, minAR = minAR))
    parts0 <- lapply(lev, function(x) {
      y <- cn[as.numeric(strsplit(x, "")[[1]])>0]
      if(length(y)==0) y <- NULL
      return(y)
    })
    names(parts0) <- lev
  } else { # case with an external factor
    disElw <- disElAll(LR, groups, alpha = alpha, rwts = TRUE)
    disEl_Prm <- lapply(1:R, function(x) disElAllPerm(LR, groups, alpha = alpha))
  }
  disEl <- disElw$disEA
  wts <- disElw$wts
  disOv <- sum(disEl^p, na.rm = TRUE)
  disOv_Prm <- sapply(disEl_Prm, function(x) sum(x^p, na.rm = TRUE))
  pvalOv <- mean(disOv_Prm >= disOv)
  if (is.null(groups0)) {
    res = list(disOv = disOv,
               pvalOv = pvalOv, 
               wts = wts,
               disEl = disEl, 
               disEl_Prm = disEl_Prm,
               parts0 = parts0)
  } else {
    res = list(disOv = disOv,
               pvalOv = pvalOv, 
               wts = wts,
               disEl = disEl, 
               disEl_Prm = disEl_Prm)
  }
  class(res) = "printres"
  
  # --- Optional Post-Hoc Analyses ---
  if (res$pvalOv < alpha) {
    if (posthoc.g) {
      res$posthoc.groups <- perlogPHg(res, p = p, alpha = alpha, mAdj = mAdj)
    }
    if (posthoc.lr) {
      res$posthoc.logratios <- perlogPHlr(res, p = p, alpha = alpha, mAdj = mAdj)
    }
  }
  else {print(res)
        stop(paste("No overall difference between groups concluded at significance level alpha =", alpha))}
  
  print(res)
}

