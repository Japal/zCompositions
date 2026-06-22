#' Test of differences in group means for compositional data
#'
#' @description A nonparametric permutation test to assess the hypothesis of equality of means between subsets of observations according to an externally or internally defined factor variable. If any, zero patterns are considered as default internal grouping factor.
#'
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param groups Factor variable indicating the grouping structure. If NULL (default), any zero patterns in the data will be used as internal grouping factor.
#' 
#' Note that if a grouping factor is set by the user, then any zeros in the data must be previously dealt with, e.g. by imputation.
#' @param p Power parameter selected in overall dissimilarity test statistic, either automatically (default = "auto") or manually fixed.
#' @param alpha Significance level parameter (default = 0.05).
#' @param R Number of permutation resamples (default = 1000).
#' @param posthoc.g Logical. If TRUE, performs post-hoc analysis for pairs of groups (default = FALSE).
#' @param posthoc.lr Logical. If TRUE, performs post-hoc analysis for logratios (default = FALSE).
#' @param mAdj Adjustment of p-values for multiple post-hoc testing (see \code{\link{p.adjust}}). Default is Benjamini and Hochberg's FDR method (default = "BH").
#'
#' @return A list object of class "perLog.output" containing summaries of results:
#' \item{disOv}{Overall dissimilarity test statistic.}
#' \item{pvalOv}{Overall permutation p-value.}
#' \item{p}{Power parameter used in overall dissimilarity test statistic.}
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
#' 
#' The power parameter p can be either automatically selected or manually fixed. For automatic selection, a simple conservative strategy is implemented starting with p = 10, as a liberal reference, and then successively setting p = \{2, 3, ..., 9\} until less significant differences are no longer obtained at the overall and group comparison levels.
#' @references Štefelová N, Palarea-Albaladejo J, Martín-Fernández JA. A permutation test of differences between externally or internally defined groupings in compositional data sets. Statistical Methods in Medical Research. 2026;0(0). <doi:10.1177/09622802251413737>
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

perLog <- function(X, groups = NULL, p = "auto", alpha = 0.05, R = 1000, 
                    posthoc.g = FALSE, posthoc.lr = FALSE, mAdj = "BH") {
  
  # Auxiliary functions
  
  fast_stats <- function(mat) {
    valid <- !is.na(mat)
    n <- colSums(valid)
    
    m <- rep(NA_real_, ncol(mat))
    v <- rep(NA_real_, ncol(mat))
    
    ok_mean <- n > 0
    m[ok_mean] <- colSums(mat[, ok_mean, drop = FALSE], na.rm = TRUE) / n[ok_mean]
    
    ok_var <- n > 1
    if (any(ok_var)) {
      diffs <- t(t(mat[, ok_var, drop = FALSE]) - m[ok_var])
      v[ok_var] <- colSums(diffs^2, na.rm = TRUE) / (n[ok_var] - 1)
    }
    
    list(n = n, m = m, v = v)
  }
  
  disElPair <- function(LR1, LR2, alpha = 0.05, rwts = FALSE, pgNAlri = NULL) {
    Dd <- ncol(LR1)
    
    if (length(pgNAlri) != 0) {
      LR1 <- LR1[, -pgNAlri, drop = FALSE]
      LR2 <- LR2[, -pgNAlri, drop = FALSE]
    }
    
    if (ncol(LR1) == 0) {
      disEP <- wts <- rep(NA_real_, Dd)
    } else {
      s1 <- fast_stats(LR1)
      s2 <- fast_stats(LR2)
      
      nmr <- s1$m - s2$m
      var1_n1 <- s1$v / s1$n
      var2_n2 <- s2$v / s2$n
      dnm <- sqrt(var1_n1 + var2_n2)
      
      wts_calc <- nmr / dnm
      
      df_calc <- dnm^4 / (
        (1 / (s1$n - 1)) * var1_n1^2 +
          (1 / (s2$n - 1)) * var2_n2^2
      )
      
      bad <- !is.finite(wts_calc) |
        !is.finite(df_calc) |
        df_calc <= 0 |
        !is.finite(dnm) |
        dnm <= 0
      
      wts_calc[bad] <- NA_real_
      df_calc[bad] <- NA_real_
      
      awts <- abs(wts_calc)
      qN <- qnorm(pt(awts, df_calc))
      qN[is.infinite(qN)] <- qnorm(0.9999999999999999)
      qN[!is.finite(qN)] <- NA_real_
      
      qNs <- qnorm(1 - alpha / 2)
      disEP_calc <- qN / qNs
      
      if (length(pgNAlri) != 0) {
        disEP <- rep(NA_real_, Dd)
        disEP[-pgNAlri] <- disEP_calc
        if (rwts) {
          wts <- rep(NA_real_, Dd)
          wts[-pgNAlri] <- wts_calc
        }
      } else {
        disEP <- disEP_calc
        if (rwts) wts <- wts_calc
      }
    }
    
    if (rwts) return(list(disEP = disEP, wts = wts))
    return(disEP)
  }
  
  # Elemental dissimilarities for all pairs
  disElAll <- function(LR, groups, alpha = 0.05, rwts = FALSE, pgNAlri = NULL) {
    Dd <- ncol(LR)
    cnlr <- colnames(LR)
    cmbg <- combn(levels(groups), 2, simplify = FALSE)
    Gg <- length(cmbg)
    grcmb <- sapply(cmbg, function(x) paste(x, collapse = " vs "))
    
    if (length(pgNAlri) != 0) {
      res_list <- lapply(1:Gg, function(x) {
        disElPair(LR[groups == cmbg[[x]][1], , drop = FALSE], 
                  LR[groups == cmbg[[x]][2], , drop = FALSE],
                  alpha = alpha, rwts = rwts, pgNAlri = pgNAlri[[x]])
      })
    } else {
      res_list <- lapply(cmbg, function(x) {
        disElPair(LR[groups == x[1], , drop = FALSE], 
                  LR[groups == x[2], , drop = FALSE],
                  alpha = alpha, rwts = rwts)
      })
    }
    
    if (rwts) {
      disEA <- matrix(unlist(lapply(res_list, `[[`, 1)), Gg, Dd, byrow = TRUE)
      wts <- matrix(unlist(lapply(res_list, `[[`, 2)), Gg, Dd, byrow = TRUE)
      rownames(disEA) <- rownames(wts) <- grcmb
      colnames(disEA) <- colnames(wts) <- cnlr
      return(list(disEA = disEA, wts = wts))
    } else {
      disEA <- matrix(unlist(res_list), Gg, Dd, byrow = TRUE)
      rownames(disEA) <- grcmb
      colnames(disEA) <- cnlr
      return(disEA)
    }
  }
  
  # Permutation function
  disElAllPerm <- function(LR, groups, alpha = 0.05, pgNAlri = NULL, gRlri = NULL, minAR = 5) {
    lev <- levels(groups)
    n <- nrow(LR)
    indR <- sample(n)
    LRR <- LR[indR, , drop = FALSE]
    
    if (length(gRlri) != 0) {
      attempts <- 0
      repeat {
        minA <- min(vapply(seq_along(lev), function(x) {
          idx <- gRlri[[x]]
          
          if (length(idx) == 0) {
            return(Inf)
          }
          
          min(colSums(!is.na(LRR[groups == lev[x], idx, drop = FALSE])))
        }, numeric(1)))
        
        if (minA >= minAR || attempts > 100) break
        indR <- sample(n)
        LRR <- LR[indR, , drop = FALSE]
        attempts <- attempts + 1
      }
      return(disElAll(LRR, groups, alpha = alpha, pgNAlri = pgNAlri))
    }
    return(disElAll(LRR, groups, alpha = alpha))
  }
  
  # Custom print helper
  print_summary <- function(x) {
    cat("\nOverall dissimilarity:\n\nTest statistic E = ", round(x$disOv, 4), 
        ", p-value = ", round(x$pvalOv, 4), "\n\n", sep = "")
    
    if (!is.null(x$posthoc.groups)) {
      cat("Post-hoc analysis by pairs of groups:\n\n")
      print(x$posthoc.groups$summary)
      if (!is.null(x$posthoc.groups$parts0)) {
        cat("\nZero patterns and zero parts therein:\n\n")
        print(x$posthoc.groups$parts0)
      }
    }
    if (!is.null(x$posthoc.logratios)) {
      cat("\nPost-hoc analysis for pairwise logratios:\n")
      for (i in seq_along(x$posthoc.logratios$summary)) {
        cat("\nPair:", names(x$posthoc.logratios$summary)[i], "\n")
        print(round(x$posthoc.logratios$summary[[i]], 4))
      }
    }
  }
  
  # Post-hoc groups
  perlogPHg <- function(perlogR, alpha = 0.05, mAdj = "BH") {
    p <- perlogR$p
    disEl <- perlogR$disEl
    Gg <- nrow(disEl)
    
    if (Gg == 1) {
      summary_mat <- c("Bet-Grp diss" = perlogR$disOv, "%Ctb overall diss" = 100, "Adj p-value" = perlogR$pvalOv)
    } else {
      disBg <- rowSums(disEl^p, na.rm = TRUE)
      rcBgOv <- 100 * disBg / sum(disBg, na.rm = TRUE)
      disBg_Prm <- t(sapply(perlogR$disEl_Prm, function(x) rowSums(x^p, na.rm = TRUE)))
      disBg_All <- rbind(disBg, disBg_Prm)
      pvalBg <- apply(disBg_All, 2, function(x) mean(x[-1] >= x[1]))
      pvalBgAdj <- p.adjust(pvalBg, method = mAdj)
      summary_mat <- cbind("Bet-Grp diss" = round(disBg, 4), "%Ctb overall diss" = round(rcBgOv, 4), "Adj p-value" = round(pvalBgAdj, 4))
    }
    
    res <- list(summary = summary_mat)
    if (!is.null(perlogR$parts0)) res$parts0 <- perlogR$parts0
    return(res)
  }
  
  # Post-hoc pairwise
  perlogPHlr <- function(perlogR, alpha = 0.05, mAdj = "BH") {
    p <- perlogR$p
    disEl <- perlogR$disEl
    pr <- rownames(disEl)
    npr <- nrow(disEl)
    
    if (npr == 1) {
      wts <- perlogR$wts[1, ]
      disEl_row <- disEl[1, ]
      disEl.p <- disEl_row^p
      rcElBg <- rcElOv <- 100 * disEl.p / perlogR$disOv
      disEl.p_Prm <- t(sapply(perlogR$disEl_Prm, function(x) x^p))
      disEl.p_All <- rbind(disEl.p, disEl.p_Prm)
      pvalEl <- apply(disEl.p_All, 2, function(x) mean(x[-1] >= x[1]))
      pvalElAdj <- p.adjust(pvalEl, method = mAdj)
      
      summary_list <- list(cbind("wts" = wts, "Elem diss" = disEl_row, "%Ctb bet-grp diss" = rcElBg, 
                                 "%Ctb overall diss" = rcElOv, "Adj p-value" = pvalElAdj))
      names(summary_list) <- pr
    } else {
      perlogPHgR <- perlogPHg(perlogR, alpha = alpha, mAdj = mAdj)
      iS <- which(perlogPHgR$summary[, 3] < alpha)
      npr <- length(iS)
      if (npr == 0) stop("No significantly different pair.")
      
      pr <- pr[iS]
      wts <- perlogR$wts[iS, , drop = FALSE]
      disEl <- disEl[iS, , drop = FALSE]
      disEl.p <- disEl^p
      rcElBg <- t(100 * t(disEl.p) / perlogPHgR$summary[iS, 1])
      rcElOv <- 100 * disEl.p / perlogR$disOv
      disEl_Prm <- lapply(perlogR$disEl_Prm, function(x) x[iS, , drop = FALSE])
      
      if (npr == 1) {
        disEl.p_Prm <- t(sapply(disEl_Prm, function(x) x^p))
        disEl.p_All <- rbind(disEl.p, disEl.p_Prm)
        pvalEl <- apply(disEl.p_All, 2, function(x) mean(x[-1] >= x[1]))
        pvalElAdj <- p.adjust(pvalEl, method = mAdj)
        summary_list <- list(cbind("wts" = wts[1,], "Elem diss" = disEl[1,], "%Ctb bet-grp diss" = rcElBg[1,], 
                                   "%Ctb overall diss" = rcElOv[1,], "Adj p-value" = pvalElAdj))
        names(summary_list) <- pr
      } else {
        disEl.p_Prm <- lapply(disEl_Prm, function(x) x^p)
        disEl.p_All <- c(list(disEl.p), disEl.p_Prm)
        disEl.p_All <- lapply(1:npr, function(i) t(sapply(disEl.p_All, function(x) x[i, ])))
        pvalEl <- t(sapply(disEl.p_All, function(y) apply(y, 2, function(x) mean(x[-1] >= x[1]))))
        rownames(pvalEl) <- pr
        pvalElAdj <- pvalEl
        pvalElAdj[] <- p.adjust(c(pvalEl), method = mAdj)
        summary_list <- lapply(1:npr, function(i) cbind("wts" = wts[i, ], "Elem diss" = disEl[i, ], 
                                                        "%Ctb bet-grp diss" = rcElBg[i, ], "%Ctb overall diss" = rcElOv[i, ], 
                                                        "Adj p-value" = pvalElAdj[i, ]))
        names(summary_list) <- pr
      }
    }
    
    significance <- lapply(summary_list, function(x) {
      cnlr <- rownames(x)
      Dd <- length(cnlr)
      D <- (1 + sqrt(1 + 8 * Dd)) / 2
      cn <- unique(unlist(strsplit(cnlr, "/")))
      mt <- matrix(NA_real_, D, D)
      sgn <- rep(0, Dd)
      sgn[which(x[, 5] < alpha)] <- 1
      sgn <- sign(x[, 1]) * sgn
      mt[lower.tri(mt)] <- -sgn
      mtt <- -t(mt)
      mt[upper.tri(mt)] <- mtt[upper.tri(mtt)]
      rownames(mt) <- colnames(mt) <- cn
      return(mt)
    })
    
    res <- list(summary = summary_list, significance = significance)
    if (!is.null(perlogR$parts0)) res$parts0 <- perlogR$parts0
    return(res)
  }
  
  # Main block
  
  
  if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
    stop("X must be a matrix or data.frame with at least two rows and two columns")
  }
  
  
  
  cn0 <- colnames(X)
  
  X <- as.data.frame(X)
  X <- as.data.frame(apply(X, 2, as.numeric))
  
  if (!is.null(cn0)) {
    colnames(X) <- cn0
  } else {
    colnames(X) <- NULL
  }
  
  if (any(X < 0, na.rm = TRUE)) {
    stop("X contains negative values")
  }
  
  if (any(is.na(X))) {
    stop("X contains missing values")
  }
  
  if (any(rowSums(X, na.rm = TRUE) <= 0)) {
    stop("At least one row has zero total")
  }
  
  if (any(X == 0, na.rm = TRUE) && !is.null(groups)) {
    stop("User-defined factor set but zero values found in the data set")
  }
  
  if (!is.null(groups)) {
    if (length(groups) != nrow(X)) {
      stop("Length of groups must be equal to the number of rows in X")
    }
    
    if (any(is.na(groups))) {
      stop("groups contains missing values")
    }
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in the interval (0, 1)")
  }
  
  if (!is.numeric(R) || length(R) != 1 || is.na(R) ||
      R < 1 || R != as.integer(R)) {
    stop("R must be a positive integer")
  }
  
  if (!is.logical(posthoc.g) || length(posthoc.g) != 1 || is.na(posthoc.g)) {
    stop("posthoc.g must be TRUE or FALSE")
  }
  
  if (!is.logical(posthoc.lr) || length(posthoc.lr) != 1 || is.na(posthoc.lr)) {
    stop("posthoc.lr must be TRUE or FALSE")
  }
  
  if (!is.character(mAdj) || length(mAdj) != 1 ||
      !(mAdj %in% p.adjust.methods)) {
    stop("mAdj must be one of p.adjust.methods")
  }
  
  groups0 <- groups
  if (is.null(groups0)) {
    zr <- 1 * (X == 0)
    groups <- apply(zr, 1, paste, collapse = "")
    zrg <- apply(do.call(rbind, strsplit(sort(unique(groups)), split = "")), 2, as.integer)
  }
  
  groups <- as.factor(groups)
  nk <- table(groups)
  
  if (length(nk) < 2) {
    stop("At least two groups are required")
  }
  
  if (is.character(p)) {
    if (length(p) != 1 || p != "auto") {
      stop("p must be a positive numeric value or 'auto'")
    }
  } else {
    if (!is.numeric(p) || length(p) != 1 || is.na(p) || p <= 0) {
      stop("p must be a positive numeric value or 'auto'")
    }
  }
  
  
  minnk <- min(nk)
  minAR <- min(minnk, 5)
  if (minnk < 5) warning("At least 1 group has less than 5 observations")
  
  lev <- levels(groups)
  G <- length(lev)
  cmbg <- combn(lev, 2, simplify = FALSE)
  Gg <- length(cmbg)
  
  cn <- colnames(X)
  D <- ncol(X)
  if (is.null(cn)) {
    cn <- paste0("P", 1:D)
    colnames(X) <- cn
  }
  
  cmbp <- combn(cn, 2, simplify = FALSE)
  cnlr <- sapply(cmbp, paste, collapse = "/")
  
  X_mat <- as.matrix(X)
  LR <- log(X_mat[, sapply(cmbp, `[`, 1)] / X_mat[, sapply(cmbp, `[`, 2)])
  LR[is.infinite(LR)] <- NA
  LR[is.nan(LR)] <- NA
  colnames(LR) <- cnlr
  
  if (is.null(groups0)) {
    Dd <- length(cnlr)
    lri <- matrix(unlist(lapply(cn, function(j) which(apply(do.call(rbind, cmbp) == j, 1, any)))), D, D - 1, byrow = TRUE)
    pgNAp <- lapply(cmbg, function(x) as.integer(unlist(strsplit(x[1], split = ""))) + as.integer(unlist(strsplit(x[2], split = ""))))
    pgNAlri <- lapply(pgNAp, function(x) sort(unique(c(lri[which(x > 0), ]))))
    gNAlri <- lapply(1:G, function(x) sort(unique(c(lri[which(zrg[x, ] > 0), ]))))
    
    gNRlri <- lapply(seq_along(gNAlri), function(i) {
      oth <- gNAlri[-i]
      com <- Reduce(intersect, oth)
      mis <- setdiff(com, gNAlri[[i]])
      c(gNAlri[[i]], mis)
    })
    
    empt <- which(lengths(gNRlri) == 0)
    if (length(empt) != 0) gNRlri[[empt]] <- Dd + 1
    gRlri <- lapply(gNRlri, function(x) c(1:Dd)[-x])
    
    disElw <- disElAll(LR, groups, alpha = alpha, rwts = TRUE, pgNAlri = pgNAlri)
    disEl_Prm <- lapply(1:R, function(x) disElAllPerm(LR, groups, alpha = alpha, pgNAlri = pgNAlri, gRlri = gRlri, minAR = minAR))
    
    parts0 <- lapply(lev, function(x) {
      y <- cn[as.numeric(strsplit(x, "")[[1]]) > 0]
      if (length(y) == 0) NULL else y
    })
    names(parts0) <- lev
  } else {
    disElw <- disElAll(LR, groups, alpha = alpha, rwts = TRUE)
    disEl_Prm <- lapply(1:R, function(x) disElAllPerm(LR, groups, alpha = alpha))
  }
  
  disEl <- disElw$disEA
  wts <- disElw$wts
  
  if (all(is.na(disEl))) {
    stop("No valid pairwise logratio dissimilarities could be computed")
  }
  
  if (p != "auto") {
    disOv <- sum(disEl^p, na.rm = TRUE)
    disOv_Prm <- sapply(disEl_Prm, function(x) sum(x^p, na.rm = TRUE))
    pvalOv <- mean(disOv_Prm >= disOv)
  } else {
    pp <- c(10, 2:9)
    nsig10 <- Inf
    nsigp <- -Inf
    ppi <- 1
    while (nsig10 > nsigp && ppi <= length(pp)) {
      p_val <- pp[ppi]
      disOv <- sum(disEl^p_val, na.rm = TRUE)
      disOv_Prm <- sapply(disEl_Prm, function(x) sum(x^p_val, na.rm = TRUE))
      pvalOv <- mean(disOv_Prm >= disOv)
      
      if (pvalOv > alpha) {
        nsig <- 0
      } else {
        disBg <- rowSums(disEl^p_val, na.rm = TRUE)
        disBg_Prm <- t(sapply(disEl_Prm, function(x) rowSums(x^p_val, na.rm = TRUE)))
        disBg_All <- rbind(disBg, disBg_Prm)
        pvalBg <- apply(disBg_All, 2, function(x) mean(x[-1] >= x[1]))
        nsig <- 1 + sum(p.adjust(pvalBg, method = mAdj) < alpha)
      }
      
      if (ppi == 1) nsig10 <- nsig else nsigp <- nsig
      ppi <- ppi + 1
      p <- p_val
    }
  }  
  
  res <- list(disOv = disOv, pvalOv = pvalOv, p = p, wts = wts, disEl = disEl, disEl_Prm = disEl_Prm)
  if (is.null(groups0)) res$parts0 <- parts0
  
  # Post-hoc & output
  
  if (res$pvalOv < alpha) {
    if (posthoc.g) res$posthoc.groups <- perlogPHg(res, alpha = alpha, mAdj = mAdj)
    if (posthoc.lr) res$posthoc.logratios <- perlogPHlr(res, alpha = alpha, mAdj = mAdj)
  } else {
    print_summary(res)
    stop(paste("No overall difference between groups concluded at significance level alpha =", alpha))
  }
  
  class(res) <- "perLog.output"
  
  print_summary(res)
  return(invisible(res))
}