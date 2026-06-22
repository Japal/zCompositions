lrSVD <- function(X, label = NULL, dl = NULL, frac = 0.65, ncp = 2,
                  imp.missing = FALSE, beta = 0.5, method = c("ridge", "EM"),
                  row.w = NULL, coeff.ridge = 1, threshold = 1e-4, seed = NULL, nb.init = 1,
                  max.iter = 1000, z.warning = 0.8, z.delete = TRUE, ...) {
  
  if (any(X < 0, na.rm = TRUE)) stop("X contains negative values")
  
  if (is.null(dim(X)) || nrow(X) < 2 || ncol(X) < 2) {
    stop("X must be a data matrix with at least two rows and two columns")
  }
  
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)) {
    if (!any(X == label, na.rm = TRUE))
      stop(paste("Label", label, "was not found in the data set"))
    if (label != 0 & any(X == 0, na.rm = TRUE))
      stop("Zero values not labelled as censored values were found in the data set")
    if (any(is.na(X)))
      stop(paste("NA values not labelled as censored values were found in the data set"))
  }
  if (is.na(label)) {
    if (any(X == 0, na.rm = TRUE))
      stop("Zero values not labelled as censored values were found in the data set")
    if (!any(is.na(X), na.rm = TRUE))
      stop(paste("Label", label, "was not found in the data set"))
  }
  
  if (!is.null(row.w)) {
    if (!is.numeric(row.w) || length(row.w) != nrow(X) ||
        any(!is.finite(row.w)) || any(row.w < 0) || sum(row.w) <= 0) {
      stop("row.w must be a numeric vector of finite non-negative weights with length equal to nrow(X)")
    }
  }
  
  if (imp.missing == FALSE){
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.null(dl)){ 
      
      if (is.na(label)) {
        dl <- apply(X, 2, function(x) min(x[!is.na(x)], na.rm = TRUE))
      } else {
        dl <- apply(X, 2, function(x) min(x[x != label], na.rm = TRUE))
      }
      
      warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
    }
    if (is.vector(dl)) dl <- matrix(dl, nrow = 1)
    dl <- as.matrix(dl)
    
    if (!is.numeric(dl)) {
      stop("dl must be numeric")
    }
    
    if (any(dl < 0, na.rm = TRUE) || any(!is.finite(dl), na.rm = TRUE)) {
      stop("Detection limits in dl must be strictly positive and finite")
    }
    
  }
  
  if (imp.missing == FALSE) {
    if (ncol(dl) != ncol(X)) stop("The number of columns in X and dl do not agree")
    
    if (nrow(dl) > 1 && nrow(dl) != nrow(X)) {
      stop("The number of rows in X and dl do not agree")
    }
    
  }

  if (!is.numeric(frac) || length(frac) != 1 || is.na(frac) || frac <= 0 || frac > 1) {
    stop("frac must be a single numeric value in the interval (0, 1]")
  }
  
  if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta < 0 || beta > 1) {
    stop("beta must be a single numeric value between 0 and 1")
  }
  
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold <= 0) {
    stop("threshold must be a positive numeric value")
  }
  
  if (!is.numeric(max.iter) || length(max.iter) != 1 || is.na(max.iter) ||
      max.iter < 1 || max.iter != as.integer(max.iter)) {
    stop("max.iter must be a positive integer")
  }
  
  if (!is.numeric(nb.init) || length(nb.init) != 1 || is.na(nb.init) ||
      nb.init < 1 || nb.init != as.integer(nb.init)) {
    stop("nb.init must be a positive integer")
  }
  
  if (!is.numeric(coeff.ridge) || length(coeff.ridge) != 1 ||
      is.na(coeff.ridge) || coeff.ridge <= 0) {
    stop("coeff.ridge must be a positive numeric value")
  }
  
  if (!is.numeric(z.warning) || length(z.warning) != 1 || is.na(z.warning) ||
      z.warning < 0 || z.warning > 1) {
    stop("z.warning must be a single numeric value between 0 and 1")
  }
  
  method <- match.arg(method)
  
  svd.triplet <- function(X, row.w = NULL, col.w = NULL, ncp = Inf) 
  {
    tryCatch.W.E <- function(expr) {
      W <- NULL
      w.handler <- function(w) {
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), 
                                       warning = w.handler), warning = W)
    }
    if (is.null(row.w)) row.w <- rep(1/nrow(X), nrow(X))
    if (is.null(col.w)) col.w <- rep(1, ncol(X))
    ncp <- min(ncp, nrow(X) - 1, ncol(X))
    row.w <- row.w/sum(row.w)
    X <- t(t(X) * sqrt(col.w)) * sqrt(row.w)
    if (ncol(X) < nrow(X)) {
      svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$value
      if (inherits(svd.usuelle, "error")) {
        svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$value
        if (!inherits(svd.usuelle, "error")) {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        } else {
          bb <- eigen(crossprod(X, X), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d <- bb$values
          svd.usuelle$d[svd.usuelle$d < 0] <- 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vectors[, 1:ncp, drop = FALSE]
          denom <- pmax(svd.usuelle$d[1:ncp], .Machine$double.eps)
          svd.usuelle$u <- t(t(crossprod(t(X), svd.usuelle$v)) / denom)
        }
      }
      U <- svd.usuelle$u
      V <- svd.usuelle$v
      if (ncp > 1) {
        mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
        mult[mult == 0] <- 1
        U <- t(t(U) * mult)
        V <- t(t(V) * mult)
      }
      U <- U/sqrt(row.w)
      V <- V/sqrt(col.w)
    } else {
      svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$value
      if (inherits(svd.usuelle, "error")) {
        svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$value
        if (!inherits(svd.usuelle, "error")) {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        } else {
          bb <- eigen(crossprod(t(X), t(X)), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d <- bb$values
          svd.usuelle$d[svd.usuelle$d < 0] <- 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vectors[, 1:ncp, drop = FALSE]
          denom <- pmax(svd.usuelle$d[1:ncp], .Machine$double.eps)
          svd.usuelle$u <- t(t(crossprod(X, svd.usuelle$v)) / denom)
        }
      }
      
      U <- svd.usuelle$v
      V <- svd.usuelle$u
      mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
      mult[mult == 0] <- 1
      V <- t(t(V) * mult)/sqrt(col.w)
      U <- t(t(U) * mult)/sqrt(row.w)
    }
    vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
    num <- which(vs[1:ncp] < 1e-15)
    if (length(num) == 1) {
      U[, num] <- U[, num, drop = FALSE] * vs[num]
      V[, num] <- V[, num, drop = FALSE] * vs[num]
    }
    if (length(num) > 1) {
      U[, num] <- t(t(U[, num]) * vs[num])
      V[, num] <- t(t(V[, num]) * vs[num])
    }
    res <- list(vs = vs, U = U, V = V)
    return(res)
  }
  
  impute <- function(X = NULL, dl = NULL, bal = NULL, frac = 0.65, ncp = 2, beta = 0.5, method=c("ridge","EM"),
                     row.w = NULL, coeff.ridge = 1, threshold = 1e-4, seed = NULL, max.iter = 1000, init = 1, ...) {
    
    moy.p <- function(V, poids) {
      res <- sum(V * poids, na.rm = TRUE) / sum(poids[!is.na(V)])
    }
    
    gm <- function(x, na.rm = TRUE) {
      n_valid <- if(na.rm) sum(!is.na(x)) else length(x)
      exp(sum(log(x), na.rm = na.rm) / n_valid)
    }
    
    nb.iter <- 1
    old <- Inf
    objective <- 0
    if (!is.null(seed)) {set.seed(seed)} 
    
    missRaw <- which(is.na(X))
    obsRaw <- which(!is.na(X))
    X <- as.matrix(X)
    Xaux <- X 
    
    caux <- rowSums(Xaux, na.rm = TRUE)
    
    if (imp.missing == FALSE) {
      if (init == 1) {X[missRaw] <- frac * dl[missRaw]} 
      else {X[missRaw] <- runif(1, 0.50, 0.8) * dl[missRaw]}
    } else {
      gmeans <- apply(Xaux, 2, function(x) gm(x[!is.na(x) & x > 0], na.rm = TRUE))
      
      if (any(!is.finite(gmeans))) {
        stop("Cannot compute finite geometric means; at least one column has no positive observed values")
      }
      
      nn <- nrow(X)
      
      gmeans <- matrix(gmeans, nrow = nn, ncol = ncol(X), byrow = TRUE)
      X[missRaw] <- gmeans[missRaw]
    }
    
    
    if (any(X <= 0, na.rm = TRUE)) {
      stop("X contains non-positive values before log-ratio transformation")
    }
    
    
    Xhat <- t(bal %*% t(log(X)))
    ncp <- min(ncp, ncol(Xhat), nrow(Xhat) - 1)
    mean.p.or = mean.p <- apply(Xhat, 2, moy.p, row.w)
    Xhat <- t(t(Xhat) - mean.p)
    X <- exp(t(t(bal) %*% t(Xhat)))
    
    X <- X / rowSums(X) 
    fittedX <- fittedXus <- Xhat
    fittedXRaw <- fittedXusRaw <- X
    
    if (ncp == 0) {nb.iter <- 0}
    
    while (nb.iter > 0) {
      X[missRaw] <- fittedXRaw[missRaw]
      Xhat <- t(bal %*% t(log(X)))
      Xhat <- t(t(Xhat) + mean.p)
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      X <- exp(t(t(bal) %*% t(Xhat)))
      
      X <- X / rowSums(X)
      Xaux2 <- X * caux
      viol <- which(Xaux2 > dl)
      Xaux2[viol] <- dl[viol]
      
      
      if (any(Xaux2 <= 0, na.rm = TRUE)) {
        stop("Xaux2 contains non-positive values before log-ratio transformation")
      }
      
      
      Xhat <- t(bal %*% t(log(Xaux2)))
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      X <- exp(t(t(bal) %*% t(Xhat)))
      
      X <- X / rowSums(X)
      
      fittedXusRC <- t(t(fittedXus) + mean.p)
      fittedXusRCRaw <- exp(t(t(bal) %*% t(fittedXusRC)))
      
      fittedXusRCRaw <- fittedXusRCRaw / rowSums(fittedXusRCRaw)
      
      X[obsRaw] <- ((fittedXusRCRaw[obsRaw]) ^ (1 - beta)) * ((X[obsRaw]) ^ beta)
      
      X <- X / rowSums(X)
      Xaux2 <- X * caux
      viol <- which(Xaux2 > dl)
      Xaux2[viol] <- dl[viol]
      
      
      if (any(Xaux2 <= 0, na.rm = TRUE)) {
        stop("Xaux2 contains non-positive values before log-ratio transformation")
      }
      
      Xhat <- t(bal %*% t(log(Xaux2)))
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      Xhat <- t(t(Xhat) - mean.p)
      X <- exp(t(t(bal) %*% t(Xhat)))

      X <- X / rowSums(X)
      svd.res <- svd.triplet(Xhat, row.w = row.w, ncp = ncp)
      
      sigma2 <- nrow(Xhat) * ncol(Xhat) / min(ncol(Xhat), nrow(Xhat) - 1) * sum((svd.res$vs[-c(1:ncp)] ^ 2) / 
                                                                                  ((nrow(Xhat) - 1) * ncol(Xhat) - (nrow(Xhat) - 1) * ncp - ncol(Xhat) * ncp + ncp ^ 2))
      
      
      if ((ncp + 1) <= length(svd.res$vs)) {
        sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1] ^ 2)
      } else {
        sigma2 <- sigma2 * coeff.ridge
      }
      
      if (method == "EM") sigma2 <- 0
      
      lambda.us <- svd.res$vs[1:ncp]
      fittedXus <- tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * row.w) * lambda.us), svd.res$V[, 1:ncp, drop = FALSE])
      fittedXus <- fittedXus / row.w
      
      lambda.shrinked <- (svd.res$vs[1:ncp] ^ 2 - sigma2) / svd.res$vs[1:ncp]
      fittedX <- tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * row.w) * lambda.shrinked), svd.res$V[, 1:ncp, drop = FALSE])
      fittedX <- fittedX / row.w
      
      fittedXRaw <- exp(t(t(bal) %*% t(fittedX)))
      
      fittedXRaw <- fittedXRaw / rowSums(fittedXRaw)
      
      diffRaw <- X / fittedXRaw
      diffRaw[missRaw] <- 1
      
      if (any(diffRaw <= 0, na.rm = TRUE)) {
        stop("diffRaw contains non-positive values before log-ratio transformation")
      }
      
      diff <- t(bal %*% t(log(diffRaw)))
      objective <- sum(diff^2 * row.w)
      criterion <- abs(1 - objective / old)
      old <- objective
      nb.iter <- nb.iter + 1
      
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
        if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
      }
      if (nb.iter > max.iter) {
        nb.iter <- 0
        warning(paste("Stopped after ", max.iter, " iterations"))
      }
    }
    
    Xhat <- t(t(Xhat) + mean.p.or)
    X <- exp(t(t(bal) %*% t(Xhat)))
    
    X <- X / rowSums(X)
    
    
    if (any(rowSums(Xaux, na.rm = TRUE) <= 0)) {
      stop("At least one row has zero observed sum; log-ratio imputation cannot proceed")
    }
    
    
    completeObs <- Xaux / rowSums(Xaux, na.rm = TRUE)
    completeObs[missRaw] <- X[missRaw]
    completeObs <- completeObs / rowSums(completeObs)
    completeObs <- completeObs * caux
    viol <- which(completeObs > dl)
    completeObs[viol] <- dl[viol]
    completeObs[obsRaw] <- Xaux[obsRaw]
    
    completeObs <- completeObs / rowSums(completeObs)
    fittedX <- t(t(fittedX) + mean.p)
    fittedXRaw <- exp(t(t(bal) %*% t(fittedX)))
    
    fittedXRaw <- fittedXRaw / rowSums(fittedXRaw)
    
    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedXRaw
    return(result) 
  }
  
  X <- as.data.frame(X)
  nn <- nrow(X)
  D <- ncol(X)
  
  if (!is.na(label)) {
    X[X == label] <- NA
  }
  
  X <- as.data.frame(apply(X, 2, as.numeric))
  
  c <- rowSums(X, na.rm = TRUE)
  
  checkNumZerosCol <- colSums(is.na(X))
  
  if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
    cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (ncol(X)-2)) {
        stop(paste("Almost all columns contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).", sep=""))
      }      
      
      X <- X[, -cases, drop = FALSE]
      
      if (imp.missing == FALSE) {
        dl <- dl[, -cases, drop = FALSE]
      }
      
      action <- "deleted"
      warning(paste("Column no. ", cases, " containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n", sep=""))
    } else {      
      action <- "found"      
      warning(paste("Column no. ", cases, " containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete. Check out with zPatterns()).\n", sep=""))      
    }
  }
  
  checkNumZerosRow <- rowSums(is.na(X))  
  if (any(checkNumZerosRow/ncol(X) > z.warning)) {    
    cases <- which(checkNumZerosRow/ncol(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (nrow(X)-2)) {
        stop(paste("Almost all rows contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).", sep=""))
      }
      
      X <- X[-cases, , drop = FALSE]
      
      if (imp.missing == FALSE && nrow(dl) > 1) {
        dl <- dl[-cases, , drop = FALSE]
      }
      
      if (!is.null(row.w)) {
        row.w <- row.w[-cases]
      }
      
      action <- "deleted"
      warning(paste("Row no. ", cases, " containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n", sep=""))
    } else {      
      action <- "found"      
      warning(paste("Row no. ", cases, " containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete. Check out with zPatterns()).\n", sep=""))      
    }
  }
  
  # Update nn, D, and c after potential row/column deletions
  nn <- nrow(X)
  D <- ncol(X)
  c <- rowSums(X, na.rm = TRUE)
  
  
  if (!is.numeric(ncp) || length(ncp) != 1 || is.na(ncp) ||
      ncp < 1 || ncp != as.integer(ncp)) {
    stop("ncp must be a positive integer")
  }
  
  if (ncp > min(nn - 2, D - 2)) {
    stop("ncp is too large for the size of the data matrix")
  }
  
  
  if (is.null(row.w)) {
    row.w <- rep(1, nn) / nn
  } else {
    
    if (!is.numeric(row.w) || length(row.w) != nn ||
        any(!is.finite(row.w)) || any(row.w < 0) || sum(row.w) <= 0) {
      stop("row.w must be a numeric vector of finite non-negative weights with length equal to nrow(X)")
    }
    
    row.w <- row.w / sum(row.w)
  }
  
  
  
  closed <- 0
  if (all(abs(c - mean(c)) < .Machine$double.eps^0.3)) closed <- 1
  
  Xaux <- as.matrix(X)
  
  
  if (any(rowSums(Xaux, na.rm = TRUE) <= 0)) {
    stop("At least one row has zero observed sum; log-ratio imputation cannot proceed")
  }
  
  
  XauxClosed <- Xaux / rowSums(Xaux, na.rm = TRUE) 
  Xna <- X 
  
  pz <- colSums(is.na(Xna)) / nn
  X <- Xna[, order(-pz), drop = FALSE]
  
  if (imp.missing == FALSE) {
    
    if (nrow(dl) == 1) {
      dl <- matrix(dl[, order(-pz), drop = FALSE], nrow = 1)
    } else {
      dl <- dl[, order(-pz), drop = FALSE]
    }
  }
  
  Smat <- diag(rep(1,D))
  Smat[upper.tri(Smat)] <- -1
  Smat <- Smat[-D,]
  bal <- Smat
  numsbp <- dim(Smat)[1]
  for (f in 1:numsbp) {
    den <- sum(bal[f, ] == -1)
    num <- sum(bal[f, ] == 1)
    bal[f, bal[f, ] == 1] <- sqrt(den / ((den + num) * num))
    bal[f, bal[f, ] == -1] <- -sqrt(num / ((den + num) * den))
  }
  
  if (imp.missing == FALSE) {
    if (nrow(dl) == 1) {
      dl <- matrix(dl, nrow = nn, ncol = ncol(X), byrow = TRUE)
    }
  } else {
    dl <- matrix(0, nrow = nn, ncol = D)
  } 
  
  observedRaw <- which(!is.na(X))
  missingRaw <- which(is.na(X))
  Xmax <- matrix(apply(X, 2, max, na.rm = TRUE), nrow = nn, ncol = ncol(X), byrow = TRUE)
  
  if (any(!is.finite(Xmax))) {
    stop("At least one column has no finite observed values")
  }
  
  dl[observedRaw] <- Xmax[observedRaw]
  if (imp.missing == TRUE) {
    dl[missingRaw] <- Xmax[missingRaw]
  }
  colnames(dl) <- colnames(X)
  
  for (i in 1:nb.init) {
    
    if (!any(is.na(X))) {
      X <- X[, colnames(Xna), drop = FALSE]
      return(as.data.frame(X))
    }
    
    res.impute <- impute(X = X, dl = dl, bal = bal, frac = frac, ncp = ncp, beta = beta, method = method,
                         row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold, 
                         seed = if (!is.null(seed)) {
                           seed + i - 1
                         } else {
                           NULL
                         }
                         , max.iter = max.iter, init = i)
  }
  
  Y <- res.impute$completeObs
  XauxClosed <- XauxClosed[, order(-pz), drop = FALSE]
  XauxClosed[missingRaw] <- Y[missingRaw]
  X <- XauxClosed * c
  
  if (closed == 1) {
    X <- (X / rowSums(X)) * c[1]
  }
  
  X <- X[, colnames(Xna), drop = FALSE]
  return(as.data.frame(X))
}