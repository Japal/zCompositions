#' @title Log-ratio SVD algorithm
#' @description This function implements an algorithm to impute left-censored data (e.g. values below detection limit,
#' rounded zeros) based on the singular value decomposition (SVD) of a compositional data set.
#' @param X Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
#' @param label Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved values in \code{X}.
#' @param dl Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}.
#' @param frac Parameter for initial multiplicative simple replacement of left-censored data (see \code{\link{multRepl}}) (default = 0.65).
#' @param ncp PARAM_DESCRIPTION, Default: 2
#' @param beta PARAM_DESCRIPTION, Default: 0.5
#' @param scale PARAM_DESCRIPTION, Default: FALSE
#' @param method PARAM_DESCRIPTION, Default: 'Regularized'
#' @param row.w PARAM_DESCRIPTION, Default: NULL
#' @param coeff.ridge PARAM_DESCRIPTION, Default: 1
#' @param threshold PARAM_DESCRIPTION, Default: 1e-04
#' @param seed PARAM_DESCRIPTION, Default: NULL
#' @param nb.init PARAM_DESCRIPTION, Default: 1
#' @param maxiter PARAM_DESCRIPTION, Default: 1000
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 

lrSVD <- function(X,label=NULL,dl=NULL,frac=0.65,ncp=2,beta=0.5,scale=FALSE,method="Regularized",row.w=NULL,
                  coeff.ridge=1,threshold=1e-4,seed=NULL,nb.init=1,maxiter=1000,...){
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)
  dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
  }
  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  
  if (ncp>min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large for the size of the data matrix")
  if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X) # Equal weight for all rows
  
  svd.triplet <- function (X, row.w = NULL, col.w = NULL, ncp = Inf) # From FactoMineR package
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
    if (is.null(row.w)) 
      row.w <- rep(1/nrow(X), nrow(X))
    if (is.null(col.w)) 
      col.w <- rep(1, ncol(X))
    ncp <- min(ncp, nrow(X) - 1, ncol(X))
    row.w <- row.w/sum(row.w)
    X <- t(t(X) * sqrt(col.w)) * sqrt(row.w)
    if (ncol(X) < nrow(X)) {
      svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "message") {
        svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
        if (names(svd.usuelle)[[1]] == "d") {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        }
        else {
          bb <- eigen(crossprod(X, X), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d[svd.usuelle$d < 0] = 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vec[, 1:ncp]
          svd.usuelle$u <- t(t(crossprod(t(X), svd.usuelle$v))/svd.usuelle$d[1:ncp])
        }
      }
      U <- svd.usuelle$u
      V <- svd.usuelle$v
      if (ncp > 1) {
        mult <- sign(as.vector(crossprod(rep(1, nrow(V)), 
                                         as.matrix(V))))
        mult[mult == 0] <- 1
        U <- t(t(U) * mult)
        V <- t(t(V) * mult)
      }
      U <- U/sqrt(row.w)
      V <- V/sqrt(col.w)
    }
    else {
      svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "message") {
        svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
        if (names(svd.usuelle)[[1]] == "d") {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        }
        else {
          bb <- eigen(crossprod(t(X), t(X)), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d[svd.usuelle$d < 0] = 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vec[, 1:ncp]
          svd.usuelle$u <- t(t(crossprod(X, svd.usuelle$v))/svd.usuelle$d[1:ncp])
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
  
  impute <- function (X,dl,bal,frac=frac,ncp=4,beta=0.5,scale=FALSE,method=NULL,threshold = 1e-6,
                      seed = NULL,init=1,maxiter=1000,row.w=NULL,coeff.ridge=1,...){
    # X == CoDa DATA
    # dl ==  matrix of DL
    # bal == MATRIX balances (rows) (D-1 x D)
    # OLR-coordinates: olr.Y <- t(bal%*%t(log(Y)))
    # INV.OLR-coord: Y.closed <-exp(t(t(bal)%*%t(olr.Y)))
    #                Y.closed <- Y.closed/apply(Y.closed,1,sum)
    # 
    ################
    #
    ##### moy.p = weighted AVERAGE = moyenne poids
    moy.p <- function(V, poids) {
      res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
    }
    ##### end
    ##### ec = weighted average SQUARE ERROR 
    ec <- function(V, poids) {
      res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
    }
    ##### end 
    
    nb.iter <- 1
    old <- Inf
    objective <- 0
    if (!is.null(seed)){set.seed(seed)} # fix seed to have same results
    #
    # OLR of initial DATA MATRIX
    # missing pattern
    missRaw <- which(is.na(X))
    obsRaw <- which(!is.na(X))
    X <- as.matrix(X)
    Xaux <- X # copy original raw data with NA
    #
    caux <- apply(Xaux, 1, sum, na.rm = TRUE)
    #
    # impute random normal values to the missing if init > 1 (= 1 DEFAULT value)
    if (init==1){
      X[missRaw] <- frac*dl[missRaw]}
    else{
      X[missRaw]<-runif(1,0.50,0.8)*dl[missRaw] #random
    }
    # DL closed
    dlclosed <- dl/apply(dl,1,sum)
    # Xhat: OLR-coordinates
    Xhat <- t(bal%*%t(log(X)))
    #
    # number of components
    ncp <- min(ncp,ncol(Xhat),nrow(Xhat)-1)
    # weighted column mean  
    mean.p <- apply(Xhat, 2, moy.p,row.w)
    # matrix centering
    Xhat <- t(t(Xhat)-mean.p)
    # weighted column variance
    et <- apply(Xhat, 2, ec,row.w)
    # scaling
    if (scale) Xhat <- t(t(Xhat)/et)
    # update X: olr.inv(Xhat)
    X <-exp(t(t(bal)%*%t(Xhat)))
    X <- X/apply(X,1,sum)
    # aux data matrix for observed and non-observed data
    fittedX <- fittedXus <- Xhat
    fittedXRaw <- fittedXusRaw <- X
    if (ncp==0) {nb.iter <- 0}
    
    while (nb.iter > 0) {
      # update data matrix  
      X[missRaw] <- fittedXRaw[missRaw]
      # Xhat: OLR-coordinates
      Xhat <- t(bal%*%t(log(X)))
      # scaling if proceed
      if (scale) {Xhat <- t(t(Xhat)*et)}
      # recover the center
      Xhat <- t(t(Xhat)+mean.p)
      # update X: olr.inv(Xhat)
      X <- exp(t(t(bal)%*%t(Xhat)))
      X <- X/apply(X,1,sum)
      # impute observed values
      fittedXusRC <- t(t(fittedXus)+mean.p)# recover center
      # RAW values
      # INV-OLR # 
      fittedXusRCRaw <- exp(t(t(bal)%*%t(fittedXusRC)))
      fittedXusRCRaw <- fittedXusRCRaw/apply(fittedXusRCRaw,1,sum)
      #
      X[obsRaw] <- ((fittedXusRCRaw[obsRaw])^(1-beta))*((X[obsRaw])^beta)
      # check the detection limit DL
      X <- X/apply(X,1,sum)
      #
      Xaux2 <- X*caux # original units
      viol <- which(Xaux2>dl)
      Xaux2[viol] <- dl[viol]
      # Xhat: OLR-coordinates
      Xhat <- t(bal%*%t(log(Xaux2)))
      # update mean
      mean.p <- apply(Xhat, 2, moy.p,row.w)
      # centering
      Xhat <- t(t(Xhat)-mean.p)
      # scaling, if proceed
      et <- apply(Xhat, 2, ec,row.w)
      if (scale) Xhat <- t(t(Xhat)/et)
      # update X
      # INV-OLR 
      X <- exp(t(t(bal)%*%t(Xhat)))
      X <- X/apply(X,1,sum)
      # SVD calculation WEIGHTED by row.w and rank ncp
      svd.res <- svd.triplet(Xhat,row.w=row.w,ncp=ncp)
      # discarted formula for sigma (appears in 2012 paper)
      #       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
      # adopted formula for sigma (appears in 2016 package paper)
      sigma2  <- nrow(Xhat)*ncol(Xhat)/min(ncol(Xhat),nrow(Xhat)-1)* sum((svd.res$vs[-c(1:ncp)]^2)/((nrow(Xhat)-1) * ncol(Xhat) - (nrow(Xhat)-1) * ncp - ncol(Xhat) * ncp + ncp^2))
      sigma2 <- min(sigma2*coeff.ridge,svd.res$vs[ncp+1]^2)
      if (method=="em") sigma2 <- 0
      # usual lambda 
      lambda.us <- svd.res$vs[1:ncp]
      # calculate the usual new matrix
      fittedXus <- tcrossprod(t(t(svd.res$U[,1:ncp,drop=FALSE]*row.w)*lambda.us),svd.res$V[,1:ncp,drop=FALSE])
      fittedXus <- fittedXus/row.w
      # lambda for regularization
      lambda.shrinked <- (svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
      # calculate the new matrix for regularization
      fittedX <- tcrossprod(t(t(svd.res$U[,1:ncp,drop=FALSE]*row.w)*lambda.shrinked),svd.res$V[,1:ncp,drop=FALSE])
      fittedX <- fittedX/row.w
      # calculate the frobenious norm of the difference between iterations (convergence)
      # INV-OLR 
      fittedXRaw <-exp(t(t(bal)%*%t(fittedX)))
      fittedXRaw <- fittedXRaw/apply(fittedXRaw,1,sum)
      #
      diffRaw <- X/fittedXRaw
      diffRaw[missRaw] <- 1
      # OLR-coordinates
      diff <- t(bal%*%t(log(diffRaw)))
      #
      objective <- sum(diff^2*row.w)
      #      objective <- mean((Xhat[-missing]-fittedX[-missing])^2)
      # CONVERGENCE
      criterion <- abs(1 - objective/old)
      old <- objective
      nb.iter <- nb.iter + 1
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
        if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
      }
      if (nb.iter > maxiter) {
        nb.iter <- 0
        warning(paste("Stopped after ",maxiter," iterations"))
      }
    }
    # END LOOP WHILE
    
    # preparing the result
    if (scale) Xhat <- t(t(Xhat)*et)
    Xhat <- t(t(Xhat)+mean.p)
    # update X
    # INV-OLR
    X <-exp(t(t(bal)%*%t(Xhat)))
    X <- X/apply(X,1,sum)
    # completeObs
    completeObs<-Xaux/apply(Xaux,1,sum,na.rm=TRUE)
    completeObs[missRaw]<-X[missRaw]
    completeObs<-completeObs/apply(completeObs,1,sum)
    #
    if (scale) fittedX <- t(t(fittedX)*et)
    fittedX <- t(t(fittedX)+mean.p)
    # 
    # INV-OLR
    fittedXRaw <- exp(t(t(bal)%*%t(fittedX)))
    fittedXRaw <- fittedXRaw/apply(fittedXRaw,1,sum)
    # return complete matrix and imputed matrix
    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedXRaw
    return(result) 
  }
  
  method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
  method <- tolower(method)
  
  ## Preliminaries ----  
  
  obj <- Inf
  
  X <- as.data.frame(X,stringsAsFactors=TRUE)
  nn <- nrow(X); D <- ncol(X)
  
  X[X==label] <- NA
  X <- as.data.frame(apply(X,2,as.numeric),stringsAsFactors=TRUE)
  c <- apply(X,1,sum,na.rm=TRUE)
  
  # Check for closure
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
  
  # Sort columns decreasingly according to observed cells
  Xaux <- as.matrix(X)
  XauxClosed <- Xaux/apply(Xaux,1,sum,na.rm=TRUE) #as.matrix(X) # copy original data set
  # replace Zeros by NA
  Xna <- X # copy of original data
  #colnames(Xna)
  pz <- apply(apply(Xna,2,is.na),2,sum)/nn
  # ordered decreasengly by number of zeros (NA)
  X <- Xna[,order(-pz)]
  if (nrow(dl)==1) {dl <- matrix(dl[,order(-pz)],nrow=1)}
  else{
    dl <- dl[,order(-pz)]
  }

  # Balance matrix for olr
  Smat <- diag(rep(1,D))
  Smat[upper.tri(Smat)] <- -1
  Smat <- Smat[-D,]
  bal <- Smat
  numsbp <- dim(Smat)[1]
  for (f in 1:numsbp) {
    den <- sum(bal[f,]==-1)
    num <- sum(bal[f,]==1)  
    bal[f,bal[f,]==1] <- sqrt(den/((den+num)*num))
    bal[f,bal[f,]==-1] <- -sqrt(num/((den+num)*den))
  }
  
  # Build dl matrix for SVD imputation
  if (nrow(dl) == 1)  dl <- matrix(rep(1, nn), ncol = 1)%*%dl
  # set observ data as upper bound for estimates of observed values
  observedRaw <- which(!is.na(X))
  missingRaw <- which(is.na(X))
  Xaux2 <- as.matrix(X)
  # DL observed values = maximum value = observed = infinity upper bound
  #dl[observedRaw]<-Xaux2[observedRaw] 
  Xmax <- apply(X,2,max,na.rm=TRUE)
  Xmax <- matrix(rep(1, nn), ncol = 1)%*%Xmax
  dl[observedRaw] <- Xmax[observedRaw]
  colnames(dl) <- colnames(X)
  
  ## Imputation ---
  
  for (i in 1:nb.init){
    if (!any(is.na(X))) return(X)

    res.impute <- impute(X,dl,bal,ncp=ncp,scale=scale,method=method,threshold=threshold,
                         seed=if(!is.null(seed)){(seed*(i-1))}else{NULL},
                         init=i,maxiter=maxiter,row.w=row.w,coeff.ridge=coeff.ridge)
    diffRaw <- as.matrix(XauxClosed/res.impute$fittedX)
    diffRaw[missingRaw] <- 1
    # OLR-coordinates
    diff <-t(bal%*%t(log(diffRaw)))
    #
    res <- res.impute
    obj <- mean((diff)^2)
  }
  
  ## Final section ---
  
  # Re-scale to original units
   Y <- res.impute$completeObs
   XauxClosed <- XauxClosed[,order(-pz)]
   XauxClosed[missingRaw] <- Y[missingRaw]
   X <- XauxClosed*c
   
   # for (i in 1:nn) {
   #   if (any(is.na(X[i, ]))) {
   #     vbdl <- which(is.na(X[i, ]))
   #     # position with obs value
   #     Nvbdl<-which(!is.na(X[i, ]))[1]
   #     X[i, vbdl] <- (X[i, Nvbdl]/Y[i, Nvbdl]) * Y[i, vbdl]
   #   }
   # }
   
   
  if (closed == 1) {
    (X/apply(X,1,sum))*c[1]
  }
  
  # Original order
  X <- X[,colnames(Xna)]
  
  return(as.data.frame(X,stringsAsFactors=TRUE))
}


