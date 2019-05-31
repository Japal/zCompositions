multReplus <- function(X, dl = NULL, delta = 0.65, suppress.print = FALSE,
                       closure = NULL){
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)

  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")

  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  
  if (any(is.na(X))==FALSE) stop("No missing data were found in the data set")
  if (any(X==0, na.rm=T)==FALSE) stop("No zeros were found in the data set")
  
  
  gm <- function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x[!is.na(x)]))
  }
  
  nam <- NULL
  if (!is.null(names(X))) nam <- names(X)
  if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)))
  
  ## Preliminaries ----

  X <- as.data.frame(X)
  nn <- nrow(X); D <- ncol(X)
  X <- as.data.frame(apply(X,2,as.numeric))
  c <- apply(X,1,sum,na.rm=TRUE)

  if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl

  # Check for closure
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1

  Y <- X
  
  if (!is.null(closure)){
    if (closed == 1) {stop("closure: The data are already closed to ",c[1])}
    resid <- apply(X,1, function(x) closure-sum(x, na.rm = TRUE))
    Xresid <- cbind(X,resid)
    c <- rep(closure,nn)
    Y <- Xresid
  }
  
  ## Imputation of missing (ignoring 0s in column if any) ----
  
  gms <- apply(X,2,function(x) gm(x[x!=0]))
  for (i in 1:nn){
    if (any(is.na(X[i,]))){
      z <- which(is.na(X[i,]))
      Y[i,z] <- gms[z]
      if (!is.null(closure)){
        Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(Xresid[i,-z]))*Xresid[i,-z]
        tmp <- Y[i,-(D+1)]
        nz_idx <- which(tmp[-z]!=0)[1] # Use a non-zero part for adjustment
        X[i,z] <- as.numeric((X[i,-z][nz_idx]/tmp[-z][nz_idx]))*Y[i,z]
      }
      else{
        Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(X[i,-z]))*X[i,-z]
        nz_idx <- which(Y[i,-z]!=0)[1] # Use a non-zero part for adjustment
        X[i,z] <- as.numeric((X[i,-z][nz_idx]/Y[i,-z][nz_idx]))*Y[i,z]
      }
    }
  } 
  
  if (closed==1){
    X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
  } 
  
  if (any(X < 0, na.rm = T)) stop("multRepl: negative imputed values were generated (please check out help for advice)")
  
  ## Imputation of zeros ----
  
  X <- multRepl(X,label=0,dl=dl,delta=delta,closure=closure)
  
  ## Final section ----
  
  if (!is.null(nam)) names(X) <- nam
  
  return(as.data.frame(X))
  
}