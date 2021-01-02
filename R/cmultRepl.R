cmultRepl <- function(X,label=0,method=c("GBM","SQ","BL","CZM","user"),output=c("prop", "p-counts"),
                      frac=0.65,threshold=0.5,correct=TRUE,t=NULL,s=NULL,z.warning=0.8,
                      suppress.print=FALSE,delta=NULL)
{
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.vector(X) | is.character(X) | (nrow(X)==1)) stop("X must be a data matrix")
  if (!is.na(label)){
    if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as count zeros were found in the data set")
    if (any(is.na(X))) stop(paste("NA values not labelled as count zeros were found in the data set"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=T)) stop("Zero values not labelled as count zeros were found in the data set")
    if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
  }

  if (!missing("delta")){
    warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
    frac <- delta
  }
  
  X <- as.data.frame(X,stringsAsFactors=TRUE)
  X[X==label] <- NA
  
  checkNumZerosCol <- apply(X,2,function(x) sum(is.na(x)))
  if (any(checkNumZerosCol/nrow(X) == 1)) {
    stop(paste("Column(s) containing all zeros/unobserved values were found (check it out using zPatterns).",sep=""))
  }
  else{
    if (any(checkNumZerosCol/nrow(X) > z.warning)) {
      warning(paste("Column(s) containing more than ",z.warning*100,"% zeros/unobserved values were found (check it out using zPatterns).
                    (You can use the z.warning argument to modify the warning threshold).",sep=""))
    }
  }

  checkNumZerosRow <- apply(X,1,function(x) sum(is.na(x)))
  if (any(checkNumZerosRow/ncol(X) == 1)) {
    stop(paste("Row(s) containing all zeros/unobserved values were found (check it out using zPatterns).",sep=""))
  }
  else{
    if (any(checkNumZerosRow/ncol(X) > z.warning)) {
      warning(paste("Row(s) containing more than ",z.warning*100,"% zeros/unobserved values were found (check it out using zPatterns).
                  (You can use the z.warning argument to modify the warning threshold).",sep=""))
    }
  }
  
  N <- nrow(X); D <- ncol(X)
  n <- apply(X,1,sum,na.rm=TRUE)
  
  # Determining t and s
  
  method <- match.arg(method)
  output <- match.arg(output)
  
  if (method!="CZM"){
    if (method=="user") {t <- t}
    else {
      alpha <- matrix(0,nrow=N,ncol=D)
      for (i in 1:N){
        alpha[i,] <- apply(X,2,function(x) sum(x[-i],na.rm=T))
        }
      t <- alpha/rowSums(alpha)
      if ((method=="GBM") && (any(t==0))) {stop("GBM method: not enough information to compute t hyper-parameter,
                                                probably there are columns with < 2 positive values.")}
    }
    s <- switch(method,
                GBM = 1/apply(t,1,function(x) exp(mean(log(x)))),
                SQ = sqrt(n),
                BL = D,
                user = s)  
    repl <- t*(s/(n+s))
  }
  
  if (method=="CZM"){
    repl <- frac*matrix(1,ncol=D,nrow=N)*(threshold/n)
  }
  
  # Multiplicative replacement on the closed data
  
  X2 <- t(apply(X,1,function(x) x/sum(x,na.rm=T)))
  colmins <- apply(X2,2,function(x) min(x,na.rm=T))
  corrected <- 0
  
  for (i in 1:N){
    if (any(is.na(X2[i,]))){
      z <- which(is.na(X2[i,]))
      X2[i,z] <- repl[i,z]
      if (correct==TRUE){
        if (any(X2[i,z] > colmins[z])){
          f <- which(X2[i,z] > colmins[z])
          X2[i,z][f] <- frac*colmins[z][f]
          corrected <- corrected + length(f)
        }
      }
      X2[i,-z] <- (1-(sum(X2[i,z])))*X2[i,-z]
    }
  }        
  
  # Rescale to p-counts
  
  if (output=="p-counts"){
    for (i in 1:N){
      
      if (any(is.na(X[i,]))){
        zero <- which(is.na(X[i,]))
        pos <- setdiff(1:D,zero)[1]
        X[i,zero] <- (X[i,pos]/X2[i,pos])*X2[i,zero]
      }
    }
  }
  
  res <- ifelse(output=="prop", X2, X)
  if (suppress.print == FALSE){
    if ((correct==TRUE) & (corrected > 0)) {cat(paste("No. corrected values: ",corrected,"\n"))}
  }
  return(as.data.frame(res,stringsAsFactors=TRUE))
}
