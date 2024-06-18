lrDA <-
  function(X,label=NULL,dl=NULL,ini.cov=c("lrEM","complete.obs","multRepl"),frac=0.65,
           imp.missing=FALSE,n.iters=1000,m=1,store.mi=FALSE,closure=NULL,z.warning=0.8,
           z.delete=TRUE,delta=NULL){
    
    if (any(X<0, na.rm=T)) stop("X contains negative values")

    if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
      if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored or missing values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored or missing values were found in the data set")
      if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    }
    
    if (imp.missing==FALSE){
      if (is.character(dl)) stop("dl must be a numeric vector or matrix")
      if (is.null(dl)){ # If dl not given use min per column
        dl <- apply(X,2, function(x) min(x[x!=label]))
        warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
      }
      if (is.vector(dl)) dl <- matrix(dl,nrow=1)
      dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
    }
    
    if (imp.missing==FALSE){
      if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
      if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
    }
    
    if ((store.mi==TRUE) & (m==1)) store.mi <- FALSE
    
    if (!missing("delta")){
      warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
      frac <- delta
    }
    
    lm.sweep <- function(M,C,varobs){
      
      sweep.matrix <- function(A,ind){ 
        
        nn <- nrow(A); p <- ncol(A)
        S <- A
        
        for (j in ind){                       
          S[j,j] <- -1/A[j,j]     
          for (i in 1:p) {           
            if (i != j){
              S[i,j] <- -A[i,j]*S[j,j]
              S[j,i] <- S[i,j]
            }
          }
          for (i in 1:p){           
            if (i != j){
              for (k  in 1:p){
                if (k != j){
                  S[i,k] <- A[i,k] - S[i,j]*A[j,k]
                  S[k,i] <- S[i,k]
                }
              }
            }
          }
          A <- S 
        }
        return(A)
      }  
      
      p <- length(M)                      
      q <- length(varobs)                     
      i <- rep(1,p)                  
      i[varobs] <- i[varobs]-1
      dep <- which(i!=0)
      ndep <- length(dep)                 
      
      A <- matrix(0,p+1,p+1)             
      A[1,1] <- -1                    
      A[1,2:(p+1)] <- M                
      A[2:(p+1),1] <- matrix(M,ncol=1)
      A[2:(p+1),2:(p+1)] <- C             
      
      reor <- c(1,varobs+1,dep+1)     
      A <- A[reor,reor]              
      A <- sweep.matrix(A,1:(q+1))             
      
      B <- A[1:(q+1),(q+2):(p+1)]       
      CR <- A[(q+2):(p+1),(q+2):(p+1)]       
      
      return(list(betas=B,resid=CR))
    }
    
    inv.raw <- function(Y,X,pos,closed,nn,c){
      
      inv.alr <- function(x,pos){
        
        ad<-1/(rowSums(exp(x))+1)
        ax<-exp(x)*ad
        if(pos==1) {
          a<-cbind(ad,ax,stringsAsFactors=TRUE)
        }
        else { 
          if (dim(x)[2] < pos){
            a<-cbind(ax,ad,stringsAsFactors=TRUE)
          }   
          else {
            a<-cbind(ax[,1:(pos-1)],ad,ax[,pos:(dim(x)[2])],stringsAsFactors=TRUE)
          }
        }
        return(a)
      }
      
      Y <- inv.alr(Y,pos)
      
      for (i in 1:nn){
        if (any(is.na(X[i,]))){
          vbdl <- which(is.na(X[i,]))
          X[i,vbdl] <- (X[i,pos]/Y[i,pos])*Y[i,vbdl]
        }
      }    
      if (closed==1){
        X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
      }
      return(as.data.frame(X,stringsAsFactors=TRUE))
    }
    
    riwish <- function(v,S){ # From ratematrix package
      S <- solve(S)
      if (!is.matrix(S)) S <- matrix(S)
      if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
      }
      p <- nrow(S)
      CC <- chol(S)
      Z <- matrix(0, p, p)
      diag(Z) <- sqrt(stats::rchisq(p, v:(v - p + 1)))
      if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- stats::rnorm(p *(p - 1)/2)
      }
      out <- crossprod(Z %*% CC)
      return(solve(out))
    }
    
    ini.cov <- match.arg(ini.cov)
    
    X <- as.data.frame(X,stringsAsFactors=TRUE)
    nn <- nrow(X); p <- ncol(X)
    if (nn <= p) stop("The lrDA algorithm works on regular data sets (no. rows > no. columns). You can consider lrSVD for wide dat sets.")
    
    rnames <- rownames(X)
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    rownames(X) <- rnames
    c <- apply(X,1,sum,na.rm=TRUE)
    
    checkNumZerosCol <- apply(X, 2, function(x) sum(is.na(x)))
    
    if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
      cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
      if (z.delete == TRUE) {
        if (length(cases) > (ncol(X)-2)) {
          stop(paste("Almost all columns contain >", z.warning*100,
                     "% zeros/unobserved values (see arguments z.warning and z.delete).",
                     sep=""))
        }      
        X <- X[,-cases]      
        action <- "deleted"
        
        warning(paste("Column no. ",cases," containing >", z.warning*100,
                      "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                      sep=""))
      } else {      
        action <- "found"      
        warning(paste("Column no. ",cases," containing >", z.warning*100,
                      "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                      sep=""))      
      }
    }
    
    checkNumZerosRow <- apply(X, 1, function(x) sum(is.na(x)))  
    if (any(checkNumZerosRow/ncol(X) > z.warning)) {    
      cases <- which(checkNumZerosRow/ncol(X) > z.warning)    
      if (z.delete == TRUE) {
        if (length(cases) > (nrow(X)-2)) {
          stop(paste("Almost all rows contain >", z.warning*100,
                     "% zeros/unobserved values (see arguments z.warning and z.delete).",
                     sep=""))
        }
        X <- X[-cases,]      
        action <- "deleted"
        
        warning(paste("Row no. ",cases," containing >", z.warning*100,
                      "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                      sep=""))
      } else {      
        action <- "found"      
        warning(paste("Row no. ", cases," containing >", z.warning*100,
                      "% zeros/unobserved values ", action,
                      " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                      sep=""))      
      }
    }
    
    if (imp.missing==FALSE) {if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl}
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1
    
    pos <- which(!is.na(colSums(X)))[1]
    if (is.na(pos)) stop("lrDA requires at least one fully observed column")
    
    if (imp.missing==FALSE){
      cpoints <- log(dl)-log(X[,pos])-.Machine$double.eps
      cpoints <- cpoints[,-pos]
    }
    
    X_alr <- log(X)-log(X[,pos]); X_alr <- as.matrix(X_alr[,-pos])
    nn <- nrow(X_alr); p <- ncol(X_alr)
    
    if (ini.cov == "complete.obs"){
      if (inherits(try(solve(cov(X,use=ini.cov)),silent=TRUE),"try-error"))
        stop("ini.cov: too few complete cases for using 'complete.obs'")
      M <- matrix(colMeans(X_alr,na.rm=T),ncol=1)
      C <- cov(X_alr,use=ini.cov)}
    if (ini.cov == "multRepl"){
      X.mr <- multRepl(X,label=NA,dl=dl,frac=frac,imp.missing=imp.missing,closure=closure,z.warning=z.warning,z.delete=z.delete)
      X.mr_alr <- t(apply(X.mr,1,function(x) log(x)-log(x[pos])))[,-pos]
      M <- matrix(colMeans(X.mr_alr,na.rm=T),ncol=1)
      C <- cov(X.mr_alr)}
    if (ini.cov == "lrEM"){
      X.em <- lrEM(X,label=NA,dl=dl,ini.cov="multRepl",frac=frac,imp.missing=imp.missing,closure=closure,suppress.print=TRUE,
                   z.warning=z.warning,z.delete=z.delete)
      X.em_alr <- t(apply(X.em,1,function(x) log(x)-log(x[pos])))[,-pos]
      M <- matrix(colMeans(X.em_alr,na.rm=T),ncol=1)
      C <- cov(X.em_alr)}
    
    misspat <- as.data.frame(is.na(X)*1,stringsAsFactors=TRUE)
    misspat <- as.factor(do.call(paste,c(misspat,sep="")))
    levels(misspat) <- 1:(length(levels(misspat)))
    
    t <- 0
    k <- 0
    runs <- 0
    alt.in <- FALSE
    alt.pat <- 0
    alt.mr <- 0
    
    if (m > 1){
      imputed <- matrix(0,nrow=m,ncol=sum(is.na(X_alr)))
      if (store.mi==TRUE) mi.list <- vector(mode="list",m)
    }
    
    while (t <= n.iters*m){
      
      Y <- X_alr                              
      runs <- runs + 1
      
      # I-step
      
      for (npat in 1:length(levels(misspat))){                     
        i <- which(misspat==npat) 
        varmiss <- which(is.na(X_alr[i[1],]))
        if (length(varmiss) == 0) {next} # Skip first pattern if all obs
        varobs <- which(!is.na(X_alr[i[1],]))
        if (length(varobs) == 0){
          alt.in <- TRUE
          temp <- multRepl(X[i,,drop=FALSE],label=NA,dl=dl[i,,drop=FALSE],frac=frac,imp.missing=imp.missing,closure=closure,z.warning=z.warning,z.delete=z.delete)
          Y[i,] <- t(apply(temp,1,function(x) log(x)-log(x[pos])))[,-pos]
          if (runs == 1){
            alt.pat <- c(alt.pat,npat)
            alt.mr <- list(alt.mr,i)
          }
          break
        }
        sigmas <- matrix(0,ncol=p)
        B <- matrix(lm.sweep(M,C,varobs)[[1]],ncol=length(varmiss))
        CR <- lm.sweep(M,C,varobs)[[2]]
        Y[i,varmiss] <- matrix(1,nrow=length(i))%*%B[1,] + X_alr[i, varobs, drop=FALSE]%*%B[2:(length(varobs)+1),]
        sigmas[varmiss] <- sqrt(diag(as.matrix(CR)))
        if (imp.missing==FALSE){
          for (j in 1:length(varmiss)){                                
            sigma <- sigmas[varmiss[j]]
            Y[i,varmiss[j]] <- rtruncnorm(1,-Inf,cpoints[i,varmiss[j]],Y[i,varmiss[j]],sigma)
          }
        }
        if (imp.missing==TRUE){
          for (j in 1:length(varmiss)){                                
            sigma <- sigmas[varmiss[j]]
            Y[i,varmiss[j]] <- rnorm(1,Y[i,varmiss[j]],sigma)
          }
        }
      }
      
      if ((t%in%((1:m)*n.iters)) & (m > 1)){
        k <- k + 1
        imputed[k,] <- Y[which(is.na(X_alr))]
        if (store.mi==TRUE){
         mi.list[[k]] <- Y
        }
      }

      # P-step
      
      C <- riwish(nn-1,nn*cov(Y))
      M <- mvrnorm(1,colMeans(Y),(1/nn)*C)
      
    t <- t + 1

    }
    
    if ((m > 1) & (store.mi == FALSE)) Y[which(is.na(X_alr))] <- colMeans(imputed) # MI estimates
    
    if (store.mi==FALSE) X <- inv.raw(Y,X,pos,closed,nn,c)
    
    if (store.mi==TRUE) X <- lapply(mi.list,FUN=function(x) inv.raw(x,X,pos,closed,nn,c))
    
    if (alt.in) {
      if (imp.missing==FALSE){
      cat("Warning: samples with only one observed component were found \n")
        for (i in 2:length(alt.pat)){
          cat(paste("  Pattern no.",alt.pat[i],"was imputed using multiplicative simple replacement \n"))
          cat("   Affected samples id: "); cat(alt.mr[[i]]); cat("\n\n")
        }
      }
    }
    
    return(X)
  }
