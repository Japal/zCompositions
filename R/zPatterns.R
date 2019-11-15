zPatterns <- function(X,label=NULL,plot=TRUE,
                      axis.labels=c("Component","Pattern ID"),
                      bar.ordered=as.character(c(FALSE,FALSE)),
                      bar.colors=c("red3","red3"), bar.labels=FALSE, 
                      show.means=FALSE,round.means=2, cex.means=1,
                      type.means=c("cgm","am"),
                      cell.colors=c("dodgerblue","white"), cell.labels=c(label,paste("No",label)),
                      cex.axis=1.1, grid.color="black",
                      grid.lty="dotted", legend=TRUE, suppress.print=FALSE, ...){
  
  cgm <- function(X, round.means = round.means)
  {
    
    ms <- apply(X,2,function(x){
      if (all(is.na(x))) {x <- NA} else {exp(mean(log(x),na.rm=T))}
    })
    ms[is.na(ms)] <- 0
    round(ms/sum(ms)*100,round.means)
    
  }
  
  am <- function(X, round.means = round.means)
  {
    
    ms <- apply(X,2,function(x){
      if (all(is.na(x))) {x <- NA} else {mean(x,na.rm=T)}
    })
    round(ms,round.means)
    
  }
  
  plot.patterns <- function(a, show.means = show.means, round.means = round.means,
                            cex.means=cex.means, X = X, pat = pat.ID, ...)
  {
    
    zones <- matrix(c(2,4,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1.5/5), heights=c(2/5,3.5/5))
    par(mar=c(3,3,0.5,0.5))
    a <- as.matrix(a[rev(rownames(a)),])
    image(1:ncol(a),1:nrow(a),t(a),col=rev(cell.colors),axes=F)
    mtext(side=1,text=axis.labels[1],line=1.75)
    mtext(side=2,text=axis.labels[2],line=1.75)
    par(mgp=c(3, .3, 0))
    axis(side = 1,at = seq(1,ncol(a),by=1),labels=colnames(a),tck=0,cex.axis=cex.axis)
    axis(side = 2,at = seq(1,nrow(a),by=1),labels=rownames(a),las=2,tck=0,cex.axis=cex.axis)
    box()
    grid(ncol(a),nrow(a),col=grid.color,lty=grid.lty)
    if (show.means == TRUE){
      if (type.means == "cgm"){
        cgmp <- by(X,pat.ID,cgm,round.means=round.means)
        cgmp <- cgmp[pat.ID2] # Reorder as in patterns table
        cgmp <- do.call(rbind,cgmp)
        cgmp <- cgmp[nrow(a):1,]
        cgmp[cgmp==0] <- NA
        for (i in 1:nrow(a)){
          for (j in 1:ncol(a)){
            text(j,i,label=cgmp[i,j],cex=cex.means)
          }
        }
      }
      if (type.means == "am"){
        amp <- by(X,pat.ID,am,round.means=round.means)
        amp <- amp[pat.ID2]
        amp <- do.call(rbind,amp)
        amp <- amp[nrow(a):1,]
        
        for (i in 1:nrow(a)){
          for (j in 1:ncol(a)){
            text(j,i,label=amp[i,j],cex=cex.means)
          }
        }
      }
    }
    par(mar=c(0,3.25,1,0.75))
    a <- barplot(as.vector(prop.col),axes=F,col=bar.colors[1],xaxs="i",
                 ylim=c(0,max(as.vector(prop.col)+0.2*max(as.vector(prop.col)))))
    if (bar.labels==TRUE)
      text(a,as.vector(prop.col),labels=as.vector(prop.col),cex=0.85,pos=3)
    par(mar=c(3.25,0,0.75,0.75))
    a <- barplot(rev(as.vector(pat.freq)),horiz=T,axes=F,col=bar.colors[2],yaxs="i",
                 xlim=c(0,max(as.vector(pat.freq)+0.3*max(as.vector(pat.freq)))))
    if (bar.labels==TRUE)
      text(rev(as.vector(pat.freq)),a,labels=rev(as.vector(pat.freq)),cex=0.85,pos=4)
    par(mar=c(0,0,3,0))
    plot.new()
    if (legend==TRUE){
      if (any(is.na(cell.labels))) cell.labels[is.na(cell.labels)] <- "NA" 
      legend("topleft",cell.labels,pch=c(22,22),bty="n",
             pt.bg=cell.colors,pt.cex=2,cex=1.1)}
  } 
  
  type.means <- match.arg(type.means)
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")  
  if (is.vector(X)) stop("X must be a matrix or data.frame class object")
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=T))
      warning("Unidentified zero values were found and will be ignored")
    if (any(is.na(X))) warning(paste("Unidentified NA values were found in the data set and will be ignored"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=T)) warning("Unidentified zero values were found in the data set and will be ignored")
    if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
  }
  
  X <- as.data.frame(X)
  
  n <- nrow(X); p <- ncol(X)
  
  if (is.na(label)) miss <- as.data.frame(is.na(X)*1)
  else miss <- as.data.frame((X==label)*1)
  
  miss[is.na(miss)] <- 0 # Ignore any unlabelled NAs/zeros to graph patterns
  miss <- cbind(miss,pat=do.call(paste,c(miss,sep="")))
  tmp <-data.frame(Pattern=names(table(miss$pat)),ID=1:nlevels(miss$pat))
  pat.ID <- as.factor(tmp$ID[match(miss$pat,tmp$Pattern)]) # IDs in original row order
  miss <- miss[order(miss$pat),] # Order according to patterns

  pat.freq <- round(as.vector(table(miss$pat))/n*100,2) # Rel freq of ordered patterns
  prop.col <- round(colSums(miss[,1:p])/n*100,2)
  prop <- round(sum(miss[,1:p])/(n*p)*100,2)
  
  tab <- miss[!duplicated(miss),1:p]
  pat.ID2 <- levels(pat.ID)
  rownames(tab) <- pat.ID2
  
  if (bar.ordered[1]==TRUE){
    tab <- tab[order(pat.freq,decreasing=T),]
    pat.ID2 <- pat.ID2[order(pat.freq,decreasing=T)]
    pat.freq <- pat.freq[order(pat.freq,decreasing=T)]
  }
  
  if (bar.ordered[2]==TRUE){
    tab <- tab[,order(prop.col,decreasing=T)]
    X <- X[,order(prop.col,decreasing=T)]
    prop.col <- prop.col[order(prop.col,decreasing=T)]
  }
  
  tab.num <- tab
  tab[tab==1] <- "+"
  tab[tab==0] <- "-"
  
  # Summary
  
  X[X==label] <- NA; X[X==0] <- NA # Ignore labelled/NAs/zeros for summaries
  
  if (plot==TRUE) plot.patterns(tab.num, show.means=show.means,round.means=round.means,
                                cex.means=cex.means,X = X, pat = pat.ID2, ...)
  tab <- cbind(Patt.ID=pat.ID2,
               tab,
               No.Unobs=rowSums(tab[,1:p]=="+"),
               Patt.Perc=pat.freq)
  
  if (suppress.print==FALSE){
    cat("Patterns ('+' means ",cell.labels[1],", '-' means ",cell.labels[2],") \n\n",sep="")
    print(tab,row.names=FALSE)
    cat("\n")
    cat("Percentage cells by component \n")
    print(prop.col)
    cat("\n")
    cat(paste("Overall percentage cells: ",prop,"% \n",sep=""))
  }
  invisible(pat.ID)
}