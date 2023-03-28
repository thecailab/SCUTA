#' plot function, which can be used to show the trajectory of outcome given an independent variable holding other variable as 0s for a gene
#' param counts, count matrix, in which the rows represent genes, while columns represent cells
#' param coldata, clinical information of each cell, e.g., group, individual, time, etc.
#' param fit.res, model fitting resutls by using SCUTA
#' param fit.res.coef, coefficient name in fit.res that will be the independent variable in the plot
#' param ind.var, independent variable name shown in the coldata which will be used as the x axis in the plot
#' param genenames, a vector contain the gene names, the log2 count of which will be plot
SCUTA.plot.facotr<-function(counts,
                            coldata,
                            fit.res,
                            fit.res.coef,
                            ind.var,
                            genenames=NULL
){
  # if the genenames is not specified
  if (length(genenames)<1){stop("Gene names are not specified!")}
  
  # if there gene names not in the fit.res
  if (  length(which( is.na(match(genenames,rownames(coef(fit.res)))))) >=1  ){
    tmp.gene.names<-genenames[which( is.na(match(genenames,rownames(coef(fit.res)))))]
    tmp.gene.names<-unique(tmp.gene.names)
    tmp.gene.names.list<-tmp.gene.names[1]
    if (length(tmp.gene.names)>=2){
      for (i in 2:length(tmp.gene.names)){
        tmp.gene.names.list<-paste0(tmp.gene.names.list," ,",tmp.gene.names[i])
      }
    }
    stop(paste0( tmp.gene.names.list, " are not found in the results"))
  }
  
  coefs<-coef(fit.res)
  if ( is.na(match("(Intercept)",colnames(coefs) )) ){
    coefs<-cbind("(Intercept)"=rep(0,nrow(coefs)),coefs)
  }
  coefs<-coefs[,c("(Intercept)",fit.res.coef)]
  
  value<-substr(fit.res.coef,nchar(ind.var)+1,nchar(fit.res.coef))
  x<-ifelse(coldata[,ind.var]==value,1,0)
  
  if (length(genenames)==1){
    subset.count<-as.matrix(counts[which(rownames(counts) %in% genenames),])
    y=log2(subset.count+1)
    colnames(y)<-genenames
    inte.data<-cbind(coldata,y)
    
    y.max<-max(y)
    col.list<-seq(1,length(genenames))
    plot(x,inte.data[,genenames[1]],col=col.list,ylim=c(0,y.max),ylab="log2 Count",xlab=fit.res.coef)
    subset.coef<-coefs[genenames[1],]
    abline(coef = subset.coef,lwd=2,col=col.list)
    legend("topright",legend=genenames,col = col.list, lty = 1,lwd=2)
  }
  
  if (length(genenames)>1){
    subset.count<-t(counts[which(rownames(counts) %in% genenames),])
    y=log2(subset.count+1)
    colnames(y)<-genenames
    inte.data<-cbind(coldata,y)
    
    col.list<-seq(1,length(genenames))
    y.max<-max(y)
    plot(x,inte.data[,genenames[1]],col=col.list[1],ylim=c(0,y.max),ylab="log2 Count",xlab=fit.res.coef)
    subset.coef<-coefs[genenames[1],]
    abline(coef = subset.coef,lwd=2,col=col.list[1])
    # par(new=TRUE)
    for (i in 2:length(genenames)){
      points(x,inte.data[,genenames[i]],col=col.list[i])
      subset.coef<-coefs[genenames[i],]
      abline(coef = subset.coef,lwd=2,col=col.list[i])
    }
    
    if (length(genenames)>10){
      warning("the number of genes is greater than 10, legend will only show the first 10 genes")
      legend("topright",legend=genenames[1:10],col = col.list[1:10], lty = 1,lwd=2)
    }else {
      legend("topright",legend=genenames,col = col.list, lty = 1,lwd=2)
    }
    
  }
  
  
}

#SCUTA.plot.facotr(counts=counts,
#                  coldata=coldata,
#                 fit.res=fit.res,
#                  fit.res.coef="condition2",
#                  ind.var="condition",
#                  genenames=c("Gene1","Gene2"))

