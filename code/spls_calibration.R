#################### Functions ####################

SetToCoordinate=function(matrix, rown, coln, value){
  matrix[rown+nrow(matrix)*(coln-1)]=value
  return(matrix)
}


PlotCalibsPLS=function(res, sparsity='X', ncomp=1, main=NULL){
  if (sparsity%in%c('X', 'Y')){
    plot(res$MSEP[res$NComp==ncomp],
         type="l", xaxt="n", xlab=paste0("Number of variables in ", sparsity), ylab="MSEP",
         main=ifelse(is.null(main), yes=paste0("Calibration of the number of variables in ", sparsity), no=main),
         sub=paste("Number of components =", ncomp), las=1)
    if (sparsity=="X"){
      axis(side = 1, at = 1:sum(res$NComp==ncomp), labels = res$NVarX[res$NComp==ncomp])
    }
    if (sparsity=="Y"){
      axis(side = 1, at = 1:sum(res$NComp==ncomp), labels = res$NVarY[res$NComp==ncomp])
    }
    points(which.min(res$MSEP[res$NComp==ncomp]),
           res$MSEP[res$NComp==ncomp][which.min(res$MSEP[res$NComp==ncomp])],
           pch=19, col="tomato")
  } else {
    Tmp=matrix(0, nrow=length(unique(res$NVarX)), ncol=length(unique(res$NVarY)))
    Tmp=SetToCoordinate(Tmp, res$NVarX, res$NVarY, res$MSEP)
    
    rownames(Tmp)=unique(res$NVarX)
    colnames(Tmp)=unique(res$NVarY)
    pheatmap(Tmp, cluster_rows = F, cluster_cols = F)
  }
}


UnivariatePLS <- function(X, Y, nrepeat=1, 
                          stepX=1, stepY=1, MinVarX=1, MaxVarX=NULL){
  # get max and min var
  MaxVarX <- ifelse(is.null(MaxVarX), yes=ncol(X), no=MaxVarX)
  MinVarY <- 1
  MaxVarY <- 1

  v = expand.grid(unique(c(MinVarY,seq(MinVarY-1, MaxVarY, by=stepY)[-1])),
                  unique(c(MinVarX, seq(MinVarX-1, MaxVarX, by=stepX)[-1])),
                  seq(1, ncol(Y), 1))

  no_cores <- detectCores()
  clust <- makeCluster(no_cores, type='FORK')

  Summary=t(parApply(clust, v, MARGIN=1, FUN=function(a){
    j=a[1];i=a[2];h=a[3]
    TmpsPLS <- spls(X,Y[,h],keepX=i, keepY=j,ncomp=1,mode='regression')
    TmpMSEP=TmpQ2.total=NULL
    for (k in 1:nrepeat){
      TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1,progressBar=FALSE)
      TmpMSEP=cbind(TmpMSEP, TmpPerf$MSEP[,1])
      TmpQ2.total=cbind(TmpQ2.total, TmpPerf$Q2.total[1, ])
    }
    return(c(h, i, j, mean(TmpMSEP), mean(TmpQ2.total)))
  }))

  stopCluster(clust)
  Summary=as.data.frame(Summary)
  colnames(Summary)=c('y', 'NVarX', 'NVarY', 'MSEP', 'Q2')
  return(Summary)
}


MultivariatePLS <- function(X, Y, nrepeat=1, 
                            stepX=1, stepY=1, MinVarX=1, MaxVarX=NULL){
  # get max and min var
  MaxVarX <- ifelse(is.null(MaxVarX), yes=ncol(X), no=MaxVarX)
  MinVarY <- ncol(Y)
  MaxVarY <- ncol(Y)
  
  v = expand.grid(unique(c(MinVarY,seq(MinVarY-1, MaxVarY, by=stepY)[-1])),
                  unique(c(MinVarX, seq(MinVarX-1, MaxVarX, by=stepX)[-1])))
  
  no_cores <- detectCores()
  clust <- makeCluster(no_cores, type='FORK')
  
  Summary=t(parApply(clust, v, MARGIN=1, FUN=function(a){
    j=a[1];i=a[2];
    TmpsPLS <- spls(X,Y,keepX=i, keepY=j,ncomp=1,mode='regression')
    TmpMSEP=TmpQ2.total=NULL
    for (k in 1:nrepeat){
      TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1,progressBar=FALSE)
      TmpMSEP=cbind(TmpMSEP, TmpPerf$MSEP[,1])
      TmpQ2.total=cbind(TmpQ2.total, TmpPerf$Q2.total[1, ])
    }
    return(c(i, j, mean(TmpMSEP), mean(TmpQ2.total)))
  }))
  
  stopCluster(clust)
  Summary=as.data.frame(Summary)
  colnames(Summary)=c('NVarX', 'NVarY', 'MSEP', 'Q2')
  return(Summary)
}


StabilityPLS<-function(X, Y, nrepeat=1, stepX=1, stepY=1, MinVarX=1, MaxVarX=NULL){

  # get max and min var
  MaxVarX <- ifelse(is.null(MaxVarX), yes=ncol(X), no=MaxVarX)
  MinVarY <- ncol(Y)
  MaxVarY <- ncol(Y)

  no_cores <- detectCores()
  clust <- makeCluster(no_cores, type='FORK')

  v=expand.grid(unique(c(MinVarY,seq(MinVarY-1, MaxVarY, by=stepY)[-1])), 
                unique(c(MinVarX, seq(MinVarX-1, MaxVarX, by=stepX)[-1])))

  Summary=t(parApply(clust, v, MARGIN=1, FUN=function(k){
    j=k[1];i=k[2]
    selected_X <- matrix(rep(0, len=ncol(X)), ncol = 1, dimnames = list(colnames(X), "comp1"))
    for (k in 1:nrepeat){
      sample_rows <- sample(nrow(X),size=floor(nrow(X) * 0.8),replace=FALSE)
      TmpsPLS <- spls(X[sample_rows,],Y[sample_rows,],keepX=c(i), keepY=c(j),ncomp=1,mode='regression')
      X_loadings <- TmpsPLS$loadings$X
      X_loadings[X_loadings != 0] = 1
      selected_X = selected_X + X_loadings
    }
    selected_X <- selected_X / nrepeat
    return(c(i, t(selected_X)))
  }))
  
  stopCluster(clust)
  Summary=as.data.frame(Summary)
  colnames(Summary)=c('NVarX', colnames(X))
  return(Summary)
}
