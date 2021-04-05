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


ComponentPLS <- function(X, Y, comp, nrepeat, stepX, stepY,
                         MinVarX, MinVarY, MaxVarX, MaxVarY,
                         SelectedX, SelectedY){
  ## Run the PLS for the nth component
  
  no_cores <- detectCores()-1
  clust <- makeCluster(no_cores, type='FORK')
  
  v=expand.grid(unique(c(MinVarY,seq(MinVarY-1, MaxVarY, by=stepY)[-1])), 
                unique(c(MinVarX, seq(MinVarX-1, MaxVarX, by=stepX)[-1])))
  Summary=t(parApply(clust, v, MARGIN=1, FUN=function(k){
    j=k[1];i=k[2]
    TmpKeepX <- c(SelectedX,i)
    TmpKeepY <- c(SelectedY,j)
    TmpsPLS <- spls(X,Y,keepX=TmpKeepX, keepY=TmpKeepY,ncomp=comp,mode='regression')
    TmpMSEP=TmpQ2.total=NULL
    for (k in 1:nrepeat){
      TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1,progressBar=FALSE)
      TmpMSEP=cbind(TmpMSEP, TmpPerf$MSEP[,comp])
      TmpQ2.total=cbind(TmpQ2.total, TmpPerf$Q2.total[comp, ])
    }
    return(c(comp, i, j, mean(TmpMSEP), mean(TmpQ2.total)))
  }))
  
  stopCluster(clust)
  Summary=as.data.frame(Summary)
  colnames(Summary)=c('NComp', 'NVarX', 'NVarY', 'MSEP', 'Q2')
  return(Summary)
}

NarrowPLS <- function(Summary, X, Y, comp, nrepeat, sparsity, 
                      stepX, stepY, MinVarX, MaxVarX,
                      MinVarY, MaxVarY, narrow_stepX, narrow_stepY){
  # Extension to run another search on a narrower range for around the previously
  # identified lowest MSEP variable
  print(paste0(comp, "_narrow"))
  
  Xrank <- Summary$NVarX + rank(Summary$MSEP)
  previousX=as.numeric(Summary[which.min(Xrank),]$NVarX)
  
  Yrank <- Summary$NVarY + rank(Summary$MSEP)
  previousY=as.numeric(Summary[which.min(Yrank),]$NVarY)
  
  if (!is.null(narrow_stepX)){
    narrowX = c(pmax(MinVarX, previousX - (stepX/2)),
                pmin(MaxVarX, previousX + (stepX/2)),
                narrow_stepX)
  }else{
    narrowX = c(previousX, previousX, 1)
  }
  
  if (!is.null(narrow_stepY)){
    narrowY = c(pmax(MinVarY, previousY - (stepY/2)),
                pmin(MaxVarY, previousY + (stepY/2)),
                narrow_stepY)
  }else{
    narrowY = c(previousY, previousY, 1)
  }
  
  Summary <- ComponentPLS(X, Y, comp, nrepeat,
                          stepX=narrowX[3], stepY=narrowY[3],
                          MinVarX=narrowX[1], MinVarY=narrowY[1],
                          MaxVarX=narrowX[2], MaxVarY=narrowY[2],
                          SelectedX, SelectedY)
  return(Summary)
}


select_num_vars <- function(Summary, var){
  score <- scale(Summary[,var]) + scale(Summary$MSEP)
  return(as.numeric(res[which.min(score),][,var]))
}

CalibratesPLS<-function(X, Y, nComp=1, nrepeat=1, sparsity='X', 
                        stepX=1, stepY=1, MinVarX=1, MaxVarX=NULL,
                        MinVarY=1, MaxVarY=NULL, narrow_stepX=NULL, narrow_stepY=NULL){
  if(!sparsity%in%c('X', 'Y', 'XY')){
    stop('sparsity should be one of "X", "Y" and "XY"')
  }
  
  SelectedX<-NULL
  SelectedY<-NULL
  
  MaxVarX <- ifelse(is.null(MaxVarX), yes=ncol(X), no=MaxVarX)
  MaxVarY <- ifelse(is.null(MaxVarY), yes=ncol(Y), no=MaxVarY)
  
  MinVarX=ifelse(sparsity%in%c('X', 'XY'), yes=MinVarX, no=MaxVarX)
  MinVarY=ifelse(sparsity%in%c('Y', 'XY'), yes=MaxVarX, no=MaxVarY)
  
  Summary <- data.frame(NComp=NULL, NVarX=NULL, NVarY=NULL, MSEP=NULL, Q2=NULL)
  full_summary <- data.frame()
  
  for(comp in c(1:nComp)){
    
    print(comp)
    
    if (comp>1){
      if (sparsity%in%c('X', 'XY')){
        previousX=select_num_vars(Summary, NVarX)
        SelectedX=c(SelectedX, previousX)
      }
      if (sparsity%in%c('Y', 'XY')){
        previousY=select_num_vars(Summary, NVarY)
        SelectedY=c(SelectedY, previousY)
      }
    }

    # Run the general search
    Summary <- ComponentPLS(X, Y, comp, nrepeat, stepX, stepY,
                            MinVarX, MinVarY, MaxVarX, MaxVarY,
                            SelectedX, SelectedY)
    Summary$search <- "general"
    full_summary <- rbind(full_summary, Summary)

    # Run an additional narrow search on a smaller step size centered on the
    # results from the above summary
    if (!is.null(narrow_stepX) | !is.null(narrow_stepY)){
      Summary <- NarrowPLS(Summary, X, Y, comp, nrepeat, sparsity, 
                           stepX, stepY, MinVarX, MaxVarX,
                           MinVarY, MaxVarY, narrow_stepX, narrow_stepY)
      Summary$search <- "narrow"
      full_summary <- rbind(full_summary, Summary)
    }
  }

  if (sparsity%in%c('X', 'XY')){
    previousX=select_num_vars(Summary, NVarX)
    SelectedX=c(SelectedX, previousX)
  }

  if (sparsity%in%c('Y', 'XY')){
    previousY=select_num_vars(Summary, NVarY)
    SelectedY=c(SelectedY, previousY)
  }
  
  if (sparsity=='X'){
    return(list(Summary=full_summary, keepX=SelectedX))
  }else if (sparsity=='Y'){
    return(list(Summary=full_summary, keepY=SelectedY))
  }else if (sparsity=='XY'){
    return(list(Summary=full_summary, keepX=SelectedX, keepY=SelectedY))
  }else{
    return(full_summary)
  }
}
