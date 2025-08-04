#' TFBvarImp
#'
#' @param fm.h5 The feature matrix, either in-memory or as HDF5Array
#' @param preds The table of predictions, with the different models as columns,
#'   as well as the "label_bin" column.
#' @param mods A named list of lightgbm models.
#' @param groups An optional (named) vector of feature assignment to groups. If 
#'   omitted, groups will be established by using the part of the feature name
#'   before the first underscore "_".
#' @param nthreads The number of threads to use.
#' @param subsampNeg The approximate number of negatives to sample. If set to 
#'   NULL, the entire feature matrix will be used.
#' @param seed The random seed for the subsampling.
#'
#' @returns A named list, including tables of feature-level and group-level 
#'   importance.
#'   
#' @importFrom lightgbm lgb.importance
#' @importFrom dplyr bind_row
#' @author Pierre-Luc Germain
#' @export
TFBvarImp <- function(fm.h5, preds, mods, groups=NULL,
                      nthreads=10, subsampNeg=1e5, seed=123){
  stopifnot(all(sapply(mods, is, class2="lgb.Booster")))
  if(!is.null(groups)) stopifnot(is.character(groups) || is.factor(groups))
  labels <- preds[,TFBlearner:::BINLABELNAME,with=FALSE]
  
  # subsetting items
  if(!is.null(subsampNeg) && subsampNeg>1000){
    p <- rowMeans(preds[,paste(TFBlearner:::PREDPREFIX, names(mods), sep="_"),with=FALSE])
    th <- .findThres(p, labels)
    w <- which(labels>0L | p>=th)
    i <- sample.int(length(labels), min(length(labels),subsampNeg+length(w)),
                    replace=FALSE)
    i <- sort(union(i,w))
    fm <- as.matrix(fm.h5[i,])
  }else{
    fm <- as.matrix(fm.h5)
  }
  # if necessary, exclude the non-feature columns:
  fm <- fm[,setdiff(colnames(fm), c(TFBlearner:::LABELCOLNAME,
                                    TFBlearner:::WIDTHFEATNAME))]
  
  # if necessary, exclude columns with only NAs:
  #rmCol <- which(colSums(is.na(fm))==nrow(fm))
  #fm <- fm[,setdiff(1:ncol(fm), rmCol)]
  #if(!is.null(groups)) groups <- groups[setdiff(1:length(groups), rmCol)]
  #print(length(groups))
  
  if(is.null(groups)){
    groups <- gsub("_.+","",colnames(fm))
  }else{
    if(is.null(names(groups))){
      stopifnot(length(groups)==ncol(fm))
    }else{
      groups <- groups[colnames(fm)]
      if(any(is.na(groups))) warning("Some features have an undefined group.")
    }
  }
  
  res <- lapply(setNames(names(mods), names(mods)), \(m){
    # get cheap variable importance
    fi1 <- lgb.importance(mods[[m]], percentage=TRUE)
    fi1 <- as.data.frame(fi1[,-1], row.names=fi1$Feature)
    
    # get shapley values
    if(m!=TFBlearner:::MODELALLNAME){
      fmSub <- fm[,setdiff(colnames(fm), 
                           paste(TFBlearner:::TFFEAT,
                                 TFBlearner:::CSCOREFEATNAME,sep="_"))]
      subGroups <- groups[which(groups!=TFBlearner:::CSCOREFEATNAME)]
    }
    else{
      fmSub <- fm
      subGroups <- groups
    }
    
    shap <- predict(mods[[m]], newdata=fmSub, type="contrib",
                    params=list(num_threads=nthreads))
    base <- shap[1,ncol(shap)]
    shap <- shap[,-ncol(shap)]
    wTP <- which(labels[i]>0L)
    
    fi2 <- data.frame(row.names=colnames(fmSub),
                      group=subGroups,
                      meanAbsShap=colMeans(abs(shap)),
                      meanPosShap=colMeans(shap[wTP,]),
                      meanNegShap=colMeans(shap[-wTP,]))
    
    m <- merge(fi1, fi2, by="row.names", all=TRUE)
    colnames(m)[1] <- "Feature"
    
    # summing SHAPs across groups
    fi3 <- aggregate(data.frame(sumAbsShap=fi2$meanAbsShap),
                     by=data.frame(Group=subGroups),
                     FUN=sum)
    list(features=m, groups=fi3, baseSHAP=base)
  })
  
  list(features=rbindlist(lapply(res, \(x) as.data.table(x[[1]],  keep.rownames=TRUE)), 
                          idcol="submodel"),
       groups=rbindlist(lapply(res, \(x) as.data.table(x[[2]], keep.rownames=TRUE)), 
                        idcol="submodel"),
       baseSHAPs=lapply(res, \(x) x[[3]]))
}

.findThres <- function(p, label){
  o <- order(-p)
  d <- data.frame(p=p[o],
                  prec=1-cumsum(1L-label[o])/seq_along(o),
                  recall=cumsum(label[o])/sum(label))
  d$dist <- sqrt(rowSums((1-as.matrix(d[,-1]))^2))
  d$p[which.min(d$dist)]
}