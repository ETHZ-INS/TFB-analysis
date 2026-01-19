assembleContextPreds <- function(nc, predsPath="../data/predictions", verbose=TRUE){
  lf <- list.files(predsPath, pattern=paste0("pred_",nc,"_.+\\.tsv\\.gz$"),
                   full=TRUE, recursive = TRUE)
  names(lf) <- gsub("\\..+|.+_","",basename(lf))
  
  tmpf <- \(x){
    x <- fread(x)$pred_stacked
    w <- which(x>0)
    data.frame(w=w, p=as.integer(1000*x[w]))
  }
  if(verbose){
    pr <- pblapply(lf, tmpf)
    message("Assembling...")
  }else{
    pr <- lapply(lf, tmpf)
  }
  pr <- sparseMatrix(i=unlist(lapply(pr, \(x) x$w)),
                     j=rep(seq_along(pr), sapply(pr, nrow)),
                     x=(unlist(lapply(pr, \(x) x$p)))/1000)
  colnames(pr) <- names(lf)
  pr
}

# reduce prediction matrix to a given set of regions
getPredsForPeaks <- function(peaks, preds, predRegions, keepSparse=FALSE, agfn=max){
  if(inherits(peaks, "SummarizedExperiment")) peaks <- rowRanges(peaks)
  o <- findOverlaps(peaks, predRegions)
  o <- o[!duplicated(subjectHits(o))]
  tmp <- as(preds[subjectHits(o),], "TsparseMatrix")
  tmp <- data.table(row=tmp@i+1L, col=tmp@j+1L, val=tmp@x, by=queryHits(o)[tmp@i+1L])
  tmp <- tmp[, .(maxpred = agfn(val)), by = .(col, by)]
  preds2 <- sparseMatrix(i=tmp$by, j=tmp$col, x=tmp$maxpred, dims=c(length(peaks), ncol(preds)))
  if(!keepSparse || (length(preds2@x)/length(preds2))>0.3 ) preds2 <- as.matrix(preds2)
  colnames(preds2) <- colnames(preds)
  preds2
}

# reweight predictions by regressing out the median prediction per locus
reweighPreds <- function(preds2){
  med <- rowMedians(preds2)
  preds3 <- apply(preds2, 2, \(x){
    co <- max(0.1,lm.fit(as.matrix(med), x)$coefficients[[1]])
    med2 <- pmin(pmax(0,co*med),1)
    pmax(0,x-med2)
  })
  row.names(preds3) <- row.names(preds2)
  preds3
}


getPredRes <- function(preds, counts){
  rs <- log1p(rowSums(assay(counts)))
  fits <- eBayes(lmFit(gtools::logit(t(preds)), model.matrix(~log1p(rs))))
  co <- fits$coefficients[,2]
  co[is.na(co)] <- 0
  tmp <- as.matrix(rs) %*% co
  gtools::inv.logit(gtools::logit(preds)-tmp)
}

# binarize weighted preds
binweighPreds <- function(preds2, reweighFirst=FALSE){
  if(reweighFirst) preds2 <- reweighPreds(preds2)
  cs <- colSums(preds2>0.5)
  out <- rep(FALSE, nrow(preds2))
  preds <- apply(preds2,2,\(x){
    r <- rank(-x, ties.method="first")
    th <- pmin(pmax(sum(x>=0.5), 500), length(x)/2)
    r<=th
  })
  if(!reweighFirst) preds <- preds*(ncol(preds)-rowSums(preds))^2
  as(preds, "sparseMatrix")
}



# convert sparse matrix to regulon format expected by viper
mat2regulon <- function(mat, minTargets=5L){
  mat <- as(drop0(mat),"TsparseMatrix")
  regulon <- lapply(seq_len(ncol(mat)), \(j){
    w <- which(mat@j==(j-1L))
    g <- row.names(mat)[mat@i[w]+1L]
    list(tfmode=setNames(rep(1L,length(g)),g),
         likelihood=setNames(mat@x[w],g))
  })
  names(regulon) <- colnames(mat)
  regulon[sapply(regulon, \(x) length(x[[1]])>=minTargets)]
}

# convert sparse matrix to network format expected by decoupleR
mat2network <- function(mat, minTargets=5L){
  mat <- mat[,colSums(mat)>=minTargets]
  mat <- as(drop0(mat),"TsparseMatrix")
  data.frame( source=factor(mat@j+1L, seq_len(ncol(mat)), colnames(mat)),
              target=factor(mat@i+1L, seq_len(nrow(mat)), row.names(mat)),
              mor=mat@x
  )
}


#' remove colinear regulons from regulon matrix
#'
#' @param net_mat Regulon matrix, with TFs as columns
#' @param tryKeep A vector of TFs that should ideally not be removed
#' @param cor.th A correlation threshold above which to remove one of the pair
#'
#' @returns A filtered version of `net_mat`
removeColinearity <- function(net_mat, tryKeep=c(), cor.th=0.97){
  cc2 <- cc <- cor(as.matrix(net_mat))
  cc[upper.tri(cc, diag=TRUE)] <- NA_real_
  removed <- vector("character")
  while(nrow(w <- which(cc>=cor.th, arr.ind=TRUE))>0){
    cat(".")
    # preferentially remove TFs we're less interested in
    w <- row.names(w)
    toRemove <- setdiff(w, tryKeep)
    if(length(toRemove)==0){
      toRemove <- unique(unlist(lapply(w, \(i){
        colnames(cc)[which(cc2[i,colnames(cc)]>=cor.th & colnames(cc)!=i)]
      })))
    }
    removed <- c(removed, toRemove)
    keep <- setdiff(row.names(cc),toRemove)
    cc <- cc[keep,][,keep]
  }
  cat("\n")
  message(paste0(ncol(cc),"/",ncol(cc2)," regulons kept."))
  if(length(ir <- intersect(removed, tryKeep))>0)
    warning("Some of the removed regulons are among those that should have been kept:",
            paste(ir, collapse=", "))
  net_mat[,setdiff(colnames(net_mat), removed)]
}

getMetaTfs <- function(m, cut.height=0.3){
  m <- as.matrix(m[,colSums(m)>0])
  cc <- cor(m)
  cl <- cutree(hclust(as.dist(1-cc)), h=cut.height)
  return(cl)
}

#' fastMLM
#'
#' @param accmat A matrix of peak (rows) accessibility per sample/cell (columns)
#' @param annotation A matrix of motif (columns) matches per peak (rows)
#' @param useIntercept Logical; whether to use an intercept in the mode
#' @param poisson Logical; whether to use poisson regression (assumes `accmat`
#'   is a count matrix).
#' @param minMatches The minimum number of matches for a motif to be considered
#' @param BPPARAM BiocParallel param for multithreading.
#'
#' @return A matrix (of `ncol(annotation)` rows and `ncol(accmat)` columns) with
#'   the activity scores (model t-values) of each motif in each sample.
#'   
#' @details 
#' Regresses each column of `accmat` on `annotation`, and uses the coefficients'
#' t-values as activity scores.
#' 
#' @export
fastMLM <- function(accmat, annotation, useIntercept=TRUE, poisson=FALSE,
                    minMatches=5L,BPPARAM=BiocParallel::SerialParam()){
  stopifnot(is.matrix(accmat) && is.matrix(annotation))
  if(useIntercept) annotation <- cbind(rep(1L,nrow(annotation)), annotation)
  annotation <- annotation[,colSums(annotation!=0)>=minMatches]
  res <- bplapply(seq_len(ncol(accmat)), BPPARAM=BPPARAM, FUN=function(i){
    if(!isTRUE(poisson)){
      mod <- RcppArmadillo::fastLmPure(annotation, accmat[,i])
      tvals <- mod$coefficients/mod$stderr
    }else{
      if(FALSE && require("Rfast", quietly=TRUE, include.only="glm_poisson")){
        mod <- glm_poisson(a, y[,1], full=TRUE)$info
      }else{
        mod <- glm(accmat[,i]~0+annotation, family="poisson")
        mod <- coef(summary(mod))
      }
      tvals <- mod[,1]/mod[,2]
    }
    tvals
  })
  res <- matrix(unlist(res), nrow=ncol(annotation))
  row.names(res) <- colnames(annotation)
  if(useIntercept) res <- res[-1,,drop=FALSE]
  colnames(res) <- colnames(accmat)
  res
}

stagedMLM <- function(mat, annotation, BP=BiocParallel::SerialParam(), ...){
  annotation <- as.matrix(annotation[,colSums(annotation)>0])
  clusters <- getMetaTfs(annotation, ...)
  clAnno <- sapply(split(names(clusters), clusters), \(x){
    rowMeans(m[,x,drop=FALSE]) })
  res1 <- fastMLM(mat, clAnno, BPPARAM=BP)
  w <- which(rowMedians(abs(res1))>1)
  
}