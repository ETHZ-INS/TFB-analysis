#' getBackgroundSet
#'
#' Select a subset of `x` so that it approximates the distribution provided by
#'   `targetDist`.
#'
#' @param x A vector of values from which the set will be sampled
#' @param targetDist A vector of values of the distribution to be approximated
#' @param nQ The number of quantiles to use to approximate the distribution. If
#'   NULL (default), it will be automatically estimated.
#' @param sameSizeBg Logical; whether to return a background set of the same 
#'   length as `targetDist`.
#' @param right Logical; whether intervals should be closed to the right. Set to
#'   FALSE by default in order to explicitly have 0s and non-0s split.
#' @param replace Logical; whether to sample with replacement. Default TRUE if
#'   `sameSizeBg=FALSE`, because this increases sampling time for large
#'   datasets.
#' @param doPlot Logical; whether to do a boxplot of the two distributions as
#'   control. If the distributions are different, the most likely explanation 
#'   is that nQ is too low.
#' @author Pierre-Luc Germain
#'
#' @returns Returns the names (if `x` has names) or indices of `x` sampled to
#'   approximate the `targetDist` distribution.
getBackgroundSet <- function(x, targetDist, nQ=NULL, sameSizeBg=FALSE,
                             right=FALSE, replace=!sameSizeBg, doPlot=TRUE){
  xn <- names(x)
  if(is.null(xn)) xn <- seq_along(x)
  # first split your genes into quantile bins
  if(is.null(nQ)){
    breaks <- c(-Inf,Inf)
    # try to find a number of quantiles that leads to the targetDist getting 
    # spread into multiple quantiles, while still having a decent count per Q
    nQ <- min(ceiling(length(unique(x))/10), 50)
    minTargetQ <- pmin(length(unique(targetDist))-1L, 4L)
    while(sum((qb_in_g <- table(cut(targetDist, breaks=breaks, right=right,
                                    include.lowest=TRUE)))>0)<minTargetQ){
      breaks <- unique(quantile(x, (0:nQ)/nQ))
      nQ <- nQ*2
    }
    if(doPlot) message("Using nQ=", nQ)
  }else{
    breaks <- unique(quantile(x, (0:nQ)/nQ))
    qb_in_g <- table(cut(targetDist, breaks=breaks, include.lowest=TRUE,
                         right=right))
  }
  qb <- cut(x, breaks = breaks, include.lowest=TRUE, right=right)
  origProb <- table(qb)
  if(any(origProb<10)) warning("Too low quantile density, reduce `nQ`.")
  # select a background that have similar values to targetDist:
  if(sum(qb_in_g>0)==1)
    warning("Target distribution all falls within a single ",
    "quantile of the sampling distribution.")
  
  probv <- (qb_in_g/origProb)[as.integer(qb)]
  if(!sameSizeBg){
    # find the maximum nb of genes we can select for background
    qb_in_g <- qb_in_g[qb_in_g>0]
    n <- ceiling(min(origProb[names(qb_in_g)]/qb_in_g)*length(targetDist))
  }else{
    n <- length(targetDist)
  }
  out <- sample(xn, n, replace = replace, prob=probv/sum(probv))
  if(doPlot){
    boxplot(list(targetDist=targetDist, sampledDist=x[out]))
  }
  out
}
