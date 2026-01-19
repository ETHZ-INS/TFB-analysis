#' removeOverlappingRanges
#'
#' Removes elements from a GRanges that overlap (or are within a given distance
#' of) other elements higher up in the list (i.e. assumes that the ranges are
#' sorted in order of priority). The function handles overlaps between more than
#' two ranges by successively removing those that overlap higher-priority ones.
#' (This function was originally made public in the scanMiR package)
#'
#' @param x A GRanges, sorted by (decreasing) importance.
#' @param minDist Minimum distance between ranges.
#' @param ignore.strand Logical. Whether the strand of the input ranges should
#' be ignored or not.
#'
#' @return A filtered GRanges
removeOverlappingRanges <- function(x, minDist=0L, ignore.strand=TRUE){
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist,
                               ignore.strand=ignore.strand)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0){
    return(x)
  }
  i <- seq_along(x)
  toRemove <- c()
  while(length(red)>0){
    ## for each overlap set, we flag the index (relative to i) of the maximum
    ## (i.e. lowest in the list)
    top <- min(red) ## indexes of the top entry per overlap set, relative to i
    ## overlap of non-top entries to the top entries:
    o <- IRanges::overlapsAny(x[i[-top]],x[i[top]],maxgap=minDist)
    torem <- i[-top][which(o)] ## entries to remove, relative to x
    toRemove <- c(toRemove, torem) ## relative to x
    i <- setdiff(i,torem)
    ## and check again overlaps among this subset (revmap ind are relative to i)
    red <- GenomicRanges::reduce(x[i], with.revmap=TRUE, min.gapwidth=minDist,
                                 ignore.strand=ignore.strand)$revmap
    red <- red[lengths(red)>1]
  }
  if(length(toRemove)>0) x <- x[-toRemove]
  x
}
