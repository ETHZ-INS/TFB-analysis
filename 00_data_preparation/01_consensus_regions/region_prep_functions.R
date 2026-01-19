#' Merge peaks, re-splitting large merges using local overlap minima
#'
#' @param peaks A list of peaks (or unlisted)
#' @param softMaxSize The (merged) peak size below which re-splitting will be attempted
#' @param relTroughDepth The minimum depth of local minima, as a fraction of the maximum.
#'   E.g. with a maxima of 12 peaks, the default of 1/4 would require the minima to be
#'   below or equal to 9.
#' @param minTroughDepth The absolute minimum depth of local minima, in number of peaks
#'   below the maxima.
#' @param minTroughWidth The minimum width of the local minima.
#' @param minDistFromBoundary The minimum distance of the local minima from the peak
#'   border.
#' @param minPeakSize The minimum final peak size.
#' @param BPPARAM BiocParallel Param object for multithreading. If set, chromosomes are
#'   split into threads.
#'
#' @return A reduced `GRanges` of non-overlapping peaks.
reduceWithResplit <- function(peaks, softMaxSize=500L, relTroughDepth=1/4, 
                              minTroughDepth=2L, minTroughWidth=1L,
                              minDistFromBoundary=150L, minPeakSize=100L,
                              BPPARAM=BiocParallel::SerialParam()){
  library(BiocParallel)
  library(GenomicRanges)
  if(is.list(peaks) || is(peaks, "GRangesList"))
    peaks <- sort(unlist(GRangesList(peaks)))
  p <- reduce(peaks, with.revmap=TRUE)
  p1 <- p[width(p)<=softMaxSize | lengths(p$revmap)==1]
  p1$revmap <- NULL
  p <- p[width(p)>softMaxSize & lengths(p$revmap)>1]
  v <- Views(coverage(peaks), p)
  gaps <- bplapply(v, BPPARAM=BPPARAM, FUN=function(v){
    unlist(IRangesList(lapply(seq_along(v), FUN=function(i){
      x <- v[[i]]
      newR <- IRanges(1L,length(x))
      prevgaps <- gaps <- IRanges()
      curgaps <- 1L
      while(length(curgaps)>0 && any(width(newR)>softMaxSize)){
        curgaps <- unlist(IRangesList(lapply(seq_along(newR), FUN=function(j){
          x <- x[start(newR[j]):end(newR[j])]
          cth <- max(0,min(max(x)-ceiling(max(x)*relTroughDepth), max(x)-minTroughDepth))
          
          x <- .doSplitIR(x, cth, minDistFromBoundary=minDistFromBoundary,
                          minTroughWidth=minTroughWidth, minPeakSize=minPeakSize)
          shift(x, start(newR[j])-1L)
        })))
        gaps <- reduce(sort(c(gaps, curgaps)))
        if(identical(gaps,prevgaps)) curgaps <- NULL
        prevgaps <- gaps
        newR <- setdiff(IRanges(1L,length(x)), gaps)
      }
      shift(gaps, start(v)[i]-1L)
    })))
  })
  gaps <- as(IRangesList(gaps), "GRanges")
  sort(c(p1,setdiff(p, gaps)))
}

.doSplitIR <- function(x, cth, minDistFromBoundary=100L, minTroughWidth=10L, minPeakSize=100L){
  vs <- slice(x, upper=cth)
  ir <- as(vs, "IRanges")
  mcols(ir)$min <- min(vs)
  w <- start(ir)>=minDistFromBoundary & end(ir)<=length(x)-minDistFromBoundary
  ir <- ir[width(ir)>=minTroughWidth & 
             (!w | (length(x)-width(ir))>=minPeakSize)]
  w <- start(ir)>=minDistFromBoundary & end(ir)<=length(x)-minDistFromBoundary
  if(length(ir)==0) return(ir)
  if(length(which(w))>0){
    # there is a non-boundary window, use that
    # if many, use the deepest trough
    ir <- ir[mcols(ir)$min==min(mcols(ir)$min)]
    # if many, use the most central one
    w <- which.min(abs(length(x)/2-(start(ir)+width(ir)/2)))
    ir2 <- ir[w]
  }else{
    # boundary window
    ir <- ir[which.max(width(ir))]
    ir <- resize(ir, floor(width(ir)/2),
                 fix=ifelse(start(ir)>1,"start","end"))
    x2 <- x[start(ir):end(ir)]
    vs <- slice(x2, upper=min(x2))
    ir2 <- shift(as(vs, "IRanges"), start(ir)-1L)
    mcols(ir2)$min <- min(vs)
    w <- start(ir2)>=minDistFromBoundary & end(ir2)<=(length(x)-minDistFromBoundary)
    ir2 <- ir2[w]
    if(length(ir2)==0) return(ir2)
    ir2 <- ir2[which(mcols(ir2)$min==min(mcols(ir2)$min))]
    ir2 <- ir2[which.max(width(ir2))]
  }
  ir2 <- resize(ir2,pmin(width(ir2),minTroughWidth),fix="center")
  # last double-check
  ir2 <- ir2[start(ir2)>=minDistFromBoundary & end(ir2)<=(length(x)-minDistFromBoundary)]
  while(any(w <- (c(start(ir2),length(x))-c(0,end(ir2)[-length(ir2)]))<minPeakSize)){
    ir2 <- ir2[-which(w)[1]]
  }
  ir2
}


#' resizeToNonOverlapping
#' 
#' Tries to extend GRanges without creating overlaps. At the moment, this 
#' function can only handle regions that overlap one or two other regions, not
#' more.
#'
#' @param g The GRanges
#' @param desiredWidth The minimum size to (try to) extend to
#'
#' @return The modified GRanges
resizeToNonOverlapping <- function(g, desiredWidth=100L){
  if(length(reduce(g, min.gapwidth=0L))!=length(g)){
    stop("The original ranges (before extension) are overlapping.")
  }
  g2 <- resize(g, pmax(desiredWidth, width(g)), fix="center")
  g2 <- .removePairsOverlaps(g,g2,desiredWidth=desiredWidth)
  g2 <- .removePairsOverlaps(g,g2,desiredWidth=desiredWidth,reextend=FALSE)
  
}

# not used... expects xo and x to be respectively the original and extended 
# GRanges each of length exactly 3
.removeTrioOverlaps <- function(xo,x){
  stopifnot(length(x)==3)
  # overlap sizes left and right:
  d <- width(disjoin(x, with.revmap=TRUE)[c(2,4)])
  # which ranges were increased and by how much:
  wIncrease <- width(x)-width(xo)
  # if the left ranges was increased, first try to move it left by the overlap
  # or by half of it's increase in size, whichever is the smallest
  if(wIncrease[1]>0){
    do.shift <- min(d[1],round(wIncrease[1]/2))
    x[1] <- IRanges::shift(x[1], shift=-do.shift)
    d[1] <- d[1]-do.shift
  }
  # if an overlaps remains, crop right until max the initial size
  if(d[1]>0){
    end(x)[1] <- max(end(xo)[1], end(x)[1]-d[1])
  }
  # same thing for the right side:
  if(wIncrease[3]>0){
    do.shift <- min(d[2],round(wIncrease[3]/2))
    x[3] <- IRanges::shift(x[3], shift=do.shift)
    d[2] <- d[2]-do.shift
  }
  if(d[2]>0){
    start(x)[3] <- min(start(xo)[3], start(x)[3]+d[2])
  }
  # trim the central object
  if(sum(d)>0){
    start(x)[2] <- start(x)[2]+d[1]
    end(x)[2] <- end(x)[2]-d[2]
  }
  x
}

.removePairsOverlaps <- function(g, g2, desiredWidth=100L, reextend=FALSE){
  o <- findOverlaps(g2,g2,minoverlap=1L)
  o <- o[from(o)<to(o),]
  oc <- table(c(from(o),to(o)))
  oc2 <- oc[which(oc==1L)]
  o2 <- as.integer(names(oc2))
  o2 <- o[which(from(o) %in% o2 | to(o) %in% o2)]
  # for each pair, number of overlapping nucleotides:
  ov <- width(pintersect(g2[from(o2)],g2[to(o2)]))
  # constraints on the extent to which each can be trimmed, based on original coords
  minEnd <- pmax(end(g)[from(o2)], end(g2)[from(o2)]-ov)
  maxStart <- pmin(start(g)[to(o2)], start(g2)[to(o2)]+ov)
  # this is the range in which we'd respect original limits:
  buffer <- IRanges(minEnd, maxStart)
  # we use the middle point
  buffer <- resize(buffer, width=2L, fix="center")
  end(g2)[from(o2)] <- start(buffer)
  start(g2)[to(o2)] <- end(buffer)
  # check if we can, we extend the regions left and right:
  if(reextend){
    g2[from(o2)] <- setdiff(resize(g2[from(o2)], width=desiredWidth, fix="end"), g2[from(o2)-1L])
    maxI <- ifelse(to(o2)==length(g2), 0L, to(o2)+1L)
    g2[to(o2)] <- setdiff(resize(g2[to(o2)], width=desiredWidth, fix="start"), g2[maxI])
  }
  g2
}

#' resplitRegions
#'
#' Uses areas devoid of motifs/subregions (`mo`) to break large `regions` into 
#' smaller components.
#'
#' @param regions GRanges of Regions 
#' @param mo GRanges of sub-regions (e.g. motif instances) to find gaps where
#'   to break `regions`
#' @param minSize The minimum size of regions
#' @param gap_width_priority The gap width above which the size of the gap 
#'   matters less than its positioning
#' @param min_gap_width The minimum gap width for a a break
#' @param verbose 
resplitRegions <- function(regions, mo, minSize=100L, gap_width_priority=15L,
                           min_gap_width=6L, verbose=TRUE){
  mo4 <- disjoin(sort(c(regions,mo)), with.revmap=TRUE)
  mo4$motifs <- lengths(mo4$revmap)-1L
  mo4$revmap <- NULL
  mgaps <- reduce(mo4[which(mo4$motifs==0L)])
  mgaps <- mgaps[width(mgaps)>5L]
  mgaps$gap_width <- width(mgaps)
  mgaps <- resize(mgaps, width=1L, fix="center")
  
  .oneSplit <- function(regions){
    o <- findOverlaps(regions, mgaps)
    x <- mgaps[to(o)]
    x$gl <- from(o)
    x$left <- start(x)-start(regions)[x$gl]
    x$right <- end(regions)[x$gl]-end(x)
    x <- x[which(x$left>minSize & x$right>minSize & x$gap_width>=min_gap_width)]
    x <- x[order(x$gap_width<gap_width_priority, abs(x$left-x$right))]
    # select breaks
    xa <- x[which(x$gap_width>=gap_width_priority)]
    xb <- x[which(!x$gl %in% unique(xa$gl))]
    xa <- xa[!duplicated(xa$gl)]
    xb <- xb[order(-xb$gap_width)]
    xb <- xb[!duplicated(xb$gl)]
    if(verbose) message("Found ", length(xa)+length(xb), " split points")
    unlist(subtract(regions, c(xa,xb)))
  }
  
  prevLen <- 0L
  while(length(regions)!=prevLen){
    prevLen <- length(regions)
    regions <- .oneSplit(regions)
  }
  regions
}
