maxNoZero <- function(x){
  m <- max(x)
  if(m==0){
    m <- min(x)
  }
  return(m)
}

rowMaxNoZero <- function(m){
  m1 <- rowMaxs(m)
  m2 <- rowMins(m)
  m <- fifelse(m1==0, m2, m1)
  return(m)
}

getLabels <- function(combs, metaData, refCoords,
                      subset=TRUE, annotate=FALSE, margin=30){

  repStats <- lapply(combs, function(comb){
    metaSub <- subset(metaData, combination_id==comb)
    metaSub <- unique(metaSub, by="chip_path_re")
    chIPDt <- lapply(metaSub$chip_path_re,readRDS)
    names(chIPDt) <- paste(1:length(chIPDt), metaSub$chip_experiment_id,
                           metaSub$caller, sep="_")
    chIPDt <- rbindlist(chIPDt,
                        use.names=TRUE,
                        fill=TRUE, idcol="caller")

    chIPDt[,c("file", "exp", "caller"):=tstrsplit(caller, split="_")]
    if("signal" %in% colnames(chIPDt))
    {
      if("qValue" %in% colnames(chIPDt)){
        chIPDt[,qValue_ext:=fifelse(is.na(qValue), signal, qValue)]}
      else{
        chIPDt[,qValue_ext:=signal]
      }
    }
    else{
      chIPDt[,qValue_ext:=qValue]
    }

    if("is_uncertain" %in% colnames(chIPDt)){
      chIPDt[,is_uncertain:=fifelse(caller=="PICS", FALSE, is_uncertain)] # this gets handled later
    }
    else
    {
      chIPDt[,is_uncertain:=FALSE]
    }

    chIPDt[,qValue_ext:=fifelse(is_uncertain, -qValue_ext, qValue_ext)]
    chIPDt[,width:=end-start]
    chIPDt[,margin_start:=center-margin]
    chIPDt[,margin_end:=center+margin]
    chIPDt[,full_start:=fifelse(width<60, margin_start, start)]
    chIPDt[,full_end:=fifelse(width<60, margin_end, end)]

    chIPDt$start <- chIPDt$end <- NULL
    setnames(chIPDt, c("margin_start", "margin_end"), c("start", "end"))
    chIPMat <- TFBlearner::genomicRangesMapping(refRanges=refCoords,
                                     assayTable=chIPDt,
                                     byCols=c("exp", "caller"),
                                     scoreCol="qValue_ext",
                                     aggregationFun=maxNoZero)

    chIPDt$start <- chIPDt$end <- NULL
    setnames(chIPDt, c("full_start", "full_end"), c("start", "end"))
    flankMat <- TFBlearner::genomicRangesMapping(refRanges=refCoords,
                                      assayTable=chIPDt,
                                      byCols=c("exp", "caller"),
                                      scoreCol="qValue_ext",
                                      aggregationFun=maxNoZero)
    flankMat <- lapply(1:length(chIPMat), function(i){
      (flankMat[[i]]!=0)-(chIPMat[[i]]!=0)})
    flankMat <- Reduce(cbind, flankMat[-1], flankMat[[1]])
    chIPMatWide <- Reduce(cbind, chIPMat[-1], chIPMat[[1]])

    # isFlank if all experiments and callers find it as non-overlapping or flanking
    isFlank <- rowSums(flankMat>0)>0 & rowSums(chIPMatWide!=0)==0

    res <- lapply(chIPMat, function(mat){
      callers <- intersect(c("MACS2", "GEM", "SISSRS"), colnames(mat))

      if(length(setdiff(colnames(mat), "PICS"))>0){
        isPeak <- rowSums(mat!=0)>0
        isPeakOther <- rowSums(mat[,c(callers), drop=FALSE]!=0)>0

        signalMat <- rowMaxNoZero(mat[,callers, drop=FALSE])
        if("PICS" %in% colnames(mat)){
          isPeakPICS <- mat[,"PICS", drop=TRUE]>0
          isUncertainPICS <- (!isPeakOther & isPeakPICS) #|  (isPeak & rowSums(mat[,callers, drop=FALSE]<0)==ncol(mat[,callers, drop=FALSE]))
          corsPICS <- lapply(callers, function(caller){
            cor(mat[isPeakPICS,"PICS", drop=TRUE],
                mat[isPeakPICS,caller,drop=TRUE], method="spearman")})
          corsPICS <- unlist(corsPICS)
          names(corsPICS) <- callers

          qValPICS <- scales::rescale(mat[,"PICS", drop=TRUE],
                                      to=c(-log10(0.05), -log10(1e-4)))
          qValPICS[mat[,"PICS", drop=TRUE]==0]  <- 0
          signalMatAll <- rowMaxNoZero(cbind(mat[,callers, drop=FALSE],
                                             Matrix(qValPICS, ncol=1)))
        }
        else{
          isUncertainPICS <- rep(FALSE, length(isPeak))
          corsPICS <- rep(NA, length(callers))
          names(corsPICS) <- callers
          signalMatAll <- signalMat
        }

        isUncertain <- signalMat<=0
        picsOnly <- FALSE
      }
      else{
        isPeak <- rowSums(mat[,"PICS", drop=FALSE]>0)>0
        qVal <- scales::rescale(mat[,"PICS", drop=TRUE],
                                to=c(-log10(0.05), -log10(1e-4)))
        qVal[mat[,"PICS", drop=TRUE]==0]  <- 0
        isUncertainPICS <- rep(TRUE, length(isPeak))
        signalMat <- qVal
        signalMatAll <- signalMat
        isUncertain <- signalMat<=0
        corsPICS <- rep(NA, length(callers))
        picsOnly <- TRUE
      }
      list(isPeakMat=Matrix(isPeak, ncol=1),
           isUncertain=Matrix(isUncertain, ncol=1),
           isUncertainPICS=Matrix(isUncertainPICS, ncol=1),
           signalMat=Matrix(signalMat, ncol=1),
           signalMatAll=Matrix(signalMatAll, ncol=1),
           corsPICS=corsPICS,
           picsOnly=picsOnly)
    })

    corsPICS <- unlist(lapply(res, function(e) e$corsPICS))
    picsOnly <- unlist(lapply(res, function(e) e$picsOnly))

    chIPBinMat <- lapply(res, function(e) e$isPeakMat)
    chIPBinMat <- Reduce(cbind, chIPBinMat[-1], chIPBinMat[[1]])

    uncertainMat <- lapply(res, function(e) e$isUncertain)
    uncertainMat <- Reduce(cbind, uncertainMat[-1], uncertainMat[[1]])

    onlyPICSMat <- lapply(res, function(e) e$isUncertainPICS)
    onlyPICSMat <- Reduce(cbind, onlyPICSMat[-1], onlyPICSMat[[1]])

    signalMat <- lapply(res, function(e) e$signalMat)
    signalMat <- Reduce(cbind, signalMat[-1], signalMat[[1]])

    signalMatAll <- lapply(res, function(e) e$signalMatAll)
    signalMatAll <- Reduce(cbind, signalMatAll[-1], signalMatAll[[1]])

    isPeak <- which(rowSums(chIPBinMat, na.rm=TRUE)>0)
    # check if peaks are not found or found only uncertain in all experiments
    certainMat <- signalMatAll>0
    isUncertain <- rowSums(uncertainMat[isPeak,,drop=FALSE])>0 &
                   rowSums(certainMat[isPeak,,drop=FALSE])==0

    # check if peak is found to be replicated (irrespective of caller)
    isRep <- rowSums(chIPBinMat[isPeak,,drop=FALSE], na.rm=TRUE)>1

    # choose signal from other callers than PICS if available
    if(sum(picsOnly)==length(picsOnly))
    {
      all_pics <- TRUE
      isOnlyPICS <- rep(TRUE, length(isPeak))

      qVal <- rowMaxNoZero(signalMat[isPeak,, drop=FALSE])
    }
    else{
      all_pics <- FALSE
      # check if peaks have been only found by PICS
      certainMat <- signalMat>0
      isOnlyPICS <- rowSums(onlyPICSMat[isPeak,,drop=FALSE])>0 &
                    rowSums(certainMat[isPeak,!picsOnly,drop=FALSE])==0
      qVal <- rowMaxNoZero(signalMat[isPeak, !picsOnly, drop=FALSE])
      qVal[isOnlyPICS] <- rowMaxNoZero(signalMatAll[isPeak,,drop=FALSE])[isOnlyPICS]
    }

    # save also is uncertain because of PICS
    repDt <- data.table(qValue=abs(qVal),
                        label_rep=isRep,
                        is_uncertain=isUncertain,
                        is_only_PICS=isOnlyPICS)

    if(subset){
      repDt <- repDt[sample(1:nrow(repDt), min(1e4, nrow(repDt))),]
    }
    if(annotate){
      coordsDt <- as.data.table(refCoords[isPeak])
      repDt <- cbind(coordsDt, repDt)
      repDt$is_flank <- FALSE

      flankDt <- as.data.table(refCoords[isFlank])
      flankDt$is_flank <- TRUE
      repDt <- rbind(repDt, flankDt, fill=TRUE, use.names=TRUE)
      repDt[,is_uncertain:=fifelse(is_uncertain | is_flank | is_only_PICS,
                                   TRUE, FALSE)]
    }

    return(list(repDt=repDt, corsPICS=corsPICS))
  })
}
