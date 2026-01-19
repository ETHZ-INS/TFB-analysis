#' getNonRedundantMotifs
#'
#' Fetches best motif per gene from MotifDb
#'
#' @param format Output format of the motifs
#' @param species Species, either "Hsapiens" or "Mmusculus", or "both" (gives
#'   priority to Hsapiens)
#' @param priority The data sources to include, in decreasing order of priority.
#'   priority="core" can be used to retrieve only the core HOCOMOCOv11 motifs.
#'
#' @return An object of the specified format
getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
                                  species=c("Hsapiens","Mmusculus","both"),
                                  priority=NULL){
  format <- match.arg(format)
  if(is.null(priority)){
    priority <-
      c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C",
        "HOCOMOCOv10", "jaspar2018", "JASPAR_CORE", "SwissRegulon", "hPDI",
        "cisbp_1.02", "HOCOMOCOv11-secondary-A", "HOCOMOCOv11-secondary-B",
        "HOCOMOCOv11-secondary-C", "HOCOMOCOv11-secondary-D",
        "jaspar2016", "JASPAR_2014",  "jolma2013", "stamlab", "UniPROBE")
  }else if(length(priority)==1 && priority %in% c("core","hocomoco-core")){
    priority <- c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C")
  }
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb,
                           ifelse(species=="both","Hsapiens|Mmusculus",species))
  m <- mcols(motifs)
  m <- m[m$dataSource %in% priority,]
  m$dataSource <- factor(m$dataSource, priority)
  m$organism <- factor(m$organism,
                       c(species, setdiff(sort(unique(m$organism)), species)))
  if(species=="both")  m$geneSymbol <- toupper(m$geneSymbol)
  m <- m[order(m$geneSymbol, m$organism, as.integer(m$dataSource)),]
  m <- m[!duplicated(m$geneSymbol) & !is.na(m$geneSymbol),]
  motifs <- motifs[row.names(m)]
  if(format=="universal")
    return(setNames(universalmotif::convert_motifs(motifs), m$geneSymbol))
  motifs <- switch(format,
    PFMatrix=universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
    PWMatrix=universalmotif::convert_motifs(motifs, class="TFBSTools-PWMatrix"))
  motifs <- mapply(mo=motifs, n=m$geneSymbol, FUN=function(mo,n){
    mo@name <- n
    mo
  })
  names(motifs) <- m$geneSymbol
  switch(format,
         PWMatrix=do.call(TFBSTools::PWMatrixList,motifs),
         PFMatrix=do.call(TFBSTools::PFMatrixList,motifs))
}



#' get Hocomoco v12 motifs
#'
#' @param format Output format of the motifs
#' @param topPerGene Logical; whether to keep only the top motif per gene (will 
#'   also rename the motif with gene symbols)
#'
#' @return An object of the specified format
getHocomoco <- function(format=c("PFMatrix","universal","PWMatrix"),
                        topPerGene=TRUE){
  library(TFBSTools)
  format <- match.arg(format)
  mo <- readRDS("/reference/H12CORE.motifs.rds")
  if(topPerGene){
    gs <- gsub("\\..+$","",names(mo))
    w <- which(!duplicated(gs))
    mo <- mo[w]
    gs <- gs[w]
    for(i in seq_along(gs)){
      x <- mo[[i]]@altname
      mo[[i]]@altname <- mo[[i]]@name
      mo[[i]]@name <- mo[[i]]@altname
    }
    names(mo) <- gs
  }
  if(format=="universal") return(mo)
  mo <- switch(format,
               PFMatrix=universalmotif::convert_motifs(mo, class="TFBSTools-PFMatrix"),
               PWMatrix=universalmotif::convert_motifs(mo, class="TFBSTools-PWMatrix"))
  
  switch(format,
         PWMatrix=do.call(PWMatrixList,mo),
         PFMatrix=do.call(PFMatrixList,mo))
}
