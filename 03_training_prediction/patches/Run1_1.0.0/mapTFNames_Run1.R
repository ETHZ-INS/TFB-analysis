library(data.table)
library(TFBlearner)

annoCol <- "context"
procTfs <- dirname(list.files("../../../data/predictions",
                              recursive=TRUE,
                              pattern="pred_286"))

cofactors <- readRDS("../../../data/meta_data/topInteractors.rds")
cofactors <- lapply(cofactors, \(x) x[!is.na(x)])
mae <- readRDS("../../../data/04_maeAll_conFeat_sub.rds")
cols <- lapply(experiments(mae), function(n) {
  colnames(n)[colnames(n) %in% unique(subset(sampleMap(mae),
                                             !isTesting)$colname)]})
mae <- subsetByColumn(mae, cols)

for(tf in procTfs){
  tfName <- tf

  # get the tf-cofactors
  tfCofactors <- unique(cofactors[[tf]])
  names(tfCofactors) <- paste(TFBlearner:::TFCOFACTORMOTIFPREFIX,
                              1:length(tfCofactors),sep="_")

  # subset & get the chIP-Matrix
  whichCol <- which(mae[[TFBlearner:::CHIPEXP]][[TFBlearner:::TFNAMECOL]]!=tf)
  chIPMat <- as(assays(mae[[TFBlearner:::CHIPEXP]])[[TFBlearner:::PEAKASSAY]][,whichCol],"CsparseMatrix")
  colnames(chIPMat) <- paste(colData(mae[[TFBlearner:::CHIPEXP]])[whichCol,annoCol],
                             colData(mae[[TFBlearner:::CHIPEXP]])[whichCol,TFBlearner:::TFNAMECOL],
                             sep="_")

  # co-bind features
  tfCols <- unlist(tstrsplit(colnames(chIPMat), split="_", keep=2))
  namesSub <- names(tfCofactors)[which(tfCofactors %in% tfCols)]
  subTfCofactors <- intersect(tfCols, tfCofactors)
  names(subTfCofactors) <- namesSub
  coBindDt <- data.table(generic_name=names(subTfCofactors),
                         full_name=subTfCofactors)
  coBindDt[,feature_name:=paste("tfFeat_coBind", gsub("Motif", "", generic_name), sep="_")]
  coBindDt[,explicit_feature_name:=paste("tfFeat_coBind", full_name, sep="_")]

  # co-count features
  # tfs <- c(tf, tfCofactors)
  # tfs <- intersect(colData(mae[[TFBlearner:::MOTIFEXP]])[[TFBlearner:::MOTIFNAMECOL]],tfs)
  # namesCoCounts <- lapply(tfs, function(tf){
  #   if(tf==tfName){name <- paste(TFBlearner:::TFMOTIFPREFIX, 1:length(tf), sep="_")}
  #   else{name <- names(tfCofactors)[which(tfCofactors==tf)]}
  # })
  coCountDt <- data.table(generic_name=c("tfMotif_1", names(tfCofactors)),
                          full_name=c(tfName, tfCofactors))
  coCountDt[,feature_name:=paste("tfFeat_coMotifCount", generic_name, sep="_")]
  coCountDt[,explicit_feature_name:=paste("tfFeat_coMotifCount", full_name, sep="_")]

  # other features
  preSelDt <- fread(file.path("../../../data/predictions",
                              tf, paste0(paste("selMotifs", tf, sep="_"), ".tsv")))
  preSelMatchDt <- subset(preSelDt, feat=="motifMatch")
  assocMotifDt <- data.table(generic_name=preSelDt$motif_class,
                             full_name=preSelDt$motif_name)
  assocMotifDt[,feature_name:=paste("tfFeat_motifMatch", generic_name, sep="_")]
  assocMotifDt[,explicit_feature_name:=paste("tfFeat_motifMatch", full_name, sep="_")]
  assocMotifDt[,feature_name:=gsub("selectedMotifco", "selectedMotif.co", feature_name)]
  assocMotifDt[,feature_name:=gsub("selectedMotifex", "selectedMotif.ex", feature_name)]

  preSelActDt <- subset(preSelDt, feat=="assocationMotifActivity")
  assocActDt <- data.table(generic_name=preSelActDt$motif_class,
                           full_name=preSelActDt$motif_name)
  assocActDt[,feature_name:=paste("tfFeat_assocationMotifActivity", generic_name, sep="_")]
  assocActDt[,explicit_feature_name:=paste("tfFeat_assocationMotifActivity", full_name, sep="_")]

  actMotifDt <- data.table(generic_name=preSelActDt$motif_class,
                           full_name=preSelActDt$motif_name)
  actMotifDt[,feature_name:=paste("contextTfFeat_chromVAR.Activity", generic_name, sep="_")]
  actMotifDt[,feature_name:=fifelse(grepl("selectedMotif", feature_name),
                                    paste(feature_name, "normed", sep="_"),
                                    feature_name)]
  actMotifDt[,explicit_feature_name:=paste("contextTfFeat_chromVAR.Activity", full_name, sep="_")]


  # insert features
  insertsDt <- data.table(generic_name=rep("tfMotif_1", 4),
                          full_name=rep(tfName, 4))
  insertsDt$feature_name <- c("contextTfFeat_weightedInserts.margin_tfMotif_1_normed",
                              "contextTfFeat_weightedInserts.within_tfMotif_1_normed",
                              "contextTfFeat_inserts.margin_tfMotif_1_normed",
                              "contextTfFeat_inserts.within_tfMotif_1_normed")
  insertsDt[,explicit_feature_name:= gsub("tfMotif_1", paste("motifMatch", tfName, sep="_"), feature_name)]

  featNameDt <- rbindlist(list(coBindDt, coCountDt, assocMotifDt,
                               assocActDt, actMotifDt, insertsDt), use.names=TRUE)
  write.table(featNameDt, sep="\t",
              file=file.path("../../../data/predictions", tf,
                             "featureNameMap.tsv"),
              quote=FALSE, row.names=FALSE)
}












