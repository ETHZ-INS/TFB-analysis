library(rmarkdown)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
batch <- args[1]

# Run with docker image v1.0.0 (TFBlearner)
outBaseDir <- "data/predictions"

seed <- 43

# prioritised TFs
prioTfs <- c("KLF1", "GATA2", "MYC", "RUNX1", "RUNX2", "NR1H4", "NR1H3",
             "BANP", "ESR1", "NR3C1", "GATA1")

# load meta-data
metaData <- readRDS("data/01_maeAll_meta.rds")
metaData <- subset(metaData, is_matched)
metaAg <- metaData[,.(n_con=length(unique(full_id))), by=tf_name]

# load cofactors
cofactors <- readRDS("data/meta_data/topInteractors.rds")
cofactors <- lapply(cofactors, \(x) x[!is.na(x)])
cofactors <- cofactors[lengths(cofactors)>0]

# get version numbers
versionFiles <- list.files("data/predictions", recursive=TRUE,
                           pattern="version.tsv", full.names=TRUE)
verDts <- lapply(versionFiles, function(versionFile){
  verDt <- fread(versionFile)
  verDt$tf <- basename(dirname(versionFile))
  verDt
})
verDt <- rbindlist(verDts, fill=TRUE, use.names=TRUE)

# determine the batch of TFs to predict
verTfs <- subset(verDt, run_version=="1.0.0" & TFBlearner_version=="0.1.0")$tf
procTfs <- dirname(list.files("data/predictions", recursive=TRUE,
                              pattern="pred_282_*"))
patchedTfs <- dirname(list.files("data/predictions", recursive=TRUE,
                                 pattern="version_pred_patched*"))

tfNamesBatch <- intersect(verTfs, procTfs)
tfNamesBatch <- c(intersect(tfNamesBatch, prioTfs),
                  setdiff(tfNamesBatch, prioTfs))
batchSize <- ceiling(length(tfNamesBatch)/7)
batches <- split(tfNamesBatch, ceiling(seq_along(tfNamesBatch)/batchSize))
tfNamesBatch <- batches[[batch]]
tfNamesBatch <- setdiff(tfNamesBatch, patchedTfs)

maePath <- "../../../data/04_maeAll_conFeat_sub.rds"
predContexts <- readRDS("data/prediction_meta.rds")
predContexts <- unique(predContexts$full_id)

for(tf in tfNamesBatch){
  print(tf)
  outDir <- file.path(outBaseDir, tf)
  outDir2 <- file.path("../../..", outBaseDir, tf)

  outFilePred1 <- file.path(outDir, paste(tf, "patched_prediction.html", sep="_"))
  outFilePred2 <- file.path(outDir2, paste(tf, "patched_prediction.html", sep="_"))
  featsToRemovePath <- file.path(outDir2, "featsToRemove.rds")

  modelPath <- file.path(outDir2,
                         paste0(paste("model", tf, sep="_"), ".txt"))

  gc()
  rmarkdown::render("./03_training_prediction/patches/Run1_1.0.0/03.5_patch_prediction.Rmd",
                    params=list(tfName=tf,
                                cofactors=cofactors[[tf]],
                                maePath=maePath,
                                modelPath=modelPath,
                                featsToRemovePath=featsToRemovePath,
                                outDir=outDir2,
                                predContext=predContexts,
                                serial=TRUE,
                                seed=seed),
                    output_file=outFilePred2)
}
