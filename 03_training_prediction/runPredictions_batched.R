library(rmarkdown)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
nContexts <- args[1]
batch <- args[2]

outBaseDir <- "data/predictions" # TODO: switch to original directory
seed <- 43

# load meta-data
metaData <- readRDS("data/05_maeAll_meta.rds")
metaData <- subset(metaData, is_matched)
metaAg <- metaData[,.(n_con=length(unique(full_id))), by=tf_name]

# load cofactors
cofactors <- readRDS("data/meta_data/topInteractors.rds")
cofactors <- lapply(cofactors, \(x) x[!is.na(x) & x!="ZNF37A"])
cofactors <- cofactors[lengths(cofactors)>0]

# determine batch of TFs
tfNamesBatch <- unique(subset(metaAg, n_con==nContexts)$tf_name)
tfNamesBatch <- intersect(tfNamesBatch, names(cofactors))
prioTFs <- c("GATA1", "GATA2", "KLF1", "RUNX2", "RUNX1", "MYC", 
             "NR1H4", "NR1H3", "BANP", "ESR1", "NR3C1", "CEBPB", 
             "MAZ", "ZNF143", "CTCF")
prioTFs <- intersect(prioTFs, tfNamesBatch)

procTfs1 <- dirname(list.files("data/predictions", recursive=TRUE,
                              pattern="version_pred_patched*"))
versionFiles <- list.files("data/predictions", recursive=TRUE,
                           pattern="version.tsv", full.names=TRUE)
verDts <- lapply(versionFiles, function(versionFile){
  verDt <- fread(versionFile)
  verDt$tf <- basename(dirname(versionFile))
  verDt
})
verDt <- rbindlist(verDts, fill=TRUE, use.names=TRUE)

procTfs2 <- subset(verDt, run_version=="1.0.2")$tf
prioTFs <- setdiff(prioTFs, procTfs2)
procTfs <- c(procTfs1, procTfs2)


tfNamesBatch <- setdiff(tfNamesBatch, procTfs)

tfNamesBatch <- setdiff(tfNamesBatch, prioTFs)
tfNamesBatch <- c(prioTFs, tfNamesBatch)

message(paste("this set of TFs is of size", length(tfNamesBatch)))
if(length(tfNamesBatch)>8){
  batchSize <- ceiling(length(tfNamesBatch)/4)
  batches <- split(tfNamesBatch, ceiling(seq_along(tfNamesBatch)/batchSize))
  tfNamesBatch <- batches[[batch]]
}

maePath <- "../data/05_maeAll_conFeat_sub.rds"
predContexts <- readRDS("data/prediction_meta.rds")
predContexts <- unique(predContexts$full_id)

for(tf in tfNamesBatch){
  print(tf)
  outDir <- file.path(outBaseDir, tf)
  outDir2 <- file.path("..", outBaseDir, tf)

  if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)

  outFileFeat1 <- file.path(outDir, paste(tf, "features.html", sep="_"))
  outFileFeat2 <- file.path(outDir2, paste(tf, "features.html", sep="_"))

  outFileTrain1 <- file.path(outDir, paste(tf, "training.html", sep="_"))
  outFileTrain2 <- file.path(outDir2, paste(tf, "training.html", sep="_"))

  outFilePred1 <- file.path(outDir, paste(tf, "prediction.html", sep="_"))
  outFilePred2 <- file.path(outDir2, paste(tf, "prediction.html", sep="_"))

  modelPath <- file.path(outDir2,
                         paste0(paste("model", tf, sep="_"), ".txt"))
  maeTfPath <- file.path(outDir2, "06_mae_conTfFeat_sub.rds")
  featsToRemovePath <- file.path(outDir2, "featsToRemove.rds")

  # get the profiles path
  profilesPath <- file.path("../data/insertions", paste(tf, "median.rds", sep="_"))
  if((file.exists(outFileFeat1) &&
      !any(grepl("Error in", readLines(outFileFeat1)))) &&
     (file.exists(outFileTrain1) && 
      !any(grepl("Error in", readLines(outFileTrain1)))) &&
     (file.exists(outFilePred1) && 
      !any(grepl("Error in", readLines(outFilePred1)))) && 
     !(tf %in% prioTFs)){
    next
  }
  else{
    rmarkdown::render("./03_training_prediction/03.5_conTf_features.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]],
                                  maePath=maePath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  crossValidate=TRUE,
                                  seed=seed),
                      output_file=outFileFeat2)

    if(nContexts>4){
      trainRmdDir <- "./03_training_prediction/03.6_training_n_fold.Rmd"
    }else{
      trainRmdDir <- "./03_training_prediction/03.6_training.Rmd"
    }
    rmarkdown::render(trainRmdDir,
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]],
                                  maePath=maeTfPath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  crossValidate=TRUE,
                                  seed=seed),
                      output_file=outFileTrain2)
    gc()
    rmarkdown::render("./03_training_prediction/03.7_prediction.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]],
                                  maePath=maeTfPath,
                                  modelPath=modelPath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  serial=TRUE,
                                  seed=seed),
                      output_file=outFilePred2)
  }
}
