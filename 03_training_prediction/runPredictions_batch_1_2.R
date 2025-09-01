library(rmarkdown)
library(data.table)

seed <- 43

metaData <- readRDS("data/01_maeAll_meta.rds")
metaData <- subset(metaData, is_matched)
metaAg <- metaData[,.(n_con=length(unique(full_id))), by=tf_name]
tfNamesBatch <- unique(subset(metaAg, n_con==1)$tf_name)
tfNamesBatch <- tfNamesBatch[(floor(length(tfNamesBatch)/2)+1):length(tfNamesBatch)]

outBaseDir <- "data/predictions"
cofactors <- readRDS("data/meta_data/topInteractors.rds")
cofactors <- lapply(cofactors, \(x) x[!is.na(x) & x!="ZNF37A"])
cofactors <- cofactors[lengths(cofactors)>0]
tfNamesBatch <- intersect(tfNamesBatch, names(cofactors))

maePath <- "../data/04_maeAll_conFeat_sub.rds"
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
  maeTfPath <- file.path(outDir2, "05_mae_conTfFeat_sub.rds")
  featsToRemovePath <- file.path(outDir2, "featsToRemove.rds")

  # get the profiles path
  profilesPath <- file.path("../data/insertions", paste(tf, "median.rds", sep="_"))
  if((file.exists(outFileFeat1) &&
      !any(grepl("Error in", readLines(outFileFeat1)))) &&
     (file.exists(outFileTrain1) && 
      !any(grepl("Error in", readLines(outFileTrain1)))) &&
     (file.exists(outFilePred1) && 
      !any(grepl("Error in", readLines(outFilePred1))))){
    next
  }
  else{
    rmarkdown::render("./03_training_prediction/03.3_conTf_features.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]],
                                  maePath=maePath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  crossValidate=TRUE,
                                  seed=seed),
                      output_file=outFileFeat2)

    rmarkdown::render("./03_training_prediction/03.4_training.Rmd",
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
    rmarkdown::render("./03_training_prediction/03.5_prediction.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]],
                                  maePath=maeTfPath,
                                  modelPath=modelPath,
                                  featsToRemovePath=featsToRemovePath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  serial=FALSE,
                                  seed=seed),
                      output_file=outFilePred2)
  }
}
