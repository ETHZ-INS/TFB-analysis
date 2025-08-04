library(rmarkdown)
library(data.table)

seed <- 43

metaData <- readRDS("data/01_maeAll_meta.rds")
metaData <- subset(metaData, is_matched)
metaAg <- metaData[,.(n_con=length(unique(full_id))), by=tf_name]
tfNamesBatch <- unique(subset(metaAg, n_con==2)$tf_name)

outBaseDir <- "data/predictions"
cofactors <- readRDS("data/meta_data/topInteractors.rds")

maePath <- "../data/04_maeAll_conFeat_sub.rds"
predContexts <- readRDS("data/prediction_meta.rds")
predContexts <- unique(predContexts$full_id)

for(tf in tfNamesBatch){
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
  if(file.exists(outFileFeat1)){
    if(!file.exists(outFileTrain1)){
      rmarkdown::render("./03_training_prediction/03.4_training.Rmd",
                        params=list(tfName=tf,
                                    cofactors=cofactors[[tf]], # provide them as an array
                                    maePath=maeTfPath,
                                    outDir=outDir2,
                                    profilesPath=profilesPath,
                                    predContext=predContexts,
                                    crossValidate=TRUE,
                                    seed=seed),
                        output_file=outFileTrain2)
    }

    if(!file.exists(outFilePred1)){
      rmarkdown::render("./03_training_prediction/03.5_prediction.Rmd",
                        params=list(tfName=tf,
                                    cofactors=cofactors[[tf]],
                                    maePath=maeTfPath,
                                    modelPath=modelPath,
                                    #featsToRemovePath=featsToRemovePath,
                                    outDir=outDir2,
                                    profilesPath=profilesPath,
                                    predContext=predContexts,
                                    serial=TRUE,
                                    seed=seed),
                        output_file=outFilePred2)
    }
    else{
      next
    }
  }
  else{
    rmarkdown::render("./03_training_prediction/03.3_conTf_features.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]], # provide them as an array
                                  maePath=maePath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  crossValidate=TRUE,
                                  seed=seed),
                      output_file=outFileFeat2)

    rmarkdown::render("./03_training_prediction/03.4_training.Rmd",
                      params=list(tfName=tf,
                                  cofactors=cofactors[[tf]], # provide them as an array
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
                                  #featsToRemovePath=featsToRemovePath,
                                  outDir=outDir2,
                                  profilesPath=profilesPath,
                                  predContext=predContexts,
                                  serial=TRUE,
                                  seed=seed),
                      output_file=outFilePred2)
  }
}
