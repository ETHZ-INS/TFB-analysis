library(data.table)
library(TFBlearner)

annoCol <- "context"
procTfs <- dirname(list.files("../../../data/predictions",
                              recursive=TRUE,
                              pattern="pred_286"))
for(tf in procTfs){
  tfName <- tf
  versionDt <- data.table(run_version="1.0.0",
                          TFBlearner_version="0.1.0",
                          TFBlearner_altVersion="0.0.1.1001")
  write.table(versionDt, sep="\t",
              file=file.path("../../../data/predictions", tf,
                             "version.tsv"),
              quote=FALSE, row.names=FALSE)
}
