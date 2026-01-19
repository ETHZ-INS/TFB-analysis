# make sure to load with:
# source("path/to/load_all.R", chdir=TRUE)
rfiles <- list.files(pattern="\\.R$")
sapply(rfiles[!grepl("load_all", rfiles)], source)
