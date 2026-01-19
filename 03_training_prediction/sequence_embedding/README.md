# Preparation of the sequence embeddings

First download in this folder the archetype motif matches from :
https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz
Accessed on 2022.05.24

Then run the markdowns in the broad and local subfolders.
These scripts assume that Starspace is installed and available at `/common/Starspace/` ; eventually replace this path in the markdowns.
The location of the ChIP data also needs to be adapted (see top of each markdown)

Finally, merge the two embeddings:
```
e1 <- readRDS("local/embedding10.rds")
e2 <- readRDS("broad/embedding10.rds")
e <- cbind(e1,e2)
saveRDS(e, file="embeddings.rds")
```
