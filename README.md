# TFBlearner Analysis

This repository contains scripts used for generating TF-Binding predictions, figures etc. using the [TFBlearner package](https://github.com/ETHZ-INS/TFBlearner).

The benchmarking scripts are saved in a separate repository: [Benchmark repo](link)

# Data Availability

* Generated TF-Binding predictions: 

* The [MultiAssayExperiment](https://www.bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) object based on which
all models were trained and predictions were generated is provided at:
[Zenodo-link](a1)

* Pretrained models can be obtained from: 
[Zenodo-link](a2)

* Precomputed insertion/footprint profiles are collected at: 
[Zenodo-link](a5)
       
* All raw data used can be found at: 
[Zenodo-link](a3)

# Generating Predictions

## Using a Docker Container

The docker image can be pulled with.
```console
 docker pull ghcr.io/emsonder/predTrainImage:latest
```
Alternatively the Dockerfile is also available at [./Dockerfiles/Dockerfile](Dockerfiles/Dockerfile)

Training & prediction on selected cellular-contexts can be performed as follows.      
* `cofactors`: Example data for the cofactors can be found at [link to this repo](a4).        
* `maePath`: The MultiAssayExperiment (`mae`) object used for generating the predictions from above can be found at [Zenodo link](a1).         
* `profilesPath`: Insertion profiles can optionally be provided via the `profilesPath` (otherwise they are computed from scratch).
Examples can be found at [Link to this repo](a5) and the full collection used in the manuscript can be found at [Zenodo link](a5).         


```console
docker run --cpus=4  -v /your_directory_with_data/data:/data  ghcr.io/emsonder/predTrainImage:latest \
  "rmarkdown::render('/home/03.3_training_prediction.Rmd', \ 
  params=list(tfName="ATF3", \ 
              cofactors="/data/cofactorMap.tsv", \ 
              predContext="A549", \
              maePath="/data/mae.rds", \
              profilesPath="/data/ATF3_profile.tsv", \ 
              outDir="/predictions")"
```        

         
If models should not be trained from scratch, but rather a pretrained model should be 
used (as provided here: [link](a2)), the prediction script can be run providing a path to a pretrained model via the `modelPath` argument
as follows.

```console
docker run --cpus=4  -v /your_directory_with_data/data:/data  ghcr.io/emsonder/predTrainImage:latest \
  "rmarkdown::render('/home/03.3_training_prediction.Rmd', \ 
  params=list(tfName="ATF3", \ 
              cofactors="/data/cofactorMap.tsv", \ 
              predContext="A549", \
              maePath="/data/mae.rds", \
              modelPath="/data/ATF3_model.txt", \ 
              profilesPath="/data/ATF3_profile.tsv", \ 
              outDir="/predictions")"
```


## Without container

The predictions can also be generated using the bare Rmarkdown scripts if all dependencies of TFBlearner are installed.

```console
rmarkdown::render('03_training_prediction/03.3_training_prediction.Rmd', \ 
  params=list(tfName="ATF3", \ 
              cofactors="./data/cofactorMap.tsv", \ 
              predContext="A549", \
              maePath="./data/mae.rds", \
              profilesPath="./data/ATF3_profile.tsv", \ 
              outDir="./predictions")
```       
        
        
Similarly as above when using a pretrained model predictions can be generated as follows. 

```console
rmarkdown::render('03_training_prediction/03.3_training_prediction.Rmd', \ 
  params=list(tfName="ATF3", \ 
              cofactors="./data/cofactorMap.tsv", \ 
              predContext="A549", \
              maePath="./data/mae.rds", \
              modelPath="./data/ATF3_model.txt", \ 
              profilesPath="./data/ATF3_profile.tsv", \ 
              outDir="./predictions")
```
