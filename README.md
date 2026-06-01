# TFBlearner Analysis

This repository contains scripts used for generating TF-Binding predictions, figures, and analyses using the [TFBlearner package](https://github.com/ETHZ-INS/TFBlearner).

The benchmarking scripts are saved in a separate repository: [TFB Prediction Benchmark](https://github.com/ETHZ-INS/TFB_Prediction_Benchmark).

# Data Availability

Data including the following can be accessed from the Zenodo record: [10.5281/zenodo.18198234](https://doi.org/10.5281/zenodo.18198234).

- Generated TF-Binding predictions
- The [MultiAssayExperiment](https://www.bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) (MAE) object used for model training and prediction
- Pretrained models
- Precomputed insertion/footprint profiles
- All raw data used

Further data is visualized under: [TFB Platform](https://www.ethz-ins.org/TFBPlatform/).

# Repository Contents

- `00_data_preparation/` - input aggregation, motif scans, obtaining TF-cofactor pairs
- `01_meta_data/` - training & testing metadata curation, QC summaries, reports
- `02_preprocessing/` - labeling, MAE construction, preprocessing outputs
- `03_training_prediction/` - model training and prediction runs
- `04_ablations/` - ablation experiments
- `06_TF_meta/` - TF annotations and metadata
- `07_DTFA/` - differential Transcription Factor analyses
- `08_crowdedness/` - crowdedness analyses
- `figures/` - scripts for generating figures of the manuscript (except the scripts for the benchmark figures which can be found on the [respective repository](https://github.com/ETHZ-INS/TFB_Prediction_Benchmark)) 

# Docker Usage

The Docker image which was used to obtain the predictions on the [TFB Platform](https://www.ethz-ins.org/TFBPlatform/) can be pulled with:

```bash
docker pull ghcr.io/ethz-ins/pred_train_tfblearner:v0.1.2
```

The version tag refers to the [TFBlearner package](https://github.com/ETHZ-INS/TFBlearner) version.

To re-generate predictions using the Zenodo data (after downloading and unpacking the record into the repository so expected data paths resolve):

```bash
docker run --rm -v "$(pwd)":"$(pwd)" -w "$(pwd)" ghcr.io/ethz-ins/pred_train_tfblearner:v0.1.2 \
	Rscript 03_training_prediction/runPredictions_batched.R 5 1
```

*Note: In the `runPredictions_batched.R` command above, the argument `5` specifies `nContexts` (the number of training cellular contexts) and `1` specifies the `batch` index of TFs to run (TFs are internally split up in batches for training).*
