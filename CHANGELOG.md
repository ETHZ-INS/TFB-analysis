## [1.0.1] - 29.08.2025

### Changes
  - Saving of explicit feature names (featureNameMap.tsv)
  - Saving of version numbers (version.tsv)

#### Data changes
  - Removed ZNF37A from cofactors (zero positive labels in the one ChIP-experiment it has been profiled)

#### Procedure changes
  - Removed contexts with zero positive labels from training
  - Increased number of subsampled regions to 2e5 for associated motif selection (TFBlearner v0.1.1)
  - Removed `Promoter_Association` from tf-features in case no promoters for requested TF in object.
  - Ensured that a TF is not amongst its cofactors specified (only the case for AHRR).

### TFBlearner version
0.1.1

### Description
  - second run generating predictions for missing TFs & TFs the first run failed. 


## [1.0.0] - 04.08.2025

### Changes

#### Data changes

#### Procedure changes

### TFBlearner version
0.1.0 (identical to version 0.0.1.1001, earlier version numbering)

### Description
  - First run generating predictions
  - Applies to all predictions generated between: 04.08.2025 & 28.08.2025
