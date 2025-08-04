import os.path as op
import pandas as pd

SEED = 42

rule all:
    input:
      "data/04_maeAll_conFeat_sub.rds",
      #"03.2_context_features.html"
      #"data/meta_data/profileMap.tsv",
      #"data/best_component_combination.tsv",
      #"data/04_maeAll_conFeat_sub.rds",
      #"data/meta/profileMap.tsv",
      #"03.2_context_features.html"
      
rule metadata:
    input:
        ref_coords="data/hybrid.custom.hg38.resized.resplit.rds"
    params:
      atac_dir="../data/atac",
      chIP_dir="../data/chIP/raw",
      out_dir="data/meta_data"
    output: 
      meta_matched="{out_dir}/01_fullMetaAll.rds",
      meta_all="{out_dir}/01_fullMetaMatched.rds"
    shell:
      """
        Rscript -e 'rmarkdown::render("01_meta_data/01.1_meta_data.Rmd", 
                                      "html_document",
                                       params=list(chIPFilesDir="{params.chIP_dir}",
                                                   atacFilesDir="{params.atac_dir}",
                                                   outDir="../{params.out_dir}"))'
      """
      
rule reformat_ChIP:
  input:
    meta_all="{out_dir_meta}/01_fullMetaAll.rds",
    meta_matched="{out_dir_meta}/01_fullMetaMatched.rds"
  params:
    out_dir_meta="data/meta_data",
    out_dir_chIP="data/chIP/processed"
  output:
    meta_re_chIP="{out_dir_meta}/02_fullMetaChIPAll.rds",
    meta_re_all="{out_dir_meta}/02_fullReMetaAll.rds"
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.1_reformat_ChIP.Rmd", 
                                    "html_document",
                                     params=list(metaDataPath="../{input.meta_all}",
                                                 outDirMeta="../{params.out_dir_meta}",
                                                 outDirChIP="../{params.out_dir_chIP}"))'
    """
    
rule train_test_split:
  input:
    meta_re_all="{out_dir}/02_fullReMetaAll.rds"
  params:
    out_dir="data/meta_data",
  output:
    meta_train="{out_dir}/03_trainMetaAll.rds",
    meta_test="{out_dir}/03_testMetaAll.rds"
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.2_train_test_split.Rmd",
                                    "html_document",
                                    params=list(metaDataPath="../{input.meta_re_all}",
                                                outDir="../{params.out_dir}"))'
    """
    
rule label_train:
  input:
    ref_coords="data/hybrid.custom.hg38.resized.resplit.rds",
    meta_train="{out_dir_meta}/03_trainMetaAll.rds",
    meta_re_chIP="{out_dir_meta}/02_fullMetaChIPAll.rds"
  params:
    out_dir_meta="data/meta_data",
    out_dir_chIP="data/chIP/labels"
  output:
    meta_labels="{out_dir_meta}/04_trainMetaChIPAll.rds",
    gam_labels="{out_dir_meta}/gam_rep.rds"
  shell:
     """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.3_labeling_train.Rmd",
                                    "html_document",
                                     params=list(refCoordsPath="../{input.ref_coords}",
                                                 metaTrainPath="../{input.meta_train}",
                                                 metaChIPPath="../{input.meta_re_chIP}",
                                                 outDirChIP="../{params.out_dir_chIP}",
                                                 outDirMeta="../{params.out_dir_meta}"))'
     """

rule label_test:
  input:
    ref_coords="data/hybrid.custom.hg38.resized.resplit.rds",
    meta_re_chIP="{out_dir_meta}/02_fullMetaChIPAll.rds",
    meta_test="{out_dir_meta}/03_testMetaAll.rds",
    gam_labels="{out_dir_meta}/gam_rep.rds"
  params:
    out_dir_meta="data/meta_data",
    out_dir_chIP="data/chIP/labels"
  output:
    meta_labels="{out_dir_meta}/04_testMetaChIPAll.rds",
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.4_labeling_test.Rmd",
                                    "html_document",
                                     params=list(refCoordsPath="../{input.ref_coords}",
                                                 gamPath="../{input.gam_labels}",
                                                 metaTestPath="../{input.meta_test}",
                                                 metaChIPPath="../{input.meta_re_chIP}",
                                                 outDirChIP="../{params.out_dir_chIP}",
                                                 outDirMeta="../{params.out_dir_meta}"))'
    """

rule process_motif:
  input:
    ref_coords="data/hybrid.custom.hg38.resized.resplit.rds",
    meta_re_chIP="{out_dir_meta}/02_fullMetaChIPAll.rds"
  params:
    motif_dir="data/motifs",
    out_dir_motifs="data/motifs",
    out_dir_meta="data/meta_data"
  output:
    tfs_list="{out_dir_meta}/allTfs.rds",
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.5_motif_processing.Rmd",
                                    "html_document",
                                     params=list(refCoordsPath="../{input.ref_coords}",
                                                 metaChIPPath="../{input.meta_re_chIP}",
                                                 motifDir="../{params.motif_dir}",
                                                 outDirMotifs="../{params.out_dir_motifs}",
                                                 outDirMeta="../{params.out_dir_meta}"))'
    """
    
rule get_cofactors:
  input:
    tfs_list="{out_dir}/allTfs.rds",
  params:
    out_dir="data/meta_data",
  output:
    cofactor_map="{out_dir}/cofactorMaps.tsv",
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.6_cofactors.Rmd",
                                    "html_document",
                                     params=list(tfListpath="../{input.tfs_list}",
                                                 outDir="../{params.out_dir}"))'
    """

rule construct_object:
  input:
    ref_coords="{out_dir}/hybrid.custom.hg38.resized.resplit.rds",
    meta_train="{out_dir}/meta_data/03_trainMetaAll.rds",
    meta_chIP_train="{out_dir}/meta_data/04_trainMetaChIPAll.rds",
    meta_test="{out_dir}/meta_data/03_testMetaAll.rds",
    meta_chIP_test="{out_dir}/meta_data/04_testMetaChIPAll.rds",
    tfs_list="{out_dir}/meta_data/allTfs.rds",
    cofactor_map="{out_dir}/meta_data/cofactorMaps.tsv",
  params:
    motif_dir="data/motifs",
    meta_dir="data/meta_data",
    out_dir="data"
  output:
    mae_object="{out_dir}/01_maeAll.rds",
    meta_object="{out_dir}/01_maeAll_meta.rds",
    meta_pred="{out_dir}/prediction_meta.rds",
  shell:
     """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.7_object_construction.Rmd",
                                    "html_document",
                                     params=list(refCoordsPath="../{input.ref_coords}",
                                                 metaTrainPath="../{input.meta_train}",
                                                 metaChIPTrainPath="../{input.meta_chIP_train}",
                                                 metaTestPath="../{input.meta_test}",
                                                 metaChIPTestPath="../{input.meta_chIP_test}",
                                                 tfsListPath="../{input.tfs_list}",
                                                 cofactorMapPath="../{input.cofactor_map}",
                                                 metaDir="../{params.meta_dir}",
                                                 motifDir="../{params.motif_dir}",
                                                 outDir="../{params.out_dir}"))'
    """
    
rule get_insertion_profiles:
  input:
    mae_object="data/01_maeAll.rds",
    meta_object="data/01_maeAll_meta.rds",
    tfs_list="data/meta_data/allTfs.rds",
  params:
    motif_dir="data/motifs/motifs",
    meta_dir="data/meta_data",
    out_dir="data/insertions",
    seed=SEED,
  output:
    profile_map="{meta_dir}/profileMap.tsv",
  shell:
    """
      Rscript -e 'rmarkdown::render("02_preprocessing/02.8_insertion_profiles.Rmd",
                                    "html_document",
                                     params=list(maePath="../{input.mae_object}",
                                                 metaObjectPath="../{input.meta_object}",
                                                 seed="{params.seed}",
                                                 tfsListPath="../{input.tfs_list}",
                                                 metaDir="../{params.meta_dir}",
                                                 motifDir="../{params.motif_dir}",
                                                 outDir="../{params.out_dir}"))'
    """
    
rule site_features:
  input:
    mae_object="{out_dir}/01_maeAll.rds",
    embeddings="{out_dir}/embeddings.rds",
    footprints="{out_dir}/consensus_footprints_and_collapsed_motifs_hg38.bed",
  params:
    out_dir="data",
  output:
    mae_object="{out_dir}/02_maeAll_siteFeat.rds",
  shell:
     """
      Rscript -e 'rmarkdown::render("03_training_prediction/03.1_site_features.Rmd",
                                    "html_document",
                                     params=list(maePath="../{input.mae_object}",
                                                 footprintPath="../{input.footprints}",
                                                 embeddingsPath="../{input.embeddings}",
                                                 outDir="../{params.out_dir}"))'
    """
    
rule context_features:
  input:
    mae_object="{out_dir}/02_maeAll_siteFeat.rds",
  params:
    out_dir="data",
    seed=SEED,
  output:
     mae_object="{out_dir}/03_maeAll_conFeat.rds",
  shell:
    """
      Rscript -e 'rmarkdown::render("03_training_prediction/03.2_context_features.Rmd",
                                    "html_document",
                                     params=list(maePath="../{input.mae_object}",
                                                 seed="{params.seed}",
                                                 outDir="../{params.out_dir}"))'
    """
    
rule subsetting:
  input:
    mae_object="{out_dir}/03_maeAll_conFeat.rds",
  params:
    out_dir="data",
    out_dir_atac="data/atac_sub",
    seed=SEED,
  output:
     mae_object="{out_dir}/04_maeAll_conFeat_sub.rds",
  shell:
    """
      Rscript -e 'rmarkdown::render("03_training_prediction/03.3_subsetting.Rmd",
                                    "html_document",
                                     params=list(maePath="../{input.mae_object}",
                                                 seed="{params.seed}",
                                                 outDirAtac="../{params.out_dir_atac}",
                                                 outDir="../{params.out_dir}"))'
    """

'''
rule model_component_ablations:
  input:
    mae_object="{out_dir}/04_maeAll_conFeat_sub.rds",
    profile_map="{out_dir}/meta_data/profileMap.tsv",
  params:
    out_dir="data",
    seed=SEED,
  output:
     component_ablation_performance="{out_dir}/performance_component_ablation.tsv",
     best_components="{out_dir}/best_component_combination.tsv",
     component_ablation_hyperparameters="{out_dir}/hyperparameters_component_ablation.tsv"
  shell:
    """
      Rscript -e 'rmarkdown::render("04_ablations/04.1_model_component_ablations.Rmd",
                                    "html_document",
                                     params=list(maePath="../{input.mae_object}",
                                                 profileMap="../{input.profile_map}",
                                                 seed="{params.seed}",
                                                 outDir="../{params.out_dir}"))'
    """


# add model_component_ablations
  # input & subsetted + profile
# add evt extract best performing component ablation in a rule
# add feature_ablations
  # input & subsetted + profile
  # inject best performing component ablation in feature_ablations
'''
