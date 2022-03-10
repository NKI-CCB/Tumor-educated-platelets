#We're using a .py suffix so that we get some text highlighting in RStudio
#Invoke via -s flag, example:
#snakemake -np gliph -s Snakefile.py -j 8


#Basic metadata exploration and validation
rule new_metadata:
    input:
        rmd="00_new_metadata.Rmd",
        script="src/utils/rmarkdown.R",
        dgeUpdated = "dataset/BrCa-dataset-201218_updatedClinicalInfo080821.RData",
        dgePrevious = "dataset/BrCa-dataset-201218.RData"
    output:
        html="00_new_metadata.html"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Assign samples to fixed data partitions for 1:1 comparison PSO vs enet
rule sample_labels:
    input:
        rmd="01_sample_label_allocation.Rmd",
        script="src/utils/rmarkdown.R",
        dgeIncludedSamples = "dataset/BrCa-dataset-201218_updatedClinicalInfo080821.RData",
    output:
        html = "01_sample_label_allocation.html",
        dgeFiltered = "Rds/01_dgeFiltered.Rds",
        sample_labels = "01_sample_labels.csv"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#This step is slow, usually >1 week of runtime on most servers
#We therefore do not include it as part of the Snakemake pipeline, to avoid any inadvertant reruns that delete the output data
#However, the commented-out version is included here for the sake of good documentation

#Run the particle swarm, by necessity on fixed data partitions
#The PSO results are in ./pso-enhanced-thromboseq, which is hardcoded in the ThromboSeq code.
#New runs will overwrite the results. To prevent this from accidentally occuring, rename the dir.
#rule pso:
#    input:
#        script="02_swarm_brca_vs_hc.R",
#        thrombo_anova = "bin/thromboSeqTools_ANOVA.R",
#        thrombo_preproc = "bin/thromboSeqTools_PreProcessing_2.R",
#        thrombo_PSO = "bin/thromboSeqTools_PSO.R",
#        dgeFiltered = "Rds/01_dgeFiltered.Rds"
#    output:
#        dir("pso-enhanced-thromboseq1")
#    shell:
#        "Rscript {input.script}"
#        "mv -v ./pso-enhanced-thromboseq ./pso-enhanced-thromboseq1"


#Assess PSO performance on internal validation set
rule pso_results:
    input: 
        rmd = "03_pso_results_vs_hc.Rmd",
        script="src/utils/rmarkdown.R",
        dgeFiltered = "Rds/01_dgeFiltered.Rds",
        thrombo_anova = "bin/thromboSeqTools_ANOVA.R",
        thrombo_preproc = "bin/thromboSeqTools_PreProcessing_2.R",
        thrombo_PSO = "bin/thromboSeqTools_PSO.R",
    output:
        html = "03_pso_results_vs_hc.html",
        dgeParticle = "Rds/03_dgeParticle.Rds",
        pso_preds = "03_pso_sample_predictions.csv",
        pso_summary = "03_pso_performance_summary.csv",
        detailed_pso = "03_pso_detailed_stats.csv"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Train elastic net on fixed data partitions, and assess performance
#Also compare with pso results
rule enet:
    input:
        rmd = "04_enet_vs_hc.Rmd",
        script = "src/utils/rmarkdown.R",
        tmm = "bin/tmm_training_norm.R",
        thrombo_anova = "bin/thromboSeqTools_ANOVA.R",
        thrombo_preproc = "bin/thromboSeqTools_PreProcessing_2.R",
        thrombo_PSO = "bin/thromboSeqTools_PSO.R",
        thromboPSO = "Rds/02_thromboPSO.Rds", 
        pso_summary = "03_pso_performance_summary.csv",
        detailed_pso = "03_pso_detailed_stats.csv",
        dgeFiltered = "Rds/01_dgeFiltered.Rds"
    output:
        glmlambda = "Rds/04_model_glmlambda.Rds",
        altlambda = "Rds/04_model_altlambda.Rds",
        enet_pred = "04_enet_sample_predictions.csv",
        detailed_enet = "04_enet_detailed_performance.csv",
        pso_vs_enet_summary = "04_enet_vs_psosvm_summary.csv",
        detailed_comparison = "04_enet_vs_psosvm_detailed.csv",
        html = "04_enet_vs_hc.html"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Graph ROCs for enet and pso together
rule combine_rocs:
    input:
        rmd = "04b_combine_ROCs.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        html = "04b_combine_ROCs.html"
    shell:
       "Rscript {input.script} {input.rmd} $PWD/{output.html}"
       
#Train elastic net using leave-one-out cross validation
rule enet_xval:
    input:
        script = "05_enet_xval.R",
        tmm = "bin/tmm_training_norm.R",
        dgeFiltered = "Rds/01_dgeFiltered.Rds"
    output:
        loocv_enet = "Rds/05_LOOCV_enet.Rds"
    shell:
        "echo $PWD\n"
        "Rscript {input.script}"

#Report on enet LOOCV performance
rule report_enet_xval:
    input:
        dgeFiltered = "Rds/01_dgeFiltered.Rds",
        loocv_enet = "Rds/05_LOOCV_enet.Rds", #Very large file that will take a while to load
        rmd = "05_enet_Xval.Rmd",
        script = "src/utils/rmarkdown.R"
    output:
        loocv_summary = "05_enet_LOOCV_performance.csv",
        loocv_predictions = "05_enet_LOOCV_predictions.csv",
        html = "05_enet_Xval.html"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Make predictions on blind validation set originating from NKI
rule predict_blindval:
    input:
        #dir(dataset/NKI-htseq), #loose htseq.ssv files from the blind validation set
        dgeFiltered = "Rds/01_dgeFiltered.Rds",
        best_pso="Rds/02_thromboPSO.Rds",
        best_enet = "Rds/04_model_altlambda.Rds",
        tmm = "bin/tmm_training_norm.R",
        rmd = "06_predict_blindval.Rmd",
        script = "src/utils/rmarkdown.R"
    output:
        dgeBlindVal = "dataset/blindVal.Rds", #data matrix produced by collect.read.counts()
        dgeVal = "Rds/06_dgeVal.Rds", #dgeList containing blind validation samples
        dgeRUV = "Rds/06_dgeRUV.Rds", #RUV normalized validation samples for PSO predictions
        preds_rds = "Rds/06_predictions.Rds",
        preds_excel = "06_predictions.xlsx",
        html = "06_predict_blindval.html"
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Check the performance of both classifiers on the now-unblinded validation samples
#Also cross-reference with performance reported by independent third party
rule blindval_performance:
    input:
        perf_summary = "dataset/third_party_val/Blind_tep_performance_summary.csv", #Performance reported by the third party
        detailed_enet = "04_enet_detailed_performance.csv",
        detailed_pso = "03_pso_detailed_stats.csv",
        dgeFiltered = "Rds/01_dgeFiltered.Rds",
        cc_third = "dataset/third_party_val/06_predictions_NKI.xlsx", #Sample predictions and true labels by third party
        preds_excel = "06_predictions.xlsx", #Own predictions
        sample_keys = "dataset/third_party_val/sample_keys/Overview_NKI_breast_6-8-2020.xlsx", #For blinded samples
        sample_labels = "dataset/third_party_val/sample_keys/202007_Validationset_caco_charact_anonymous.xlsx", #For blinded samples
        dgeVal = "Rds/06_dgeVal.Rds", #For matching names from NKI with those received from VUMC
        rmd = "07_blind_val_performance.Rmd",
        script = "src/utils/rmarkdown.R"
    output:
        html = "07_blind_val_performance.html",
        sample_dict = "07_validated_sample_dictionary.csv", #Unblinded cases/controls
        perf_blindval = "07_blind_val_performance.csv" #Previous and current performance of both classifiers on (un)blinded data
    shell:
        "echo $PWD\n"
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Create a DGEList that contains both the original data and the now-unblinded validation samples
rule dgeAll:
    input:
        dgeVal = "Rds/06_dgeVal.Rds",
        sample_dict = "07_validated_sample_dictionary.csv",
        pellet_qual = "dataset/third_party_val/Classify BRCA NKI july 2020 plus ery data 24-7-2020_VU.xlsx", #Pellet quality info
        dgeFiltered = "Rds/01_dgeFiltered.Rds", 
        script = "07b_dgeAll.R",
    output:
        dgeVal = "Rds/07b_dgeVal.Rds", #validation sample DGEList with metadata
        dgeAll = "Rds/07b_dgeAll.Rds" #All samples, original & blindval, with metadata
    shell:
        "Rscript {input.script}"

#Post hoc analysis: dimensionality plots (PCA, t-SNE etc, and a few heatmaps)
rule dimensionality:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        enet = "Rds/04_model_altlambda.Rds",
        pso = "Rds/02_thromboPSO.Rds",
        rmd = "08_dimensionality.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        html = "08_dimensionality.html",
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"    

#Post hoc analysis: differential expression
rule post_hoc_diffex:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        enet = "Rds/04_model_altlambda.Rds",
        pso = "Rds/02_thromboPSO.Rds",
        rmd = "09_post_hoc_diffex.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        #martFile = "Rds/martFile.Rds",
        html = "09_post_hoc_diffex.html",
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Post hoc: interhospital classifiers
#Train on NKI, predict on other centers
rule interhosp:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        enet = "Rds/04_model_altlambda.Rds",
        loocv="05_enet_LOOCV_predictions.csv",
        rmd = "10_interhosp_classifiers.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        kappaNKImodel = "Rds/10_interhosp_kappa.Rds",
        dsNKImodel = "Rds/10_interhosp_ds.Rds",
        interhosp_res = "10_interhospital_performance_summary.csv",
        html = "10_interhosp_classifiers.html"
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"        
        
#Variance partition
rule varpar:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        enet = "Rds/04_model_altlambda.Rds",
        pso = "Rds/02_thromboPSO.Rds",
        rmd = "11_variancepartition.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        html = "11_variancepartition.html"
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Feature plots
rule features:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        enet = "Rds/04_model_altlambda.Rds",
        pso = "Rds/02_thromboPSO.Rds",
        rmd = "12_features.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        html = "12_features.html"
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Exploration of the efficacy of batch correction methods
rule batchcor:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        tmmnorm = "bin/tmm_training_norm.R",
        enet_preds = "04_enet_sample_predictions.csv",
        blindval_preds = "06_predictions.xlsx",
        dgeRUV = "Rds/01_dgeQC.Rds",
        dgeBlindRUV = "Rds/06_dgeRUV.Rds",
        rmd = "13_batch_correction.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        enet_combat = "Rds/13_combat.Rds",
        enet_ruv = "Rds/13_RUV.Rds",
        enet_blindcombat = "Rds/13_blind_combat.Rds",
        enet_blindruv = "Rds/13_blind_RUV.Rds",
        batchres = "13_batch_summary.csv",
        html = "13_batch_correction.html"
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
        
#Assess potential lymphocyte & erythrocyte contamination
rule contam:
    input:
        dgeAll = "Rds/07b_dgeAll.Rds",
        erydf = "dataset/third_party_val/20201217_NKI_ValSet_VU_ID_merged_redness_complete.xlsx",
        enet_preds = "04_enet_sample_predictions.csv",
        blindval_preds = "06_predictions.xlsx",
        sampledict = "07_validated_sample_dictionary.csv",
        rmd = "14_ery_lymph_contam.Rmd",
        script = "src/utils/rmarkdown.R",
    output:
        dgePred = "Rds/14_dgePred.Rds", #As dgeAll but includes predictions in the metadata
        html = "14_ery_lymph_contam.html"
    shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Manuscript figures
rule figures:
  input:
        rmd = "15_manuscript_figures.Rmd",
        script = "src/utils/rmarkdown.R",
  output:
        html = "15_manuscript_figures.html",
        metadata_tbl = "15_metadata_table.csv"
  shell:
        "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    
