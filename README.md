# Tumor Educated Platelets (TEPs) in Breast Cancer

## Project Intro/Objective

Developing a blood-based classifier for early detection of breast cancer based on tumor educated platelets. We examine two types of classifiers: a particle swarm-optimized support vector machine (PSO-SVM) and an elastic net (EN).

### Approach

* Train the PSO-SVM and EN on identical data partitions with age/stage matching
* Additional nested cross-validation approaches for EN
* Differential expression analysis to further investigate cancer-associated transcripts and batch-effect-associated transcripts.

### Software

* R, esp `glmnet` and `caret` for EN training, and `tidyverse` for data wrangling 
* For PSO-SVM: [ThromboSeq](https://github.com/MyronBest/thromboSeq_source_code)
* Snakemake

## Contact persons

* Esther Lips (Corresponding author): e.lips@nki.nl
* Kat Moore (Postdoc/bioinformatician): k.moore@nki.nl
* Marte Liefaard (MD/PhD candidate): m.liefaard@nki.nl
* Lodewyk Wessels (Group leader): l.wessels@nki.nl

### Warning

PSO-SVM input code and output directory locations are HARDCODED. Do not change these unless you want to change the source code. Be warned also that repeated runs of the PSO-SVM will overwrite the previous data unless the output directory is renamed.

Project Organization
------------


    ├── Snakefile             <- Snakemake file with rules to run this project. 
    |
    ├── renv.lock             <- renv lockfile containing package versioning
    │
    ├── README.md             <- The top-level README for developers using this project.
    │
    ├── dataset               <- Directory for data storage (excluded from git repo)
    │
    ├── bin                   <- Location (hardcoded) for PSO-SVM code, and adapted TMM normalization.
    |
    ├── figureOutputFolder    <- A (hardcoded) figure folder for PSO-SVM
    │    
    ├── pso-enhanced-thromboseq1  <- PSO-SVM output
    │
    ├── Rds                   <- Interim data files associated with notebooks of the same number.
    │    
    └── src                   <- Source code for use in this project.
        │
        └── utils             <- Utility scripts, such as rmarkdown.R



--------

Note: Not all directories within the project organization hierarchy are part of the github repository. The omitted directories typically contain raw data, older analyses, or figures.
