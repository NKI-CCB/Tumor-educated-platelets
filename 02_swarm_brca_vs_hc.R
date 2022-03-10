library(doMC)
library(e1071)
library(foreach)
library(ggplot2)
library(methods)
library(pROC)
library(ppso)
library(reshape)
library(ROCR)
library(RUVSeq)
library(parallel)
library(dendroextras)


#### Set up ----
#Low count filtering and QC done in previous notebook
dgeFiltered <- readRDS("~/TEP/Rds/01_dgeFiltered_v2.Rds") #Compatible with R 3.4.4
training_samples <- rownames(dgeFiltered$samples[dgeFiltered$samples$Training==T,])
evaluation_samples <- rownames(dgeFiltered$samples[dgeFiltered$samples$Evaluation==T,])
validation_samples <- rownames(dgeFiltered$samples[dgeFiltered$samples$Validation==T,])

#There should be no overlap
stopifnot(length(Reduce(intersect, list(training_samples,evaluation_samples,validation_samples))) == 0)

#Sanity check
stopifnot(all(c(training_samples, evaluation_samples, validation_samples) %in% colnames(dgeFiltered)))

#Location of these scripts is HARDCODED: Must be in ./bin relative to working directory.
source('~/TEP/bin/thromboSeqTools_PreProcessing_2.R')
source('~/TEP/bin/thromboSeqTools_ANOVA.R')
source('~/TEP/bin/thromboSeqTools_PSO.R')


#### ThromboSeqPSO ----
# For this, first select the training and evaluation series. Next, perform PSO-optimization. Finally, summarize data
# and output the PSO progression plot.
#
# Args:
#   dge: DGEList with dataset count table, sample info and gene info tables.
#   percentage.for.training: Numeric value indicating the percentage of samples per group to be
#                            assigned to the training series.
#   percentage.for.evaluation: Numeric value indicating the percentage of samples per group to be
#                            assigned to the evaluation series.
#   training.samples.provided: Vector with specified column names of samples that have to be assigned to the 
#                           training series.
#   evaluation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
#                           evaluation series.
#   swarm.parameters: Vector with parameters to be optimized by swarm intelligence. By default the parameters FDR,
#                     stable genes based on lib size correlation, correlated genes, and ranked genes are included.
#                     Additional clinical parameters may be added by adding the sample info column names to this vector
#                     plus the 'variable.to.assess'-vector and a default value to the 'variable.threshold'-vector.
#   swarm.boundaries: Vector with lower and upper boundaries per variable in swarm parameters, separated
#                     by a comma. Number of input values should match with the number of swarm parameters provided.
#   k.variables: Number of (expected) variables/axes in the dataset.
#   variable.to.assess: Vector containing the column names of the sample info
#                       that are potential confounding factors. Of note, for the columns
#                       'age' or 'Age' the transcripts with a correlation coefficient below
#                       and over the provided variable.thresholds are included (see Best et al.
#                       Cancer Cell 2017).
#   variable.threshold: Vector containing manually set thresholds for the potential
#                       confounding factors in the same order as for variable.to.assess.
#   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
#                                be reached by the correlation between the counts and  any variable
#                                in order to bypass the wanted variable 'group' to be selected.
#   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
#                                be reached by the correlation between the counts and the specific 
#                                variable in order to this variable to be assigned to the RUVg axis.
#   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
#                         or PValue (FALSE) statistics.
#   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
#                                         biomarker panel.
#   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
#   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
#   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
#   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
#   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
#   n.particles: Numeric-value with number of PSO particles per iteration for classifier development.
#   n.iterations: Numeric-value with number of PSO iterations in total for classifier development.
#   figureDir: String with directory in which figures can be outputted.
#   number.cores: Vector indicating number of computational cores to be used.
#   verbose: Whether to print output or not (TRUE/FALSE).

thromboPSOfile <- "~/TEP/Rds/02_thromboPSO.Rds"


thromboPSO <- thromboSeqPSO(dge = dgeFiltered,
                            k.variables = 4,
                            variable.to.assess = c("Age","lib.size", "isolationlocation"),
                            variable.threshold = c(0.2,0.6,0.2),
                            n.particles = 100,
                            n.iterations = 10,
                            number.cores = 8,
                            training.samples.provided = training_samples,
                            evaluation.samples.provided = evaluation_samples,
                            verbose=T)


save.image(paste(Sys.Date(),"swarm_brca_vs_hc.RData",sep="_"))
saveRDS(object = thromboPSO, file = thromboPSOfile)
