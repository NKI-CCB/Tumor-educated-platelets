library(edgeR)
library(glmnet)
library(here)
library(foreach)
library(doParallel)
library(caret)
library(tidyverse)

#### Load data ----
dgeFiltered <- readRDS(here("Rds/01_dgeFiltered.Rds"))

#The original cross-correlation analysis left us with some unwanted list items.
#Clean these up
dgeFiltered <- edgeR::DGEList(counts = dgeFiltered$counts,
                              samples = dgeFiltered$samples,
                              genes = dgeFiltered$genes)
#Remove dup columns
dgeFiltered$samples$group.1 <- NULL
dgeFiltered$samples$lib.size.1 <- NULL
dgeFiltered$samples$norm.factors.1 <- NULL

#### Tuning grid

lambdas <- 10^seq(-4, 2, length = 100)
#lambdas %>% signif(2)

alphas <- seq(0,1, by=0.1)

tuning_grid <- expand.grid(alpha = alphas, lambda = lambdas)

#### Normalization ----

#Allow TMM normalization that will never use the leave-one-out sample as the reference

source(here("bin/tmm_training_norm.R")) 

#### Training function ----
enet_train <- function(
  dge = dgeFiltered,
  val_sample, #A string matching a single column name in dge$counts
  grid = tuning_grid, #An object created by expand.grid
  verboseIter = T #Whether to print progress by fold
){
  
  #Ensure that sample names don't get shuffled
  dge$samples$sample_name <- colnames(dge)
  
  #Training samples
  train_samples <- colnames(dge)[colnames(dge) != val_sample]
  
  #TMM normalization, ensure left-out sample does not serve as TMM reference
  #Followed by cpm correction and log transformation
  refsample <- getTMMref(dge, samples.for.training = train_samples)$col.index
  
  counts <- normalize.tmm(dge = dge,
                          tmm.ref.from.training = T,
                          training.samples = train_samples,
                          refCol = refsample)
  #return(counts)
  
  #Remove the LOO sample from the training counts
  train <- counts[, train_samples]
  val <- counts[, val_sample, drop = F]
  
  #return(list(train, val))
  
  #Retrieve true classes
  #Fort training set
  train_true <- dge$samples %>%
    filter(sample_name %in% colnames(train)) %>%
    select(sample_name, group)
  
  #For validation set
  val_true <- dge$samples %>%
    filter(sample_name %in% colnames(val)) %>%
    select(sample_name, group)
  
  #return(list(train_true, val_true))
  
  #Ensure that column/sample names don't get shuffled
  train_true <- train_true[order(match(train_true$sample_name,colnames(train))),]
  stopifnot(all(train_true$sample_name == colnames(train)))
  
  val_true <- val_true[order(match(val_true$sample_name, colnames(val))),]
  stopifnot(all(val_true$sample_name == colnames(val)))
  
  model <- caret::train(
    x = t(train),
    y = train_true$group,
    method = "glmnet",
    metric = "ROC",
    tuneGrid = grid,
    #tuneLength = 10,
    #trControl = trainControl("cv", number = 10),
    trControl = trainControl(
      method = "cv", number = 10,
      verboseIter = verboseIter,
      classProbs=TRUE,
      summaryFunction = twoClassSummary
    ),
  )
  
  #Return fit, training data/labels, validation sample/labels
  list(fit = model,
       train = list(data = t(train),
                    labels = train_true),
       test = list(data = t(val),
                   labels = val_true)
  )
}

#Some testing

#{
#  start <- Sys.time()
#  test = enet_train(dge = dgeFiltered, val_sample = colnames(dgeFiltered)[1])
#  end <- Sys.time()
#  end - start  
#}

#which(lambdas==test$fit$bestTune$lambda)
#tibble(lambda = lambdas, 
#       best_lambda = (lambda == test$fit$bestTune$lambda)) %>%
#  ggplot(., aes(y = lambda, x = "", color = best_lambda)) +
#  geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
  #scale_y_log10() +
#  ylim(c(-0, .25)) +
#  ggthemes::scale_color_colorblind() +
#  ggtitle("Lambda tuning")
#modelList <- foreach::foreach(a = colnames(dgeFiltered)[1:2]) %do% (
#  enet_train(dge = dgeFiltered, val_sample = a, grid = tuning_grid, verboseIter = T)
#)


#We use a leave one out cross validation strategy
{
  
  doMC::registerDoMC(cores = 8)
  set.seed(123)
  
  start <- Sys.time()
  
  modelList <- foreach::foreach(a = colnames(dgeFiltered)) %do% (
    enet_train(dge = dgeFiltered, val_sample = a, grid = tuning_grid, verboseIter = T)
  )
  
  saveRDS(modelList, file = here("Rds/05_LOOCV_enet.Rds"))
  
  end <- Sys.time()
  
  end - start
}
