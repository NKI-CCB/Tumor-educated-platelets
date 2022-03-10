# Create a dgeList object that contains both the samples used to train the classifier, 
# the internal validation set that was part of the same dataset,
# and the now-unblinded validation samples

library(here)
library(edgeR)
library(tidyverse)

## Validation data ----

#Load the validation samples
dgeVal <- readRDS(here("Rds/06_dgeVal.Rds"))

head(dgeVal$samples)

#The case-control status of the blind validation set:
  
dict_df <- read_csv(here("07_validated_sample_dictionary.csv"))
head(dict_df)

## Pellet quality ----

#In the results sheet from the VUMC, there is some additional information regarding 
#ery contamination and platelet pellet quality.
#Let's add this to the dictionary.

dict_df <- readxl::read_excel(here("dataset/third_party_val",
                        "Classify BRCA NKI july 2020 plus ery data 24-7-2020_VU.xlsx"),
                   sheet = "RNAs detected") %>%
  select(sample, ery_contam = `Ery contamination`, pellet_qual = `Platelet pellet`) %>%
  left_join(dict_df, ., by = "sample")

#The rest are NA
table(dict_df$ery_contam)
table(dict_df$pellet_qual)

## Case-control status ----

#Add the case-control status to the (no longer) blind validation set:
dgeVal$samples <- dgeVal$samples %>%
  rownames_to_column("sample") %>%
  left_join(., dict_df, by = "sample") %>%
  mutate(
    group = factor(
      if_else(Status_description == "Case", "breastCancer", "healthyControl"),
      levels = c("healthyControl", "breastCancer")),
    stage_detailed = Stage,
    stage = as.integer(substr(stage_detailed, 1, 1))
    )

#Tibbles removes row names
rownames(dgeVal$samples) <- dgeVal$samples$sample

#Sanity checks
stopifnot(all(rownames(dgeVal$samples) == colnames(dgeVal)))
stopifnot(all((dgeVal$samples$group == "breastCancer") == (dgeVal$samples$Status_description == "Case")))
stopifnot(all((dgeVal$samples$group == "healthyControl") == (dgeVal$samples$Status_description == "Control")))

#Add a Label column for later
dgeVal$samples$Label <- "Blind Validation"

head(dgeVal$samples)

## Stage ----

#Recode the stage column so it matches the same factor levels we had before:
dgeVal$samples$stage <- dgeVal$samples$stage %>%
  replace_na("healthyControl") %>%
  dplyr::recode_factor(
    "healthyControl" = "healthyControl",
    "1" = "I", "2" = "II" ,
    "3" = "III", "4" = "IV"
    )
levels(dgeVal$samples$stage)

#Save the validation data as its own DGEList including metadata
saveRDS(dgeVal, here("Rds/07b_dgeVal.Rds"))

## Training data ----

#The count matrix used to train both models:
dgeFiltered <- readRDS(here("Rds/01_dgeFiltered.Rds")) 

stopifnot(all(levels(dgeFiltered$samples$stage) == levels(dgeVal$samples$stage)))

## Combine val set with training set ----

#Create a new DGEList that contains both the samples used to train the 
#models and the blind validation set:

dgeAll <- edgeR::DGEList(
  counts = cbind(
    dgeFiltered$counts, #Samples used to train the model
    dgeVal$counts #New samples
    ),
  samples = bind_rows(
    #Training samples
    #Keep original data partition for reference
    #Change label to "training"
    dgeFiltered$samples %>%
      rownames_to_column("sample") %>%
      mutate(Original_Label = Label) %>%
      mutate(Dataset = "Original",
             Label = "Training") %>%
      select(group, sample, stage, Dataset, Label, Age,
             Original_Label, isolationlocation),
    #Blind validation samples
    dgeVal$samples %>%
      mutate(Dataset = "blindVal", isolationlocation = "blindNKI",
             Label = "Validation") %>%
      select(group, sample, stage, Dataset, Label, Age,
             isolationlocation, ery_contam, pellet_qual)
  ),
  genes = dgeFiltered$genes
)

#head(dgeAll$samples)

#sanity checks
stopifnot(all(colnames(dgeAll) == rownames(dgeAll$samples)))
stopifnot(all(colnames(dgeAll) == dgeAll$samples$sample))
stopifnot(all(rowSums(is.na(dgeAll$counts)) == 0))


#Remove the extra group columns
dgeAll$samples$group.1 <- NULL
head(dgeAll$samples)

## Recode hospital ---

#Some hospitals barely contributed any samples
dgeAll$samples$isolationlocation %>% table()

#Put low count hospitals into an "other" category
dgeAll$samples$hosp <- ifelse(dgeAll$samples$isolationlocation %in% c("AMC", "UMCU", "VIENNA"),
                              "other", dgeAll$samples$isolationlocation)

dgeAll$samples$hosp <- factor(dgeAll$samples$hosp,
                              levels = c("other", "VUMC", "MGH", "NKI", "blindNKI"))
table(dgeAll$samples$isolationlocation, dgeAll$samples$hosp)

## Write DGEList ----
saveRDS(object = dgeAll, file = here("Rds/07b_dgeAll.Rds"))
