###############################################################################

## Set WD and load packages 

setwd("/work/tfs3/gsAI")
library(tidyverse)
library(readxl)
library(vegan)
library(RColorBrewer)
library(broom)
library(car)
library(rstatix)
library(FSA)
library(viridis)
library(emmeans)
library(tibble)

###############################################################################

### MAF 01 ALL
MLmaf01all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.01AllOld_fold_metrics.csv"
Rmaf01all <- "/work/tfs3/gsAI/Rmodels/gebvs/Dec06_CrossValFoldMetrics_allMAF01QC.csv"
MLmaf01all <- read.csv(MLmaf01all) 
Rmaf01all <- read.csv(Rmaf01all)
names(Rmaf01all) <- tolower(names(Rmaf01all))
names(MLmaf01all) <- tolower(names(MLmaf01all))
MLmaf01all <- MLmaf01all %>%
  rename(correlation = pearsonr)
Rmaf01all <- Rmaf01all %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf01all <- MLmaf01all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF01ALL <- rbind(MLmaf01all, Rmaf01all)

MAF01ALL <- MAF01ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF01ALL$MAF <- 0.01
MAF01ALL$gen <- "all"
MAF01ALL <- MAF01ALL %>%
  select(-auc_iter)


###############################################################################

### MAF 05 ALL
MLmaf05all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.05AllOld_fold_metrics.csv"
Rmaf05all <- "/work/tfs3/gsAI/Rmodels/gebvs/Nov30_CrossValFoldMetrics_allMAF05QC.csv"
MLmaf05all <- read.csv(MLmaf05all) 
Rmaf05all <- read.csv(Rmaf05all)
names(Rmaf05all) <- tolower(names(Rmaf05all))
names(MLmaf05all) <- tolower(names(MLmaf05all))
MLmaf05all <- MLmaf05all %>%
  rename(correlation = pearsonr)
Rmaf05all <- Rmaf05all %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf05all <- MLmaf05all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF05ALL <- rbind(MLmaf05all, Rmaf05all)

MAF05ALL <- MAF05ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF05ALL$MAF <- 0.05
MAF05ALL$gen <- "all"
MAF05ALL <- MAF05ALL %>%
  select(-auc_iter)

###############################################################################

## MAF 05 F2 
MLmaf05F2 <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.05F2Old_fold_metrics.csv"
Rmaf05F2 <- "/work/tfs3/gsAI/Rmodels/gebvs/Nov26_CrossValFoldMetrics_F2MAF05QC.csv"
MLmaf05F2 <- read.csv(MLmaf05F2) 
Rmaf05F2 <- read.csv(Rmaf05F2)
names(Rmaf05F2) <- tolower(names(Rmaf05F2))
names(MLmaf05F2) <- tolower(names(MLmaf05F2))
MLmaf05F2 <- MLmaf05F2 %>%
  rename(correlation = pearsonr)
Rmaf05F2 <- Rmaf05F2 %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf05F2 <- MLmaf05F2 %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF05F2 <- rbind(MLmaf05F2, Rmaf05F2)

MAF05F2 <- MAF05F2 %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF05F2$MAF <- 0.05
MAF05F2$gen <- "F2"
MAF05F2 <- MAF05F2 %>%
  select(-auc_iter)

###############################################################################

## MAF 01 F2

MLmaf01F2 <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.01F2Old_fold_metrics.csv"
Rmaf01F2 <- "/work/tfs3/gsAI/Rmodels/gebvs/Nov27_CrossValFoldMetrics_F2MAF01QC.csv"
MLmaf01F2 <- read.csv(MLmaf01F2) 
Rmaf01F2 <- read.csv(Rmaf01F2)
names(Rmaf01F2) <- tolower(names(Rmaf01F2))
names(MLmaf01F2) <- tolower(names(MLmaf01F2))
MLmaf01F2 <- MLmaf01F2 %>%
  rename(correlation = pearsonr)
Rmaf01F2 <- Rmaf01F2 %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf01F2 <- MLmaf01F2 %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF01F2 <- rbind(MLmaf01F2, Rmaf01F2)

MAF01F2 <- MAF01F2 %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF01F2$MAF <- 0.01
MAF01F2$gen <- "F2"
MAF01F2 <- MAF01F2 %>%
  select(-auc_iter)

###############################################################################

### TBA MAF005 All
# 
# MLmaf005all <- "X"
# Rmaf005all <- "Y"
# MLmaf005all <- read.csv(MLmaf005all) 
# Rmaf005all <- read.csv(Rmaf005all)
# names(Rmaf005all) <- tolower(names(Rmaf005all))
# names(MLmaf005all) <- tolower(names(MLmaf005all))
# MLmaf005all <- MLmaf005all %>%
#   rename(correlation = pearsonr)
# Rmaf005all <- Rmaf005all %>%
#   mutate(model=dplyr::recode(model, "BL" = "LASSO"))
# 
# MLmaf005all <- MLmaf005all %>%
#   group_by(model) %>%
#   mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
#   ungroup()
# 
# MAF005ALL <- rbind(MLmaf005all, Rmaf005all)
# 
# MAF005ALL <- MAF005ALL %>%
#   group_by(model, iteration) %>%
#   summarise(
#     auc_iter = mean(auc, na.rm = TRUE),
#     corr_iter = mean(correlation, na.rm = TRUE),
#     .groups = "drop"
#   )
# MAF005ALL$MAF <- 0.005
# MAF005ALL$gen <- "all"
# MAF005ALL <- MAF005ALL %>%
#   select(-auc_iter)

###############################################################################
### TBA MAF005 All

MLmaf005f2 <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.005F2Old_fold_metrics.csv"
Rmaf005f2 <- "/work/tfs3/gsAI/Rmodels/gebvs/Nov29_CrossValFoldMetrics_F2MAF005QC.csv"
MLmaf005f2 <- read.csv(MLmaf005f2) 
Rmaf005f2 <- read.csv(Rmaf005f2)
names(Rmaf005f2) <- tolower(names(Rmaf005f2))
names(MLmaf005f2) <- tolower(names(MLmaf005f2))
MLmaf005f2 <- MLmaf005f2 %>%
  rename(correlation = pearsonr)
Rmaf005f2 <- Rmaf005f2 %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf005f2 <- MLmaf005f2 %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF005F2 <- rbind(MLmaf005f2, Rmaf005f2)

MAF005F2 <- MAF005F2 %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF005F2$MAF <- 0.005
MAF005F2$gen <- "F2"
MAF005F2 <- MAF005F2 %>%
  select(-auc_iter)

###############################################################################

all_foldmets <- rbind(MAF005F2, MAF005ALL, MAF01F2, MAF01ALL, MAF05F2, MAF05ALL)
write.csv(all_foldmets, "/work/tfs3/gsAI/data/combined_fold_metrics.csv", row.names = FALSE)
