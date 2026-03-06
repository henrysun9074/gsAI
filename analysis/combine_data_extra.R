###############################################################################

## Set WD and load packages 

setwd("/work/tfs3/gsAI")
library(tidyverse)
library(readxl)
library(readr)
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
MLmaf01all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.01AllExtra_fold_metrics.csv"
Rmaf01all <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar05_CrossValFoldMetrics_extraMAF01QC.csv"
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
MLmaf05all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.05AllExtra_fold_metrics.csv"
Rmaf05all <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar04_CrossValFoldMetrics_extraMAF05QC.csv"
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

## MAF005 All

MLmaf005all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.005AllExtra_fold_metrics.csv"
Rmaf005all <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar03_CrossValFoldMetrics_extraMAF005QC.csv"
MLmaf005all <- read.csv(MLmaf005all)
Rmaf005all <- read.csv(Rmaf005all)
names(Rmaf005all) <- tolower(names(Rmaf005all))
names(MLmaf005all) <- tolower(names(MLmaf005all))
MLmaf005all <- MLmaf005all %>%
  rename(correlation = pearsonr)
Rmaf005all <- Rmaf005all %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLmaf005all <- MLmaf005all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

MAF005ALL <- rbind(MLmaf005all, Rmaf005all)

MAF005ALL <- MAF005ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF005ALL$MAF <- 0.005
MAF005ALL$gen <- "all"
MAF005ALL <- MAF005ALL %>%
  select(-auc_iter)

########################

all_foldmets <- rbind(MAF005ALL, MAF01ALL, MAF05ALL)
write.csv(all_foldmets, "/work/tfs3/gsAI/data/combined_fold_metrics_extra.csv", row.names = FALSE)


###############################################################################
###############################################################################
###############################################################################

## combine GEBV data

MLMAF05gebv <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.05AllExtra_GEBVs_10foldCV.csv"
RMLMAF05gebv <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar04_CrossValDarpaGebv_extraMAF05QC.xlsx"
DARPAGEBV <- read.csv(MLMAF05gebv)
DARPAGEBV2 <- read_xlsx(RMLMAF05gebv)
DARPAGEBV2 <- DARPAGEBV2 %>%
  rename_with(~ str_remove(., "_Mean$"), ends_with("_Mean"))
df <- inner_join(DARPAGEBV, DARPAGEBV2, by = "ID")
df <- df %>%
  dplyr::rename(Status = Status.x) %>%
  dplyr::select(-Status.y)


MLMAF01gebv <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.01AllExtra_GEBVs_10foldCV.csv"
RMLMAF01gebv <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar05_CrossValDarpaGebv_extraMAF01QC.xlsx"
DARPAGEBV01 <- read.csv(MLMAF01gebv)
DARPAGEBV01_2 <- read_xlsx(RMLMAF01gebv)
DARPAGEBV01_2 <- DARPAGEBV01_2 %>%
  rename_with(~ str_remove(., "_Mean$"), ends_with("_Mean"))
df2 <- inner_join(DARPAGEBV01, DARPAGEBV01_2, by = "ID")
df2 <- df2 %>%
  dplyr::rename(Status = Status.x) %>%
  dplyr::select(-Status.y)

MLMAF005gebv <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.005AllExtra_GEBVs_10foldCV.csv"
RMLMAF005gebv <- "/work/tfs3/gsAI/Rmodels/gebvs/Mar03_CrossValDarpaGebv_extraMAF005QC.xlsx"
DARPAGEBV005 <- read.csv(MLMAF005gebv)
DARPAGEBV005_2 <- read_xlsx(RMLMAF005gebv)
DARPAGEBV005_2 <- DARPAGEBV005_2 %>%
  rename_with(~ str_remove(., "_Mean$"), ends_with("_Mean"))
df3 <- inner_join(DARPAGEBV005, DARPAGEBV005_2,by = "ID")
df3 <- df3 %>%
  dplyr::rename(Status = Status.x) %>%
  dplyr::select(-Status.y)

df$MAF <- 0.05
df2$MAF <- 0.01
df3$MAF <- 0.005
combined_df <- rbind(df, df2, df3)
write.csv(combined_df, file = "/work/tfs3/gsAI/data/combined_gebvs_extra.csv", row.names = FALSE)

gebvdf<-read.csv("/work/tfs3/gsAI/data/combined_gebvs_extra.csv")