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

###############################################################################

## load data for AI vs R comparison

ml_models_path <- "/work/tfs3/gsAI/gebvs/nocal/Oct22_MAF01_GEBVs_10foldCV.csv"
ml_metrics_path <- "/work/tfs3/gsAI/gebvs/nocal/Oct22_MAF01_fold_metrics.csv"

r_models_path <- "/work/tfs3/gsAI/4paul/gebvs/Nov08_CrossValDarpaGebv_allMAF01QC.xlsx"
r_metrics_path <- "/work/tfs3/gsAI/4paul/gebvs/Nov08_CrossValFoldMetrics_allMAF01QC.csv"

DARPAGEBV <- read.csv(ml_models_path) 
DARPAGEBV2 <- read_xlsx(r_models_path) 
df <- merge(DARPAGEBV, DARPAGEBV2)

MLMETRICS <- read.csv(ml_metrics_path) 
RMETRICS <- read.csv(r_metrics_path)
names(RMETRICS) <- tolower(names(RMETRICS))
names(MLMETRICS) <- tolower(names(MLMETRICS))
MLMETRICS <- MLMETRICS %>%
  rename(correlation = pearsonr)
RMETRICS <- RMETRICS %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

MLMETRICS <- MLMETRICS %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

allMetrics <- rbind(MLMETRICS, RMETRICS)

per_iter <- allMetrics %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

#### ANOVA assumptions check

# Normality
models <- unique(per_iter$model)
models <- models[models != "GBLUP"]
for (m in models) {
  cat("Model:", m)
  data_m <- per_iter %>% filter(model == m)
  qqPlot(data_m$corr_iter, main = paste("QQ Plot -", m))
  shapiro_res <- shapiro.test(data_m$corr_iter)
  print(shapiro_res)
}
# all pass

# Variance
levene_res <- car::leveneTest(corr_iter ~ model, data = per_iter)
print(levene_res)


############# Hypothesis Testing

models <- unique(per_iter$model)

## Visualize results
ggplot(per_iter, aes(x = model, y = corr_iter, fill = model)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
  geom_jitter(color = "black", alpha = 0.7, size = 2) +  
  labs(x = "Model", y = "Correlation") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")


#### ANOVA
anova_res <- aov(corr_iter ~ model, data = per_iter)
anova_posthoc <- TukeyHSD(anova_res)
print(anova_posthoc)
emmeans(anova_res, pairwise ~ model)

tukey_df <- as.data.frame(anova_posthoc$model)
tukey_df <- tukey_df %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
tukey_df$comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL
nonsig_tukey <- tukey_df %>%
  dplyr::filter(`p adj` > 0.05)
print(nonsig_tukey)

#### KW tests
kruskal_res <- kruskal.test(corr_iter ~ model, data = per_iter)
print(kruskal_res)
dunn_res <- dunnTest(corr_iter ~ model, data = per_iter, method="none")
print(dunn_res)

dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
nonsig <- dunn_table %>%
  dplyr::filter(P.adj > 0.05)
print(nonsig)


###############################################################################

## load data for F2 vs all comparison for AI models

f2_ai <- "/work/tfs3/gsAI/gebvs/nocal/Oct27_MAF01_F2_fold_metrics.csv"
all_ai <- "/work/tfs3/gsAI/gebvs/nocal/Oct22_MAF01_fold_metrics.csv"

f2_ai <- read.csv(f2_ai) 
all_ai <- read.csv(all_ai)
names(f2_ai) <- tolower(names(f2_ai))
names(all_ai) <- tolower(names(all_ai))

f2_ai <- f2_ai %>%
  mutate(model = dplyr::recode(model,
                               "LR" = "LR_F2",
                               "RF" = "RF_F2",
                               "GB" = "GB_F2"))

all_ai <- all_ai %>%
  mutate(model = dplyr::recode(model,
                               "LR" = "LR_all",
                               "RF" = "RF_all",
                               "GB" = "GB_all"))

allMetrics <- rbind(f2_ai, all_ai)
allMetrics <- allMetrics %>%
  rename(correlation = pearsonr)

allMetrics <- allMetrics %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

per_iter <- allMetrics %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

######### Assumption checks
models <- unique(per_iter$model)
for (m in models) {
  cat("Model:", m)
  data_m <- per_iter %>% filter(model == m)
  qqPlot(data_m$corr_iter, main = paste("QQ Plot -", m))
  shapiro_res <- shapiro.test(data_m$corr_iter)
  print(shapiro_res)
}
# all pass

# Variance
levene_res <- car::leveneTest(corr_iter ~ model, data = per_iter)
print(levene_res)
#fail

## Visualize results
ggplot(per_iter, aes(x = model, y = corr_iter, fill = model)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
  geom_jitter(color = "black", alpha = 0.7, size = 2) +  
  labs(x = "Model", y = "Correlation") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

compare_models <- function(model_prefix) {
  f2  <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_F2")]
  all <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_all")]
  tibble(
    model = model_prefix,
    t_p = t.test(f2, all, paired = TRUE)$p.value,
    wilcox_p = wilcox.test(f2, all, paired = TRUE)$p.value,
    mean_diff = mean(f2 - all)
  )
}

results <- bind_rows(
  compare_models("RF"),
  compare_models("GB"),
  compare_models("LR")
)

results

###############################################################################

## load data for MAF 0.01 vs MAF 0.05 comparison for AI models

MAF01_ai <- "/work/tfs3/gsAI/gebvs/nocal/Oct22_MAF01_fold_metrics.csv"
MAF05_ai <- "/work/tfs3/gsAI/gebvs/nocal/Oct10_allnewQC_fold_metrics.csv"

MAF01_ai <- read.csv(MAF01_ai) 
MAF05_ai <- read.csv(MAF05_ai)
names(MAF01_ai) <- tolower(names(MAF01_ai))
names(MAF05_ai) <- tolower(names(MAF05_ai))

MAF01_ai <- MAF01_ai %>%
  mutate(model = dplyr::recode(model,
                               "LR" = "LR_MAF01",
                               "RF" = "RF_MAF01",
                               "GB" = "GB_MAF01"))

MAF05_ai <- MAF05_ai %>%
  mutate(model = dplyr::recode(model,
                               "LR" = "LR_MAF05",
                               "RF" = "RF_MAF05",
                               "GB" = "GB_MAF05"))

allMetrics <- rbind(MAF01_ai, MAF05_ai)
allMetrics <- allMetrics %>%
  rename(correlation = pearsonr)

allMetrics <- allMetrics %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

per_iter <- allMetrics %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

######### Assumption checks
models <- unique(per_iter$model)
for (m in models) {
  cat("Model:", m)
  data_m <- per_iter %>% filter(model == m)
  qqPlot(data_m$corr_iter, main = paste("QQ Plot -", m))
  shapiro_res <- shapiro.test(data_m$corr_iter)
  print(shapiro_res)
}
# all pass

# Variance
levene_res <- car::leveneTest(corr_iter ~ model, data = per_iter)
print(levene_res)
#pass

## Visualize results
ggplot(per_iter, aes(x = model, y = corr_iter, fill = model)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
  geom_jitter(color = "black", alpha = 0.7, size = 2) +  
  labs(x = "Model", y = "Correlation") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

compare_models <- function(model_prefix) {
  maf01  <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_MAF01")]
  maf05 <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_MAF05")]
  tibble(
    model = model_prefix,
    t_p = t.test(maf01, maf05, paired = TRUE)$p.value,
    wilcox_p = wilcox.test(maf01, maf05, paired = TRUE)$p.value,
    mean_diff = mean(maf01 - maf05)
  )
}

results <- bind_rows(
  compare_models("RF"),
  compare_models("GB"),
  compare_models("LR")
)

results

###############################################################################

## load data for F2 vs all comparison for R models

f2_r <- "/work/tfs3/gsAI/4paul/gebvs/Nov10_CrossValFoldMetrics_F2MAF01QC.csv"
all_r <- "/work/tfs3/gsAI/4paul/gebvs/Nov08_CrossValFoldMetrics_allMAF01QC.csv"

f2_r <- read.csv(f2_r) 
all_r <- read.csv(all_r)

names(f2_r) <- tolower(names(f2_r))
names(all_r) <- tolower(names(all_r))
f2_r <- f2_r %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))
all_r <- all_r %>%
  mutate(model=dplyr::recode(model, "BL" = "LASSO"))

f2_r <- f2_r %>%
  mutate(model = dplyr::recode(model,
                               "GBLUP" = "GBLUP_F2",
                               "LASSO" = "LASSO_F2",
                               "EGBLUP" = "EGBLUP_F2",
                               "RKHS" = "RKHS_F2",
                               "BB" = "BB_F2",
                               "BRR" = "BRR_F2"))

all_r <- all_r %>%
  mutate(model = dplyr::recode(model,
                               "GBLUP" = "GBLUP_all",
                               "LASSO" = "LASSO_all",
                               "EGBLUP" = "EGBLUP_all",
                               "RKHS" = "RKHS_all",
                               "BB" = "BB_all",
                               "BRR" = "BRR_all"))

allMetrics <- rbind(f2_r, all_r)

allMetrics <- allMetrics %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()

per_iter <- allMetrics %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

######### Assumption checks
models <- unique(per_iter$model)
for (m in models) {
  cat("Model:", m)
  data_m <- per_iter %>% filter(model == m)
  qqPlot(data_m$corr_iter, main = paste("QQ Plot -", m))
  shapiro_res <- shapiro.test(data_m$corr_iter)
  print(shapiro_res)
}
# all pass except GBLUP

# Variance
levene_res <- car::leveneTest(corr_iter ~ model, data = per_iter)
print(levene_res)
#fail

## Visualize results
ggplot(per_iter, aes(x = model, y = corr_iter, fill = model)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
  geom_jitter(color = "black", alpha = 0.7, size = 2) +  
  labs(x = "Model", y = "Correlation") +
  theme_classic(base_size = 10) + 
  theme(axis.text.x = element_text(size = 6)) +
  theme(legend.position = "none")

compare_models <- function(model_prefix) {
  f2  <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_F2")]
  all <- per_iter$corr_iter[per_iter$model == paste0(model_prefix, "_all")]
  
  tibble(
    model = model_prefix,
    mean_F2 = mean(f2, na.rm = TRUE),
    mean_all = mean(all, na.rm = TRUE),
    mean_diff = mean(f2 - all, na.rm = TRUE),
    t_p = t.test(f2, all, paired = TRUE)$p.value,
    wilcox_p = wilcox.test(f2, all, paired = TRUE)$p.value
  )
}

# List of model prefixes
models <- c("BB", "BRR","EGBLUP", "RKHS", "LASSO", "GBLUP")
# remove GBLUP throws error

# Run all comparisons and bind into a single dataframe
results <- bind_rows(lapply(models, compare_models))

# View neatly formatted results
print(results)

