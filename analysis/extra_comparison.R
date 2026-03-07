setwd("/work/tfs3/gsAI")
library(tidyverse)
library(svglite)
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
library(FSA)
library(ggsignif)
library(ggthemes)
library(multcompView)
library(rcompanion)

model_names <- c("GB", "LR", "RF", "BayesB", "BRR", "EGBLUP","GBLUP", 
                 "LASSO","RKHS")
# hex_codes <- turbo(n = 9, alpha = 0.8)
hex_codes <-c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", 
              "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF")
model_color_palette <- setNames(hex_codes, model_names)

df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics.csv")
df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))
df$extra <- 0
df <- df %>%
  filter(gen == 'all')

extra_df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics_extra.csv")
extra_df$MAF <- factor(extra_df$MAF, levels = unique(sort(extra_df$MAF)))
extra_df$extra <- 1

df <- rbind(df, extra_df)

################################################################################

MLmaf01all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.01AllExtra_fold_metrics.csv"
MLmaf01all <- read.csv(MLmaf01all) 
names(MLmaf01all) <- tolower(names(MLmaf01all))
MLmaf01all <- MLmaf01all %>%
  rename(correlation = pearsonr)
MAF01ALL <- MLmaf01all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()
MAF01ALL <- MAF01ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF01ALL$MAF <- 0.01
MAF01ALL$gen <- "all"
MAF01ALL$extra <- 1
MAF01ALL <- MAF01ALL %>%
  select(-auc_iter)

MLmaf05all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.05AllExtra_fold_metrics.csv"
MLmaf05all <- read.csv(MLmaf05all) 
names(MLmaf05all) <- tolower(names(MLmaf05all))
MLmaf05all <- MLmaf05all %>%
  rename(correlation = pearsonr)
MAF05ALL <- MLmaf05all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()
MAF05ALL <- MAF05ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF05ALL$MAF <- 0.05
MAF05ALL$gen <- "all"
MAF05ALL$extra <- 1
MAF05ALL <- MAF05ALL %>%
  select(-auc_iter)

MLmaf005all <- "/work/tfs3/gsAI/MLmodels/gebvs/MAF0.005AllExtra_fold_metrics.csv"
MLmaf005all <- read.csv(MLmaf005all) 
names(MLmaf005all) <- tolower(names(MLmaf005all))
MLmaf005all <- MLmaf005all %>%
  rename(correlation = pearsonr)
MAF005ALL <- MLmaf005all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()
MAF005ALL <- MAF005ALL %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )
MAF005ALL$MAF <- 0.005
MAF005ALL$gen <- "all"
MAF005ALL$extra <- 1
MAF005ALL <- MAF005ALL %>%
  select(-auc_iter)

df <- rbind(df, MAF01ALL, MAF05ALL, MAF005ALL)
df$extra <- factor(df$extra, levels = unique(sort(df$extra)))

################################################################################

bp <- ggplot(df[df$MAF == '0.05',], aes(x = extra, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 1, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c(levels(df$extra)[1], levels(df$extra)[2])),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Dataset", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
bp


bp2 <- ggplot(df[df$extra==1,], aes(x = MAF, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 1, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "t.test",# Placeholder comparison list
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05), # Define significance levels
    step_increase = 0.1
  ) + 
  labs(x = "MAF", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
bp2


