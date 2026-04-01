###############################################################################

## Set WD and load packages + data

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
library(cowplot)
library(ggsignif)
library(ggthemes)
library(multcompView)
library(rcompanion)

model_names <- c("GB", "LR", "RF", "BayesB", "BRR", "EGBLUP","GBLUP", 
                 "LASSO","RKHS")
hex_codes <-c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", 
              "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF")
model_color_palette <- setNames(hex_codes, model_names)

df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics.csv")
df <- df %>%
  filter(model %in% c("GB", "LR", "RF")) 
df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))

################################################################################

# The below chunk plots comparisons between the 100 iteration shuffled bayesian optimized
# hyperparameter tuning and models trained without any hyperparameter tuning.

df$hpt <- 1
df <- df %>%
  filter(gen == 'all')

MLmaf01all <- "/work/tfs3/gsAI/MLmodels/gebvs/defaultHP0.01_fold_metrics.csv"
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
MAF01ALL$hpt <- 0
MAF01ALL <- MAF01ALL %>%
  select(-auc_iter)

MLmaf05all <- "/work/tfs3/gsAI/MLmodels/gebvs/defaultHP0.05_fold_metrics.csv"
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
MAF05ALL$hpt <- 0
MAF05ALL <- MAF05ALL %>%
  select(-auc_iter)

MLmaf005all <- "/work/tfs3/gsAI/MLmodels/gebvs/defaultHP0.005_fold_metrics.csv"
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
MAF005ALL$hpt <- 0
MAF005ALL <- MAF005ALL %>%
  select(-auc_iter)

df <- rbind(df, MAF01ALL, MAF05ALL, MAF005ALL)
df$hpt <- factor(df$hpt, levels = unique(sort(df$hpt)))
df <- df %>%
  mutate(hpt = factor(hpt, levels = c(0, 1), labels = c("No", "Yes")))


maf05bp <- ggplot(df[df$MAF == '0.05', ], aes(x = hpt, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model), alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(aes(fill = model),
              shape = 21,
              color = "transparent",
              alpha = 0.4,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  geom_jitter(aes(color = model),
              shape = 21,
              fill = NA,
              stroke = 0.8,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  scale_color_manual(values = model_color_palette) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model) +
  geom_signif(
    comparisons = list(c(levels(df$hpt)[1], levels(df$hpt)[2])),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Hyperparameter Tuning", y = "Correlation Accuracy") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
maf05bp

maf01bp <- ggplot(df[df$MAF == '0.01', ], aes(x = hpt, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model), alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(aes(fill = model),
              shape = 21,
              color = "transparent",
              alpha = 0.4,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  geom_jitter(aes(color = model),
              shape = 21,
              fill = NA,
              stroke = 0.8,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  scale_color_manual(values = model_color_palette) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model) +
  geom_signif(
    comparisons = list(c(levels(df$hpt)[1], levels(df$hpt)[2])),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Hyperparameter Tuning", y = "Correlation Accuracy") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
maf01bp

maf005bp <- ggplot(df[df$MAF == '0.005', ], aes(x = hpt, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model), alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(aes(fill = model),
              shape = 21,
              color = "transparent",
              alpha = 0.4,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  geom_jitter(aes(color = model),
              shape = 21,
              fill = NA,
              stroke = 0.8,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  facet_wrap(~model) +
  scale_color_manual(values = model_color_palette) + 
  scale_fill_manual(values = model_color_palette) +
  geom_signif(
    comparisons = list(c(levels(df$hpt)[1], levels(df$hpt)[2])),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Hyperparameter Tuning", y = "Correlation Accuracy") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
maf005bp

combined_hpt <- plot_grid(
  maf05bp, maf01bp, maf005bp, 
  ncol = 3, 
  labels = c("A", "B", "C"),     
  label_size = 18,               
  rel_widths = c(1, 1, 1),        
  align = 'h',                  
  axis = 'bt' 
)
combined_hpt
ggsave("/work/tfs3/gsAI/analysis/pdfs/HPTinfluence.pdf", combined_hpt, width = 11, height = 6, dpi = 500)

################################################################################

# The below chunk tests influence of nested hyperparameter tuning on model performance.

MLmaf005all <- "/work/tfs3/gsAI/MLmodels/gebvs/unused/MAF0.005All/MAF0.005All_nested_metrics.csv"
MLmaf005all <- read.csv(MLmaf005all)
names(MLmaf005all) <- tolower(names(MLmaf005all))
MLmaf005all <- MLmaf005all %>%
  rename(correlation = pearsonr)
MLmaf005all <- MLmaf005all %>%
  group_by(model) %>%
  mutate(iteration = rep(1:10, each = 5, length.out = n())) %>%
  ungroup()
MLmaf005all <- MLmaf005all %>%
  group_by(model, iteration) %>%
  summarise(
    auc_iter = mean(auc, na.rm = TRUE),
    corr_iter = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  )

#### For matching and bind with main df
MLmaf005all$MAF <- 0.005
MLmaf005all$gen <- "all"
MLmaf005all$nested <- 1
MLmaf005all <- MLmaf005all %>%
  select(-auc_iter)

df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics.csv")
df <- df %>%
  filter(model %in% c("GB", "LR", "RF")) 
df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))
df$nested <- 0
df <- rbind(df,MLmaf005all)
df$nested <- factor(df$nested, levels = unique(sort(df$nested)))


maf005bp <- ggplot(df[df$MAF == '0.005' & df$gen=='all', ], aes(x = nested, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model), alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(fill = model),
              shape = 21,
              color = "transparent",
              alpha = 0.4,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  geom_jitter(aes(color = model),
              shape = 21,
              fill = NA,
              stroke = 0.8,
              size = 2.7,
              position = position_jitter(width = 0.2, seed = 123)) +
  facet_wrap(~model, scales = "free") +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  geom_signif(
    comparisons = list(c(levels(df$nested)[1], levels(df$nested)[2])),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  labs(x = "Nested Hyperparameter Tuning", y = "Correlation Accuracy") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none")
maf005bp
ggsave("/work/tfs3/gsAI/analysis/pdfs/NestedHPT.pdf", maf005bp, width = 10, height = 7, dpi = 500)
