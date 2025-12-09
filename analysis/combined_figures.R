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
library(ggsignif)
library(ggthemes)
library(multcompView)

df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics.csv")
df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))
model_names <- c("BayesB", "BRR", "EGBLUP", "GB", "GBLUP", 
                 "LASSO", "LR", "RF", "RKHS")
hex_codes <- viridis(n = 9, alpha = 0.6)
model_color_palette <- setNames(hex_codes, model_names)


################################################################################

##### plot boxplot all generations faceted by MAF
my_plot<-ggplot(df[df$gen == 'all', ], aes(x = corr_iter, y = model, fill = model)) +
    geom_boxplot(orientation = "y") +
    facet_wrap(~ MAF,ncol=1,scales = "free_y") + 
    scale_fill_manual(values = model_color_palette) +
    geom_jitter(color = "black", alpha = 0.3, size = 1) +
    labs(x = "Correlation", y = "Model", fill="Model") +
    theme_classic(base_size = 12) +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/misc/boxplot.png", my_plot, width = 6, height = 8, units = "in")

###### point plot all generations faceted by MAF 
summary_by_model <- df %>%
  pivot_longer(cols = c(corr_iter), names_to = "metric", values_to = "iter_mean") %>%
  group_by(MAF, gen, model, metric) %>%
  summarise(
    mean_of_iters = mean(iter_mean, na.rm = TRUE),    # single point to plot
    sd_of_iters   = sd(iter_mean, na.rm = TRUE),      # used for error bars
    n_iters       = n(),
    .groups = "drop"
  )
new_labels <- c("0.005" = "MAF 0.005", "0.05" = "MAF 0.05", "0.01" = "MAF 0.01")

my_plot2<-ggplot(summary_by_model[summary_by_model$gen == 'all', ], aes(x = model, y = mean_of_iters, color = model)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
                width = 0.2, color="black",linewidth = 0.5) +
  facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(y = "Correlation",color="Model") +
  theme(strip.text = element_text(size = 12)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("/work/tfs3/gsAI/analysis/misc/point.svg", my_plot2, width = 8, height = 5, units = "in")

my_plot3<-ggplot(summary_by_model[summary_by_model$gen == 'all', ], aes(x = model, y = mean_of_iters, color = model)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
                width = 0.2, color="black",linewidth = 0.5) +
  facet_wrap(~ MAF, labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(size = 12)) + 
  labs(y = "Correlation",color="Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("/work/tfs3/gsAI/analysis/misc/point2.svg", my_plot3, width = 8, height = 5, units = "in")

###### point plot all generations MAF0.01
MAF01df <- summary_by_model[summary_by_model$gen == 'all' & summary_by_model$MAF == 0.01,]
ggplot(MAF01df, aes(x = model, y = mean_of_iters, color = model)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
                width = 0.2, color="black", linewidth = 0.5) +
  scale_color_manual(values = model_color_palette) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(y = "Correlation",color="Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())


##### plot all MAF facet by generation
maf05bp <- ggplot(df[df$MAF == '0.05', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/misc/MAF05boxplot.png", maf05bp, width = 7, height = 5, units = "in")

maf01bp <- ggplot(df[df$MAF == '0.01', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/misc/MAF01boxplot.png", maf01bp, width = 7, height = 5, units = "in")

maf005bp <- ggplot(df[df$MAF == '0.005', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/misc/MAF005boxplot.png", maf005bp, width = 7, height = 5, units = "in")


##### plot all ML models facet by generation
ML_models <- c("GB", "LR", "RF")
ML_df <- df[df$gen == 'all' & df$model %in% ML_models, ]
annotation_data <- ML_df %>% filter(model == "RF") %>% distinct(model)
ggplot(ML_df, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(alpha = 0.5, aes(fill = model), outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  geom_jitter(color = "black", alpha = 0.7, size = 2) +
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "t.test",# Placeholder comparison list
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05), # Define significance levels
    step_increase = 0.1
  ) +
  labs(x = "Minor Allele Frequency", y = "Correlation") +
  theme_classic(base_size = 10) + 
  geom_text(
    data = annotation_data, 
    inherit.aes = FALSE,
    x = Inf, 
    y = 0.23, 
    label = "*** = p < 0.001\n** = p < 0.01\n* = p < 0.05", 
    size = 3.5, 
    hjust = 1.05,
    vjust = -0.5
  ) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(strip.text = element_text(size = 12)) +
  theme(legend.position = "none") + 
  facet_wrap(~model,scale="free")

R_models <- c("GBLUP", "LASSO", "EGBLUP", "BayesB", "BRR", "RKHS")
R_models <- df[df$gen == 'all' & df$model %in% R_models, ]
annotation_data <- R_models %>% filter(model == "GBLUP") %>% distinct(model)
ggplot(R_models, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(alpha = 0.5, aes(fill = model), outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  geom_jitter(color = "black", alpha = 0.7, size = 0.5) +
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  labs(x = "Minor Allele Frequency", y = "Correlation") +
  theme_classic(base_size = 10) + 
  geom_text(
    data = annotation_data, 
    inherit.aes = FALSE,
    x = 0.67, 
    y = 0.235, 
    label = "*** = p < 0.001\n** = p < 0.01\n* = p < 0.05", 
    size = 3, 
    hjust = 0,
    vjust = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(strip.text = element_text(size = 10)) +
  theme(legend.position = "none") + 
  facet_wrap(~model,scale="free")


################################################################################
MAF01df <- df[df$MAF == '0.01', ]
MAF005df <- df[df$MAF == '0.005', ]
MAF05df <- df[df$MAF == '0.05', ]


kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
nonsig <- dunn_table %>%
  dplyr::filter(P.adj > 0.05)
print(nonsig)

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
nonsig <- dunn_table %>%
  dplyr::filter(P.adj > 0.05)
print(nonsig)

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
nonsig <- dunn_table %>%
  dplyr::filter(P.adj > 0.05)
print(nonsig)
