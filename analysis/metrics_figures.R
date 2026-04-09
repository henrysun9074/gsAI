###############################################################################

## Set WD and load packages + data

setwd("/work/tfs3/gsAI")
library(tidyverse)
library(svglite)
library(readxl)
library(vegan)
library(ggpubr)
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
library(ggsci)
library(lsr)
library(scales)
library(cowplot)
library(rstatix)
library(tidyplots)
library(ggrepel)
library(ggstats)
library(ggalign)

df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics.csv")
df <- df %>%
  mutate(model = dplyr::recode(model,
                        "BB" = "BayesB"))
df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))
model_names <- c("GB", "LR", "RF", "BayesB", "BRR", "EGBLUP","GBLUP", 
                 "LASSO","RKHS")
hex_codes <-c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", 
              "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF")
model_color_palette <- setNames(hex_codes, model_names)

priority_models <- c("GB", "LR", "RF")
all_models <- unique(as.character(df$model))
other_models <- setdiff(all_models, priority_models)
new_model_order <- c(priority_models, other_models)
df$model <- factor(df$model, levels = new_model_order)

extra_df <- read.csv("/work/tfs3/gsAI/data/combined_fold_metrics_extra.csv")
extra_df <- extra_df %>%
  mutate(model = dplyr::recode(model,
                               "BB" = "BayesB"))
extra_df$MAF <- factor(extra_df$MAF, levels = unique(sort(extra_df$MAF)))
extra_df$MAF <- factor(extra_df$MAF, levels = unique(sort(extra_df$MAF)))
extra_df$extra <- 1
extra_df$model <- factor(extra_df$model, levels = new_model_order)

new_labels <- c("0.005" = "MAF 0.005", "0.05" = "MAF 0.05", "0.01" = "MAF 0.01")

################################################################################

# boxplot all generations
df$gen <- factor(df$gen, levels = c("F2", "all"))
gen_colors <- c("F2" = "#D2A6B4FF", "all" = "#8E2043FF")

bpF2vAll <- ggplot(df, aes(x = MAF, y = corr_iter, fill = gen, color = gen)) +
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15, 
      dodge.width = 0.8, 
      seed = 123
    ),
    shape = 21,
    size = 2,
    alpha = 0.4,
    stroke = 0.5
  ) +
  facet_wrap(~model, scales = "free_y") +
  stat_compare_means(
    aes(group = gen), 
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = FALSE,
    label.y.npc = "top",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                        symbols = c("***", "**", "*", "ns"))
  ) +
  scale_fill_manual(values = gen_colors) +
  scale_color_manual(values = gen_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr() +
  labs(
    x = "Minor Allele Frequency",
    y = "Correlation Accuracy",
    fill = "Generation",
    color = "Generation"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size =14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/AllMAFsF2vsAll.pdf", bpF2vAll, width = 12, height = 10, dpi = 300)
# *p<.05, **p<.01, ***p<0.001

################################################################################

### point plots for F2 generation only

MAF01df <- df[df$gen == 'F2' & df$MAF == '0.01', ]
MAF005df <- df[df$gen == 'F2' & df$MAF == '0.005',]
MAF05df <- df[df$gen == 'F2' & df$MAF == '0.05', ]

## kruskal and posthoc dunn tests for each MAF
kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
# Kruskal-Wallis chi-squared = 77.601, df = 8, p-value = 1.484e-13
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="none")
dunn_table <- dunn_res$res
model_order_index <- setNames(seq_along(new_model_order), new_model_order)
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
F2_dunn_MAF005 <- dunn_table_ordered

CLD1 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD1

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
# Kruskal-Wallis chi-squared = 72.596, df = 8, p-value = 1.492e-12
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
F2_dunn_MAF01 <- dunn_table_ordered

CLD2 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD2

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
# Kruskal-Wallis chi-squared = 59.59, df = 8, p-value = 5.61e-10
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
F2_dunn_MAF05 <- dunn_table_ordered

CLD3 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD3

################################
CLD_list <- list(
  list(cld_result = CLD1, MAF = '0.005'),
  list(cld_result = CLD2, MAF = '0.01'),
  list(cld_result = CLD3, MAF = '0.05')
)
all_CLD_df <- bind_rows(
  lapply(CLD_list, function(cld_item) {
    as.data.frame(cld_item$cld_result) %>%
      dplyr::rename(model = Group) %>%
      dplyr::mutate(
        MAF = cld_item$MAF,
        model = factor(model, levels = new_model_order)
      )
  })
)

###### revert back to df for generation comparison
summary_by_model <- df %>%
  pivot_longer(cols = c(corr_iter), names_to = "metric", values_to = "iter_mean") %>%
  group_by(MAF, gen, model, metric) %>%
  summarise(
    mean_of_iters = mean(iter_mean, na.rm = TRUE),    # single point to plot
    sd_of_iters   = sd(iter_mean, na.rm = TRUE),      # used for error bars
    n_iters       = n(),
    .groups = "drop"
  )
plot_data_with_cld_all <- summary_by_model[summary_by_model$gen == 'all', ] %>%
  dplyr::left_join(all_CLD_df, by = c("model", "MAF")) %>%
  dplyr::mutate(
    CLD_y_pos = mean_of_iters + sd_of_iters + 0.005
  )

### free y scale
# cld_point <- ggplot(
#   summary_by_model[summary_by_model$gen == 'all', ],
#   aes(x = model, y = mean_of_iters, color = model)
# ) +
#   geom_point(size = 5) +
#   geom_errorbar(
#     aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
#     width = 0.2, color = "black", linewidth = 0.5
#   ) +
#   facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
#   scale_color_manual(values = model_color_palette) +
#   theme_pubr() +
#   theme(strip.text = element_text(size = 12)) +
#   labs(y = "Correlation Accuracy", color = "Model") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   geom_text(
#     data = plot_data_with_cld_all, 
#     aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
#     inherit.aes = FALSE,
#     color = "black",
#     size = 4,
#     vjust = 0
#   )
# ggsave("/work/tfs3/gsAI/analysis/misc/point_f2_cld_all_facets.png", cld_point, width = 8, height = 5, units = "in")

### fixed y scale
cld_point2 <- ggplot(
  summary_by_model[summary_by_model$gen == 'all', ],
  aes(x = model, y = mean_of_iters, color = model)
) +
  geom_point(size = 5) +
  geom_errorbar(
    aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
    width = 0.2, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ MAF, labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_pubr() +
  theme(strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(y = "Correlation Accuracy", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
    data = plot_data_with_cld_all, 
    aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = 0
  )
cld_point2 <- ggdraw(cld_point2) +
  draw_label("KW p < 0.001", x = 0.9, y = 0.04, hjust = 0.5, vjust = 0)
ggsave("/work/tfs3/gsAI/analysis/pdfs/F2atallMAFs.pdf", cld_point2, width = 8, height = 5, units = "in")


################################################################################

##### plot F2 vs all performance at different MAFs

maf05bp <- ggplot(df[df$MAF == '0.05', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.8, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
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
  scale_color_manual(values = model_color_palette, guide = "none") +
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF05F2vsAll.pdf", maf05bp, width = 7, height = 6, units = "in")

maf01bp <- ggplot(df[df$MAF == '0.01', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
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
  scale_color_manual(values = model_color_palette, guide = "none") +
  facet_wrap(~model, scale="free") +
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01F2vsAll.pdf", maf01bp, width = 7, height = 6, units = "in")

maf005bp <- ggplot(df[df$MAF == '0.005', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
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
  scale_color_manual(values = model_color_palette, guide = "none") +
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005F2vsAll.pdf", maf005bp, width = 7, height = 6, units = "in")

################################################################################

##### plot all ML models facet by MAF
ML_models <- c("GB", "LR", "RF")
ML_df <- df[df$gen == 'all' & df$model %in% ML_models, ]
annotation_data <- ML_df %>% filter(model == "RF") %>% distinct(model)
ML_all_bp <- ggplot(ML_df, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
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
  scale_color_manual(values = model_color_palette, guide = "none") +
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05), # Define significance levels
    step_increase = 0.1
  ) +
  labs(x = "Minor Allele Frequency", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) + 
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "none") + 
  facet_wrap(~model,scale="free")
ML_all_bp

##### plot all R models facet by MAF
R_models <- c("GBLUP", "LASSO", "EGBLUP", "BayesB", "BRR", "RKHS")
R_models <- df[df$gen == 'all' & df$model %in% R_models, ]
annotation_data <- R_models %>% filter(model == "GBLUP") %>% distinct(model)
R_all_bp <- ggplot(R_models, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
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
  scale_color_manual(values = model_color_palette, guide = "none") +
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  labs(x = "Minor Allele Frequency", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.position = "none") + 
  facet_wrap(~model,scale="free")
R_all_bp


## combine ML_all_bp and R_all_bp 
combined_plot <- plot_grid(
  ML_all_bp, 
  R_all_bp,
  labels = c("A", "B"),
  label_size = 18,
  ncol = 2)
ggsave("/work/tfs3/gsAI/analysis/pdfs/AllModelsMAFboxplot.pdf",
       combined_plot,
       width = 11,
       height = 7,
       units = "in",
       dpi = 300)

################################################################################

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

MAF01df <- df[df$gen == 'all' & df$MAF == '0.01', ]
MAF005df <- df[df$gen == 'all' & df$MAF == '0.005',]
MAF05df <- df[df$gen == 'all' & df$MAF == '0.05', ]

## calculate kruskal-wallis and dunn test results for 'all' generations at each MAF
kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
# Kruskal-Wallis chi-squared = 82.683, df = 8, p-value = 1.408e-14
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="none")
dunn_table <- dunn_res$res
model_order_index <- setNames(seq_along(new_model_order), new_model_order)
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
    dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
all_dunn_MAF005 <- dunn_table_ordered

CLD1 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD1

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
# Kruskal-Wallis chi-squared = 84.215, df = 8, p-value = 6.904e-15
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
all_dunn_MAF01 <- dunn_table_ordered

CLD2 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD2

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
# Kruskal-Wallis chi-squared = 77.059, df = 8, p-value = 1.907e-13
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
all_dunn_MAF05 <- dunn_table_ordered

CLD3 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD3

CLD_list <- list(
  list(cld_result = CLD1, MAF = '0.005'),
  list(cld_result = CLD2, MAF = '0.01'),
  list(cld_result = CLD3, MAF = '0.05')
)
all_CLD_df <- bind_rows(
  lapply(CLD_list, function(cld_item) {
    as.data.frame(cld_item$cld_result) %>%
      dplyr::rename(model = Group) %>%
      dplyr::mutate(
        MAF = cld_item$MAF,
        model = factor(model, levels = new_model_order)
      )
  })
)
plot_data_with_cld_all <- summary_by_model[summary_by_model$gen == 'all', ] %>%
  dplyr::left_join(all_CLD_df, by = c("model", "MAF")) %>%
  dplyr::mutate(
    CLD_y_pos = mean_of_iters + sd_of_iters + 0.005
  )

################################

# this draws the plot with variable y-axis scale for each facet
# cld_point <- ggplot(
#   summary_by_model[summary_by_model$gen == 'all', ],
#   aes(x = model, y = mean_of_iters, color = model)
# ) +
#   geom_point(size = 5) +
#   geom_errorbar(
#     aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
#     width = 0.2, color = "black", linewidth = 0.5
#   ) +
#   facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
#   scale_color_manual(values = model_color_palette) +
#   theme_pubr() +
#   theme(strip.text = element_text(size = 12)) +
#   labs(y = "Correlation Accuracy", color = "Model") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   geom_text(
#     data = plot_data_with_cld_all, 
#     aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
#     inherit.aes = FALSE,
#     color = "black",
#     size = 4,
#     vjust = 0
#   ) 
# cld_point <- ggdraw(cld_point) +
#   draw_label("KW p < 0.001", x = 0.9, y = 0.05, hjust = 0.5, vjust = 0) +
#   theme_cowplot()
# ggsave("/work/tfs3/gsAI/analysis/misc/point_cld_all_facets.png", cld_point, width = 8, height = 5, units = "in")

# this fixes the y-axis scale for all facets
cld_point2_all <- ggplot(
  summary_by_model[summary_by_model$gen == 'all', ],
  aes(x = model, y = mean_of_iters, color = model)
) +
  geom_point(size = 5) +
  geom_errorbar(
    aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
    width = 0.2, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ MAF, labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_pubr() +
  theme(strip.text = element_text(size = 12),
        legend.position = "top" # turn this off if just plotting this plot due to ggdraw
        ) +
  labs(y = "Correlation Accuracy", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
  data = plot_data_with_cld_all,
  aes(x = model, y = CLD_y_pos,
      label = Letter,
      group = MAF),
  inherit.aes = FALSE,
  color = "black",
  size = 4,
  vjust = 0
)
legend_combined <- get_legend(cld_point2_all)
cld_point2_all <- cld_point2_all + theme(legend.position = "none")
cld_point2_all <- ggdraw(cld_point2_all) +
draw_label("KW p < 0.001", x = 0.9, y = 0.05, hjust = 0.5, vjust = 0)
ggsave("/work/tfs3/gsAI/analysis/pdfs/AllModelsMAFpointplot.pdf", cld_point2_all, width = 8, height = 5, units = "in")

################################################################################

combined_pointplot <- plot_grid(
  legend_combined,                                    
  cld_point2 + theme(legend.position = "none"),   
  cld_point2_all + theme(legend.position = "none"),
  ncol = 1,
  labels = c("", "A", "B"), 
  label_size = 18,
  rel_heights = c(0.3, 1, 1),
  align = "v",      
  axis = "lr"       
)
ggsave("/work/tfs3/gsAI/analysis/pdfs/CombinedPointPlotF2All.pdf", combined_pointplot,
       width = 10, height = 8, dpi = 300)

################################################################################
################################################################################

### point plot with correlations per MAF but for extra alleles with imputation

summary_by_model <- extra_df %>%
  pivot_longer(cols = c(corr_iter), names_to = "metric", values_to = "iter_mean") %>%
  group_by(MAF, gen, model, metric) %>%
  summarise(
    mean_of_iters = mean(iter_mean, na.rm = TRUE),    
    sd_of_iters   = sd(iter_mean, na.rm = TRUE), 
    n_iters       = n(),
    .groups = "drop"
  )

MAF01df <- extra_df[extra_df$gen == 'all' & extra_df$MAF == '0.01', ]
MAF005df <- extra_df[extra_df$gen == 'all' & extra_df$MAF == '0.005',]
MAF05df <- extra_df[extra_df$gen == 'all' & extra_df$MAF == '0.05', ]

kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
# Kruskal-Wallis chi-squared = 87.739, df = 8, p-value = 1.337e-15
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="none")
dunn_table <- dunn_res$res
model_order_index <- setNames(seq_along(new_model_order), new_model_order)
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
imputed_dunn_MAF005 <- dunn_table_ordered

CLD1 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD1

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
# Kruskal-Wallis chi-squared = 86.875, df = 8, p-value = 2e-15
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
imputed_dunn_MAF01 <- dunn_table_ordered


CLD2 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD2

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
# Kruskal-Wallis chi-squared = 87.974, df = 8, p-value = 1.198e-15
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="none")
dunn_table <- dunn_res$res
dunn_table_ordered <- dunn_table %>%
  tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
  dplyr::mutate(
    Model1_Rank = model_order_index[Model1],
    Model2_Rank = model_order_index[Model2]
  ) %>%
  dplyr::arrange(
    Model1_Rank,
    Model2_Rank
  ) %>% dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
imputed_dunn_MAF05 <- dunn_table_ordered

CLD3 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD3

CLD_list <- list(
  list(cld_result = CLD1, MAF = '0.005'),
  list(cld_result = CLD2, MAF = '0.01'),
  list(cld_result = CLD3, MAF = '0.05')
)
all_CLD_df <- bind_rows(
  lapply(CLD_list, function(cld_item) {
    as.data.frame(cld_item$cld_result) %>%
      dplyr::rename(model = Group) %>%
      dplyr::mutate(
        MAF = cld_item$MAF,
        model = factor(model, levels = new_model_order)
      )
  })
)
plot_data_with_cld_all <- summary_by_model[summary_by_model$gen == 'all', ] %>%
  dplyr::left_join(all_CLD_df, by = c("model", "MAF")) %>%
  dplyr::mutate(
    CLD_y_pos = mean_of_iters + sd_of_iters + 0.005
  )

################################

# free y scale for facets
# cld_point <- ggplot(
#   summary_by_model[summary_by_model$gen == 'all', ],
#   aes(x = model, y = mean_of_iters, color = model)
# ) +
#   geom_point(size = 5) +
#   geom_errorbar(
#     aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
#     width = 0.2, color = "black", linewidth = 0.5
#   ) +
#   facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
#   scale_color_manual(values = model_color_palette) +
#   theme_pubr() +
#   theme(strip.text = element_text(size = 12)) +
#   labs(y = "Correlation Accuracy", color = "Model") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   geom_text(
#     data = plot_data_with_cld_all, 
#     aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
#     inherit.aes = FALSE,
#     color = "black",
#     size = 4,
#     vjust = 0
#   ) 
# cld_point <- ggdraw(cld_point) +
#   draw_label("KW p < 0.001", x = 0.9, y = 0.05, hjust = 0.5, vjust = 0) +
#   theme_cowplot()
# ggsave("/work/tfs3/gsAI/analysis/pdfs/ImputedSNPsMAFpointplot.pdf", cld_point, width = 8, height = 5, units = "in")

# fixed y scale for facets
cld_point2 <- ggplot(
  summary_by_model[summary_by_model$gen == 'all', ],
  aes(x = model, y = mean_of_iters, color = model)
) +
  geom_point(size = 5) +
  geom_errorbar(
    aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
    width = 0.2, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ MAF, labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_pubr() +
  theme(strip.text = element_text(size = 12)) +
  labs(y = "Correlation Accuracy", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
    data = plot_data_with_cld_all, 
    aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = 0
  )
cld_point2 <- ggdraw(cld_point2) +
  draw_label("KW p < 0.001", x = 0.17, y = 0.04, hjust = 0.5, vjust = 0) +
  theme_cowplot()
ggsave("/work/tfs3/gsAI/analysis/pdfs/ImputedSNPsMAFpointplot.pdf", cld_point2, width = 8, height = 5, units = "in")

################################################################################
################################################################################

## combined between 
df$extra <- 0
combined_df_all <- rbind(df,extra_df)

################################################################################

# box plot all generations, at each MAF for default vs GSM markers
combined_df_allgens <- combined_df_all %>%
  filter(gen == "all") %>%
  mutate(extra = factor(extra, levels = c(0, 1), labels = c("Default", "With GSM")))

extra_colors <- c("Default" = "#BAB97DFF", "With GSM" = "#426737FF")

bpDefaultvsGSM <- ggplot(combined_df_allgens, aes(x = MAF, y = corr_iter, fill = extra, color = extra)) +
  geom_boxplot(
    position = position_dodge(width = 0.8), 
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15, 
      dodge.width = 0.8, 
      seed = 123
    ),
    shape = 21,
    size = 2,
    alpha = 0.4,
    stroke = 0.5
  ) +
  facet_wrap(~model, scales = "free_y") +
  stat_compare_means(
    aes(group = extra), 
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = FALSE,
    label.y.npc = "top",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                       symbols = c("***", "**", "*", "ns"))
  ) +
  scale_fill_manual(values = extra_colors) +
  scale_color_manual(values = extra_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_pubr() +
  labs(
    x = "Minor Allele Frequency",
    y = "Correlation Accuracy",
    fill = "Dataset",
    color = "Dataset"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size =14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/AllMAFsDefaultvsImputed.pdf", bpDefaultvsGSM, width = 12, height = 10, dpi = 300)
# *p<.05, **p<.01, ***p<0.001

#############################
# MAF05
MAF05_combined_df <- rbind(df, extra_df) %>%
  filter(gen == "all", MAF == 0.05) %>% ### modify MAF here
  mutate(extra = factor(extra, levels = c(0, 1), labels = c("Default", "With GSM")))

MAF05_extra_plot <- ggplot(MAF05_combined_df, aes(x = extra, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(
    x = "Dataset", 
    y = "Correlation Accuracy"
  ) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size=12)
  ) +
  geom_signif(
    comparisons = list(c("Default", "With GSM")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "N.S."=2),
    step_increase = 0.1
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF05Imputed.pdf", MAF05_extra_plot, width = 8, height = 5, units = "in", dpi = 300)

#############################
# MAF01
MAF01_combined_df <- rbind(df, extra_df) %>%
  filter(gen == "all", MAF == 0.01) %>% ### modify MAF here
  mutate(extra = factor(extra, levels = c(0, 1), labels = c("Default", "With GSM")))

MAF01_extra_plot <- ggplot(MAF01_combined_df, aes(x = extra, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(
    x = "Dataset", 
    y = "Correlation Accuracy"
  ) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size=12)
  ) +
  geom_signif(
    comparisons = list(c("Default", "With GSM")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "N.S."=2),
    step_increase = 0.1
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01Imputed.pdf", MAF01_extra_plot, width = 8, height = 5, units = "in", dpi = 300)


#############################
# MAF005
MAF005_combined_df <- rbind(df, extra_df) %>%
  filter(gen == "all", MAF == 0.005) %>% ### modify MAF here
  mutate(extra = factor(extra, levels = c(0, 1), labels = c("Default", "With GSM")))

MAF005_extra_plot <- ggplot(MAF005_combined_df, aes(x = extra, y = corr_iter)) +
  geom_boxplot(aes(fill = model, color = model),alpha = 0.6, outlier.shape = NA) + 
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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(
    x = "Dataset", 
    y = "Correlation Accuracy"
  ) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size=12)
  ) +
  geom_signif(
    comparisons = list(c("Default", "With GSM")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "N.S."=2),
    step_increase = 0.1
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005Imputed.pdf", MAF005_extra_plot, width = 8, height = 5, units = "in", dpi = 300)


################################################################################ 
##### plot all models facet by MAF with extra dataset
Extra_models <- c("GBLUP", "LASSO", "EGBLUP", "BayesB", "BRR", "RKHS")
Extra_models <- extra_df[extra_df$gen == 'all', ]
annotation_data <- Extra_models %>% filter(model == "GBLUP") %>% distinct(model)
Extra_models_bp <- ggplot(Extra_models, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(aes(color = model, fill = model), alpha = 0.6, outlier.shape = NA) + 
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
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  labs(x = "Minor Allele Frequency", y = "Correlation Accuracy") +
  theme_pubr(base_size = 12) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(strip.text = element_text(size = 10)) +
  theme(legend.position = "none") + 
  facet_wrap(~model,scale="free")
Extra_models_bp
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAFAllImputedboxplot.pdf", Extra_models_bp, width = 8, height = 6, units = "in", dpi = 300)

combined_GSM_plot <- plot_grid(cld_point2, Extra_models_bp,
                             labels = c("A", "B"),
                             label_size = 18,
                             ncol = 1)
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAFAllImputed.pdf", combined_GSM_plot, 
       width = 8, height = 10, units = "in", dpi = 300)

################################################################################

# Table with average corr_iter at F2 vs all and at each MAF
model_summary_table <- df %>%
  group_by(MAF, model, gen) %>%
  summarise(
    mean_corr = mean(corr_iter, na.rm = TRUE),
    sd_corr   = sd(corr_iter, na.rm = TRUE),
    .groups = "drop"
  )
# write.csv(model_summary_table, "/work/tfs3/gsAI/analysis/stats/model_summary.csv")

# Table with average corr_iter at imputed and not imputed at each MAF
extra_model_summary_table <- extra_df %>%
  group_by(MAF, model, gen) %>%
  summarise(
    mean_corr = mean(corr_iter, na.rm = TRUE),
    sd_corr   = sd(corr_iter, na.rm = TRUE),
    .groups = "drop"
  )
# write.csv(extra_model_summary_table, "/work/tfs3/gsAI/analysis/stats/imputed_model_summary.csv")

################################################################################
# P-values KW/dunn between models at each MAF on all generation
all_dunn_MAF05$MAF <- 0.05
all_dunn_MAF01$ MAF <- 0.01
all_dunn_MAF005$MAF <- 0.005
dunn_results_all <- rbind(all_dunn_MAF05, all_dunn_MAF01, all_dunn_MAF005)
# write.csv(dunn_results_all, "/work/tfs3/gsAI/analysis/stats/dunn_results_allgen.csv")

# P-values KW/dunn between models at each MAF on F2 generation
F2_dunn_MAF05$MAF <- 0.05
F2_dunn_MAF01$MAF <- 0.01
F2_dunn_MAF005$MAF <- 0.005
dunn_results_f2 <- rbind(F2_dunn_MAF05, F2_dunn_MAF01, F2_dunn_MAF005)
# write.csv(dunn_results_f2, "/work/tfs3/gsAI/analysis/stats/dunn_results_f2gen.csv")

# P-values KW/dunn between models at each MAF on imputed data
imputed_dunn_MAF005$MAF <- 0.05
imputed_dunn_MAF05$MAF <- 0.01
imputed_dunn_MAF01$MAF <- 0.005
imputed_dunn <- rbind(imputed_dunn_MAF05, imputed_dunn_MAF01, imputed_dunn_MAF005)
# write.csv(dunn_results_f2, "/work/tfs3/gsAI/analysis/stats/dunn_results_imputed.csv")


################################################################################
# P-values wilcoxon between MAF for each model

maf_pvals <- df %>%
  filter(gen == "all") %>%
  group_by(model) %>%
  pairwise_wilcox_test(corr_iter ~ MAF) %>%
  ungroup()
maf_pvals
write.csv(maf_pvals, "/work/tfs3/gsAI/analysis/stats/MAF_wilcoxon.csv")

# P-values wilcoxon at each MAF between F2 and all for each model 
f2_vs_all_pvals <- df %>%
  group_by(MAF, model) %>%
  wilcox_test(corr_iter ~ gen) %>%
  add_significance()
f2_vs_all_pvals
write.csv(f2_vs_all_pvals, "/work/tfs3/gsAI/analysis/stats/F2vAll_wilcoxon.csv")

# P-values wilcoxon at each MAF between default and imputed for each model
final_df <- rbind(MAF005_combined_df, MAF01_combined_df, MAF05_combined_df)
imputed_pvals <- final_df %>%
  group_by(MAF, model) %>%
  wilcox_test(corr_iter ~ extra) %>%
  add_significance()
imputed_pvals
write.csv(imputed_pvals, "/work/tfs3/gsAI/analysis/stats/Imputed_wilcoxon.csv")


################################################################################
# Effect sizes

# Compute cohen's d
df_all <- df[df$gen == 'all',]
cohens_d_all <- df_all %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
write.csv(cohens_d_all, "/work/tfs3/gsAI/analysis/stats/allgens_cohensd.csv")

## Compute cohen's d for models after imputation
extra_df_all <- extra_df[extra_df$gen == 'all',]
cohens_d_extra <- extra_df_all %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
write.csv(cohens_d_extra, "/work/tfs3/gsAI/analysis/stats/imputation_cohensd.csv")

## Compute cohen's d for F2 models
df_f2 <- df[df$gen == 'F2',]
cohens_d_F2 <- df_f2 %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
write.csv(cohens_d_F2, "/work/tfs3/gsAI/analysis/stats/f2_cohensd.csv")

## saved effect sizes - GB has large effect in all pairwise comparisons