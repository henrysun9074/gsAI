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
# hex_codes <- turbo(n = 9, alpha = 0.8)
hex_codes <-c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", 
              "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF")
model_color_palette <- setNames(hex_codes, model_names)

priority_models <- c("GB", "LR", "RF")
all_models <- unique(as.character(df$model))
other_models <- setdiff(all_models, priority_models)
new_model_order <- c(priority_models, other_models)
df$model <- factor(df$model, levels = new_model_order)


################################################################################
new_labels <- c("0.005" = "MAF 0.005", "0.05" = "MAF 0.05", "0.01" = "MAF 0.01")

##### plot boxplot all generations faceted by MAF
my_plot<-ggplot(df[df$gen == 'all', ], aes(x = model, y = corr_iter, fill = model)) +
    geom_boxplot() +
    facet_wrap(~ MAF,ncol=3,scales = "free_y",labeller = as_labeller(new_labels)) + 
    scale_fill_manual(values = model_color_palette) +
    geom_jitter(color = "black", alpha = 0.8, size = 1) +
    labs(x = "Model", y = "Correlation", fill="Model") +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(size = 12)) + 
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")))  # Adjusted vjust, hjust, and added negative bottom margin to pull labels closer to the axis.
ggsave("/work/tfs3/gsAI/analysis/misc/boxplot.png", my_plot, width = 10, height = 8, units = "in")

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

##### plot all MAF facet by generation
maf05bp <- ggplot(df[df$MAF == '0.05', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 1, outlier.shape = NA) + 
  scale_fill_manual(values = model_color_palette) +
  facet_wrap(~model, scale="free") +
  geom_jitter(color = "black", alpha = 0.7, size = 0.7) + 
  geom_signif(
    comparisons = list(c("F2", "all")),
    test = "t.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05," "=2),
    tip_length=0,
    step_increase = 0.3
  ) + 
  labs(x = "Generation", y = "Correlation") +
  theme_classic(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme(legend.position = "none") 
ggsave("/work/tfs3/gsAI/analysis/misc/MAF05boxplot.png", maf05bp, width = 7, height = 5, units = "in")

maf01bp <- ggplot(df[df$MAF == '0.01', ], aes(x = gen, y = corr_iter)) +
  geom_boxplot(aes(fill = model),alpha = 1, outlier.shape = NA) + 
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
  geom_boxplot(aes(fill = model),alpha = 1, outlier.shape = NA) + 
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


##### plot all ML models facet by MAF
ML_models <- c("GB", "LR", "RF")
ML_df <- df[df$gen == 'all' & df$model %in% ML_models, ]
annotation_data <- ML_df %>% filter(model == "RF") %>% distinct(model)
ML_all_bp <- ggplot(ML_df, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(alpha = 1, aes(fill = model), outlier.shape = NA) + 
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
    y = 0.22, 
    label = "*** = p < 0.001\n** = p < 0.01\n* = p < 0.05\nt-test", 
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
ML_all_bp

##### plot all R models facet by MAF
R_models <- c("GBLUP", "LASSO", "EGBLUP", "BayesB", "BRR", "RKHS")
R_models <- df[df$gen == 'all' & df$model %in% R_models, ]
annotation_data <- R_models %>% filter(model == "GBLUP") %>% distinct(model)
R_all_bp <- ggplot(R_models, aes(x = MAF, y = corr_iter)) +
  geom_boxplot(alpha = 1, aes(fill = model), outlier.shape = NA) + 
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
    y = 0.2335, 
    label = "*** = p < 0.001\n** = p < 0.01\n* = p < 0.05\nt-test", 
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
R_all_bp

################################################################################
################################################################################

MAF01df <- df[df$gen == 'all' & df$MAF == '0.01', ]
MAF005df <- df[df$gen == 'all' & df$MAF == '0.005',]
MAF05df <- df[df$gen == 'all' & df$MAF == '0.05', ]

kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
CLD1 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD1

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
CLD2 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD2

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
cld_point <- ggplot(
  summary_by_model[summary_by_model$gen == 'all', ],
  aes(x = model, y = mean_of_iters, color = model)
) +
  geom_point(size = 5) +
  geom_errorbar(
    aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
    width = 0.2, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_classic() +
  theme(strip.text = element_text(size = 12)) +
  labs(y = "Correlation", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
    data = plot_data_with_cld_all, 
    aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = 0
  ) 
cld_point <- ggdraw(cld_point) +
  draw_label("KW p < 0.001", x = 0.9, y = 0.05, hjust = 0.5, vjust = 0) +
  theme_cowplot()
ggsave("/work/tfs3/gsAI/analysis/misc/point_cld_all_facets.png", cld_point, width = 8, height = 5, units = "in")

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
  theme_classic() +
  theme(strip.text = element_text(size = 12)) +
  labs(y = "Correlation", color = "Model") +
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
  draw_label("KW p < 0.001", x = 0.9, y = 0.05, hjust = 0.5, vjust = 0) +
  theme_cowplot()
ggsave("/work/tfs3/gsAI/analysis/misc/point2_cld_all_facets.png", cld_point2, width = 8, height = 5, units = "in")

##grid
test <- align_plots(cld_point,cld_point2,nrow = 2)
ggsave("/work/tfs3/gsAI/analysis/misc/test.png", test, width = 12, height = 8, units = "in")


################################################################################
################################################################################

MAF01df <- df[df$gen == 'F2' & df$MAF == '0.01', ]
MAF005df <- df[df$gen == 'F2' & df$MAF == '0.005',]
MAF05df <- df[df$gen == 'F2' & df$MAF == '0.05', ]

kruskal <- kruskal.test(corr_iter ~ model, data = MAF005df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF005df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
CLD1 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD1

kruskal <- kruskal.test(corr_iter ~ model, data = MAF01df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF01df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
CLD2 = cldList(P.adj ~ Comparison, data=dunn_table_ordered)
CLD2

kruskal <- kruskal.test(corr_iter ~ model, data = MAF05df)
dunn_res <- dunnTest(corr_iter ~ model, data = MAF05df, method="bh")
dunn_table <- dunn_res$res
dunn_table <- dunn_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
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
plot_data_with_cld_all <- summary_by_model[summary_by_model$gen == 'all', ] %>%
  dplyr::left_join(all_CLD_df, by = c("model", "MAF")) %>%
  dplyr::mutate(
    CLD_y_pos = mean_of_iters + sd_of_iters + 0.005
  )

cld_point <- ggplot(
  summary_by_model[summary_by_model$gen == 'all', ],
  aes(x = model, y = mean_of_iters, color = model)
) +
  geom_point(size = 5) +
  geom_errorbar(
    aes(ymin = mean_of_iters - sd_of_iters, ymax = mean_of_iters + sd_of_iters),
    width = 0.2, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ MAF, scales = "free_y", labeller = as_labeller(new_labels)) +
  scale_color_manual(values = model_color_palette) +
  theme_classic() +
  theme(strip.text = element_text(size = 12)) +
  labs(y = "Correlation", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
    data = plot_data_with_cld_all, 
    aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = 0
  )
ggsave("/work/tfs3/gsAI/analysis/misc/point_f2_cld_all_facets.png", cld_point, width = 8, height = 5, units = "in")

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
  theme_classic() +
  theme(strip.text = element_text(size = 12)) +
  labs(y = "Correlation", color = "Model") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_text(
    data = plot_data_with_cld_all, 
    aes(x = model, y = CLD_y_pos, label = Letter, group = MAF),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = 0
  )
ggsave("/work/tfs3/gsAI/analysis/misc/point2_f2_cld_all_facets.png", cld_point2, width = 8, height = 5, units = "in")


################################################################################

## Compute cohen's d for all models
df_all <- df[df$gen == 'all',]
cohens_d_all <- df_all %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)

## Compute cohen's d for F2 models
df_f2 <- df[df$gen == 'F2',]
cohens_d_F2 <- df_f2 %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)

## saved effect sizes - GB has large effect in all pairwise comparisons