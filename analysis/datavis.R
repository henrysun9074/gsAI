library(tidyverse)
library(svglite)
library(readxl)
library(vegan)
library(ggpubr)
library(RColorBrewer)
library(broom)
library(car)
library(rstatix)
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
library(tidyplots)
library(ggrepel)
library(ggstats)
library(ggalign)

################################################################################

setwd("/work/tfs3/gsAI/analysis/v2")

df <- readRDS("masterdf.rds")

df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))

model_names <- c("GB", "LR", "RF", "BayesB", "BRR", "EGBLUP","GBLUP", 
                 "LASSO","RKHS")
# hex_codes <- turbo(n = 9, alpha = 0.8)
hex_codes <-c("#5E81ACFF", "#8FA87AFF", "#BF616AFF", "#E7D202FF", "#7D5329FF", 
              "#F49538FF", "#66CDAAFF", "#D070B9FF", "#98FB98FF", "#FCA3B7FF")
model_color_palette <- setNames(hex_codes, model_names)

priority_models <- c("GB", "LR", "RF")
all_models <- c("BayesB", "BRR", "EGBLUP", "GBLUP", "LASSO", "RKHS")
other_models <- setdiff(all_models, priority_models)
new_model_order <- c(priority_models, other_models)
df$model <- factor(df$model, levels = new_model_order)

new_labels <- c("0.005" = "MAF 0.005", "0.05" = "MAF 0.05", "0.01" = "MAF 0.01")


################################################################################

### Figure S2 - gsm vs all
nogsm_df <- df %>%
  filter(gsm == 0, hpt == "100iter") %>%
  mutate(gen = dplyr::recode(gen, "all" = "All"))

nogsm_df$gen <- factor(nogsm_df$gen, levels = c("F2", "All"))
gen_colors <- c("F2" = "#D2A6B4FF", "All" = "#8E2043FF")

bpF2vAll <- ggplot(nogsm_df, aes(x = MAF, y = corr_iter, fill = gen, color = gen)) +
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
    strip.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("pdfs/FigureS2.png", bpF2vAll, width = 12, height = 10, dpi = 300)

################################################################################

# F2 and all boxplots at all MAF
df_plot <- df %>%
  filter(
    gsm == 0,
    hpt == "100iter",
    model %in% new_model_order
  ) %>%
  mutate(
    MAF = factor(
      as.character(MAF),
      levels = c("0.05", "0.01", "0.005")
    ),
    model = factor(
      model,
      levels = new_model_order
    ),
    gen = factor(
      gen,
      levels = c("F2", "all")
    )
  )

get_cld <- function(data, generation) {
  
  generation_data <- data %>%
    filter(gen == generation)
  
  cld_list <- lapply(
    levels(droplevels(generation_data$MAF)),
    function(maf_level) {
      
      tmp <- generation_data %>%
        filter(MAF == maf_level)
      kw <- kruskal.test(
        corr_iter ~ model,
        data = tmp
      )
      dunn_res <- FSA::dunnTest(
        corr_iter ~ model,
        data = tmp,
        method = "bh"
      )
      dunn_table <- dunn_res$res
      model_order_index <- setNames(
        seq_along(new_model_order),
        new_model_order
      )
      dunn_table_ordered <- dunn_table %>%
        tidyr::separate(
          Comparison,
          into = c("Model1", "Model2"),
          sep = " - ",
          remove = FALSE
        ) %>%
        mutate(
          Model1_Rank = model_order_index[Model1],
          Model2_Rank = model_order_index[Model2]
        ) %>%
        arrange(
          Model1_Rank,
          Model2_Rank
        )
      
      # Compact letter display
      cld_res <- rcompanion::cldList(
        P.adj ~ Comparison,
        data = dunn_table_ordered
      )
      
      cld_res <- as.data.frame(cld_res) %>%
        rename(model = Group) %>%
        mutate(
          MAF = maf_level,
          KW_p = kw$p.value
        )
      
      cld_res
    }
  )
  
  bind_rows(cld_list)
}

cld_F2 <- get_cld(df_plot, "F2")
cld_all <- get_cld(df_plot, "all")

make_generation_plot <- function(data, cld_data, generation) {
  
  plot_data <- data %>%
    filter(gen == generation)
  
  cld_positions <- plot_data %>%
    group_by(MAF, model) %>%
    summarise(
      model_max = max(corr_iter, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      cld_data %>% select(model, MAF, Letter),
      by = c("model", "MAF")
    ) %>%
    group_by(MAF) %>%
    mutate(
      facet_range = max(model_max, na.rm = TRUE) -
        min(model_max, na.rm = TRUE),
      facet_range = ifelse(
        facet_range == 0 | is.na(facet_range),
        0.01,
        facet_range
      ),
      CLD_y_pos = model_max + 0.08 * facet_range
    ) %>%
    ungroup()
  
  ggplot(
    plot_data,
    aes(
      x = model,
      y = corr_iter,
      color = model
    )
  ) +
    
    geom_boxplot(
      aes(fill = model),
      width = 0.65,
      alpha = 0.25,
      linewidth = 0.6,
      outlier.shape = NA
    ) +
    
    geom_point(
      position = position_jitter(
        width = 0.15,
        height = 0,
        seed = 123
      ),
      size = 1.7,
      alpha = 0.65
    ) +
    
    geom_text(
      data = cld_positions,
      aes(
        x = model,
        y = CLD_y_pos,
        label = Letter
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 4,
      vjust = 0
    ) +
    
    facet_wrap(
      ~ MAF,
      labeller = as_labeller(new_labels)
    ) +
    
    scale_color_manual(
      values = model_color_palette
    ) +
    
    scale_fill_manual(
      values = model_color_palette
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    
    labs(
      x = NULL,
      y = "Correlation Accuracy",
      color = "Model",
      fill = "Model"
    ) +
    
    theme_pubr() +
    
    theme(
      strip.text = element_text(size = 12),
      strip.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
}

cld_boxplot_F2 <- make_generation_plot(
  df_plot,
  cld_F2,
  "F2"
)

cld_boxplot_F2 <- ggdraw(cld_boxplot_F2) +
  draw_label(
    "KW p < 0.001",
    x = 0.9,
    y = 0.04,
    hjust = 0.5,
    vjust = 0
  )

cld_boxplot_all <- make_generation_plot(
  df_plot,
  cld_all,
  "all"
)

cld_boxplot_all <- ggdraw(cld_boxplot_all) +
  draw_label(
    "KW p < 0.001",
    x = 0.9,
    y = 0.04,
    hjust = 0.5,
    vjust = 0
  )
cld_boxplot_all

legend_plot <- ggplot(
  df_plot,
  aes(
    x = model,
    y = corr_iter,
    color = model,
    fill = model
  )
) +
  geom_point(size = 4) +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  guides(
    color = guide_legend(
      title = "Model",
      nrow = 1,
      override.aes = list(alpha = 1)
    ),
    fill = "none"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

legend_combined <- cowplot::get_legend(legend_plot)

combined_boxplot <- cowplot::plot_grid(
  legend_combined,
  cld_boxplot_F2,
  cld_boxplot_all,
  ncol = 1,
  labels = c("", "A", "B"),
  label_size = 18,
  rel_heights = c(0.18, 1, 1),
  align = "v"
)
combined_boxplot

ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/Figure1.pdf",
  combined_boxplot,
  width = 10,
  height = 8,
  dpi = 300
)

################################################################################

ML_models <- c("GB", "LR", "RF")

ML_df <- df %>%
  filter(
    gen == "all",
    model %in% ML_models,
    gsm == 0,
    hpt == "100iter"
  ) %>%
  mutate(
    model = factor(model, levels = ML_models)
  )

ML_all_bp <- ggplot(
  ML_df,
  aes(x = MAF, y = corr_iter)
) +
  geom_boxplot(
    aes(fill = model, color = model),
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(fill = model),
    shape = 21,
    color = "transparent",
    alpha = 0.4,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) + geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  geom_jitter(
    aes(color = model),
    shape = 21,
    fill = NA,
    stroke = 0.8,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) +
  scale_fill_manual(values = model_color_palette) +
  scale_color_manual(
    values = model_color_palette,
    guide = "none"
  ) +
  facet_wrap(~model, scales = "free_y") +
  labs(
    x = "Minor Allele Frequency",
    y = "Correlation Accuracy"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  )

ML_all_bp

R_models <- c(
  "GBLUP",
  "LASSO",
  "EGBLUP",
  "BayesB",
  "BRR",
  "RKHS"
)

R_df <- df %>%
  filter(
    gen == "all",
    model %in% R_models,
    gsm == 0,
    hpt == "100iter"
  ) %>%
  mutate(
    model = factor(model, levels = R_models)
  )

R_all_bp <- ggplot(
  R_df,
  aes(x = MAF, y = corr_iter)
) +
  geom_boxplot(
    aes(fill = model, color = model),
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(fill = model),
    shape = 21,
    color = "transparent",
    alpha = 0.4,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) +
  geom_jitter(
    aes(color = model),
    shape = 21,
    fill = NA,
    stroke = 0.8,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) + geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    step_increase = 0.3
  ) +
  scale_fill_manual(values = model_color_palette) +
  scale_color_manual(
    values = model_color_palette,
    guide = "none"
  ) +
  facet_wrap(~model, scales = "free_y") +
  labs(
    x = "Minor Allele Frequency",
    y = "Correlation Accuracy"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  )
R_all_bp

combined_plot <- plot_grid(
  ML_all_bp,
  R_all_bp,
  labels = c("A", "B"),
  label_size = 18,
  ncol = 2,
  align = "hv"
)
combined_plot

ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/Figure2.pdf",
  combined_plot,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)

################################################################################

# box plot all generations, at each MAF for default vs GSM markers
combined_df_allgens <- df %>%
  filter(gen == "all", hpt == "100iter") %>%
  mutate(gsm = factor(gsm, levels = c(0, 1), labels = c("Default", "With GSM")))

extra_colors <- c("Default" = "#BAB97DFF", "With GSM" = "#426737FF")

bpDefaultvsGSM <- ggplot(combined_df_allgens, aes(x = MAF, y = corr_iter, fill = gsm, color = gsm)) +
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
    aes(group = gsm), 
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
    strip.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("/work/tfs3/gsAI/analysis/v2/pdfs/Figure3.pdf", bpDefaultvsGSM, width = 12, height = 10, dpi = 300)

################################################################################

plot_data_all_gsm <- df %>%
  filter(
    gen == "all",
    gsm == 1,
    hpt == "100iter",
    model %in% new_model_order
  ) %>%
  mutate(
    MAF = factor(
      as.character(MAF),
      levels = c("0.005", "0.01", "0.05")
    ),
    model = factor(
      model,
      levels = new_model_order
    )
  )

### box plot with correlations per MAF but with GSM inclusion
MAF005_gsm <- plot_data_all_gsm %>% filter(MAF == "0.005")
MAF01_gsm  <- plot_data_all_gsm %>% filter(MAF == "0.01")
MAF05_gsm  <- plot_data_all_gsm %>% filter(MAF == "0.05")

model_order_index <- setNames(seq_along(new_model_order), new_model_order)

get_ordered_cld <- function(df_subset) {
  dunn_res <- dunnTest(corr_iter ~ model, data = df_subset, method = "bh")$res
  
  dunn_ordered <- dunn_res %>%
    tidyr::separate(Comparison, into = c("Model1", "Model2"), sep = " - ", remove = FALSE) %>%
    dplyr::mutate(
      Model1_Rank = model_order_index[Model1],
      Model2_Rank = model_order_index[Model2]
    ) %>%
    dplyr::arrange(Model1_Rank, Model2_Rank) %>%
    dplyr::select(-Model1, -Model2, -Model1_Rank, -Model2_Rank)
  
  cldList(P.adj ~ Comparison, data = dunn_ordered)
}

# MAF 0.005
kw_005 <- kruskal.test(corr_iter ~ model, data = MAF005_gsm)
CLD1   <- get_ordered_cld(MAF005_gsm)

# MAF 0.01
kw_01 <- kruskal.test(corr_iter ~ model, data = MAF01_gsm)
CLD2  <- get_ordered_cld(MAF01_gsm)

# MAF 0.05
kw_05 <- kruskal.test(corr_iter ~ model, data = MAF05_gsm)
CLD3  <- get_ordered_cld(MAF05_gsm)

CLD_list <- list(
  list(cld_result = CLD1, MAF = "0.005"),
  list(cld_result = CLD2, MAF = "0.01"),
  list(cld_result = CLD3, MAF = "0.05")
)

all_CLD_df <- bind_rows(
  lapply(CLD_list, function(cld_item) {
    as.data.frame(cld_item$cld_result) %>%
      dplyr::rename(model = Group) %>%
      dplyr::mutate(
        MAF = factor(cld_item$MAF, levels = c("0.05", "0.01", "0.005")),
        model = factor(model, levels = new_model_order)
      )
  })
)

plot_data_with_cld_gsm <- plot_data_all_gsm %>%
  dplyr::group_by(MAF, model) %>%
  dplyr::summarise(
    model_max = max(corr_iter, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(all_CLD_df, by = c("model", "MAF")) %>%
  dplyr::group_by(MAF) %>%
  dplyr::mutate(
    facet_range = diff(range(plot_data_all_gsm$corr_iter[
      plot_data_all_gsm$MAF == dplyr::first(MAF)
    ], na.rm = TRUE)),
    CLD_y_pos = model_max + 0.06 * facet_range
  ) %>%
  dplyr::ungroup()

cld_boxplot_all_gsm <- ggplot(
  plot_data_all_gsm,
  aes(x = model, y = corr_iter, color = model)
) +
  geom_boxplot(
    aes(fill = model),
    width = 0.65,
    alpha = 0.25,
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = 123
    ),
    size = 1.7,
    alpha = 0.65
  ) +
  geom_text(
    data = plot_data_with_cld_gsm,
    aes(
      x = model,
      y = CLD_y_pos,
      label = Letter
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3.8,
    vjust = 0
  ) +
  facet_wrap(
    ~ MAF,
    labeller = as_labeller(new_labels)
  ) +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.14))
  ) +
  labs(
    x = NULL,
    y = "Correlation Accuracy",
    color = "Model",
    fill = "Model"
  ) +
  theme_pubr() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  )

cld_boxplot_all_gsm <- ggdraw(cld_boxplot_all_gsm) +
  draw_label(
    "KW p < 0.001",
    x = 0.17,
    y = 0.04,
    hjust = 0.5,
    vjust = 0
  )

cld_boxplot_all_gsm

gsm_all_plot <- ggplot(
  plot_data_all_gsm,
  aes(x = MAF, y = corr_iter)
) +
  geom_boxplot(
    aes(color = model, fill = model),
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(fill = model),
    shape = 21,
    color = "transparent",
    alpha = 0.4,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) +
  geom_jitter(
    aes(color = model),
    shape = 21,
    fill = NA,
    stroke = 0.8,
    size = 2.7,
    position = position_jitter(width = 0.2, seed = 123)
  ) +
  facet_wrap(~model, scales = "free") +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  geom_signif(
    comparisons = list(c("0.005", "0.01"), c("0.005", "0.05"), c("0.05", "0.01")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = 0.3
  ) +
  labs(
    x = "Minor Allele Frequency",
    y = "Correlation Accuracy"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.2))
  ) +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill="white"),
    legend.position = "none"
  )
gsm_all_plot

model_legend_plot <- ggplot(
  plot_data_all_gsm, 
  aes(x = model, y = corr_iter, color = model, fill = model)
) +
  geom_point(shape = 21, size = 5, stroke = 1) +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  labs(color = "Model", fill = "Model") +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(size = 5, shape = 21, stroke = 1)
    ),
    fill = guide_legend(nrow = 1)
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key = element_blank()
  )

shared_model_legend <- cowplot::get_legend(model_legend_plot)

combined_GSM_boxplot <- plot_grid(
  shared_model_legend,
  cld_boxplot_all_gsm,
  gsm_all_plot,
  labels = c("", "A", "B"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(0.08, 0.6, 1.2)
)

ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/FigureS4.jpg",
  combined_GSM_boxplot, 
  width = 8,
  height = 9,
  units = "in",
  dpi = 300
)

################################################################################

# NEW plot of hpt influence on RF and GB performance

hpt_data <- df %>%
  filter(
    gen == "all",
    model %in% c("RF", "GB"),
    hpt %in% c("none", "100iter", "nested"),
    MAF %in% c(0.05, 0.01, 0.005)
  ) %>%
  mutate(
    hpt = factor(
      dplyr::recode(
        hpt,
        "none"    = "None",
        "100iter" = "100 Iteration",
        "nested"  = "Nested"
      ),
      levels = c("None", "100 Iteration", "Nested")
    ),
    MAF = factor(as.character(MAF), levels = c("0.05", "0.01", "0.005")),
    gsm_label = ifelse(gsm == 1 | gsm == "yes", "GSM (+)", "GSM (-)"),
    model = factor(model, levels = c("RF", "GB")),
    model_hpt = factor(
      paste(model, hpt, sep = "_"),
      levels = c(
        "RF_None", "RF_100 Iteration", "RF_Nested",
        "GB_None", "GB_100 Iteration", "GB_Nested"
      )
    )
  )

best_conv <- df %>%
  filter(
    gen == "all",
    ml == 0,
    MAF %in% c(0.05, 0.01, 0.005)
  ) %>%
  mutate(
    MAF = factor(as.character(MAF), levels = c("0.05", "0.01", "0.005")),
    gsm_label = ifelse(gsm == 1 | gsm == "yes", "GSM (+)", "GSM (-)")
  ) %>%
  group_by(MAF, gsm_label, model) %>%
  summarise(mean_acc = mean(corr_iter, na.rm = TRUE), .groups = "drop") %>%
  group_by(MAF, gsm_label) %>%
  slice_max(order_by = mean_acc, n = 1, with_ties = FALSE) %>%
  ungroup()

build_hpt_panel <- function(curr_maf, curr_gsm) {
  sub_data <- hpt_data %>%
    filter(MAF == curr_maf, gsm_label == curr_gsm)
  
  bench <- best_conv %>%
    filter(MAF == curr_maf, gsm_label == curr_gsm)
  
  hpt_comparisons <- list(
    c("RF_None", "RF_100 Iteration"),
    c("RF_100 Iteration", "RF_Nested"),
    c("RF_None", "RF_Nested"),
    c("GB_None", "GB_100 Iteration"),
    c("GB_100 Iteration", "GB_Nested"),
    c("GB_None", "GB_Nested")
  )
  
  p <- ggplot(sub_data, aes(x = model_hpt, y = corr_iter)) +
    geom_hline(
      yintercept = bench$mean_acc,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.75
    ) +
    geom_boxplot(
      aes(color = model, fill = model),
      alpha = 0.6,
      outlier.shape = NA,
      width = 0.65
    ) +
    geom_jitter(
      aes(fill = model),
      shape = 21,
      color = "transparent",
      alpha = 0.4,
      size = 2.0,
      position = position_jitter(width = 0.18, seed = 123)
    ) +
    geom_jitter(
      aes(color = model),
      shape = 21,
      fill = NA,
      stroke = 0.7,
      size = 2.0,
      position = position_jitter(width = 0.18, seed = 123)
    ) +
    # Significance tests
    geom_signif(
      comparisons = hpt_comparisons,
      test = "wilcox.test",
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      step_increase = 0.12,
      tip_length = 0.015,
      textsize = 3.2,
      color = "black"
    ) +
    scale_x_discrete(
      labels = c(
        "RF_None"          = "None",
        "RF_100 Iteration" = "100 Iterations",
        "RF_Nested"        = "Nested",
        "GB_None"          = "None",
        "GB_100 Iteration" = "100 Iterations",
        "GB_Nested"        = "Nested"
      )
    ) +
    scale_color_manual(values = model_color_palette) +
    scale_fill_manual(values = model_color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      x = NULL,
      y = "Correlation Accuracy",
      fill = "Model",
      color = "Model"
    ) +
    theme_pubr(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 8.5),
      axis.title.y = element_text(size = 10),
      legend.position = "none"
    )
  
  return(p)
}

pA <- build_hpt_panel("0.05",  "GSM (-)")
pB <- build_hpt_panel("0.01",  "GSM (-)")
pC <- build_hpt_panel("0.005", "GSM (-)")
pD <- build_hpt_panel("0.05",  "GSM (+)")
pE <- build_hpt_panel("0.01",  "GSM (+)")
pF <- build_hpt_panel("0.005", "GSM (+)")

legend_plot <- ggplot(hpt_data, aes(x = model, y = corr_iter, fill = model, color = model)) +
  geom_point(alpha = 1, size =5) +
  scale_color_manual(values = model_color_palette) +
  scale_fill_manual(values = model_color_palette) +
  labs(fill = "Model", color = "Model") +
  theme_pubr() +
  theme(legend.position = "bottom")

shared_legend <- get_legend(legend_plot)

panels_grid <- plot_grid(
  pA, pB, pC,
  pD, pE, pF,
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 3,
  align = "hv"
)

combined_hpt_figure <- plot_grid(
  shared_legend,
  panels_grid,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

combined_hpt_figure

ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/Figure4.pdf",
  plot = combined_hpt_figure,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300
)

################################################################################

### stats for paper

summary_stats <- df %>%
  group_by(model, gen, MAF, gsm, hpt) %>%
  summarise(
    n = n(),
    mean_corr   = round(mean(corr_iter, na.rm = TRUE), 3),
    sd_corr     = round(sd(corr_iter, na.rm = TRUE), 3),
    median_corr = median(corr_iter, na.rm = TRUE),
    IQR_corr    = IQR(corr_iter, na.rm = TRUE),
    .groups = "drop"
  )


## base subset for all following stats
df_100iter <- df %>% 
  filter(hpt == "100iter")

#############

run_kw_by_maf <- function(data_subset, dataset_name = "") {
  maf_levels <- c("0.05", "0.01", "0.005")
  
  lapply(maf_levels, function(m) {
    df_m <- data_subset %>% filter(as.character(MAF) == m)
    
    # Run Kruskal-Wallis test
    kw <- kruskal.test(corr_iter ~ model, data = df_m)
    
    data.frame(
      Dataset   = dataset_name,
      MAF       = as.numeric(m),
      statistic = round(kw$statistic, 3),
      df        = kw$parameter,
      p.value   = kw$p.value,
      p_format  = ifelse(kw$p.value < 0.001, "p < 0.001", sprintf("p = %.4f", kw$p.value))
    )
  }) %>% bind_rows()
}

# F2 generation (no GSM)
kw_f2_no_gsm <- df_100iter %>%
  filter(gen == "F2", gsm == 0) %>%
  run_kw_by_maf("F2 (No GSM)")

# All generations (no GSM)
kw_all_no_gsm <- df_100iter %>%
  filter(gen == "all", gsm == 0) %>%
  run_kw_by_maf("All (No GSM)")

# All generations (with GSM)
kw_all_gsm <- df_100iter %>%
  filter(gen == "all", gsm == 1) %>%
  run_kw_by_maf("All (With GSM)")

print(kw_f2_no_gsm)
print(kw_all_no_gsm)
print(kw_all_gsm)

# Helper function to compute Dunn's tests per MAF level
run_dunn_by_maf <- function(data_subset) {
  maf_levels <- c("0.05", "0.01", "0.005")
  
  lapply(maf_levels, function(m) {
    df_m <- data_subset %>% filter(as.character(MAF) == m)
    res <- dunnTest(corr_iter ~ model, data = df_m, method = "bh")$res
    res$MAF <- as.numeric(m)
    res
  }) %>% bind_rows()
}


### All generations with GSM included (gsm == 1)
dunn_results_all <- df_100iter %>%
  filter(gen == "all", gsm == 1) %>%
  run_dunn_by_maf()
dunn_results_all

# F2 generation no GSM
dunn_results_f2 <- df_100iter %>%
  filter(gen == "F2", gsm == 0) %>%
  run_dunn_by_maf()
dunn_results_f2

# All generations no GSM (gsm == 0)
dunn_results_no_gsm <- df_100iter %>%
  filter(gen == "all", gsm == 0) %>%
  run_dunn_by_maf()
dunn_results_no_gsm

################# wilcoxon tests

# Pairwise between MAF levels for each model (all generations, gsm == 1)
maf_pvals <- df_100iter %>%
  filter(gen == "all", gsm == 1) %>%
  group_by(model) %>%
  wilcox_test(corr_iter ~ MAF, p.adjust.method = "fdr") %>%
  ungroup()
maf_pvals

# F2 vs All generations at each MAF for each model
f2_vs_all_pvals <- df_100iter %>%
  filter(gen %in% c("all", "F2"), gsm == 0) %>%
  group_by(MAF, model) %>%
  wilcox_test(corr_iter ~ gen) %>%
  add_significance() %>%
  ungroup()
f2_vs_all_pvals

# GSM Effect: Without GSM (gsm = 0) vs With GSM (gsm = 1) at each MAF
gsm_pvals <- df_100iter %>%
  filter(gen == "all", !is.na(gsm)) %>%
  group_by(MAF, model) %>%
  wilcox_test(corr_iter ~ gsm) %>%
  add_significance() %>%
  ungroup()
gsm_pvals

####################################
# % change between MAF levels split by GSM status 
pct_change_maf <- df_100iter %>%
  filter(
    gen == "all",
    !is.na(gsm),
    MAF %in% c(0.05, 0.01, 0.005)
  ) %>%
  mutate(MAF = paste0("MAF_", MAF)) %>%
  select(gsm, model, iteration, MAF, corr_iter) %>%
  pivot_wider(
    id_cols = c(gsm, model, iteration),
    names_from = MAF,
    values_from = corr_iter
  ) %>%
  mutate(
    pct_05_to_01   = 100 * (MAF_0.01 - MAF_0.05) / MAF_0.05,
    pct_05_to_005  = 100 * (MAF_0.005 - MAF_0.05) / MAF_0.05,
    pct_01_to_005  = 100 * (MAF_0.005 - MAF_0.01) / MAF_0.01
  ) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "comparison",
    names_prefix = "pct_",
    values_to = "pct_change"
  ) %>%
  mutate(
    comparison = dplyr::recode(
      comparison,
      "05_to_01"  = "0.05 -> 0.01",
      "05_to_005" = "0.05 -> 0.005",
      "01_to_005" = "0.01 -> 0.005"
    ),
    gsm_label = ifelse(gsm == 1, "GSM (+)", "GSM (-)")
  )

summary_pct_maf <- pct_change_maf %>%
  group_by(gsm_label, model, comparison) %>%
  summarise(
    mean_pct_change = round(mean(pct_change, na.rm = TRUE), 1),
    sd_pct_change   = round(sd(pct_change, na.rm = TRUE), 1),
    n               = sum(!is.na(pct_change)),
    .groups         = "drop"
  ) %>%
  arrange(gsm_label, model, comparison)

summary_pct_maf

####################################
# % inc after GSMs
pct_change_df <- df_100iter %>%
  filter(gen == "all", !is.na(gsm)) %>%
  select(gen, MAF, model, iteration, gsm, corr_iter) %>%
  pivot_wider(
    id_cols = c(gen, MAF, model, iteration),
    names_from = gsm,
    values_from = corr_iter
  ) %>%
  rename(
    without_gsm = `0`,
    with_gsm    = `1`
  ) %>%
  mutate(
    pct_change = 100 * (with_gsm - without_gsm) / without_gsm
  )

summary_df <- pct_change_df %>%
  group_by(model, MAF) %>%
  summarise(
    mean_pct_change = round(mean(pct_change, na.rm = TRUE), 1),
    sd_pct_change   = round(sd(pct_change, na.rm = TRUE), 1),
    n               = sum(!is.na(pct_change)),
    .groups         = "drop"
  ) %>%
  arrange(model, MAF)
summary_df

####################################
pct_change_gen <- df_100iter %>%
  filter(gsm == 0, gen %in% c("F2", "all")) %>%
  select(MAF, model, iteration, gen, corr_iter) %>%
  pivot_wider(
    id_cols = c(MAF, model, iteration),
    names_from = gen,
    values_from = corr_iter
  ) %>%
  rename(
    f2_gen  = F2,
    all_gen = all
  ) %>%
  mutate(
    pct_change = 100 * (all_gen - f2_gen) / f2_gen
  )

summary_pct_gen <- pct_change_gen %>%
  group_by(model, MAF) %>%
  summarise(
    mean_pct_change = round(mean(pct_change, na.rm = TRUE), 1),
    sd_pct_change   = round(sd(pct_change, na.rm = TRUE), 1),
    n               = sum(!is.na(pct_change)),
    .groups         = "drop"
  ) %>%
  arrange(model, MAF)

summary_pct_gen


####################################
pct_change_hpt <- df %>%
  filter(
    gen == "all",
    model %in% c("RF", "GB"),
    hpt %in% c("none", "100iter", "nested")
  ) %>%
  select(gsm, MAF, model, iteration, hpt, corr_iter) %>%
  pivot_wider(
    id_cols = c(gsm, MAF, model, iteration),
    names_from = hpt,
    values_from = corr_iter
  ) %>%
  mutate(
    pct_none_to_100    = 100 * (`100iter` - none) / none,
    pct_none_to_nested = 100 * (nested - none) / none,
    pct_100_to_nested  = 100 * (nested - `100iter`) / `100iter`
  ) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "comparison",
    names_prefix = "pct_",
    values_to = "pct_change"
  ) %>%
  mutate(
    comparison = dplyr::recode(
      comparison,
      "none_to_100"    = "None -> 100iter",
      "none_to_nested" = "None -> Nested",
      "100_to_nested"  = "100iter -> Nested"
    )
  )

summary_pct_hpt <- pct_change_hpt %>%
  group_by(gsm, model, MAF, comparison) %>%
  summarise(
    mean_pct_change = round(mean(pct_change, na.rm = TRUE), 1),
    sd_pct_change   = round(sd(pct_change, na.rm = TRUE), 1),
    n               = sum(!is.na(pct_change)),
    .groups         = "drop"
  ) %>%
  arrange(gsm, model, MAF, comparison)
summary_pct_hpt


# effect sizes for diff models 

# All generations with GSM (gsm == 1)
cohens_d_all_gsm <- df_100iter %>%
  filter(gen == "all", gsm == 1) %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
cohens_d_all_gsm

# All generations without GSM (gsm == 0)
cohens_d_all_no_gsm <- df_100iter %>%
  filter(gen == "all", gsm == 0) %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
cohens_d_all_no_gsm

# F2 generation
cohens_d_F2 <- df_100iter %>%
  filter(gen == "F2") %>%
  group_by(MAF) %>%
  cohens_d(corr_iter ~ model)
cohens_d_F2

# wilcox test hpt strategy
hpt_pvals <- df %>%
  filter(
    gen == "all",
    model %in% c("RF", "GB"),
    hpt %in% c("none", "100iter", "nested"),
    MAF %in% c(0.05, 0.01, 0.005)
  ) %>%
  mutate(
    gsm_label = ifelse(gsm == 1 | gsm == "yes", "GSM (+)", "GSM (-)"),
    hpt = factor(hpt, levels = c("none", "100iter", "nested"))
  ) %>%
  group_by(gsm_label, MAF, model) %>%
  wilcox_test(
    corr_iter ~ hpt,
    p.adjust.method = "none"
  ) %>%
  add_significance() %>%
  ungroup() %>%
  arrange(gsm_label, desc(MAF), model)

hpt_pvals

################################################################################

# save r env
save.image(file = "metrics_results.RData")

