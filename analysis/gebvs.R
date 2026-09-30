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
library(corrplot)
library(FSA)
library(ggsignif)
library(ggridges)
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

setwd("/work/tfs3/gsAI/analysis/v2")

################################################################################

df <- readRDS("gebvdf.rds")

df$MAF <- factor(df$MAF, levels = unique(sort(df$MAF)))
df$Status <- factor(df$Status, levels = unique(sort(df$Status)))

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

#### check that values for non ML models are identical
# conv_df <- df %>% select(-LR, -LR_SD, -RF, -RF_SD, -GB, -GB_SD, -prefix)
# conv_df_a <- conv_df %>% filter(hpt == "none" & MAF == 0.05)
# conv_df_b <- conv_df %>% filter(hpt == "100iter" & MAF == 0.05)
# conv_df_a <- conv_df_a %>% select(-hpt)
# conv_df_b <- conv_df_b %>% select(-hpt)
# identical(conv_df_a, conv_df_b) # TRUE, good = proceed

#### dont need other hpt results for this analysis
df <- df %>% filter(hpt == "100iter") %>% select(-prefix)
unique(df$hpt)
nrow(df)

colnames(df)

################################################################################

gebv_models <- c(
  "LR", "RF", "GB",
  "GBLUP", "LASSO", "RKHS",
  "EGBLUP", "BRR", "BayesB"
)

# ks 
ks_results <- df %>%
  select(
    MAF,
    gsm,
    ID,
    Status,
    all_of(gebv_models)
  ) %>%
  pivot_longer(
    cols = all_of(gebv_models),
    names_to = "model",
    values_to = "GEBV"
  ) %>%
  group_by(model, MAF) %>%
  summarise(
    KS = list(
      ks.test(
        GEBV[gsm == "0"],
        GEBV[gsm == "1"]
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    statistic = map_dbl(KS, ~ .x$statistic),
    p_value = map_dbl(KS, ~ .x$p.value)
  ) %>%
  select(-KS)

################################################################################

# (gsm == 0)
MAF05df  <- df %>% filter(MAF == "0.05"  & (gsm == 0 | gsm == "0"))
MAF01df  <- df %>% filter(MAF == "0.01"  & (gsm == 0 | gsm == "0"))
MAF005df <- df %>% filter(MAF == "0.005" & (gsm == 0 | gsm == "0"))

# (gsm == 1)
MAF05extradf  <- df %>% filter(MAF == "0.05"  & (gsm == 1 | gsm == "1"))
MAF01extradf  <- df %>% filter(MAF == "0.01"  & (gsm == 1 | gsm == "1"))
MAF005extradf <- df %>% filter(MAF == "0.005" & (gsm == 1 | gsm == "1"))

get_cor_matrix <- function(data_subset, model_order = new_model_order) {
  mat <- cor(data_subset[, model_order], use = "complete.obs")
  mat[model_order, model_order]
}

cor_05  <- get_cor_matrix(MAF05df)
cor_01  <- get_cor_matrix(MAF01df)
cor_005 <- get_cor_matrix(MAF005df)
cor_05_extra  <- get_cor_matrix(MAF05extradf)
cor_01_extra  <- get_cor_matrix(MAF01extradf)
cor_005_extra <- get_cor_matrix(MAF005extradf)

pdf("/work/tfs3/gsAI/analysis/v2/pdfs/CorrelationMatrixAllMAF_Combined.pdf", 
    width = 14, height = 9.5)

par(mfrow = c(2, 3), mar = c(1, 1, 2.5, 1))

corrplot(cor_05, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.05", line = 0.8, cex.main = 1.3)

corrplot(cor_01, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.01", line = 0.8, cex.main = 1.3)

corrplot(cor_005, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.005", line = 0.8, cex.main = 1.3)

corrplot(cor_05_extra, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.05 (GSM)", line = 0.8, cex.main = 1.3)

corrplot(cor_01_extra, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.01 (GSM)", line = 0.8, cex.main = 1.3)

corrplot(cor_005_extra, method = "color", tl.cex = 1, tl.col = "black", col.lim = c(0, 1),
         diag = TRUE, type = "upper", order = "original", addCoef.col = "white", mar = c(0, 0, 0, 0))
title("MAF = 0.005 (GSM)", line = 0.8, cex.main = 1.3)

dev.off()

################################################################################
# Ridgeline Plots

# helper to reshape data for ridgelines
gebv_cols <- c("GBLUP", "EGBLUP", "BRR", "BayesB",
               "LASSO", "RKHS", "LR", "RF", "GB")

prepare_ridgeline_long <- function(data_subset, target_maf, target_gsm) {
  data_subset %>%
    filter(
      as.character(MAF) == as.character(target_maf),
      gsm == target_gsm
    ) %>%
    dplyr::select(all_of(gebv_cols), Status) %>%
    tidyr::pivot_longer(
      cols = all_of(gebv_cols),
      names_to = "Model",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Value = as.numeric(Value),
      Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive")),
      Model = factor(Model, levels = new_model_order)
    )
}

MAF05df_long       <- prepare_ridgeline_long(df, "0.05", 0)
MAF01df_long       <- prepare_ridgeline_long(df, "0.01", 0)
MAF005df_long      <- prepare_ridgeline_long(df, "0.005", 0)

MAF05extradf_long  <- prepare_ridgeline_long(df, "0.05", 1)
MAF01extradf_long  <- prepare_ridgeline_long(df, "0.01", 1)
MAF005extradf_long <- prepare_ridgeline_long(df, "0.005", 1)

# helper for default no GSM 
create_default_ridge <- function(df_long) {
  p <- ggplot(df_long, aes(x = Value, y = Model, fill = Status_Label)) +
    geom_density_ridges(alpha = 0.7) +
    theme_ridges() +
    theme(legend.position = "top") +
    labs(
      x = "Breeding Value",
      y = "Model",
      fill = "Status"
    ) +
    theme(panel.grid.major = element_blank()) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1))
  
  ggpar(p, palette = "startrek")
}

# helper for GSM
create_imputed_ridge <- function(df_long) {
  ggplot(df_long, aes(x = Value, y = Model, fill = Status_Label)) +
    geom_density_ridges(alpha = 0.7) +
    theme_ridges() +
    scale_fill_manual(values = c(
      "Dead" = "palevioletred3",
      "Alive" = "dodgerblue3"
    )) +
    theme(legend.position = "top") +
    labs(
      x = "Breeding Value",
      y = "Model",
      fill = "Status"
    ) +
    theme(panel.grid.major = element_blank()) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1))
}

# Generate individual plots
MAF05_p_overlay_status        <- create_default_ridge(MAF05df_long)
MAF01_p_overlay_status        <- create_default_ridge(MAF01df_long)
MAF005_p_overlay_status       <- create_default_ridge(MAF005df_long)

MAF05_p_overlay_status_extra  <- create_imputed_ridge(MAF05extradf_long)
MAF01_p_overlay_status_extra  <- create_imputed_ridge(MAF01extradf_long)
MAF005_p_overlay_status_extra <- create_imputed_ridge(MAF005extradf_long)

ridgeline05plot <- plot_grid(
  MAF05_p_overlay_status,
  MAF05_p_overlay_status_extra,
  nrow = 2,
  align = "v",
  labels = c("A", "B"),
  label_size = 18
)
ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/MAF05ridgelinesCombined.pdf",
  ridgeline05plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

ridgeline01plot <- plot_grid(
  MAF01_p_overlay_status,
  MAF01_p_overlay_status_extra,
  nrow = 2,
  align = "v",
  labels = c("A", "B"),
  label_size = 18
)
ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/MAF01ridgelinesCombined.pdf",
  ridgeline01plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

ridgeline005plot <- plot_grid(
  MAF005_p_overlay_status,
  MAF005_p_overlay_status_extra,
  nrow = 2,
  align = "v",
  labels = c("A", "B"),
  label_size = 18
)
ggsave(
  "/work/tfs3/gsAI/analysis/v2/pdfs/MAF005ridgelinesCombined.pdf",
  ridgeline005plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

################################################################################

# save image
save.image(file = "gebv_results.RData")

