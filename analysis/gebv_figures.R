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
library(cowplot)
library(corrplot)
library(rstatix)
library(tidyplots)
library(ggrepel)
library(ggstats)
library(ggalign)
library(ggridges)
library(GGally)
library(ggplotify)
library(patchwork)
library(gridGraphics)

# extra SNPs from GSM 
extra_df <- read.csv("/work/tfs3/gsAI/data/combined_gebvs_extra.csv")

# default
df <- read.csv("/work/tfs3/gsAI/data/combined_gebvs.csv")

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

MAF05df <- df[df$MAF == 0.05,]
MAF01df <- df[df$MAF == 0.01,]
MAF005df <- df[df$MAF == 0.005,]

MAF05extradf <- extra_df[extra_df$MAF == 0.05,]
MAF01extradf <- extra_df[extra_df$MAF == 0.01,]
MAF005extradf <- extra_df[extra_df$MAF == 0.005,]

################################################################################
## correlation plots for default dataset

# MAF05
Pairwise <- MAF05df[, c("GBLUP", "LASSO", "RKHS",
                   "EGBLUP", "BRR", "BayesB",
                   "LR", "RF", "GB")]
colnames(Pairwise) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")
cor_matrix <- cor(Pairwise, use = "complete.obs")
# pdf(file = "/work/tfs3/gsAI/analysis/pdfs/CorrelationMatrixMAF05.pdf", width = 7, height = 6)
corrplot(cor_matrix,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         tl.pos = 'l',
         diag = FALSE,
         type = 'upper',
         order = 'hclust',
         addCoef.col = "white",
         # number.cex = 0.8,
         mar = c(0, 0, 0, 0))
# dev.off()

# MAF01
Pairwise01 <- MAF01df[, c("GBLUP", "LASSO", "RKHS",
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
colnames(Pairwise01) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")
cor_matrix_01 <- cor(Pairwise01, use = "complete.obs")
# pdf(file = "/work/tfs3/gsAI/analysis/pdfs/CorrelationMatrixMAF01.pdf", width = 7, height = 6)
corrplot(cor_matrix_01,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
# dev.off()

## MAF005
Pairwise005 <- MAF005df[, c("GBLUP", "LASSO", "RKHS",
                          "EGBLUP", "BRR", "BayesB",
                          "LR", "RF", "GB")]

colnames(Pairwise005) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                          "LR", "RF", "GB")
cor_matrix005 <- cor(Pairwise005, use = "complete.obs")
cor_matrix_ordered005 <- cor_matrix005[new_model_order, new_model_order]
# pdf(file = "/work/tfs3/gsAI/analysis/pdfs/CorrelationMatrixMAF005.pdf", width = 7, height = 6)
corrplot(cor_matrix_ordered005,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
# dev.off()


pdf("/work/tfs3/gsAI/analysis/pdfs/CorrelationMatrixAllMAF.pdf",
    width = 12, height = 5)

par(mfrow = c(1, 3),   
    mar = c(1,1,1.5,1)) 

corrplot(cor_matrix,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.05", line = 0.5)

corrplot(cor_matrix_01,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.01", line = 0.5)

corrplot(cor_matrix_ordered005,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.005", line = 0.5)

dev.off()

################################################################################

# correlation plots for GSMs

Pairwise <- MAF05extradf[, c("GBLUP", "LASSO", "RKHS",
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
colnames(Pairwise) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")
cor_matrix_extra <- cor(Pairwise, use = "complete.obs")
corrplot(cor_matrix_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
# MAF01
Pairwise01 <- MAF01extradf[, c("GBLUP", "LASSO", "RKHS",
                          "EGBLUP", "BRR", "BayesB",
                          "LR", "RF", "GB")]
colnames(Pairwise01) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                          "LR", "RF", "GB")
cor_matrix_01_extra <- cor(Pairwise01, use = "complete.obs")
corrplot(cor_matrix_01_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))

## MAF005
Pairwise005 <- MAF005extradf[, c("GBLUP", "LASSO", "RKHS",
                            "EGBLUP", "BRR", "BayesB",
                            "LR", "RF", "GB")]

colnames(Pairwise005) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                           "LR", "RF", "GB")
cor_matrix005_extra <- cor(Pairwise005, use = "complete.obs")
cor_matrix_ordered005_extra <- cor_matrix005_extra[new_model_order, new_model_order]
corrplot(cor_matrix_ordered005_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))



pdf("/work/tfs3/gsAI/analysis/pdfs/CorrelationMatrixAllMAFImputed.pdf",
    width = 12, height = 5)

par(mfrow = c(1, 3),  
    mar = c(1,1,1.5,1)) 

corrplot(cor_matrix_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.05", line = 0.5)

corrplot(cor_matrix_01_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.01", line = 0.5)

corrplot(cor_matrix_ordered005_extra,
         method = "color",
         tl.cex = 1,
         tl.col = "black",
         col.lim = c(0, 1),
         diag = TRUE,
         type = "upper",
         order = "hclust",
         addCoef.col = "white",
         mar = c(0, 0, 0, 0))
title("MAF = 0.005", line = 0.5)

dev.off()

################################################################################
# SPLOM of model GEBVs for normal dataset

lim_min <- 0
lim_max <- 1
diag_label <- function(data, mapping, ...) {
  label_text <- rlang::as_label(mapping$x)
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, 
             size = 5, fontface = "bold") +
    theme_void()
}

# this helper function generates the scatterplot
custom_scatter <- function(data, mapping, ...) {
  x_val <- eval_data_col(data, mapping$x)
  y_val <- eval_data_col(data, mapping$y)
  correlation_rho <- round(cor(x_val, y_val, use = "complete.obs"), 3)
  rho_label <- paste0("rho*paste(' = ',", correlation_rho, ")")
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.4, size = 0.6, color = "midnightblue") +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.5) +
    scale_x_continuous(limits = c(lim_min, lim_max), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(lim_min, lim_max), breaks = seq(0, 1, 0.25)) +
    annotate("text", x = lim_min, y = lim_max, 
             label = rho_label, parse = TRUE,
             hjust = -0.1, vjust = 1.5, size = 3.5) +
    theme_pubr() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) 
}

# MAF05
Pairwise <- MAF05df[, c("GBLUP", "LASSO", "RKHS", 
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
p05 <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF05SPLOM.pdf", p05, width = 15, height = 10, units = "in")


# MAF01
Pairwise <- MAF01df[, c("GBLUP", "LASSO", "RKHS", 
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
p01 <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01SPLOM.pdf", p01, width = 15, height = 10, units = "in")


# MAF005
Pairwise <- MAF005df[, c("GBLUP", "LASSO", "RKHS", 
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
p005 <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005SPLOM.pdf", p005, width = 15, height = 10, units = "in")

################################################################################

lim_min <- 0
lim_max <- 1
diag_label <- function(data, mapping, ...) {
  label_text <- rlang::as_label(mapping$x)
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, 
             size = 5, fontface = "bold") +
    theme_void()
}

# this helper function generates the scatterplot
custom_scatter <- function(data, mapping, ...) {
  x_val <- eval_data_col(data, mapping$x)
  y_val <- eval_data_col(data, mapping$y)
  correlation_rho <- round(cor(x_val, y_val, use = "complete.obs"), 3)
  rho_label <- paste0("rho*paste(' = ',", correlation_rho, ")")
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.4, size = 0.6, color = "skyblue4") +
    geom_smooth(method = "lm", se = FALSE, color = "lightcoral", linewidth = 0.5) +
    scale_x_continuous(limits = c(lim_min, lim_max), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(lim_min, lim_max), breaks = seq(0, 1, 0.25)) +
    annotate("text", x = lim_min, y = lim_max, 
             label = rho_label, parse = TRUE,
             hjust = -0.1, vjust = 1.5, size = 3.5) +
    theme_pubr() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) 
}

# MAF05
Pairwise <- MAF05extradf[, c("GBLUP", "LASSO", "RKHS", 
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
p05_extra <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF05SPLOMImputed.pdf", p05_extra, width = 15, height = 10, units = "in")


# MAF01
Pairwise <- MAF01extradf[, c("GBLUP", "LASSO", "RKHS", 
                        "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")]
p01_extra <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01SPLOMImputed.pdf", p01_extra, width = 15, height = 10, units = "in")


# MAF005
Pairwise <- MAF005extradf[, c("GBLUP", "LASSO", "RKHS", 
                         "EGBLUP", "BRR", "BayesB",
                         "LR", "RF", "GB")]
p005_extra <- ggpairs(
  Pairwise,
  columns = 1:ncol(Pairwise),
  upper = "blank", 
  diag = list(continuous = diag_label),
  lower = list(continuous = custom_scatter),
  xlab = "Breeding Value"
) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 5))
  )
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005SPLOMImputed.pdf", p005_extra, width = 15, height = 10, units = "in")

################################################################################

# ridgeline plots for standard dataset

gebv_cols <- c("GBLUP", "EGBLUP", "BRR", "BayesB",
               "LASSO", "RKHS", "LR", "RF", "GB")
MAF05df_long <- MAF05df %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )
MAF01df_long <- MAF01df %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )
MAF005df_long <- MAF005df %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )

MAF05df_long$Model <- factor(MAF05df_long$Model, levels = new_model_order)
MAF01df_long$Model <- factor(MAF01df_long$Model, levels = new_model_order)
MAF005df_long$Model <- factor(MAF005df_long$Model, levels = new_model_order)

MAF05_p_overlay_status <- ggplot(MAF05df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
MAF05_p_overlay_status <- ggpar(MAF05_p_overlay_status,
      palette = "startrek")
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF05ridgeline.pdf", MAF05_p_overlay_status, width=8,height=7,dpi=300)

MAF01_p_overlay_status <- ggplot(MAF01df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # mean survival = 0.53
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
MAF01_p_overlay_status <- ggpar(MAF01_p_overlay_status,
      palette = "startrek")
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01ridgeline.pdf", MAF01_p_overlay_status, width=10,height=6,dpi=300)

MAF005_p_overlay_status <- ggplot(MAF005df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
MAF005_p_overlay_status <- ggpar(MAF005_p_overlay_status,
      palette = "startrek")
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005ridgeline.pdf", MAF005_p_overlay_status, width=10,height=6,dpi=300)

################################################################################

# ridgeline plots for imputed SNPs

gebv_cols <- c("GBLUP", "EGBLUP", "BRR", "BayesB",
               "LASSO", "RKHS", "LR", "RF", "GB")
MAF05extradf_long <- MAF05extradf %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )
MAF01extradf_long <- MAF01extradf %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )
MAF005extradf_long <- MAF005extradf %>%
  dplyr::select(gebv_cols, "Status") %>% # Include 'Status'
  tidyr::pivot_longer(cols = all_of(gebv_cols), 
                      names_to = "Model",
                      values_to = "Value") %>%
  dplyr::mutate(
    Value = as.numeric(Value),
    Status_Label = factor(Status, levels = c(0, 1), labels = c("Dead", "Alive"))
  )

MAF05extradf_long$Model <- factor(MAF05extradf_long$Model, levels = new_model_order)
MAF01extradf_long$Model <- factor(MAF01extradf_long$Model, levels = new_model_order)
MAF005extradf_long$Model <- factor(MAF005extradf_long$Model, levels = new_model_order)

MAF05_p_overlay_status_extra <- ggplot(MAF05extradf_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  scale_fill_manual(values = c(
    "Dead" = "palevioletred3",   
    "Alive" = "dodgerblue3"   
  )) + 
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggsave("/work/tfs3/gsAI/analysis/pdfs/ImputedMAF05ridgeline.pdf", MAF05_p_overlay_status_extra, 
       width=8,height=7,dpi=300)

MAF01_p_overlay_status_extra <- ggplot(MAF01extradf_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  scale_fill_manual(values = c(
    "Dead" = "palevioletred3",   
    "Alive" = "dodgerblue3"   
  )) +
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggsave("/work/tfs3/gsAI/analysis/pdfs/ImputedMAF01ridgeline.pdf", MAF01_p_overlay_status_extra, 
       width=10,height=6,dpi=300)

MAF005_p_overlay_status_extra <- ggplot(MAF005extradf_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  scale_fill_manual(values = c(
    "Dead" = "palevioletred3",   
    "Alive" = "dodgerblue3"   
  )) +
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggsave("/work/tfs3/gsAI/analysis/pdfs/ImputedMAF005ridgeline.pdf", MAF005_p_overlay_status_extra, width=10,height=6,dpi=300)


################################################################################

# combine multipanel ridgeline plots

# MAF 05
MAF01_p_overlay_status <- MAF01_p_overlay_status + theme(legend.position = "none")
MAF01_p_overlay_status_extra <- MAF01_p_overlay_status_extra + theme(legend.position = "none")
MAF005_p_overlay_status <- MAF005_p_overlay_status + theme(legend.position = "none")
MAF005_p_overlay_status_extra <- MAF005_p_overlay_status_extra + theme(legend.position = "none")

margin_adj <- theme(plot.margin = margin(t = 10, l = 5, r = 10, b =5))
ridgeline05plot <- plot_grid(
  MAF05_p_overlay_status + margin_adj, MAF05_p_overlay_status_extra + margin_adj,
  MAF01_p_overlay_status + margin_adj, MAF01_p_overlay_status_extra + margin_adj,
  MAF005_p_overlay_status + margin_adj, MAF005_p_overlay_status_extra + margin_adj,
  nrow = 3, ncol = 2,
  align = "hv",
  labels = c("A","B","C","D","E","F"), 
  label_size = 18
)
ggsave("/work/tfs3/gsAI/analysis/pdfs/MAFAllRidgelinesCombined.pdf",
       ridgeline05plot,
       width = 12,
       height = 11,
       units = "in",
       dpi = 300)

# MAF 01
# ridgeline01plot <- plot_grid(
#   MAF01_p_overlay_status, MAF01_p_overlay_status_extra,
#   nrow = 2,
#   labels = c("A","B"), 
#   label_size = 18
# )
# ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF01ridgelinesCombined.pdf",
#        ridgeline01plot,
#        width = 7,
#        height = 9,
#        units = "in",
#        dpi = 300)
# 
# # MAF005
# ridgeline005plot <- plot_grid(
#   MAF005_p_overlay_status, MAF005_p_overlay_status_extra,
#   nrow = 2,
#   labels = c("A","B"), 
#   label_size = 18
# )
# ggsave("/work/tfs3/gsAI/analysis/pdfs/MAF005ridgelinesCombined.pdf",
#        ridgeline005plot,
#        width = 7,
#        height = 9,
#        units = "in",
#        dpi = 300)

################################################################################

# KS tests and comparison of breeding value shifts following GSM imputation
models <- c("LR", "RF", "GB", "GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB")

df_long <- df %>% mutate(Source = "Before") %>% select(Source, MAF, all_of(models))
extra_long <- extra_df %>% mutate(Source = "After") %>% select(Source, MAF, all_of(models))

df_long <- df %>% 
  mutate(
    Source = "Before",
    MAF = as.numeric(as.character(MAF)) # as.character handles factors safely
  ) %>% 
  select(Source, MAF, all_of(models))

extra_long <- extra_df %>% 
  mutate(
    Source = "After",
    MAF = as.numeric(as.character(MAF))
  ) %>% 
  select(Source, MAF, all_of(models))

combined_data <- bind_rows(df_long, extra_long) %>%
  pivot_longer(cols = all_of(models), names_to = "Model", values_to = "GEBV")

summary_stats <- combined_data %>%
  group_by(Model, MAF, Source) %>%
  summarise(
    Mean   = mean(GEBV, na.rm = TRUE),
    Median = median(GEBV, na.rm = TRUE),
    Q10     = quantile(GEBV, 0.1, na.rm = TRUE),
    Q90     = quantile(GEBV, 0.9, na.rm = TRUE),
    IQR    = IQR(GEBV, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

quantile_diffs <- summary_stats %>%
  select(Model, MAF, Source, Q10, Q90, Mean, Median) %>%
  pivot_wider(
    names_from = Source, 
    values_from = c(Q10, Q90, Mean, Median)
  ) %>%
  mutate(
    diff_Q10 = Q10_After - Q10_Before,
    diff_Q90 = Q90_After - Q90_Before,
    diff_mean = Mean_After - Mean_Before,
    diff_median = Median_After - Median_Before
  )
print(quantile_diffs)

quantile_diffs_ordered <- quantile_diffs %>%
  arrange((abs(diff_Q10)))
print(quantile_diffs_ordered)

ks_results <- combined_data %>%
  group_by(Model, MAF) %>%
  summarise(
    ks_p_value = ks.test(
      GEBV[Source == "Before"], 
      GEBV[Source == "After"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(Significant = ifelse(ks_p_value < 0.05, "*", ""))
print(ks_results)
# write.csv(ks_results, "ksresults.csv")
# write.csv(quantile_diffs_ordered, "quantilediffs.csv")

## Shapiro Wilk test
shapiro_results <- combined_data %>%
  group_by(Model, MAF, Source) %>%
  summarise(
    shapiro_p_value = if(n() >= 3 & n() <= 5000 & var(GEBV, na.rm = TRUE) > 0) {
      shapiro.test(GEBV)$p.value
    } else if (n() > 5000 & var(GEBV, na.rm = TRUE) > 0) {
      # subsample if N > 5000, as shapiro.test fails otherwise
      shapiro.test(sample(GEBV, 5000))$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    Is_Normal = ifelse(shapiro_p_value > 0.05, "Yes", "No"),
    Significant_Departure = ifelse(shapiro_p_value < 0.05, "*", "")
  )
print(shapiro_results)