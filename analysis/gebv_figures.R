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

################################################################################
## correlation

## select MAF
Pairwise <- MAF01df[, c("GBLUP", "LASSO", "RKHS", 
                   "EGBLUP", "BRR", "BayesB",
                   "LR", "RF", "GB")]

colnames(Pairwise) <- c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB",
                        "LR", "RF", "GB")
cor_matrix <- cor(Pairwise, use = "complete.obs")
cor_matrix_ordered <- cor_matrix[new_model_order, new_model_order]

# Plot correlation heatmap
# my_colors <- colorRampPalette(c("#80cdc1", "#dfc27d"))(200)
corrplot(cor_matrix, 
         method = "square", 
         tl.cex = 1,
         tl.col = "black",  
         col.lim = c(0, 1),
         tl.pos = 'l',
         # addCoef.col = "white",       
         # number.cex = 0.8,
         mar = c(0, 0, 0, 0)) %>% 
  corrRect(c(1, 6), col = "green", lwd = 4) %>%  
  corrRect(c(7, 9), col = "red", lwd = 4)


################################################################################

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

p_overlay_status <- ggplot(MAF05df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  # geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggpar(p_overlay_status,
      palette = "startrek")

p_overlay_status <- ggplot(MAF01df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  # geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggpar(p_overlay_status,
      palette = "startrek")

p_overlay_status <- ggplot(MAF005df_long, aes(x = Value, y = Model, fill = Status_Label)) +
  geom_density_ridges(alpha=0.7) +
  theme_ridges() + 
  theme(legend.position = "top") +
  labs(
    x = "Breeding Value",
    y = "Model",
    fill = "Status") +
  theme(panel.grid.major = element_blank())+
  # geom_vline(xintercept = 0.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0.5328808, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1))
ggpar(p_overlay_status,
      palette = "startrek")

################################################################################

ggplot(MAF01df_long, aes(x = Model, y = Value, fill = Model)) +
  geom_violin() +
  scale_fill_manual(values = model_color_palette) +
  labs(x = "Model", y = "Breeding Value", fill="Model") +
  # geom_jitter(aes(color = Model, 
  #                 shape = as.factor(Status_Label)), 
  #             width = 0.2, size = 2, alpha = 0.7) + 
  scale_color_manual(values = model_color_palette) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(size = 12)) + 
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0.5, l = 0, unit = "pt")))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.75, hjust = 1, 
                                   margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")))