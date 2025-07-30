.libPaths("/work/tfs3/gsAI/4paul/rlib")
cat("Library path:\n")
print(.libPaths())

# library(rlang)
library(tidyverse)
library(readxl)

setwd("/work/tfs3/gsAI/4paul")
df <- read_xlsx("CrossValDarpaGebv.xlsx")
DARPAGEBV2 <- df

cat("Correlations with Status:\n")
cat("GBLUP: ", cor(DARPAGEBV2$GBLUP_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")
cat("LASSO: ", cor(DARPAGEBV2$LASSO_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")
cat("RKHS: ",  cor(DARPAGEBV2$RKHS_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")
cat("EGBLUP:", cor(DARPAGEBV2$EGBLUP_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")
cat("BRR:   ", cor(DARPAGEBV2$BRR_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")
cat("BayesB:", cor(DARPAGEBV2$BayesB_Mean, DARPAGEBV2$Status, use = "complete.obs"), "\n")

cor_df <- data.frame(
  Model = c("GBLUP", "LASSO", "RKHS", "EGBLUP", "BRR", "BayesB"),
  Correlation = c(
    cor(DARPAGEBV2$GBLUP_Mean,  DARPAGEBV2$Status, use = "complete.obs"),
    cor(DARPAGEBV2$LASSO_Mean,  DARPAGEBV2$Status, use = "complete.obs"),
    cor(DARPAGEBV2$RKHS_Mean,   DARPAGEBV2$Status, use = "complete.obs"),
    cor(DARPAGEBV2$EGBLUP_Mean, DARPAGEBV2$Status, use = "complete.obs"),
    cor(DARPAGEBV2$BRR_Mean,    DARPAGEBV2$Status, use = "complete.obs"),
    cor(DARPAGEBV2$BayesB_Mean, DARPAGEBV2$Status, use = "complete.obs")
  ),
  SD = c(
    sd(DARPAGEBV2$GBLUP_Mean,  na.rm = TRUE),
    sd(DARPAGEBV2$LASSO_Mean,  na.rm = TRUE),
    sd(DARPAGEBV2$RKHS_Mean,   na.rm = TRUE),
    sd(DARPAGEBV2$EGBLUP_Mean, na.rm = TRUE),
    sd(DARPAGEBV2$BRR_Mean,    na.rm = TRUE),
    sd(DARPAGEBV2$BayesB_Mean, na.rm = TRUE)
  )
)

# Barchart
ggplot(cor_df, aes(x = Model, y = Correlation)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    y = "Correlation (± SD)",
    x = "Model"
  ) +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

### must figure out - how to get correlation SDs rather than GEBV SDs
### have to compute raw correlation for each run rather than taking mean

# geom_errorbar(aes(ymin = Correlation - SD, ymax = Correlation + SD),
#               width = 0.2, color = "black") +