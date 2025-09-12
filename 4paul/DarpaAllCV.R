cat("==> R script has started and loaded ==>\n")
.libPaths("/work/tfs3/gsAI/4paul/rlib")
cat("Library path:\n")
print(.libPaths())

library(rlang)
library(tidyverse)
library(rrBLUP)
library(BGLR)
# library(SNPRelate)
library(ggplot2)
# library(randomForest)
library(BWGS)
library(qqman)
library(readxl)
library(dplyr)
library(ggplot2)
library(genetics)
# library(snpStats)
library(corrplot)
library(writexl)

## modified file paths
cat("R script started\n")
setwd("/work/tfs3/gsAI/4paul")
phenTrain<-read_xlsx("phenTrainDarpa.xlsx")
train_base<-read_csv("train_baseDarpa.csv")

train_base <- column_to_rownames(train_base, var = colnames(train_base)[1])

n_samples <- nrow(train_base)
set.seed(123)
all_indices <- 1:n_samples
shuffled_indices <- sample(all_indices)
test_size<-472
split_indices <- split(shuffled_indices, ceiling(seq_along(shuffled_indices) / test_size))
length(split_indices)

GBLUP_all_preds2 <- list()
LASSO_all_preds2 <- list()
RKHS_all_preds2 <- list()
EGBLUP_all_preds2 <- list()
BRR_all_preds2 <- list()
BayesB_all_preds2 <- list()

###############################################################
cat("==> Entering main prediction loop...\n")

# Precompute once
pheno_train_vec1 <- phenTrain$Status
names(pheno_train_vec1) <- phenTrain$ID

# Initialize final output lists
GBLUP_all_preds2 <- vector("list", 10)
LASSO_all_preds2 <- vector("list", 10)
RKHS_all_preds2 <- vector("list", 10)
EGBLUP_all_preds2 <- vector("list", 10)
BRR_all_preds2 <- vector("list", 10)
BayesB_all_preds2 <- vector("list", 10)

for (i in 1:10) {
  cat("Starting iteration i =", i, "\n")
  
  # Store fold-level predictions in lists
  pred_GBLUP <- vector("list", 5)
  pred_LASSO <- vector("list", 5)
  pred_RKHS <- vector("list", 5)
  pred_EGBLUP <- vector("list", 5)
  pred_BRR <- vector("list", 5)
  pred_BayesB <- vector("list", 5)
  
  for (j in 1:5) {
    cat("  Fold j =", j, "\n")
    
    Statusgeno_target <- train_base[split_indices[[j]], ]
    Statusgeno_train <- train_base[-split_indices[[j]], ]
    Statuspheno_train_vec_j <- pheno_train_vec1[rownames(Statusgeno_train)]
    
    common_args <- list(
      geno_train = Statusgeno_train,
      pheno_train = Statuspheno_train_vec_j,
      geno_target = Statusgeno_target,
      MAXNA = 1, MAF = 0,
      geno.reduct.method = "NULL", reduct.size = "NULL",
      r2 = "NULL", pval = "NULL", MAP = "NULL",
      geno.impute.method = "MNI"
    )
    
    pred_GBLUP[[j]]  <- do.call(bwgs.predict, c(common_args, list(predict.method = "GBLUP")))
    pred_LASSO[[j]]  <- do.call(bwgs.predict, c(common_args, list(predict.method = "BL")))
    pred_RKHS[[j]]   <- do.call(bwgs.predict, c(common_args, list(predict.method = "RKHS")))
    pred_EGBLUP[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "EGBLUP")))
    pred_BRR[[j]]    <- do.call(bwgs.predict, c(common_args, list(predict.method = "BRR")))
    pred_BayesB[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "BB")))
  }
  
  # Combine predictions for current i
  GBLUP_all_preds2[[i]]  <- do.call(rbind, pred_GBLUP)
  LASSO_all_preds2[[i]]  <- do.call(rbind, pred_LASSO)
  RKHS_all_preds2[[i]]   <- do.call(rbind, pred_RKHS)
  EGBLUP_all_preds2[[i]] <- do.call(rbind, pred_EGBLUP)
  BRR_all_preds2[[i]]    <- do.call(rbind, pred_BRR)
  BayesB_all_preds2[[i]] <- do.call(rbind, pred_BayesB)
  
  cat("Completed iteration i =", i, "\n")
}

#################################################################

average_GEBV2 <- function(preds_list2) {
  preds_matrix2 <- do.call(cbind, preds_list2)
  row_means2 <- apply(preds_matrix2, 1, mean)
  row_sd2 <- apply(preds_matrix2, 1, sd)
  return(list(mean = row_means2, sd = row_sd2))
}

model_lists <- list(
  GBLUP2 = GBLUP_all_preds2,
  LASSO2 = LASSO_all_preds2,
  RKHS2 = RKHS_all_preds2,
  EGBLUP2 = EGBLUP_all_preds2,
  BRR2 = BRR_all_preds2,
  BayesB2 = BayesB_all_preds2
)

mean_results2 <- lapply(model_lists, function(pred_list2) {
  predicts2 <- sapply(pred_list2, function(x) x[, 1])
  stds2     <- sapply(pred_list2, function(x) x[, 2])
  data.frame(
    mean_predict2 = rowMeans(predicts2),
    mean_std2     = rowMeans(stds2)
  )
})

GBLUP <- sapply(GBLUP_all_preds2, function(x) x[, 1])
GBLUP_means <- rowMeans(GBLUP)
GBLUP_sds <- apply(GBLUP, 1, sd)
GBLUP_GEBVS <- cbind(mean = GBLUP_means, sd = GBLUP_sds)

LASSO <- sapply(LASSO_all_preds2, function(x) x[, 1])
LASSO_means <- rowMeans(LASSO)
LASSO_sds <- apply(LASSO, 1, sd)
LASSO_GEBVS <- cbind(mean = LASSO_means, sd = LASSO_sds)

RKHS <- sapply(RKHS_all_preds2, function(x) x[, 1])
RKHS_means <- rowMeans(RKHS)
RKHS_sds <- apply(RKHS, 1, sd)
RKHS_GEBVS <- cbind(mean = RKHS_means, sd = RKHS_sds)

EGBLUP <- sapply(EGBLUP_all_preds2, function(x) x[, 1])
EGBLUP_means <- rowMeans(EGBLUP)
EGBLUP_sds <- apply(EGBLUP, 1, sd)
EGBLUP_GEBVS <- cbind(mean = EGBLUP_means, sd = EGBLUP_sds)

BRR <- sapply(BRR_all_preds2, function(x) x[, 1])
BRR_means <- rowMeans(BRR)
BRR_sds <- apply(BRR, 1, sd)
BRR_GEBVS <- cbind(mean = BRR_means, sd = BRR_sds)

BayesB <- sapply(BayesB_all_preds2, function(x) x[, 1])
BayesB_means <- rowMeans(BayesB)
BayesB_sds <- apply(BayesB, 1, sd)
BayesB_GEBVS <- cbind(mean = BayesB_means, sd = BayesB_sds)

DARPAGEBV <- cbind(GBLUP_GEBVS, LASSO_GEBVS, RKHS_GEBVS, EGBLUP_GEBVS, BRR_GEBVS, BayesB_GEBVS)
DARPAGEBV <- data.frame(RowName = rownames(DARPAGEBV), DARPAGEBV, row.names = NULL)
colnames(DARPAGEBV) <- c("ID","GBLUP_Mean", "GBLUP_SD", "LASSO_Mean", "LASSO_SD","RKHS_Mean", "RKHS_SD","EGBLUP_Mean", "EGBLUP_SD","BRR_Mean", "BRR_SD","BayesB_Mean", "BayesB_SD")

DARPAGEBV2 <- merge(DARPAGEBV, phenTrain[, c("ID", "Status")],
                   by = "ID", all.x = TRUE)


write_xlsx(DARPAGEBV2, "CrossValDarpaGebv.xlsx")

cor(DARPAGEBV2$GBLUP_Mean, DARPAGEBV2$Status, use = "complete.obs")
cor(DARPAGEBV2$LASSO_Mean,DARPAGEBV2$Status)
cor(DARPAGEBV2$RKHS_Mean,DARPAGEBV2$Status)
cor(DARPAGEBV2$EGBLUP_Mean,DARPAGEBV2$Status)
cor(DARPAGEBV2$BRR_Mean,DARPAGEBV2$Status)
cor(DARPAGEBV2$BayesB_Mean,DARPAGEBV2$Status)

Pairwise <- cbind(GBLUP_means, LASSO_means, RKHS_means, EGBLUP_means, BRR_means, BayesB_means)
colnames(Pairwise) <- c("GBLUP", "LASSO","RKHS","EGBLUP","BRR","BayesB")
Pairwise<-as.data.frame(Pairwise)
pairs(Pairwise[1:6], main = "Comparision of Genomic Prediction Models for Status", cex.labels = 1.5, upper.panel = NULL, cex = 0.75)


################################################################################################

# run for data viz

.libPaths("/work/tfs3/gsAI/4paul/rlib")
cat("Library path:\n")
print(.libPaths())

# library(rlang)
library(tidyverse)
library(readxl)

setwd("/work/tfs3/gsAI/4paul")
load("cv_parallel.RData")

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
