library(rlang)
library(tidyverse)
library(rrBLUP)
library(BGLR)
library(SNPRelate)
library(ggplot2)
library(randomForest)
library(BWGS)
library(qqman)
library(readxl)
library(dplyr)
library(ggplot2)
library(genetics)
library(snpStats)
library(corrplot)
library(writexl)

phenTrain<-read_xlsx("/Users/Paul/Darpa All Gen/phenTrainDarpa.xlsx")
train_base<-read_csv("/Users/Paul/Darpa All Gen/train_baseDarpa.csv")

train_base <- column_to_rownames(train_base, var = colnames(train_base)[1])

n_samples <- nrow(train_base)
set.seed(123)
all_indices <- 1:n_samples
shuffled_indices <- sample(all_indices)
test_size<-472
split_indices <- split(shuffled_indices, ceiling(seq_along(shuffled_indices) / test_size))
length(split_indices)

#############################

GBLUP_all_preds2 <- list()
LASSO_all_preds2 <- list()
RKHS_all_preds2 <- list()
EGBLUP_all_preds2 <- list()
BRR_all_preds2 <- list()
BayesB_all_preds2 <- list()

for (i in 1:10) {
  testPREDICT_GBLUP2_i <- NULL
  testPREDICT_LASSO2_i <- NULL
  testPREDICT_RKHS2_i <- NULL
  testPREDICT_EGBLUP2_i <- NULL
  testPREDICT_BRR2_i <- NULL
  testPREDICT_BayesB2_i <- NULL
  
  pheno_train_vec1 <- phenTrain$Status
  names(pheno_train_vec1) <- phenTrain$ID
  class(pheno_train_vec1)
  
  
  for (j in 1:5) {
    Statusgeno_target <- train_base[split_indices[[j]], ]
    Statusgeno_train <- train_base[-split_indices[[j]], ]
    Statuspheno_train_vec_j <- pheno_train_vec1[rownames(Statusgeno_train)]
    
    # Predict for each method using the current split
    testPREDICT_GBLUP2_i <- rbind(testPREDICT_GBLUP2_i, 
                                  bwgs.predict(geno_train = Statusgeno_train,
                                               pheno_train = Statuspheno_train_vec_j,
                                               geno_target = Statusgeno_target,
                                               MAXNA = 1, MAF = 0,
                                               geno.reduct.method = "NULL", reduct.size = "NULL",
                                               r2 = "NULL", pval = "NULL", MAP = "NULL",
                                               geno.impute.method = "MNI", predict.method = "GBLUP"))
    
    testPREDICT_LASSO2_i <- rbind(testPREDICT_LASSO2_i, 
                                  bwgs.predict(geno_train = Statusgeno_train,
                                               pheno_train = Statuspheno_train_vec_j,
                                               geno_target = Statusgeno_target,
                                               MAXNA = 1, MAF = 0,
                                               geno.reduct.method = "NULL", reduct.size = "NULL",
                                               r2 = "NULL", pval = "NULL", MAP = "NULL",
                                               geno.impute.method = "MNI", predict.method = "BL"))
    
    testPREDICT_RKHS2_i <- rbind(testPREDICT_RKHS2_i, 
                                 bwgs.predict(geno_train = Statusgeno_train,
                                              pheno_train = Statuspheno_train_vec_j,
                                              geno_target = Statusgeno_target,
                                              MAXNA = 1, MAF = 0,
                                              geno.reduct.method = "NULL", reduct.size = "NULL",
                                              r2 = "NULL", pval = "NULL", MAP = "NULL",
                                              geno.impute.method = "MNI", predict.method = "RKHS"))
    
    testPREDICT_EGBLUP2_i <- rbind(testPREDICT_EGBLUP2_i, 
                                   bwgs.predict(geno_train = Statusgeno_train,
                                                pheno_train = Statuspheno_train_vec_j,
                                                geno_target = Statusgeno_target,
                                                MAXNA = 1, MAF = 0,
                                                geno.reduct.method = "NULL", reduct.size = "NULL",
                                                r2 = "NULL", pval = "NULL", MAP = "NULL",
                                                geno.impute.method = "MNI", predict.method = "EGBLUP"))
    
    testPREDICT_BRR2_i <- rbind(testPREDICT_BRR2_i, 
                                bwgs.predict(geno_train = Statusgeno_train,
                                             pheno_train = Statuspheno_train_vec_j,
                                             geno_target = Statusgeno_target,
                                             MAXNA = 1, MAF = 0,
                                             geno.reduct.method = "NULL", reduct.size = "NULL",
                                             r2 = "NULL", pval = "NULL", MAP = "NULL",
                                             geno.impute.method = "MNI", predict.method = "BRR"))
    
    testPREDICT_BayesB2_i <- rbind(testPREDICT_BayesB2_i, 
                                   bwgs.predict(geno_train = Statusgeno_train,
                                                pheno_train = Statuspheno_train_vec_j,
                                                geno_target = Statusgeno_target,
                                                MAXNA = 1, MAF = 0,
                                                geno.reduct.method = "NULL", reduct.size = "NULL",
                                                r2 = "NULL", pval = "NULL", MAP = "NULL",
                                                geno.impute.method = "MNI", predict.method = "BB"))
  }
  
  GBLUP_all_preds2[[i]] <- testPREDICT_GBLUP2_i
  LASSO_all_preds2[[i]] <- testPREDICT_LASSO2_i
  RKHS_all_preds2[[i]] <- testPREDICT_RKHS2_i
  EGBLUP_all_preds2[[i]] <- testPREDICT_EGBLUP2_i
  BRR_all_preds2[[i]] <- testPREDICT_BRR2_i
  BayesB_all_preds2[[i]] <- testPREDICT_BayesB2_i
}

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
