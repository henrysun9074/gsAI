args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript DarpaCV.R <phenTrain_path> <geno_path> <output_prefix> <split_fraction>")
}

# Rscript DarpaCV.R "phenTrainDarpa.xlsx" "/work/tfs3/gsAI/MAF01_DarpaQCFiltered.csv" "allMAF01QC" 0.2

phenTrain_path <- args[1]      # e.g., "phenTrainDarpa.xlsx"
geno_path      <- args[2]      # e.g., "/work/tfs3/gsAI/data/MAF0.01.csv"
output_prefix  <- args[3]      # e.g., "allMAF01QC"
split_fraction <- as.numeric(args[4]) # e.g., 0.2 for 20% test size

cat("Arguments parsed:\n")
cat("  Phenotype file:", phenTrain_path, "\n")
cat("  Genotype file: ", geno_path, "\n")
cat("  Output prefix: ", output_prefix, "\n")
cat("  Split fraction:", split_fraction, "\n")

cat("==> R script has started and loaded ==>\n")
.libPaths("/work/tfs3/gsAI/rlib")
cat("Library path:\n")
print(.libPaths())

setwd("/work/tfs3/gsAI/Rmodels")
library(future)
library(furrr)
library(pryr)
library(tidyverse)
library(rrBLUP)
library(BGLR)
library(BWGS)
library(qqman)
library(readxl)
library(genetics)
library(corrplot)
library(writexl)
library(pROC)
cat("R script started\n")

phenTrain <- read_xlsx(phenTrain_path)
train_base <- read_csv(geno_path)

train_base <- column_to_rownames(train_base, var = "ID")
train_base <- train_base[, grepl("^AX", colnames(train_base))]

n_samples <- nrow(train_base)
folds <- 5
set.seed(123)
all_indices <- 1:n_samples
shuffled_indices <- sample(all_indices)
split_indices <- split(shuffled_indices, rep(1:folds, length.out = n_samples))
length(split_indices)

pheno_train_vec1 <- phenTrain$Status
names(pheno_train_vec1) <- phenTrain$ID



#############################CROSS VALIDATION################################

cat("==> Entering main prediction loop...\n")
results <- list()

for (i in 1:10) {
  message(sprintf("Starting iteration i = %d", i))
  
  pred_GBLUP <- vector("list", 5)
  pred_LASSO <- vector("list", 5)
  pred_RKHS <- vector("list", 5)
  pred_EGBLUP <- vector("list", 5)
  pred_BRR <- vector("list", 5)
  pred_BayesB <- vector("list", 5)
  
  metrics_list <- list()
  
  run_and_eval <- function(method_name, train_args) {
    cat("Running", method_name, "...\n")
    res <- do.call(bwgs.predict, c(train_args, list(predict.method = method_name)))
    preds <- res[, 1]
    obs <- as.numeric(y_test)
    
    # ROC AUC
    auc_val <- tryCatch({
      roc_obj <- pROC::roc(obs, preds, quiet = TRUE)
      as.numeric(pROC::auc(roc_obj))
    }, error = function(e) NA)
    
    # Pearson correlation
    corr_val <- suppressWarnings(cor(preds, obs, use = "complete.obs"))
    
    # Brier score
    brier_val <- mean((obs - preds)^2, na.rm = TRUE)
    
    # Log loss (with clipping)
    eps <- 1e-15
    preds_clipped <- pmin(pmax(preds, eps), 1 - eps)
    logloss_val <- -mean(obs * log(preds_clipped) + (1 - obs) * log(1 - preds_clipped))
    
    metrics_list[[length(metrics_list)+1]] <<- data.frame(
      Iteration   = i,
      Fold        = j,
      Model       = method_name,
      AUC         = auc_val,
      Correlation = corr_val,
      Brier       = brier_val,
      LogLoss     = logloss_val
    )
    return(res)
  }
  
  for (j in 1:5) {
    cat(sprintf("-- Fold %d --\n", j))
    
    Statusgeno_target <- train_base[split_indices[[j]], ]
    Statusgeno_train <- train_base[-split_indices[[j]], ]
    Statuspheno_train_vec_j <- pheno_train_vec1[rownames(Statusgeno_train)]
    y_test <- pheno_train_vec1[rownames(Statusgeno_target)]
    
    common_args <- list(
      geno_train = Statusgeno_train,
      pheno_train = Statuspheno_train_vec_j,
      geno_target = Statusgeno_target,
      MAXNA = 1, MAF = 0,
      geno.reduct.method = "NULL", reduct.size = "NULL",
      r2 = "NULL", pval = "NULL", MAP = "NULL",
      geno.impute.method = "MNI"
    )
    
    pred_GBLUP[[j]]  <- run_and_eval("GBLUP",  common_args)
    pred_LASSO[[j]]  <- run_and_eval("BL",     common_args)
    pred_RKHS[[j]]   <- run_and_eval("RKHS",   common_args)
    pred_EGBLUP[[j]] <- run_and_eval("EGBLUP", common_args)
    pred_BRR[[j]]    <- run_and_eval("BRR",    common_args)
    pred_BayesB[[j]] <- run_and_eval("BB",     common_args)
  }
  
  results[[i]] <- list(
    preds = list(
      GBLUP  = do.call(rbind, pred_GBLUP),
      LASSO  = do.call(rbind, pred_LASSO),
      RKHS   = do.call(rbind, pred_RKHS),
      EGBLUP = do.call(rbind, pred_EGBLUP),
      BRR    = do.call(rbind, pred_BRR),
      BayesB = do.call(rbind, pred_BayesB)
    ),
    metrics = do.call(rbind, metrics_list)
  )
}

dir.create("gebvs", showWarnings = FALSE)

cat("==> Collecting and saving fold metrics...\n")
all_metrics_df <- do.call(rbind, lapply(results, function(x) x$metrics))
metrics_outfile <- sprintf("gebvs/%s_CrossValFoldMetrics_%s.csv", format(Sys.Date(), "%b%d"), output_prefix)
write.csv(all_metrics_df, metrics_outfile, row.names = FALSE)
cat("Saved fold metrics to", metrics_outfile, "\n")


#############################GEBV Computation################################

GBLUP_all_preds2  <- lapply(results, function(x) x$preds$GBLUP)
LASSO_all_preds2  <- lapply(results, function(x) x$preds$LASSO)
RKHS_all_preds2   <- lapply(results, function(x) x$preds$RKHS)
EGBLUP_all_preds2 <- lapply(results, function(x) x$preds$EGBLUP)
BRR_all_preds2    <- lapply(results, function(x) x$preds$BRR)
BayesB_all_preds2 <- lapply(results, function(x) x$preds$BayesB)

# Post-processing
cat("==> Aggregating results...\n")

# Function to combine across repetitions/folds
combine_summary <- function(preds) {
  m <- sapply(preds, function(x) x[, 1])
  s <- sapply(preds, function(x) x[, 2])
  cbind(mean = rowMeans(m), sd = apply(m, 1, sd))
}

GBLUP_GEBVS  <- combine_summary(GBLUP_all_preds2)
LASSO_GEBVS  <- combine_summary(LASSO_all_preds2)
RKHS_GEBVS   <- combine_summary(RKHS_all_preds2)
EGBLUP_GEBVS <- combine_summary(EGBLUP_all_preds2)
BRR_GEBVS    <- combine_summary(BRR_all_preds2)
BayesB_GEBVS <- combine_summary(BayesB_all_preds2)

DARPAGEBV <- cbind(
  GBLUP_GEBVS, LASSO_GEBVS, RKHS_GEBVS, EGBLUP_GEBVS, BRR_GEBVS, BayesB_GEBVS
)
DARPAGEBV <- data.frame(ID = rownames(DARPAGEBV), DARPAGEBV, row.names = NULL)
colnames(DARPAGEBV) <- c(
  "ID",
  "GBLUP_Mean", "GBLUP_SD",
  "LASSO_Mean", "LASSO_SD",
  "RKHS_Mean", "RKHS_SD",
  "EGBLUP_Mean", "EGBLUP_SD",
  "BRR_Mean", "BRR_SD",
  "BayesB_Mean", "BayesB_SD"
)

DARPAGEBV2 <- merge(DARPAGEBV, phenTrain[, c("ID", "Status")], by = "ID", all.x = TRUE)

## Modify to save correct name
outfile <- sprintf("gebvs/%s_CrossValDarpaGebv_%s.xlsx", format(Sys.Date(), "%b%d"), output_prefix)
write_xlsx(DARPAGEBV2, outfile)
cat("Saved GEBVs to", outfile, "\n")
cat("==> Script completed at", Sys.time(), "\n")
