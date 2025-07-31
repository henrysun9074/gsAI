# ================================
# DarpaAllParallel.R - Parallel Genomic Prediction
# ================================

cat("==> R script has started and loaded ==>\n")
.libPaths("/work/tfs3/gsAI/4paul/rlib")
cat("Library path:\n")
print(.libPaths())

# Load Libraries
# library(rlang)
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

# Set working directory
setwd("/work/tfs3/gsAI/4paul")

# cat("Reading phenotype and genotype files...\n")
# phenTrain <- read_xlsx("phenTrainDarpa.xlsx")
# train_base <- read_csv("train_baseDarpa.csv")
# train_base <- column_to_rownames(train_base, var = colnames(train_base)[1])
# save(phenTrain, train_base, file = "darpa_input_data.RData")
# cat("Saved preloaded data to darpa_input_data.RData\n")

## for repeat loads
cat("Loading pre-saved data from RData file...\n")
load("darpa_input_data.RData")
cat("Loaded data from darpa_input_data.RData\n")

############################# Preprocessing and parallelization ###########################

# Split setup
n_samples <- nrow(train_base)
set.seed(123)
shuffled_indices <- sample(1:n_samples)
test_size <- 472
split_indices <- split(shuffled_indices, ceiling(seq_along(shuffled_indices) / test_size))
cat("Number of folds:", length(split_indices), "\n")

# Phenotype vector
pheno_train_vec1 <- phenTrain$Status
names(pheno_train_vec1) <- phenTrain$ID

# Parallel setup
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB per worker
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1
cat("Using", n_cores, "cores via future::multicore\n")
plan(multicore, workers = n_cores)
dir.create("logs", showWarnings = FALSE)

############################# CROSS VALIDATION ################################

cat("==> Entering main prediction loop...\n")
results <- future_map(1:10, function(i) {
  message(sprintf("Starting iteration i = %d", i))
  
  log_file <- sprintf("logs/worker_%02d.log", i)
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  on.exit({ sink(); close(log_con) }, add = TRUE)
  cat(sprintf("==> Worker %d started at %s\n", i, Sys.time()))
  
  pred_GBLUP <- vector("list", 5)
  pred_LASSO <- vector("list", 5)
  pred_RKHS <- vector("list", 5)
  pred_EGBLUP <- vector("list", 5)
  pred_BRR <- vector("list", 5)
  pred_BayesB <- vector("list", 5)
  
  for (j in 1:5) {
    cat(sprintf("-- Fold %d --\n", j))
    
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
    
    cat("Running GBLUP...\n")
    t1 <- system.time(pred_GBLUP[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "GBLUP"))))
    print(t1)
    
    cat("Running LASSO...\n")
    t2 <- system.time(pred_LASSO[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "BL"))))
    print(t2)
    
    cat("Running RKHS...\n")
    t3 <- system.time(pred_RKHS[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "RKHS"))))
    print(t3)
    
    cat("Running EGBLUP...\n")
    t4 <- system.time(pred_EGBLUP[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "EGBLUP"))))
    print(t4)
    
    cat("Running BRR...\n")
    t5 <- system.time(pred_BRR[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "BRR"))))
    print(t5)
    
    cat("Running BayesB...\n")
    t6 <- system.time(pred_BayesB[[j]] <- do.call(bwgs.predict, c(common_args, list(predict.method = "BB"))))
    print(t6)
  }
  
  cat(sprintf("==> Worker %d finished at %s\n", i, Sys.time()))
  mem <- pryr::mem_used()
  message(sprintf("Worker %d memory used: %.2f GB", i, as.numeric(mem) / 1e9))
  
  list(
    GBLUP  = do.call(rbind, pred_GBLUP),
    LASSO  = do.call(rbind, pred_LASSO),
    RKHS   = do.call(rbind, pred_RKHS),
    EGBLUP = do.call(rbind, pred_EGBLUP),
    BRR    = do.call(rbind, pred_BRR),
    BayesB = do.call(rbind, pred_BayesB)
  )
}, .options = furrr_options(seed = TRUE))
save.image("cv_parallel.RData")  # optional full workspace save


############################# GEBV Computation ################################

cat("==> Collecting model predictions from future_map results...\n")
GBLUP_all_preds2  <- lapply(results, `[[`, "GBLUP")
LASSO_all_preds2  <- lapply(results, `[[`, "LASSO")
RKHS_all_preds2   <- lapply(results, `[[`, "RKHS")
EGBLUP_all_preds2 <- lapply(results, `[[`, "EGBLUP")
BRR_all_preds2    <- lapply(results, `[[`, "BRR")
BayesB_all_preds2 <- lapply(results, `[[`, "BayesB")

cat("==> Saving cross-validation results to disk...\n")
saveRDS(list(
  GBLUP = GBLUP_all_preds2,
  LASSO = LASSO_all_preds2,
  RKHS = RKHS_all_preds2,
  EGBLUP = EGBLUP_all_preds2,
  BRR = BRR_all_preds2,
  BayesB = BayesB_all_preds2
), file = "cv_results_checkpoint.rds")

# Post-processing
cat("==> Aggregating results...\n")

average_GEBV2 <- function(preds_list2) {
  preds_matrix2 <- do.call(cbind, preds_list2)
  row_means2 <- apply(preds_matrix2, 1, mean)
  row_sd2 <- apply(preds_matrix2, 1, sd)
  list(mean = row_means2, sd = row_sd2)
}

# Summary per model
summarize_model <- function(pred_list2) {
  predicts2 <- sapply(pred_list2, function(x) x[, 1])
  stds2     <- sapply(pred_list2, function(x) x[, 2])
  data.frame(
    mean_predict2 = rowMeans(predicts2),
    mean_std2     = rowMeans(stds2)
  )
}

model_lists <- list(
  GBLUP2  = GBLUP_all_preds2,
  LASSO2  = LASSO_all_preds2,
  RKHS2   = RKHS_all_preds2,
  EGBLUP2 = EGBLUP_all_preds2,
  BRR2    = BRR_all_preds2,
  BayesB2 = BayesB_all_preds2
)

mean_results2 <- lapply(model_lists, summarize_model)

# Combine to final table
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

DARPAGEBV <- cbind(GBLUP_GEBVS, LASSO_GEBVS, RKHS_GEBVS, EGBLUP_GEBVS, BRR_GEBVS, BayesB_GEBVS)
DARPAGEBV <- data.frame(ID = rownames(DARPAGEBV), DARPAGEBV, row.names = NULL)
colnames(DARPAGEBV) <- c("ID","GBLUP_Mean", "GBLUP_SD", "LASSO_Mean", "LASSO_SD","RKHS_Mean", "RKHS_SD","EGBLUP_Mean", "EGBLUP_SD","BRR_Mean", "BRR_SD","BayesB_Mean", "BayesB_SD")

# Merge with phenotype
DARPAGEBV2 <- merge(DARPAGEBV, phenTrain[, c("ID", "Status")], by = "ID", all.x = TRUE)

# Write output
write_xlsx(DARPAGEBV2, "CrossValDarpaGebv_Parallel.xlsx")

# Correlations
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

# TODO: modify the below section so it automatically computes SD of correlation rather than
# SD of GEBVs

# Barchart
ggplot(cor_df, aes(x = Model, y = Correlation)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Correlation - SD, ymax = Correlation + SD), 
                width = 0.2, color = "black") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    y = "Correlation (± SD)",
    x = "Model"
  ) +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# Pairwise scatterplot
Pairwise <- DARPAGEBV2 %>%
  select(GBLUP_Mean, LASSO_Mean, RKHS_Mean, EGBLUP_Mean, BRR_Mean, BayesB_Mean)
pairs(Pairwise, main = "Comparison of Genomic Prediction Models for Status", cex.labels = 1.5, upper.panel = NULL, cex = 0.75)

cat("==> Script completed at", Sys.time(), "\n")
