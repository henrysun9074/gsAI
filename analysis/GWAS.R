rm(list = ls())
gc()                                
graphics.off()  

library(rlang)
library(asreml)
library(ASRgenomics)
library(tidyverse)
library(rrBLUP)
library(BGLR)
library(patchwork)
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
library(ASRgwas)
library(corrplot)
library(writexl)
library(impute)
library(ggbreak)

theme_set(
  theme_minimal(base_family = "Arial"))

getwd()
setwd("/Users/Paul/Darpa All Gen")

#Loading Genotypes, Phenotypes, SNP Map, and Labels
phen <- read_xlsx("/Users/Paul/Darpa All Gen/F0F1F2_phenotype.xlsx")
geno <- read_ped("/Users/Paul/Darpa All Gen/fixed_output.ped")
genotable <-read.table("/Users/Paul/Darpa All Gen/fixed_output.ped")
map <- read.table("/Users/Paul/Darpa All Gen/F0F1F2.map")
map <- map[,-3]
colnames(map) <- c("Chr","Marker","Pos")
genoID <- read_xlsx("/Users/Paul/Darpa All Gen/F0F1F2_phenotype.xlsx")
genoID <- genoID[genoID$Pop == "Training", ]
genoID <- genoID[, 1]

#Adjusting Labels and Reformating Genotypes
genoID$ID[1:1450] <- sub("_.*", "", genoID$ID[1:1450])
genotable$V1[1:1523] <- sub("_.*", "", genotable$V1[1:1523])
genotable$V1[1524:2288] <- sub("^([^_]+_[^_]+)_.*", "\\1", genotable$V1[1524:2288])
genotable$V1[2289:3058] <- sub("_.*", "", genotable$V1[2289:3058])
genotable$V1[3059:3816] <- sub("^([^_]+_[^_]+)_.*", "\\1", genotable$V1[3059:3816])
nonmatchingrows <- which(!genotable$V1 %in% genoID$ID)
p = geno$p 
n = geno$n 
geno = geno$x 
geno[geno == 2] <- NA 
geno[geno == 3] <- 2 
Geno <- matrix(geno, nrow = p, ncol = n, byrow = TRUE) 
Geno <- t(Geno) 
Geno <- Geno[-nonmatchingrows, ]
colnames(Geno) <- map$Marker 
genotable<-genotable[-nonmatchingrows,]
rownames(Geno) <- genotable$V1
Geno<-as.data.frame(Geno)
Geno<-as.matrix(Geno)

#Cleaning Phenotypes
phenTrain <- phen[phen$Pop == "Training", ]
phenTrain$ID[1:1450] <- sub("_.*", "", phenTrain$ID[1:1450])

#Assigning Generations for Nearest Neighbor Imputation to pull from
F0_ids <- phenTrain$ID[phenTrain$Generation == "F0"]
GenoF0 <- Geno[rownames(Geno) %in% F0_ids, ]
F1_ids <- phenTrain$ID[phenTrain$Generation == "F1"]
GenoF1 <- Geno[rownames(Geno) %in% F1_ids, ]
F2_ids <- phenTrain$ID[phenTrain$Generation == "F2"]
GenoF2 <- Geno[rownames(Geno) %in% F2_ids, ]
groups <- list(
  F0 = F0_ids,
  F1 = F1_ids,
  F2 = F2_ids
)

###Explorative GWAS
GenoAll<-qc.filtering(Geno, marker.callrate = 1,
                      ind.callrate = 0.1,
                      maf = 0,
                      impute=FALSE)

#Mean Numeric Imputation per Generation
GenoImputedAll <- as.matrix(GenoAll$M.clean)  
overall_col_means <- colMeans(as.matrix(GenoAll$M.clean), na.rm = TRUE)
for (grp in names(groups)) {
  rows <- groups[[grp]]
  rows_clean <- gsub("_", "", rows)
  rownames_clean <- gsub("_", "", rownames(GenoAll$M.clean))
  idx <- which(rownames_clean %in% rows_clean)
  subset_mat <- as.matrix(GenoAll$M.clean[idx, , drop = FALSE])
  col_means <- colMeans(subset_mat, na.rm = TRUE)
  for (j in seq_len(ncol(subset_mat))) {
    na_rows <- is.na(subset_mat[, j])
    subset_mat[na_rows, j] <- col_means[j]
  }
  
  GenoImputedAll[idx, ] <- subset_mat
}

#Imputation with Nearest Neighbor, k = 10
GenoImputedAll <- as.matrix(GenoAll$M.clean)  
for (grp in names(groups)) {
  rows <- groups[[grp]]
  rows_clean <- gsub("_", "", rows)
  rownames_clean <- gsub("_", "", rownames(GenoAll$M.clean))
  idx <- which(rownames_clean %in% rows_clean)
  subset_mat <- as.matrix(GenoAll$M.clean[idx, , drop = FALSE])
  imputed <- impute.knn(subset_mat)$data
  GenoImputedAll[idx, ] <- imputed
}

#More Phenotype Cleaning
Phen_to_remove<-which(!phenTrain$`ID` %in% rownames(GenoImputedAll))
phenTrain<-phenTrain[-Phen_to_remove,]

#PreGWAS, GWAS
GenoRoundedAll<-round(GenoImputedAll)

###Memory Issues###
rm(list = setdiff(ls(), c("phenTrain", "map", "GenoRoundedAll")))
###

Map <- map[map$Marker%in%colnames(GenoRoundedAll), ]
pre.gwas <- pre.gwas(
  pheno.data = phenTrain, indiv = "ID", resp = "Status", 
  geno.data = GenoRoundedAll, Q.method = "K",
  method = "VanRaden",
  maf = 0, marker.callrate = 1, ind.callrate = 1,
  impute = FALSE)
###Graphics###
variance<-pre.gwas$plot.scree 
variance+ggtitle("Genetic Explained Variance")
###

gwasS <- gwas.asreml(
  pheno.data = pre.gwas$pheno.data, resp = "Status", gen = "ID",
  Kinv = pre.gwas$Kinv, Q = pre.gwas$Q, npc = 7,family = c("gaussian"),
  geno.data = pre.gwas$geno.data, map.data = pre.gwas$map.data,
  pvalue.thr = 0.05, bonferroni = TRUE, threads = 1,
  P3D = TRUE)


#GWAS Statistics
pvals <- gwasS$gwas.all$p.value
pvals <- pvals[!is.na(pvals)]
chisq_obs <- qchisq(pvals, df = 1, lower.tail = FALSE)
lambda <- median(chisq_obs) / qchisq(0.5, df = 1, lower.tail = FALSE)
gwasS$heritability
gwasSall<-gwasS$gwas.all
gwasSall$marker <- gsub("_", "-", gwasSall$marker)
gwasSall$chrom <- Map$Chr[match(gwasSall$marker, Map$Marker)]
gwasSall$chrom <- as.numeric(as.character(gwasSall$chrom))
gwasSall <- gwasSall[order(gwasSall$chrom), ]
gwasSall <- gwasSall[gwasSall$chrom != 11, ]
gwasSalltemp <- gwasSall %>%
  group_by(chrom) %>%
  mutate(pos = rank(pos, na.last = "keep"))
gwasSalltemp <- gwasSalltemp[-((nrow(gwasSalltemp)-1):nrow(gwasSalltemp)), ]

###Graphics###
ASRgwas::qq.plot(gwas.table = gwasSalltemp) +
  geom_point(size = 3, color = "#5B9BD5") +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16))
###

###Graphics###
p <- manhattan.plot(gwas.table = gwasSalltemp, pvalue.thr = 0.05/6586, point.size = 1)
p + 
  geom_hline(yintercept = -log10(0.05/65860), color = "black", linetype = "dashed") +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))
###

###Handling Candidate Genes for Recovery [bot few GSMs can vary depending on NNI seed]
sigall<-gwasSalltemp

sigall <- sigall %>%
  filter(!is.na(`p.value`), `p.value` <= 7.59e-6) 
markers <- sigall[[1]]
markers <- markers[markers %in% colnames(Geno)] #need to reload if removed Geno for memory

missing_percent <- sapply(markers, function(m) {
  mean(is.na(Geno[, m])) * 100
})

result <- data.frame(
  Marker = markers,
  PercentMissing = missing_percent
)

good_markers <- result$Marker[result$PercentMissing < 15] 
sigallvector<-good_markers
sigall_clean <- trimws(sigallvector)
geno_names   <- trimws(colnames(Geno))
sigall_valid <- intersect(sigallvector, colnames(Geno))

#if sig_all already found
sigmat<-read_xlsx("GSMs.xlsx")
sigall_valid<-sigmat$SNP

#Building the EX Genotypes
F0_ids <- phenTrain$ID[phenTrain$Generation == "F0"]
F1_ids <- phenTrain$ID[phenTrain$Generation == "F1"]
F2_ids <- phenTrain$ID[phenTrain$Generation == "F2"]

groups <- list(
  F0 = F0_ids,
  F1 = F1_ids,
  F2 = F2_ids
)

#ImputeMNI
GenoSelectImpute <- as.matrix(Geno)
for (grp in names(groups)) {
  ids <- groups[[grp]]
  idx <- which(rownames(GenoSelectImpute) %in% ids)
  if (length(idx) > 1) {
    subset_mat <- GenoSelectImpute[idx, sigall_valid, drop = FALSE]
    col_means <- colMeans(subset_mat, na.rm = TRUE)
    for (j in seq_len(ncol(subset_mat))) {
      na_rows <- is.na(subset_mat[, j])
      subset_mat[na_rows, j] <- col_means[j]
    }
    GenoSelectImpute[idx, sigall_valid] <- subset_mat
  }
}

#Impute KNN
GenoSelectImpute <- as.matrix(Geno)
for (grp in names(groups)) {
  ids <- groups[[grp]]
  idx <- which(rownames(GenoSelectImpute) %in% ids)
  if (length(idx) > 1) {
    subset_mat <- GenoSelectImpute[idx, sigall_valid, drop = FALSE]
    imputed <- impute.knn(subset_mat)$data
    GenoSelectImpute[idx, sigall_valid] <- imputed
  }
}

GenoSelectImpute<-as.data.frame(GenoSelectImpute)
GenoSelectImpute<-as.matrix(GenoSelectImpute)
GenoSelectImpute <- round(GenoSelectImpute)

########


###GWAS for EX
rm(list = setdiff(ls(), c("GenoSelectImpute", "sigall_valid", "phenTrain", "map")))
GenoEX1<-qc.filtering(GenoSelectImpute, marker.callrate = 0.05, 
                      ind.callrate = 0.10,
                      maf = 0.01,
                      impute=FALSE)

GenoImputedEX1 <- as.matrix(GenoEX1$M.clean)  
#train_base<-GenoRoundedAll

sum(sigall_valid %in% colnames(GenoImputedEX1)) #check the recovered SNPs were passed

write.csv(GenoImputedEX1, "TrueEX0.01Final.csv", row.names = TRUE)

##############################################################
Map <- map[map$Marker%in%colnames(GenoImputedEX1), ]
pre.gwas <- pre.gwas(
  pheno.data = phenTrain, indiv = "ID", resp = "Status", 
  geno.data = GenoImputedEX1, Q.method = "K",
  method = "VanRaden",
  maf = 0, marker.callrate = 1, ind.callrate = 1,
  impute = FALSE)

###Graphics###
variance<-pre.gwas$plot.scree 
variance+ggtitle("Genetic Explained Variance")
###

gwasS <- gwas.asreml(
  pheno.data = pre.gwas$pheno.data, resp = "Status", gen = "ID",
  Kinv = pre.gwas$Kinv, Q = pre.gwas$Q, npc = 7,family = c("gaussian"),
  geno.data = pre.gwas$geno.data, map.data = pre.gwas$map.data,
  pvalue.thr = 0.05, bonferroni = TRUE, threads = 1,
  P3D = TRUE)


###

#GWAS Stats
pvals <- gwasS$gwas.all$p.value
pvals <- pvals[!is.na(pvals)]
chisq_obs <- qchisq(pvals, df = 1, lower.tail = FALSE)
lambda1 <- median(chisq_obs) / qchisq(0.5, df = 1, lower.tail = FALSE)

gwasS$heritability
gwasSall<-gwasS$gwas.all
gwasSall$marker <- gsub("_", "-", gwasSall$marker)
gwasSall$chrom <- map$Chr[match(gwasSall$marker,map$Marker)]
gwasSall$chrom <- as.numeric(as.character(gwasSall$chrom))
gwasSall <- gwasSall[order(gwasSall$chrom), ]
gwasSall <- gwasSall[gwasSall$chrom != 11, ]
gwasSalltemp <- gwasSall %>%
  group_by(chrom) %>%
  mutate(pos = rank(pos, na.last = "keep"))
gwasSalltemp <- gwasSalltemp[-((nrow(gwasSalltemp)-1):nrow(gwasSalltemp)), ]
gwasSalltemp$pos <- Map$Pos[match(gwasSalltemp$marker, Map$Marker)]


chr_colours <- c("red", "blue")

gwasSalltemp <- gwasSalltemp %>%
  mutate(point_col = chr_colours[(chrom %% 2) + 1])
###Graphics###
p <- manhattan.plot(gwas.table = gwasSalltemp, pvalue.thr = NULL, point.size = 0, point.alpha = 0.8)

p_manhattan <- p + 
  geom_hline(yintercept = -log10(0.05/48166), color = "black", linetype = "solid") +
  geom_hline(yintercept = -log10(0.05/4816), color = "black", linetype = "dashed") +
  labs(x = NULL, y = expression(-log[10](p-value))) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank()) +
  scale_colour_manual(values = rep(c("black", "gray60"), 5)) +
  geom_point(size = 3, alpha = 0.8, shape = 16)+
  scale_y_continuous(limits = c(0, 44)) +
  scale_y_break(c(20, 42), scales = 0.1, space = 0, ticklabels = c(42, 43, 44), symbol = 'slash')

###Graphics###
p_qq <- ASRgwas::qq.plot(gwas.table = gwasSalltemp) +
  geom_point(size = 3, color = "black") +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        text = element_text(family = "Arial"))
###
patchworkplot <- p_qq | p_manhattan
patchworkplot +
  plot_annotation(title = "Inclusive Expanded QC: QQ and Manhattan Plots", theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Arial")))
