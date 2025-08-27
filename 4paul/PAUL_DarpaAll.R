rm(list = ls())

library(rlang)
# library(asreml)
library(ASRgenomics)
library(tidyverse)
library(rrBLUP)
library(BGLR)
# library(SNPRelate)
# library(BWGS)
# library(qqman)
library(readxl)
library(dplyr)
library(ggplot2)
library(genetics)
# library(snpStats)
# library(ASRgwas)
library(corrplot)
library(writexl)

theme_set(
  theme_minimal(base_family = "Arial"))

getwd()
setwd("/Users/henrysun_1/Desktop/Duke/PhD 2025-2026/genomic selection/DARPA data 3 generations combined")
phen <- read_xlsx("/Users/henrysun_1/Desktop/Duke/PhD 2025-2026/genomic selection/DARPA data 3 generations combined/F0F1F2_phenotype.xlsx")
filter_phen <- read_csv("/Users/henrysun_1/Desktop/Duke/PhD 2025-2026/genomic selection/DARPA data 3 generations combined/FilteredPhenotypeData.csv")

genotable <- read.table("F0F1F2.ped")
n <- nrow(genotable)
meta <- data.frame(
 IID = genotable[[1]],            # Individual ID (assumed column 1)
 FID = paste0("FAM", seq_len(n)), # Family ID
 PID = 0,
 MID = 0,
 SEX = 0,
 PHENO = -9)

geno_only <- genotable[, -1]
ped_fixed_no273 <- cbind(meta, geno_only)

write.table(ped_fixed_no273, "fixed_output.ped", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
file.exists("fixed_output.ped")

#geno <- read_ped("/Users/Paul/Darpa All Gen/F0F1F2.ped")
#geno <- read_ped("/Users/Paul/Darpa All Gen/F0F1F2_no273.ped")
geno <- read_ped("F0F1F2.ped")
genotable <-read.table("F0F1F2.ped")


map <- read.table("/Users/Paul/Darpa All Gen/F0F1F2.map")
map <- map[,-3]
colnames(map) <- c("Chr","Marker","Pos")


genoID <- read_xlsx("/Users/Paul/Darpa All Gen/F0F1F2_phenotype.xlsx")
genoID <- genoID[genoID$Pop == "Training", ]
genoID <- genoID[, 1]
genoID$ID <- sub("_.*", "", genoID$ID)

genotable$V1 <- sub("_.*", "", genotable$V1)
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

write_xlsx(Geno, "DarpaGeno.xlsx")
library(readr)
write_csv(as.data.frame(Geno), "DarpaGeno.csv")
test<-as.data.frame(Geno)

#############

MaxQC <- qc.filtering(Geno, marker.callrate = 0.05, 
                      ind.callrate = 0.10,
                      maf = 0.05,
                      heterozygosity = 0.9,
                      impute=TRUE)

MinQC <- qc.filtering(Geno, marker.callrate = 0.05, 
                      ind.callrate = 0.10,
                      maf = 0.05,
                      heterozygosity = 0.1,
                      impute=TRUE)

MinCol <- colnames(MinQC$M.clean)
MaxQC$M.clean <- MaxQC$M.clean[, !colnames(MaxQC$M.clean) %in% MinCol]
MaxMinQC<-MaxQC$M.clean #2856 ID and 36848 SNPs

####### Phenotype Cleaning

phenTrain <- phen[phen$Pop == "Training", ]
phenTrain$ID <- sub("_.*", "", phenTrain$ID)

Phen_to_remove<-which(!phenTrain$`ID` %in% rownames(MaxMinQC))
phenTrain<-phenTrain[-Phen_to_remove,]

######## GWAS

MaxMinQC <- round(MaxMinQC, digits=0)
train_base <- MaxMinQC
train_base <- train_base[complete.cases(train_base), ]
dim(train_base)

Map <- map[map$Marker%in%colnames(MaxMinQC), ]
Map <- Map[, c(2, 1, 3:ncol(MaxMinQC))]
colnames(Map) <- c("marker","chrom","pos")

phenTrain <- phenTrain[ , -5]
pheno.data = phenTrain
sum(is.na(phenTrain))

phenTrain[[1]] <- seq_len(nrow(phenTrain))
rownames(train_base) <- seq_len(nrow(train_base))


Base.pre.gwas <- pre.gwas(
  pheno.data = pheno.data, indiv = "ID", resp = "Status", 
  geno.data = train_base, Q.method = "K",
  method = "VanRaden",
  maf = 0.05, marker.callrate = 0.05, ind.callrate = 0.05,
  impute = FALSE)

variance<-Base.pre.gwas$plot.scree 
variance+ggtitle("Status Base QC Explained Variance")

gwasS <- gwas.asreml(
  pheno.data = Base.pre.gwas$pheno.data, resp = "Status", gen = "ID",
  Kinv = Base.pre.gwas$Kinv, Q = Base.pre.gwas$Q, npc = 7,
  geno.data = Base.pre.gwas$geno.data, map.data = Base.pre.gwas$map.data,
  pvalue.thr = 0.0005, bonferroni = FALSE, 
  P3D = TRUE)############### 48 markers

gwasS$heritability  ################################# vc =  SE =  // pev = 
gwasSall<-gwasS$gwas.all
gwasSall$marker <- gsub("_", "-", gwasSall$marker)
gwasSall$chrom <- Map$marker[match(gwasSall$marker, Map$chrom)]
gwasSall$chrom <- as.numeric(as.character(gwasSall$chrom))
gwasSall <- gwasSall[order(gwasSall$chrom), ]
gwasSalltemp <- gwasSall %>%
  group_by(chrom) %>%
  mutate(pos = rank(pos, na.last = "keep"))
manhattan.plot(gwas.table = gwasSalltemp, pvalue.thr = 0.05/36000, point.size = 1)
manhattan(gwasSalltemp, chr = "chrom", bp = "pos", p = "p.value", snp = "marker")


