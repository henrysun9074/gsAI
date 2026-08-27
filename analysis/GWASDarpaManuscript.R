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
library(ggrepel)


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

###Explorative GWAS####
GenoAll<-qc.filtering(Geno, marker.callrate = 0.2,
                      ind.callrate = 0.1,
                      maf = 0,
                      impute=FALSE)
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

impute_geno_knn <- function(geno, k = 10, maxp = 5000) {
  
  geno_t <- t(as.matrix(geno))  
  all_na <- rowSums(!is.na(geno_t)) == 0
  invariant <- apply(geno_t, 1, function(x) {
    ux <- unique(x[!is.na(x)])
    length(ux) <= 1
  })
  keep <- !(all_na | invariant)
  message("Keeping ", sum(keep), " markers; temporarily removing ", sum(!keep))
  geno_t_imp <- impute.knn(
    geno_t[keep, , drop = FALSE],
    k = k,
    maxp = maxp
  )$data
  geno_t[keep, ] <- geno_t_imp
  for (r in which(!keep)) {
    obs <- geno_t[r, !is.na(geno_t[r, ])]
    if (length(obs) > 0) {
      geno_t[r, is.na(geno_t[r, ])] <- obs[1]
    }
  }
  t(geno_t)
}
GenoF0_imp <- impute_geno_knn(GenoF0, k = 10, maxp = 5000)
GenoF1_imp <- impute_geno_knn(GenoF1, k = 10, maxp = 5000)
GenoF2_imp <- impute_geno_knn(GenoF2, k = 10, maxp = 5000)
GenoImputedAll <- rbind(GenoF0_imp, GenoF1_imp, GenoF2_imp)
GenoAll_t <- t(GenoImputedAll)
GenoAll_t<-impute.knn(GenoAll_t, k =10)$data
GenoImputedAll <- t(GenoAll_t)

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
  pvalue.thr = 0.05, bonferroni = TRUE, threads = 1, workspace = "4Gb",
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

sigall <- sigall[sigall$marker %in% sigall_valid, ]

#if sig_all already found
sigmat<-read_xlsx("GSMs.xlsx")


sigall_valid<-sigmat$marker

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
GenoSelectImpute <- t(as.matrix(Geno))  

for (grp in names(groups)) {
  ids <- groups[[grp]]
  ids_clean <- gsub("_", "", ids)
  colnames_clean <- gsub("_", "", colnames(GenoSelectImpute))
  idx <- which(colnames_clean %in% ids_clean)
  if (length(idx) > 1) {
    subset_mat <- GenoSelectImpute[sigall_valid, idx, drop = FALSE]  
    imputed <- impute.knn(subset_mat)$data
    GenoSelectImpute[sigall_valid, idx] <- imputed
  }
}

GenoSelectImpute <- t(GenoSelectImpute)  
GenoSelectImpute<-as.data.frame(GenoSelectImpute)
GenoSelectImpute<-as.matrix(GenoSelectImpute)
GenoSelectImpute <- round(GenoSelectImpute)

########

###GWAS for EX####
rm(list = setdiff(ls(), c("GenoSelectImpute", "sigall_valid", "phenTrain", "map")))
GenoEX1<-qc.filtering(GenoSelectImpute, marker.callrate = 0.05, 
                      ind.callrate = 0.10,
                      maf = 0.01,
                      impute=FALSE)

GenoImputedEX1 <- as.matrix(GenoEX1$M.clean)  
#train_base<-GenoRoundedAll

sum(sigall_valid %in% colnames(GenoImputedEX1)) #check the recovered SNPs were passed

write.csv(GenoImputedEX1, "SecondTrueEX0.01Final.csv", row.names = TRUE)

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

gwasSalltemp<-read_xlsx("EX0.01NNFinalPls.xlsx")

chr_colours <- c("red", "blue")

gwasSalltemp <- gwasSalltemp %>%
  mutate(point_col = chr_colours[(chrom %% 2) + 1])
###Graphics###
gwasSalltemp$marker <- paste0(gwasSalltemp$chrom, ":", gwasSalltemp$pos)

tag.table <- data.frame(
  marker = gwasSalltemp$marker[!is.na(gwasSalltemp$Annotation)],
  tag    = gwasSalltemp$Annotation[!is.na(gwasSalltemp$Annotation)]
)
tag.table$nudge.x <- -10000   
tag.table$nudge.y <- 0.5

tag.table$nudge.x[9] <- -500000
tag.table$nudge.x[10] <- 10000000
tag.table$nudge.y[10] <- 1.2

p <- manhattan.plot(gwas.table = gwasSalltemp, pvalue.thr = NULL, point.size = 0, point.alpha = 0.8, tag.table = tag.table, tag.repel = FALSE)

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
  geom_point(size = 3, alpha = 0.8, shape = 16)+
  scale_y_continuous(limits = c(0, 44)) +
  scale_y_break(c(20, 42), scales = 0.1, space = 0, ticklabels = c(42, 43, 44), symbol = 'slash')

###Graphics###
p_qq <- ASRgwas::qq.plot(gwas.table = gwasSalltemp) +
  geom_point(size = 3, color = "purple4") +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        text = element_text(family = "Helvetica"))
###
patchworkplot <- p_qq | p_manhattan
patchworkplot +
  plot_annotation(title = "Inclusive Expanded QC: QQ and Manhattan Plots", theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Helvetica")))


##############

##############Getting MAF Deltas#############
Graphs <- read_xlsx("/Users/Paul/Darpa All Gen/AllMarkersGeneSearch.xlsx")
DesiredMarkers<-read_xlsx("GSMs.xlsx")
#DRMS
DesiredMarkers<-c("AX-563280343", "AX-564181349", "AX-564439151", "AX-563911687",
                  "AX-570458254", "AX-576810971", "AX-567945231", "AX-563866890",
                  "AX-563859187", "AX-570229026", "AX-574123461", "AX-574557543",
                  "AX-563453997", "AX-563703370", "AX-562975825", "AX-570890079",
                  "AX-576402328", "AX-564118820", "AX-564085983", "AX-563298068",
                  "AX-574116952", "AX-568235396", "AX-576902130", "AX-563318564",
                  "AX-567276575", "AX-564193914", "AX-574389243", "AX-570642792",
                  "AX-567116016", "AX-576895722", "AX-570332532", "AX-563564797",
                  "AX-574408724", "AX-567965988", "AX-576610402", "AX-567679698",
                  "AX-570494413", "AX-571052524", "AX-574115859", "AX-575695122",
                  "AX-575122210", "AX-574092264", "AX-575653763", "AX-576896226",
                  "AX-570465218", "AX-564057487", "AX-563408780", "AX-575684868")
DesiredMarkers<-as.vector(DesiredMarkers)

Graphs <- Graphs[Graphs$marker %in% DesiredMarkers, ]
colnames(Graphs) <- gsub("_", "-", colnames(Graphs))

Graphs$Delta<-Graphs$`All-MAF-Dead`-Graphs$`All-MAF-Live`

Graphs$`Delta-F0` <- Graphs$`F0-MAF-Dead` - Graphs$`F0-MAF-Live`
Graphs$`Delta-F1` <- Graphs$`F1-MAF-Dead` - Graphs$`F1-MAF-Live`
Graphs$`Delta-F2` <- Graphs$`F2-MAF-Dead` - Graphs$`F2-MAF-Live`

Graphs$`Avg-Delta` <- rowMeans(Graphs[, c("Delta-F0","Delta-F1","Delta-F2")], na.rm = TRUE)
`marker-summary` <- Graphs[, c("marker", "Avg-Delta")]

######Graphics#########
####For each individual marker#######
for (m in unique(Graphs$marker)) {
  p <- Graphs %>%
    filter(marker == m) %>%
    pivot_longer(
      cols = c(`F0-MAF-Live`, `F0-MAF-Dead`, 
               `F1-MAF-Live`, `F1-MAF-Dead`, 
               `F2-MAF-Live`, `F2-MAF-Dead`),
      names_to = "Group",
      values_to = "MAF"
    ) %>%
    mutate(
      Generation = case_when(
        grepl("F0", Group) ~ "F0",
        grepl("F1", Group) ~ "F1",
        grepl("F2", Group) ~ "F2"
      ),
      Status = case_when(
        grepl("Live", Group) ~ "Live",
        grepl("Dead", Group) ~ "Dead"
      ),
      Generation = factor(Generation, levels = c("F0", "F1", "F2")),
      Status = factor(Status, levels = c("Live", "Dead"))
    ) %>%
    ggplot(aes(x = Generation, y = MAF, fill = Status)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Live" = "firebrick", "Dead" = "grey70")) +
    labs(title = paste0(m, ": \u0394 MAF Across Generations")) +
    theme(
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      plot.title = element_text(size = 20, family = "Helvetica", hjust = 0.5)
    )
  
  print(p)
}

#######For all markers together##########

fgen<-p


fgenp<-Graphs %>%
  pivot_longer(
    cols = c(`F0-MAF-Live`, `F0-MAF-Dead`, 
             `F1-MAF-Live`, `F1-MAF-Dead`, 
             `F2-MAF-Live`, `F2-MAF-Dead`),
    names_to = "Group",
    values_to = "MAF"
  ) %>%
  mutate(
    Generation = case_when(
      grepl("F0", Group) ~ "F0",
      grepl("F1", Group) ~ "F1",
      grepl("F2", Group) ~ "F2"
    ),
    Status = case_when(
      grepl("Live", Group) ~ "Live",
      grepl("Dead", Group) ~ "Dead"
    ),
    Generation = factor(Generation, levels = c("F0", "F1", "F2")),
    Status = factor(Status, levels = c("Live", "Dead"))
  ) %>%
  group_by(Generation, Status) %>%
  summarise(MAF = mean(MAF, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Generation, y = MAF, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Live" = "firebrick", "Dead" = "grey")) +
  labs(title = "Average \u0394 MAF Across Generations") +
  theme(
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    plot.title = element_text(size = 24, family = "Helvetica", hjust = 0.5),
    legend.position = "none"
  )

genop <- Graphs %>%
  filter(marker %in% DesiredMarkers) %>%
  pivot_longer(
    cols = c(`All-AA-Dead`, `All-AA-Live`,
             `All-AB-Dead`, `All-AB-Live`,
             `All-BB-Dead`, `All-BB-Live`),
    names_to = "Group",
    values_to = "Count"
  ) %>%
  mutate(
    Genotype = case_when(
      grepl("AA", Group) ~ "AA",
      grepl("AB", Group) ~ "AB",
      grepl("BB", Group) ~ "BB"
    ),
    
    Status = case_when(
      grepl("Dead", Group) ~ "Live",
      grepl("Live", Group) ~ "Dead"
    )
  ) %>%
  group_by(marker, Status) %>%
  mutate(
    Percent = Count / sum(Count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  mutate(
    Genotype = factor(Genotype, levels = c("AA", "AB", "BB")),
    Status = factor(Status, levels = c("Live", "Dead"))
  ) %>%
  
  group_by(marker, Genotype) %>%
  mutate(overall_freq = mean(Percent, na.rm = TRUE)) %>%
  ungroup() %>%
  
  group_by(marker) %>%
  mutate(
    AA_freq = mean(overall_freq[Genotype == "AA"], na.rm = TRUE),
    BB_freq = mean(overall_freq[Genotype == "BB"], na.rm = TRUE),
    
    Genotype = case_when(
      AA_freq < BB_freq & Genotype == "AA" ~ "BB",
      AA_freq < BB_freq & Genotype == "BB" ~ "AA",
      TRUE ~ as.character(Genotype)
    ),
    
    Genotype = factor(Genotype, levels = c("AA", "AB", "BB"))
  ) %>%
  ungroup() %>%
  
  group_by(Genotype, Status) %>%
  summarise(mean_Percent = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
  
  ggplot(aes(x = Genotype, y = mean_Percent, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Live" = "firebrick", "Dead" = "grey")) +
  labs(y = "Frequency", x = "Genotype",
       title = "Average Genotype Distribution") +
  theme(
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    plot.title = element_text(size = 24, family = "Helvetica", hjust = 0.5),
    legend.text = element_text(size = 20),   
    legend.title = element_text(size = 22),
    legend.key.size = unit(1.5, "cm")
  )

fgenp | genop
wrap_plots(fgenp, genop) +
  plot_annotation(title = "All GSMs and NonGSMs",
                  theme = theme(plot.title = element_text(size = 24,
                                                          face = "bold",
                                                          hjust = 0.5,
                                                          family = "Helvetica")))

##### Mucin 5-ac specific####
marker_id <- "AX-563280343"
p_maf <- Graphs %>%
  filter(marker == marker_id) %>%
  pivot_longer(
    cols = c(`F0-MAF-Live`, `F0-MAF-Dead`, 
             `F1-MAF-Live`, `F1-MAF-Dead`, 
             `F2-MAF-Live`, `F2-MAF-Dead`),
    names_to = "Group",
    values_to = "MAF"
  ) %>%
  mutate(
    Generation = case_when(
      grepl("F0", Group) ~ "F0",
      grepl("F1", Group) ~ "F1",
      grepl("F2", Group) ~ "F2"
    ),
    Status = case_when(
      grepl("Dead", Group) ~ "Dead",
      grepl("Live", Group) ~ "Live"
    ),
    Generation = factor(Generation, levels = c("F0", "F1", "F2")),
    Status = factor(Status, levels = c("Live", "Dead"))
  ) %>%
  ggplot(aes(x = Generation, y = MAF, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Live" = "firebrick", "Dead" = "grey")) +
  labs(y = "MAF", x = "Generation", fill = "Status", title = "\u0394 MAF Across Generations") +
  theme(
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    plot.title = element_text(size = 24, hjust = 0.5)
  )

p_geno <- Graphs %>%
  filter(marker == marker_id) %>%
  pivot_longer(
    cols = c(`All-AA-Dead`, `All-AA-Live`,
             `All-AB-Dead`, `All-AB-Live`,
             `All-BB-Dead`, `All-BB-Live`),
    names_to = "Group",
    values_to = "Count"
  ) %>%
  mutate(
    Genotype = case_when(
      grepl("AA", Group) ~ "AA",
      grepl("AB", Group) ~ "AB",
      grepl("BB", Group) ~ "BB"
    ),
    
    Status = case_when(
      grepl("Dead", Group) ~ "Live",
      grepl("Live", Group) ~ "Dead"
    )
  ) %>%
  group_by(marker, Status) %>%
  mutate(
    Percent = Count / sum(Count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  
  mutate(
    Genotype = factor(Genotype, levels = c("AA", "AB", "BB")),
    Status = factor(Status, levels = c("Live", "Dead"))
  ) %>%
  
  group_by(marker, Genotype) %>%
  mutate(overall_freq = mean(Percent, na.rm = TRUE)) %>%
  ungroup() %>%
  
  group_by(marker) %>%
  mutate(
    AA_freq = mean(overall_freq[Genotype == "AA"], na.rm = TRUE),
    BB_freq = mean(overall_freq[Genotype == "BB"], na.rm = TRUE),
    
    Genotype = case_when(
      AA_freq < BB_freq & Genotype == "AA" ~ "BB",
      AA_freq < BB_freq & Genotype == "BB" ~ "AA",
      TRUE ~ as.character(Genotype)
    ),
    
    Genotype = factor(Genotype, levels = c("AA", "AB", "BB"))
  ) %>%
  ungroup() %>%
  
  group_by(Genotype, Status) %>%
  summarise(mean_Percent = mean(Percent, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_Percent = pmin(mean_Percent, 100)) %>%
  
  ggplot(aes(x = Genotype, y = mean_Percent, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Live" = "firebrick", "Dead" = "grey")) +
  labs(y = "Frequency", x = "Genotype",
       title = "Genotype Distribution") +
  theme(
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    plot.title = element_text(size = 24, family = "Helvetica", hjust = 0.5),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(1.5, "cm")
  )
print(p_maf)
print(p_geno)


p_maf | p_geno
wrap_plots(p_maf, p_geno) +
  plot_annotation(title = "AX-563280343 / Mucin-5AC-Like",
                  theme = theme(plot.title = element_text(size = 24,
                                                          face = "bold",
                                                          hjust = 0.5,
                                                          family = "Helvetica")))

##############

########Annotation (much was done manually in native GUI, but redone in iterations below#####
setwd("/Users/Paul/Darpa All Gen")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("GenomicRanges", "rtracklayer", "IRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

for (pkg in c("dplyr", "readr", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(dplyr)
library(readr)
library(stringr)
library(readxl)

df<-read_xlsx("GSMs.xlsx")
map<-map <- read.table("/Users/Paul/Darpa All Gen/F0F1F2.map")
map <- map[,-3]
colnames(map) <- c("Chr","Marker","Pos")

library(dplyr)

df <- df %>%
  left_join(
    map %>%
      transmute(
        marker = as.character(Marker),
        new_pos = Pos
      ),
    by = "marker"
  ) %>%
  mutate(
    pos = coalesce(new_pos, pos)
  ) %>%
  select(-new_pos)

stopifnot(all(c("marker", "chrom", "pos") %in% names(df)))

df_markers <- df %>%
  select(marker, chrom, pos) %>%
  filter(!is.na(marker), !is.na(chrom), !is.na(pos)) %>%
  mutate(
    marker = as.character(marker),
    chrom = as.character(chrom),
    pos = as.integer(pos)
  ) %>%
  distinct(marker, .keep_all = TRUE)

dup_check <- df %>%
  select(marker, chrom, pos) %>%
  filter(!is.na(marker), !is.na(chrom), !is.na(pos)) %>%
  distinct() %>%
  count(marker) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  warning("Some marker IDs have multiple chrom/pos combinations. Keeping first occurrence per marker.")
}

###Downloading 2017 Assembly###

gff_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/",
  "GCF_002022765.2_C_virginica-3.0/",
  "GCF_002022765.2_C_virginica-3.0_genomic.gff.gz"
)

gff_file <- "GCF_002022765.2_C_virginica-3.0_genomic.gff.gz"

if (!file.exists(gff_file)) {
  message("Downloading GCF_002022765.2 annotation GFF...")
  download.file(gff_url, destfile = gff_file, mode = "wb")
}

###Getting Annotations###

gff <- import("genomic.gff")
genes <- gff[gff$type == "gene"]
gene_df <- as.data.frame(genes)
safe_col <- function(x, nm) {
  if (nm %in% names(x)) as.character(x[[nm]]) else rep(NA_character_, nrow(x))
}
gene_annot <- tibble(
  seqid = as.character(seqnames(genes)),
  gene_start = start(genes),
  gene_end = end(genes),
  gene_strand = as.character(strand(genes)),
  gene_id_raw = safe_col(gene_df, "ID"),
  gene_name = dplyr::coalesce(
    safe_col(gene_df, "Name"),
    safe_col(gene_df, "gene"),
    safe_col(gene_df, "locus_tag")
  ),
  gene_biotype = dplyr::coalesce(
    safe_col(gene_df, "gene_biotype"),
    safe_col(gene_df, "gbkey")
  ),
  description = dplyr::coalesce(
    safe_col(gene_df, "description"),
    safe_col(gene_df, "product")
  ),
  dbxref = safe_col(gene_df, "Dbxref")
) %>%
  mutate(
    gene_id = str_remove(gene_id_raw, "^gene-"),
    ncbi_geneid = str_extract(dbxref, "GeneID:[0-9]+"),
    ncbi_geneid = str_remove(ncbi_geneid, "GeneID:")
  )

###Renaming Chromosomes to match assembly####

gff_seqids <- unique(gene_annot$seqid)
df_chroms <- unique(df_markers$chrom)

chr_map <- data.frame(
  chrom = as.character(1:10),
  seqid = c(
    "NC_035780.1",
    "NC_035781.1",
    "NC_035782.1",
    "NC_035783.1",
    "NC_035784.1",
    "NC_035785.1",
    "NC_035786.1",
    "NC_035787.1",
    "NC_035788.1",
    "NC_035789.1"
  )
)

df_markers2 <- df_markers %>%
  left_join(chr_map, by = "chrom")

df_markers2 <- df_markers2 %>%
  filter(!is.na(seqid))

###Build marker GRanges and find nearest gene ###
marker_gr <- GRanges(
  seqnames = df_markers2$seqid,
  ranges = IRanges(start = df_markers2$pos, end = df_markers2$pos),
  marker = df_markers2$marker,
  chrom_original = df_markers2$chrom,
  pos = df_markers2$pos
)

gene_gr <- GRanges(
  seqnames = gene_annot$seqid,
  ranges = IRanges(start = gene_annot$gene_start, end = gene_annot$gene_end),
  strand = gene_annot$gene_strand
)

nearest_idx <- nearest(marker_gr, gene_gr, ignore.strand = TRUE)

out <- tibble(
  marker = mcols(marker_gr)$marker,
  chrom = mcols(marker_gr)$chrom_original,
  pos = mcols(marker_gr)$pos,
  seqid = as.character(seqnames(marker_gr)),
  nearest_gene_index = nearest_idx
) %>%
  mutate(
    nearest_gene_id = gene_annot$gene_id[nearest_gene_index],
    nearest_gene_name = gene_annot$gene_name[nearest_gene_index],
    nearest_ncbi_geneid = gene_annot$ncbi_geneid[nearest_gene_index],
    nearest_gene_biotype = gene_annot$gene_biotype[nearest_gene_index],
    nearest_description = gene_annot$description[nearest_gene_index],
    nearest_gene_start = gene_annot$gene_start[nearest_gene_index],
    nearest_gene_end = gene_annot$gene_end[nearest_gene_index],
    nearest_gene_strand = gene_annot$gene_strand[nearest_gene_index],
    
    distance_to_gene = case_when(
      pos >= nearest_gene_start & pos <= nearest_gene_end ~ 0L,
      pos < nearest_gene_start ~ nearest_gene_start - pos,
      pos > nearest_gene_end ~ pos - nearest_gene_end,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(
    marker, chrom, pos, seqid,
    nearest_gene_id,
    nearest_gene_name,
    nearest_ncbi_geneid,
    nearest_description,
    nearest_gene_biotype,
    nearest_gene_start,
    nearest_gene_end,
    nearest_gene_strand,
    distance_to_gene
  )

###Making Df###

df_nearest_gene_annotated <- df %>%
  left_join(out, by = "marker")

library(dplyr)
library(readr)
library(httr)
library(stringr)

gene_ids <- df_nearest_gene_annotated %>%
  pull(nearest_gene_id) %>%
  unique() %>%
  na.omit() %>%
  as.character()
gene_ids <- gsub("^gene-", "", gene_ids)

head(gene_ids)
length(gene_ids)

lookup_uniprot_gene <- function(gene_id, organism_id = 6565) {
  query <- paste0(
    "(gene_exact:", gene_id, ") AND (organism_id:", organism_id, ")"
  )
  
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=", URLencode(query, reserved = TRUE),
    "&fields=accession,reviewed,protein_name,gene_names,organism_name",
    "&format=tsv",
    "&size=5"
  )
  
  res <- httr::GET(url)
  
  if (httr::status_code(res) != 200) {
    return(tibble(
      nearest_gene_id_clean = gene_id,
      uniprot_accession = NA_character_,
      uniprot_reviewed = NA_character_,
      uniprot_protein_name = NA_character_,
      uniprot_gene_names = NA_character_,
      uniprot_organism = NA_character_
    ))
  }
  
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\n")[[1]]
  
  if (length(lines) <= 1) {
    return(tibble(
      nearest_gene_id_clean = gene_id,
      uniprot_accession = NA_character_,
      uniprot_reviewed = NA_character_,
      uniprot_protein_name = NA_character_,
      uniprot_gene_names = NA_character_,
      uniprot_organism = NA_character_
    ))
  }
  
  dat <- readr::read_tsv(I(txt), show_col_types = FALSE)
  
  dat %>%
    slice(1) %>%   
    transmute(
      nearest_gene_id_clean = gene_id,
      uniprot_accession = Entry,
      uniprot_reviewed = Reviewed,
      uniprot_protein_name = `Protein names`,
      uniprot_gene_names = `Gene Names`,
      uniprot_organism = Organism
    )
}
###UniProt Search###
uniprot_lookup <- bind_rows(lapply(gene_ids, lookup_uniprot_gene))

head(uniprot_lookup)

df_nearest_gene_annotated <- df_nearest_gene_annotated %>%
  mutate(
    nearest_gene_id_clean = gsub("^gene-", "", nearest_gene_id)
  ) %>%
  left_join(uniprot_lookup, by = "nearest_gene_id_clean")

df_nearest_gene_annotated <- df_nearest_gene_annotated %>%
  mutate(
    uniprot_accession_protein = ifelse(
      is.na(uniprot_accession),
      NA_character_,
      paste(uniprot_accession, uniprot_protein_name, sep = " | ")
    )
  )

df_nearest_gene_annotated %>%
  select(
    marker,
    nearest_gene_id,
    nearest_gene_name,
    nearest_description,
    uniprot_accession,
    uniprot_protein_name,
    uniprot_accession_protein,
    distance_to_gene
  ) %>%
  head(30)

##################Mutations##############

rm(list = ls())
gc()
graphics.off()


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "VariantAnnotation",
  "GenomicFeatures",
  "txdbmaker",
  "GenomicRanges",
  "IRanges",
  "Rsamtools",
  "Biostrings",
  "GenomeInfoDb"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(txdbmaker)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(Biostrings)
library(GenomeInfoDb)


setwd("/Users/Paul/DarpaMutations/")

fasta_file <- paste0(
  "/Users/Paul/DarpaMutations/",
  "GCF_002022765.2_C_virginica-3.0_genomic.fna"
)

gff_file <- paste0(
  "/Users/Paul/DarpaMutations/",
  "genomic.gff"
)

annotation_file <- paste0(
  "/Users/Paul/DarpaMutations/",
  "OysterCv SNP annotation.csv"
)

output_file <- paste0(
  "/Users/Paul/Desktop/",
  "SNP_coding_consequences.csv"
)

all_output_file <- paste0(
  "/Users/Paul/Desktop/",
  "All_SNP_coding_status.csv"
)

mismatch_output_file <- paste0(
  "/Users/Paul/Desktop/",
  "SNP_reference_mismatches.csv"
)


snps <- read.csv(
  annotation_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Input annotation rows:", nrow(snps), "\n")

required_columns <- c(
  "probeset_id",
  "custchr",
  "custpos",
  "Ref_Allele",
  "Alt_Allele"
)

missing_columns <- setdiff(
  required_columns,
  names(snps)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


if (!file.exists(paste0(fasta_file, ".fai"))) {
  Rsamtools::indexFa(fasta_file)
}

reference_genome <- Rsamtools::FaFile(
  fasta_file
)

fasta_index <- Rsamtools::scanFaIndex(
  reference_genome
)

fasta_seqnames <- unique(
  as.character(
    GenomicRanges::seqnames(fasta_index)
  )
)

fasta_lengths <- setNames(
  width(fasta_index),
  as.character(
    GenomicRanges::seqnames(fasta_index)
  )
)

cat("\nFASTA sequence names:\n")
print(fasta_seqnames)


txdb <- txdbmaker::makeTxDbFromGFF(
  file = gff_file,
  format = "gff3"
)

gff_seqnames <- GenomeInfoDb::seqlevels(
  txdb
)

cat("\nGFF sequence names:\n")
print(gff_seqnames)

cat(
  "\nNumber of transcripts:",
  length(GenomicFeatures::transcripts(txdb)),
  "\n"
)

cat(
  "Number of CDS entries:",
  length(GenomicFeatures::cds(txdb)),
  "\n"
)

if (length(GenomicFeatures::cds(txdb)) == 0) {
  stop(
    "The GFF transcript database contains no CDS features."
  )
}


chr_map <- c(
  "1"  = "NC_035780.1",
  "2"  = "NC_035781.1",
  "3"  = "NC_035782.1",
  "4"  = "NC_035783.1",
  "5"  = "NC_035784.1",
  "6"  = "NC_035785.1",
  "7"  = "NC_035786.1",
  "8"  = "NC_035787.1",
  "9"  = "NC_035788.1",
  "10" = "NC_035789.1"
)

snps <- snps %>%
  dplyr::mutate(
    probeset_id = as.character(
      probeset_id
    ),
    
    custchr_original = trimws(
      as.character(custchr)
    ),
    
    custchr_original = sub(
      "\\.0$",
      "",
      custchr_original
    ),
    
    custchr = unname(
      chr_map[custchr_original]
    ),
    
    custpos = suppressWarnings(
      as.integer(
        as.character(custpos)
      )
    ),
    
    Ref_Allele = toupper(
      trimws(
        as.character(Ref_Allele)
      )
    ),
    
    Alt_Allele = toupper(
      trimws(
        as.character(Alt_Allele)
      )
    )
  )

cat("\nChromosome conversion:\n")

print(
  table(
    original = snps$custchr_original,
    mapped = snps$custchr,
    useNA = "ifany"
  )
)


snp_seqnames <- unique(
  na.omit(snps$custchr)
)

common_seqnames <- Reduce(
  intersect,
  list(
    snp_seqnames,
    gff_seqnames,
    fasta_seqnames
  )
)

cat("\nShared SNP/GFF/FASTA chromosomes:\n")
print(common_seqnames)

if (length(common_seqnames) == 0) {
  stop(
    paste(
      "No chromosome names match"
    )
  )
}

txdb <- GenomeInfoDb::keepSeqlevels(
  txdb,
  value = common_seqnames,
  pruning.mode = "coarse"
)


snps_clean <- snps %>%
  dplyr::filter(
    !is.na(probeset_id),
    probeset_id != "",
    !is.na(custchr),
    custchr %in% common_seqnames,
    !is.na(custpos),
    custpos > 0,
    Ref_Allele %in% c("A", "C", "G", "T"),
    Alt_Allele %in% c("A", "C", "G", "T"),
    Ref_Allele != Alt_Allele
  )

cat("\nOriginal SNP rows:", nrow(snps), "\n")
cat("Rows retained:", nrow(snps_clean), "\n")
cat(
  "Rows removed:",
  nrow(snps) - nrow(snps_clean),
  "\n"
)


snp_ranges <- GenomicRanges::GRanges(
  seqnames = snps_clean$custchr,
  
  ranges = IRanges::IRanges(
    start = snps_clean$custpos,
    width = 1
  ),
  
  strand = "*"
)

names(snp_ranges) <- snps_clean$probeset_id

cat(
  "\nSNP ranges created:",
  length(snp_ranges),
  "\n"
)


range_chr_lengths <- unname(
  fasta_lengths[
    as.character(
      GenomicRanges::seqnames(snp_ranges)
    )
  ]
)

position_in_bounds <- (
  GenomicRanges::start(snp_ranges) >= 1 &
    GenomicRanges::start(snp_ranges) <= range_chr_lengths
)

cat("\nOriginal SNP positions inside chromosome boundaries:\n")

print(
  table(
    position_in_bounds,
    useNA = "ifany"
  )
)

position_in_bounds[
  is.na(position_in_bounds)
] <- FALSE

snp_ranges <- snp_ranges[
  position_in_bounds
]

snps_clean <- snps_clean[
  position_in_bounds,
  ,
  drop = FALSE
]


fasta_ref <- as.character(
  Rsamtools::scanFa(
    reference_genome,
    param = snp_ranges
  )
)

ref_match <- (
  toupper(snps_clean$Ref_Allele) ==
    toupper(fasta_ref)
)

ref_match[
  is.na(ref_match)
] <- FALSE

cat("\nReference allele comparison:\n")

print(
  table(
    ref_match,
    useNA = "ifany"
  )
)

cat(
  "Reference matches:",
  sum(ref_match),
  "\n"
)

cat(
  "Reference mismatches:",
  sum(!ref_match),
  "\n"
)

ref_mismatches <- snps_clean[
  !ref_match,
  ,
  drop = FALSE
]

if (nrow(ref_mismatches) > 0) {
  ref_mismatches$FASTA_REF <- fasta_ref[
    !ref_match
  ]
  
  write.csv(
    ref_mismatches,
    mismatch_output_file,
    row.names = FALSE,
    na = ""
  )
  
  cat(
    "\nSaved reference mismatches to:\n",
    mismatch_output_file,
    "\n"
  )
}


snp_ranges_valid <- snp_ranges[
  ref_match
]

snps_valid <- snps_clean[
  ref_match,
  ,
  drop = FALSE
]

var_alleles_valid <- Biostrings::DNAStringSet(
  snps_valid$Alt_Allele
)

cat("\nObject checks before predictCoding():\n")

cat(
  "Query class:",
  class(snp_ranges_valid)[1],
  "\n"
)

cat(
  "ALT class:",
  class(var_alleles_valid)[1],
  "\n"
)

cat(
  "Query length:",
  length(snp_ranges_valid),
  "\n"
)

cat(
  "ALT length:",
  length(var_alleles_valid),
  "\n"
)

stopifnot(
  inherits(snp_ranges_valid, "GRanges"),
  inherits(var_alleles_valid, "DNAStringSet"),
  length(snp_ranges_valid) ==
    length(var_alleles_valid)
)


seqinfo_valid <- GenomeInfoDb::Seqinfo(
  seqnames = common_seqnames,
  seqlengths = fasta_lengths[common_seqnames]
)

GenomeInfoDb::seqinfo(
  snp_ranges_valid
) <- seqinfo_valid


coding <- VariantAnnotation::predictCoding(
  query = snp_ranges_valid,
  subject = txdb,
  seqSource = reference_genome,
  varAllele = var_alleles_valid
)

cat(
  "\nCoding transcript consequences:",
  length(coding),
  "\n"
)

if (length(coding) == 0) {
  stop(
    paste(
      "predictCoding returned no coding consequences."
    )
  )
}


coding_marker <- names(coding)

names(coding) <- NULL

coding_df <- as.data.frame(
  coding,
  row.names = NULL
)

coding_df$marker <- coding_marker

cat(
  "\nCoding result rows:",
  nrow(coding_df),
  "\n"
)

cat(
  "Unique coding SNP markers:",
  length(unique(coding_df$marker)),
  "\n"
)


sequence_columns <- c(
  "REFCODON",
  "VARCODON",
  "REFAA",
  "VARAA",
  "varAllele"
)

for (
  column_name in intersect(
    sequence_columns,
    names(coding_df)
  )
) {
  coding_df[[column_name]] <- as.character(
    coding_df[[column_name]]
  )
}

if ("varAllele" %in% names(coding_df)) {
  coding_df$ALT <- coding_df$varAllele
}


coding_df$coding_class <- dplyr::case_when(
  coding_df$CONSEQUENCE == "synonymous" ~
    "Synonymous",
  
  coding_df$CONSEQUENCE == "nonsense" ~
    "Stop gained",
  
  coding_df$CONSEQUENCE == "nonsynonymous" &
    coding_df$REFAA == "*" ~
    "Stop lost",
  
  coding_df$CONSEQUENCE == "nonsynonymous" ~
    "Nonsynonymous",
  
  TRUE ~ as.character(
    coding_df$CONSEQUENCE
  )
)

cat("\nTranscript-level coding classifications:\n")

print(
  table(
    coding_df$coding_class,
    useNA = "ifany"
  )
)


result <- dplyr::left_join(
  coding_df,
  snps_valid,
  by = c(
    "marker" = "probeset_id"
  )
)


first_columns <- c(
  "marker",
  "custchr_original",
  "custchr",
  "custpos",
  "Ref_Allele",
  "Alt_Allele",
  "GENEID",
  "TXID",
  "TXNAME",
  "CDSID",
  "CDSLOC",
  "PROTEINLOC",
  "REFCODON",
  "VARCODON",
  "REFAA",
  "VARAA",
  "CONSEQUENCE",
  "coding_class"
)

first_columns <- intersect(
  first_columns,
  names(result)
)

result <- dplyr::select(
  result,
  dplyr::all_of(first_columns),
  dplyr::everything()
)


list_columns <- names(result)[
  vapply(
    result,
    is.list,
    logical(1)
  )
]

cat("\nList-like result columns being converted:\n")
print(list_columns)

for (column_name in list_columns) {
  result[[column_name]] <- vapply(
    result[[column_name]],
    
    function(x) {
      if (length(x) == 0) {
        return("")
      }
      
      x <- as.character(x)
      x <- x[!is.na(x)]
      
      if (length(x) == 0) {
        return("")
      }
      
      paste(
        x,
        collapse = ";"
      )
    },
    
    character(1)
  )
}

for (column_name in names(result)) {
  column_object <- result[[column_name]]
  
  if (
    inherits(column_object, "DNAStringSet") ||
    inherits(column_object, "AAStringSet") ||
    inherits(column_object, "CharacterList") ||
    inherits(column_object, "IntegerList") ||
    inherits(column_object, "CompressedList")
  ) {
    result[[column_name]] <- as.character(
      column_object
    )
  }
}

bad_columns <- names(result)[
  !vapply(
    result,
    
    function(x) {
      is.atomic(x) || is.factor(x)
    },
    
    logical(1)
  )
]

cat("\nColumns still unsafe for CSV export:\n")
print(bad_columns)

if (length(bad_columns) > 0) {
  stop(
    "These columns still cannot be exported: ",
    paste(
      bad_columns,
      collapse = ", "
    )
  )
}


write.csv(
  result,
  output_file,
  row.names = FALSE,
  na = ""
)

cat("\nSaved transcript-level coding results to:\n")
cat(output_file, "\n")

cat("\nTranscript-level classifications:\n")

print(
  table(
    result$coding_class,
    useNA = "ifany"
  )
)


coding_summary <- result %>%
  dplyr::group_by(marker) %>%
  dplyr::summarise(
    coding_class = paste(
      sort(
        unique(
          coding_class[
            !is.na(coding_class)
          ]
        )
      ),
      collapse = "; "
    ),
    
    gene_ids = paste(
      sort(
        unique(
          as.character(
            GENEID[
              !is.na(GENEID)
            ]
          )
        )
      ),
      collapse = "; "
    ),
    
    transcript_ids = paste(
      sort(
        unique(
          as.character(
            TXID[
              !is.na(TXID)
            ]
          )
        )
      ),
      collapse = "; "
    ),
    
    amino_acid_changes = paste(
      sort(
        unique(
          paste0(
            as.character(REFAA),
            as.character(PROTEINLOC),
            as.character(VARAA)
          )
        )
      ),
      collapse = "; "
    ),
    
    codon_changes = paste(
      sort(
        unique(
          paste0(
            as.character(REFCODON),
            ">",
            as.character(VARCODON)
          )
        )
      ),
      collapse = "; "
    ),
    
    .groups = "drop"
  )


all_snps_result <- dplyr::left_join(
  snps_clean,
  coding_summary,
  by = c(
    "probeset_id" = "marker"
  )
)

all_snps_result <- all_snps_result %>%
  dplyr::mutate(
    coding_class = dplyr::if_else(
      is.na(coding_class) |
        coding_class == "",
      "Not in annotated CDS",
      coding_class
    )
  )


all_list_columns <- names(all_snps_result)[
  vapply(
    all_snps_result,
    is.list,
    logical(1)
  )
]

for (column_name in all_list_columns) {
  all_snps_result[[column_name]] <- vapply(
    all_snps_result[[column_name]],
    
    function(x) {
      if (length(x) == 0) {
        return("")
      }
      
      x <- as.character(x)
      x <- x[!is.na(x)]
      
      if (length(x) == 0) {
        return("")
      }
      
      paste(
        x,
        collapse = ";"
      )
    },
    
    character(1)
  )
}


write.csv(
  all_snps_result,
  all_output_file,
  row.names = FALSE,
  na = ""
)

cat("\nSaved one-row-per-SNP results to:\n")
cat(all_output_file, "\n")

cat("\nOne-row-per-SNP classifications:\n")

print(
  table(
    all_snps_result$coding_class,
    useNA = "ifany"
  )
)

cat("\nFinished successfully.\n")

library(readxl)
library(dplyr)


gsm_file <- "/Users/Paul/DarpaMutations/FinalGSMTable copy.xlsx"

syno_file <- "/Users/Paul/Desktop/All_SNP_coding_status.csv"

output_file <- "/Users/Paul/Desktop/FinalGSMTable_SNP_coding_status.csv"


gsm <- readxl::read_excel(
  gsm_file
)

syno <- read.csv(
  syno_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


names(gsm)[1] <- "SNP"

if (!"SNP" %in% names(gsm)) {
  stop(
    "The FinalGSMTable file does not contain a column named SNP."
  )
}

if (!"probeset_id" %in% names(syno)) {
  stop(
    "The SNP coding-status file does not contain a column named probeset_id."
  )
}


gsm <- gsm %>%
  mutate(
    SNP = trimws(as.character(SNP))
  )

syno <- syno %>%
  mutate(
    probeset_id = trimws(as.character(probeset_id))
  )


cat("Rows in FinalGSMTable:", nrow(gsm), "\n")
cat("Unique SNPs in FinalGSMTable:", length(unique(gsm$SNP)), "\n")

cat(
  "SNPs matching coding annotation:",
  sum(gsm$SNP %in% syno$probeset_id),
  "\n"
)

cat(
  "SNPs not matching coding annotation:",
  sum(!gsm$SNP %in% syno$probeset_id),
  "\n"
)

unmatched_snps <- gsm %>%
  filter(
    is.na(SNP) |
      SNP == "" |
      !SNP %in% syno$probeset_id
  ) %>%
  distinct(SNP)

cat("\nUnmatched SNPs:\n")
print(unmatched_snps)


syno_subset <- syno %>%
  filter(
    probeset_id %in% gsm$SNP
  )

write.csv(
  syno_subset,
  "/Users/Paul/Desktop/FinalGSM_SNP_annotation_subset.csv",
  row.names = FALSE,
  na = ""
)


gsm_annotated <- gsm %>%
  left_join(
    syno,
    by = c("SNP" = "probeset_id")
  )

write.csv(
  gsm_annotated,
  output_file,
  row.names = FALSE,
  na = ""
)


cat("\nCoding classifications in FinalGSMTable:\n")

print(
  table(
    gsm_annotated$coding_class,
    useNA = "ifany"
  )
)

cat("\nSubset annotation saved to:\n")
cat("/Users/Paul/Desktop/FinalGSM_SNP_annotation_subset.csv\n")

cat("\nFinalGSMTable with coding annotations saved to:\n")
cat(output_file, "\n")


test_snp <- snps_clean %>%
  dplyr::filter(probeset_id == "AX-563280343")

test_range <- GenomicRanges::GRanges(
  seqnames = test_snp$custchr,
  ranges = IRanges::IRanges(
    start = test_snp$custpos,
    width = 1
  )
)

fasta_base <- as.character(
  Rsamtools::scanFa(
    reference_genome,
    param = test_range
  )
)

data.frame(
  marker = test_snp$probeset_id,
  chromosome = test_snp$custchr,
  position = test_snp$custpos,
  annotation_REF = test_snp$Ref_Allele,
  annotation_ALT = test_snp$Alt_Allele,
  FASTA_base = fasta_base,
  REF_matches_FASTA =
    test_snp$Ref_Allele == fasta_base
)





#################################

###GO Scoring###

library(readxl)
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(GO.db)
library(writexl)

setwd("/Users/Paul/Darpa All Gen")

###File with SwissProt GOs###
df <- readxl::read_xlsx(
  "/Users/Paul/Darpa All Gen/SwissProt222.xlsx"
)
colnames(df) <- c("pvalue", "GO")


###Chr.P threshold###
threshold <- 1.0381e-05
log_threshold <- -log10(threshold)

df <- df %>%
  dplyr::mutate(
    SNP_ID = dplyr::row_number(),
    pvalue = as.numeric(pvalue)
  ) %>%
  tidyr::separate_rows(
    GO,
    sep = ";"
  ) %>%
  dplyr::mutate(
    GO = trimws(GO)
  ) %>%
  dplyr::filter(
    !is.na(GO),
    GO != "",
    !is.na(pvalue),
    pvalue > 0
  ) %>%
  dplyr::distinct(
    SNP_ID,
    GO,
    .keep_all = TRUE
  )

df_ready <- df %>%
  dplyr::mutate(
    neg_log_p = -log10(pvalue),
    score_weighted = neg_log_p / log_threshold
  )


go_terms <- AnnotationDbi::select(
  GO.db,
  keys = unique(df_ready$GO),
  columns = c("GOID", "ONTOLOGY"),
  keytype = "GOID"
) %>%
  dplyr::distinct(
    GOID,
    .keep_all = TRUE
  )

###Keep only BP GO###

df_ready <- df_ready %>%
  dplyr::left_join(
    go_terms,
    by = c("GO" = "GOID")
  ) %>%
  dplyr::filter(
    ONTOLOGY == "BP"
  )


df_ready$label <- AnnotationDbi::mapIds(
  GO.db,
  keys = df_ready$GO,
  column = "TERM",
  keytype = "GOID",
  multiVals = "first"
)

df_ready <- df_ready %>%
  dplyr::filter(
    !is.na(label)
  )

###Grouping GO Terms

df_ready <- df_ready %>%
  dplyr::mutate(
    cluster = dplyr::case_when(
      

      
      grepl(
        paste0(
          "immune|innate immune|defense response|inflammat|",
          "pathogen|bacteri|viral|virus|antigen|microbiota|",
          "host.microbe|host pathogen|mucus|mucin|mucosal|",
          "maintenance of gastrointestinal epithelium|",
          "epithelial barrier|epithelium maintenance|",
          "glycan|glycosylation|O-linked|lectin|complement|",
          "phagocyt|hemocyte|wounding|response to toxic|",
          "detoxification|oxidative stress|reactive oxygen|",
          "NF-kappaB|Toll|pattern recognition|microbial|",
          "platelet activation|T cell proliferation|",
          "leukocyte|cytokine|interferon"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Immune Functions",
      
      grepl(
        paste0(
          "biomineral|mineralization|calcification|carbonate|",
          "shell|extracellular matrix organization|",
          "extracellular structure organization|",
          "extracellular matrix assembly|collagen fibril|",
          "chitin|skeletal|odont|dentin|osteoblast|",
          "odontoblast|matrix organization"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Biomineralization",
      
      
      ###Apoptosis later grouped with Cell Cycle###
      grepl(
        paste0(
          "apoptotic|apoptosis|programmed cell death|",
          "cell death|cell survival|anti-apopt|caspase"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Apoptosis",
      
      
      grepl(
        paste0(
          "proteolysis|protein catabolic|ubiquitin|",
          "deubiquitin|proteasome|protein folding|",
          "protein stability|protein quality|",
          "protein modification|protein processing|",
          "unfolded protein|chaperone|autophagy"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Proteostasis",
      
      
      grepl(
        paste0(
          "transcription|RNA polymerase|",
          "regulation of gene expression|chromatin|histone|",
          "mRNA processing|RNA processing|RNA splicing|",
          "RNA 3'|RNA destabil|miRNA|microRNA|pre-miRNA|",
          "primary miRNA|genomic imprinting|",
          "RNA modification|RNA localization|",
          "RNA stability|translation regulation"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Transcription & RNA Regulation",
      
      
      grepl(
        paste0(
          "mitotic|meiotic|cell cycle|",
          "chromosome segregation|chromosome organization|",
          "spindle|kinetochore|chromatid|metaphase|",
          "cytokinesis|DNA replication|DNA repair|",
          "sister chromatid|G1|G2|DNA damage|",
          "telomere|centromere"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Cell Cycle",
      
      
      grepl(
        paste0(
          "neuron|neural|axon|dendrite|synapse|synaptic|",
          "presynaptic|postsynaptic|glial|neurotransmitter|",
          "sensory|mechanosensory|chemosensory|photoreceptor|",
          "taste receptor|gustatory|cerebell|cerebral cortex|",
          "brain development|spinal cord|neural tube|",
          "neural crest|floor plate|commissure|",
          "substantia nigra|cochlear nucleus|neurofilament|",
          "neuromuscular|behavior|learning|memory|",
          "response to stimulus"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Neurosensory",
      
      
      grepl(
        paste0(
          "mitochondrial|mitochondrion|",
          "oxidative phosphorylation|ATP synthesis|",
          "electron transport|respiratory chain|",
          "proton motive|energy derivation|",
          "mitochondrion organization|",
          "aerobic respiration|cellular respiration|",
          "ATP metabolic process"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Mitochondria & Energy",
      
      
      grepl(
        paste0(
          "signaling|signalling|signal transduction|",
          "cell communication|receptor signaling|",
          "receptor signalling|G protein|Wnt|SMAD|MAPK|",
          "JNK cascade|JAK|STAT|integrin|adenylate|",
          "phospholipase|smoothened|plexin|inositol|",
          "protein phosphorylation|protein autophosphorylation|",
          "response to retinoic acid|hormone|thyroid|",
          "endocrine|neuroendocrine|norepinephrine secretion|",
          "catecholamine secretion|second messenger|",
          "receptor activity regulation"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Cell Signalling",
      
      
      grepl(
        paste0(
          "actin|cytoskeleton|actomyosin|microtubule|",
          "stress fiber|stress fibre|cell motility|",
          "cell migration|cell adhesion|focal adhesion|",
          "cilium|cilia|flagell|muscle contraction|",
          "contractile|heart contraction|cell shape|",
          "cell projection|microfilament"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Cell Structure",
      
      
      grepl(
        paste0(
          "transport|localization|localisation|",
          "ion homeostasis|ion transmembrane|solute|",
          "amino acid transport|metal ion transport|",
          "calcium ion transport|proton transport|",
          "sodium|potassium|chloride|osmoreg|osmotic|",
          "water homeostasis|plasma membrane fusion|",
          "secretion|exocytosis|endocytosis|",
          "vesicle-mediated|vesicle transport|Golgi|",
          "secretory|protein transport|nuclear import|",
          "nuclear export|membrane trafficking"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Transport",
      
      
      grepl(
        paste0(
          "cell differentiation|cell development|cell fate|",
          "development|developmental|morphogenesis|formation|",
          "pattern formation|pattern specification|",
          "proximal.distal pattern|midline development|",
          "mesoderm|epidermis development|",
          "epithelial cell proliferation|",
          "endothelial cell proliferation|",
          "cell population proliferation|",
          "multicellular organismal process|",
          "multicellular organismal development|",
          "tissue development|tissue remodeling|",
          "organ development|tube development|branching|",
          "embryonic|larval|epithelial differentiation|",
          "epidermal differentiation|",
          "digestive system development|growth|",
          "angiogenesis|vasculogenesis|",
          "keratinocyte differentiation"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Cell Development",
      
      
      grepl(
        paste0(
          "metabolic|metabolism|biosynthetic|catabolic|",
          "lipid|fatty acid|amino acid|carbohydrate|",
          "glycogen|glucose|cholesterol|triglyceride|",
          "glycerolipid|nucleoside|nucleotide|xenobiotic|",
          "acyl-CoA|phosphatidic|taurine|monocarboxylic|",
          "steroid|lipoprotein|glycolysis|gluconeogenesis|",
          "oxidation|reduction|redox|cofactor|vitamin"
        ),
        label,
        ignore.case = TRUE
      ) ~ "Metabolism",
      
      TRUE ~ "Other"
    )
  )


snp_cluster_counts <- df_ready %>%
  dplyr::distinct(
    SNP_ID,
    cluster
  ) %>%
  dplyr::count(
    SNP_ID,
    name = "number_of_clusters"
  )

###Scoring only once per cluster per SNP###

cluster_snp_scores <- df_ready %>%
  dplyr::group_by(
    cluster,
    SNP_ID
  ) %>%
  dplyr::summarise(
    pvalue = dplyr::first(pvalue),
    raw_SNP_score = dplyr::first(score_weighted),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    snp_cluster_counts,
    by = "SNP_ID"
  ) %>%
  dplyr::mutate(
    SNP_score = raw_SNP_score / number_of_clusters
  )

###Adding terms to df####

cluster_go_terms <- df_ready %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::summarise(
    unique_GO_count = dplyr::n_distinct(GO),
    
    contributing_GO_IDs = paste(
      sort(unique(GO)),
      collapse = "; "
    ),
    
    contributing_GO_terms = paste(
      sort(unique(label)),
      collapse = "; "
    ),
    
    .groups = "drop"
  )

cluster_summary <- cluster_snp_scores %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::summarise(
    total_score = sum(
      SNP_score,
      na.rm = TRUE
    ),
    
    unique_SNP_count = dplyr::n_distinct(
      SNP_ID
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    cluster_go_terms,
    by = "cluster"
  ) %>%
  dplyr::arrange(
    dplyr::desc(total_score)
  )
