
#Define iterations and fractioning
n_iter <- 50
train_frac <- 0.80

common_ids <- intersect(
  as.character(phenTrain$ID),
  rownames(GenoRoundedAll)
)

length(common_ids)

# Set Seed
set.seed(123)


# Iterative Loop
for (i in 1:n_iter) {
  
  cat("\n============================\n")
  cat("Iteration:", i, "of", n_iter, "\n")
  cat("============================\n")
  

  # Fractionate
  ids80 <- sample(
    common_ids,
    size = floor(length(common_ids) * train_frac),
    replace = FALSE
  )
  

  # Subset phen

  phen80 <- phenTrain[
    as.character(phenTrain$ID) %in% ids80,
    ,
    drop = FALSE
  ]
  
  # Reorder
  phen80 <- phen80[
    match(ids80, as.character(phen80$ID)),
    ,
    drop = FALSE
  ]
  

  # Subset Geno

  geno80 <- GenoRoundedAll[
    ids80,
    ,
    drop = FALSE
  ]
  

  # Pre-Gwas

  pre80 <- pre.gwas(
    pheno.data = phen80,
    indiv = "ID",
    resp = "Status",
    geno.data = geno80,
    Q.method = "K",
    method = "VanRaden",
    maf = 0,
    marker.callrate = 1,
    ind.callrate = 1,
    impute = FALSE
  )
  

  # GWAS

  gwas80 <- gwas.asreml(
    pheno.data = pre80$pheno.data,
    resp = "Status",
    gen = "ID",
    Kinv = pre80$Kinv,
    Q = pre80$Q,
    npc = 7,
    family = c("gaussian"),
    geno.data = pre80$geno.data,
    map.data = pre80$map.data,
    pvalue.thr = 0.05,
    bonferroni = TRUE,
    threads = 1,
    workspace = "4Gb",
    P3D = TRUE
  )
  

  # Save all but NAs

  all80 <- gwas80$gwas.all[
    !is.na(gwas80$gwas.all$p.value),
    ,
    drop = FALSE
  ]
  
  all80$iteration <- i
  
  cat(
    "Markers saved:",
    nrow(all80),
    "\n"
  )
  

  # Save iteration as an R obj

  saveRDS(
    all80,
    paste0("GWAS80_iteration_", i, ".rds")
  )
  

  # Clean and repeat

  rm(
    ids80,
    phen80,
    geno80,
    pre80,
    gwas80,
    all80
  )
  
  gc()
}

#Outside of Loop


# Combine Iterations

all_results <- lapply(
  1:50,
  function(i) {
    readRDS(
      paste0("GWAS80_iteration_", i, ".rds")
    )
  }
)


# Remove NA Markers

all_results80 <- do.call(
  rbind,
  lapply(all_results, as.data.frame)
)


# Remove NA pvalues

all_results80 <- all_results80[
  !is.na(all_results80$p.value),
  ,
  drop = FALSE
]


# Save R obj

saveRDS(
  all_results80,
  "GWAS_80percent_50rep_ALL_combined.rds"
)


# Save xlsx

wb <- createWorkbook()

for (i in 1:length(all_results)) {
  
  sheet_name <- paste0("Iteration_", i)
  
  addWorksheet(
    wb,
    sheet_name
  )
  
  writeData(
    wb,
    sheet = sheet_name,
    x = all_results[[i]],
    rowNames = FALSE
  )
}

saveWorkbook(
  wb,
  "GWAS_80percent_50rep_ALL_results.xlsx",
  overwrite = TRUE
)

cat("\nTotal rows combined:", nrow(all_results80), "\n")

cat(
  "Any NA p-values remaining:",
  any(is.na(all_results80$p.value)),
  "\n"
)

cat(
  "Iterations present:",
  paste(
    sort(unique(all_results80$iteration)),
    collapse = ", "
  ),
  "\n"
)

# Read R obj back

rep80 <- readRDS(
  "GWAS_80percent_50rep_ALL_combined.rds"
)

# Read in or Create 'True' Full Test Set GWAS
true_gwas <- read.xlsx(
  "DARPA_AllMarkers_Final.xlsx"
)

true_gwas<-gwasS$gwas.all

# Matching SNP names
true_gwas$marker <- gsub("-", "_", true_gwas$marker)
rep80$marker <- gsub("-", "_", rep80$marker)



# Rescaling

rep_summary <- rep80 %>%
  filter(
    !is.na(p.value),
    p.value > 0
  ) %>%
  mutate(
    logp = -log10(p.value)
  ) %>%
  group_by(marker) %>%
  summarise(
    mean_logp = mean(logp),
    sd_logp = sd(logp),
    n_iterations = n_distinct(iteration),
    .groups = "drop"
  )



# Rescaling

true_snps <- true_gwas %>%
  filter(
    !is.na(p.value),
    p.value > 0
  ) %>%
  transmute(
    marker,
    true_logp = -log10(p.value)
  )



# Combining iterations and True set

corr_data <- inner_join(
  true_snps,
  rep_summary,
  by = "marker"
)


table(corr_data$n_iterations)

dim(corr_data)



# Pearson Corr

cor_value <- cor(
  corr_data$true_logp,
  corr_data$mean_logp,
  method = "pearson",
  use = "complete.obs"
)

cor_value

#Define GSMs
GSM_markers <- c(
  # insert GSMs here
)
  
  # Define only GSM set
  GSM_data <- corr_data %>%
    filter(marker %in% GSM_markers)
  
  # Define Full Set
  all_data <- corr_data
  
  
  # Develop Theoretical Relationships (ie 1:1 and 0.8:1) and Fxn for RMSE
  agreement_stats <- function(dat) {
    
    residual_1to1 <- dat$mean_logp - dat$true_logp
    residual_08   <- dat$mean_logp - (0.8 * dat$true_logp)
    
    data.frame(
      n = nrow(dat),
      
      # Distance from theoretical 1:1 line
      RMSE_1to1 = sqrt(mean(residual_1to1^2, na.rm = TRUE)),
      MeanResidual_1to1 = mean(residual_1to1, na.rm = TRUE),
      SDResidual_1to1 = sd(residual_1to1, na.rm = TRUE),
      
      # Distance from theoretical 0.8:1 line
      RMSE_08 = sqrt(mean(residual_08^2, na.rm = TRUE)),
      MeanResidual_08 = mean(residual_08, na.rm = TRUE),
      SDResidual_08 = sd(residual_08, na.rm = TRUE)
    )
  }
  
  
  # GSMS
  GSM_stats <- agreement_stats(GSM_data)
  
  # Whole Dataset
  all_stats <- agreement_stats(all_data)
  
  agreement_results <- bind_rows(
    cbind(Group = "GSM SNPs", GSM_stats),
    cbind(Group = "All markers", all_stats)
  )
  
  print(agreement_results)
  cat(
    "\nGSM SNPs\n",
    "1:1 mean residual = ",
    round(GSM_stats$MeanResidual_1to1, 6),
    " ± ",
    round(GSM_stats$SDResidual_1to1, 6),
    "; RMSE = ",
    round(GSM_stats$RMSE_1to1, 6),
    "\n",
    "0.8:1 mean residual = ",
    round(GSM_stats$MeanResidual_08, 6),
    " ± ",
    round(GSM_stats$SDResidual_08, 6),
    "; RMSE = ",
    round(GSM_stats$RMSE_08, 6),
    
    "\n\nAll markers\n",
    "1:1 mean residual = ",
    round(all_stats$MeanResidual_1to1, 6),
    " ± ",
    round(all_stats$SDResidual_1to1, 6),
    "; RMSE = ",
    round(all_stats$RMSE_1to1, 6),
    "\n",
    "0.8:1 mean residual = ",
    round(all_stats$MeanResidual_08, 6),
    " ± ",
    round(all_stats$SDResidual_08, 6),
    "; RMSE = ",
    round(all_stats$RMSE_08, 6),
    "\n"
  )
  
  # Pearson Corr for GSMs
  cor_GSM <- cor.test(
    GSM_data$true_logp,
    GSM_data$mean_logp,
    method = "pearson"
  )
  
  # Pearson for All
  cor_all <- cor.test(
    all_data$true_logp,
    all_data$mean_logp,
    method = "pearson"
  )
  
  cat(
    "GSM SNPs: r =",
    round(unname(cor_GSM$estimate), 6),
    "; 95% CI =",
    round(cor_GSM$conf.int[1], 6),
    "to",
    round(cor_GSM$conf.int[2], 6),
    "\n",
    "All markers: r =",
    round(unname(cor_all$estimate), 6),
    "; 95% CI =",
    round(cor_all$conf.int[1], 6),
    "to",
    round(cor_all$conf.int[2], 6),
    "\n"
  )
  
  GSM_pm <- (cor_GSM$conf.int[2] - cor_GSM$conf.int[1]) / 2
  all_pm <- (cor_all$conf.int[2] - cor_all$conf.int[1]) / 2
  
  cat(
    "GSM SNPs: r =",
    round(unname(cor_GSM$estimate), 6),
    "±",
    round(GSM_pm, 6),
    "\n",
    "All markers: r =",
    round(unname(cor_all$estimate), 6),
    "±",
    round(all_pm, 6),
    "\n"
  )
  

# Graphics

p_corr_all <- ggplot(
  corr_data,
  aes(
    x = true_logp,
    y = mean_logp
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.25
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0(
      "Pearson r = ",
      round(cor_value, 2)
    ),
    hjust = 1.1,
    vjust = -1,
    size = 5
  ) +
  labs(
    x = expression(-log[10](P)~"from Full GWAS"),
    y = expression("Mean " * -log[10](P)~"Across 80% Replicates")
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_corr_all

corr_data <- corr_data %>%
  arrange(true_logp)

p_corr_all <- ggplot(
  corr_data,
  aes(
    x = true_logp,
    y = mean_logp
  )
) +
  geom_ribbon(
    aes(
      ymin = mean_logp - sd_logp,
      ymax = mean_logp + sd_logp
    ),
    alpha = 0.15
  ) +
  geom_abline(
    slope = 0.8,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.25
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0(
      "Pearson r = ",
      round(cor_value, 2)
    ),
    hjust = 1.1,
    vjust = -1,
    size = 5
  ) +
  labs(
    x = expression(-log[10](P)~"from Full GWAS"),
    y = expression("Mean " * -log[10](P)~"Across 80% Replicates")
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

########################################## PCA and LD Pruning Stuff ##########
# Repeat for 0.01 / 0.05 / All, etc.

#Remove MT and align for kbps
keep_markers <- Map$Marker[Map$Chr != 11]
GenoRoundedAll <- GenoRoundedAll[
  ,
  colnames(GenoRoundedAll) %in% keep_markers
]

Map <- Map[Map$Chr != 11, ]
geno_mat <- as.matrix(GenoRoundedAll)
map_match <- Map[match(colnames(GenoRoundedAll), Map$Marker), ]

# GDS Create

snpgdsCreateGeno(
  gds.fn = "GenoRounded0.01.gds",
  genmat = geno_mat,
  sample.id = rownames(geno_mat),
  snp.id = map_match$Marker,
  snp.chromosome = as.integer(map_match$Chr),
  snp.position = as.integer(map_match$Pos),
  snpfirstdim = FALSE
)

# GDS open
genofile <- snpgdsOpen("GenoRounded0.01.gds")
snpgdsSummary("GenoRounded0.01.gds")

# GDS Pruning, tested multiple slide windows, minimal change 50kbp = middle
ld_pruned <- snpgdsLDpruning(
  genofile,
  autosome.only = FALSE,
  maf = 0,
  missing.rate = 1,
  method = "corr",
  slide.max.bp = 50000,
  ld.threshold = sqrt(0.2),
  start.pos = "first",
  num.thread = 4
)

pruned_snps <- unlist(ld_pruned, use.names = FALSE)
pruned0.01 <- pruned_snps

snpgdsClose(genofile)
 
# PCA
GenoPCA <- scale(Geno_pruned0.01,center=T,scale=F) 
PCAResults <- svd(GenoPCA) 
GraphPCAResults <- as.data.frame(GenoPCA%*%PCAResults$v)

# PCs
PCA1 <- 100*round((PCAResults$d[1])^2/sum((PCAResults$d)^2),d=3); PCA1
PCA2 <- 100*round((PCAResults$d[2])^2/sum((PCAResults$d)^2),d=3); PCA2

GraphPCAResults<-PCAResults$u
GraphPCAResults<-as.data.frame(GraphPCAResults)

# Add names back
rownames(GraphPCAResults) <- rownames(Geno_pruned0.01)
GraphPCAResults$ID <- rownames(GraphPCAResults)

GraphPCAResults$Status <- phenTrain$Status[
  match(GraphPCAResults$ID, phenTrain$ID)
]

# Generation Assignments
GraphPCAResults$Generation <- NA

GraphPCAResults$Generation[
  GraphPCAResults$ID %in% F0_ids
] <- "F0"

GraphPCAResults$Generation[
  GraphPCAResults$ID %in% F1_ids
] <- "F1"

GraphPCAResults$Generation[
  GraphPCAResults$ID %in% F2_ids
] <- "F2"

GraphPCAResults$Generation <- factor(
  GraphPCAResults$Generation,
  levels = c("F0", "F1", "F2")
)

# Status as Factor
GraphPCAResults$Status <- factor(GraphPCAResults$Status)
pc_cols <- setdiff(
  names(GraphPCAResults),
  c("ID", "Status", "Generation")
)

names(GraphPCAResults)[match(pc_cols, names(GraphPCAResults))] <-
  paste0("PC", seq_along(pc_cols))

# Make F0 pop better
GraphPCAResults <- GraphPCAResults[
  order(GraphPCAResults$Generation == "F0"),
]

# Graphics
p1 <- ggplot(
  GraphPCAResults,
  aes(
    x = PC1,
    y = PC2,
    color = Generation,
    shape = Status
  )
) +
  geom_point(size = 3) +
  scale_shape_manual(
    values = c(
      "0" = 4,
      "1" = 1
    ),
    labels = c(
      "0" = "Dead",
      "1" = "Live"
    )
  ) +
  scale_color_manual(
    values = c(
      "F0" = "#0072B2",
      "F1" = "#D55E00",
      "F2" = "#A23BEC"
    )
  ) +
  labs(
    x = paste0("PC1 (", round(PCA1, 2), "%)"),
    y = paste0("PC2 (", round(PCA2, 2), "%)"),
    color = "Generation",
    shape = "Status"
  ) +
  theme_classic(
    base_family = "",
    base_size = 14
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5)
    ),
    shape = guide_legend(
      override.aes = list(size = 5)
    )
  )

p1

ggsave(
  filename = "0.01PCAGenerationStatus.pdf",
  plot = p1,
  width = 7,
  height = 5.5,
  units = "in"
)

p_combined <- (p + theme(legend.position = "none")) + p1 +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

p_combined

ggsave(
  "AandBPCAGenerationStatus.pdf",
  plot = p_combined,
  width = 14,
  height = 5.5,
  units = "in"
)

p_combined
