rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("linc2.RData")

library(dplyr)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(vsn) # to plot shrinkage variance
library(gridExtra) # to have multiple plots in grids
library(dendextend) # to plot dendrogram
library(car) #Levene's test

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)

## Load data -------------
sz_count <- read.csv("C:/Users/chuot/Documents/Skoltech Study/Research/lncRNA/sz_human.counts_with.biotypes.csv", header = T)
# Filter lincRNA genes only
lincRNA <- sz_count[sz_count$Biotype == "lincRNA",]
lincRNA <- lincRNA[,1:ncol(lincRNA)-1]
# Move X to row names 
rownames(lincRNA) <- lincRNA$X
lincRNA$X <- NULL

# Meta data
description <- read.csv("table_sample_description.csv", header = T, row.names = 1)
rownames(description) <- description$sample

## Filter non-zero genes ---------------------

filter_gene <- function(counts, max_zero, description) {
  # Each region has 4 samples
  regions <- split(description$sample, description$sample_region)
  # Remove genes wih maximal 2 zero counts
  for (region in regions) {
    counts %>%
      dplyr::select(all_of(region)) == 0 -> x
    counts[rowSums(x) <= max_zero, ] -> counts
  }
  return(counts)
}

filtered_lincRNA <- filter_gene(lincRNA, 2, description)
filtered_lincRNA <- as.matrix(filtered_lincRNA)

## Transform + normalize for size factor + variance shrinkage --------------------
# Form DESeq Data Set (dds)
dds <- DESeqDataSetFromMatrix(countData = filtered_lincRNA, 
                              colData = description,
                              design = ~ condition + Acronym)

# Log 2 transformation of read counts including variance shrinkage
rlog_transformed <- rlog(dds,blind = TRUE)
rlog_transformed <- assay(rlog_transformed)

msd_plot <- meanSdPlot(rlog_transformed,
                       ranks = FALSE,
                       plot = FALSE)
msd_plot$gg + ylab("standard deviation")

pca_plot <- function(df, description, group) {
  pca <- prcomp(t(df))
  plot <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, col = description[,group])) + 
    geom_point() + labs(col = group)
  return(plot)
}

(pca_raw <- pca_plot(rlog_transformed, description, "condition"))

## Residuals of RIN normalization ------
# Matrix of size 0 x samples

take_residuals <- function(df, x_var) {
  residuals <- matrix(nrow = 0, ncol = ncol(df))
  colnames(residuals) <- colnames(df)
  # Keep residuals of regression for all genes
  for (i in 1:nrow(df)) {
    model <- lm(df[i,] ~ description[,x_var])
    residuals <- rbind(residuals, model$residuals)
  }
  rownames(residuals) <- rownames(df)
  return(residuals)
}

RIN_norm <- take_residuals(rlog_transformed, "RIN")

## Individual value ------------------

difference_to_brain_mean <- function (matrix, description) {
  # Merge info
  df <- as.data.frame(t(matrix))
  df <- merge(df, description[,c('condition', 'Acronym', 'individuum')], by = "row.names", all = T)
  
  genes <- rownames(matrix)
  # Calculate mean for each brains
  df %>%
    group_by(individuum) %>%
    dplyr::summarise_at(.vars = genes, .funs =  function(i) {i - mean(i)}) -> differ_from_mean
  # Reformat
  indi_norm <- t(differ_from_mean[,genes])
  colnames(indi_norm) <- description$sample
  return(indi_norm)
}
individual_norm <- difference_to_brain_mean(RIN_norm, description)

(pca_indi <- pca_plot(individual_norm, description, "condition"))
grid.arrange(pca_raw,pca_indi, nrow=1)

## Clustering -------------

difference_by_regions <- function(matrix, description) {
  # Combine matrix with description
  df_indi <- as.data.frame(t(matrix))
  df_indi <- merge(df_indi, description, by = "row.names", all = T)  
  genes <- rownames(matrix)
  # Mean of condition_region
  df_indi %>%
    group_by(condition, Acronym) %>%
    summarise_if(.predicate = is.numeric, .funs = mean) -> mean_regions
  # Difference in each region
  dif_by_regions <- mean_regions[mean_regions$condition == "X", genes] - mean_regions[mean_regions$condition == "H", genes]
  rownames(dif_by_regions) <- mean_regions$Acronym[mean_regions$condition == "H"]
  return(dif_by_regions)
}

make_dhc <- function(matrix, k) {
  # Pearson correlation
  pearson.cor <- cor(matrix, method = "pearson")
  distance <- as.dist(1 - pearson.cor)
  # Get descending hierarchical clustering
  dhc <- as.dendrogram(hclust(distance))
  dhc <- color_branches(dhc, k = k)
  dhc <- color_labels(dhc, k = k)
  plot <- (ggplot(as.ggdend(dhc), horiz = TRUE))
  print(plot)
  return(dhc)
}

mean_regions <- difference_by_regions(individual_norm, description)
dhc <- make_dhc(t(mean_regions),description)

dend_list <- get_subdendrograms(dhc, k = 4)
dend_list <- lapply(dend_list, labels)
names(dend_list) <- as.roman(1:length(dend_list))

dend_list <- melt(dend_list)
colnames(dend_list) <- c("Acronym", "Cluster")

metadat_clust <- merge(x = description,  y = dend_list, by = "Acronym")
metadat_clust <- metadat_clust[order(metadat_clust$sample),]

pca_plot(individual_norm, metadat_clust, "Cluster")

## Levene's Test-------------
# Prepare df
df_indi <- melt(individual_norm)
names(df_indi) <- c("gene", "sample", "value")
df_indi <- merge(df_indi, description, by = "sample")

#Actual test
lapply(split(df_indi, df_indi$gene), function(i) {
  a <- leveneTest(value ~ Acronym*condition, data = i)
  levene_p <- a$`Pr(>F)`[1]
  return(levene_p)
}) -> levene_pvals

levene_pvals <- data.frame(unlist(levene_pvals))
sum(levene_pvals > 0.05)

## ANOVA new -----------

# aaaa <- df_indi[df_indi$gene == df_indi$gene[2], c("gene","Acronym","condition", "value")]
# write.csv(aaaa, "scz_data.csv",row.names = F)
# 
# ggplot(aaaa, aes(x = Acronym, y = value, col = condition)) + geom_boxplot() + labs(x = "Region", y = "Expression level", col = "Condition")
# 
# summary(lm(value ~ condition + Acronym + Acronym:condition, data = df_indi[df_indi$gene == df_indi$gene[2],]))
# a <- aov(value ~ Acronym:condition + condition + Acronym, data = df_indi[df_indi$gene == df_indi$gene[2],])
# 
# aa <- TukeyHSD(a)
# names(aa$`Acronym:condition`)
# aaa <- as.data.frame(aa$`Acronym:condition`)
# aaa[paste(unique(description$Acronym),":X-", unique(description$Acronym),":H", sep =""),]
# 
# lapply(split(df_indi, df_indi$gene), function(i) {
#   a <- summary(lm(value ~ Acronym + condition + Acronym:condition, data = i))
#   p_val <- a$coefficients[,4]
#   return(p_val)
# }) -> all_lm
# 
# b <- matrix(unlist(all_lm), ncol = length(all_lm[[1]]), byrow = TRUE, dimnames = list(names(all_lm), names(all_lm[[1]])))
# b <- b[levene_pvals > 0.05, 37:70]
# sum(apply(b, 1, function(i) any(i < 0.05)))
# 
# corrected_lm <- p.adjust(b, method = "BH")
# corrected_lm <- matrix(corrected_lm, ncol = ncol(b), byrow = TRUE, dimnames = dimnames(b))
# aa <- apply(corrected_lm, 1, function(i) any(i < 0.05))
# length(names(aa)[aa])
# sum(names(aa)[aa] %in% candidates)
# candidates <- names(aa)[aa]
# 
# head(b)
# all_lm$ENSG00000132832
# names(all_lm)
# length(unlist(all_lm))
# corrected_lm <- p.adjust(unlist(all_lm), method = "BH")
# corrected_lm <- matrix(corrected_lm, ncol = length(all_lm[[1]]), byrow = TRUE, dimnames = list(names(all_lm), names(all_lm[[1]])))
# 
# aa <- apply(corrected_lm[,37:70], 1, function(i) any(i < 0.05))
# length(names(aa)[aa])
# sum(names(aa)[aa] %in% candidates)
# 
# 
## ANOVA old --------------
lapply(split(df_indi, df_indi$gene), function(i) {
  a <- summary(aov(value ~ Acronym + condition + Acronym:condition, data = i))
  p_val <- a[[1]]$`Pr(>F)`[3]
  return(p_val)
}) -> all_anovas
raw_p <- data.frame(unlist(all_anovas[levene_pvals > 0.05]))
sum(raw_p < 0.05)

# Multiple testing correction
correct_pvals <- p.adjust(all_anovas[levene_pvals > 0.05], method = "BH")
correct_pvals <- unlist(correct_pvals)
candidates <- names(correct_pvals)[correct_pvals < 0.05]

## Wilcoxon test for each gene_region, between H and X -------------------
# Subset the candidate genes only
inf.candidate <- subset(df_indi, gene %in% candidates)
inf.candidate <- droplevels(inf.candidate)

# Wilconxon test
## Define pairwise test for each gene,  with correction
pairwise_test_by_regions <- function(gene_extract) {
  lapply(split(gene_extract, gene_extract$Acronym), function(i) {
    wilcox.test(value ~ condition, data = i)$p.value
  })
}
# Apply for all genes
lapply(split(inf.candidate, inf.candidate$gene), function(extract) {
  pairwise_test_by_regions(extract)
}) -> signi_pairs
signi_pairs <- melt(signi_pairs)
colnames(signi_pairs) <- c("pval", "region", "gene")

# Which combination is significant?
sum(signi_pairs$pval < 0.05)
table(signi_pairs[signi_pairs$pval < 0.05, ]$region)
ggplot(data = signi_pairs, 
       mapping = aes(x = region, y = gene, fill = (pval<0.05))) +
  geom_tile() + scale_fill_manual(values = c("grey", "blue")) +
  theme(axis.text.y = element_blank()) +
  labs(fill = "Significant")


## Dysregulation counts and directions -------------------
up_or_down <- melt(t(mean_regions), varnames = c("gene", "region"))
up_or_down$upregulate <- up_or_down$value > 0 

final <- merge(signi_pairs, up_or_down, by = c("gene", "region"))
final <- final[final$pval < 0.05, ]
final$upregulate[final$upregulate == "TRUE"] <- "Up"
final$upregulate[final$upregulate == "FALSE"] <- "Down"
colnames(final) <- c("Gene", "Region", "pval", "Value", "Regulation")

updown_count <- final %>%
  group_by(Region, Regulation) %>%
  tally()
updown_count <- merge(updown_count, dend_list, by.x = "Region", by.y = "Acronym")
updown_count$n[updown_count$Regulation == "Down"] <- updown_count$n[updown_count$Regulation == "Down"]*(-1)
#x11()
ggplot(data = updown_count, aes(x = Region, y = n, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_abline(slope = 0, color = "white") +
  coord_flip() +
  labs(y = "Number of lincRNA genes")

# Combine with info about cluster
final <- merge(final, dend_list, by.x = "Region", by.y = "Acronym")
final$Region <- factor(final$Region, levels = dend_list$Acronym)

## Modules of genes ------------
df_profile <- up_or_down[up_or_down$gene %in% unique(final$Gene),]
profile <- dcast(df_profile, gene ~ region, value.var = "value")
rownames(profile) <- profile$gene
profile$gene <- NULL

modules_0 <- get_subdendrograms(make_dhc(t(individual_norm[candidates,]), 4), k = 4)
modules_0 <- lapply(modules_0, labels)
names(modules_0) <- as.roman(1:length(modules_0))
modules_0 <- melt(modules_0)
colnames(modules_0) <- c("gene", "module")

modules <- get_subdendrograms(make_dhc(t(profile), 6), k = 6)
modules <- lapply(modules, labels)
names(modules) <- as.roman(1:length(modules))

modules <- melt(modules)
colnames(modules) <- c("gene", "module")
table(modules$module)

# If use lm
df_profile <- up_or_down[up_or_down$gene %in% candidates,]
profile <- dcast(df_profile, gene ~ region, value.var = "value")
rownames(profile) <- profile$gene
profile$gene <- NULL

modules_0 <- get_subdendrograms(make_dhc(t(individual_norm[candidates,]), 4), k = 4)
modules_0 <- lapply(modules_0, labels)
names(modules_0) <- as.roman(1:length(modules_0))
modules_0 <- melt(modules_0)
colnames(modules_0) <- c("gene", "module")

modules <- get_subdendrograms(make_dhc(t(profile), 6), k = 6)
modules <- lapply(modules, labels)
names(modules) <- as.roman(1:length(modules))

modules <- melt(modules)
colnames(modules) <- c("gene", "module")
table(modules$module)

## DGE mRNA -----------
# Filter mRNA
mRNA <- sz_count[sz_count$Biotype == "protein_coding",]
rownames(mRNA) <- mRNA$X
mRNA <- mRNA[,description$sample]
filtered_mRNA <- filter_gene(mRNA, 2, description)
filtered_mRNA <- as.matrix(filtered_mRNA)

# Find DGE mRNA
mRNA_dds <- DESeqDataSetFromMatrix(countData = filtered_mRNA, 
                                   colData = description,
                                   design = ~ condition + Acronym)

mRNA_transformed <- varianceStabilizingTransformation(mRNA_dds,blind = TRUE)
mRNA_transformed <- assay(mRNA_transformed)

mRNA_RIN <- take_residuals(mRNA_transformed, "RIN")
mRNA_indi <- difference_to_brain_mean(mRNA_RIN, description)

# Test ...
# Prepare df
df_dge_mRNA <- melt(mRNA_indi)
names(df_dge_mRNA) <- c("gene", "sample", "value")
df_dge_mRNA <- merge(df_dge_mRNA, description, by = "sample")

## Levene test
lapply(split(df_dge_mRNA, df_dge_mRNA$gene), function(i) {
  a <- leveneTest(value ~ Acronym*condition, data = i)
  levene_p <- a$`Pr(>F)`[1]
  return(levene_p)
}) -> mRNA_levene_pvals

mRNA_levene_pvals <- data.frame(unlist(mRNA_levene_pvals))
sum(mRNA_levene_pvals > 0.05)

## New ANOVA
lapply(split(df_dge_mRNA, df_dge_mRNA$gene), function(i) {
  a <- summary(lm(value ~ Acronym + condition + Acronym:condition, data = i))
  p_val <- a$coefficients[,4]
  return(p_val)
}) -> mRNA_lm

b <- matrix(unlist(mRNA_lm), ncol = length(mRNA_lm[[1]]), byrow = TRUE, dimnames = list(names(mRNA_lm), names(mRNA_lm[[1]])))
b <- b[levene_pvals > 0.05, 37:70]
sum(apply(b, 1, function(i) any(i < 0.05)))

mRNA_pvals <- p.adjust(b, method = "BH")
mRNA_pvals <- matrix(mRNA_pvals, ncol = ncol(b), byrow = TRUE, dimnames = dimnames(b))
sum(apply(mRNA_pvals, 1, function(i) any(i < 0.05)))
length(names(aa)[aa])
sum(names(aa)[aa] %in% candidates)
candidates <- names(aa)[aa]

## ANOVA
lapply(split(df_dge_mRNA, df_dge_mRNA$gene), function(i) {
  a <- summary(aov(value ~ Acronym + condition + Acronym:condition, data = i))
  p_val <- a[[1]]$`Pr(>F)`[3]
  return(p_val)
}) -> mRNA_anovas

# Multiple testing correction
mRNA_pvals <- p.adjust(mRNA_anovas[mRNA_levene_pvals > 0.05], method = "BH")
mRNA_pvals <- unlist(mRNA_pvals)
dge_mRNA_list <- names(mRNA_pvals)[mRNA_pvals < 0.05]
rm(list = c("mRNA_pvals", "mRNA_levene_pvals", "mRNA_RIN", "mRNA_dds", "mRNA_transformed", "mRNA_anovas"))
# saveRDS(dge_mRNA_list, file = "mRNA_list.rds")
# 

dge_mRNA_list <- readRDS(file = "mRNA_list.rds")

dge_mRNA <- filtered_mRNA[dge_mRNA_list,]
dge_lincRNA <- filtered_lincRNA[candidates, ]
dim(dge_mRNA)

## mRNA Correlations ----------
## Without mean ---------------
associations <- matrix(nrow = length(candidates), ncol = length(dge_mRNA_list))
#correlation_pvals <- matrix(nrow = length(candidates), ncol = length(dge_mRNA_list))
rownames(associations) <- candidates
colnames(associations) <- dge_mRNA_list
candidates

for (i in 1:length(candidates)) {
  for (j in 1:length(dge_mRNA_list)) {
    correlation <- cor.test(individual_norm[candidates[i],], mRNA_indi[dge_mRNA_list[j],])
    associations[i,j] <- correlation$estimate
  }
}
#rownames(associations) <- candidates
#colnames(associations) <- dge_mRNA_list
#sum(p.adjust(p = correlation_pvals, method = "BH") < 0.01)
#length(correlation_pvals)

hist(associations)
# moduleI <- associations[modules$gene[modules$module == "I"],]
# hist(colMeans(moduleI))
dim(associations)
associations[trans_target[order(trans_target$trans_target),]]

trans_target <- apply(associations, 1, function(i) {sum(i > 0.8)})
trans_target <- as.data.frame(trans_target)
trans_target$gene <- rownames(trans_target)
trans_target <- trans_target[order(trans_target$trans_target, decreasing = T),]

top_target <- trans_target[1:5,] 

merge(top_target, modules, by.x = "row.names", by.y = "gene")
correlated_with_top_target <- colnames(associations)[(associations[top_target$gene[2],] > 0.8)]

# sum(correlated_with_top_target %in% mRNA_in_modules$IV)
# sum(mRNA_in_modules$IV %in% correlated_with_top_target)
# length(mRNA_in_modules$IV)
# length(correlated_with_top_target)

# Test one gene
ego <- enrichGO(gene          = correlated_with_top_target,
                universe      = dge_mRNA_list,
                keyType = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH")


ego <- enrichGO(gene          = rownames(filtered_mRNA),
                universe      = rownames(filtered_mRNA),
                keyType = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH")

dotplot(simplify(ego), showCategory = 20)
## With mean --------------
module_means <- aggregate(individual_norm[modules$gene,], list(modules$module), mean)
rownames(module_means) <- module_means$Group.1
module_means$Group.1 <- NULL

associations_1 <- matrix(nrow = nrow(module_means), ncol = length(dge_mRNA_list))
rownames(associations_1) <- rownames(module_means)
colnames(associations_1) <- dge_mRNA_list


for (i in 1:nrow(module_means)) {
  for (j in 1:length(dge_mRNA_list)) {
    correlation <- cor.test(unlist(module_means[i,]), mRNA_indi[dge_mRNA_list[j],])
    associations_1[i,j] <- correlation$estimate
  }
}


apply(associations_1, 1, function(i) sum(i > 0.8))

mRNA_in_modules <- apply(associations_1, 1, function(i) names(i)[i > 0.8])
enrichGO <- compareCluster(geneCluster = mRNA_in_modules, 
                           universe = dge_mRNA_list, 
                           fun = "enrichGO", 
                           keyType = "ENSEMBL", 
                           OrgDb = org.Hs.eg.db,
                           ont = "CC")
dotplot(enrichGO, showCategory = 8)
clusters_GO <- simplify(enrichGO)#, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(clusters_GO, showCategory = 7) + ggtitle("CC GO, dge genes")



# With all genes
associations_2 <- matrix(nrow = nrow(module_means), ncol = nrow(filtered_mRNA))
rownames(associations_2) <- rownames(module_means)
colnames(associations_2) <- rownames(filtered_mRNA)


for (i in 1:nrow(module_means)) {
  for (j in 1:nrow(filtered_mRNA)) {
    correlation <- cor.test(unlist(module_means[i,]), mRNA_indi[j,])
    associations_2[i,j] <- correlation$estimate
  }
}
apply(associations_2, 1, function(i) sum(i > 0.8))
hist(associations_2)
mRNA_in_modules_2 <- apply(associations_2, 1, function(i) names(i)[i > 0.8])
enrichGO2 <- compareCluster(geneCluster = mRNA_in_modules_2, 
                           universe = rownames(filtered_mRNA), 
                           fun = "enrichGO", 
                           keyType = "ENSEMBL", 
                           OrgDb = org.Hs.eg.db,
                           ont = "CC")
dotplot(enrichGO2) + ggtitle("CC GO")
dotplot(simplify(enrichGO2), showCategory = 8) + ggtitle("CC GO, all genes")

a <- enrichGO(mRNA_in_modules_2$III, universe = rownames(filtered_mRNA), OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC")
barplot(a)

# KEGG
mRNA_in_modules_converted <- lapply(mRNA_in_modules_2, function(i) {
  convert <- bitr(i, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(convert$ENTREZID)
  })
universe_converted <- bitr(rownames(filtered_mRNA), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cluster_KEGG <- compareCluster(geneCluster = mRNA_in_modules_converted, 
                            universe = universe_converted, 
                            fun = "enrichKEGG")

dotplot(cluster_KEGG)

cluster_Path <- compareCluster(geneCluster = mRNA_in_modules_converted, 
                               universe = universe_converted, 
                               fun = "enrichPathway")

dotplot(cluster_Path)


#ggplot(melt(individual_norm), aes(x = Var2, y = value)) + geom_boxplot(fill = factor(rep(1:8, each= 35)))

## Cell type markers ------------------
markers <- readRDS("markers_combined_cor_075_correct.rds")

table(markers$cell_type)

find_enriched_marker <- function(dge_genes, markers_list) {
  no.dge_markers <- sum(dge_genes %in% markers_list)
  no.markers <- length(markers_list)
  no.non_markers <- nrow(filtered_mRNA) - no.markers
  no.dge <- length(dge_genes)
  print(c(no.dge_markers, no.dge, no.markers, no.non_markers))
  prob <- dhyper(no.dge_markers, no.markers, no.non_markers, no.dge) + phyper(no.dge_markers, no.markers, no.non_markers, no.dge, lower.tail = F)
  return(prob)
}


cellmarker <- matrix(nrow = length(unique(markers$cell_type)), ncol = length(mRNA_in_modules_2))
colnames(cellmarker) <- names(mRNA_in_modules_2)
rownames(cellmarker) <- unique(markers$cell_type)

for (i in 1:length(unique(markers$cell_type))) {
  cell_type <- unique(markers$cell_type)[i]
  markers_list <- markers$gene[markers$cell_type == cell_type]
  for (j in 1:length(mRNA_in_modules_2)) {
    cellmarker[i,j] <- find_enriched_marker(mRNA_in_modules_2[[j]], markers_list)
  }
}

cellmarker_pval <- matrix(p.adjust(cellmarker, method = "BH"), nrow = nrow(cellmarker), ncol = ncol(cellmarker), dimnames = dimnames(cellmarker))
cellmarker_pval <- melt(cellmarker_pval)
names(cellmarker_pval) <- c("Type", "Module", "pvalue")
cellmarker_pval <- cellmarker_pval#[cellmarker_pval$Module %in% c("I","III","IV","VI"),]
ggplot(cellmarker_pval, aes(x = Module, y = Type, fill = pvalue)) + geom_tile()
cellmarker_pval

## With 4k dge genes only





## End ------------------
save.image("linc2.RData")
