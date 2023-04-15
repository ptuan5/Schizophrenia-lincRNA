rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("linc2.RData")

library(grid)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(vsn) # to plot shrinkage variance
library(gridExtra) # to have multiple plots in grids
library(cowplot)
library(dendextend) # to plot dendrogram
library(car) #Levene's test
library(showtext)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(EnsDb.Hsapiens.v79)
library(ggh4x)
library(ggdendro)
library(umap)

font_add_google("Open Sans")
showtext_auto()
basesize = 15 # Size for export figure
# library(plotly)
# library(EDASeq)

## Load data -------------
sz_count <- read.csv("C:/Users/chuot/Documents/Skoltech Study/Research/lncRNA/sz_human.counts_with.biotypes.csv", header = T)
# Filter lincRNA genes only
lincRNA <- sz_count[sz_count$Biotype == "lincRNA",]
lincRNA <- lincRNA[,1:ncol(lincRNA)-1]
# Move X to row names 
rownames(lincRNA) <- lincRNA$X
lincRNA$X <- NULL

mRNA <- sz_count[sz_count$Biotype %in% c("protein_coding", "lincRNA"),]
rownames(mRNA) <- mRNA$X
mRNA$X <- NULL
mRNA$CPM = log2(rowMeans(mRNA[,1:280])+1)
table(mRNA$Biotype)
ggplot(mRNA, aes(x = CPM, fill = Biotype)) + geom_density(alpha = .5)
gene_length <- getGeneLengthAndGCContent(rownames(lincRNA), "hsa")


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

## Transform + normalize for size factor + variance shrinkage --------------------
# Form DESeq Data Set (dds)
dds <- DESeqDataSetFromMatrix(countData = filtered_lincRNA, 
                              colData = description,
                              design = ~ condition + Acronym + condition:Acronym)

# Log 2 transformation of read counts including variance shrinkage
# rlog_data <- rlog(dds,blind = FALSE)
# saveRDS(object = rlog_data, file = "rlog_transformed.rds")
rlog_data <- readRDS("rlog_transformed.rds")
rlog_transformed <- assay(rlog_data)

msd_plot <- meanSdPlot(rlog_transformed,
                       ranks = FALSE,
                       plot = FALSE)
msd_plot$gg + labs(x = "Mean", y = "SD", title = "rlog-transformed") + theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))

pca <- prcomp(t(rlog_transformed))
pca



pca_plot <- function(df, description, group) {
  pca <- prcomp(t(df))
  var1 <- round((pca$sdev[1])^2/sum(pca$sdev^2)*100, 1)
  var2 <- round((pca$sdev[2])^2/sum(pca$sdev^2)*100, 1)
  plot <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, fill = description[,group])) +
    geom_point(shape = 21, colour = "black", size = 1.5) + #stroke = 1.5, size = 9) +
    labs(col = group, x = paste("PC1: ",var1, "% variance", sep = ""), y = paste("PC2: ", var2, "% variance", sep = "")) +
    theme_bw(base_size = 15, base_family = "Open Sans")
  return(plot)
}


# pca_plot(individual_norm, description, "big_area")

description[,"condition"]
pca_raw_condi <- pca_plot(rlog_transformed, description, "condition") +
  scale_fill_discrete(name = "Condition", labels= c("HC", "SZ")) + theme_bw()
pca_raw_region <- pca_plot(rlog_transformed, description, "big_area") + 
  scale_fill_discrete(name = "Brain Area") + theme_bw()
pca <- prcomp(t(rlog_transformed))
 
grid.arrange(pca_raw_condi, pca_raw_region, nrow = 1)

## Residuals of RIN normalization ------
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

grid.arrange(pca_plot(RIN_norm, description, "condition") + scale_fill_discrete(name = "Condition", labels= c("HC", "SZ")) + theme_bw(),pca_plot(RIN_norm, description, "big_area") + scale_fill_discrete(name = "Brain Area") + theme_bw(), nrow = 1)

## Individual value ------------------
difference_to_brain_mean <- function (matrix, description) {
  # Merge info
  df <- as.data.frame(t(matrix))
  df <- merge(df, description[,c('condition', 'Acronym', 'individuum')], by = "row.names")
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

# Reformat
indi_norm <- t(differ_from_mean[,genes])

individual_norm <- difference_to_brain_mean(RIN_norm, description)
rowSums(individual_norm)
# PCA
pca_indi_condi <- pca_plot(individual_norm, description, "condition") +
  scale_fill_discrete(name = "Condition", labels= c("HC", "SZ")) + labs(tag = "A") + 
  theme(legend.key.size = unit(1, "mm"), plot.tag = element_text(face = "bold"),  aspect.ratio = 1)


pca_indi_region <- pca_plot(individual_norm, description, "big_area") + 
  scale_fill_manual(name = "Section", values = colorspace::qualitative_hcl(10)) + labs(tag = "B") + 
  theme(legend.key.size = unit(1, "mm"), plot.tag = element_text(face = "bold"), aspect.ratio = 1)

grid.arrange(pca_indi_condi,pca_indi_region, nrow=1)

# a <- prcomp(t(individual_norm))
# plot_ly(as.data.frame(a$x), x = ~PC1, y = ~PC2, z = ~PC3, color = description[,"big_area"], type = "scatter3d", mode = "markers")

# description$check <- description$big_area %in% c("Frontal Lobe", "Occipital Lobe", "Parietal Lobe", "Temporal Lobe", "Insula")
# pca_plot(individual_norm, description, "check")

# Data frame
indi_df <- as.data.frame(t(individual_norm))
indi_df <- merge(indi_df, description, by = "row.names", all = T)


## Calculate fold change -------------
difference_by_regions <- function(matrix, description) {
  # Combine matrix with description
  df_indi <- as.data.frame(t(matrix))
  df_indi <- merge(df_indi, description, by = "row.names", all = T)  
  genes <- rownames(matrix)
  # Mean of condition_region
  df_indi %>%
    group_by(condition, Acronym) %>%
    dplyr::summarise_if(.predicate = is.numeric, .funs = mean) -> mean_regions
  # Difference in each region
  dif_by_regions <- mean_regions[mean_regions$condition == "X", genes] - mean_regions[mean_regions$condition == "H", genes]
  rownames(dif_by_regions) <- mean_regions$Acronym[mean_regions$condition == "H"]
  return(dif_by_regions)
}

foldchange_regions <- difference_by_regions(individual_norm, description)

## Hierarchical clustering ---------------
# Calculate mean by regions
indi_df %>%
  group_by(Acronym) %>%
  dplyr::summarise_if(.predicate = is.numeric, .funs = mean) -> mean_regions

mean_regions <- as.data.frame(mean_regions)
rownames(mean_regions) <- mean_regions$Acronym
mean_regions <- mean_regions[,rownames(individual_norm)]

# Create dendrogram
make_dhc <- function(matrix, k) {
  # Pearson correlation
  pearson.cor <- cor(matrix, method = "pearson")
  distance <- as.dist(1 - pearson.cor)
  # Get descending hierarchical clustering
  hc <- hclust(distance, method = "ward.D2")
  hc$height <- hc$height/5
  dhc <- as.dendrogram(hc)
  dhc <- color_branches(dhc, k)
  return(dhc)
}

# hclust alone -------
pearson.cor <- cor(t(mean_regions), method = "pearson")
distance <- as.dist(1 - pearson.cor)
# Get descending hierarchical clustering
hc <- hclust(distance, method = "ward.D2")
hc$height <- hc$height/5

labels(hc) <- maditr::vlookup(lookup_value = labels(hc), dict = abbre_region, result_column = 4, lookup_column = 3)

ggplot(abbre_region, aes(x = 1, y = area, fill = big_area, width = 0.01)) + 
  geom_tile() + scale_y_dendrogram(hclust = hc) + theme_minimal(base_size = basesize)

#-----
abc

dhc <- make_dhc(t(mean_regions), 3)

labels(dhc) <- maditr::vlookup(lookup_value = labels(dhc), dict = abbre_region, result_column = 4, lookup_column = 3)
ggplot(dhc, horiz = TRUE)

data <- dendro_data(dhc, type = "rectangle")
x <- ggplot(segment(data)) + geom_segment(aes(x=x, y=y, xend = xend, yend = yend)) + coord_flip() +
  scale_y_reverse(expand = c(0.2, 0))

# Rename
# abbre_region <- description[description$individuum == "HA", c("big_area","area","Acronym")]
abbre_region$area_character <- as.character(abbre_region$area)
labels(dhc) <- maditr::vlookup(lookup_value = labels(dhc), dict = abbre_region, result_column = 4, lookup_column = 3)
rownames(abbre_region) <- NULL
col_big_area <- factor(maditr::vlookup(lookup_value = labels(dhc), dict = abbre_region, result_column = 1, lookup_column = 2))

region_labels <- maditr::vlookup(lookup_value = labels(dhc), dict = abbre_region, result_column = 2, lookup_column = 3)
region_labels <- data.frame(labels(dhc),region_labels)
region_labels <- edit(region_labels)

abbre_region$area


par(mfrow=c(1,1))
# png(filename = "Global pattern.png", width = 1000, height = 1000)

layout(matrix(c(1,2, 3,3), nrow = 2, ncol = 2, byrow = TRUE))
pca_indi_condi
pca_indi_region
plot(dhc, horiz = TRUE, axes = FALSE, xlim = c(1,-0.8))
legend("topleft",title = "Anatomical Divisions", legend = levels(col_big_area), fill = colorspace::qualitative_hcl(10), cex = 2.5)  


abbre_region$area <- factor(abbre_region$area, levels = labels(dhc))
regions <- ggplot(abbre_region, aes(x = 0, y = area, fill = big_area)) + 
  geom_tile() + scale_y_dendrogram(hclust = dhc) 
grid.arrange(regions, dendro, nrow = 1)


labels(dhc) 
labels(dhc) <- region_labels$region_labels
col_big_area

labels_colors(dhc) <- colorspace::qualitative_hcl(10)[col_big_area]

labels_cex(dhc) <- .7
ggplot(dhc, horiz = TRUE)
dendro <- ggplot(dhc, horiz = TRUE) + theme_void(base_size = 18, base_family = "Open Sans") + ylim(1.1, -1.3) +
  theme(plot.tag = element_text(face="bold")) # axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
dendro
dendro + annotate(geom = "label", x = 30, y = 5, label = paste(unlist(levels(col_big_area)), collapse = "\n"), col = paste(unlist(colorspace::qualitative_hcl(10)), collapse = "\n"))

dendro #+ geom_text(data = label(dhc), aes(colour = col_big_area))

hcdata <- dendro_data(hc, type = "rectangle")
hcdata$labels <- merge(x = hcdata$labels, y = abbre_region, by.x = "label", by.y = "area")

jpeg("Paper/Global.jpg", width = 3000, height = 1500, quality = 750)
#grid.arrange(pca_indi_condi, pca_indi_region, dendro, layout_matrix = rbind(c(1,3), c(2,3)))
abc <- plot_grid(pca_indi_condi, pca_indi_region, ncol = 1,align = "v")
grid.arrange(abc,  dendro, ncol = 2)
dev.off()


tiff("Paper/Global.tiff", width = 272, height = 80, units = "mm", res = 200)
#grid.arrange(pca_indi_condi, pca_indi_region, dendro, layout_matrix = rbind(c(1,3), c(2,3)))
# abc <- plot_grid(pca_indi_condi, pca_indi_region, ncol = 2,align = "h")
abc
# plot_grid(abc,  dendro, ncol = 2, align = "h")
dev.off()

tiff("Paper/Dendro.tiff", width = 180, height = 180, units = "mm", res = 200)
#grid.arrange(pca_indi_condi, pca_indi_region, dendro, layout_matrix = rbind(c(1,3), c(2,3)))
# dendro
plot(dhc, horiz = TRUE, axes = FALSE, xlim = c(1,-0.6))
legend("topleft",title = "Anatomical Divisions", legend = levels(col_big_area), fill = colorspace::qualitative_hcl(10), cex = 0.7)

dev.off()




# Extract Cluster
dend_list <- get_subdendrograms(make_dhc(t(mean_regions), 3), k = 3)
dend_list <- lapply(dend_list, labels)
names(dend_list) <- as.roman(1:length(dend_list))

dend_list <- melt(dend_list)
colnames(dend_list) <- c("Acronym", "Cluster")
metadat_clust <- merge(x = description,  y = dend_list, by = "Acronym")
metadat_clust <- metadat_clust[order(metadat_clust$sample),]

tiff("Paper/Fig S1. Intersection of clusters and regions.tif", width = 170, height = 57, unit = "mm", res = 200)
ggplot(metadat_clust, aes(x = big_area, y = Cluster)) + geom_tile() + theme_bw(base_size = basesize, base_family = "Open Sans") + labs(x = "Anatomical Area")
dev.off()

## Levene's Test-------------
# Prepare df
df_indi <- reshape2::melt(individual_norm)
names(df_indi) <- c("gene", "sample", "value")
df_indi <- merge(df_indi, description, by = "sample")

#Actual test
lapply(split(df_indi, df_indi$gene), function(i) {
  a <- leveneTest(value ~ Acronym*condition, data = i)
  levene_p <- a$`Pr(>F)`[1]
  return(levene_p)
}) -> levene_pvals

levene_pvals <- data.frame(unlist(levene_pvals))
sum(levene_pvals > 0.05) # 754 genes
## ANOVA -----------
# Perform ANOVA test for all genes
lapply(split(df_indi, df_indi$gene), function(i) {
  a <- summary(aov(value ~ Acronym + condition + Acronym:condition, data = i))
  p_val <- a[[1]]$`Pr(>F)`[3]
  return(p_val)
}) -> all_anovas
raw_p <- data.frame(unlist(all_anovas[levene_pvals > 0.05]))
sum(raw_p < 0.05) # 234 genes

# Multiple testing correction
correct_pvals <- p.adjust(all_anovas[levene_pvals > 0.05], method = "BH")
correct_pvals <- unlist(correct_pvals)
significant_anova <- correct_pvals[correct_pvals < 0.05]
length(significant_anova) # 135 genes

names(significant_anova)
# Find pairs with Tukey
region_anova <- function(matrix) {
  a <- aov(value ~ Acronym:condition + condition + Acronym, data = matrix)
  aa <- TukeyHSD(a)
  names(aa$`Acronym:condition`)
  aaa <- as.data.frame(aa$`Acronym:condition`)
  aaa <- aaa[paste(unique(description$Acronym),":X-", unique(description$Acronym),":H", sep =""),]
  aaa$region <- unique(description$Acronym)
  aaa <- aaa[,c("region", "p adj")]
  rownames(aaa) <- NULL
  return(aaa)  
}
lapply(split(df_indi, df_indi$gene)[names(significant_anova)], region_anova) -> all_Tukey
df_Tukey <- reshape2::melt(all_Tukey, id = c("region", "p adj"))

signi_pairs <- df_Tukey[df_Tukey$`p adj` < 0.05,]
colnames(signi_pairs) <- c("region", "pval", "gene")
candidates <- unique(signi_pairs$gene)
table(signi_pairs$region)

# Which combination is significant?
table(signi_pairs$region)
ggplot(data = signi_pairs, 
       mapping = aes(x = region, y = gene, fill = (pval<0.05))) +
  geom_tile() + scale_fill_manual(values = c("grey", "blue")) +
  #theme(axis.text.y = element_blank()) +
  labs(fill = "Significant")

## Dysregulation counts and directions -------------------
# Decide up or down
foldchange_df <- melt(t(foldchange_regions), varnames = c("gene", "region"))
foldchange_df$up_down[foldchange_df$value > 0] <- "Up"
foldchange_df$up_down[foldchange_df$value < 0] <- "Down"
table(foldchange_df$up_down)

gene_symbols <- ensembldb::select(EnsDb.Hsapiens.v79, keys= as.character(foldchange_df$gene), keytype = "GENEID", columns = c("GENEID","SYMBOL"))
foldchange_df <- merge(foldchange_df, gene_symbols, by.x = "gene", by.y = "GENEID", all = TRUE)
foldchange_df$SYMBOL[is.na(foldchange_df$SYMBOL)] <- as.character(foldchange_df$gene[is.na(foldchange_df$SYMBOL)])

saveRDS(foldchange_regions, file = "foldchange.rds")


# Merge all info
final <- merge(signi_pairs, foldchange_df, by = c("gene", "region"))
final <- merge(final, dend_list, by.x = "region", by.y = "Acronym")
final$region <- factor(final$region, levels = dend_list$Acronym)
write.csv(x = final, file = "DE lincRNA.csv")

table(final$gene[final$value > -1 & final$value < 1])

max_foldchanges <- apply(X = t(foldchange_regions), MARGIN = 1, FUN = function(x) {x[which.max(abs(x))]})

volcano_data <- merge(reshape2::melt(correct_pvals), reshape2::melt(max_foldchanges), by = "row.names")
colnames(volcano_data) <- c("gene", "pval", "log2FC")
volcano_data$pval <- -log(volcano_data$pval, 10)

volcano_data$note[volcano_data$gene %in% candidates] <- "DEL"
volcano_data$note[!(volcano_data$gene %in% candidates)] <- "Not DEL"

volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = pval, col = note)) + geom_point(size = 1.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = "blue", lty = "dashed") + scale_x_continuous(limits = c(-4.1,4.1)) +
  labs(x = expression(paste("Maximum ", log[2], "(fold-change)")), y = expression(paste("\u2013", log[10], "(p\u00advalue)")), col = NULL) +
  theme_bw(base_family = "Open Sans", base_size = 15) + labs(tag = "A") +
  theme(legend.position = "top", plot.tag = element_text(face = "bold")) 
volcano_plot
 

merge(volcano_data, max_foldchanges, by = "row.names") 
# Count
updown_count <- final %>%
  group_by(region, up_down) %>%
  tally()
updown_count
updown_count <- merge(updown_count, dend_list, by.x = "region", by.y = "Acronym")
updown_count$n[updown_count$up_down == "Down"] <- updown_count$n[updown_count$up_down == "Down"]*(-1)
#x11()
updown_plot <- ggplot(data = updown_count, aes(x = region, y = n, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_abline(slope = 0, color = "black") +
  coord_flip() + scale_y_continuous(labels = abs, limits = c(-15, 15)) +
  labs(x = "Region", y = "Number of DELs") + 
  annotate(geom="label", x = 10, y = -10, label = "Down") + #, size = 18) +
  annotate(geom="label", x = 10, y = 10, label = "Up") + #, size = 18) +
  theme_bw(base_family = "Open Sans", base_size = 15) + labs(tag = "B") + 
  theme(plot.tag = element_text(face = "bold"))
updown_plot

font_add_google(name = "Press Start 2P")
grid.arrange(volcano_plot, updown_plot, nrow = 1)
## Plot individual profile ----
indi_df$Acronym <- factor(indi_df$Acronym, levels = dend_list$Acronym)
indi_df <- merge(indi_df, dend_list, by = "Acronym")

foldchange_df$region <- factor(foldchange_df$region, levels = dend_list$Acronym)
foldchange_df <- merge(foldchange_df, dend_list, by.x = "region", by.y = "Acronym")

 
chosen_candidates <- c("MEG3", "MEG9", "MALAT1", "LINC01252", "RP11-247L20.4")

five_prof <- ggplot(foldchange_df[foldchange_df$SYMBOL %in% chosen_candidates, ], 
       aes(x = region, y = value, group = gene, col = SYMBOL)) + geom_line() +
  # facet_grid(~Cluster, scales = "free_x", space = "free_x")
  geom_vline(xintercept = c(20.5,26.5), linetype = "dashed") + #, size = 2.4) + 
  geom_hline(yintercept = 0, col = "grey") + #, size = 1.8) +
  theme_bw(base_family = "Open Sans", base_size = 15) + 
  labs(x = "Region", y = expression(paste(log[2], "(SZ/HC)")), col = "Gene") +
  theme(legend.position = c(0.1, 0.77),
        legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(face = "italic"), 
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.3), 
        plot.tag = element_text(face = "bold")) + labs(tag = "C")

five_prof


tiff("Paper/DEL.tiff", width = 170, height = 153, units = "mm", res = 200)
grid.arrange(volcano_plot, updown_plot, five_prof, 
             layout_matrix = rbind(c(1,2),c(3,3)))                    
dev.off()


ggplot(foldchange_df[foldchange_df$SYMBOL %in% chosen_candidates, ], 
       aes(x = region, y = value, group = gene, col = SYMBOL))
table(dend_list$Cluster)
## mRNA Correlations ----------
both_indi <- readRDS("both_indi.rds")

associations <- matrix(nrow = length(candidates), ncol = nrow(filtered_mRNA))
rownames(associations) <- candidates
colnames(associations) <- rownames(filtered_mRNA)

for (i in 1:length(candidates)) {
  for (j in 1:nrow(filtered_mRNA)) {
    associations[i,j] <- cor(both_indi[candidates[i],], both_indi[rownames(filtered_mRNA)[j],])
  }
}

get_top_targets <- function(associations, min_r, no_target) {
  trans_target <- apply(associations, 1, function(i) {sum(i > min_r) + sum(i < -min_r)})
  trans_target <- as.data.frame(trans_target)
  trans_target$gene <- rownames(trans_target)
  trans_target <- trans_target[order(trans_target$trans_target, decreasing = T),]
  
  top_target_list <- rownames(trans_target)[1:no_target]
  top_targets <- associations[top_target_list,]
  
  lapply(split(top_targets, top_target_list), function(i) {
    genes <- colnames(top_targets)[i > min_r | i < -min_r]
    return(genes)
  }) -> targets
  return(targets)
}


top_targets <- get_top_targets(associations,0.85,3)
top_targets
names(top_targets) <- maditr::vlookup(lookup_value = names(top_targets), dict = gene_symbols, result_column = 2, lookup_column = 1)

chosen <- c("ENSG00000260328", "ENSG00000260804", "ENSG00000277200")
chosen <- t(both_indi[chosen,])
typeof(chosen)
cor(chosen)
chosen <- as.data.frame(chosen)

sum(top_targets$PKI55 %in% top_targets$`RP11-74E22.8`)
length(top_targets)
names(top_targets)[2] <- "LINC01963"

cluster_GO_BP <- compareCluster(geneClusters = top_targets, 
                             universe = rownames(filtered_mRNA), 
                             fun = "enrichGO",
                             keyType = "ENSEMBL",
                             OrgDb = org.Hs.eg.db,
                             ont = "BP")

dotplot(simplify(cluster_GO_BP)) + labs(x = "Differentially Expressed lincRNAs", col = "Adjusted P-value", size = "Gene Ratio") #+
  ggtitle("GO Cellular Compartment Terms")
0
table(cluster_GO_BP@compareClusterResult$ID) == 3
View(cluster_GO_CC@compareClusterResult)

dge_mRNA <- readRDS("mRNA_dge.rds")


lapply(top_targets, function(i) i[i %in% dge_mRNA]) -> dge_targets ## so low :((
dge_targets

cluster_GO <- compareCluster(geneClusters = dge_targets,
                             universe = dge_mRNA,
                             keyType = "ENSEMBL",
                             fun = "enrichGO",
                             OrgDb = org.Hs.eg.db,
                             ont = "MF")
dotplot(simplify(cluster_GO))
length(unique(cluster_GO@compareClusterResult$ID))

# KEGG
entrezid_converted <- lapply(top_targets, function(i) {
  convert <- bitr(i, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(convert$ENTREZID)
})
universe_converted <- bitr(rownames(filtered_mRNA), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cluster_KEGG <- compareCluster(geneCluster = entrezid_converted, 
                               universe = universe_converted$ENTREZID, 
                               fun = "enrichKEGG")

c <- dotplot(cluster_KEGG) + ggtitle("KEGG Pathways")

# Reactome
cluster_Path <- compareCluster(geneCluster = entrezid_converted, 
                               universe = universe_converted$ENTREZID, 
                               fun = "enrichPathway")

d <- dotplot(cluster_Path) + ggtitle("Reactome")
d
a
## Modules of genes ------------
modules <- get_subdendrograms(make_dhc(foldchange_regions[,candidates], 6), k = 6)
modules <- lapply(modules, labels)
names(modules) <- as.roman(1:length(modules))

modules <- reshape2::melt(modules)
colnames(modules) <- c("gene", "module")
table(modules$module)

# Module Means
module_means <- aggregate(both_indi[modules$gene,], list(modules$module), mean)
rownames(module_means) <- module_means$Group.1
module_means$Group.1 <- NULL

# Find associations
associations_modules <- matrix(nrow = nrow(module_means), ncol = nrow(filtered_mRNA))
rownames(associations_modules) <- rownames(module_means)
colnames(associations_modules) <- rownames(filtered_mRNA)

for (i in 1:nrow(module_means)) {
  for (j in 1:nrow(filtered_mRNA)) {
    associations_modules[i,j] <- cor(unlist(module_means[i,]), both_indi[rownames(filtered_mRNA)[j],])
  }
}

apply(associations_modules, 1, function(i) sum(i > 0.85))
mRNA_in_modules <- apply(associations_modules, 1, function(i) names(i)[i > 0.85])

# Enrichment analysis
enrichGO <- enrichGO(gene = mRNA_in_modules$III,
                     universe = rownames(filtered_mRNA),
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     ont = "CC")
dotplot(simplify(enrichGO))
barplot(enrichGO)

# KEGG
mRNA_in_modules_converted <- lapply(mRNA_in_modules, function(i) {
  convert <- bitr(i, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(convert$ENTREZID)
  })

modules_KEGG <- enrichKEGG(gene = mRNA_in_modules_converted$III,
                           universe = universe_converted$ENTREZID)
barplot(modules_KEGG)

modules_Path <- enrichPathway(gene = mRNA_in_modules_converted$III,
                              universe = universe_converted$ENTREZID)

barplot(modules_Path)

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


cellmarker <- matrix(nrow = length(unique(markers$cell_type)), ncol = length(top_targets))
colnames(cellmarker) <- names(top_targets)
rownames(cellmarker) <- unique(markers$cell_type)

for (i in 1:length(unique(markers$cell_type))) {
  cell_type <- unique(markers$cell_type)[i]
  markers_list <- markers$gene[markers$cell_type == cell_type]
  for (j in 1:length(top_targets)) {
    cellmarker[i,j] <- find_enriched_marker(top_targets[[j]], markers_list)
  }
}


top_targets$`RP11-74E22.8`[top_targets$`RP11-74E22.8` %in% markers$gene[markers$cell_type == "In"]]

cellmarker_pval <- matrix(p.adjust(cellmarker, method = "BH"), nrow = nrow(cellmarker), ncol = ncol(cellmarker), dimnames = dimnames(cellmarker))
cellmarker_pval <- reshape2::melt(cellmarker_pval)
names(cellmarker_pval) <- c("Type", "Module", "pvalue")
cellmarker_pval <- cellmarker_pval#[cellmarker_pval$Module %in% c("I","III","IV","VI"),]
levels(cellmarker_pval$Module) <- c("RP11-74E22.8", "LINC01963", "RP11-416I2.1")


GO_BP_plot <- dotplot(simplify(cluster_GO_BP),) + 
  theme_bw(base_family = "Open Sans", base_size = 11) + 
  scale_size_area(max_size = 5) +
  labs(x = "DE lincRNAs", col = "Adjusted P-value", size = "Gene Ratio", tag = "A") + 
  theme(legend.key.size = unit(3, 'mm'), axis.text = element_text(lineheight = 0.8),
        plot.tag = element_text(face="bold"))
GO_BP_plot
cell_marker_plot <- ggplot(cellmarker_pval, aes(x = Type, y = Module)) + geom_tile(aes(fill = pvalue)) + 
  geom_text(aes(label = round(pvalue, 2)), size = 5) + 
  scale_fill_gradient(low = "orangered", high = "gray97") + 
  labs(x = "Cell Type", y = "DE lincRNAs", fill = "P-value", tag = "B") + 
  theme_bw(base_family = "Open Sans", base_size = basesize) + 
  theme(legend.key.size = unit(3, 'mm'), plot.tag = element_text(face = "bold")) 

cell_marker_plot

tiff("Paper/Enrichment.tiff", width = 170, height = 136, units = "mm", res = 200)
grid.arrange(GO_BP_plot, cell_marker_plot, heights = c(0.7, 0.3))
dev.off()
# basesize = 10
# maxsize = 3
GO_CC_plot <- dotplot(simplify(cluster_GO_CC),) + 
  theme_bw(base_family = "Open Sans", base_size = basesize) + 
  scale_size_area(max_size = maxsize) +
  labs(x = "DE lincRNAs", col = "Adjusted P-value", size = "Gene Ratio", tag = "A") + 
  theme(legend.key.size = unit(3, 'mm'), axis.text = element_text(lineheight = 0.4),
        plot.tag = element_text(face="bold"))
GO_MF_plot <- dotplot(simplify(cluster_GO_MF),) + 
  theme_bw(base_family = "Open Sans", base_size = basesize) + 
  scale_size_area(max_size = maxsize) +
  labs(x = "DE lincRNAs", col = "Adjusted P-value", size = "Gene Ratio", tag = "B") + 
  theme(legend.key.size = unit(3, 'mm'), axis.text = element_text(lineheight = 0.4),
        plot.tag = element_text(face="bold"))
KEGG_plot <- dotplot(cluster_KEGG) + 
  theme_bw(base_family = "Open Sans", base_size = basesize) + 
  scale_size_area(max_size = maxsize) +
  labs(x = "DE lincRNAs", col = "Adjusted P-value", size = "Gene Ratio", tag = "C") + 
  theme(legend.key.size = unit(3, 'mm'), axis.text = element_text(lineheight = 0.4),
        plot.tag = element_text(face="bold"))
Path_plot <- dotplot(cluster_Path) + 
  theme_bw(base_family = "Open Sans", base_size = basesize) + 
  scale_size_area(max_size = maxsize) +
  labs(x = "DE lincRNAs", col = "Adjusted P-value", size = "Gene Ratio", tag = "D") + 
  theme(legend.key.size = unit(3, 'mm'), axis.text = element_text(lineheight = 0.4),
        plot.tag = element_text(face="bold"))

tiff("Paper/Fig S3. Enrich analyses.tiff", width = 170, height = 120, units = "mm", res = 200)
plot_grid(GO_CC_plot, GO_MF_plot, KEGG_plot, Path_plot, ncol = 2, align = "hv")
dev.off()


saveRDS(object = candidates, file = "lincRNA_DE.rds")

## Reviewer's questions-----------------
mapped_stat <- read.csv("contr_sz_meta.csv")
p1 <- ggplot(mapped_stat, aes(x = characteristics..individual, y = total.reads.count, fill = description)) + geom_boxplot() +
  labs(x = "Individual", y = "Count") +
  scale_fill_discrete(name = "Condition", labels= c("HC", "SZ")) +
  scale_x_discrete(labels = c("H1", "H2", "H3", "H4", "X1", "X2", "X3", "X4")) +
  theme_bw()
(p2 <- ggplot(mapped_stat, aes(x = Acronym, y = total.reads.count, fill = description)) + geom_boxplot() +
  labs(x = "Region", y = "Count") +
  scale_fill_discrete(name = "Condition", labels= c("HC", "SZ")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.3)))
grid.arrange(p1,p2, ncol = 1, heights = c(.4, .6), top = "Total reads count")

ggplot(mapped_stat, aes(x = alignment.rate, fill = characteristics..individual)) + geom_density(alpha = .5) +
labs(x = "Alignment Rate (%)", y = "Density", title = "Mapping alignment rate") +
scale_fill_discrete(name = "Individual", labels= c("H1", "H2", "H3", "H4", "X1", "X2", "X3", "X4")) +
theme_bw()

#UMAP
set.seed(1234)
indi.umap <- t(individual_norm) %>%
  scale() %>%
  umap()

indi.umap$layout %>%
  as.data.frame() %>%
  rename(UMAP1 = "V1", UMAP2 = "V2") -> indi.umap.df


grid.arrange(ggplot(indi.umap.df, aes(x = UMAP1, y =UMAP2, fill = description$condition)) + geom_point(shape = 21, colour = "black",stroke = 0.5, size = 3) + scale_fill_discrete(name = "Brain Area") + theme_bw(),
             ggplot(indi.umap.df, aes(x = UMAP1, y =UMAP2, fill = description$big_area)) + geom_point(shape = 21, colour = "black",stroke = 0.5, size = 3) + scale_fill_discrete(name = "Brain Area") + theme_bw(),
             ncol = 2, top = "After removing batch effect")

## End ------------------
save.image("linc2.RData")
