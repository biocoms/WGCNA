###Install gene expression data set from TGCA

library("WGCNA")
library("tidyverse")
library("dendextend")
library("gplots")
library("ggplot2")
library("ggpubr")
library("VennDiagram")
library("dplyr")
library("GO.db")
library("DESeq2")
library("genefilter")
library("clusterProfiler")
library("org.Hs.eg.db")

## Load tumor and normal data
tumor_data <- read.csv("GeneExpression/OSCC_TCGA_gene_expression_337.csv", row.names = 1)
normal_data <- read.csv("GeneExpression/OSCC_TCGA_gene_expression_32.csv", row.names = 1)

# View data for 5 rows and 10 columns of tumor data
tumor_data[1:5,1:10] 
# View data for 5 rows and 10 columns of normal data
normal_data[1:5,1:10]

# ## remove the series number in gene ID
# rownames(tumor_data) <- sapply(strsplit(rownames(tumor_data), "\\."), '[',1)
# rownames(normal_data) <- sapply(strsplit(rownames(normal_data), "\\."), '[',1)
# 
# ##### Remove expression estimates with counts in less than 20% of case to robust the analysis
# tumor_data = tumor_data[apply(tumor_data,1,function(x) sum(x==0))<ncol(tumor_data)*0.8,]
# normal_data = normal_data[apply(normal_data,1,function(x) sum(x==0))<ncol(normal_data)*0.8,]

### Normalize expression raw count with DESeq ##########
# Prepare metadata
metadata_tumor <- data.frame(Sample = colnames(tumor_data),
                             Condition = rep(c('Tumor'),337))
metadata_normal <- data.frame(Sample = colnames(normal_data),
                             Condition = rep(c('Normal'),32))

###############Normalize Data#########################
#Tumor
# 1. Construct the DESeqDataSet object for the tumor data
dds_tumor <- DESeqDataSetFromMatrix(countData = round(tumor_data), colData = metadata_tumor,design = ~1 )
# 2. Run the DESeq function to estimate size factors
dds_tumor <- DESeq(dds_tumor)
# 3. Extract the normalized counts matrix
normalized_counts_tumor <- counts(dds_tumor, normalized = TRUE)

# Repeat the entire process for normal_data
#Normal
dds_normal <- DESeqDataSetFromMatrix(countData = round(normal_data), colData = metadata_normal, design = ~1)
dds_normal <- DESeq(dds_normal)
normalized_counts_normal <- counts(dds_normal, normalized = TRUE)


########Filter the low variance#######
#Variance Stabilization
vsd_tumor <- varianceStabilizingTransformation(dds_tumor)
vsd_normal <- varianceStabilizingTransformation(dds_normal)

# Retain genes with high variance
rv_tumor <- rowVars(assay(vsd_tumor))
rv_normal <- rowVars(assay(vsd_normal))
q95_tumor <- quantile(rv_tumor, 0.95)
q95_normal <- quantile(rv_normal, 0.95)
filtered_tumor <- assay(vsd_tumor)[rv_tumor > q95_tumor, ]
filtered_normal <- assay(vsd_normal)[rv_normal > q95_normal, ]

# ##Violin plot for normalized and 95 quantile Expression
# normalized_counts_tumor_long <- as.data.frame(normalized_counts_tumor[,1:30]) %>%
#   rownames_to_column(var = "Gene") %>%  # If genes are row names
#   pivot_longer(
#     cols = -Gene,                      # All columns except 'Gene' (sample IDs)
#     names_to = "name",                 # Column for sample names
#     values_to = "value"                # Column for expression values
#   )
# 
# ggplot(normalized_counts_tumor_long, aes(x = name, y = value)) +
#   geom_violin() +
#   geom_point(alpha = 0.5) +  # Optional: Add transparency to points
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1)  # Adjust axis labels
#   ) +
#   labs(
#     title = "Normalized and 95 Quantile Expression",
#     x = "Treatment",
#     y = "Normalized Expression"
#   )

################################# Perform WGCNA ########################################## 
##Choose softthreshold
#Enable multithreads
allowWGCNAThreads()

#  Choose soft threshold parameter
# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 10, to=20, by=2))

# Tumor
sft_tumor <- pickSoftThreshold(t(filtered_tumor), powerVector = powers, verbose = 5)
#print sft_tumor
sft_tumor

# Normal
sft_normal <- pickSoftThreshold(t(filtered_normal), powerVector = powers, verbose = 5)
#print sft_normal
sft_normal

# Choose power values (e.g., based on R^2 > 0.9)
power_tumor <- 3
power_normal <- 14
# Scale-free topology fit index as a function of the soft-thresholding power
#Plotting the results
par(mfrow = c(1,2))
cex1 = 0.9

#Index the scale free topology adjust as a function of the power soft thresholding - Tumor
plot(sft_tumor$fitIndices[,1], -sign(sft_tumor$fitIndices[,3])*sft_tumor$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft_tumor$fitIndices[,1], -sign(sft_tumor$fitIndices[,3])*sft_tumor$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.9,col="red")

#Connectivity mean as a function of soft power thresholding - Tumor
plot(sft_tumor$fitIndices[,1], sft_tumor$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_tumor$fitIndices[,1], sft_tumor$fitIndices[,5], labels=powers, cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.9,col="red")

#Index the scale free topology adjust as a function of the power soft thresholding - Normal
plot(sft_normal$fitIndices[,1], -sign(sft_normal$fitIndices[,3])*sft_normal$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft_normal$fitIndices[,1], -sign(sft_normal$fitIndices[,3])*sft_normal$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.9,col="red")

#Connectivity mean as a function of soft power thresholding -Normal
plot(sft_normal$fitIndices[,1], sft_normal$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_normal$fitIndices[,1], sft_normal$fitIndices[,5], labels=powers, cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.9,col="red")

#2.	Construct the co-expression networks
power_tumor = 3
power_normal = 14

# set cor to ensure cor() are used from WGCNA instead of the cor() generic in incase you get error of unused argument error
cor <- WGCNA::cor 

# Tumor
netwk_tumor <- blockwiseModules(
  t(filtered_tumor),
  power = power_tumor,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "tumor",
  verbose = 3
)
mergedColors = labels2colors(netwk_tumor$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk_tumor$dendrograms[[1]],
  mergedColors[netwk_tumor$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# Normal
netwk_normal <- blockwiseModules(
   t(filtered_normal),
  power = power_normal,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal",
  verbose = 3
)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk_normal$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk_normal$dendrograms[[1]],
  mergedColors[netwk_normal$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# check number of genes in each module
#Tumor
table(netwk_tumor$colors)
#Normal
table(netwk_normal$colors)

#
#===============================================================================
# Tumor eigengenes
MEs_tumor <- moduleEigengenes(t(filtered_tumor), colors = labels2colors(netwk_tumor$colors))$eigengenes
MEs_tumor <- orderMEs(MEs_tumor)
colnames(MEs_tumor) = names(MEs_tumor) %>% gsub("ME","", .)

# Normal eigengenes
MEs_normal <- moduleEigengenes(t(filtered_normal), colors = labels2colors(netwk_normal$colors))$eigengenes
MEs_normal <- orderMEs(MEs_normal)
colnames(MEs_normal) = names(MEs_normal) %>% gsub("ME","", .)

# Correlation heatmap of module eigengenes for both conditions
cor_tumor <- cor(MEs_tumor)
cor_normal <- cor(MEs_normal)

# png(
#   filename = paste0("image/tumor/module_eigen_cor.png"),
#   width = 200,
#   height = 200,
#   res = 200
# )
heatmap.2(cor_tumor,
          main = "Tumor Module Eigengene Correlation",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50), # Color gradient
          key = TRUE,  # Adds a color key (legend)
          key.title = "Correlation",
          key.xlab = "Value",
          density.info = "none",  # Disable histogram in the legend
          denscol = NA,
          cexRow = 1.0,         # Adjust text size for rows
          cexCol = 1.0)           # Remove histogram color)
# dev.off()

# Normal heatmap with legend
heatmap.2(cor_normal,
          main = "Normal Module Eigengene Correlation",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50), # Color gradient
          key = TRUE,  # Adds a color key (legend)
          key.title = "Correlation",
          key.xlab = "Value",
          density.info = "none",  # Disable histogram in the legend
          denscol = NA,
          cexRow = 1.0,         # Adjust text size for rows
          cexCol = 1.0)           # Remove histogram color)

#Generate edge list
# Tumor TOM and edge list
TOM_tumor <- TOMsimilarityFromExpr(t(filtered_tumor), power = power_tumor)
row.names(TOM_tumor) <- row.names(filtered_tumor)
colnames(TOM_tumor) <- row.names(filtered_tumor)

# Convert the TOM matrix into a long-format edge list
edge_list_tumor <- data.frame(TOM_tumor) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
  filter(gene1 != gene2) %>%
  filter(correlation > 0.1)

# Remove duplicate edges by considering them undirected (keeping only one direction)
edge_list_tumor <- edge_list_tumor %>%
  group_by(gene1) %>%
  slice_max(order_by = correlation, n = 5)

library(org.Hs.eg.db)
#Convert to symbol for gene_pair data and expresion data
edge_list_tumor$gene1.name <- mapIds(org.Hs.eg.db, keys = edge_list_tumor$gene1, column = "SYMBOL", keytype = "ENSEMBL")
edge_list_tumor$gene2.name <- mapIds(org.Hs.eg.db, keys = edge_list_tumor$gene2, column = "SYMBOL", keytype = "ENSEMBL")

edge_list_tumor.cyto <- data.frame(gene1 = edge_list_tumor$gene1.name, gene2 = edge_list_tumor$gene2.name , value = edge_list_tumor$correlation)
edge_list_tumor.cyto <- na.omit(edge_list_tumor.cyto)
write.table(edge_list_tumor.cyto, "edge_list_tumor.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

module_tumor <- data.frame(gene = names(netwk_tumor$colors), color = labels2colors(netwk_tumor$colors))
module_tumor$genename <- mapIds(org.Hs.eg.db, keys = module_tumor$gene, column = "SYMBOL", keytype = "ENSEMBL")
module_tumor <- module_tumor[module_tumor$genename %in% c(unique(edge_list_tumor.cyto$gene1), unique(edge_list_tumor.cyto$gene2)),]

write.table(data.frame(genename = module_tumor$genename, value = module_tumor$color), "node_list_tumor.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)


# Normal TOM and edge list
TOM_normal <- TOMsimilarityFromExpr(t(filtered_normal), power = power_normal)
row.names(TOM_normal) <- row.names(filtered_normal)
colnames(TOM_normal) <- row.names(filtered_normal)

edge_list_normal <- data.frame(TOM_normal) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
  filter(gene1 != gene2) %>%
  filter(correlation > 0.3)

edge_list_normal <- edge_list_normal %>%
  group_by(gene1) %>%
  slice_max(order_by = correlation, n = 5)

edge_list_normal$gene1.name <- mapIds(org.Hs.eg.db, keys = edge_list_normal$gene1, column = "SYMBOL", keytype = "ENSEMBL")
edge_list_normal$gene2.name <- mapIds(org.Hs.eg.db, keys = edge_list_normal$gene2, column = "SYMBOL", keytype = "ENSEMBL")

edge_list_normal.cyto <- data.frame(gene1 = edge_list_normal$gene1.name, gene2 = edge_list_normal$gene2.name , value = edge_list_normal$correlation)
edge_list_normal.cyto <- na.omit(edge_list_normal.cyto)
write.table(edge_list_normal.cyto, "edge_list_normal.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

module_normal <- data.frame(gene = names(netwk_normal$colors), color = labels2colors(netwk_normal$colors))
module_normal$genename <- mapIds(org.Hs.eg.db, keys = module_normal$gene, column = "SYMBOL", keytype = "ENSEMBL")
module_normal <- module_normal[module_normal$genename %in% c(unique(edge_list_normal.cyto$gene1), unique(edge_list_normal.cyto$gene2)),]

write.table(data.frame(genename = module_normal$genename, value = module_normal$color), "node_list_normal.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)


################ Module Preservation and Reproducibility ##################################

# Extract module assignments and colors
tumor_modules <- netwk_tumor$colors
tumor_module_colors <- labels2colors(tumor_modules)
names(tumor_module_colors) <- names(netwk_tumor$colors)

normal_module_colors <- labels2colors(netwk_normal$colors)
names(normal_module_colors) <- names(netwk_normal$colors)

# Prepare input data for module preservation
multiData <- list(
  Tumor = list(data = t(filtered_tumor)),   # Transpose to make samples columns
  Normal = list(data = t(filtered_normal)) # Transpose to make samples columns
)

# Extract module assignments and colors
tumor_modules <- netwk_tumor$colors
tumor_module_colors <- labels2colors(tumor_modules)
names(tumor_module_colors) <- names(netwk_tumor$colors)

normal_module_colors <- labels2colors(netwk_normal$colors)
names(normal_module_colors) <- names(netwk_normal$colors)

# Module colors for tumor and normal datasets
multiColor <- list(
  Tumor = tumor_module_colors,   # Named vector of module colors for tumor
  Normal = normal_module_colors  # Named vector of module colors for normal
)
all(names(multiColor$Tumor) %in% rownames(filtered_tumor))  # Should return TRUE
all(names(multiColor$Normal) %in% rownames(filtered_normal))  # Should return TRUE

preservation_results <- modulePreservation(
  multiData = multiData,          # List of datasets
  multiColor = multiColor,         # List of module colors
  referenceNetworks = 1,          # Use tumor as reference (index in multiData)
  nPermutations = 100,            # Number of permutations (higher for real analysis, e.g., 200)
  randomSeed = 12345,             # For reproducibility
  verbose = 3                     # Verbose output
)

# Preservation statistics for modules
preservation_stats <- preservation_results$preservation$Z$ref.Tumor$inColumnsAlsoPresentIn.Normal
#The presence of the "gold" module in your module preservation analysis is a feature of the modulePreservation function in WGCNA. 
#The "gold" module is not a biological module but rather an artificial module included as a reference and it should not be included in the functional enrichment analysis.
preservation_stats <- preservation_stats[rownames(preservation_stats) != "gold", ]

# Plot Z-summary statistics
mod_colors <- rownames(preservation_stats)  # Module colors
Z_summary <- preservation_stats$Zsummary.pres

barplot(
  Z_summary, names.arg = mod_colors,
  col = mod_colors, las = 2,
  ylab = "Preservation Z-summary",
  main = "Module Preservation Statistics"
)
abline(h = 2, col = "blue", lty = 2)  # Moderate preservation threshold
abline(h = 10, col = "red", lty = 2)  # High preservation threshold'

# Identify highly preserved modules
highly_preserved <- rownames(preservation_stats)[preservation_stats$Zsummary.pres > 10]

print(highly_preserved)  # List of module colors
module_names <- unique(tumor_module_colors)

# Identify moderate preserved modules
moderately_preserved <- rownames(preservation_stats)[preservation_stats$Zsummary.pres > 2 & preservation_stats$Zsummary.pres < 10]
print(moderately_preserved)

# Identify low preserved modules
low_preserved<- rownames(preservation_stats)[preservation_stats$Zsummary.pres < 2]
print(low_preserved)

### make a point plot
mod_colors <- rownames(preservation_stats)  # Module colors
Z_summary <- preservation_stats$Zsummary.pres
preservation_data <- data.frame(Module = mod_colors, Z_summary = Z_summary)


# Create the point plot
ggplot(preservation_data, aes(x = Module, y = Z_summary, color = Module)) +
  geom_point(size = 8) +  # Larger points
  scale_color_manual(values = mod_colors) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 1) +  # Threshold line
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue", size = 1) +  # Strong preservation line
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),  # Larger x-axis tick text
    axis.text.y = element_text(size = 14),  # Larger y-axis tick text
    axis.title.x = element_text(size = 14, face = "bold"),  # Larger x-axis label
    axis.title.y = element_text(size = 14, face = "bold"),  # Larger y-axis label
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Title formatting
  ) +
  labs(
    title = "Module Preservation Statistics",
    x = "Module Colors",
    y = "Preservation Z-summary"
  )


# Extract genes for each  module
genes_in_modules <- lapply(unique(tumor_module_colors), function(module) {
  names(tumor_module_colors[tumor_module_colors == module])
})
names(genes_in_modules) <- module_names
# View genes for a specific module
head(genes_in_modules$`turquoise`)
# Perform Functional Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
# GO enrichment for the "turquoise" module
blue_genes <- genes_in_modules$`blue`

# Perform GO enrichment
go_results <- enrichGO(
  gene          = blue_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",  # Adjust to your gene ID type
  ont           = "BP",       # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# Examine first 6 rows of GO Enrichment result
head(go_results)

#Create enrich path folder with sub folder GO_T_N to store results
# Use file.path() for OS-independent paths
out_dir <- file.path("enrich", "GO_T_N")

# Create (only if absent)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(go_results, "enrich/GO_T_N/GO_BP_blue.csv", row.names = TRUE)

dotplot(go_results,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Biological Process Enrichment of blue module")

# GO enrichment for the "turquoise" module
blue_genes <- genes_in_modules$`blue`
# Define a function for GO enrichment and CSV writing
perform_go_enrichment <- function(gene_list, ontology, output_path) {
  go_results <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL", 
    ont           = ontology,  # Ontology (BP, CC, MF)
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  # Write results to CSV
  write.csv(go_results, file = output_path, row.names = TRUE)
  return(go_results)
}

# Perform GO enrichment for BP, CC, and MF
go_results_BP <- perform_go_enrichment(blue_genes, "BP", "enrich/GO_T_N/GO_BP_blue.csv")
go_results_CC <- perform_go_enrichment(blue_genes, "CC", "enrich/GO_T_N/GO_CC_blue.csv")
go_results_MF <- perform_go_enrichment(blue_genes, "MF", "enrich/GO_T_N/GO_MF_blue.csv")

#print GO BP for checking
head(go_results_CC)
str(go_results_CC)

# Generate dotplots for each ontology
# Because there is no term found in green module, we skip dot plot of green module for BP 
dotplot_BP <- dotplot(go_results_BP, showCategory=20, font.size=10, label_format=70) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() +
  ggtitle("GO Enrichment - Biological Process (BP) - Blue module")

dotplot_CC <- dotplot(go_results_CC, showCategory=20, font.size=10, label_format=70) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() +
  ggtitle("GO Enrichment - Cellular Component (CC) - Blue module")

dotplot_MF <- dotplot(go_results_MF, showCategory=20, font.size=10, label_format=70) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() +
  ggtitle("GO Enrichment - Molecular Function (MF) - Blue module")

# Combine dotplots into a single image
combined_plot <- ggarrange(
                           dotplot_BP,
                           dotplot_CC, 
                           dotplot_MF, ncol=1, nrow=3)

print(combined_plot)

# Save the combined plot
ggsave("enrich/GO_T_N/GO_combined_dotplot_blue_module.png", combined_plot, width=10, height=15)

library(topGO)
library(GO.db)
library(dplyr)

# Define directories and file patterns
input_dir <- "enrich/GO_T_N/"
file_pattern <- "GO_(BP|CC|MF)_(.*).csv"  # Matches files like go_BP_module1.csv

# Define module preservation levels
module_preservation <- data.frame(
  Module = c("blue", "brown", "green", "grey", "turquoise", "yellow"),  # Replace with your module names
  PreservationLevel = c("High", "Moderate", "Low","Moderate", "High", "High" )  # Corresponding preservation levels
)

# Function to process files
process_go_files <- function(file_path, ontology, module) {
  # Read the GO results file
  go_data <- read.csv(file_path)
  
  # Add ontology and module columns
  go_data <- go_data %>%
    mutate(
      Ontology = ontology,
      Module = module
    )
  
  return(go_data)
}

# List and process all files
go_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)

# Combine all GO results
combined_go_results <- do.call(rbind, lapply(go_files, function(file) {
  # Extract ontology and module from the filename
  match <- regmatches(basename(file), regexec(file_pattern, basename(file)))
  ontology <- match[[1]][2]
  module <- match[[1]][3]
  
  # Process the file
  process_go_files(file, ontology, module)
}))

# Add PreservationLevel to the combined results
combined_go_results <- combined_go_results %>%
  left_join(module_preservation, by = c("Module"))

# View combined data
head(combined_go_results)

# Save combined results for each ontology
for (ontology in c("BP", "CC", "MF")) {
  ontology_data <- combined_go_results %>% filter(Ontology == ontology)
  write.csv(ontology_data, file = paste0("enrich/combined_go_", ontology, ".csv"), row.names = FALSE)
}

go_results_BP <- read.csv("enrich/combined_go_BP.csv")
go_results_CC <- read.csv("enrich/combined_go_CC.csv")
go_results_MF <- read.csv("enrich/combined_go_MF.csv")

# Split by preservation level
high_preservation <- go_results_CC[go_results_MF$PreservationLevel == "High", ]
moderate_preservation <- go_results_CC[go_results_MF$PreservationLevel == "Moderate", ]
low_preservation <- go_results_CC[go_results_MF$PreservationLevel == "Low", ]

high_terms <- unique(high_preservation$Description)
moderate_terms <- unique(moderate_preservation$Description)
low_terms <- unique(low_preservation$Description)


# Overlap analysis
library(VennDiagram)

venn.diagram(
  x = list(
    High = high_terms,
    Moderate = moderate_terms,
    Low = low_terms
  ),
  filename = "enrich/GO_comparison_venn_MF.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5
)

#Find the common genes in high preservation module and low preservation module
intersect_CC <- intersect(high_terms, low_terms)

intersect_MF_mod_high <- intersect(high_terms, moderate_terms)
write.csv(intersect_MF_mod_high,"enrich/intersect_MF_mod_high.csv")
intersect_MF_mod_low <- intersect(low_terms, moderate_terms)
write.csv(intersect_MF_mod_low,"enrich/intersect_MF_mod_low.csv")


###### KEGG/REACTOME pathway####
high_preservation_genes <- c(genes_in_modules$blue, genes_in_modules$yellow, genes_in_modules$turquoise)    # High preservation
high_preservation_genes.entrez_ids <- bitr(high_preservation_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

moderate_preservation_genes <- c(genes_in_modules$brown, genes_in_modules$grey) # Moderate preservation
moderate_preservation_genes.entrez_ids <- bitr(moderate_preservation_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
low_preservation_genes <- c(genes_in_modules$green)    # Low preservation
low_preservation_genes.entrez_ids <- bitr(low_preservation_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment for High preservation genes
kegg_high <- enrichKEGG(
  gene = high_preservation_genes.entrez_ids$ENTREZID,
  organism = 'hsa',  # Human
  pvalueCutoff = 0.05
)

# Perform KEGG enrichment for Moderate preservation genes
kegg_moderate <- enrichKEGG(
  gene = moderate_preservation_genes.entrez_ids$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05
)

# Perform KEGG enrichment for Low preservation genes
kegg_low <- enrichKEGG(
  gene = low_preservation_genes.entrez_ids$ENTREZID,
  organism = 'hsa',
  pvalueCutoff = 0.05
)

# Visualize KEGG enrichment
library(ggplot2)

dotplot(kegg_high, title = "KEGG Pathways - High Preservation") +
  theme_minimal()

dotplot(kegg_moderate, title = "KEGG Pathways - Moderate Preservation") +
  theme_minimal()
dotplot(kegg_low, title = "KEGG Pathways - Low Preservation") +
  theme_minimal()

# Save results to CSV
write.csv(as.data.frame(kegg_high), "kegg_high_preservation.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_moderate), "kegg_moderate_preservation.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_low), "kegg_low_preservation.csv", row.names = FALSE)