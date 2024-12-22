

################ Module Preservation and Reproducibility ##################################
# seurat_ref <- filtered_tumor  # Tumor dataset
# seurat_query <- filtered_normal      # Normal dataset
# 
# Extract module assignments and colors
tumor_modules <- netwk_tumor$colors
tumor_module_colors <- labels2colors(tumor_modules)
names(tumor_module_colors) <- names(netwk_tumor$colors)

normal_module_colors <- labels2colors(netwk_normal$colors)
names(normal_module_colors) <- names(netwk_normal$colors)
# 
# #Map the module from the ref onto query dataset
# # Convert dense matrix to sparse matrix
# seurat_query <- as(seurat_query, "dgCMatrix")
# #Create seurat object
# seurat_query <- CreateSeuratObject(
#   counts = seurat_query,  # The gene expression matrix
#   assay = "RNA"           # Specify the assay name
# )
# 
# load("normal-block.1.RData")
# 
# seurat_query@misc$normal_network <- list(
#   module_colors = normal_module_colors,          # Gene-to-module assignments
#   module_eigengenes = netwk_normal$MEs,          # Module eigengenes
#   TOM = TOM_normal                      # Topological overlap matrix
# )
# 
# 
# seurat_ref <- as(seurat_ref, "dgCMatrix")
# #Create seurat object
# seurat_ref <- CreateSeuratObject(
#   counts = seurat_ref,  # The gene expression matrix
#   assay = "RNA"           # Specify the assay name
# )
# 
# #add WGCNA to seurat
# seurat_ref@misc$tumor_network <- list(
#   module_colors = tumor_module_colors,          # Gene-to-module assignments
#   module_eigengenes = netwk_tumor$MEs,          # Module eigengenes
#   TOM = TOM_tumor                         # Topological overlap matrix
# )
# 
# 
# 
# # Project modules from tumor (reference) to normal (test)
# seurat_query <- ProjectModules(
#   seurat_obj = seurat_query,
#   seurat_ref = seurat_ref,
#   wgcna_name = "tumor_network",  # Reference WGCNA name
#   wgcna_name_proj = "normal_network",  # Query WGCNA name
#   assay = "RNA"  # Specify the assay
# )
# 
# # Extract results for visualization
# preservation_stats <- preservation$preservation$Z[[2]][, -1]

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
data <- data.frame(Module = mod_colors, Z_summary = Z_summary)


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

head(go_results)

write.csv(go_results, "enrich/GO_T_N/GO_BP_blue.csv", row.names = TRUE)

dotplot(go_results,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Biological Process Enrichment of blue module")


# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
# GO enrichment for the "turquoise" module
green_genes <- genes_in_modules$`green`
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
go_results_BP <- perform_go_enrichment(green_genes, "BP", "enrich/GO_T_N/GO_BP_green.csv")
go_results_CC <- perform_go_enrichment(green_genes, "CC", "enrich/GO_T_N/GO_CC_green.csv")
go_results_MF <- perform_go_enrichment(green_genes, "MF", "enrich/GO_T_N/GO_MF_green.csv")

#print GO BP for checking
head(go_results_CC)
str(go_results_CC)

# Generate dotplots for each ontology
# dotplot_BP <- dotplot(go_results_BP, showCategory=20, font.size=10, label_format=70) +
#   scale_size_continuous(range=c(1, 7)) +
#   theme_minimal() +
#   ggtitle("GO Enrichment - Biological Process (BP) - Green module")

dotplot_CC <- dotplot(go_results_CC, showCategory=20, font.size=10, label_format=70) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() +
  ggtitle("GO Enrichment - Cellular Component (CC) - Green module")

dotplot_MF <- dotplot(go_results_MF, showCategory=20, font.size=10, label_format=70) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() +
  ggtitle("GO Enrichment - Molecular Function (MF) - Green module")

# Combine dotplots into a single image
combined_plot <- ggarrange(
  # dotplot_BP,
  dotplot_CC, 
  dotplot_MF, ncol=1, nrow=3)

print(combined_plot)

# Save the combined plot
ggsave("enrich/GO_T_N/GO_combined_dotplot_green_module.png", combined_plot, width=10, height=15)

library(topGO)
library(GO.db)

# Guangchuang Yu. Gene Ontology Semantic Similarity Analysis Using GOSemSim. In: Kidder B. (eds) Stem Cell
# Transcriptional Networks. Methods in Molecular Biology. 2020, 2117:207-215. Humana, New York, NY.
go_results_MF
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