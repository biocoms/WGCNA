# WGCNA and Module Preservation
This repository demonstrates a complete **Weighted Gene-Co-expression Network Analysis (WGCNA)** pipeline, including module preservation across datasets.

---

## 1  Dataset

1. **Download** `GeneExpression.zip`.
2. **Unzip** it in the project root – you should end up with:
```text
GeneExpression/
├── OSCC_TCGA_gene_expression_337.csv
├── OSCC_TCGA_gene_expression_32.csv
```

## 2  Code

### `wgcna.R`

1. Pre-processes the expression matrix  
2. Builds the WGCNA network  
3. Calculates module eigengenes  
4. Performs module-preservation analysis  

## 3  How to run
1. Place `wgcna.R` in the same directory **alongside** the `GeneExpression/` folder.  
2. Open the script in **RStudio** to run the code

## Citation

If you use this repository or the pipeline in your research, please cite:

Nguyen, P. and Zeng, E. (2025). A Protocol for Weighted Gene Co-expression Network Analysis With Module Preservation and Functional Enrichment Analysis for Tumor and Normal Transcriptomic Data. Bio-protocol, 15(18): e5447. https://doi.org/10.21769/BioProtoc.5447
