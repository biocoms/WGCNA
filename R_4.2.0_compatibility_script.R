# For R version 4.2.0 on Windows
# Check your R version
R.version.string

# If using R 4.2.0, follow these special instructions:

# 1. First, manually install specific versions of problematic packages
install.packages("https://cran.r-project.org/src/contrib/Archive/xfun/xfun_0.39.tar.gz", 
                 repos = NULL, type = "source")

install.packages("https://cran.r-project.org/src/contrib/Archive/knitr/knitr_1.42.tar.gz", 
                 repos = NULL, type = "source")

install.packages("https://cran.r-project.org/src/contrib/Archive/htmlTable/htmlTable_2.4.0.tar.gz", 
                 repos = NULL, type = "source")

# 2. Install Matrix package compatible with R 4.2.0
install.packages("Matrix", dependencies = TRUE)

# 3. Install BiocManager for R 4.2.0
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.16")

# 4. Install Bioconductor packages compatible with R 4.2.0
BiocManager::install(version = "3.16")
BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"), version = "3.16")

# 5. Install WGCNA with specific method for Windows
install.packages("WGCNA", type = "binary", dependencies = TRUE)

# 6. Install other required packages
required_packages <- c("tidyverse", "dendextend", "gplots", "ggplot2", 
                       "VennDiagram", "dplyr")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, type = "binary")
  }
}

# 7. Install remaining Bioconductor packages
BiocManager::install(c("DESeq2", "genefilter", "clusterProfiler", "org.Hs.eg.db"), 
                     version = "3.16", ask = FALSE)

# 8. For ggpubr (which was problematic in the logs)
install.packages("car", dependencies = TRUE)  # Install car first
install.packages("ggpubr", dependencies = TRUE)
