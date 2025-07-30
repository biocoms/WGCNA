# For R version 4.2.0 on Windows
# Check your R version
R.version.string

# If using R 4.2.0, follow these special instructions:

# 1. First, manually install specific versions of problematic packages
install.packages("https://cran.r-project.org/src/contrib/Archive/xfun/xfun_0.39.tar.gz", 
                 repos = NULL, type = "source")

install.packages("https://cran.r-project.org/src/contrib/Archive/knitr/knitr_1.42.tar.gz", 
                 repos = NULL, type = "source")

#Caution: You might encounter error when installing htmlTable due to dependency. The common dependency is "checkmate"
install.packages("https://cran.r-project.org/src/contrib/Archive/checkmate/checkmate_2.1.0.tar.gz", 
                 repos = NULL, type = "source")

install.packages("https://cran.r-project.org/src/contrib/Archive/htmlTable/htmlTable_2.4.0.tar.gz", 
                 repos = NULL, type = "source")

# 2. Install Matrix package compatible with R 4.2.0
install.packages("Matrix", dependencies = TRUE)
#if above does not works, manually install older version of Matrix package
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.tar.gz", 
                 repos = NULL, type = "source")

#Or
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.4-1.tar.gz", 
                 repos = NULL, type = "source")

# 3. Install BiocManager for R 4.2.0
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.16")

# 4. Install Bioconductor packages compatible with R 4.2.0
BiocManager::install(version = "3.16")
BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"), version = "3.16")

# If you get warning message:
#Warning message:
#package(s) not installed when version(s) same as or greater than current; use `force = TRUE` to
#  re-install: 'AnnotationDbi' 

BiocManager::install("AnnotationDbi", force = TRUE)


# 5. Install WGCNA with specific method for Windows
install.packages("WGCNA", type = "binary", dependencies = TRUE)
#If using Linux, remove type = "binary"
install.packages("WGCNA", dependencies = TRUE)

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

# 8. For ggpubr
install.packages("car", dependencies = TRUE)  # Install car first
install.packages("ggpubr", dependencies = TRUE)
