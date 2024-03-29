
install.packages("ggplot2",repos = "http://cran.us.r-project.org")
install.packages("tidyr",repos = "http://cran.us.r-project.org")
install.packages("data.table",repos = "http://cran.us.r-project.org")
install.packages("feather",repos="http://cran.us.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",version="3.10")

BiocManager::install(version="3.1.0")
BiocManager::install("DropletUtils")
BiocManager::install("scds")
BiocManager::install("scran")
BiocManager::install("limma")

