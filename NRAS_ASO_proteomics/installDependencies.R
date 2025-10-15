## install dependencies if needed
# from CRAN
install.packages(c("plyr","dplyr","readxl","plotly","stringr","tidyr","knitr",
                   "pheatmap","ggplot2","tidyverse","reticulate","patchwork",
                   "data.table","scales","tibble","reshape2"))

# from Bioconductor
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

# from other sources
if (!require("amlresistancenetworks", quietly = TRUE)) {
  remotes::install_github("https://github.com/PNNL-CompBio/amlresistancenetworks")
}

if (!require(panSEA, quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  
  devtools::install_github("BelindaBGarana/panSEA", "add-tie-handling")
}

if (!require(synapser, quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Install the latest version of synapser (automatically installs compatible dependency versions)
  remotes::install_cran("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))
}

if (!require(PNNL.DMS.utils, quietly = TRUE)) {
  if(!require("remotes", quietly = T)) install.packages("remotes")
  remotes::install_github("PNNL-Comp-Mass-Spec/PNNL.DMS.utils")
}

if (!require(PlexedPiper, quietly = TRUE)) {
  if(!require("remotes", quietly = T)) install.packages("remotes")
  remotes::install_github("PNNL-Comp-Mass-Spec/PlexedPiper")
}

if (!require(MSnSet.utils, quietly = TRUE)) {
  devtools::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
}

if (!require(PCSF, quietly = TRUE)) {
  devtools::install_github("sgosline/PCSF")
}

if (!require(DMEA, quietly = TRUE)) {
  devtools::install_github("BelindaBGarana/DMEA", build=FALSE)
  # build=FALSE option is for Windows users, not sure if it will cause
  # issues for Mac users
}