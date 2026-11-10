options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
pk <- c("remotes","data.table","lavaan","gdata","stringr","doParallel","foreach","R.utils","plyr","readr","Rcpp","jsonlite")
inst <- pk[!pk %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst)
remotes::install_github("GenomicSEM/GenomicSEM", upgrade="never")
library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "\n")
