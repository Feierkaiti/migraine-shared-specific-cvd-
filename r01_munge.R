# Munge HapMap3-filtered standardised files (from p04_hm3.py) with GenomicSEM.
# Usage: Rscript r01_munge.R <root> <name> [<name> ...]
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; nms <- args[-1]
suppressMessages(library(GenomicSEM))
hm3 <- file.path(root, "data", "ref", "w_hm3.snplist")
setwd(file.path(root, "data", "munged"))
for (nm in nms) {
  g <- file.path(root, "data", "tmp", paste0(nm, ".hm3.tsv"))
  munge(files = g, hm3 = hm3, trait.names = nm, N = NA, info.filter = 0.9, maf.filter = 0.01,
        log.name = file.path(root, "logs", paste0("munge_", nm)))
}
