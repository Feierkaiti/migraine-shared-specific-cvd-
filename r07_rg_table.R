# Export the genetic correlation matrix (with SEs) from a GenomicSEM ldsc run as a long table.
# Usage: Rscript r07_rg_table.R <root> <ldsc_rds> <out.csv>
args <- commandArgs(trailingOnly = TRUE)
L <- readRDS(args[2]); nm <- colnames(L$S); k <- length(nm)
se <- matrix(0, k, k); se[lower.tri(se, diag = TRUE)] <- sqrt(diag(L$V_Stand)); se <- se + t(se) - diag(diag(se))
out <- expand.grid(trait1 = nm, trait2 = nm, stringsAsFactors = FALSE)
out$rg <- as.vector(L$S_Stand); out$se <- as.vector(se); out$p <- 2 * pnorm(-abs(out$rg / out$se))
out$intercept <- as.vector(L$I)
out <- out[out$trait1 != out$trait2, ]
write.csv(out, args[3], row.names = FALSE)
mig <- out[out$trait1 == "MIG", c("trait2", "rg", "se")]; names(mig) <- c("outcome", "rg", "se")
write.csv(mig, file.path(args[1], "results", "rg_mig_se.csv"), row.names = FALSE)
