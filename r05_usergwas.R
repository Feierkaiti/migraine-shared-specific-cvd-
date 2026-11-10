# SNP-level decomposition GWAS (GenomicSEM sumstats + userGWAS): effects of each SNP on S and M.
# Usage: Rscript r05_usergwas.R <root> <tag> <ldsc_rds> <mig> <shared,...> <cores> [<n_test_snps>|0] [Q_SNP] [stage1|stage2]
# stage1 = HapMap3 variants plus any variant with P < 1e-4 in an input GWAS (all candidate loci; LDSC set);
# stage2 = all remaining variants. stage1 + stage2 = the full genome-wide run.
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; tag <- args[2]; rds <- args[3]; mig <- args[4]
shared <- strsplit(args[5], ",")[[1]]; cores <- as.integer(args[6])
ntest <- if (length(args) >= 7) as.integer(args[7]) else NA; if (!is.na(ntest) && ntest == 0) ntest <- NA
QSNP <- if (length(args) >= 8) as.logical(args[8]) else FALSE
stage <- if (length(args) >= 9) args[9] else NA
suppressMessages({library(GenomicSEM); library(data.table)})
traits <- c(mig, shared)
indir <- file.path(root, "data", "gwas_in", tag)
ref <- file.path(root, "data", "ref", "ref_1000G_EUR_maf005.txt")
outdir <- file.path(root, "results", "gwas", tag); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ss_file <- file.path(outdir, "sumstats.rds")
if (!file.exists(ss_file)) {
  ss <- sumstats(files = file.path(indir, paste0(traits, ".tsv")), ref = ref, trait.names = traits,
                 se.logit = rep(TRUE, length(traits)), OLS = NULL, linprob = NULL,
                 info.filter = 0.6, maf.filter = 0.01, keep.indel = FALSE, parallel = TRUE, cores = cores)
  saveRDS(ss, ss_file)
} else ss <- readRDS(ss_file)
cat("sumstats SNPs:", nrow(ss), "\n")

L <- readRDS(rds); nm <- colnames(L$S); i <- match(traits, nm); k <- length(nm)
idx <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE); kv <- which(idx[, 1] %in% i & idx[, 2] %in% i)
L2 <- L; L2$S <- L$S[i, i]; L2$V <- L$V[kv, kv]; L2$I <- L$I[i, i]; L2$S_Stand <- L$S_Stand[i, i]; L2$V_Stand <- L$V_Stand[kv, kv]
colnames(L2$S) <- traits

model <- paste0("S =~ NA*", mig, " + ", paste(shared, collapse = " + "), "\n",
                "M =~ NA*", mig, "\n", "S ~~ 1*S\nM ~~ 1*M\nS ~~ 0*M\n", mig, " ~~ 0*", mig, "\n",
                if (length(shared) == 1) paste0(shared, " ~~ 0*", shared, "\n") else "",
                "S ~ SNP\nM ~ SNP\n")
if (!is.na(ntest)) { set.seed(20261105); ss <- ss[sort(sample(nrow(ss), ntest)), ] }
if (!is.na(stage)) {
  hm3 <- fread(file.path(root, "data", "ref", "w_hm3.snplist"))$SNP
  minp <- do.call(pmin, lapply(traits, function(t) 2 * pnorm(-abs(ss[[paste0("beta.", t)]] / ss[[paste0("se.", t)]]))))
  in1 <- ss$SNP %in% hm3 | minp < 1e-4
  cat("stage1 SNPs:", sum(in1), " stage2 SNPs:", sum(!in1), "\n")
  ss <- if (stage == "stage1") ss[in1, ] else ss[!in1, ]
}
t0 <- Sys.time()
res <- userGWAS(covstruc = L2, SNPs = ss, estimation = "DWLS", model = model, sub = c("S~SNP", "M~SNP"),
                cores = cores, parallel = TRUE, fix_measurement = TRUE, Q_SNP = QSNP, printwarn = TRUE)  # printwarn = FALSE breaks column naming in GenomicSEM 0.0.5
cat("userGWAS time:", format(Sys.time() - t0), "for", nrow(ss), "SNPs\n")
sfx <- if (!is.na(stage)) paste0("_", stage) else if (is.na(ntest)) "" else paste0("_test", ntest)
fwrite(res[[1]], file.path(outdir, paste0("S_SNP", sfx, ".tsv.gz")), sep = "\t")
fwrite(res[[2]], file.path(outdir, paste0("M_SNP", sfx, ".tsv.gz")), sep = "\t")
cat("done\n")
