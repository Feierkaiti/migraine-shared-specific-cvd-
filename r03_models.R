# Shared vs migraine-specific decomposition and outcome models (GenomicSEM usermodel, DWLS).
# Usage: Rscript r03_models.R <root> <ldsc_rds> <tag> <mig> <shared,...> <outcomes,...|NONE> [<covariates,...>]
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; rds <- args[2]; tag <- args[3]; mig <- args[4]
shared <- strsplit(args[5], ",")[[1]]
outs <- if (args[6] == "NONE") character(0) else strsplit(args[6], ",")[[1]]
covs <- if (length(args) >= 7) strsplit(args[7], ",")[[1]] else character(0)
suppressMessages(library(GenomicSEM))
L <- readRDS(rds)

sub_cov <- function(L, keep) {
  nm <- colnames(L$S); i <- match(keep, nm)
  if (any(is.na(i))) stop("not in LDSC output: ", paste(keep[is.na(i)], collapse = ","))
  k <- length(nm); idx <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  keepv <- which(idx[, 1] %in% i & idx[, 2] %in% i)
  out <- L
  out$S <- L$S[i, i]; out$V <- L$V[keepv, keepv]; out$I <- L$I[i, i]
  if (!is.null(L$S_Stand)) { out$S_Stand <- L$S_Stand[i, i]; out$V_Stand <- L$V_Stand[keepv, keepv] }
  dimnames(out$S) <- list(NULL, keep); colnames(out$S) <- keep
  out
}
base <- paste0("S =~ NA*lS*", mig, " + ", paste(shared, collapse = " + "), "\n",
               "M =~ NA*lM*", mig, "\n",
               "S ~~ 1*S\nM ~~ 1*M\nS ~~ 0*M\n", mig, " ~~ 0*", mig, "\n",
               # one shared indicator: S is that trait's genetic factor (Cholesky), so its residual is 0
               if (length(shared) == 1) paste0(shared, " ~~ 0*", shared, "\n") else "",
               "share := lS^2/(lS^2 + lM^2)\n")
get <- function(res, lhs, op, rhs) res[res$lhs == lhs & res$op == op & res$rhs == rhs, ][1, ]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

rows <- list()
fit0 <- usermodel(sub_cov(L, c(mig, shared)), estimation = "DWLS", model = base, std.lv = FALSE)
r0 <- fit0$results
mf <- fit0$modelfit
lab <- function(r, lhs, op, rhs) {
  x <- get(r, lhs, op, rhs)
  c(est = num(x$STD_Genotype), se = num(x$STD_Genotype_SE), unstd = num(x$Unstand_Est), unstd_se = num(x$Unstand_SE))
}
loadS <- lab(r0, "S", "=~", mig); loadM <- lab(r0, "M", "=~", mig)
sh <- get(r0, "share", ":=", "lS^2/(lS^2+lM^2)")
if (nrow(sh) == 0 || is.na(sh$lhs)) sh <- r0[r0$op == ":=", ][1, ]
base_row <- data.frame(tag = tag, outcome = "BASE", mig = mig, shared = paste(shared, collapse = "+"),
  loadS = loadS["est"], loadS_se = loadS["se"], loadM = loadM["est"], loadM_se = loadM["se"],
  share = num(sh$Unstand_Est), share_se = num(sh$Unstand_SE),
  chisq = num(mf$chisq), df = num(mf$df), CFI = num(mf$CFI), SRMR = num(mf$SRMR))
others <- r0[r0$op == "=~" & r0$lhs == "S" & r0$rhs != mig, c("rhs", "STD_Genotype", "STD_Genotype_SE")]
write.csv(others, file.path(root, "results", paste0("loadings_", tag, ".csv")), row.names = FALSE)

for (y in outs) {
  nmS <- colnames(L$S); rg_y <- L$S_Stand[match(mig, nmS), match(y, nmS)]
  model <- paste0(base, y, " ~ bS*S + bM*M\n", "d := bS - bM\n")
  if (length(covs)) {
    model <- paste0(base, y, " ~ bS*S + bM*M + ", paste(covs, collapse = " + "), "\n",
                    paste0("S ~~ ", paste(covs, collapse = " + "), "\n", "M ~~ ", paste(covs, collapse = " + "), "\n"),
                    if (length(covs) > 1) paste0(covs[1], " ~~ ", paste(covs[-1], collapse = " + "), "\n") else "",
                    "d := bS - bM\n")
  }
  fit <- tryCatch(usermodel(sub_cov(L, c(mig, shared, y, covs)), estimation = "DWLS", model = model, std.lv = FALSE),
                  error = function(e) { message(y, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  r <- fit$results
  bS <- lab(r, y, "~", "S"); bM <- lab(r, y, "~", "M")
  dd <- r[r$op == ":=" & r$lhs == "d", ][1, ]
  # difference on the standardised scale: both paths share the same scaling (factor variances 1)
  sdy <- num(bS["unstd"]) / num(bS["est"])
  rows[[y]] <- data.frame(tag = tag, outcome = y, mig = mig, shared = paste(shared, collapse = "+"),
    covariates = paste(covs, collapse = "+"),
    rg_mig = rg_y, bS = bS["est"], bS_se = bS["se"], bM = bM["est"], bM_se = bM["se"],
    diff = num(dd$Unstand_Est) / sdy, diff_se = num(dd$Unstand_SE) / abs(sdy),
    CFI = num(fit$modelfit$CFI), SRMR = num(fit$modelfit$SRMR))
}
out <- if (length(rows)) do.call(rbind, rows) else NULL
dir.create(file.path(root, "results", "models"), showWarnings = FALSE)
write.csv(base_row, file.path(root, "results", "models", paste0("base_", tag, ".csv")), row.names = FALSE)
if (!is.null(out)) {
  out$bS_p <- 2 * pnorm(-abs(out$bS / out$bS_se)); out$bM_p <- 2 * pnorm(-abs(out$bM / out$bM_se))
  out$diff_p <- 2 * pnorm(-abs(out$diff / out$diff_se))
  write.csv(out, file.path(root, "results", "models", paste0("outcomes_", tag, ".csv")), row.names = FALSE)
  print(out[, c("outcome", "rg_mig", "bS", "bS_se", "bM", "bM_se", "diff", "diff_se", "diff_p")], digits = 3)
}
print(base_row, digits = 3)
