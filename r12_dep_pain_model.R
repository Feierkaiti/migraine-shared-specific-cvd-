# Three orthogonal components of migraine liability: depression (D), chronic pain beyond depression (P), migraine-specific (M).
# Usage: Rscript r12_dep_pain_model.R <root> <outcomes,...>
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; outs <- strsplit(args[2], ",")[[1]]
suppressMessages(library(GenomicSEM))
L <- readRDS(file.path(root, "results", "ldsc_full.rds"))
sub_cov <- function(L, keep) {
  nm <- colnames(L$S); i <- match(keep, nm)
  k <- length(nm); idx <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  keepv <- which(idx[, 1] %in% i & idx[, 2] %in% i)
  out <- L
  out$S <- L$S[i, i]; out$V <- L$V[keepv, keepv]; out$I <- L$I[i, i]
  if (!is.null(L$S_Stand)) { out$S_Stand <- L$S_Stand[i, i]; out$V_Stand <- L$V_Stand[keepv, keepv] }
  dimnames(out$S) <- list(NULL, keep); colnames(out$S) <- keep
  out
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
base <- paste("D =~ NA*lD*MIG + MDD25 + BACK + FIB",
              "P =~ NA*lP*MIG + BACK + FIB",
              "M =~ NA*lM*MIG",
              "D ~~ 1*D", "P ~~ 1*P", "M ~~ 1*M",
              "D ~~ 0*P", "D ~~ 0*M", "P ~~ 0*M",
              "MIG ~~ 0*MIG", "MDD25 ~~ 0*MDD25",
              "BACK ~~ a*BACK", "FIB ~~ b*FIB", "a > .001", "b > .001",
              "shareD := lD^2/(lD^2 + lP^2 + lM^2)",
              "shareP := lP^2/(lD^2 + lP^2 + lM^2)",
              "shareM := lM^2/(lD^2 + lP^2 + lM^2)", sep = "\n")
rows <- list(); shares <- NULL
for (y in outs) {
  model <- paste0(base, "\n", y, " ~ bD*D + bP*P + bM*M\n", "dDM := bD - bM\n", "dPM := bP - bM\n", "dDP := bD - bP\n")
  fit <- tryCatch(usermodel(sub_cov(L, c("MIG", "MDD25", "BACK", "FIB", y)), estimation = "DWLS", model = model, std.lv = FALSE),
                  error = function(e) { message(y, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  r <- fit$results
  g <- function(l, o, rr) r[r$lhs == l & r$op == o & r$rhs == rr, ][1, ]
  dfn <- function(l) r[r$op == ":=" & r$lhs == l, ][1, ]
  if (is.null(shares)) {
    shares <- data.frame(component = c("D", "P", "M"),
      loading = c(num(g("D", "=~", "MIG")$STD_Genotype), num(g("P", "=~", "MIG")$STD_Genotype), num(g("M", "=~", "MIG")$STD_Genotype)),
      loading_se = c(num(g("D", "=~", "MIG")$STD_Genotype_SE), num(g("P", "=~", "MIG")$STD_Genotype_SE), num(g("M", "=~", "MIG")$STD_Genotype_SE)),
      share = c(num(dfn("shareD")$Unstand_Est), num(dfn("shareP")$Unstand_Est), num(dfn("shareM")$Unstand_Est)),
      share_se = c(num(dfn("shareD")$Unstand_SE), num(dfn("shareP")$Unstand_SE), num(dfn("shareM")$Unstand_SE)))
    lo <- r[r$op == "=~", c("lhs", "rhs", "STD_Genotype", "STD_Genotype_SE")]
    write.csv(lo, file.path(root, "results", "models", "loadings_dep_pain.csv"), row.names = FALSE)
  }
  bD <- g(y, "~", "D"); bP <- g(y, "~", "P"); bM <- g(y, "~", "M")
  sdy <- num(bD$Unstand_Est) / num(bD$STD_Genotype)
  rows[[y]] <- data.frame(outcome = y,
    bD = num(bD$STD_Genotype), bD_se = num(bD$STD_Genotype_SE),
    bP = num(bP$STD_Genotype), bP_se = num(bP$STD_Genotype_SE),
    bM = num(bM$STD_Genotype), bM_se = num(bM$STD_Genotype_SE),
    dPM = num(dfn("dPM")$Unstand_Est) / sdy, dPM_se = num(dfn("dPM")$Unstand_SE) / abs(sdy),
    dDM = num(dfn("dDM")$Unstand_Est) / sdy, dDM_se = num(dfn("dDM")$Unstand_SE) / abs(sdy),
    dDP = num(dfn("dDP")$Unstand_Est) / sdy, dDP_se = num(dfn("dDP")$Unstand_SE) / abs(sdy),
    CFI = num(fit$modelfit$CFI), SRMR = num(fit$modelfit$SRMR))
}
out <- do.call(rbind, rows)
for (x in c("bD", "bP", "bM", "dPM", "dDM", "dDP")) out[[paste0(x, "_p")]] <- 2 * pnorm(-abs(out[[x]] / out[[paste0(x, "_se")]]))
write.csv(out, file.path(root, "results", "models", "outcomes_dep_pain.csv"), row.names = FALSE)
write.csv(shares, file.path(root, "results", "models", "base_dep_pain.csv"), row.names = FALSE)
print(shares, digits = 3); print(out, digits = 3)
