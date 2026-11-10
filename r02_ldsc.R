# Multivariable LDSC with GenomicSEM for one trait set.
# Usage: Rscript r02_ldsc.R <root> <set_name> <trait,trait,...>
# Binary traits were munged with effective N, so sample prevalence is 0.5 and the
# population prevalences below only rescale h2 to the liability scale.
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; set_name <- args[2]; traits <- strsplit(args[3], ",")[[1]]
suppressMessages(library(GenomicSEM))
pop <- c(MIG = 0.15, FG_MIG = 0.15, MVP_MIG = 0.15, FG_MIGTRIP = 0.15, FG_MA = 0.05, FG_MO = 0.10,
         MVP_MA = 0.05, CHOQ = 0.15, MDD25 = 0.15, MDD25_CLIN = 0.15, BACK = 0.25, FIB = 0.02,
         FG_ANX = 0.10, FG_IBS = 0.10, CAD = 0.06, IS = 0.03, AS = 0.04, LAS = 0.005, CES = 0.008,
         SVS = 0.006, HF = 0.02, AF = 0.025, UKB_T2D = 0.08, UKB_IHD = 0.08, FG_CHD = 0.06,
         FG_T2D = 0.08, FG_HTN = 0.25, FG_STR = 0.04, EDU = 0.32, SMK = 0.60)
cont <- c("SBP", "DBP", "HDL", "TG", "BMI", "COMP_S", "COMP_M")
sp <- ifelse(traits %in% cont, NA, 0.5)
pp <- ifelse(traits %in% cont, NA, pop[traits])
if (any(is.na(pp) & !traits %in% cont)) stop("missing population prevalence: ", paste(traits[is.na(pp) & !traits %in% cont], collapse = ","))
files <- file.path(root, "data", "munged", paste0(traits, ".sumstats.gz"))
ld <- paste0(file.path(root, "data", "ref", "eur_w_ld_chr"), "/")
setwd(file.path(root, "results"))
out <- ldsc(traits = files, sample.prev = sp, population.prev = pp, ld = ld, wld = ld,
            trait.names = traits, stand = TRUE, ldsc.log = file.path(root, "logs", paste0("ldsc_", set_name)))
saveRDS(out, file.path(root, "results", paste0("ldsc_", set_name, ".rds")))
k <- length(traits)
h2 <- diag(out$S); se <- sqrt(diag(out$V)[cumsum(c(1, k:2))])
write.csv(data.frame(trait = traits, h2_liab = h2, se = se, z = h2 / se, intercept = diag(out$I)),
          file.path(root, "results", paste0("h2_", set_name, ".csv")), row.names = FALSE)
cat("done", set_name, "\n")
