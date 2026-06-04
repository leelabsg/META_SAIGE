# Synthesize a binary-trait SAIGE step2 GWAS file from the continuous LDL file.
# We are NOT trying to reproduce the LDL p-values — we just need a structurally valid
# binary input so we can confirm the binary code path runs end-to-end without errors.
suppressMessages({
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
in_path  <- args[1]
out_path <- args[2]

gwas <- fread(in_path)
n_total <- gwas$N[1]
n_case  <- round(n_total / 2)
n_ctrl  <- n_total - n_case

bin <- data.table(
  CHR        = gwas$CHR,
  POS        = gwas$POS,
  MarkerID   = gwas$MarkerID,
  Allele1    = gwas$Allele1,
  Allele2    = gwas$Allele2,
  AC_Allele2 = gwas$AC_Allele2,
  AF_Allele2 = gwas$AF_Allele2,
  MissingRate= gwas$MissingRate,
  BETA       = gwas$BETA,
  SE         = gwas$SE,
  Tstat      = gwas$Tstat,
  var        = gwas$var,
  p.value    = gwas$p.value,
  p.value.NA = gwas$p.value,
  Is.SPA     = FALSE,
  N_case     = n_case,
  N_ctrl     = n_ctrl,
  N_case_hom = 0L,
  N_ctrl_hom = 0L,
  N_case_het = pmax(0L, as.integer(round(gwas$AC_Allele2 * n_case / n_total))),
  N_ctrl_het = pmax(0L, as.integer(round(gwas$AC_Allele2 * n_ctrl / n_total))),
  N          = gwas$N
)

fwrite(bin, out_path, sep = "\t")
cat("Wrote binary fixture:", out_path, "rows:", nrow(bin), "\n")
