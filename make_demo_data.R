## make_demo_data.R -----------------------------------------------------------
##
## Simulate a small two-species SNP panel with a KNOWN answer, so the demo file
## shipped with the app doubles as a correctness test: the true individual for
## every sample is written into the file, and any method worth using should
## recover it.
##
## Nothing here is derived from real data.
## Run from genoID/:  Rscript make_demo_data.R
## ---------------------------------------------------------------------------
set.seed(20260806)

N_LOCI      <- 40     # neutral loci for identification
N_DIAG      <- 2      # species-diagnostic loci
N_IND_A     <- 22     # true individuals, species A
N_IND_B     <- 6      # true individuals, species B
DROPOUT     <- 0.030  # per-consensus-genotype allelic dropout
FALSE_ALLE  <- 0.010  # per-consensus-genotype false allele
MISSING     <- 0.070  # per-cell missingness
BASES       <- c("A", "C", "G", "T")

## ---- allele frequencies differ between species, as they would in reality ---
locus_alleles <- replicate(N_LOCI, sample(BASES, 2), simplify = FALSE)
freq_A <- runif(N_LOCI, 0.28, 0.72)   # informative loci, as a real panel would be
freq_B <- pmin(pmax(freq_A + runif(N_LOCI, -0.35, 0.35), 0.05), 0.95)

draw_ind <- function(freq) {
  vapply(seq_len(N_LOCI), function(l) {
    a <- locus_alleles[[l]]
    paste(sort(c(sample(a, 1, prob = c(freq[l], 1 - freq[l])),
                 sample(a, 1, prob = c(freq[l], 1 - freq[l])))), collapse = "")
  }, "")
}

## Species B individuals are drawn as a family group, which is the realistic
## hard case: relatives that a weak panel cannot tell apart from recaptures.
founders <- lapply(1:2, function(i) draw_ind(freq_B))
## one allele transmitted from each parent, Mendel-style
draw_offspring <- function() {
  pick <- function(g) { i <- sample(1:2, 1); substr(g, i, i) }
  vapply(seq_len(N_LOCI), function(l)
    paste(sort(c(pick(founders[[1]][l]), pick(founders[[2]][l]))), collapse = ""), "")
}

true_ind <- c(
  lapply(seq_len(N_IND_A), function(i) list(sp = "alpha", gt = draw_ind(freq_A))),
  lapply(seq_len(2),       function(i) list(sp = "beta",  gt = founders[[i]])),
  lapply(seq_len(N_IND_B - 2), function(i) list(sp = "beta", gt = draw_offspring())))
names(true_ind) <- sprintf("TRUE_%02d", seq_along(true_ind))

## ---- how many times each individual was sampled ---------------------------
n_times <- c(rep(1, 14), 2, 2, 2, 3, 3, 4, 6, 9,      # alpha: one heavily resampled
             1, 1, 2, 2, 3, 1)                        # beta
stopifnot(length(n_times) == length(true_ind))

## ---- observation process: dropout, false alleles, missing data ------------
observe <- function(g, l) {
  a <- strsplit(g, "")[[1]]
  if (runif(1) < MISSING) return("00")
  if (a[1] != a[2] && runif(1) < DROPOUT) a <- rep(sample(a, 1), 2)
  if (runif(1) < FALSE_ALLE) a[sample(1:2, 1)] <- sample(locus_alleles[[l]], 1)
  paste(sort(a), collapse = "")
}

rows <- list(); k <- 0
for (i in seq_along(true_ind)) {
  ind <- true_ind[[i]]
  for (rep_i in seq_len(n_times[i])) {
    k <- k + 1
    obs <- vapply(seq_len(N_LOCI), function(l) observe(ind$gt[l], l), "")
    ## species-diagnostic loci: near-fixed differences, occasionally missing
    diag_gt <- vapply(seq_len(N_DIAG), function(d)
      if (runif(1) < 0.05) "00" else if (ind$sp == "alpha") "CC" else "TT", "")
    rows[[k]] <- c(
      SampleID   = sprintf("DEMO_%03d", k),
      Plate      = sprintf("P%d", 1 + (k %% 3)),
      CollDate   = format(as.Date("2025-03-01") + (i * 7 + rep_i * 3) %% 210),
      setNames(as.list(diag_gt), sprintf("DIAG%02d", seq_len(N_DIAG))),
      setNames(as.list(obs), sprintf("LOC%02d", seq_len(N_LOCI))),
      SexMarker  = if (i %% 3 == 0) "XY" else "XX",
      TrueIndividual = names(true_ind)[i],
      TrueSpecies    = ind$sp)
  }
}

df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
df <- df[sample(nrow(df)), ]                 # shuffle: recaptures are not adjacent
rownames(df) <- NULL

dir.create("app/demo", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "app/demo/demo_genotypes.csv", row.names = FALSE)

cat(sprintf("%d samples from %d true individuals (%d alpha, %d beta)\n",
            nrow(df), length(true_ind), N_IND_A, N_IND_B))
cat(sprintf("largest true cluster: %d samples;  missing cells: %.1f%%\n",
            max(n_times), 100 * mean(as.matrix(df[, grep("^LOC", names(df))]) == "00")))

## ---- the demo must actually be recoverable, or it is a bad demo -----------

source("app/genoID_core.R")
loci <- sprintf("LOC%02d", seq_len(N_LOCI))
gt   <- gid_matrix(df, "SampleID", loci)
res  <- gid_by_group(gt, df$TrueSpecies, gid_method_lr, dropout = DROPOUT,
                     false_allele = FALSE_ALLE, kinship = "full_sib",
                     post_cut = 0.999, min_loci = 12)
cat(sprintf("recovered %d individuals (truth: %d);  ARI vs truth = %.4f\n",
            length(unique(res$assignment$individual)), length(true_ind),
            gid_ari(res$assignment$individual, df$TrueIndividual)))

## also check it survives a mis-specified error rate, since that is what a new
## user will hit before they have measured their own
res2 <- gid_by_group(gt, df$TrueSpecies, gid_method_lr, dropout = 0.005,
                     false_allele = 0.002, kinship = "full_sib",
                     post_cut = 0.999, min_loci = 12)
cat(sprintf("at app defaults (dropout 0.005): %d individuals, ARI = %.4f\n",
            length(unique(res2$assignment$individual)),
            gid_ari(res2$assignment$individual, df$TrueIndividual)))

ex <- gid_by_group(gt, df$TrueSpecies, gid_method_exact, min_loci = 12)
cat(sprintf("exact match for contrast:      %d individuals, ARI = %.4f\n",
            length(unique(ex$assignment$individual)),
            gid_ari(ex$assignment$individual, df$TrueIndividual)))
