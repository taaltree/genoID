## compare_sample_quality.R ---------------------------------------------------
##
## compare_reps_vs_consensus.R asks whether every PCR replicate beats a
## consensus when every sample is equally good. Scat is not like that: most
## samples are fine and a minority -- old, wet, sun-baked -- drop alleles many
## times as often. This asks what that does to each route, and whether giving
## every sample its own error rate fixes it.
##
## Six routes, each scored against the known truth:
##   consensus first, Sethi and likelihood ratio, at the residual rates that
##     survive the multi-tube rule
##   every replicate, one error rate for the whole dataset
##   every replicate, each sample with its own error scale (genoID's default)
##
## Simulated families, so the competing relatives are real: about 40 animals,
## about 100 samples, three reactions each, 80 SNPs. The error rates are measured
## from the simulated replicates exactly as the app measures them, never given
## the truth.
##
## Run from genoID/:  Rscript compare_sample_quality.R [datasets per scenario]
## Takes about half a minute per dataset.
## ---------------------------------------------------------------------------
source("app/genoID_core.R")

N_SIM    <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[1]) else 8 }
N_LOCI   <- 80
MIN_LOCI <- 30
REL      <- c("unrelated", "full_sib", "parent_offspring")   # the app's default

## Each scenario: the share of poor samples, and (dropout, false allele,
## missing) per reaction for good and poor samples.
SCENARIOS <- list(
  "Uniform quality"         = list(share = 0,    good = c(0.09, 0.013, 0.08), poor = c(0.09, 0.013, 0.08)),
  "30% noisy samples"       = list(share = 0.30, good = c(0.02, 0.005, 0.03), poor = c(0.30, 0.03, 0.20)),
  "15% very noisy samples"  = list(share = 0.15, good = c(0.02, 0.005, 0.03), poor = c(0.45, 0.06, 0.10)),
  "10% with false alleles"  = list(share = 0.10, good = c(0.02, 0.005, 0.03), poor = c(0.04, 0.25, 0.10)))

simulate_one <- function(sc, seed) {
  set.seed(seed)
  p <- runif(N_LOCI, 0.1, 0.5)                          # minor allele frequency
  loci <- sprintf("L%02d", seq_len(N_LOCI))
  founder <- function() rbinom(N_LOCI, 2, p)            # copies of the minor allele
  kid <- function(m, f) rbinom(N_LOCI, 1, m / 2) + rbinom(N_LOCI, 1, f / 2)
  G <- list()
  for (fam in 1:12) {
    m <- founder(); f <- founder()
    if (runif(1) < 0.5) G[[length(G) + 1]] <- m         # a parent is sometimes sampled
    for (k in seq_len(sample(2:4, 1))) G[[length(G) + 1]] <- kid(m, f)
  }
  for (k in 1:6) G[[length(G) + 1]] <- founder()        # unrelated animals
  G <- do.call(rbind, G)
  truth <- rep(sprintf("A%02d", seq_len(nrow(G))), 1 + rpois(nrow(G), 1.6))
  samples <- sprintf("s%03d", seq_along(truth))
  poor <- runif(length(samples)) < sc$share
  rows <- list(); rs <- character(0)
  for (i in seq_along(samples)) {
    e <- if (poor[i]) sc$poor else sc$good
    g <- G[match(truth[i], sprintf("A%02d", seq_len(nrow(G)))), ]
    for (r in 1:3) {
      a1 <- ifelse(runif(N_LOCI) < g / 2, "T", "C")       # one copy per allele slot
      a2 <- ifelse(g == 2, "T", ifelse(g == 0, "C", ifelse(a1 == "T", "C", "T")))
      drop <- a1 != a2 & runif(N_LOCI) < e[1]             # dropout: one allele lost
      who <- runif(N_LOCI) < 0.5
      a1[drop & who] <- a2[drop & who]; a2[drop & !who] <- a1[drop & !who]
      fa <- runif(N_LOCI) < e[2]                          # false allele: a slot redrawn
      new <- ifelse(runif(N_LOCI) < p, "T", "C"); slot <- runif(N_LOCI) < 0.5
      a1[fa & slot] <- new[fa & slot]; a2[fa & !slot] <- new[fa & !slot]
      call <- ifelse(a1 < a2, paste0(a1, "/", a2), paste0(a2, "/", a1))
      call[runif(N_LOCI) < e[3]] <- NA
      rows[[length(rows) + 1]] <- call; rs <- c(rs, samples[i])
    }
  }
  R <- do.call(rbind, rows); colnames(R) <- loci
  list(R = R, sample = rs, truth = setNames(truth, samples))
}

score <- function(asg, truth) {
  a <- setNames(asg$individual, asg$sample); s <- names(a); t <- truth[s]
  same_t <- outer(t, t, "=="); same_a <- outer(a, a, "=="); ut <- upper.tri(same_t)
  c(individuals = length(unique(a)), true_individuals = length(unique(t)),
    ARI = gid_ari(t, a), missed = mean(!same_a[ut][same_t[ut]]),
    joined = sum(same_a[ut] & !same_t[ut]))
}

run_one <- function(nm, seed) {
  z <- simulate_one(SCENARIOS[[nm]], seed)
  C <- gid_consensus_from_reps(z$R, z$sample)
  C <- gid_filter(C, min_locus_call = 0.25, min_sample_call = 0.5)$gt
  keep <- z$sample %in% rownames(C)
  reps <- list(gt = z$R[keep, colnames(C), drop = FALSE], sample = z$sample[keep])
  e  <- gid_estimate_error(C, reps)
  pe <- gid_propagate_error(e$dropout, e$false_allele, n_rep = 3, rule = "taberlet",
                            hom_n = 3, het_n = 2, per_rep_missing = mean(is.na(reps$gt)))
  pe <- pmax(pe, 1e-4)                                  # the app's floor
  fr <- gid_freq_cover(gid_allele_freq(C), reps$gt)
  sethi <- function(d, f, ...) gid_method_sethi(C, freqs = fr, dropout = d, false_allele = f,
                                                relationships = REL, min_loci = MIN_LOCI, ...)
  lr <- function(d, f, ...) gid_method_lr(C, freqs = fr, dropout = d, false_allele = f,
                                          min_loci = MIN_LOCI, ...)
  res <- list(
    "Consensus first: Sethi"     = sethi(pe[["dropout"]], pe[["false_allele"]]),
    "Consensus first: LR"        = lr(pe[["dropout"]], pe[["false_allele"]]),
    "Every replicate, one rate: Sethi" = sethi(e$dropout, e$false_allele, reps = reps, sample_error = FALSE),
    "Every replicate, one rate: LR"    = lr(e$dropout, e$false_allele, reps = reps, sample_error = FALSE),
    "Every replicate, own rates: Sethi" = sethi(e$dropout, e$false_allele, reps = reps),
    "Every replicate, own rates: LR"    = lr(e$dropout, e$false_allele, reps = reps))
  do.call(rbind, lapply(names(res), function(r)
    data.frame(scenario = nm, seed = seed, route = r, t(score(res[[r]]$assignment, z$truth)),
               check.names = FALSE)))
}

all <- do.call(rbind, lapply(names(SCENARIOS), function(nm) {
  cat(nm, "\n")
  do.call(rbind, lapply(seq_len(N_SIM), function(s) run_one(nm, 500 + s)))
}))
agg <- aggregate(cbind(individuals, true_individuals, ARI, missed, joined) ~ scenario + route, all, mean)
agg <- agg[order(match(agg$scenario, names(SCENARIOS)), agg$route), ]
agg$missed <- sprintf("%.1f%%", 100 * agg$missed)
agg$ARI <- round(agg$ARI, 3); agg$individuals <- round(agg$individuals, 1)
agg$true_individuals <- round(agg$true_individuals, 1); agg$joined <- round(agg$joined, 2)
names(agg)[names(agg) == "missed"] <- "recapture pairs missed"
names(agg)[names(agg) == "joined"] <- "pairs wrongly joined"
cat(sprintf("\n%d datasets per scenario\n", N_SIM))
print(agg, row.names = FALSE)
