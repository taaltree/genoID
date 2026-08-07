## compare_reps_vs_consensus.R ------------------------------------------------
##
## Does conditioning on every PCR replicate beat collapsing them to a consensus
## first? On real data you cannot tell, because you do not know the answer. So
## simulate: build individuals, sample them, generate replicate genotypes with
## known dropout and false-allele rates, then run both routes and score each
## against the truth.
##
## Both routes are given the best error rates available to them: the consensus
## route gets the residual rates that survive its own multi-tube rule, the
## replicate route gets the per-reaction rates. Neither is handicapped.
##
## Run from genoID/:  Rscript compare_reps_vs_consensus.R
## ---------------------------------------------------------------------------
source("app/genoID_core.R")

N_IND      <- 30
N_REPS     <- 3          # PCR replicates per sample
## Replicate failures are NOT independent: in the AITRC data 9.2% of loci fail
## in all three reps, 41x more than independence would predict, because it is
## the sample-by-locus combination that is bad rather than the reaction. Model
## that explicitly, or the consensus route is handicapped by missingness it
## would never actually suffer. Calibrated to the real file: 13.1% of replicate
## cells missing, 17.3% of consensus cells missing.
P_HARD     <- 0.09       # locus fails in every replicate
REP_MISS   <- 0.045      # additional independent per-replicate failure
N_SIM      <- 20         # simulations per condition
GRID       <- expand.grid(n_loci = c(20, 30, 40),
                          dropout = c(0.03, 0.10, 0.25),
                          stringsAsFactors = FALSE)
POST_CUT   <- 0.99       # reachable by all these panels; 0.999 is not, and a
                         # threshold the panel cannot meet makes both routes
                         # return every sample as its own individual, which
                         # measures the threshold rather than the method
FALSE_ALLE <- 0.01
BASES      <- c("A", "C", "G", "T")

## how many times each individual was sampled: a realistic skew
TIMES <- c(rep(1, 16), 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 7, 8)
stopifnot(length(TIMES) == N_IND)

simulate_one <- function(n_loci, dropout, seed) {
  set.seed(seed)
  alle <- replicate(n_loci, sample(BASES, 2), simplify = FALSE)
  p    <- runif(n_loci, 0.25, 0.75)

  draw <- function() vapply(seq_len(n_loci), function(l)
    paste(sort(c(sample(alle[[l]], 1, prob = c(p[l], 1 - p[l])),
                 sample(alle[[l]], 1, prob = c(p[l], 1 - p[l])))), collapse = "/"), "")
  truth <- lapply(seq_len(N_IND), function(i) draw())

  ## observation process, one PCR replicate
  obs_one <- function(g, l) {
    a <- strsplit(g, "/")[[1]]
    if (runif(1) < REP_MISS) return(NA_character_)
    if (a[1] != a[2] && runif(1) < dropout) a <- rep(sample(a, 1), 2)
    if (runif(1) < FALSE_ALLE) a[sample(1:2, 1)] <- sample(alle[[l]], 1)
    paste(sort(a), collapse = "/")
  }

  ## Taberlet multi-tube consensus: heterozygote on 2 replicates, homozygote
  ## needs all 3. The same rule the AITRC data were called with.
  consensus_of <- function(v) {
    nz <- v[!is.na(v)]
    if (!length(nz)) return(NA_character_)
    tb   <- table(nz)
    hets <- tb[vapply(names(tb), function(z) {
      s <- strsplit(z, "/")[[1]]; s[1] != s[2] }, TRUE)]
    if (length(hets) && max(hets) >= 2) return(names(hets)[which.max(hets)])
    ac <- table(unlist(lapply(nz, function(z) unique(strsplit(z, "/")[[1]]))))
    if (sum(ac >= 2) == 2) return(paste(sort(names(ac)[ac >= 2]), collapse = "/"))
    homs <- tb[!names(tb) %in% names(hets)]
    if (length(homs) && max(homs) >= 3) return(names(homs)[which.max(homs)])
    NA_character_
  }

  rep_rows <- list(); cons_rows <- list(); who <- character(0); k <- 0
  for (i in seq_len(N_IND)) for (s in seq_len(TIMES[i])) {
    k <- k + 1; id <- sprintf("S%03d", k)
    hard <- runif(n_loci) < P_HARD          # these loci fail in every replicate
    reps <- replicate(N_REPS,
      vapply(seq_len(n_loci), function(l)
        if (hard[l]) NA_character_ else obs_one(truth[[i]][l], l), ""))
    for (r in seq_len(N_REPS)) rep_rows[[length(rep_rows) + 1L]] <-
      list(id = id, gt = reps[, r])
    cons_rows[[k]] <- consensus_of_row <- vapply(seq_len(n_loci),
      function(l) consensus_of(reps[l, ]), "")
    who[k] <- sprintf("IND%02d", i)
  }

  loci <- sprintf("L%02d", seq_len(n_loci))
  cons <- do.call(rbind, cons_rows)
  dimnames(cons) <- list(sprintf("S%03d", seq_len(nrow(cons))), loci)
  rgt <- do.call(rbind, lapply(rep_rows, `[[`, "gt"))
  colnames(rgt) <- loci
  list(consensus = cons, reps = list(gt = rgt, sample = vapply(rep_rows, `[[`, "", "id")),
       truth = who)
}

## residual error left in a consensus call after the multi-tube rule
resid <- function(dropout) {
  r <- gid_propagate_error(dropout, FALSE_ALLE, n_rep = N_REPS, rule = "taberlet",
                           hom_n = 3, het_n = 2,
                           per_rep_missing = P_HARD + (1 - P_HARD) * REP_MISS,
                           nsim = 8000, seed = 7)
  c(d = max(r[["dropout"]], 1e-5), f = max(r[["false_allele"]], 1e-5))
}

out <- list(); row <- 0
for (g in seq_len(nrow(GRID))) {
  n_loci <- GRID$n_loci[g]; dd <- GRID$dropout[g]
  rr <- resid(dd)
  cat(sprintf("loci %2d  dropout %.2f  ->  consensus residual d=%.4f f=%.4f\n",
              n_loci, dd, rr[["d"]], rr[["f"]]))
  for (sim in seq_len(N_SIM)) {
    S <- simulate_one(n_loci, dd, seed = 1000 * g + sim)
    ml <- max(6, floor(n_loci * 0.4))

    a <- gid_method_lr(S$consensus, dropout = rr[["d"]], false_allele = rr[["f"]],
                       kinship = "full_sib", post_cut = POST_CUT, min_loci = ml)
    b <- gid_method_lr(S$consensus, dropout = dd, false_allele = FALSE_ALLE,
                       kinship = "full_sib", post_cut = POST_CUT, min_loci = ml,
                       reps = S$reps)

    ## Threshold-free measure of how well each route separates the two kinds of
    ## pair: the gap between the weakest true recapture and the strongest false
    ## one. Positive means a threshold exists that gets every pair right.
    sep <- function(r) {
      same <- S$truth[match(r$pairs$id1, rownames(S$consensus))] ==
              S$truth[match(r$pairs$id2, rownames(S$consensus))]
      ok <- r$pairs$n_compared >= ml
      if (!any(same & ok) || !any(!same & ok)) return(c(NA, NA, NA))
      c(mean(r$pairs$log10_LR[same & ok]),
        quantile(r$pairs$log10_LR[!same & ok], 0.999, names = FALSE),
        min(r$pairs$log10_LR[same & ok]) -
          quantile(r$pairs$log10_LR[!same & ok], 0.999, names = FALSE))
    }
    sa <- sep(a); sb <- sep(b)

    row <- row + 1
    out[[row]] <- data.frame(
      n_loci = n_loci, dropout = dd, sim = sim,
      truth_n = length(unique(S$truth)),
      cons_n = length(unique(a$assignment$individual)),
      reps_n = length(unique(b$assignment$individual)),
      cons_ari = gid_ari(a$assignment$individual, S$truth),
      reps_ari = gid_ari(b$assignment$individual, S$truth),
      cons_loci = median(a$pairs$n_compared),
      reps_loci = median(b$pairs$n_compared),
      cons_missing = mean(is.na(S$consensus)),
      cons_same_lr = sa[1], reps_same_lr = sb[1],
      cons_margin = sa[3], reps_margin = sb[3],
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, out)
dir.create("outputs", showWarnings = FALSE)
write.csv(res, "outputs/reps_vs_consensus_simulation.csv", row.names = FALSE)

agg <- aggregate(cbind(cons_ari, reps_ari, cons_n, reps_n, cons_loci, reps_loci,
                       cons_missing, cons_same_lr, reps_same_lr,
                       cons_margin, reps_margin) ~ n_loci + dropout,
                 data = res, FUN = function(z) mean(z, na.rm = TRUE))
agg$ari_gain <- agg$reps_ari - agg$cons_ari
agg$err_cons <- abs(agg$cons_n - N_IND)
agg$err_reps <- abs(agg$reps_n - N_IND)

cat("\n================ mean over", N_SIM, "simulations, truth =", N_IND, "individuals ================\n")
print(within(agg, {
  cons_ari <- round(cons_ari, 3); reps_ari <- round(reps_ari, 3)
  ari_gain <- round(ari_gain, 3); cons_n <- round(cons_n, 1); reps_n <- round(reps_n, 1)
  cons_loci <- round(cons_loci, 1); reps_loci <- round(reps_loci, 1)
  cons_missing <- round(cons_missing, 3)
})[, c("n_loci", "dropout", "cons_loci", "reps_loci",
       "cons_n", "reps_n", "cons_ari", "reps_ari", "ari_gain")], row.names = FALSE)

cat("\nEvidence strength, independent of any threshold:\n")
cat("  same_lr = mean log10 LR for pairs that really ARE the same animal\n")
cat("  margin  = weakest true recapture minus strongest false one; >0 means a\n")
cat("            threshold exists that classifies every pair correctly\n\n")
print(within(agg, {
  cons_same_lr <- round(cons_same_lr, 2); reps_same_lr <- round(reps_same_lr, 2)
  cons_margin <- round(cons_margin, 2); reps_margin <- round(reps_margin, 2)
})[, c("n_loci", "dropout", "cons_same_lr", "reps_same_lr",
       "cons_margin", "reps_margin")], row.names = FALSE)

cat("\nHow often did each route land closer to the truth?\n")
cat(sprintf("  replicates closer : %d of %d simulations\n",
            sum(abs(res$reps_n - N_IND) < abs(res$cons_n - N_IND)), nrow(res)))
cat(sprintf("  consensus closer  : %d\n",
            sum(abs(res$reps_n - N_IND) > abs(res$cons_n - N_IND))))
cat(sprintf("  tied              : %d\n",
            sum(abs(res$reps_n - N_IND) == abs(res$cons_n - N_IND))))
cat(sprintf("\nMean ARI: consensus %.4f, replicates %.4f (gain %+.4f)\n",
            mean(res$cons_ari), mean(res$reps_ari),
            mean(res$reps_ari - res$cons_ari)))
tt <- t.test(res$reps_ari, res$cons_ari, paired = TRUE)
cat(sprintf("Paired t-test on ARI: t = %.2f, df = %d, p = %.3g\n",
            tt$statistic, tt$parameter, tt$p.value))
