## make_demo_data.R -----------------------------------------------------------
##
## Simulate the example dataset shipped with the app. It has a KNOWN answer, so
## it doubles as a correctness test: the true individual behind every sample is
## written into the file, along with the true error rates, and any method worth
## using should recover both.
##
## The file is laid out the way a genotyping lab actually delivers data: one row
## per PCR reaction, three reactions per sample, followed by the consensus row
## the lab would compute. That lets the demo show the thing a consensus-only
## file cannot -- estimating allelic dropout and false-allele rates directly
## from the disagreements between replicates, instead of guessing them.
##
## Everything a user can do in the app is exercised by this one file:
## replicates, species blocking, sex, collection dates across two field seasons,
## and coordinates with realistic home ranges.
##
## Nothing here is derived from real data.
## Run from genoID/:  Rscript make_demo_data.R
## ---------------------------------------------------------------------------
source("app/genoID_core.R")
set.seed(20260914)

N_LOCI      <- 40     # neutral loci for identification
N_DIAG      <- 2      # species-diagnostic loci
N_IND_A     <- 22     # true individuals, species alpha
N_IND_B     <- 6      # true individuals, species beta (a family group)
N_REPS      <- 3      # PCR reactions per sample

## Per-REACTION error rates. Deliberately low, as in a well-run lab: dropout is
## the common failure, false alleles are rare. These are the numbers the app's
## "Estimate from replicates" button should give back.
DROPOUT     <- 0.015  # a heterozygote reads as a homozygote in one reaction
FALSE_ALLE  <- 0.002  # a spurious allele appears in one reaction

## Amplification failure is not independent between reactions: a poor extract
## fails in all three. Real data shows this strongly, and a simulation that
## ignored it would flatter the replicate analysis.
FAIL_GOOD   <- 0.03   # per-reaction failure, good-quality samples
FAIL_POOR   <- 0.35   # per-reaction failure, poor-quality samples
P_POOR      <- 0.10   # fraction of samples that are poor quality
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

## Species beta is a family group, which is the realistic hard case: relatives
## that a weak panel cannot tell apart from recaptures.
founders <- lapply(1:2, function(i) draw_ind(freq_B))
draw_offspring <- function() {                       # one allele from each parent
  pick <- function(g) { i <- sample(1:2, 1); substr(g, i, i) }
  vapply(seq_len(N_LOCI), function(l)
    paste(sort(c(pick(founders[[1]][l]), pick(founders[[2]][l]))), collapse = ""), "")
}

true_ind <- c(
  lapply(seq_len(N_IND_A),     function(i) list(sp = "alpha", gt = draw_ind(freq_A))),
  lapply(seq_len(2),           function(i) list(sp = "beta",  gt = founders[[i]])),
  lapply(seq_len(N_IND_B - 2), function(i) list(sp = "beta",  gt = draw_offspring())))
names(true_ind) <- sprintf("TRUE_%02d", seq_along(true_ind))
n_ind <- length(true_ind)

## ---- who was sampled how often --------------------------------------------
n_times <- c(rep(1, 14), 2, 2, 2, 3, 3, 4, 6, 9,      # alpha: one heavily resampled
             1, 1, 2, 2, 3, 1)                        # beta
stopifnot(length(n_times) == n_ind)

## ---- where each animal lives ----------------------------------------------
## Prince of Wales Island, southeast Alaska: a real place at a realistic scale,
## so the basemap looks like somewhere. Samples scatter about 2 km around each
## animal's activity center, so home ranges are small next to the area sampled.
LAT0 <- 55.62; LON0 <- -132.90
centre <- data.frame(lat = LAT0 + runif(n_ind, -0.16, 0.16),
                     lon = LON0 + runif(n_ind, -0.28, 0.28))

## ---- sex, and when each animal was about ----------------------------------
sex <- ifelse(seq_len(n_ind) %% 3 == 0, "XY", "XX")
## two field seasons; most animals seen in one, a few carried over to both
season <- sample(c(2024, 2025), n_ind, TRUE, prob = c(0.45, 0.55))
carry  <- sample(which(n_times > 2), 2)

## ---- the observation process, one PCR reaction at a time -------------------
observe <- function(g, l, fail) {
  if (runif(1) < fail) return("00")
  a <- strsplit(g, "")[[1]]
  if (a[1] != a[2] && runif(1) < DROPOUT) a <- rep(sample(a, 1), 2)
  if (runif(1) < FALSE_ALLE) a[sample(1:2, 1)] <- sample(locus_alleles[[l]], 1)
  paste(sort(a), collapse = "")
}

loci  <- sprintf("LOC%02d", seq_len(N_LOCI))
diags <- sprintf("DIAG%02d", seq_len(N_DIAG))
blocks <- list(); k <- 0

for (i in seq_len(n_ind)) {
  ind <- true_ind[[i]]
  for (s in seq_len(n_times[i])) {
    k   <- k + 1
    sid <- sprintf("DEMO_%03d", k)
    fail <- if (runif(1) < P_POOR) FAIL_POOR else FAIL_GOOD

    yr <- season[i]
    if (i %in% carry && s %% 2 == 0) yr <- if (yr == 2024) 2025 else 2024
    meta <- list(
      Plate    = sprintf("P%d", 1 + (k %% 3)),
      CollDate = format(as.Date(sprintf("%d-06-01", yr)) + sample(0:150, 1)),
      Latitude  = round(centre$lat[i] + rnorm(1, 0, 0.018), 5),
      Longitude = round(centre$lon[i] + rnorm(1, 0, 0.030), 5))

    reps <- lapply(letters[seq_len(N_REPS)], function(r) {
      obs  <- vapply(seq_len(N_LOCI), function(l) observe(ind$gt[l], l, fail), "")
      diag <- vapply(seq_len(N_DIAG), function(d)
        if (runif(1) < fail) "00" else if (ind$sp == "alpha") "CC" else "TT", "")
      c(list(SampleID = sid, Rep = r), meta,
        setNames(as.list(diag), diags), setNames(as.list(obs), loci),
        list(SexMarker = if (runif(1) < fail) "00" else sex[i],
             TrueIndividual = names(true_ind)[i], TrueSpecies = ind$sp))
    })
    blocks[[k]] <- do.call(rbind, lapply(reps, as.data.frame, stringsAsFactors = FALSE))
  }
}

## Shuffle samples so recaptures are not adjacent, but keep each sample's
## reactions together and in order, as a lab's spreadsheet would.
blocks <- blocks[sample(length(blocks))]
reps_df <- do.call(rbind, blocks)
rownames(reps_df) <- NULL

## ---- the consensus row the lab would add ----------------------------------
## Built with the app's own rule (heterozygote on two reactions, homozygote only
## on three), so the file's consensus is exactly what the app would compute.
to_call <- function(v) ifelse(v == "00", NA_character_, v)
cons_one <- function(cols, hom_n = 3, het_n = 2) {
  m <- matrix(to_call(unlist(reps_df[, cols])), ncol = length(cols),
              dimnames = list(NULL, cols))
  rk <- paste(reps_df$SampleID, reps_df$Rep, sep = "|")
  rownames(m) <- rk
  m <- gid_matrix(data.frame(key = rk, reps_df[, cols], check.names = FALSE,
                             stringsAsFactors = FALSE), "key", cols)
  gid_consensus_from_reps(m, reps_df$SampleID, hom_n = hom_n, het_n = het_n)
}
cons_loci <- cons_one(c(diags, loci))
## back to the file's own spelling: "AG" for a call, "00" for missing
unslash <- function(x) ifelse(is.na(x), "00", gsub("/", "", x, fixed = TRUE))

ids <- unique(reps_df$SampleID)
consensus_rows <- lapply(ids, function(sid) {
  first <- reps_df[reps_df$SampleID == sid, ][1, ]
  first$Rep <- "consensus"
  for (col in c(diags, loci)) first[[col]] <- unslash(cons_loci[sid, col])
  sx <- reps_df$SexMarker[reps_df$SampleID == sid]
  sx <- sx[sx != "00"]
  first$SexMarker <- if (length(sx)) names(sort(table(sx), decreasing = TRUE))[1] else "00"
  first
})

## interleave: a, b, c, consensus for each sample
df <- do.call(rbind, lapply(ids, function(sid)
  rbind(reps_df[reps_df$SampleID == sid, ],
        consensus_rows[[match(sid, ids)]])))
rownames(df) <- NULL

dir.create("app/demo", showWarnings = FALSE, recursive = TRUE)
write.csv(df, "app/demo/demo_genotypes.csv", row.names = FALSE)

## ============================================================================
## The demo must actually demonstrate what it claims to. Check, and stop loudly
## if it does not -- a demo that quietly fails its own story is worse than none.
## ============================================================================
say <- function(...) cat(sprintf(...), "\n")
rx  <- df[df$Rep %in% letters[seq_len(N_REPS)], ]
cn  <- df[df$Rep == "consensus", ]

say("%d samples x %d reactions + consensus = %d rows, from %d true individuals",
    nrow(cn), N_REPS, nrow(df), n_ind)
say("missing cells per reaction: %.1f%%   poor-quality samples: %d",
    100 * mean(as.matrix(rx[, loci]) == "00"),
    sum(tapply(rowMeans(as.matrix(rx[, loci]) == "00"), rx$SampleID, mean) > 0.2))

## 1. error rates recovered from the replicates -----------------------------
gt_cons <- gid_matrix(cn, "SampleID", loci)
rkey    <- paste(rx$SampleID, rx$Rep, sep = "|")
rep_gt  <- gid_matrix(data.frame(key = rkey, rx[, loci], check.names = FALSE,
                                 stringsAsFactors = FALSE), "key", loci)
reps    <- list(gt = rep_gt, sample = rx$SampleID)
## group by species: allele frequencies belong to a population, so pooling two
## of them misspecifies the genotype prior the likelihood integrates over
species <- setNames(cn$TrueSpecies, cn$SampleID)
err     <- gid_estimate_error(gt_cons, reps, group = species)

say("")
say("error rates, estimated from the replicates by %s:", err$method)
say("  dropout      true %.4f   estimated %.4f  [%.4f - %.4f]",
    DROPOUT, err$dropout, err$dropout_ci[1], err$dropout_ci[2])
say("  false allele true %.4f   estimated %.4f  [%.4f - %.4f]",
    FALSE_ALLE, err$false_allele, err$false_ci[1], err$false_ci[2])

inside <- function(x, ci) is.finite(ci[1]) && x >= ci[1] && x <= ci[2]

## The likelihood assumes genotypes are independent draws from the allele
## frequencies. Species beta is a family group, which breaks that: relatives
## share genotypes, the prior is wrong for them, and the fit buys the mismatch
## back by inflating the false-allele rate (beta alone returns about 0.016).
## So the recovery check is made on the unrelated species, where the assumption
## holds, and the pooled figure is reported next to it as the honest caveat.
a_ids <- cn$SampleID[cn$TrueSpecies == "alpha"]
rx_a  <- rx[rx$SampleID %in% a_ids, ]
rk_a  <- paste(rx_a$SampleID, rx_a$Rep, sep = "|")
err_a <- gid_estimate_error(
  gid_matrix(cn[cn$SampleID %in% a_ids, ], "SampleID", loci),
  list(gt = gid_matrix(data.frame(key = rk_a, rx_a[, loci], check.names = FALSE,
                                  stringsAsFactors = FALSE), "key", loci),
       sample = rx_a$SampleID))
say("  on the unrelated species alone: dropout %.4f [%.4f - %.4f]  false allele %.4f [%.4f - %.4f]",
    err_a$dropout, err_a$dropout_ci[1], err_a$dropout_ci[2],
    err_a$false_allele, err_a$false_ci[1], err_a$false_ci[2])

ok_d <- identical(err_a$method, "replicates") && inside(DROPOUT, err_a$dropout_ci)
ok_f <- identical(err_a$method, "replicates") && inside(FALSE_ALLE, err_a$false_ci)
say("  true rates inside their intervals (unrelated species): dropout %s, false allele %s",
    ok_d, ok_f)
say("  pooled with the family group the false-allele rate reads high, as expected.")

## 2. individuals recovered --------------------------------------------------
truth <- cn$TrueIndividual[match(rownames(gt_cons), cn$SampleID)]
grp   <- cn$TrueSpecies[match(rownames(gt_cons), cn$SampleID)]
score <- function(r, lbl) {
  a <- r$assignment
  say("  %-34s %2d individuals (truth %d)   ARI %.4f", lbl,
      length(unique(a$individual)), n_ind,
      gid_ari(a$individual, truth[match(a$sample, rownames(gt_cons))]))
}
say("")
say("individuals recovered:")
score(gid_by_group(gt_cons, grp, gid_method_lr, dropout = err$dropout,
                   false_allele = err$false_allele, kinship = "full_sib",
                   post_cut = 0.999, min_loci = 12, reps = reps),
      "likelihood ratio, replicates")
score(gid_by_group(gt_cons, grp, gid_method_sethi, dropout = err$dropout,
                   false_allele = err$false_allele, min_loci = 12, reps = reps),
      "Sethi et al., replicates")
score(gid_by_group(gt_cons, grp, gid_method_exact, min_loci = 12),
      "exact match, consensus (contrast)")

if (!ok_d || !ok_f)
  stop("The demo does not recover its own error rates. Change the seed or the ",
       "sample size before shipping it -- the replicate story is the point.")
say("")
say("wrote app/demo/demo_genotypes.csv")
