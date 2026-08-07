## run_analysis.R ------------------------------------------------------------
##
## Collapse the AITRC July 2026 consensus genotypes into unique individuals
## using four independent methods, and compare them.
##
##   Method 1  exact match                (baseline; no error tolerance)
##   Method 2  mismatch threshold + sweep (Paetkau-style; the field standard)
##   Method 3  GenAlEx Matches            (Peakall & Smouse; what many labs run)
##   Method 4  allelematch                (Galpern et al. 2012; semi-automated)
##   Method 5  likelihood ratio           (probabilistic; recommended)
##
## Run from the genoID/ directory:  Rscript run_analysis.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(igraph); library(readxl); library(ggplot2)
})
source("app/genoID_core.R")

DATA_DIR <- ".."
CSV  <- file.path(DATA_DIR, "AITRC_DietMetabarcoding_Genotyping_Results_July2026_consensus.csv")
XLSX <- file.path(DATA_DIR, "AITRC_DietMetabarcoding_Genotyping_Results_July2026.xlsx")
OUT  <- "outputs"; dir.create(OUT, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")
tsv <- function(x, f) utils::write.csv(x, file.path(OUT, f), row.names = FALSE)

## ===========================================================================
## 1. LOAD AND DEFINE THE PANEL
## ===========================================================================
raw <- gid_read(CSV)
say("Loaded %d rows x %d columns", nrow(raw), ncol(raw))

META      <- c("PlateID", "SampleID", "Rep", "Raw_Reads", "On_Target_Reads",
               "Percent_On_Target", "IFI", "Pct_success",
               "Species_check_1", "Species_check_2", "Species_check_3")
SEX_LOCUS <- "OmyY1_2SEXY"
DIAG_LOCI <- c("CL07", "CL08", "CL09")   # species-diagnostic -> excluded from ID

all_loci <- gid_detect_loci(raw, exclude = META)
id_loci  <- setdiff(all_loci, c(SEX_LOCUS, DIAG_LOCI))
say("Detected %d genotype columns; %d used for individual ID", length(all_loci), length(id_loci))
say("  excluded: %s (species-diagnostic), %s (sex)",
    paste(DIAG_LOCI, collapse = ", "), SEX_LOCUS)

## Sample IDs are not always unique: a sample re-run on a second plate appears
## twice. Keep both rows -- they are a blind positive control, since any correct
## method must group them -- but give each row a unique key.
raw$row_key <- ifelse(duplicated(raw$SampleID) | duplicated(raw$SampleID, fromLast = TRUE),
                      paste0(raw$SampleID, "|", raw$PlateID), raw$SampleID)
dup_ids <- unique(raw$SampleID[duplicated(raw$SampleID)])
say("Duplicate SampleIDs (built-in positive controls): %s",
    if (length(dup_ids)) paste(dup_ids, collapse = ", ") else "none")

gt_all <- gid_matrix(raw, "row_key", id_loci)          # ID panel
gt_dia <- gid_matrix(raw, "row_key", DIAG_LOCI)        # species panel
sex    <- setNames(gid_norm_gt(raw[[SEX_LOCUS]]), raw$row_key)

## Allele-order normalisation: how many cells would have produced a false
## mismatch if compared as raw strings?
raw_str  <- toupper(gsub("[*?!#]", "", as.matrix(raw[, all_loci])))
norm_str <- gsub("/", "", gid_matrix(raw, "row_key", all_loci))
n_reordered <- sum(raw_str != norm_str & !is.na(norm_str))
say("Genotype cells whose allele order was normalised (%s): %d",
    paste(unique(all_loci[which(colSums(raw_str != norm_str & !is.na(norm_str)) > 0)]),
          collapse = ", "), n_reordered)
say("Quality-flagged (*) genotype cells retained: %d",
    sum(grepl("[*]", as.matrix(raw[, all_loci]))))

## ===========================================================================
## 2. SPECIES ASSIGNMENT (from the three diagnostic loci)
## ===========================================================================
## The lab's own Species_check_1/2/3 columns are simply CL07/CL08/CL09 read
## as wolf vs coyote. Take the majority call across the three.
sp_call <- t(apply(as.data.frame(raw[, c("Species_check_1", "Species_check_2", "Species_check_3")]),
                   1, function(z) gsub("[*]", "", tolower(trimws(as.character(z))))))
species <- apply(sp_call, 1, function(z) {
  z <- z[z %in% c("wolf", "coyote")]
  if (!length(z)) return("unknown")
  names(sort(table(z), decreasing = TRUE))[1]
})
names(species) <- raw$row_key
say("Species: %s", paste(sprintf("%s=%d", names(table(species)), table(species)), collapse = "  "))

## ===========================================================================
## 3. QC
## ===========================================================================
lstat <- gid_locus_stats(gt_all)
sstat <- gid_sample_stats(gt_all)
sstat$species <- species[sstat$sample]
sstat$IFI     <- as.numeric(raw$IFI)[match(sstat$sample, raw$row_key)]
sstat$sex     <- sex[sstat$sample]
tsv(lstat, "qc_locus_stats.csv"); tsv(sstat, "qc_sample_stats.csv")
say("Locus call rate: min %.2f  median %.2f", min(lstat$call_rate), median(lstat$call_rate))
say("Sample call rate: min %.2f  median %.2f", min(sstat$call_rate), median(sstat$call_rate))
## Drop loci that cannot contribute: monomorphic, or barely typed.
dead <- lstat$locus[lstat$n_alleles < 2 | lstat$call_rate < 0.25]
if (length(dead)) {
  say("Dropping uninformative loci: %s", paste(dead, collapse = ", "))
  gt_all <- gt_all[, setdiff(colnames(gt_all), dead), drop = FALSE]
  id_loci <- colnames(gt_all)
  lstat <- gid_locus_stats(gt_all); sstat <- gid_sample_stats(gt_all)
  sstat$species <- species[sstat$sample]
  sstat$IFI <- as.numeric(raw$IFI)[match(sstat$sample, raw$row_key)]
  sstat$sex <- sex[sstat$sample]
  tsv(lstat, "qc_locus_stats.csv"); tsv(sstat, "qc_sample_stats.csv")
}
say("Panel used for individual ID: %d loci", ncol(gt_all))

## ===========================================================================
## 4. GENOTYPING ERROR RATES FROM THE REPLICATE DATA
## ===========================================================================
gsheet <- as.data.frame(read_excel(XLSX, sheet = "Genotyping"), check.names = FALSE)
reps   <- gsheet[gsheet$Rep %in% c("a", "b", "c"), ]
cons   <- gsheet[gsheet$Rep == "consensus", ]
cons   <- cons[!duplicated(cons$SampleID), ]

rep_gt  <- gid_matrix(reps, "SampleID", id_loci)
cons_gt <- gid_matrix(cons, "SampleID", id_loci)

## (a) Broquet & Petit estimator, replicate scored against the consensus.
err <- gid_error_rates(rep_gt, cons_gt, reps$SampleID, cons$SampleID)
pooled <- attr(err, "pooled")
tsv(err, "error_rates_per_replicate.csv")
say("Per-replicate error vs consensus (Broquet & Petit): dropout %.4f  false allele %.4f",
    pooled["dropout"], pooled["false_allele"])
say("  ^ the false-allele rate here is forced to 0 by the consensus rule itself")
say("    (a homozygous consensus call REQUIRES three identical replicates),")
say("    so it is not a usable estimate. Use the ML estimator below.")

## (b) Maximum likelihood from the replicate triplets alone. Never touches the
##     consensus, so the false-allele rate is actually estimable.
mle <- gid_error_ml(rep_gt, reps$SampleID, freqs = gid_allele_freq(cons_gt))
say("ML per-replicate error: dropout %.4f [%.4f-%.4f]  false allele %.4f [%.4f-%.4f]",
    mle["dropout"], mle["dropout_lo"], mle["dropout_hi"],
    mle["false_allele"], mle["false_lo"], mle["false_hi"])
pooled <- c(dropout = unname(mle["dropout"]), false_allele = unname(mle["false_allele"]))

## Residual error in the *consensus* genotypes, given the multi-tube rule the
## lab used (heterozygote accepted on 2 replicates, homozygote needs 3).
per_rep_missing <- mean(is.na(rep_gt))
resid <- gid_propagate_error(pooled["dropout"], pooled["false_allele"],
                             n_rep = 3, rule = "taberlet", hom_n = 3, het_n = 2,
                             per_rep_missing = per_rep_missing, nsim = 50000)
say("Per-replicate missingness: %.3f", per_rep_missing)
say("Residual error in CONSENSUS calls: dropout %.4f  false allele %.4f",
    resid["dropout"], resid["false_allele"])
D_CONS <- max(resid["dropout"], 1e-4); F_CONS <- max(resid["false_allele"], 1e-4)
write.csv(data.frame(
  level = c("per_replicate_BroquetPetit", "per_replicate_ML", "consensus_residual"),
  dropout = c(attr(err, "pooled")["dropout"], mle["dropout"], D_CONS),
  false_allele = c(attr(err, "pooled")["false_allele"], mle["false_allele"], F_CONS)),
  file.path(OUT, "error_rates_summary.csv"), row.names = FALSE)

## ===========================================================================
## 5. ALLELE FREQUENCIES AND P(ID)   (wolves only; recapture-corrected)
## ===========================================================================
wolf <- names(species)[species == "wolf"]
gt_w <- gt_all[wolf, , drop = FALSE]

## Pass 1: naive frequencies -> provisional clusters -> down-weight recaptures
## -> Pass 2 frequencies. Otherwise a heavily resampled wolf inflates its own
## alleles and P(ID) looks better than it is.
f1  <- gid_allele_freq(gt_w)
prov <- gid_method_lr(gt_w, freqs = f1, dropout = D_CONS, false_allele = F_CONS,
                      kinship = "full_sib", post_cut = 0.999, min_loci = 15)
csz  <- table(prov$assignment$individual)
wts  <- 1 / as.numeric(csz[prov$assignment$individual])
freqs <- gid_allele_freq(gt_w, weights = wts[match(rownames(gt_w), prov$assignment$sample)])

pid <- gid_pid(freqs)
pid <- pid[order(pid$pid), ]                  # most informative loci first
pid$pid_cum <- cumprod(pid$pid); pid$pid_sib_cum <- cumprod(pid$pid_sib)
tsv(pid, "pid_by_locus.csv")
say("Full panel (%d loci): P(ID) = %.3g   P(ID)sib = %.3g",
    nrow(pid), tail(pid$pid_cum, 1), tail(pid$pid_sib_cum, 1))
say("Loci needed for P(ID)sib < 0.01: %d", which(pid$pid_sib_cum < 0.01)[1])
say("Loci needed for P(ID)sib < 0.001: %d",
    ifelse(any(pid$pid_sib_cum < 0.001), which(pid$pid_sib_cum < 0.001)[1], NA))

## Per-sample resolving power, using only the loci each sample actually has
## Fis, linkage and allele frequencies are all meaningless while recaptures
## are still in the table: one wolf sampled 16 times contributes its genotype
## 16 times. Everything below uses one sample per provisional individual.
uniq_rep <- prov$assignment$sample[!duplicated(prov$assignment$individual)]
gt_u <- gt_w[uniq_rep, , drop = FALSE]

## Linked loci. Six markers are a second SNP from the same amplicon as another
## marker (X and XSNP2), so P(ID) must not multiply them as independent.
ld <- gid_ld(gt_u)
tsv(ld, "linkage_r2.csv")
say("Locus pairs with r^2 > 0.3 among unique individuals: %s",
    paste(sprintf("%s~%s(%.2f)", ld$locus1, ld$locus2, ld$r2)[ld$r2 > 0.3], collapse = " "))
say("  same-amplicon pairs: %s",
    paste(sprintf("%s~%s r2=%.2f", ld$locus1, ld$locus2, ld$r2)[
      sub("SNP2$", "", ld$locus1) == sub("SNP2$", "", ld$locus2)], collapse = "; "))

lstat_ind <- gid_locus_stats(gt_u)
tsv(lstat_ind, "qc_locus_stats_unique_individuals.csv")
say("Loci with |Fis| > 0.3 among unique individuals (n=%d): %s", length(uniq_rep),
    paste(sprintf("%s(%+.2f)", lstat_ind$locus, lstat_ind$Fis)[abs(lstat_ind$Fis) > 0.3 &
          !is.na(lstat_ind$Fis)], collapse = " "))

## P(ID) with one SNP per amplicon, since SNP2 markers are physically linked
paired <- intersect(sub("SNP2$", "", grep("SNP2$", colnames(gt_w), value = TRUE)),
                    colnames(gt_w))
pruned <- setdiff(colnames(gt_w), grep("SNP2$", colnames(gt_w), value = TRUE))
pid_pruned <- gid_pid(freqs[pruned])
say("Linkage-pruned panel (%d loci, one SNP per amplicon): P(ID) = %.3g  P(ID)sib = %.3g",
    length(pruned), tail(pid_pruned$pid_cum, 1), tail(pid_pruned$pid_sib_cum, 1))
tsv(pid_pruned, "pid_by_locus_linkage_pruned.csv")

pw <- gid_sample_power(gt_w, pid)
pw$call_rate <- sstat$call_rate[match(pw$sample, sstat$sample)]
tsv(pw, "sample_power.csv")
say("Samples with P(ID)sib > 0.01 given their own typed loci: %d",
    sum(pw$pid_sib > 0.01))

## ===========================================================================
## 6. THE FOUR METHODS   (blocked by species throughout)
## ===========================================================================
enc <- gid_encode(gt_all)
MIN_LOCI <- 15

say("\n--- Method 1: exact match ---")
m1 <- gid_by_group(gt_all, species, gid_method_exact, min_loci = MIN_LOCI)

say("--- Method 2: mismatch threshold ---")
sweep_tab <- gid_threshold_sweep(gt_w, max_k = 6, min_loci = MIN_LOCI)
tsv(sweep_tab, "threshold_sweep.csv"); print(sweep_tab)
## choose the threshold at the plateau: the largest k for which the number of
## individuals is unchanged from k-1 and no cluster is internally inconsistent
plateau <- with(sweep_tab, which(c(FALSE, diff(n_individuals) == 0) &
                                   n_conflicting_clusters == 0))
K <- if (length(plateau)) sweep_tab$max_mismatch[max(plateau)] else 1
say("Chosen mismatch threshold from the sweep: %d", K)
m2 <- gid_by_group(gt_all, species, gid_method_threshold,
                   max_mismatch = K, min_loci = MIN_LOCI, max_mismatch_2allele = 0)

say("--- Method 3: GenAlEx Matches ---")
mgx <- gid_by_group(gt_all, species, gid_method_genalex, min_loci = MIN_LOCI,
                    near_match_loci = 2)
gxw <- mgx$by_group$wolf
say("GenAlEx: %d exact-match pairs, %d near matches; P(ID)=%.3g P(ID)sib=%.3g",
    nrow(gxw$matched_pairs), nrow(gxw$near_matches),
    gxw$pid_total[["pid"]], gxw$pid_total[["pid_sib"]])
tsv(gxw$match_distribution, "genalex_match_distribution.csv")
tsv(gxw$near_matches[, c("id1", "id2", "n_compared", "n_mismatch")], "genalex_near_matches.csv")

say("--- Method 4: allelematch ---")
m3 <- gid_by_group(gt_all, species, gid_method_allelematch, min_loci = MIN_LOCI)
am <- m3$by_group$wolf
say("allelematch engine: %s  alleleMismatch = %s  matchThreshold = %s  cutHeight = %s",
    am$engine, am$settings$alleleMismatch, signif(am$settings$matchThreshold %||% NA, 3),
    signif(am$settings$cutHeight %||% NA, 3))
if (!is.null(am$profile) && is.data.frame(am$profile)) {
  cat("allelematch alleleMismatch profile (its own threshold-selection curve):\n")
  print(am$profile); tsv(am$profile, "allelematch_profile.csv")
}
say("allelematch flagged %s samples as unclassified or multiple-match",
    length(am$unclassified %||% character(0)))

say("--- Method 5: likelihood ratio ---")
m4 <- gid_by_group(gt_all, species, gid_method_lr, freqs = NULL,
                   dropout = D_CONS, false_allele = F_CONS,
                   kinship = "full_sib", post_cut = 0.999, min_loci = MIN_LOCI)
m4u <- gid_by_group(gt_all, species, gid_method_lr, freqs = NULL,
                    dropout = D_CONS, false_allele = F_CONS,
                    kinship = "unrelated", post_cut = 0.999, min_loci = MIN_LOCI)

results <- list(exact = m1, threshold = m2, genalex = mgx, allelematch = m3,
                LR_fullsib = m4, LR_unrelated = m4u)

## ===========================================================================
## 7. COMPARISON
## ===========================================================================
cmp <- gid_compare_methods(results)
print(cmp$summary); cat("\nAdjusted Rand Index between methods:\n"); print(round(cmp$ari, 3))
tsv(cmp$summary, "method_summary.csv")
tsv(as.data.frame(round(cmp$ari, 4)), "method_ari.csv")
tsv(cmp$table, "assignments_all_methods.csv")

## Positive control: the two plate runs of a duplicated sample must land together.
if (length(dup_ids)) {
  ctrl <- do.call(rbind, lapply(names(results), function(nm) {
    a <- results[[nm]]$assignment
    k <- a$individual[grepl(dup_ids[1], a$sample, fixed = TRUE)]
    data.frame(method = nm, n_rows = length(k), same_individual = length(unique(k)) == 1)
  }))
  cat("\nPositive control (", dup_ids[1], " run on two plates):\n", sep = "")
  print(ctrl); tsv(ctrl, "positive_control.csv")
}

## Sex-marker consistency inside each individual (an independent check that
## never entered the matching)
sex_check <- do.call(rbind, lapply(names(results), function(nm) {
  a <- results[[nm]]$assignment; a$sex <- sex[a$sample]
  bad <- tapply(a$sex, a$individual, function(z) {
    z <- z[!is.na(z)]; length(unique(z)) > 1 })
  data.frame(method = nm, n_individuals = length(bad),
             n_sex_conflicts = sum(bad, na.rm = TRUE))
}))
cat("\nSex-marker conflicts within inferred individuals:\n"); print(sex_check)
tsv(sex_check, "sex_consistency.csv")

## Cluster consistency: are all members of every cluster mutually matched?
## A cluster held together by single linkage but not internally complete is
## the classic chaining failure, and worth looking at by hand every time.
cat("\nInternally inconsistent clusters (single-linkage chaining):\n")
print(do.call(rbind, lapply(names(results), function(nm)
  data.frame(method = nm, n_conflicting_clusters = results[[nm]]$n_conflict))))
cf_all <- do.call(rbind, lapply(names(results), function(nm) {
  x <- results[[nm]]$conflicts; if (is.null(x)) NULL else cbind(method = nm, x) }))
if (!is.null(cf_all)) {
  cat("\nDetails:\n"); print(cf_all[, c("method", "n_samples", "completeness", "members")])
  tsv(cf_all, "cluster_conflicts.csv")
}

## Complete-linkage version of the recommended method, as a stability check
m4c <- gid_by_group(gt_all, species, gid_method_lr, freqs = NULL,
                    dropout = D_CONS, false_allele = F_CONS, kinship = "full_sib",
                    post_cut = 0.999, min_loci = MIN_LOCI, linkage = "complete")
say("LR with complete linkage instead of single: %d individuals (vs %d)",
    length(unique(m4c$assignment$individual)), length(unique(m4$assignment$individual)))
say("  ARI between the two linkages: %.4f",
    gid_ari(m4$assignment$individual, m4c$assignment$individual[match(m4$assignment$sample, m4c$assignment$sample)]))

## Where the methods disagree, sample by sample
say("\nSamples the five methods do not co-assign identically: %d of %d",
    cmp$n_disagree, nrow(cmp$table))
disagree <- cmp$table[cmp$table$disputed, ]
tsv(disagree, "method_disagreements.csv")
if (nrow(disagree)) {
  cat("\nDisputed samples (grouped by what the recommended method did):\n")
  d <- disagree
  d$call_rate <- round(sstat$call_rate[match(d$sample, sstat$sample)], 3)
  print(d[order(d$LR_fullsib), c("sample", "call_rate", "exact", "threshold",
                                 "allelematch", "LR_fullsib")], row.names = FALSE)
}

## ===========================================================================
## 8. FINAL ANSWER FROM THE RECOMMENDED METHOD
## ===========================================================================
final <- m4$assignment
final$species    <- species[final$sample]
final$sex        <- sex[final$sample]
final$call_rate  <- sstat$call_rate[match(final$sample, sstat$sample)]
final$IFI        <- sstat$IFI[match(final$sample, sstat$sample)]
final$pid_sib    <- pw$pid_sib[match(final$sample, pw$sample)]
final$n_samples_for_individual <- as.integer(table(final$individual)[final$individual])
final <- final[order(final$individual, -final$call_rate), ]
tsv(final, "FINAL_individual_assignments.csv")

pairs_w <- m4$by_group$wolf$pairs
tsv(pairs_w[order(-pairs_w$log10_LR), ][1:200, ], "top_pairs_LR.csv")
tsv(m4$matched_pairs, "matched_pairs_LR.csv")

ic <- gid_individual_consensus(gt_all, m4$assignment)
write.csv(data.frame(individual = rownames(ic$genotypes), ic$genotypes),
          file.path(OUT, "individual_consensus_genotypes.csv"), row.names = FALSE)
if (!is.null(ic$discordance)) tsv(ic$discordance, "within_individual_discordance.csv")

s <- gid_summarise_assignment(m4$assignment)
say("\n=== RECOMMENDED (LR, full-sib null, P>=0.999) ===")
say("%d samples -> %d individuals (%d wolf, %d coyote)",
    s$n_samples, s$n_individuals,
    length(unique(m4$assignment$individual[species[m4$assignment$sample] == "wolf"])),
    length(unique(m4$assignment$individual[species[m4$assignment$sample] == "coyote"])))
say("Sampled once: %d   resampled: %d   largest cluster: %d",
    s$n_singletons, s$n_individuals - s$n_singletons, s$max_cluster)
say("Recapture rate: %.1f%%", 100 * s$recapture_rate)

## ===========================================================================
## 9. FIGURES
## ===========================================================================
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold")))
ggsave2 <- function(p, f, w = 7, h = 4.4)
  suppressWarnings(ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 200))

pp <- m4$by_group$wolf$pairs
pp <- pp[pp$n_compared >= MIN_LOCI, ]
raw_mm <- gid_method_threshold(gt_w, max_mismatch = 0, min_loci = MIN_LOCI)$pairs
raw_mm <- raw_mm[raw_mm$n_compared >= MIN_LOCI, ]

ggsave2(ggplot(raw_mm, aes(n_mismatch)) +
  geom_histogram(binwidth = 1, fill = "#3b6ea5", colour = "white") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
  geom_vline(xintercept = K + 0.5, linetype = 2, colour = "#c0392b") +
  labs(title = "Pairwise locus mismatches (wolves)",
       subtitle = sprintf("Dashed line = chosen threshold (<=%d mismatches). Log-scaled count axis.", K),
       x = "Mismatching loci", y = "Number of sample pairs"),
  "fig1_mismatch_distribution.png")

ggsave2(ggplot(pp, aes(log10_LR)) +
  geom_histogram(bins = 70, fill = "#3b6ea5", colour = NA) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  labs(title = "Log10 likelihood ratio, same individual vs full siblings",
       subtitle = "Wolves only. Positive = evidence the two samples are one animal.",
       x = expression(log[10]~LR), y = "Number of sample pairs"),
  "fig2_LR_distribution.png")

pid_long <- rbind(
  data.frame(k = seq_len(nrow(pid)), value = pid$pid_cum,     stat = "P(ID)"),
  data.frame(k = seq_len(nrow(pid)), value = pid$pid_sib_cum, stat = "P(ID)sib"))
ggsave2(ggplot(pid_long, aes(k, value, colour = stat)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
  scale_y_log10() + scale_colour_manual(values = c("#3b6ea5", "#c0392b")) +
  geom_hline(yintercept = 0.01, linetype = 2, colour = "grey40") +
  labs(title = "Cumulative probability of identity",
       subtitle = "Loci added most-informative first. Dashed line = 0.01 benchmark.",
       x = "Number of loci", y = "Cumulative probability", colour = NULL),
  "fig3_pid_accumulation.png")

ggsave2(ggplot(sweep_tab, aes(max_mismatch, n_individuals)) +
  geom_line(colour = "#3b6ea5", linewidth = 0.9) +
  geom_point(aes(size = n_conflicting_clusters + 1), colour = "#3b6ea5") +
  geom_vline(xintercept = K, linetype = 2, colour = "#c0392b") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "Sensitivity of the answer to the mismatch threshold",
       subtitle = "Point size = number of internally inconsistent clusters. Flat = stable.",
       x = "Maximum mismatching loci allowed", y = "Individuals inferred"),
  "fig4_threshold_sweep.png")

ggsave2(ggplot(pw, aes(call_rate, pid_sib)) +
  geom_point(colour = "#3b6ea5", alpha = 0.75, size = 2) +
  scale_y_log10() + geom_hline(yintercept = 0.01, linetype = 2, colour = "#c0392b") +
  labs(title = "Per-sample resolving power",
       subtitle = "P(ID)sib computed from the loci each sample actually has. Above the line = not safely identifiable.",
       x = "Sample call rate", y = "P(ID)sib for this sample"),
  "fig5_sample_power.png")

saveRDS(list(results = results, cmp = cmp, pid = pid, err = err, resid = resid,
             freqs = freqs, species = species, sex = sex, gt = gt_all,
             sweep = sweep_tab, power = pw, final = final),
        file.path(OUT, "analysis_objects.rds"))
say("\nWrote %d files to %s/", length(list.files(OUT)), OUT)
