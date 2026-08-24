## methods_disagreement.R --------------------------------------------------------
##
## Where do the matching methods actually disagree, and how much of the spread is
## driven by allelematch alone?
##
## The headline comparison reports a range of individual counts, but a count is a
## weak summary: two methods can return the same number of individuals while
## disagreeing about which samples go together. This works at the level of sample
## PAIRS instead, which is where a disagreement is a fact rather than an artifact
## of how the clusters happened to fall out.
##
## Run from genoID/:  Rscript methods_disagreement.R
## -----------------------------------------------------------------------------
source("app/genoID_core.R")

OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
MIN_LOCI <- 15
POST_CUT <- 0.999

## ------------------------------------------------------------------- setup ---
o   <- readRDS(file.path(OUT, "wolf_individuals.rds"))
gt  <- o$gt; species <- o$species; err <- o$err
d   <- gid_read("../AITRC_DietMetabarcoding_Genotyping_Results_July2026_rawreps.csv")
rd  <- d[d$Rep %in% c("a", "b", "c"), ]
reps <- list(gt = gid_matrix(rd, "SampleID", colnames(gt)), sample = rd$SampleID)

## The consensus route needs the error rates that survive the 3-replicate rule,
## not the per-reaction rates -- feeding per-reaction rates to a consensus would
## double-count the error suppression the consensus already bought.
pr <- gid_propagate_error(err$dropout, err$false_allele, n_rep = 3, rule = "taberlet",
                          hom_n = 3, het_n = 2,
                          per_rep_missing = mean(is.na(reps$gt)), nsim = 30000)

methods <- list(
  "LR, replicates"  = function() gid_by_group(gt, species, gid_method_lr,
      dropout = err$dropout, false_allele = err$false_allele, kinship = "full_sib",
      post_cut = POST_CUT, min_loci = MIN_LOCI, reps = reps),
  "LR, consensus"   = function() gid_by_group(gt, species, gid_method_lr,
      dropout = max(pr[["dropout"]], 1e-4), false_allele = max(pr[["false_allele"]], 1e-4),
      kinship = "full_sib", post_cut = POST_CUT, min_loci = MIN_LOCI),
  "Exact / GenAlEx" = function() gid_by_group(gt, species, gid_method_genalex, min_loci = MIN_LOCI),
  "Sethi 2016"      = function() gid_by_group(gt, species, gid_method_sethi,
      dropout = err$dropout, false_allele = err$false_allele, lambda_cut = 1,
      min_loci = MIN_LOCI, reps = reps),
  "Threshold k=1"   = function() gid_by_group(gt, species, gid_method_threshold,
      max_mismatch = 1, min_loci = MIN_LOCI),
  "allelematch"     = function() gid_by_group(gt, species, gid_method_allelematch,
      min_loci = MIN_LOCI))

runs <- lapply(methods, function(f) f())
ids  <- rownames(gt)
## one column per method, one row per sample, holding that method's cluster label
part <- sapply(runs, function(r) r$assignment$individual[match(ids, r$assignment$sample)])
rownames(part) <- ids

CORE <- setdiff(names(methods), "allelematch")

## --------------------------------------------------- 1. pairwise agreement ---
ari <- outer(seq_along(methods), seq_along(methods),
             Vectorize(function(i, j) gid_ari(part[, i], part[, j])))
dimnames(ari) <- list(names(methods), names(methods))

say("=== pairwise agreement (Adjusted Rand Index) ===")
print(round(ari, 3))
off <- function(m) m[upper.tri(m)]
say("\nmean pairwise ARI, all 6 methods : %.4f  (range %.3f - %.3f)",
    mean(off(ari)), min(off(ari)), max(off(ari)))
say("mean pairwise ARI, without allelematch: %.4f  (range %.3f - %.3f)",
    mean(off(ari[CORE, CORE])), min(off(ari[CORE, CORE])), max(off(ari[CORE, CORE])))
say("allelematch vs the other five        : %.4f  (range %.3f - %.3f)",
    mean(ari["allelematch", CORE]), min(ari["allelematch", CORE]), max(ari["allelematch", CORE]))

## ------------------------------------------------------ 2. contested pairs ---
## For every pair of samples, count how many methods put them together. A pair
## voted together by all methods, or by none, is settled; anything between is a
## genuine disagreement about the data.
co <- function(lab) { m <- outer(lab, lab, "=="); m[is.na(m)] <- FALSE; m }
votes_all  <- Reduce(`+`, lapply(seq_along(methods), function(k) co(part[, k])))
votes_core <- Reduce(`+`, lapply(CORE, function(k) co(part[, k])))

pair_tab <- function(votes, n_methods, who) {
  ut <- which(upper.tri(votes) & votes > 0 & votes < n_methods, arr.ind = TRUE)
  if (!nrow(ut)) return(NULL)
  data.frame(
    sample_a = ids[ut[, 1]], sample_b = ids[ut[, 2]],
    n_together = votes[ut], n_methods = n_methods,
    methods_joining = apply(ut, 1, function(k)
      paste(who[part[k[1], who] == part[k[2], who]], collapse = "; ")),
    stringsAsFactors = FALSE)
}
cp_all  <- pair_tab(votes_all,  length(methods), names(methods))
cp_core <- pair_tab(votes_core, length(CORE),    CORE)

settled_all  <- sum(upper.tri(votes_all)  & (votes_all  == 0 | votes_all  == length(methods)))
settled_core <- sum(upper.tri(votes_core) & (votes_core == 0 | votes_core == length(CORE)))
n_pairs <- choose(length(ids), 2)

say("\n=== contested sample pairs (together under some methods, not others) ===")
say("all 6 methods        : %d contested of %d pairs (%.3f%% unsettled)",
    nrow(cp_all), n_pairs, 100 * nrow(cp_all) / n_pairs)
say("without allelematch  : %d contested of %d pairs (%.3f%% unsettled)",
    nrow(cp_core), n_pairs, 100 * nrow(cp_core) / n_pairs)
say("pairs contested ONLY because of allelematch: %d",
    nrow(cp_all) - nrow(cp_core))

if (!is.null(cp_core)) {
  cp_core <- cp_core[order(-cp_core$n_together, cp_core$sample_a), ]
  say("\nthe %d pairs the five core methods argue over:", nrow(cp_core))
  for (i in seq_len(nrow(cp_core))) say("  %s + %s  -- joined by %d/5: %s",
      cp_core$sample_a[i], cp_core$sample_b[i], cp_core$n_together[i], cp_core$methods_joining[i])
}

## -------------------------------------------- 3. how many samples are stable --
## A sample is stable if every method places it with exactly the same partners.
stable <- sapply(seq_along(ids), function(i)
  length(unique(apply(part[, CORE, drop = FALSE], 2,
                      function(lab) paste(sort(ids[lab == lab[i] & !is.na(lab)]), collapse = ",")))) == 1)
stable_all <- sapply(seq_along(ids), function(i)
  length(unique(apply(part, 2,
                      function(lab) paste(sort(ids[lab == lab[i] & !is.na(lab)]), collapse = ",")))) == 1)
say("\n=== sample stability ===")
say("identically grouped by all 6 methods       : %d of %d (%.1f%%)",
    sum(stable_all), length(ids), 100 * mean(stable_all))
say("identically grouped by the five core methods: %d of %d (%.1f%%)",
    sum(stable), length(ids), 100 * mean(stable))
say("samples destabilised only by allelematch    : %d", sum(stable) - sum(stable_all))

## ------------------------------------------------ 4. majority-rule consensus --
## Take the pairs a majority of the five core methods agree on, then resolve the
## resulting graph the same way any single method's output is resolved.
maj  <- votes_core >= ceiling(length(CORE) / 2)
pm   <- which(upper.tri(maj) & maj, arr.ind = TRUE)
## gid_resolve() keys on id1/id2; anything else matches to NA and silently
## returns every sample as its own individual.
cons <- gid_resolve(data.frame(id1 = ids[pm[, 1]], id2 = ids[pm[, 2]],
                               stringsAsFactors = FALSE), ids, linkage = "complete")
stopifnot(nrow(pm) > 0, length(unique(cons$assignment$individual)) < length(ids))
cw <- length(unique(cons$assignment$individual[species[cons$assignment$sample] == "wolf"]))
say("\n=== majority-rule consensus of the five core methods ===")
say("%d individuals (%d wolves), internally inconsistent clusters: %d",
    length(unique(cons$assignment$individual)), cw, cons$n_conflict)
say("agreement with the reported LR/replicates result: ARI %.4f",
    gid_ari(part[, "LR, replicates"], cons$assignment$individual[match(ids, cons$assignment$sample)]))

## ------------------------------------------------------------ 5. the counts ---
cnt <- data.frame(
  method = names(methods),
  individuals = sapply(runs, function(r) length(unique(r$assignment$individual))),
  wolves = sapply(runs, function(r) length(unique(
    r$assignment$individual[species[r$assignment$sample] == "wolf"]))),
  inconsistent_clusters = sapply(runs, function(r) r$n_conflict),
  stringsAsFactors = FALSE)
say("\n=== individual counts ===")
print(cnt, row.names = FALSE)
say("\nspread including allelematch : %d - %d individuals (%d wide)",
    min(cnt$individuals), max(cnt$individuals), diff(range(cnt$individuals)))
core_cnt <- cnt$individuals[cnt$method %in% CORE]
say("spread excluding allelematch : %d - %d individuals (%d wide)",
    min(core_cnt), max(core_cnt), diff(range(core_cnt)))

## ------------------------------------------ 6. are the singletons for real? ---
## 48 wolves seen exactly once is a lot. The question is whether they are truly
## distinct animals or whether the cutoff is splitting resampled ones apart, and
## the way to tell is to ask each singleton how close its nearest neighbour is.
## A singleton whose best rival sits at 1e-20 is not a borderline call.
cf <- o$conf
pr <- do.call(rbind, lapply(names(res_by <- o$res$by_group), function(g) res_by[[g]]$pairs))
sing <- cf$sample[cf$n_in_individual == 1 & cf$species == "wolf"]
nn <- do.call(rbind, lapply(sing, function(s) {
  r <- pr[pr$id1 == s | pr$id2 == s, ]
  r <- r[order(-r$posterior_same), ][1, ]
  other <- if (r$id1 == s) r$id2 else r$id1
  data.frame(sample = s, animal = cf$animal[cf$sample == s],
             nearest = other, nearest_animal = cf$animal[match(other, cf$sample)],
             posterior = r$posterior_same, log10_LR = r$log10_LR,
             n_loci = r$n_compared, stringsAsFactors = FALSE)
}))
nn <- nn[order(-nn$posterior), ]

say("\n=== how isolated is each of the %d singleton wolves? ===", nrow(nn))
brk <- c(-Inf, 1e-12, 1e-6, 0.001, 0.01, 0.1, 0.5, Inf)
lbl <- c("below 1e-12 (decisive)", "1e-12 to 1e-6", "1e-6 to 0.001", "0.001 to 0.01",
         "0.01 to 0.1", "0.1 to 0.5", "above 0.5 (near-match)")
print(table(`nearest-neighbour posterior` = cut(nn$posterior, brk, labels = lbl)))
say("median nearest-neighbour posterior: %.3g   median log10 LR: %.1f",
    median(nn$posterior), median(nn$log10_LR))
say("\nsingletons with a nearest neighbour worth a second look (posterior > 1e-6):")
print(format(nn[nn$posterior > 1e-6, ], digits = 3), row.names = FALSE)
write.csv(nn, file.path(OUT, "WOLF_singleton_nearest.csv"), row.names = FALSE)

## The permissive end of the method range is exactly "merge every near-match":
## check that, rather than asserting it.
near <- unique(unlist(lapply(which(nn$posterior > 0.1), function(k)
  sort(c(nn$sample[k], nn$nearest[k])))))
n_ref <- cnt$wolves[cnt$method == "LR, replicates"]
say("\nthe %d samples involved in a near-match form %d candidate merges;",
    length(near), length(near) / 2)
say("merging all of them would take the reported %d wolves to %d;",
    n_ref, n_ref - length(near) / 2)
say("Sethi and the k=1 threshold return %d and %d.",
    cnt$wolves[cnt$method == "Sethi 2016"], cnt$wolves[cnt$method == "Threshold k=1"])

write.csv(cbind(method = rownames(ari), as.data.frame(round(ari, 4))),
          file.path(OUT, "WOLF_method_ari.csv"), row.names = FALSE)
if (!is.null(cp_all)) write.csv(cp_all, file.path(OUT, "WOLF_contested_pairs.csv"), row.names = FALSE)
saveRDS(list(part = part, ari = ari, cp_all = cp_all, cp_core = cp_core,
             cons = cons, cnt = cnt, stable = stable, stable_all = stable_all),
        file.path(OUT, "wolf_disagreement.rds"))
say("\nwrote WOLF_method_ari.csv, WOLF_contested_pairs.csv")
