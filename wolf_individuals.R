## wolf_individuals.R ---------------------------------------------------------
##
## The deliverable run: assign every AITRC sample to an individual wolf, with a
## confidence measure per sample, and cross-check the answer against field
## metadata that never entered the genetics.
##
## Settings are not conventional defaults. Error rates are measured from the PCR
## replicates by maximum likelihood, the analysis conditions on every replicate
## rather than a consensus, and the posterior cutoff is chosen from the expected
## number of mistakes it would make.
##
## Run from genoID/:  Rscript wolf_individuals.R
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(readxl))
source("app/genoID_core.R")

RAW  <- "../AITRC_DietMetabarcoding_Genotyping_Results_July2026_rawreps.csv"
XLSX <- "../AITRC_DietMetabarcoding_Genotyping_Results_July2026.xlsx"
OUT  <- "outputs"; dir.create(OUT, showWarnings = FALSE)
say  <- function(...) cat(sprintf(...), "\n")
POST_CUT <- 0.999      # sits in the middle of the plateau; see calibration below
MIN_LOCI <- 15

## ---------------------------------------------------------------- 1. load ---
d    <- gid_read(RAW)
loci <- setdiff(gid_detect_loci(d), c("CL07", "CL08", "CL09"))   # species-diagnostic
rd   <- d[d$Rep %in% c("a", "b", "c"), ]
cn   <- d[d$Rep == "consensus", ]

## One sample was run on two plates and so appears twice. Give the rows unique
## keys: duplicate matrix rownames silently alias to the first match, which
## would drop the second copy from every comparison -- and that copy is the
## blind positive control, so losing it defeats the one check that matters.
cn$key <- ave(cn$SampleID, cn$SampleID, FUN = function(z)
  if (length(z) == 1) z else paste0(z, "|", seq_along(z)))
dup_key <- cn$key[cn$SampleID %in% names(which(table(cn$SampleID) > 1))]

gt   <- gid_matrix(cn, "key", loci)
## replicates are keyed by SampleID, so a duplicated sample shares them; that is
## correct, both plate runs are observations of the same DNA
reps <- list(gt = gid_matrix(rd, "SampleID", loci), sample = rd$SampleID)
f    <- gid_filter(gt, min_locus_call = 0.25, min_sample_call = 0.5)
gt   <- f$gt
say("%d samples in file; %d pass the call-rate filter; %d loci used",
    nrow(cn), nrow(gt), ncol(gt))
say("dropped for low call rate: %s", paste(f$dropped_samples, collapse = ", "))

species <- setNames(gsub("[*]", "", tolower(cn$Species_check_3)), cn$key)[rownames(gt)]
sex     <- setNames(gid_norm_gt(cn$OmyY1_2SEXY), cn$key)[rownames(gt)]
say("species: %s", paste(sprintf("%s=%d", names(table(species)), table(species)), collapse = "  "))

## --------------------------------------------------- 2. measure the error ---
err <- gid_estimate_error(gt, reps)
say("\nerror rates measured from the replicates (%s):", err$method)
say("  dropout      %.4f  [%.4f - %.4f] per reaction", err$dropout, err$dropout_ci[1], err$dropout_ci[2])
say("  false allele %.4f  [%.4f - %.4f] per reaction", err$false_allele, err$false_ci[1], err$false_ci[2])

## ------------------------------------------------------------- 3. assign ---
run <- function(cut) gid_by_group(gt, species, gid_method_lr,
                                  dropout = err$dropout, false_allele = err$false_allele,
                                  kinship = "full_sib", post_cut = cut,
                                  min_loci = MIN_LOCI, reps = reps)
res   <- run(POST_CUT)
sethi <- gid_by_group(gt, species, gid_method_sethi,
                      dropout = err$dropout, false_allele = err$false_allele,
                      relationships = c("unrelated", "full_sib", "parent_offspring"),
                      lambda_cut = 1, min_loci = MIN_LOCI, reps = reps)

say("\nlikelihood ratio, full-sib null, P >= %s : %d individuals", POST_CUT,
    length(unique(res$assignment$individual)))
say("Sethi et al. (2016), lambda > 1              : %d individuals",
    length(unique(sethi$assignment$individual)))

## ------------------------------------------- 4. was the cutoff arbitrary? ---
cal <- do.call(rbind, lapply(names(res$by_group), function(g) {
  e <- res$by_group[[g]]
  tb <- gid_calibrate_threshold(e$pairs, e$assignment$sample, min_loci = MIN_LOCI)
  if (is.null(tb)) NULL else cbind(group = g, tb)
}))
write.csv(cal, file.path(OUT, "wolf_cutoff_calibration.csv"), row.names = FALSE)
wolfcal <- cal[cal$group == "wolf", ]
plateau <- wolfcal$cutoff[wolfcal$n_individuals == wolfcal$n_individuals[wolfcal$cutoff == POST_CUT]]
say("\ncutoff plateau: %d wolf individuals for every cutoff from %s to %s",
    wolfcal$n_individuals[wolfcal$cutoff == POST_CUT], min(plateau), max(plateau))
say("at the reported cutoff: %.3f expected false merges, %.1f expected missed pairs",
    wolfcal$exp_false_merges[wolfcal$cutoff == POST_CUT],
    wolfcal$exp_missed_pairs[wolfcal$cutoff == POST_CUT])

## ---------------------------------------------------------- 5. confidence ---
freqs <- gid_allele_freq(gt[names(species)[species == "wolf"], , drop = FALSE])
pid   <- gid_pid(freqs); pid <- pid[order(pid$pid), ]
pid$pid_cum <- cumprod(pid$pid); pid$pid_sib_cum <- cumprod(pid$pid_sib)
say("\npanel: P(ID) = %.3g, P(ID)sib = %.3g over %d loci",
    tail(pid$pid_cum, 1), tail(pid$pid_sib_cum, 1), nrow(pid))

conf <- gid_sample_confidence(res, gt = gt, pid_tab = pid,
                              post_cut = POST_CUT, min_loci = MIN_LOCI)
say("confidence: %s", paste(sprintf("%s=%d", names(table(conf$status)),
                                    table(conf$status)), collapse = "  "))

## ------------------------------------------------------- 6. field metadata --
lab <- suppressMessages(as.data.frame(read_excel(XLSX, sheet = "labdata"),
                                      check.names = FALSE))
lab$SampleID <- trimws(as.character(lab$SampleID))
lab$area <- trimws(as.character(lab$UNK))

## The collection dates come through as a mix of Excel serial numbers stored as
## text and ordinary date strings, so parse both rather than trusting one.
parse_date <- function(x) {
  x <- trimws(as.character(x))
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  ser <- !is.na(x) & grepl("^[0-9]{5}$", x)
  out[ser] <- as.Date(as.numeric(x[ser]), origin = "1899-12-30")
  ## as.Date() throws on unparseable text rather than returning NA, and this
  ## column also carries free-text notes, so give it explicit formats.
  txt <- which(!is.na(x) & !ser)
  for (k in txt) for (fmt in c("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%d/%m/%Y")) {
    v <- as.Date(x[k], format = fmt)
    if (!is.na(v)) { out[k] <- v; break }
  }
  out
}
lab$date <- parse_date(lab$CollDate)
meta <- lab[!duplicated(lab$SampleID), c("SampleID", "area", "date")]

conf$species <- species[conf$sample]
conf$sex     <- sex[conf$sample]
base_id      <- sub("\\|[0-9]+$", "", conf$sample)
conf$area    <- meta$area[match(base_id, meta$SampleID)]
conf$date    <- meta$date[match(base_id, meta$SampleID)]

## ------------------------------------------------ 7. name the individuals ---
## Number them by how often they were sampled, then by first capture, so the
## most-encountered animals are W01 onwards and the ordering means something.
ord <- do.call(rbind, lapply(split(conf, conf$individual), function(x)
  data.frame(individual = x$individual[1], species = x$species[1],
             n = nrow(x), first = suppressWarnings(min(x$date, na.rm = TRUE)),
             stringsAsFactors = FALSE)))
ord$first[!is.finite(ord$first)] <- as.Date(NA)
ord <- ord[order(ord$species, -ord$n, ord$first), ]
ord$label <- ave(ord$species, ord$species, FUN = function(z)
  sprintf("%s%02d", ifelse(z[1] == "wolf", "W", "C"), seq_along(z)))
conf$animal <- ord$label[match(conf$individual, ord$individual)]

## --------------------------------------------------------- 8. the roster ---
roster <- do.call(rbind, lapply(split(conf, conf$animal), function(x) {
  dts <- x$date[!is.na(x$date)]
  data.frame(
    animal = x$animal[1], species = x$species[1],
    n_samples = nrow(x),
    sex = paste(unique(x$sex[!is.na(x$sex)]), collapse = "/"),
    area = paste(unique(x$area[!is.na(x$area)]), collapse = "/"),
    first_seen = if (length(dts)) min(dts) else as.Date(NA),
    last_seen  = if (length(dts)) max(dts) else as.Date(NA),
    span_days  = if (length(dts) > 1) as.integer(max(dts) - min(dts)) else 0L,
    weakest_margin = min(x$margin, na.rm = TRUE),
    status = as.character(x$status[which.min(x$margin)]),
    samples = paste(x$sample, collapse = ", "),
    stringsAsFactors = FALSE)
}))
roster <- roster[order(roster$species, -roster$n_samples, roster$first_seen), ]

sample_map <- conf[order(conf$animal, conf$sample),
                   c("sample", "animal", "species", "sex", "area", "date",
                     "n_in_individual", "held_by", "n_loci", "support", "rival",
                     "rival_sample", "margin", "status")]
sample_map$rival_animal <- conf$animal[match(sample_map$rival_sample, conf$sample)]

write.csv(sample_map, file.path(OUT, "WOLF_sample_map.csv"), row.names = FALSE)
write.csv(roster,     file.path(OUT, "WOLF_individual_roster.csv"), row.names = FALSE)

## ------------------------------------------------------------- 9. checks ---
say("\n--- checks that never entered the genetics ---")
if (length(dup_key)) {
  k <- conf$animal[conf$sample %in% dup_key]
  say("positive control (%s, one sample run on two plates): same animal = %s (%s)",
      sub("\\|.*", "", dup_key[1]), length(unique(k)) == 1,
      paste(unique(k), collapse = ", "))
}
sx <- tapply(conf$sex, conf$animal, function(z) length(unique(z[!is.na(z)])) > 1)
say("individuals with conflicting sex calls: %d of %d", sum(sx, na.rm = TRUE), length(sx))
ar <- tapply(conf$area, conf$animal, function(z) length(unique(z[!is.na(z)])) > 1)
say("individuals spanning more than one study area: %d", sum(ar, na.rm = TRUE))
multi <- roster[roster$n_samples > 1, ]
say("resampled animals: %d, spanning %d to %d days",
    nrow(multi), min(multi$span_days), max(multi$span_days))
say("internally inconsistent clusters: %d", res$n_conflict)

say("\n=== RESULT ===")
say("%d samples -> %d individuals (%d wolves, %d coyotes)",
    nrow(gt), nrow(roster), sum(roster$species == "wolf"), sum(roster$species == "coyote"))
s <- gid_summarise_assignment(res$assignment)
say("samples per animal: median %g, mean %.2f, max %d",
    s$median_samples, s$mean_samples, s$max_cluster)
say("seen once: %d   resampled: %d", s$n_singletons, s$n_individuals - s$n_singletons)
say("flagged for review: %d samples", sum(conf$status != "Confident"))

saveRDS(list(res = res, sethi = sethi, conf = conf, roster = roster,
             sample_map = sample_map, cal = cal, err = err, pid = pid,
             species = species, sex = sex, gt = gt),
        file.path(OUT, "wolf_individuals.rds"))
say("\nwrote WOLF_sample_map.csv, WOLF_individual_roster.csv, wolf_cutoff_calibration.csv")
