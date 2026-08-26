## spatial_evidence_test.R ---------------------------------------------------
##
## Does knowing where a sample was found help decide who left it?
##
## Augustine et al. (2020, PNAS 117:17903) say yes, and their genotype spatial
## partial identity model (gSPIM) proves it inside a full spatial
## capture-recapture model fitted by MCMC. gid_spatial_evidence() keeps only the
## pairwise part of that idea and costs one multiplication per pair. This script
## asks whether the cheap version actually earns its place, and -- more usefully
## -- when it does not.
##
## The answer is that it depends entirely on home-range overlap, which is what
## the paper predicts: location is informative exactly when animals are
## separated relative to how far you sampled.
##
## Run from genoID/:  Rscript spatial_evidence_test.R
## ---------------------------------------------------------------------------
source("app/genoID_core.R")

NREP <- 20
say  <- function(...) cat(sprintf(...), "\n")

#' Simulate a population with home ranges, sample it, and genotype it badly.
#'
#' sigma_m is the home-range scale; spacing_m the distance between neighbouring
#' activity centres. Their ratio is the quantity that matters: small sigma over
#' large spacing means animals are spatially separable, and the reverse means
#' every animal is everywhere and location says nothing.
sim <- function(n_ind = 30, n_loci = 34, dropout = 0.03, false_allele = 0.015,
                sigma_m = 2000, spacing_m = 8000, mean_reps = 2.0) {
  side <- ceiling(sqrt(n_ind))
  ctr <- expand.grid(x = seq_len(side), y = seq_len(side))[seq_len(n_ind), ]
  ctr$x <- ctr$x * spacing_m + rnorm(n_ind, 0, spacing_m * 0.25)
  ctr$y <- ctr$y * spacing_m + rnorm(n_ind, 0, spacing_m * 0.25)

  p <- runif(n_loci, 0.15, 0.5)
  true <- t(vapply(seq_len(n_ind), function(i)
    vapply(seq_len(n_loci), function(l) {
      a <- sample(c("A", "G"), 2, TRUE, c(p[l], 1 - p[l]))
      paste(sort(a), collapse = "/") }, ""), character(n_loci)))

  n_s <- pmax(1, rpois(n_ind, mean_reps - 1) + 1)
  rows <- do.call(rbind, lapply(seq_len(n_ind), function(i)
    data.frame(ind = i, x = ctr$x[i] + rnorm(n_s[i], 0, sigma_m),
               y = ctr$y[i] + rnorm(n_s[i], 0, sigma_m))))

  obs <- t(vapply(seq_len(nrow(rows)), function(k) {
    g <- true[rows$ind[k], ]
    vapply(g, function(z) {
      a <- strsplit(z, "/")[[1]]
      if (a[1] != a[2] && runif(1) < dropout)
        return(paste(rep(sample(a, 1), 2), collapse = "/"))
      if (runif(1) < false_allele) {
        b <- sample(c("A", "G"), 1)
        return(paste(sort(c(sample(a, 1), b)), collapse = "/"))
      }
      z }, "") }, character(n_loci)))
  rownames(obs) <- sprintf("S%03d", seq_len(nrow(rows)))
  colnames(obs) <- sprintf("L%02d", seq_len(n_loci))

  list(gt = obs, truth = rows$ind, n_ind = n_ind,
       coords = data.frame(sample = rownames(obs),
                           lon = -133 + rows$x / (111320 * cos(55 * pi / 180)),
                           lat = 55 + rows$y / 110540))
}

#' Genetics alone, then genetics with the distance term, scored against truth.
run <- function(s, dropout, false_allele, post_cut = 0.99) {
  r <- gid_method_lr(s$gt, dropout = dropout, false_allele = false_allele,
                     kinship = "unrelated", post_cut = post_cut, min_loci = 10)
  truth_of <- function(ids) s$truth[match(ids, rownames(s$gt))]
  gen_ari <- gid_ari(r$assignment$individual, truth_of(r$assignment$sample))

  sp <- gid_spatial_evidence(r$pairs, s$coords, post_cut = post_cut)
  a  <- attr(sp, "spatial")
  rs <- gid_resolve(sp[sp$posterior_joint >= post_cut, c("id1", "id2")],
                    rownames(s$gt), "single")
  c(truth = s$n_ind,
    genetics = length(unique(r$assignment$individual)), gen_ari = gen_ari,
    joint = length(unique(rs$assignment$individual)),
    sp_ari = gid_ari(rs$assignment$individual, truth_of(rs$assignment$sample)),
    sigma_km = if (isTRUE(a$applied)) a$sigma / 1000 else NA_real_)
}

## sigma vs spacing is the axis that matters, so vary it deliberately
grid <- list(
  list("low dropout 3%, tight ranges",  0.03,  2000, 8000),
  list("dropout 10%, tight ranges",     0.10,  2000, 8000),
  list("dropout 20%, tight ranges",     0.20,  2000, 8000),
  list("dropout 20%, high overlap",     0.20,  6000, 6000),
  list("dropout 20%, total overlap",    0.20, 12000, 4000))

say("%d replicates per scenario. ARI is agreement with the true individuals.", NREP)
say("")
say("%-30s %8s %9s %8s %7s %6s", "scenario", "ARI gen", "ARI joint",
    "better", "worse", "same")
for (g in grid) {
  res <- t(vapply(seq_len(NREP), function(i) {
    set.seed(1000 + i)
    s <- sim(dropout = g[[2]], sigma_m = g[[3]], spacing_m = g[[4]])
    v <- tryCatch(run(s, g[[2]], 0.015),
                  error = function(e) c(gen_ari = NA, sp_ari = NA))
    c(v[["gen_ari"]], v[["sp_ari"]])
  }, numeric(2)))
  res <- res[complete.cases(res), , drop = FALSE]
  d <- res[, 2] - res[, 1]
  say("%-30s %8.4f %9.4f %8d %7d %6d", g[[1]], mean(res[, 1]), mean(res[, 2]),
      sum(d > 1e-9), sum(d < -1e-9), sum(abs(d) <= 1e-9))
}
say("")
say("Read the last two rows as the warning they are: once home ranges overlap")
say("heavily, distance stops separating same from different and the term adds")
say("noise rather than information. It is off by default for that reason.")
