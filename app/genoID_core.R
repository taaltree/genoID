## genoID_core.R -------------------------------------------------------------
##
## Species-agnostic core for collapsing consensus multilocus genotypes into
## unique individuals from noninvasive samples.
##
## Nothing in this file is specific to wolves, to SNPs, or to any one panel.
## It is loaded both by run_analysis.R (batch) and by app/app.R (Shiny).
##
## Genotype representation used internally
## ---------------------------------------
##   gt   : sample x locus character matrix, alleles sorted, NA = missing
##   code : sample x locus integer matrix, 0 = missing, 1..G = locus genotype id
##   For every locus we keep a dictionary of its observed genotypes, which lets
##   every pairwise quantity be precomputed as a small G x G lookup table.
##   That is what makes the probabilistic method fast enough for a web app and
##   is also what makes it work unchanged for microsatellites.
##
## Dependencies: BASE R ONLY. allelematch is optional -- that method degrades
## to a built-in equivalent if the package is absent. Keeping the core free of
## compiled packages is what lets the whole thing run in a browser under
## WebAssembly, and lets you source this file anywhere without installing
## anything.
## ---------------------------------------------------------------------------

MISSING_CODES <- c("00", "0", "000", "NA", "N/A", "", "--", "-", "?", "00/00",
                   "0/0", "000000", "NN", "..", ".")


# =============================================================================
# 1. INPUT, NORMALIZATION, QC
# =============================================================================

#' Read a genotype table from csv/tsv/xlsx
gid_read <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("readxl is required to read Excel files.")
    if (is.null(sheet)) sheet <- 1
    df <- as.data.frame(readxl::read_excel(path, sheet = sheet), check.names = FALSE)
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    df <- utils::read.csv(path, sep = sep, stringsAsFactors = FALSE,
                          check.names = FALSE, fileEncoding = "UTF-8-BOM",
                          colClasses = "character")
  }
  names(df) <- trimws(names(df))
  df
}


#' Normalise one genotype string.
#'
#' Strips quality flags (*, ?, !), upper-cases, splits on a separator if one is
#' present, sorts the two alleles so that "TC" and "CT" are the same genotype,
#' and maps every recognised missing code to NA.
gid_norm_gt <- function(x, sep = NULL, strip_flags = TRUE) {
  x <- toupper(trimws(as.character(x)))
  if (strip_flags) x <- gsub("[*?!#]", "", x)
  x[is.na(x)] <- ""
  out <- rep(NA_character_, length(x))

  has_sep <- !is.null(sep) && nzchar(sep)
  if (has_sep) {
    parts <- strsplit(x, sep, fixed = TRUE)
    a1 <- vapply(parts, function(p) if (length(p) >= 1) p[1] else "", "")
    a2 <- vapply(parts, function(p) if (length(p) >= 2) p[2] else "", "")
    too_long <- vapply(parts, length, 1L) > 2
  } else {
    # single-character alleles concatenated, e.g. "AG"; also tolerate "A".
    # Anything longer is not a genotype -- refusing to truncate it is what
    # stops PlateID and Species_check columns being read as loci.
    a1 <- substr(x, 1, 1)
    a2 <- ifelse(nchar(x) >= 2, substr(x, 2, 2), substr(x, 1, 1))
    too_long <- nchar(x) > 2
  }
  a1 <- trimws(a1); a2 <- trimws(a2)

  bad <- a1 %in% MISSING_CODES | a2 %in% MISSING_CODES | !nzchar(a1) | !nzchar(a2) |
    x %in% MISSING_CODES | (too_long & !x %in% MISSING_CODES)
  # order-insensitive genotype
  lo <- pmin(a1, a2); hi <- pmax(a1, a2)
  out <- paste0(lo, "/", hi)
  out[bad] <- NA_character_
  out
}


GID_DNA <- c("A", "C", "G", "T", "N")


#' Guess the allele separator used in a table, if any.
gid_guess_sep <- function(df, exclude = NULL) {
  x <- unlist(lapply(setdiff(names(df), exclude), function(c) as.character(df[[c]])))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NULL)
  for (s in c("/", "|", ":", "-", " ")) {
    hit <- mean(vapply(strsplit(x, s, fixed = TRUE), length, 1L) == 2)
    if (hit > 0.4) return(s)
  }
  NULL
}


#' Guess which columns are genotype loci.
#'
#' Deliberately strict. A permissive detector will happily read a plate ID, a
#' read count or a species-call column as a locus and quietly poison every
#' downstream number, so a column has to look like a diploid genotype in
#' several independent ways before it is accepted:
#'   - every value parses to exactly two alleles (or a recognised missing code)
#'   - all alleles are the same width
#'   - the alphabet is DNA, or numeric fragment sizes when a separator is used
#'   - between 2 and max_alleles distinct alleles, and few distinct genotypes
#'   - the column is not simply a continuous number
#' Whatever it returns should still be shown to the user for confirmation.
gid_detect_loci <- function(df, exclude = NULL, sep = NULL,
                            max_genotypes = 60, max_alleles = 40,
                            min_call_rate = 0.05,
                            alphabet = if (is.null(sep)) "dna" else "any") {
  cand <- setdiff(names(df), exclude)
  keep <- character(0)
  for (cn in cand) {
    v <- df[[cn]]
    if (is.numeric(v)) next
    ch <- toupper(trimws(as.character(v)))
    ch <- ch[!is.na(ch) & nzchar(ch)]
    if (!length(ch)) next

    # a column that is a continuous measurement is not a locus
    suppressWarnings(num <- as.numeric(ch))
    if (!anyNA(num) && length(unique(num)) > 12) next

    g  <- gid_norm_gt(v, sep = sep)
    if (mean(!is.na(g)) < min_call_rate) next
    ug <- unique(g[!is.na(g)])
    if (!length(ug) || length(ug) > max_genotypes) next
    if (length(ug) > max(3, 0.5 * sum(!is.na(g)))) next

    al <- unique(unlist(strsplit(ug, "/", fixed = TRUE)))
    if (length(al) < 2 || length(al) > max_alleles) next
    if (length(unique(nchar(al))) != 1) next
    ok_alpha <- switch(alphabet,
      dna = all(al %in% GID_DNA),
      any = all(al %in% GID_DNA) || all(grepl("^[0-9]+$", al)),
      TRUE)
    if (!ok_alpha) next

    keep <- c(keep, cn)
  }
  keep
}


#' Guess which column holds the sample identifier.
#'
#' Prefers a name that looks like an identifier, then uniqueness, and penalises
#' purely numeric columns -- read counts are unique per row but are not IDs.
gid_guess_id_col <- function(df, exclude = NULL) {
  cand <- setdiff(names(df), exclude)
  if (!length(cand)) return(NULL)
  score <- vapply(cand, function(cn) {
    v <- as.character(df[[cn]])
    s <- length(unique(v)) / max(1, length(v))
    suppressWarnings(num <- as.numeric(v[!is.na(v)]))
    if (length(num) && !anyNA(num)) s <- s - 1
    if (grepl("sample|specimen|individual|indiv|barcode|(^|[^a-z])id$|^id",
              cn, ignore.case = TRUE)) s <- s + 2
    s
  }, 0)
  cand[which.max(score)]
}


#' Build the normalised sample x locus genotype matrix.
gid_matrix <- function(df, id_col, loci, sep = NULL, strip_flags = TRUE) {
  ids <- as.character(df[[id_col]])
  m <- vapply(loci, function(L) gid_norm_gt(df[[L]], sep = sep, strip_flags = strip_flags),
              character(nrow(df)))
  if (is.null(dim(m))) m <- matrix(m, nrow = nrow(df))
  dimnames(m) <- list(ids, loci)
  m
}


#' Per-locus summary: call rate, alleles, MAF, observed/expected heterozygosity.
gid_locus_stats <- function(gt) {
  do.call(rbind, lapply(colnames(gt), function(L) {
    g <- gt[, L]; g <- g[!is.na(g)]
    if (!length(g)) return(data.frame(locus = L, n = 0, call_rate = 0, n_alleles = 0,
                                      maf = NA, Ho = NA, He = NA, Fis = NA))
    a <- do.call(rbind, strsplit(g, "/", fixed = TRUE))
    af <- table(c(a[, 1], a[, 2])) / (2 * length(g))
    ho <- mean(a[, 1] != a[, 2])
    he <- 1 - sum(af^2)
    data.frame(locus = L, n = length(g), call_rate = length(g) / nrow(gt),
               n_alleles = length(af), maf = min(af), Ho = ho, He = he,
               Fis = ifelse(he > 0, 1 - ho / he, NA))
  }))
}


gid_sample_stats <- function(gt) {
  data.frame(sample = rownames(gt),
             n_typed = rowSums(!is.na(gt)),
             call_rate = rowMeans(!is.na(gt)),
             row.names = NULL, stringsAsFactors = FALSE)
}


#' Drop loci / samples that are too poorly typed to be useful.
gid_filter <- function(gt, min_locus_call = 0.5, min_sample_call = 0.6,
                       drop_loci = character(0), min_maf = 0) {
  gt <- gt[, !colnames(gt) %in% drop_loci, drop = FALSE]
  ls <- gid_locus_stats(gt)
  keepL <- ls$call_rate >= min_locus_call & (is.na(ls$maf) | ls$maf >= min_maf)
  gt <- gt[, keepL, drop = FALSE]
  keepS <- rowMeans(!is.na(gt)) >= min_sample_call
  list(gt = gt[keepS, , drop = FALSE],
       dropped_loci = ls$locus[!keepL],
       dropped_samples = rownames(gt)[!keepS])
}


# =============================================================================
# 2. ALLELE FREQUENCIES, P(ID)
# =============================================================================

#' Allele frequencies per locus.
#'
#' @param weights optional per-sample weights. Pass 1/(cluster size) to avoid
#'   letting a heavily recaptured individual dominate the frequencies -- this
#'   matters because P(ID) and the LR are both driven by these numbers.
gid_allele_freq <- function(gt, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(gt))
  lapply(setNames(colnames(gt), colnames(gt)), function(L) {
    g <- gt[, L]; w <- weights[!is.na(g)]; g <- g[!is.na(g)]
    if (!length(g)) return(c(A = 1))
    a <- do.call(rbind, strsplit(g, "/", fixed = TRUE))
    tab <- tapply(c(w, w), c(a[, 1], a[, 2]), sum)
    tab <- tab[!is.na(tab)]
    tab / sum(tab)
  })
}


#' Probability of identity, unrelated and full-sib (Waits, Luikart & Taberlet 2001).
gid_pid <- function(freqs, n = NULL) {
  res <- do.call(rbind, lapply(names(freqs), function(L) {
    p <- freqs[[L]]
    pid <- sum(p^4) + sum(outer(p, p, function(a, b) (2 * a * b)^2)[upper.tri(diag(length(p)))])
    s2 <- sum(p^2); s4 <- sum(p^4)
    pid_sib <- 0.25 + 0.5 * s2 + 0.5 * s2^2 - 0.25 * s4
    data.frame(locus = L, n_alleles = length(p), pid = pid, pid_sib = pid_sib)
  }))
  res$pid_cum     <- cumprod(res$pid)
  res$pid_sib_cum <- cumprod(res$pid_sib)
  res
}


#' Composite linkage disequilibrium (r^2) between every pair of loci.
#'
#' P(ID) multiplies across loci, which assumes they are independent. Panels
#' that carry two SNPs from one amplicon violate that badly and make the panel
#' look more powerful than it is. Allele dosage correlation on unphased
#' genotypes is the standard composite estimator.
gid_ld <- function(gt, min_pairs = 20) {
  loci <- colnames(gt)
  dose <- vapply(loci, function(L) {
    g <- gt[, L]; a <- strsplit(g, "/", fixed = TRUE)
    alle <- names(sort(table(unlist(a)), decreasing = TRUE))
    if (length(alle) < 2) return(rep(NA_real_, nrow(gt)))
    ref <- alle[1]
    vapply(a, function(z) if (is.null(z) || any(is.na(z))) NA_real_ else sum(z != ref), 0)
  }, numeric(nrow(gt)))
  colnames(dose) <- loci
  ij <- which(upper.tri(matrix(0, length(loci), length(loci))), arr.ind = TRUE)
  out <- do.call(rbind, lapply(seq_len(nrow(ij)), function(k) {
    x <- dose[, ij[k, 1]]; y <- dose[, ij[k, 2]]
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < min_pairs || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NULL)
    r <- suppressWarnings(cor(x[ok], y[ok]))
    data.frame(locus1 = loci[ij[k, 1]], locus2 = loci[ij[k, 2]],
               n = sum(ok), r2 = r^2, stringsAsFactors = FALSE)
  }))
  out[order(-out$r2), ]
}


#' Cumulative P(ID)sib for an arbitrary subset of loci (used per sample).
gid_pid_subset <- function(pid_tab, loci) {
  s <- pid_tab[pid_tab$locus %in% loci, ]
  c(pid = prod(s$pid), pid_sib = prod(s$pid_sib), n_loci = nrow(s))
}


# =============================================================================
# 3. GENOTYPING ERROR RATES FROM REPLICATES
# =============================================================================

#' Estimate allelic dropout and false-allele rates from replicate genotypes.
#'
#' Uses the consensus (or, if absent, a per-sample majority call) as the
#' reference truth, following Broquet & Petit (2004):
#'   dropout      = P(observe homozygote | true heterozygote)
#'   false allele = P(observe a non-reference allele | true homozygote)
#' Both are returned per locus and pooled. The reference is derived from the
#' same replicates, so these are lower bounds -- errors shared by all
#' replicates are invisible. That is stated in the report rather than hidden.
gid_error_rates <- function(rep_gt, ref_gt, sample_id_rep, sample_id_ref) {
  stopifnot(identical(colnames(rep_gt), colnames(ref_gt)))
  idx <- match(sample_id_rep, sample_id_ref)
  ok  <- !is.na(idx)
  rep_gt <- rep_gt[ok, , drop = FALSE]; idx <- idx[ok]
  ref <- ref_gt[idx, , drop = FALSE]

  out <- do.call(rbind, lapply(colnames(rep_gt), function(L) {
    r <- rep_gt[, L]; t <- ref[, L]
    keep <- !is.na(t) & !is.na(r)
    r <- r[keep]; t <- t[keep]
    ta <- strsplit(t, "/", fixed = TRUE); ra <- strsplit(r, "/", fixed = TRUE)
    het <- vapply(ta, function(x) x[1] != x[2], TRUE)

    n_het <- sum(het)
    # dropout: true het observed as a homozygote for one of its own alleles
    drop <- if (n_het) sum(mapply(function(tt, rr) rr[1] == rr[2] && rr[1] %in% tt,
                                  ta[het], ra[het])) else 0
    n_hom <- sum(!het)
    # false allele: true hom observed carrying an allele it should not have
    fa <- if (n_hom) sum(mapply(function(tt, rr) !all(rr %in% tt),
                                ta[!het], ra[!het])) else 0
    data.frame(locus = L, n_het = n_het, n_dropout = drop,
               dropout = ifelse(n_het > 0, drop / n_het, NA),
               n_hom = n_hom, n_false = fa,
               false_allele = ifelse(n_hom > 0, fa / n_hom, NA))
  }))
  attr(out, "pooled") <- c(
    dropout      = sum(out$n_dropout) / max(1, sum(out$n_het)),
    false_allele = sum(out$n_false)   / max(1, sum(out$n_hom)))
  out
}


#' Maximum-likelihood dropout and false-allele rates from replicate genotypes.
#'
#' The Broquet & Petit estimator above scores each replicate against a
#' consensus that was itself built from those replicates. When the consensus
#' rule demands three identical replicates to call a homozygote, the false
#' allele rate it returns is forced to zero by construction. This estimator
#' avoids that: it never looks at a consensus. For each sample x locus it
#' treats the true genotype as unknown and integrates over it,
#'
#'   L = prod_{sample,locus} sum_g P(g) prod_k P(obs_k | g, d, f)
#'
#' and maximises over the two rates. Requires >= 2 replicates per sample.
gid_error_ml <- function(rep_gt, sample_ids, freqs = NULL,
                         init = c(0.05, 0.01), min_reps = 2) {
  loci <- colnames(rep_gt)
  if (is.null(freqs)) freqs <- gid_allele_freq(rep_gt)

  # collapse to one list per locus of replicate-observation vectors
  by_locus <- lapply(loci, function(L) {
    v <- split(rep_gt[, L], sample_ids)
    v <- lapply(v, function(z) z[!is.na(z)])
    v[vapply(v, length, 1L) >= min_reps]
  })
  names(by_locus) <- loci

  nll <- function(par) {
    d <- plogis(par[1]); f <- plogis(par[2])
    tot <- 0
    for (L in loci) {
      obs <- by_locus[[L]]
      if (!length(obs)) next
      p <- freqs[[L]]; dict <- sort(unique(unlist(obs)))
      if (length(dict) < 1) next
      E  <- gid_obs_matrix(dict, p, d, f)          # E[obs, true]
      gp <- gid_geno_prob(p, dict)
      idx <- lapply(obs, function(z) match(z, dict))
      lik <- vapply(idx, function(ii) sum(gp * apply(E[ii, , drop = FALSE], 2, prod)), 0)
      tot <- tot - sum(log(pmax(lik, 1e-300)))
    }
    tot
  }
  fit <- stats::optim(qlogis(init), nll, method = "Nelder-Mead",
                      control = list(reltol = 1e-9, maxit = 800))
  # profile-likelihood 95% interval on each rate (2 log-likelihood units)
  ci <- function(k) {
    grid <- qlogis(seq(1e-4, 0.5, length.out = 120))
    ll <- vapply(grid, function(v) { par <- fit$par; par[k] <- v; -nll(par) }, 0)
    ok <- plogis(grid[ll >= max(ll) - 1.92])
    c(min(ok), max(ok))
  }
  c(dropout = plogis(fit$par[1]), false_allele = plogis(fit$par[2]),
    dropout_lo = ci(1)[1], dropout_hi = ci(1)[2],
    false_lo = ci(2)[1], false_hi = ci(2)[2],
    logLik = -fit$value, converged = fit$convergence == 0)
}


#' Propagate per-replicate error into the residual error of a consensus rule.
#'
#' The consensus genotypes the lab actually analyses have already survived a
#' multi-tube rule, so their error rate is far below the per-replicate rate.
#' Monte-Carlo the rule to get the residual rates that the LR model should use.
#'
#' @param rule "taberlet" (het accepted on het_n reps, hom needs hom_n reps),
#'   "majority", or "unanimous".
gid_propagate_error <- function(dropout, false_allele, n_rep = 3,
                                rule = c("taberlet", "majority", "unanimous"),
                                hom_n = 3, het_n = 2, per_rep_missing = 0.2,
                                nsim = 20000, seed = 1) {
  rule <- match.arg(rule)
  set.seed(seed)
  simulate_one <- function(true_het) {
    obs <- character(n_rep)
    for (k in seq_len(n_rep)) {
      if (runif(1) < per_rep_missing) { obs[k] <- "00"; next }
      if (true_het) {
        obs[k] <- if (runif(1) < dropout) sample(c("A/A", "B/B"), 1) else "A/B"
      } else {
        obs[k] <- if (runif(1) < false_allele) "A/B" else "A/A"
      }
    }
    obs
  }
  call_consensus <- function(obs) {
    nz <- obs[obs != "00"]
    if (!length(nz)) return(NA_character_)
    tb <- table(nz)
    hets <- tb[names(tb) == "A/B"]
    homs <- tb[names(tb) != "A/B"]
    if (rule == "unanimous")
      return(if (length(nz) == n_rep && length(unique(nz)) == 1) nz[1] else NA_character_)
    if (rule == "majority") {
      if (length(tb) > 1 && sort(tb, decreasing = TRUE)[1] == sort(tb, decreasing = TRUE)[2])
        return(NA_character_)
      return(if (max(tb) >= 2) names(which.max(tb)) else NA_character_)
    }
    # taberlet
    if (length(hets) && hets[1] >= het_n) return("A/B")
    al <- table(unlist(strsplit(unique(nz), "/", fixed = TRUE)))
    al <- table(unlist(lapply(nz, function(z) unique(strsplit(z, "/")[[1]]))))
    if (sum(al >= het_n) == 2) return("A/B")
    if (length(homs) && max(homs) >= hom_n) return(names(homs)[which.max(homs)])
    NA_character_
  }
  het_res <- replicate(nsim, call_consensus(simulate_one(TRUE)))
  hom_res <- replicate(nsim, call_consensus(simulate_one(FALSE)))
  c(dropout      = mean(het_res %in% c("A/A", "B/B"), na.rm = TRUE),
    false_allele = mean(hom_res == "A/B", na.rm = TRUE),
    missing_het  = mean(is.na(het_res)),
    missing_hom  = mean(is.na(hom_res)))
}


# =============================================================================
# 4. INTEGER CODING + PAIRWISE MACHINERY
# =============================================================================

#' Recode the genotype matrix to integers with a per-locus dictionary.
gid_encode <- function(gt) {
  dict <- lapply(colnames(gt), function(L) sort(unique(gt[, L][!is.na(gt[, L])])))
  names(dict) <- colnames(gt)
  code <- vapply(colnames(gt), function(L) {
    m <- match(gt[, L], dict[[L]]); m[is.na(m)] <- 0L; as.integer(m)
  }, integer(nrow(gt)))
  if (is.null(dim(code))) code <- matrix(code, nrow = nrow(gt))
  dimnames(code) <- dimnames(gt)
  list(code = code, dict = dict)
}


#' All pairs, with loci compared and mismatch counts.
#'
#' n_mismatch        loci where the two genotypes differ at all
#' n_mismatch_2allele loci where they share no allele (a "hard" mismatch that
#'                    dropout alone cannot produce)
gid_pairwise <- function(enc, gt) {
  code <- enc$code; n <- nrow(code); L <- ncol(code)
  if (n < 2) return(data.frame())
  ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]

  # share-an-allele lookup per locus
  share <- lapply(colnames(gt), function(l) {
    d <- enc$dict[[l]]; G <- length(d)
    if (!G) return(matrix(TRUE, 1, 1))
    al <- strsplit(d, "/", fixed = TRUE)
    s <- outer(seq_len(G), seq_len(G), Vectorize(function(a, b) length(intersect(al[[a]], al[[b]])) > 0))
    s
  })
  names(share) <- colnames(gt)

  n_cmp <- integer(length(i)); n_mm <- integer(length(i)); n_mm2 <- integer(length(i))
  for (l in seq_len(L)) {
    a <- code[i, l]; b <- code[j, l]
    both <- a > 0 & b > 0
    n_cmp <- n_cmp + both
    diff <- both & a != b
    n_mm <- n_mm + diff
    if (any(diff)) {
      sh <- share[[l]]
      hard <- diff & !sh[cbind(pmax(a, 1L), pmax(b, 1L))]
      n_mm2 <- n_mm2 + hard
    }
  }
  data.frame(i = i, j = j,
             id1 = rownames(code)[i], id2 = rownames(code)[j],
             n_compared = n_cmp, n_mismatch = n_mm, n_mismatch_2allele = n_mm2,
             stringsAsFactors = FALSE)
}


# =============================================================================
# 5. RESOLVING A MATCH GRAPH INTO INDIVIDUALS
# =============================================================================

#' Turn a set of matching pairs into individual assignments.
#'
#' linkage = "single"    connected components: A~B, B~C implies A,B,C are one
#'                       individual even if A and C were never matched. Fast but
#'                       chains: one bad edge merges two individuals.
#'           "complete"  every member must match every other member (maximal
#'                       cliques, greedily assigned). Conservative, splits.
#' Conflicts (open triangles) are always reported -- they are the single most
#' useful diagnostic for whether the threshold is set sensibly.
#' Connected components of a symmetric logical adjacency matrix. Base R, so the
#' core has no compiled dependency and runs anywhere R does -- including inside
#' a browser under WebAssembly.
gid_components <- function(adj) {
  n <- nrow(adj)
  comp <- integer(n); k <- 0L
  for (s in seq_len(n)) {
    if (comp[s]) next
    k <- k + 1L; stack <- s
    while (length(stack)) {
      v <- stack[length(stack)]; stack <- stack[-length(stack)]
      if (comp[v]) next
      comp[v] <- k
      stack <- c(stack, which(adj[v, ] & comp == 0L))
    }
  }
  comp
}


#' Maximal cliques among the vertices `vs`, by Bron-Kerbosch with pivoting.
#' Called one connected component at a time, which keeps the search bounded --
#' components in this setting are the recapture groups, so they are small.
gid_max_cliques <- function(adj, vs) {
  res <- list()
  nbr <- lapply(vs, function(v) vs[adj[v, vs]])
  names(nbr) <- as.character(vs)
  nb <- function(v) nbr[[as.character(v)]]

  bk <- function(R, P, X) {
    if (!length(P) && !length(X)) { res[[length(res) + 1L]] <<- R; return(invisible(NULL)) }
    PX  <- c(P, X)
    piv <- PX[which.max(vapply(PX, function(u) length(intersect(P, nb(u))), 1L))]
    for (v in setdiff(P, nb(piv))) {
      bk(c(R, v), intersect(P, nb(v)), intersect(X, nb(v)))
      P <- setdiff(P, v); X <- c(X, v)
    }
  }
  bk(integer(0), vs, integer(0))
  res
}


gid_resolve <- function(pairs_matched, all_ids, linkage = c("single", "complete")) {
  linkage <- match.arg(linkage)
  n <- length(all_ids)
  adj <- matrix(FALSE, n, n)

  if (!is.null(pairs_matched) && nrow(pairs_matched)) {
    a <- match(pairs_matched$id1, all_ids)
    b <- match(pairs_matched$id2, all_ids)
    ok <- !is.na(a) & !is.na(b) & a != b
    if (any(ok)) {
      adj[cbind(a[ok], b[ok])] <- TRUE
      adj[cbind(b[ok], a[ok])] <- TRUE
    }
  }

  comp <- gid_components(adj)

  # A cluster whose members do not all match each other is being held together
  # by a chain of intermediate samples, and may be two animals.
  n_conflict <- 0L; conflicts <- list()
  for (cid in unique(comp)) {
    v <- which(comp == cid)
    if (length(v) < 3) next
    n_edge <- sum(adj[v, v]) / 2
    dens <- n_edge / choose(length(v), 2)
    if (dens < 1) {
      n_conflict <- n_conflict + 1L
      conflicts[[length(conflicts) + 1L]] <- data.frame(
        component = cid, n_samples = length(v),
        n_edges = n_edge, n_possible = choose(length(v), 2),
        completeness = dens,
        members = paste(all_ids[v], collapse = ", "),
        stringsAsFactors = FALSE)
    }
  }

  if (linkage == "single") {
    memb <- comp
  } else {
    # greedy maximal cliques, largest first, within each component
    memb <- rep(NA_integer_, n); k <- 0L
    cl <- unlist(lapply(unique(comp), function(cid) {
      v <- which(comp == cid)
      if (length(v) == 1L) list(v) else gid_max_cliques(adj, v)
    }), recursive = FALSE)
    cl <- cl[order(-vapply(cl, length, 1L))]
    for (cq in cl) {
      v <- as.integer(cq); v <- v[is.na(memb[v])]
      if (!length(v)) next
      k <- k + 1L; memb[v] <- k
    }
    lone <- which(is.na(memb))
    if (length(lone)) memb[lone] <- seq.int(k + 1L, length.out = length(lone))
  }

  ind <- sprintf("IND_%03d", as.integer(factor(memb, levels = unique(memb[order(memb)]))))
  asg <- data.frame(sample = all_ids, individual = ind, stringsAsFactors = FALSE)
  # returned as list elements, not attributes: rbind() in gid_by_group would
  # silently drop attributes and the conflict count would read as zero
  list(assignment = asg,
       conflicts = if (length(conflicts)) do.call(rbind, conflicts) else NULL,
       n_conflict = n_conflict, adjacency = adj)
}


#' Collapse replicate genotypes to one consensus call per sample, using the
#' multi-tube rule of Taberlet et al. (1996): a heterozygote is accepted once
#' two replicates show it, a homozygote only when hom_n replicates agree,
#' because a homozygote can be manufactured by dropout and a heterozygote
#' cannot.
gid_consensus_from_reps <- function(rep_gt, rep_sample, hom_n = 3, het_n = 2) {
  samples <- unique(as.character(rep_sample))
  out <- matrix(NA_character_, length(samples), ncol(rep_gt),
                dimnames = list(samples, colnames(rep_gt)))
  idx <- split(seq_along(rep_sample), factor(as.character(rep_sample), levels = samples))
  for (s in samples) for (L in colnames(rep_gt)) {
    v <- rep_gt[idx[[s]], L]; v <- v[!is.na(v)]
    if (!length(v)) next
    tb  <- table(v)
    het <- names(tb)[vapply(names(tb), function(z) {
      a <- strsplit(z, "/", fixed = TRUE)[[1]]; a[1] != a[2] }, TRUE)]
    if (length(het) && max(tb[het]) >= het_n) {
      out[s, L] <- het[which.max(tb[het])]; next
    }
    ## two alleles each seen in het_n replicates also make a heterozygote
    ac <- table(unlist(lapply(v, function(z) unique(strsplit(z, "/", fixed = TRUE)[[1]]))))
    if (sum(ac >= het_n) == 2) {
      out[s, L] <- paste(sort(names(ac)[ac >= het_n]), collapse = "/"); next
    }
    hom <- setdiff(names(tb), het)
    if (length(hom) && max(tb[hom]) >= hom_n) out[s, L] <- hom[which.max(tb[hom])]
  }
  out
}


#' Run any of the methods separately within groups (species, study area, year)
#' and glue the answers back together with globally unique individual IDs.
#'
#' Blocking on species is not optional when a panel spans two species: allele
#' frequencies differ, and a wolf must never be able to match a coyote.
gid_by_group <- function(gt, group, fun, ...) {
  group <- as.character(group)
  stopifnot(length(group) == nrow(gt))
  parts <- split(rownames(gt), group)
  offset <- 0L; out <- list(); extras <- list()
  for (g in names(parts)) {
    sub <- gt[parts[[g]], , drop = FALSE]
    if (nrow(sub) == 1L) {
      a <- data.frame(sample = rownames(sub), individual = "IND_001",
                      stringsAsFactors = FALSE)
      r <- list(method = "single", assignment = a, n_conflict = 0L)
    } else r <- fun(sub, ...)
    a <- r$assignment
    a$individual <- sprintf("%s_%s", g, a$individual)
    a$group <- g
    out[[g]] <- a
    extras[[g]] <- r
  }
  asg <- do.call(rbind, out); rownames(asg) <- NULL
  asg <- asg[match(rownames(gt), asg$sample), ]
  cf <- do.call(rbind, lapply(names(extras), function(g) {
    x <- extras[[g]]$conflicts
    if (is.null(x)) NULL else cbind(group = g, x)
  }))
  list(method = paste0(extras[[1]]$method, "_bygroup"), assignment = asg,
       by_group = extras, conflicts = cf,
       n_conflict = sum(vapply(extras, function(e) e$n_conflict %||% 0L, 0L)),
       matched_pairs = do.call(rbind, lapply(extras, function(e) e$matched_pairs)),
       pairs = do.call(rbind, lapply(extras, function(e) e$pairs)))
}


gid_summarise_assignment <- function(asg) {
  tb <- table(asg$individual)
  list(n_samples = nrow(asg), n_individuals = length(tb),
       n_singletons = sum(tb == 1), max_cluster = max(tb),
       recapture_rate = 1 - length(tb) / nrow(asg),
       cluster_sizes = as.integer(table(tb)))
}


# =============================================================================
# 6. METHOD 1 -- EXACT MATCH
# =============================================================================

#' Samples are the same individual only if every co-typed locus is identical.
#' Baseline. With any dropout at all it over-splits; that is the point of
#' including it.
gid_method_exact <- function(gt, enc = NULL, min_loci = 1,
                             linkage = "single") {
  if (is.null(enc)) enc <- gid_encode(gt)
  p <- gid_pairwise(enc, gt)
  m <- p[p$n_compared >= min_loci & p$n_mismatch == 0, , drop = FALSE]
  rs <- gid_resolve(m, rownames(gt), linkage)
  list(method = "exact", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m, pairs = p,
       settings = list(min_loci = min_loci, linkage = linkage))
}


# =============================================================================
# 7. METHOD 2 -- MISMATCH THRESHOLD (+ SWEEP)
# =============================================================================

#' Classic threshold rule: same individual if mismatches <= max_mismatch over
#' at least min_loci co-typed loci. The threshold is normally chosen from the
#' shape of the mismatch distribution (Paetkau 2003): with a good panel there
#' is a gap between the near-zero pile (recaptures + error) and the bulk
#' (different individuals). gid_threshold_sweep() produces that diagnostic.
gid_method_threshold <- function(gt, enc = NULL, max_mismatch = 1, min_loci = 20,
                                 max_mismatch_2allele = 0, linkage = "single") {
  if (is.null(enc)) enc <- gid_encode(gt)
  p <- gid_pairwise(enc, gt)
  m <- p[p$n_compared >= min_loci &
           p$n_mismatch <= max_mismatch &
           p$n_mismatch_2allele <= max_mismatch_2allele, , drop = FALSE]
  rs <- gid_resolve(m, rownames(gt), linkage)
  list(method = "threshold", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m, pairs = p,
       settings = list(max_mismatch = max_mismatch, min_loci = min_loci,
                       max_mismatch_2allele = max_mismatch_2allele,
                       linkage = linkage))
}


#' Sweep the mismatch threshold and record how the answer moves.
#' A flat plateau across several thresholds = a well-separated, trustworthy
#' answer. A monotone slide = the panel cannot separate individuals from
#' near-matches and the threshold choice is doing the work, not the data.
gid_threshold_sweep <- function(gt, enc = NULL, max_k = 6, min_loci = 20,
                                linkage = "single") {
  if (is.null(enc)) enc <- gid_encode(gt)
  p <- gid_pairwise(enc, gt)
  do.call(rbind, lapply(0:max_k, function(k) {
    m <- p[p$n_compared >= min_loci & p$n_mismatch <= k, , drop = FALSE]
    rs <- gid_resolve(m, rownames(gt), linkage)
    s <- gid_summarise_assignment(rs$assignment)
    data.frame(max_mismatch = k, n_matched_pairs = nrow(m),
               n_individuals = s$n_individuals, max_cluster = s$max_cluster,
               n_conflicting_clusters = rs$n_conflict)
  }))
}


# =============================================================================
# 8. METHOD 3 -- ALLELEMATCH (Galpern et al. 2012)
# =============================================================================

#' allelematch scores dissimilarity as mismatching loci / co-typed loci and
#' picks the threshold that minimises samples it cannot confidently classify.
#' Its real contribution is the "unclassified" category: samples that are
#' ambiguous get flagged rather than forced into a cluster.
#'
#' If the package is unavailable we fall back to the same dissimilarity metric
#' with a threshold chosen by the same criterion, so this method always runs.
gid_method_allelematch <- function(gt, alleleMismatch = NULL, cutHeight = NULL,
                                   min_loci = 20, profile = TRUE) {
  if (!requireNamespace("allelematch", quietly = TRUE))
    return(gid_method_allelematch_fallback(gt, alleleMismatch, min_loci))

  df <- data.frame(id = rownames(gt), stringsAsFactors = FALSE)
  for (L in colnames(gt)) {
    a <- do.call(rbind, strsplit(ifelse(is.na(gt[, L]), "NA/NA", gt[, L]), "/", fixed = TRUE))
    df[[paste0(L, "a")]] <- ifelse(a[, 1] == "NA", "-99", a[, 1])
    df[[paste0(L, "b")]] <- ifelse(a[, 2] == "NA", "-99", a[, 2])
  }

  amd <- try(allelematch::amDataset(df, missingCode = "-99", indexColumn = "id"),
             silent = TRUE)
  if (inherits(amd, "try-error"))
    return(gid_method_allelematch_fallback(gt, alleleMismatch, min_loci))

  prof <- NULL
  if (is.null(alleleMismatch) && profile) {
    prof <- try(allelematch::amUniqueProfile(amd, doPlot = FALSE, verbose = FALSE),
                silent = TRUE)

    if (!inherits(prof, "try-error") && "guessOptimum" %in% names(prof) &&
        any(prof$guessOptimum)) {
      # The package computes its own optimum (the second minimum of the
      # unclassified curve). Use it -- picking the largest or smallest value
      # that merely ties the minimum is not the same thing and goes badly
      # wrong at both ends of the sweep.
      alleleMismatch <- prof$alleleMismatch[which(prof$guessOptimum)][1]

    } else {
      # amUniqueProfile() errors on some SNP panels. Reproduce its criterion:
      # sweep, and take the first value where the count of samples that are
      # unclassified or match more than one genotype returns to its minimum
      # AFTER a stretch where it was higher. Both ends of the sweep are
      # degenerate -- at 0 nothing merges, at very high values everything
      # does -- and both make the ambiguity count vanish for the wrong reason.
      kmax <- min(30L, max(6L, 2L * ncol(gt)))
      prof <- do.call(rbind, lapply(0:kmax, function(k) {
        u <- try(allelematch::amUnique(amd, alleleMismatch = k,
                                       cutHeight = cutHeight, verbose = FALSE),
                 silent = TRUE)
        if (inherits(u, "try-error")) return(NULL)
        data.frame(alleleMismatch = k, unique = length(u$unique$index),
                   unclassified = u$numUnclassified,
                   multipleMatch = u$numMultipleMatches,
                   ambiguous = u$numUnclassified + u$numMultipleMatches)
      }))
      if (is.null(prof)) return(gid_method_allelematch_fallback(gt, NULL, min_loci))
      mn   <- min(prof$ambiguous)
      rose <- which(prof$ambiguous > mn)
      cand <- if (length(rose))
        prof$alleleMismatch[prof$ambiguous == mn & seq_len(nrow(prof)) > min(rose)]
      else integer(0)
      alleleMismatch <- if (length(cand)) cand[1] else
        min(prof$alleleMismatch[prof$ambiguous == mn & prof$alleleMismatch > 0])
    }
    # alleleMismatch = 0 forces matchThreshold = 1, which no sample carrying
    # missing data can reach; allelematch then splits every incomplete sample
    # into its own cluster. Never select it automatically.
    if (is.na(alleleMismatch) || alleleMismatch == 0) alleleMismatch <- 1
  }
  if (is.null(alleleMismatch)) alleleMismatch <- 2

  res <- try(allelematch::amUnique(amd, alleleMismatch = alleleMismatch,
                                   cutHeight = cutHeight, verbose = FALSE),
             silent = TRUE)
  if (inherits(res, "try-error"))
    return(gid_method_allelematch_fallback(gt, alleleMismatch, min_loci))

  tf <- tempfile(fileext = ".csv")
  invisible(utils::capture.output(allelematch::amCSV.amUnique(res, tf)))
  cs <- utils::read.csv(tf, stringsAsFactors = FALSE, colClasses = "character")
  unlink(tf)

  asgv <- setNames(sprintf("IND_%03d", as.integer(cs$uniqueGroup)), cs$matchIndex)
  asgv <- asgv[!duplicated(names(asgv))]
  miss <- setdiff(rownames(gt), names(asgv))
  if (length(miss))
    asgv[miss] <- sprintf("IND_U%02d", seq_along(miss))
  asg <- data.frame(sample = rownames(gt), individual = unname(asgv[rownames(gt)]),
                    stringsAsFactors = FALSE)

  unclass_ids <- cs$matchIndex[cs$rowType %in% c("UNCLASSIFIED", "MULTIPLE_MATCH")]
  list(method = "allelematch", assignment = asg,
       unclassified = unique(unclass_ids), engine = "allelematch",
       profile = prof, psib = suppressWarnings(as.numeric(cs$Psib)),
       n_unclassified = res$numUnclassified, n_multiple = res$numMultipleMatches,
       settings = list(alleleMismatch = alleleMismatch,
                       cutHeight = res$cutHeight,
                       matchThreshold = res$matchThreshold))
}


#' Built-in equivalent: dissimilarity = mismatching loci, threshold chosen by
#' minimising the number of samples that fall in the ambiguous zone.
gid_method_allelematch_fallback <- function(gt, alleleMismatch = NULL, min_loci = 20) {
  enc <- gid_encode(gt); p <- gid_pairwise(enc, gt)
  p <- p[p$n_compared >= min_loci, , drop = FALSE]
  if (is.null(alleleMismatch)) {
    prof <- do.call(rbind, lapply(0:6, function(k) {
      # ambiguous = pairs sitting exactly at the boundary (k and k+1)
      amb <- sum(p$n_mismatch == k + 1)
      data.frame(alleleMismatch = k, unclassified = amb)
    }))
    alleleMismatch <- prof$alleleMismatch[which.min(prof$unclassified)]
  }
  m <- p[p$n_mismatch <= alleleMismatch, , drop = FALSE]
  rs <- gid_resolve(m, rownames(gt), "single")
  amb <- p[p$n_mismatch == alleleMismatch + 1, , drop = FALSE]
  list(method = "allelematch", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m,
       unclassified = unique(c(amb$id1, amb$id2)), engine = "builtin",
       settings = list(alleleMismatch = alleleMismatch))
}


# =============================================================================
# 8b. METHOD -- GenAlEx "Multilocus Matches" (Peakall & Smouse 2006, 2012)
# =============================================================================
#
# GenAlEx's Matches routine is what a lot of labs already run, so it belongs
# here for continuity even though it is the least powerful option. It does
# three things, and this function reproduces all three:
#
#   1. compares every pair of samples and tabulates how many loci match
#   2. reports pairs matching at every co-typed locus as the same individual,
#      and pairs one or two loci short as "near matches" to inspect by hand
#   3. reports single-locus and cumulative P(ID) and P(ID)sib as the
#      justification that the panel is powerful enough
#
# Note on estimators: this uses the standard Waits, Luikart & Taberlet (2001)
# equations, which assume Hardy-Weinberg and linkage equilibrium and take the
# observed allele frequencies at face value. GenAlEx can additionally report a
# small-sample-corrected ("unbiased") P(ID); that correction is not applied
# here, and at the sample sizes typical of these studies it is negligible
# relative to the effect of population substructure, which neither version
# accounts for. The Methods page in the app states this.
gid_method_genalex <- function(gt, enc = NULL, min_loci = 1,
                               near_match_loci = 2, linkage = "single",
                               freqs = NULL) {
  if (is.null(enc))   enc   <- gid_encode(gt)
  if (is.null(freqs)) freqs <- gid_allele_freq(gt)
  p <- gid_pairwise(enc, gt)

  # GenAlEx counts MATCHING loci, not mismatching ones
  p$n_match <- p$n_compared - p$n_mismatch

  # (1) the match-distribution table GenAlEx prints
  dist <- as.data.frame(table(mismatching_loci = p$n_mismatch[p$n_compared >= min_loci]),
                        stringsAsFactors = FALSE)
  names(dist)[2] <- "n_pairs"
  dist$mismatching_loci <- as.integer(dist$mismatching_loci)

  # (2) exact matches = same individual; near matches flagged, not merged
  m    <- p[p$n_compared >= min_loci & p$n_mismatch == 0, , drop = FALSE]
  near <- p[p$n_compared >= min_loci & p$n_mismatch > 0 &
              p$n_mismatch <= near_match_loci, , drop = FALSE]

  # (3) P(ID) and P(ID)sib, whole panel and per matched pair over shared loci
  pid <- gid_pid(freqs)
  if (nrow(m)) {
    shared <- t(vapply(seq_len(nrow(m)), function(k) {
      L <- colnames(gt)[!is.na(gt[m$id1[k], ]) & !is.na(gt[m$id2[k], ])]
      v <- gid_pid_subset(pid, L); c(v[["pid"]], v[["pid_sib"]])
    }, numeric(2)))
    m$pid <- shared[, 1]; m$pid_sib <- shared[, 2]
  }
  if (nrow(near)) {
    near$loci_matching <- near$n_compared - near$n_mismatch
  }

  rs <- gid_resolve(m, rownames(gt), linkage)
  list(method = "genalex", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m, near_matches = near,
       pairs = p, match_distribution = dist, pid = pid,
       pid_total = c(pid = prod(pid$pid), pid_sib = prod(pid$pid_sib)),
       settings = list(min_loci = min_loci, near_match_loci = near_match_loci,
                       linkage = linkage))
}


# =============================================================================
# 9. METHOD 4 -- PROBABILISTIC LIKELIHOOD RATIO
# =============================================================================
#
# For every pair of samples and every co-typed locus we compare
#   H1: one individual, observed twice through an explicit error process
#   H0: two individuals drawn from the population at a stated relatedness
# and sum log10 LR over loci. The error process is
#   true het  -> het  1 - d      ; -> either hom  d/2 each
#   true hom  -> hom  1 - f - f2 ; -> het  f      ; -> other hom  f2
# with d (dropout) and f (false allele) taken from the replicate data.
#
# The H0 relatedness matters a great deal for social carnivores: the samples
# competing to be "not the same wolf" are usually packmates, not random wolves.
# k = (k0, k1, k2) IBD coefficients: unrelated (1,0,0), half sib (.5,.5,0),
# full sib (.25,.5,.25), parent-offspring (0,1,0).
# =============================================================================

GID_KINSHIP <- list(
  unrelated        = c(1,    0,   0),
  half_sib         = c(0.5,  0.5, 0),
  full_sib         = c(0.25, 0.5, 0.25),
  parent_offspring = c(0,    1,   0)
)


#' Genotype probabilities under HWE for one locus.
gid_geno_prob <- function(p, dict) {
  vapply(dict, function(g) {
    a <- strsplit(g, "/", fixed = TRUE)[[1]]
    if (a[1] == a[2]) p[[a[1]]]^2 else 2 * p[[a[1]]] * p[[a[2]]]
  }, 0)
}


GID_KEY <- function(x, y) if (x <= y) paste0(x, "/", y) else paste0(y, "/", x)


#' Observation model P(observed genotype | true genotype), G x G.
#'
#' Two-stage generative process, so it stays coherent for any number of
#' alleles (SNPs today, microsatellites tomorrow):
#'   1. dropout   -- with prob d a true heterozygote loses one allele at random
#'                   and is seen as a homozygote
#'   2. false allele -- with prob f one of the two surviving alleles is replaced
#'                   by an allele drawn from the population
#' Mass landing on genotypes never observed at this locus is renormalised away.
gid_obs_matrix <- function(dict, p, dropout, false_allele) {
  G <- length(dict)
  if (!G) return(matrix(0, 0, 0))
  al <- strsplit(dict, "/", fixed = TRUE)
  alleles <- names(p)
  gidx <- setNames(seq_len(G), dict)
  M <- matrix(0, G, G, dimnames = list(obs = dict, true = dict))

  for (tt in seq_len(G)) {
    t_al <- al[[tt]]
    stage1 <- if (t_al[1] != t_al[2])
      list(list(a = c(t_al[1], t_al[1]), pr = dropout / 2),
           list(a = c(t_al[2], t_al[2]), pr = dropout / 2),
           list(a = t_al,                pr = 1 - dropout))
    else list(list(a = t_al, pr = 1))

    for (s in stage1) {
      gk <- gidx[GID_KEY(s$a[1], s$a[2])]
      if (!is.na(gk)) M[gk, tt] <- M[gk, tt] + s$pr * (1 - false_allele)
      for (pos in 1:2) for (x in alleles) {
        nb <- s$a; nb[pos] <- x
        gk2 <- gidx[GID_KEY(nb[1], nb[2])]
        if (!is.na(gk2))
          M[gk2, tt] <- M[gk2, tt] + s$pr * false_allele * 0.5 * p[[x]]
      }
    }
  }
  M <- M + 1e-12
  sweep(M, 2, colSums(M), "/")
}


#' P(g1, g2) for two individuals at IBD coefficients k.
gid_pair_prob <- function(p, dict, k) {
  G <- length(dict)
  gp <- gid_geno_prob(p, dict)
  alleles <- names(p)
  gidx <- setNames(seq_len(G), dict)

  P0 <- outer(gp, gp)                                        # 0 alleles IBD
  P2 <- diag(gp, nrow = G)                                   # 2 alleles IBD
  P1 <- matrix(0, G, G)                                      # exactly 1 IBD
  for (s in alleles) for (a in alleles) for (b in alleles) {
    g1 <- gidx[GID_KEY(s, a)]; g2 <- gidx[GID_KEY(s, b)]
    if (is.na(g1) || is.na(g2)) next
    P1[g1, g2] <- P1[g1, g2] + p[[s]] * p[[a]] * p[[b]]
  }
  dimnames(P0) <- dimnames(P1) <- dimnames(P2) <- list(dict, dict)
  k[1] * P0 + k[2] * P1 + k[3] * P2
}


#' Per-locus log10 likelihood of an observed genotype pair under each
#' hypothesis separately: same individual, and each relatedness state.
#'
#' Kept separate (rather than pre-divided into a ratio) so that the denominator
#' can be chosen after summing across loci -- which is what Sethi et al. (2016)
#' require and what a single fixed alternative cannot do.
gid_hypothesis_tables <- function(p, dict, dropout, false_allele,
                                  relationships = c("unrelated", "full_sib",
                                                    "parent_offspring")) {
  if (!length(dict)) return(NULL)
  E  <- gid_obs_matrix(dict, p, dropout, false_allele)
  gp <- gid_geno_prob(p, dict)
  out <- list(same_individual =
                log10(pmax((E %*% diag(gp, nrow = length(gp))) %*% t(E), 1e-300)))
  for (r in relationships)
    out[[r]] <- log10(pmax(E %*% gid_pair_prob(p, dict, GID_KINSHIP[[r]]) %*% t(E),
                           1e-300))
  out
}


#' Per-locus G x G matrix of log10 LR (same individual vs H0).
gid_lr_table <- function(p, dict, dropout, false_allele, k) {
  if (!length(dict)) return(matrix(0, 0, 0))
  E  <- gid_obs_matrix(dict, p, dropout, false_allele)   # E[obs, true]
  gp <- gid_geno_prob(p, dict)
  H1 <- (E %*% diag(gp, nrow = length(gp))) %*% t(E)     # same individual
  H0 <- E %*% gid_pair_prob(p, dict, k) %*% t(E)         # two individuals
  log10(pmax(H1, 1e-300)) - log10(pmax(H0, 1e-300))
}


#' Probabilistic individual identification.
#'
#' @param dropout,false_allele residual error rates of the *consensus*
#'   genotypes (see gid_propagate_error). Scalars or named per-locus vectors.
#' @param kinship H0 relatedness. Use "full_sib" for pack- or group-living
#'   species; it is the conservative choice.
#' @param prior_same prior probability that a random pair of samples is the
#'   same individual. Defaults to a self-consistent estimate from the data.
#' @param post_cut posterior probability required to declare a match.
# =============================================================================
# 9a. USING PCR REPLICATES DIRECTLY INSTEAD OF A CONSENSUS
# =============================================================================
#
# A consensus genotype is a decision made before the analysis starts, and every
# decision throws information away. A locus where two replicates said A/G and
# one said A/A becomes "A/G" -- the evidence of a dropout is gone. A locus the
# consensus rule refused to call becomes missing -- all three observations are
# gone. Neither loss is necessary: the likelihood already integrates over the
# unknown true genotype, so it can just as easily condition on every
# observation instead of one summary of them.
#
# For sample i at locus l with replicate observations o_1 ... o_R, define the
# genotype likelihood
#
#     Lik_il(g) = prod_r  P(o_r | g)
#
# and then, exactly as before but with Lik in place of a single observation,
#
#     P(H1)_l = sum_g   P(g) Lik_il(g) Lik_jl(g)
#     P(H0)_l = sum_g1 sum_g2  P(g1,g2) Lik_il(g1) Lik_jl(g2)
#
# A sample with no observation at a locus gets Lik = 1 everywhere, which makes
# H1 and H0 collapse to the same marginal and contribute exactly zero evidence.
# Missing data needs no special case; it simply says nothing.
#
# Note the error rates to use here are the PER-REPLICATE rates, not the much
# smaller residual rates of a consensus call.

#' Per-locus genotype likelihoods from a long table of replicate genotypes.
#'
#' @param rep_gt sample-replicate x locus genotype matrix (rows are individual
#'   PCR replicates, so one sample contributes several rows).
#' @param rep_sample the sample each row belongs to.
#' @param samples the sample order to return, defaulting to first appearance.
gid_geno_lik <- function(rep_gt, rep_sample, freqs, dropout, false_allele,
                         samples = NULL, dict = NULL) {
  loci <- colnames(rep_gt)
  rep_sample <- as.character(rep_sample)
  if (is.null(samples)) samples <- unique(rep_sample)
  idx <- match(rep_sample, samples)
  keep <- !is.na(idx)
  rep_gt <- rep_gt[keep, , drop = FALSE]; idx <- idx[keep]
  n <- length(samples)

  d <- if (length(dropout) == 1) setNames(rep(dropout, length(loci)), loci) else dropout
  f <- if (length(false_allele) == 1) setNames(rep(false_allele, length(loci)), loci) else false_allele
  if (is.null(dict))
    dict <- lapply(loci, function(L) sort(unique(rep_gt[, L][!is.na(rep_gt[, L])])))
  names(dict) <- loci

  lik <- vector("list", length(loci)); names(lik) <- loci
  observed <- matrix(FALSE, n, length(loci), dimnames = list(samples, loci))

  for (L in loci) {
    dd <- dict[[L]]; G <- length(dd)
    if (!G) { lik[[L]] <- matrix(1, n, 1); next }
    E <- gid_obs_matrix(dd, freqs[[L]], d[[L]], f[[L]])   # E[obs, true]
    m <- matrix(1, n, G, dimnames = list(samples, dd))
    o <- match(rep_gt[, L], dd)
    got <- !is.na(o)
    if (any(got)) {
      # multiply each sample's row by P(observation | g) for every replicate
      for (r in which(got)) {
        s <- idx[r]
        m[s, ] <- m[s, ] * E[o[r], ]
        observed[s, L] <- TRUE
      }
    }
    lik[[L]] <- m
  }
  list(samples = samples, loci = loci, dict = dict, lik = lik,
       observed = observed, n_reps = as.integer(table(factor(idx, levels = seq_len(n)))))
}


#' Turn a consensus genotype matrix into the same structure, so both routes
#' share one code path. Each sample simply has a single "replicate".
gid_geno_lik_from_gt <- function(gt, freqs, dropout, false_allele, dict = NULL) {
  gid_geno_lik(gt, rownames(gt), freqs, dropout, false_allele,
               samples = rownames(gt), dict = dict)
}


#' Pairwise log10 likelihood ratios from genotype likelihoods.
#' Returns a list of n x n matrices, one per hypothesis, plus the count of loci
#' both samples were observed at.
gid_pair_loglik <- function(glik, freqs, hypotheses = c("same_individual", "full_sib")) {
  n <- length(glik$samples)
  out <- lapply(hypotheses, function(h) matrix(0, n, n))
  names(out) <- hypotheses
  n_cmp <- matrix(0L, n, n)

  for (L in glik$loci) {
    dd <- glik$dict[[L]]
    if (!length(dd)) next
    M  <- glik$lik[[L]]                      # n x G
    gp <- gid_geno_prob(freqs[[L]], dd)
    ob <- glik$observed[, L]
    n_cmp <- n_cmp + outer(ob, ob, "&")

    for (h in hypotheses) {
      P <- if (h == "same_individual") diag(gp, nrow = length(gp))
           else gid_pair_prob(freqs[[L]], dd, GID_KINSHIP[[h]])
      out[[h]] <- out[[h]] + log10(pmax(M %*% P %*% t(M), 1e-300))
    }
  }
  c(out, list(n_compared = n_cmp))
}


gid_method_lr <- function(gt, freqs = NULL, enc = NULL,
                          dropout = 0.02, false_allele = 0.005,
                          kinship = "full_sib", prior_same = NULL,
                          post_cut = 0.999, min_loci = 15,
                          linkage = "single", weights = NULL, reps = NULL) {

  ## ---- replicate route: condition on every PCR observation ---------------
  if (!is.null(reps)) {
    if (is.null(freqs)) freqs <- gid_allele_freq(gt, weights)
    kin <- if (is.character(kinship)) kinship else "full_sib"
    keep <- reps$sample %in% rownames(gt)
    rgt <- reps$gt[keep, colnames(gt), drop = FALSE]
    rsm <- as.character(reps$sample[keep])
    ## A sample with no replicate rows would otherwise contribute no evidence
    ## at all and be split off as a singleton. Fall back to its consensus call
    ## as a single observation, so a partially replicated dataset still works.
    orphan <- setdiff(rownames(gt), unique(rsm))
    if (length(orphan)) {
      rgt <- rbind(rgt, gt[orphan, , drop = FALSE])
      rsm <- c(rsm, orphan)
    }
    glik <- gid_geno_lik(rgt, rsm, freqs, dropout, false_allele,
                         samples = rownames(gt),
                         dict = gid_encode(gt)$dict)
    pl <- gid_pair_loglik(glik, freqs, c("same_individual", kin))
    return(gid_lr_finish(rownames(gt),
                         pl[["same_individual"]] - pl[[kin]], pl$n_compared,
                         prior_same, post_cut, min_loci, linkage,
                         settings = list(dropout = mean(dropout),
                                         false_allele = mean(false_allele),
                                         kinship = kinship, post_cut = post_cut,
                                         min_loci = min_loci, linkage = linkage,
                                         input = "replicates",
                                         reps_per_sample = mean(glik$n_reps),
                                         n_from_consensus = length(orphan))))
  }

  if (is.null(enc))   enc   <- gid_encode(gt)
  if (is.null(freqs)) freqs <- gid_allele_freq(gt, weights)
  k <- if (is.character(kinship)) GID_KINSHIP[[kinship]] else kinship
  loci <- colnames(gt)
  d <- if (length(dropout) == 1) setNames(rep(dropout, length(loci)), loci) else dropout
  f <- if (length(false_allele) == 1) setNames(rep(false_allele, length(loci)), loci) else false_allele

  tabs <- lapply(loci, function(L)
    gid_lr_table(freqs[[L]], enc$dict[[L]], d[[L]], f[[L]], k))
  names(tabs) <- loci

  code <- enc$code; n <- nrow(code)
  ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]
  lr <- numeric(length(i)); n_cmp <- integer(length(i))
  for (L in loci) {
    a <- code[i, L]; b <- code[j, L]
    both <- a > 0 & b > 0
    if (!any(both)) next
    n_cmp <- n_cmp + both
    lr[both] <- lr[both] + tabs[[L]][cbind(a[both], b[both])]
  }

  # Prior odds that a random pair of samples is a recapture.
  #
  # Estimated self-consistently rather than from a fixed LR cutoff: start at
  # one recapture per sample, count how many pairs match at that prior, feed
  # the count back, repeat until it stops moving. A hard-coded seed threshold
  # would be silently unreachable on a panel whose maximum attainable LR falls
  # below it, and the prior would collapse to its floor.
  prior_iter <- NULL
  if (is.null(prior_same)) {
    pri <- 1 / max(n, 2)
    prior_iter <- pri
    for (it in seq_len(25)) {
      po <- pri / (1 - pri)
      nm <- sum(n_cmp >= min_loci & 1 / (1 + 10^(-(log10(po) + lr))) >= post_cut)
      new <- min(0.5, max(1 / max(n, 2), (nm + 1) / length(i)))
      prior_iter <- c(prior_iter, new)
      if (abs(new - pri) < 1e-10) break
      pri <- new
    }
    prior_same <- pri
  }
  prior_odds <- prior_same / (1 - prior_same)
  log10_post_odds <- log10(prior_odds) + lr
  post <- 1 / (1 + 10^(-log10_post_odds))

  pr <- data.frame(i = i, j = j, id1 = rownames(gt)[i], id2 = rownames(gt)[j],
                   n_compared = n_cmp, log10_LR = lr,
                   posterior_same = post, stringsAsFactors = FALSE)
  m <- pr[pr$n_compared >= min_loci & pr$posterior_same >= post_cut, , drop = FALSE]
  rs <- gid_resolve(m, rownames(gt), linkage)

  # A panel carries a finite amount of evidence, so a posterior threshold can
  # be set beyond anything the data could ever supply -- and then every pair
  # fails, which looks like a result and is not one. The theoretical ceiling
  # (best genotype pair at every locus) is a loose bound; what matters in
  # practice is the strongest evidence any real pair actually reached.
  max_lr_theoretical <- sum(vapply(tabs, function(x) if (length(x)) max(x) else 0, 0))

  list(method = "lr", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m, pairs = pr,
       lr_tables = tabs,
       max_log10_LR_theoretical = max_lr_theoretical,
       max_log10_LR_observed = if (nrow(pr)) max(pr$log10_LR[pr$n_compared >= min_loci], -Inf) else NA,
       max_posterior_observed = if (nrow(pr)) max(pr$posterior_same[pr$n_compared >= min_loci], 0) else NA,
       prior_path = prior_iter,
       settings = list(dropout = mean(d), false_allele = mean(f),
                       kinship = kinship, prior_same = prior_same,
                       post_cut = post_cut, min_loci = min_loci, linkage = linkage))
}


#' Shared tail of the likelihood-ratio methods: estimate the prior, convert
#' evidence to posteriors, threshold, and resolve the match graph. Both the
#' consensus route and the replicate route end here, so a change to the
#' decision rule cannot apply to one and not the other.
#'
#' @param lrm n x n matrix of log10 likelihood ratios, or a vector over pairs.
gid_lr_finish <- function(ids, lrm, n_cmpm, prior_same, post_cut, min_loci,
                          linkage, settings = list(), extra = list()) {
  n <- length(ids)
  ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]
  lr    <- lrm[cbind(i, j)]
  n_cmp <- n_cmpm[cbind(i, j)]

  prior_iter <- NULL
  if (is.null(prior_same)) {
    pri <- 1 / max(n, 2); prior_iter <- pri
    for (it in seq_len(25)) {
      po <- pri / (1 - pri)
      nm <- sum(n_cmp >= min_loci & 1 / (1 + 10^(-(log10(po) + lr))) >= post_cut)
      new <- min(0.5, max(1 / max(n, 2), (nm + 1) / length(i)))
      prior_iter <- c(prior_iter, new)
      if (abs(new - pri) < 1e-10) break
      pri <- new
    }
    prior_same <- pri
  }
  post <- 1 / (1 + 10^(-(log10(prior_same / (1 - prior_same)) + lr)))

  pr <- data.frame(i = i, j = j, id1 = ids[i], id2 = ids[j],
                   n_compared = n_cmp, log10_LR = lr, posterior_same = post,
                   stringsAsFactors = FALSE)
  m  <- pr[pr$n_compared >= min_loci & pr$posterior_same >= post_cut, , drop = FALSE]
  rs <- gid_resolve(m, ids, linkage)

  c(list(method = "lr", assignment = rs$assignment, conflicts = rs$conflicts,
         n_conflict = rs$n_conflict, matched_pairs = m, pairs = pr,
         max_log10_LR_observed =
           if (nrow(pr)) max(pr$log10_LR[pr$n_compared >= min_loci], -Inf) else NA,
         max_posterior_observed =
           if (nrow(pr)) max(pr$posterior_same[pr$n_compared >= min_loci], 0) else NA,
         prior_path = prior_iter,
         settings = c(settings, list(prior_same = prior_same))),
    extra)
}


#' Per-sample power: how identifiable is this sample given only the loci it
#' actually has? This is the honest answer to "can I trust a singleton?"
gid_sample_power <- function(gt, pid_tab) {
  do.call(rbind, lapply(rownames(gt), function(s) {
    L <- colnames(gt)[!is.na(gt[s, ])]
    v <- gid_pid_subset(pid_tab, L)
    data.frame(sample = s, n_loci = as.integer(v[["n_loci"]]),
               pid = v[["pid"]], pid_sib = v[["pid_sib"]],
               stringsAsFactors = FALSE)
  }))
}


# =============================================================================
# 9b. METHOD -- SETHI ET AL. (2016) ERROR-TOLERANT MATCH CALLING
# =============================================================================
#
# Same likelihood machinery as gid_method_lr(), but two deliberate differences
# that come straight from the paper:
#
#   1. The denominator is the BEST non-match explanation rather than one the
#      user picked. Sethi et al. define
#
#         Lambda = L(same individual) / max{ L(unrelated), L(full sib), L(parent-offspring) }
#
#      taken over the whole multilocus genotype, not locus by locus. A pair is
#      only called a match if it beats every competing relationship, which
#      removes the awkward "which alternative do I assume?" decision and is
#      conservative in the right direction.
#
#   2. The decision rule is Lambda > 1 -- pure strength of evidence, no prior
#      and no posterior. Simpler and assumption-free, but it does not scale
#      with how many pairs you are testing, so on a large dataset it will make
#      more false matches than a posterior rule. Both are offered.
#
# Their clustering step (compare each sample to all others, merge on a match,
# repeat until memberships stop changing) is exactly connected components of
# the match graph, which is what linkage = "single" does here.
#
# Sethi, S.A. et al. (2016) R. Soc. open sci. 3:160457.

gid_method_sethi <- function(gt, freqs = NULL, enc = NULL,
                             dropout = 0.02, false_allele = 0.005,
                             relationships = c("unrelated", "full_sib",
                                               "parent_offspring"),
                             lambda_cut = 1, min_loci = 15,
                             linkage = "single", weights = NULL, reps = NULL) {
  if (is.null(enc))   enc   <- gid_encode(gt)
  if (is.null(freqs)) freqs <- gid_allele_freq(gt, weights)

  ## ---- replicate route: condition on every PCR observation ---------------
  if (!is.null(reps)) {
    keep <- reps$sample %in% rownames(gt)
    rgt <- reps$gt[keep, colnames(gt), drop = FALSE]
    rsm <- as.character(reps$sample[keep])
    orphan <- setdiff(rownames(gt), unique(rsm))
    if (length(orphan)) { rgt <- rbind(rgt, gt[orphan, , drop = FALSE])
                          rsm <- c(rsm, orphan) }
    glik <- gid_geno_lik(rgt, rsm, freqs, dropout, false_allele,
                         samples = rownames(gt), dict = enc$dict)
    pl <- gid_pair_loglik(glik, freqs, c("same_individual", relationships))
    n  <- length(rownames(gt))
    alt <- array(unlist(pl[relationships]), dim = c(n, n, length(relationships)))
    best_ll   <- apply(alt, 1:2, max)
    best_name <- matrix(relationships[apply(alt, 1:2, which.max)], n, n)
    lam <- pl[["same_individual"]] - best_ll

    ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
    i <- ij[, 1]; j <- ij[, 2]
    pr <- data.frame(i = i, j = j, id1 = rownames(gt)[i], id2 = rownames(gt)[j],
                     n_compared = pl$n_compared[cbind(i, j)],
                     log10_lambda = lam[cbind(i, j)],
                     best_alternative = best_name[cbind(i, j)],
                     stringsAsFactors = FALSE)
    m  <- pr[pr$n_compared >= min_loci & pr$log10_lambda > log10(lambda_cut), , drop = FALSE]
    rs <- gid_resolve(m, rownames(gt), linkage)
    return(list(method = "sethi", assignment = rs$assignment,
                conflicts = rs$conflicts, n_conflict = rs$n_conflict,
                matched_pairs = m, pairs = pr,
                max_log10_lambda_observed =
                  if (nrow(pr)) max(pr$log10_lambda[pr$n_compared >= min_loci], -Inf) else NA,
                alt_used = table(m$best_alternative),
                settings = list(dropout = mean(dropout),
                                false_allele = mean(false_allele),
                                relationships = relationships,
                                lambda_cut = lambda_cut, min_loci = min_loci,
                                linkage = linkage, input = "replicates",
                                reps_per_sample = mean(glik$n_reps),
                                n_from_consensus = length(orphan))))
  }

  loci <- colnames(gt)
  d <- if (length(dropout) == 1) setNames(rep(dropout, length(loci)), loci) else dropout
  f <- if (length(false_allele) == 1) setNames(rep(false_allele, length(loci)), loci) else false_allele

  tabs <- lapply(loci, function(L)
    gid_hypothesis_tables(freqs[[L]], enc$dict[[L]], d[[L]], f[[L]], relationships))
  names(tabs) <- loci

  code <- enc$code; n <- nrow(code)
  ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i <- ij[, 1]; j <- ij[, 2]
  hyp <- c("same_individual", relationships)
  ll  <- matrix(0, length(i), length(hyp), dimnames = list(NULL, hyp))
  n_cmp <- integer(length(i))

  for (L in loci) {
    tb <- tabs[[L]]
    if (is.null(tb)) next
    a <- code[i, L]; b <- code[j, L]
    both <- a > 0 & b > 0
    if (!any(both)) next
    n_cmp <- n_cmp + both
    idx <- cbind(a[both], b[both])
    for (h in hyp) ll[both, h] <- ll[both, h] + tb[[h]][idx]
  }

  best_alt_ll   <- apply(ll[, relationships, drop = FALSE], 1, max)
  best_alt_name <- relationships[max.col(ll[, relationships, drop = FALSE], "first")]
  log10_lambda  <- ll[, "same_individual"] - best_alt_ll

  pr <- data.frame(i = i, j = j, id1 = rownames(gt)[i], id2 = rownames(gt)[j],
                   n_compared = n_cmp, log10_lambda = log10_lambda,
                   best_alternative = best_alt_name,
                   stringsAsFactors = FALSE)
  pr <- cbind(pr, setNames(as.data.frame(ll), paste0("logL_", hyp)))

  m  <- pr[pr$n_compared >= min_loci & pr$log10_lambda > log10(lambda_cut), , drop = FALSE]
  rs <- gid_resolve(m, rownames(gt), linkage)

  list(method = "sethi", assignment = rs$assignment, conflicts = rs$conflicts,
       n_conflict = rs$n_conflict, matched_pairs = m, pairs = pr,
       max_log10_lambda_observed =
         if (nrow(pr)) max(pr$log10_lambda[pr$n_compared >= min_loci], -Inf) else NA,
       alt_used = table(m$best_alternative),
       settings = list(dropout = mean(d), false_allele = mean(f),
                       relationships = relationships, lambda_cut = lambda_cut,
                       min_loci = min_loci, linkage = linkage))
}


# =============================================================================
# 10. COMPARING METHODS
# =============================================================================

#' Adjusted Rand Index between two partitions of the same samples.
gid_ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  s_ij <- sum(choose(tab, 2))
  s_i <- sum(choose(rowSums(tab), 2)); s_j <- sum(choose(colSums(tab), 2))
  exp_ <- s_i * s_j / choose(n, 2); mx <- (s_i + s_j) / 2
  if (isTRUE(all.equal(mx, exp_))) return(1)
  (s_ij - exp_) / (mx - exp_)
}


gid_compare_methods <- function(results) {
  ids <- results[[1]]$assignment$sample
  parts <- lapply(results, function(r) r$assignment$individual[match(ids, r$assignment$sample)])
  summ <- do.call(rbind, lapply(names(results), function(nm) {
    s <- gid_summarise_assignment(results[[nm]]$assignment)
    data.frame(method = nm, n_samples = s$n_samples, n_individuals = s$n_individuals,
               n_singletons = s$n_singletons, max_cluster = s$max_cluster,
               recapture_rate = s$recapture_rate,
               n_conflicting_clusters = results[[nm]]$n_conflict %||% NA_integer_)
  }))
  m <- outer(seq_along(parts), seq_along(parts),
             Vectorize(function(x, y) gid_ari(parts[[x]], parts[[y]])))
  dimnames(m) <- list(names(results), names(results))
  # Individual labels are arbitrary, so comparing them across methods is
  # meaningless. Two methods agree about a sample when they put it with the
  # same set of OTHER samples.
  cosets <- lapply(parts, function(v) split(ids, v)[v])
  disagree <- vapply(seq_along(ids), function(k)
    length(unique(lapply(cosets, function(cs) sort(setdiff(cs[[k]], ids[k]))))) > 1, TRUE)

  list(summary = summ, ari = m, n_disagree = sum(disagree),
       disagree_samples = ids[disagree],
       table = data.frame(sample = ids, disputed = disagree,
                          setNames(as.data.frame(parts), names(results)),
                          stringsAsFactors = FALSE))
}


#' Merge each individual's samples into one consensus genotype and record
#' any locus where its samples disagreed.
gid_individual_consensus <- function(gt, asg) {
  inds <- unique(asg$individual)
  out <- matrix(NA_character_, length(inds), ncol(gt),
                dimnames = list(inds, colnames(gt)))
  disc <- list()
  for (ii in inds) {
    ss <- asg$sample[asg$individual == ii]
    sub <- gt[ss, , drop = FALSE]
    for (L in colnames(gt)) {
      v <- sub[, L]; v <- v[!is.na(v)]
      if (!length(v)) next
      tb <- sort(table(v), decreasing = TRUE)
      out[ii, L] <- names(tb)[1]
      if (length(tb) > 1)
        disc[[length(disc) + 1L]] <- data.frame(individual = ii, locus = L,
          calls = paste(sprintf("%s(%d)", names(tb), as.integer(tb)), collapse = " "),
          stringsAsFactors = FALSE)
    }
  }
  list(genotypes = out,
       discordance = if (length(disc)) do.call(rbind, disc) else NULL)
}
