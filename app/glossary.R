## glossary.R -----------------------------------------------------------------
##
## One definition per piece of jargon, used in three places at once:
##   - tables show the plain-English label instead of the raw column name
##   - hovering a column header shows the definition
##   - the Methods tab prints the whole thing as a glossary
##
## If you add a column anywhere in the app, add it here too. Anything missing
## falls through and displays its raw name, so nothing breaks -- it just looks
## like code instead of English.
## ---------------------------------------------------------------------------

## key = list(label shown in tables, definition)
GID_GLOSSARY <- list(

  # ---- counts of samples and animals -------------------------------------
  n_samples = list("Samples",
    "How many samples went into this analysis, after any call-rate filtering."),
  n_individuals = list("Individuals",
    "How many distinct animals those samples came from, according to this method. This is the number you are usually after."),
  n_singletons = list("Seen once",
    "Individuals represented by exactly one sample. A high proportion usually means you have not yet sampled the population thoroughly."),
  max_cluster = list("Largest cluster",
    "The number of samples assigned to the single most-resampled animal. Normal when one animal is easy to sample repeatedly, for example if it uses a trail you walk often, but worth checking against collection dates and extraction batches in case it is a contamination artifact."),
  median_samples = list("Samples per animal (median)",
    "The typical animal was sampled this many times. 1 means most animals were seen only once. Less sensitive than the mean to one heavily resampled animal."),
  mean_samples = list("Samples per animal (mean)",
    "Total samples divided by individuals. Compare with the median: if the mean is much larger, a few animals dominate your collection."),
  median_recaptures = list("Recaptures per animal (median)",
    "Samples per animal minus one, for the typical animal. 0 means the typical animal was never recaptured."),
  mean_recaptures = list("Recaptures per animal (mean)",
    "Average number of times an animal was seen again after the first time."),
  recapture_rate = list("Recapture rate",
    "The fraction of samples that are repeats of an animal already seen: 1 - individuals/samples. 0 means every sample was a different animal."),
  n_conflicting_clusters = list("Inconsistent clusters",
    "Clusters where some members do not actually match each other, and are held together only through an intermediate sample. Each one may really be two animals. This is the single most useful warning on the page - inspect these before reporting a population size."),
  n_conflicting = list("Inconsistent clusters",
    "Clusters where some members do not match each other and are held together only through an intermediate sample."),
  n_samples_for_individual = list("Samples of this animal",
    "How many samples in total were assigned to the same animal as this one."),

  # ---- the cluster-conflict table ------------------------------------------
  component = list("Cluster",
    "An internal number identifying the cluster within its group."),
  n_edges = list("Pairs that matched",
    "How many pairs inside this cluster the method actually called a match."),
  n_possible = list("Pairs possible",
    "How many pairs there are inside this cluster in total."),
  completeness = list("Completeness",
    "Matched pairs divided by possible pairs. 1.0 means every member matched every other member, which is what you want. Below 1.0 means the cluster is held together by a chain and may be two animals."),
  members = list("Samples in this cluster", "Which samples were grouped together."),

  # ---- per-sample and per-locus quality ------------------------------------
  sample = list("Sample", "Your sample identifier."),
  individual = list("Animal",
    "The individual this sample was assigned to. These labels are arbitrary - IND_001 is simply the first animal found, and the numbers mean nothing across methods."),
  group = list("Group",
    "The species or population block this sample was analysed within. Samples in different groups are never compared."),
  species = list("Species", "Species assigned from the diagnostic loci."),
  sex = list("Sex",
    "Sex from the sex marker. Held out of matching, so it works as an independent check: two samples from one animal cannot disagree about sex."),
  call_rate = list("Call rate",
    "The fraction of loci successfully genotyped. 1.0 means every locus worked; 0.8 means a fifth failed."),
  n_typed = list("Loci typed", "How many loci were successfully genotyped in this sample."),
  n_loci = list("Loci typed", "How many loci were successfully genotyped in this sample."),
  IFI = list("IFI",
    "Individual Fuzziness Index from the GT-seq pipeline: a contamination signal. High values suggest DNA from more than one animal in the sample."),

  locus = list("Locus", "The marker name, as it appeared in your file."),
  n = list("Samples typed", "How many samples were successfully genotyped at this locus."),
  n_alleles = list("Alleles",
    "How many different allele versions were seen at this locus. A locus with one allele carries no information and is dropped."),
  maf = list("Minor allele freq.",
    "How common the rarer allele is. Near 0 means almost everyone is the same and the locus barely helps tell animals apart; near 0.5 is maximally informative."),
  Ho = list("Observed heterozygosity",
    "The fraction of animals carrying two different alleles at this locus, as observed."),
  He = list("Expected heterozygosity",
    "The fraction you would expect to carry two different alleles if the population were mating at random."),
  Fis = list("Fis",
    "How far observed heterozygosity falls from expected: 1 - Ho/He. Around 0 is normal. Strongly positive means fewer heterozygotes than expected, the classic signature of allelic dropout or a null allele at that marker. Strongly negative means more than expected, which can indicate two markers being mistaken for one."),

  # ---- linkage --------------------------------------------------------------
  locus1 = list("Locus 1", "First marker in the pair."),
  locus2 = list("Locus 2", "Second marker in the pair."),
  r2 = list("r-squared",
    "How strongly the two markers are correlated, from 0 to 1. Probability of identity assumes markers are independent, so 0 is the assumption. Values near 1 mean the second marker adds almost nothing - drop one of the pair before quoting panel power."),

  # ---- probability of identity ---------------------------------------------
  pid = list("P(ID)",
    "Probability of identity: the chance two UNRELATED animals share this genotype by coincidence. Smaller is better."),
  pid_sib = list("P(ID) siblings",
    "The same probability for two full siblings, who share far more of their genome. Always larger, and the number to quote for social or family-structured species. Below 0.01 is conventionally adequate; below 0.001 is comfortable."),
  pid_cum = list("P(ID) cumulative",
    "P(ID) multiplied across all loci up to and including this one."),
  pid_sib_cum = list("P(ID) siblings, cumulative",
    "P(ID) for siblings multiplied across all loci up to and including this one."),

  # ---- pairwise comparisons -------------------------------------------------
  id1 = list("Sample 1", "First sample in the pair."),
  id2 = list("Sample 2", "Second sample in the pair."),
  n_compared = list("Loci compared",
    "How many loci were successfully genotyped in BOTH samples. Only these can be compared, so this varies from pair to pair."),
  n_mismatch = list("Loci that differ",
    "How many of the compared loci have different genotypes in the two samples."),
  n_mismatch_2allele = list("Hard mismatches",
    "Loci where the two genotypes share no allele at all. Dropout cannot cause these, because dropout only ever removes an allele the animal already has - so even one is strong evidence for two different animals."),
  mismatching_loci = list("Loci that differ",
    "The number of loci at which a pair of samples differs."),
  n_pairs = list("Pairs", "How many sample pairs fall in this category."),
  loci_matching = list("Loci that agree", "How many compared loci had identical genotypes."),

  log10_LR = list("log10 LR",
    "The strength of evidence, on a log scale: how much better 'one animal sampled twice' explains this pair than 'two different animals'. 0 means the two explanations fit equally well. 4 means one animal fits ten thousand times better. Negative favours two animals."),
  posterior_same = list("P(same animal)",
    "The probability this pair is the same animal, combining the evidence above with a prior for how often recaptures occur in your dataset. This is what gets compared against your acceptance cutoff."),
  log10_lambda = list("log10 Lambda",
    "As log10 LR, but measured against whichever relationship best explains the pair rather than one you nominated. Above 0 means one animal is the better explanation."),
  best_alternative = list("Closest alternative",
    "Which relationship came nearest to explaining this pair without invoking one animal. If full siblings dominates, that is the null your other analyses should use."),
  logL_same_individual = list("log-likelihood: one animal", "Log-likelihood of the pair under one animal."),
  logL_unrelated = list("log-likelihood: unrelated", "Log-likelihood of the pair under two unrelated animals."),
  logL_full_sib = list("log-likelihood: full sibs", "Log-likelihood of the pair under two full siblings."),
  logL_parent_offspring = list("log-likelihood: parent-offspring", "Log-likelihood of the pair under a parent and its offspring."),

  # ---- threshold calibration -----------------------------------------------
  cutoff = list("Posterior cutoff",
    "The probability required before a pair is accepted as the same animal."),
  n_accepted = list("Pairs accepted",
    "How many sample pairs clear the cutoff and are treated as the same animal."),
  exp_false_merges = list("Expected false merges",
    "How many accepted pairs are expected to be wrong, adding up (1 - probability) across every pair you accepted. A false merge joins two animals into one and biases abundance downward."),
  exp_missed_pairs = list("Expected missed pairs",
    "How many rejected pairs are expected to be genuine recaptures, adding up the probability across every pair you rejected. A miss splits one animal into two and biases abundance upward."),
  fdr = list("False discovery rate",
    "Expected false merges divided by pairs accepted: the fraction of your matches that are expected to be wrong."),
  max_post = list("Best posterior available",
    "The highest posterior any pair in this group reaches. Set a cutoff above it and nothing can match."),

  # ---- settings and comparison ---------------------------------------------
  method = list("Method", "Which method produced this row."),
  disputed = list("Methods disagree",
    "TRUE when the methods do not all place this sample with the same set of other samples. These are the samples worth a second look."),
  kinship = list("Alternative hypothesis",
    "What 'a different animal' was taken to mean when scoring the evidence."),
  dropout = list("Assumed dropout rate",
    "The allelic dropout rate the error model was given."),
  post_cut = list("Posterior cutoff",
    "The probability required before a pair is accepted as the same animal."),
  setting = list("", "Marks the setting currently selected in the sidebar."),
  alleleMismatch = list("Alleles allowed to differ",
    "allelematch's m-hat: how many mismatching alleles it tolerates before refusing to call a match."),
  unclassified = list("Unclassified",
    "Samples allelematch could not confidently assign to any individual."),
  multipleMatch = list("Ambiguous",
    "Samples that matched more than one individual, so allelematch declined to choose."),
  ambiguous = list("Unresolved",
    "Unclassified plus ambiguous samples: the total allelematch could not settle."),
  unique = list("Unique genotypes", "How many distinct genotypes allelematch identified."),

  # ---- method names, used as column headers in the comparison table --------
  exact = list("Exact match", "Same animal only if identical at every compared locus."),
  threshold = list("Mismatch threshold", "Same animal if no more than k loci differ."),
  genalex = list("GenAlEx", "GenAlEx's Multilocus Matches routine."),
  allelematch = list("allelematch", "Galpern et al. (2012), with its own automatic threshold."),
  sethi = list("Sethi et al.", "Likelihood against the best competing relationship, accepted when Lambda > 1."),
  probabilistic = list("Likelihood ratio", "Likelihood with an explicit error model, converted to a posterior probability."),
  LR_fullsib = list("Likelihood ratio (full sib)", "Likelihood ratio tested against full siblings."),
  LR_unrelated = list("Likelihood ratio (unrelated)", "Likelihood ratio tested against unrelated animals.")
)


## internal bookkeeping columns nobody needs to see
GID_HIDE <- c("i", "j")


#' Swap raw column names for plain-English labels, remembering the definitions
#' so the table can show them on hover. Unknown columns pass through unchanged.
gid_relabel <- function(df) {
  df <- df[, setdiff(names(df), GID_HIDE), drop = FALSE]
  df <- gid_relabel_values(df)
  raw <- names(df)
  lab <- vapply(raw, function(k) {
    g <- GID_GLOSSARY[[k]]
    if (is.null(g) || !nzchar(g[[1]])) k else g[[1]]
  }, "")
  tips <- vapply(raw, function(k) {
    g <- GID_GLOSSARY[[k]]
    if (is.null(g)) "" else sprintf("%s  (column name in downloads: %s)", g[[2]], k)
  }, "")
  names(df) <- make.unique(unname(lab))
  attr(df, "gid_tips") <- unname(tips)
  df
}


## Values that are internal codes rather than data, and the words for them.
GID_VALUE_LABELS <- c(
  exact = "Exact match", threshold = "Mismatch threshold", genalex = "GenAlEx",
  allelematch = "allelematch", sethi = "Sethi et al.", probabilistic = "Likelihood ratio",
  LR_fullsib = "Likelihood ratio (full sib)", LR_unrelated = "Likelihood ratio (unrelated)",
  unrelated = "Unrelated", full_sib = "Full siblings", half_sib = "Half siblings",
  parent_offspring = "Parent-offspring", single = "Single linkage",
  complete = "Complete linkage")

## columns whose contents are those codes
GID_CODED_COLS <- c("method", "kinship", "best_alternative", "linkage")


#' Translate internal codes appearing in table cells into the same words used
#' everywhere else, so "probabilistic" does not appear next to "Likelihood
#' ratio" for the same thing.
gid_relabel_values <- function(df) {
  for (cn in intersect(names(df), GID_CODED_COLS)) {
    v <- as.character(df[[cn]])
    hit <- !is.na(v) & v %in% names(GID_VALUE_LABELS)
    v[hit] <- unname(GID_VALUE_LABELS[v[hit]])
    df[[cn]] <- v
  }
  df
}


#' The whole glossary, rendered for the Methods tab.
gid_glossary_ui <- function() {
  keys <- names(GID_GLOSSARY)
  keys <- keys[vapply(keys, function(k) nzchar(GID_GLOSSARY[[k]][[1]]), TRUE)]
  tags$table(class = "m-tbl",
    tags$thead(tags$tr(tags$th("Shown as"), tags$th("In downloads"), tags$th("Means"))),
    tags$tbody(lapply(sort(keys), function(k) {
      g <- GID_GLOSSARY[[k]]
      tags$tr(tags$td(tags$b(g[[1]])), tags$td(tags$code(k)), tags$td(g[[2]]))
    })))
}
