# genoID

**Turn consensus multilocus genotypes into unique individuals — in your browser.**

Upload a genotype table, get a list of animals. Six methods run side by side so
you can see whether your answer is a property of the data or of the settings you
chose. Species-agnostic: SNPs or microsatellites, any panel, any organism.

### → [Open the app](https://taaltree.github.io/genoID/)

The app is compiled to WebAssembly and runs **entirely inside your browser**.
There is no server. Uploaded genotypes never leave your machine, which matters
if your data are unpublished or covered by a data-sharing agreement.

---

## What it does

| | |
|---|---|
| **Data & QC** | Call rates, minor allele frequencies, F<sub>IS</sub>, linkage (r²) between loci, and a checklist of the things that silently corrupt identity analyses — reversed allele order, duplicated sample IDs, quality-flagged cells, monomorphic loci |
| **Panel power** | P(ID) and P(ID)<sub>sib</sub> per locus and cumulative, plus a per-sample version computed from the loci each sample actually has |
| **Individuals** | Pick which method to identify with. Explains what it does, what it is good for and what to watch, links to its full derivation, then shows the answer, the evidence for every pair, and any cluster whose members do not all match each other |
| **Sample map** | Which animal each sample belongs to, with a confidence measure per sample: how firmly it is held in its own cluster, what it nearly matched instead, and whether the assignment would survive a different cutoff. Downloadable as your capture history |
| **Scat map** | Every sample plotted where it was collected, coloured by which animal it belongs to. Switch the colouring between methods to see which samples change hands, click a scat for its ID, its individual and its full genotype, and link an animal's samples with a spider plot or a convex hull. Filter by individual, sex or year. Save it as vector PDF, raster JPG, or one self-contained HTML file a collaborator can open and manipulate |
| **Method comparison** | All five methods on the same data at once, adjusted Rand index between them, the samples they disagree about, and a sensitivity grid across 18 combinations of settings |
| **Methods** | Notation, a worked example carried through all five methods, the statistical model behind each one, and a glossary defining every column the app can show you |

Everything downloads as CSV, plus a zip of the whole run, the settings you used,
and an R script that reproduces it.

**Sharing a map.** *Shareable interactive map* writes the map you are looking
at to a single HTML file — points, links, colors and every popup embedded.
Send that file to a collaborator and they pan, zoom and click it in any
browser, with no genoID and no R. Put it on any web host or shared drive and it
becomes a link. Basemap tiles still come from the internet, so offline it draws
your samples on a blank background. The same file, plus every table as CSV,
is inside **Everything (zip)**.

> Exported maps carry your real sample IDs and coordinates. Do not commit them
> to a public repository — this one included.

**Save your work.** *Save project* writes one JSON file holding your data and
every setting you chose; *Open project* puts them all back. The data travels
inside the file rather than as a path, so a project still opens after the file
moves, is renamed, or lands on a colleague's laptop. Nothing is uploaded
anywhere — the file is written and read in your browser.

## The five methods

| Method | What it does | Use it for |
|---|---|---|
| **Exact match** | Identical at every co-typed locus | Baseline sanity check |
| **Mismatch threshold** | Allow up to *k* differing loci; sweeps *k* to show how sensitive the answer is | The field standard (Paetkau 2003) |
| **GenAlEx Matches** | Reproduces GenAlEx's Multilocus → Matches: match distribution, near-match list, P(ID) / P(ID)<sub>sib</sub> | Continuity with what your lab already runs |
| **Sethi et al. (2016)** | Explicit dropout / false-allele error model, divided by whichever relationship best explains the pair (unrelated, full sib, parent–offspring) — so you don't pick one — and decided on evidence alone (Λ > 1) | **The default** |
| **Likelihood ratio** | The same likelihood, tested against one alternative you name (unrelated, half sib, full sib, parent–offspring), returning a posterior probability per pair | **Recommended** |

### Optional: where the sample was found

Two scats a few hundred meters apart are more likely to be one animal than two
scats twenty kilometers apart, because animals have home ranges. Augustine et
al. ([2020, *PNAS*](https://www.pnas.org/doi/10.1073/pnas.2000247117)) turned
that into the **genotype spatial partial identity model (gSPIM)**, which
resolves identity, genotyping error and density together inside a spatial
capture-recapture model fitted by MCMC — and recovered 157 fisher samples that
had been discarded as unidentifiable.

genoID does not implement gSPIM: that needs a detector array and an SCR design,
which opportunistic scat collection does not have. What it implements is the
pairwise part of the same idea, at the cost of one multiplication per pair:

<p align="center">posterior odds = prior odds &times; LR(genetics) &times; LR(space)</p>

where LR(space) = *f*<sub>same</sub>(*d*) / *f*<sub>diff</sub>(*d*). The
numerator is a Rayleigh distribution — two draws from one home range of scale
&sigma; — and the denominator is read straight off the pairs the genetics has
already ruled out, so no assumption is made about the shape of the sampled area.
&sigma; comes from the pairs the genetics has already matched. Both are
estimated from the confident tails and applied only to the uncertain middle.

**It is off by default, and the reason is in the numbers.** Twenty simulated
replicates per row, scored against known truth:

| Scenario | ARI, genetics | ARI, joint | Better | Worse | Same |
|---|---|---|---|---|---|
| 3% dropout, tight ranges | 0.999 | 1.000 | 1 | 0 | 19 |
| 10% dropout, tight ranges | 0.999 | 1.000 | 1 | 0 | 19 |
| **20% dropout, tight ranges** | 0.979 | **0.995** | **9** | 1 | 10 |
| 20% dropout, high overlap | 0.979 | 0.980 | 3 | 3 | 14 |
| 20% dropout, total overlap | 0.979 | 0.977 | 0 | 2 | 18 |

Location earns its keep in exactly one situation: **when dropout is bad enough
to make the genetics ambiguous *and* home ranges are small next to the area you
sampled.** Where the genetics is already clean it changes nothing, and where
home ranges overlap heavily it adds noise rather than information — the last row
is a small net loss. Reproduce it with `Rscript spatial_evidence_test.R`.

The Individuals tab reports the estimated home-range scale, how many pairs it
was calibrated on, and how far apart matched and unrelated pairs actually sit,
so you can see whether your data is in the regime where this helps before you
believe it.

The last two are the ones to report. Both attach a number to every pair instead
of a yes/no, and both make the alternative hypothesis explicit — which, for
pack-, pride-, or colony-living species, matters more than any threshold you can
tune.

They differ in one place: **the decision rule**. Sethi et al. accept a match
whenever it is more likely than the best competing relationship (Λ > 1), which
needs no prior but takes no account of how many pairs you tested. The likelihood
ratio method converts evidence to a posterior using a prior estimated from the
data, which is stricter on large datasets. Run both and report the range they
bracket — on the AITRC wolf dataset that is 59 to 63 animals.

**Sethi is what the app shows unless you choose otherwise.** It asks the least
of a newcomer — no alternative to name, no prior to defend — and in the
simulations below it never split recaptures off more often than the likelihood
ratio did, with neither ever joining two animals. The likelihood ratio always
runs alongside it, and the Method comparison tab explains every pair the two
decide differently.

## Try it without your own data

Click **"or load the example dataset"** in the sidebar. It's 55 simulated
samples from 28 known individuals across two species, each amplified three
times at 1.5% dropout and 0.2% false alleles per reaction, followed by the lab's
consensus row. The true answer is in the `TrueIndividual` column, so you can
score any method against it.

Press **Identify individuals** without touching anything. The app measures the
error rates from the replicates (dropout comes back near 1.4%), gives the one
noisy sample its own error scale, and finds the 28 animals by every method and
by both routes: three reactions per sample wash most errors out. The Data & QC
tab shows how alike the samples' replicates are, and the Method comparison tab
confirms the two routes agree on every sample.

## PCR replicates

If your file has several rows per sample — one per PCR replicate — the app
detects it and **runs both routes**: every replicate reaction used directly, and
the reactions collapsed to a consensus first. Both are reported side by side;
you choose which one the results tabs show. The two likelihood methods extend
naturally: for sample *i* at locus *ℓ* with replicate observations
*o*₁…*o*<sub>R</sub>, the genotype likelihood is

```
L_i(g) = prod_r P(o_r | g)
```

and it drops straight into the same two hypotheses. Missing observations
contribute a factor of 1, which makes both hypotheses collapse to the same
marginal — so missing data says nothing rather than needing a special rule.

**Every sample gets its own error rates.** One pair of rates describes the
typical sample, and scat is not typical: a minority of old or wet samples err
many times as often as the rest. Held to the typical rate, the model reads their
disagreements with themselves as differences between animals and splits their
recaptures off as new individuals. So the replicate route gives each sample one
multiplier on both rates, measured from its own replicated loci, never below 1,
and left at 1 unless its reactions reject that. The Data & QC tab shows how much
your samples differ before anything is run.

**Does it help?** With uniform sample quality, using every replicate helps
where the panel is marginal: across 180 simulations (20–40 loci, 3–25% dropout)
it landed closer to the truth in 114 cases, the consensus route in 3, and they
tied in 63 (`Rscript compare_reps_vs_consensus.R`). With mixed quality, one
shared rate is the weak point. Share of true recapture pairs each route failed
to join, eight simulated datasets per scenario:

| Route | Uniform quality | 30% noisy samples | 15% very noisy | 10% false-allele prone |
|---|---|---|---|---|
| Consensus first, Sethi | 0.0% | 1.3% | 0.0% | 0.0% |
| Consensus first, likelihood ratio | 0.0% | 2.6% | 3.0% | 1.4% |
| Every replicate, one rate, Sethi | 0.0% | 0.0% | 1.9% | 0.0% |
| Every replicate, one rate, likelihood ratio | 0.0% | 0.8% | 7.6% | 0.7% |
| **Every replicate, own rates, Sethi** | 0.0% | 0.0% | 0.0% | 0.0% |
| **Every replicate, own rates, likelihood ratio** | 0.0% | 0.3% | 0.0% | 0.0% |

No route joined two different animals in any dataset. Reproduce it with
`Rscript compare_sample_quality.R`.

**Why both routes are reported.** They assume different things about error.
The consensus route treats every consensus call as equally reliable and refuses
to call a locus whose reactions disagree; the replicate route assumes a
sample's reactions err independently at that sample's rate, which an error
repeated in all of one extract's reactions breaks. Where they agree, the count
does not depend on the choice. Where they disagree, the Method comparison tab
lists the samples — those are the ones to check by hand before you report a
number.

> **Watch the error rates.** The replicate route wants *per-reaction* dropout
> and false-allele rates (a few percent). The consensus route wants the much
> smaller residual rates that survive the multi-tube rule. The app keeps both:
> leave the sliders alone and it measures the per-reaction rates from your
> replicates when you run the analysis, then works out the consensus rates from
> them. The sliders show the rate for the route on screen.

Only the two likelihood methods can use replicates. Exact matching, the mismatch
threshold and GenAlEx compare one genotype per sample, so they always run on the
consensus and the comparison tab stays meaningful.

## Your file format

**One column per locus.** One row per sample, or one row per PCR replicate.
CSV, TSV or Excel.

### Format A — one row per sample

Use this when replicates are already collapsed to a consensus.

| SampleID | Species | Sex | LOC01 | LOC02 | LOC03 |
|---|---|---|---|---|---|
| WFS_001 | wolf | XX | AG | CC | 00 |
| WFS_002 | wolf | XY | AA | CT | TT |
| WFS_003 | coyote | XX | AG | CC | TT |

The only column the app truly needs is the sample identifier. `Species` and
`Sex` are ordinary extra columns — name them whatever you like, or leave them
out.

### Format B — one row per PCR replicate

Repeat the sample ID once per reaction. The app notices and offers to use every
observation directly, which beats collapsing them first.

| SampleID | Rep | Species | LOC01 | LOC02 | LOC03 |
|---|---|---|---|---|---|
| WFS_001 | a | wolf | AG | CC | 00 |
| WFS_001 | b | wolf | AA | CC | CT |
| WFS_001 | c | wolf | AG | 00 | 00 |
| WFS_001 | consensus | wolf | AG | CC | 00 |
| WFS_002 | a | wolf | AA | CT | TT |

The replicate label column can be called anything; you pick it in the sidebar.
If your file also carries pre-computed `consensus` rows, list that label under
**Row labels to exclude** so it is not counted as a fourth reaction — the app
finds and pre-selects labels like `consensus` for you.

A Format B file can be analyzed **either way**, and you can switch between them
in the sidebar to see whether it changes anything.

### Writing a genotype

| Thing | Write it like this | Notes |
|---|---|---|
| SNP genotype | `AG`, `CC`, `T/C` | Two alleles in one cell; separator optional |
| Microsatellite | `120/124`, `120\|124` | Separator **required**, so the two sizes can be told apart |
| Missing | `00`, `NA`, `-`, `?`, blank | All recognized. Missing loci are skipped, not counted against a pair |
| Allele order | `AG` = `GA` | Sorted before comparison, so reversed cells never become false differences |
| Quality flags | `AG*`, `CT?` | Flag stripped, genotype kept. If a flag means "do not trust this", set those cells to missing first |

> **One column per locus, not two.** If your file has `LOC01a` and `LOC01b` as
> separate columns — the layout GenAlEx uses internally — join each pair into a
> single column first.

### Sex and date, for naming and filtering

Two more optional columns. With a sex column present, animals are named
**F1, F2, M1, M2** in order of how often each was detected, and the map can be
filtered by sex. With a date column, it can be filtered by year.

| SampleID | Sex | CollDate | Latitude | Longitude |
|---|---|---|---|---|
| WFS_001 | XX | 2024-11-19 | 55.62771 | -132.87711 |
| WFS_002 | XY | 2025-07-08 | 55.48210 | -132.81940 |

Sex is read from `XX`/`XY`, `X/X`/`X/Y`, `M`/`F` or `Male`/`Female`. Numeric
codings such as 1/2 are deliberately **not** guessed — 1 = female and 1 = male
are both in common use — so recode those to letters first. Dates can be ISO,
`MM/DD/YYYY`, `DD/MM/YYYY`, an Excel serial number, or a bare year.

Both columns are found automatically by what their values look like, and you can
override the choice in the sidebar.

### Coordinates, for the Scat map tab

Add two more columns and every sample appears on a map. They are optional —
every other tab works without them.

| SampleID | Latitude | Longitude | LOC01 | LOC02 |
|---|---|---|---|---|
| WFS_001 | 55.62771 | -132.87711 | AG | CC |
| WFS_002 | 55.48210 | -132.81940 | AA | CT |

**Decimal degrees**, with west and south written negative. Columns named
`Latitude`/`Longitude` (or `Lat`/`Lon`, `Lat`/`Long`, `Y`/`X` spelled out as
`ycoord`/`xcoord`) are found automatically; anything else you pick in the
sidebar.

**UTM works too.** Name the columns `Easting` and `Northing` and the app asks
for the zone and hemisphere, then converts to WGS84 for you.

On a Format B file — one row per replicate — put the same position on every row
of a sample; the app takes one position per sample. Samples with no position are
left off the map and reported nowhere else, so a partly-mapped file is fine.

### Failed samples: upload raw, do not pre-curate

**Give it everything, including the samples where nothing amplified.** The
call-rate filters in the sidebar are the curation step, and doing it in the app
is better than doing it in Excel first, because the app shows you what was
dropped and lets you move the line and watch the answer move.

What happens if you *don't* filter is the reason it matters. A sample with no
called loci matches nothing, so every method hands it its own individual —
turning a plate of failed extractions into a row of phantom animals and
inflating your population estimate. On one real 123-sample run, leaving the
filter off took the count from 62 individuals to 81, purely from 19 samples
that had no usable data.

The default (`Minimum sample call rate` 0.5) is a sensible line, and the reason
is visible in **Panel power &rarr; per-sample**, which reports P(ID)<sub>sib</sub>
for each sample using only the loci that sample actually has:

| Call rate | Median loci typed | Median P(ID)<sub>sib</sub> | Can it be told from a sibling? |
|---|---|---|---|
| 0&ndash;25% | 0 | 1.0 | No |
| 25&ndash;50% | 10 | 0.028 | No |
| 50&ndash;75% | 24 | 8 &times; 10<sup>&minus;5</sup> | Yes |
| 75&ndash;100% | 33 | 3 &times; 10<sup>&minus;6</sup> | Yes |

The cliff is sharp and it is a property of your panel, not a rule of thumb.
Rather than trusting a blanket threshold, check the per-sample figure: anything
below about P(ID)<sub>sib</sub> = 0.01 can be trusted to be distinguishable from
a close relative, and a marginal sample above the cutoff can be judged on its
own merits instead of being discarded by a round number.

Samples that fail the filter are listed by name, not silently dropped, and they
are excluded from allele-frequency estimation as well as from matching.

### What the app does for you

- **Finds the loci** — columns that parse as diploid genotypes with a consistent
  alphabet, then lists them in the sidebar for you to correct. Every column is
  offerable, so nothing it missed is out of reach.
- **Finds the sample column** — prefers identifier-looking names and ignores
  pure numbers like read counts.
- **Sorts alleles** so order never creates a false difference.
- **Drops dead loci** — monomorphic, or below your minimum call rate — and says
  which.

### Two things to do yourself

1. **Remove species-diagnostic and sex markers** from the locus list. They are
   near-fixed within a species, so they add no power and distort the allele
   frequencies every probability depends on.
2. **Set "analyze separately by"** to your species or population column, so
   samples from different groups are never compared and each gets its own
   frequencies.

## Running it locally

```bash
git clone https://github.com/taaltree/genoID.git
```

```bash
cd genoID && Rscript -e "shiny::runApp('app', port = 4599)"
```

Needs `shiny`, `bslib`, `DT`, `ggplot2` and `leaflet`.

## Using the R functions directly

`app/genoID_core.R` is **base R only** — no Shiny, no compiled packages, no
dependencies at all. Source it anywhere.

```r
source("app/genoID_core.R")

raw  <- gid_read("genotypes.csv")
loci <- setdiff(gid_detect_loci(raw, exclude = "SampleID"), c("SexMarker", "DIAG01"))
gt   <- gid_matrix(raw, "SampleID", loci)

res <- gid_by_group(gt, raw$Species, gid_method_sethi,
                    dropout = 0.005, false_allele = 0.002,
                    relationships = c("unrelated", "full_sib", "parent_offspring"),
                    min_loci = 15)

write.csv(res$assignment, "individuals.csv", row.names = FALSE)
```

```r
gid_estimate_error(gt, reps)            # picks the best method your data supports
gid_error_ml(rep_gt, sample_ids)        # dropout + false allele by maximum likelihood
gid_error_from_fis(gt)                  # rough dropout from a heterozygote deficit
gid_propagate_error(d, f, rule = "taberlet")   # per-reaction -> consensus residual
```

With one row per PCR reaction, pass the reactions and the methods condition on
every one, each sample at its own error scale:

```r
reps <- list(gt = rep_gt, sample = rep_sample_ids)   # one row per reaction
gt   <- gid_consensus_from_reps(rep_gt, rep_sample_ids)
e    <- gid_estimate_error(gt, reps)                  # per-reaction rates
a    <- gid_method_sethi(gt, dropout = e$dropout, false_allele = e$false_allele,
                         reps = reps)                 # every reaction
r    <- gid_propagate_error(e$dropout, e$false_allele)
b    <- gid_method_sethi(gt, dropout = max(r[["dropout"]], 1e-4),
                         false_allele = max(r[["false_allele"]], 1e-4))   # consensus
gid_compare_routes(a, b)$moved                         # samples the routes disagree on
gid_replicate_discordance(rep_gt, rep_sample_ids)      # how noisy each sample is
```

## Where do the error rates come from?

Guessing the dropout and false-allele rates is the weakest link in the analysis,
and you usually don't have to. Press **Estimate these from my data** in the
sidebar; the app picks the best method your data supports. With replicates in
the file, it measures them on its own when you run the analysis, unless you have
set the sliders yourself.

| What you have | Method | On the AITRC panel |
|---|---|---|
| **PCR replicates** | Maximum likelihood over the replicates, integrating over the unknown true genotype. No consensus involved, so nothing is circular. Returns confidence intervals. | dropout **2.9%** [2.5–3.4], false allele **1.5%** [1.3–1.7] per reaction |
| **No replicates** | Heterozygote deficit: dropout turns hets into homs, so *d* ≈ 1 − H<sub>o</sub>/H<sub>e</sub> = F<sub>IS</sub> | ≈ 6.0% — right order of magnitude, not a number to report |
| **Nothing** | Literature values, then check with the sensitivity panel | answer moves by 1 animal across a 100× change |

Two replicates on a *subset* of samples is enough — you don't have to replicate
everything, and the app uses whatever replicates exist even when you're
analyzing consensus calls.

> **The rate depends on what you're analyzing.** Raw replicates want the
> per-reaction rate (a few percent). Consensus calls want the much smaller
> residual rate that survives the multi-tube rule (a few tenths of a percent) —
> a factor of ~300 apart on this panel. The app measures the per-reaction rate,
> works out the residual rate exactly, and uses each for its own route. No
> residual rate is set below one in ten thousand: a swapped tube does not care
> how many reactions agreed.

**A trap worth knowing.** The standard Broquet & Petit estimator scores each
replicate against a consensus. If your consensus calls a homozygote only when
all replicates agree, no replicate can ever be seen disagreeing with a
homozygous consensus, and the false-allele rate is forced to **exactly zero**.
On this dataset that estimator returns 0.00% where maximum likelihood, given the
same replicates, returns 1.5%.

## Linking samples to individuals, with confidence

The **Sample map** tab is the output most people actually need: one row per
sample, the animal it belongs to, and how much to trust that. For each sample it
asks two questions and reports the gap between them.

| Column | Means |
|---|---|
| **Support** | The *weakest* link holding this sample inside its own animal. A sample joined through one marginal pair is far less secure than one matching every other member outright. Singletons get 1 |
| **Nearest rival** | The *strongest* link to a sample assigned to a **different** animal — what it nearly matched instead |
| **Margin** | Support − rival. Wide means the assignment survives any reasonable cutoff; near zero means this sample is the reason your answer depends on where you drew the line |
| **Directly matched** | `1 of 8` means it is held into an 8-mate cluster by a single link, everything else inferred through intermediates |
| **Status** | Confident · Check · Uncertain · Underpowered |

`Underpowered` is deliberately separate: a sample with too few loci to resolve
either way isn't a borderline call, it's **no call**, and treating it as
borderline hides the real problem.

On the AITRC data 13 of 104 samples are flagged — and they are exactly the
samples the six methods disagreed about, arrived at independently. The worst
sits at posterior 0.9989 against a five-sample individual, just under the 0.999
cutoff: one sample, one thousandth of a probability, and a whole animal either
appears or does not.

It works for every method: the likelihood methods are scored on the posterior
scale, the counting methods on mismatching loci. A method that clusters
internally and exposes no per-pair score borrows the mismatch counts, which are
a property of the data rather than of any method.

## Every sample coming back as a unique individual?

That is a symptom, not a result. It almost always means the posterior cutoff is
set above what your panel can supply. A panel carries a finite amount of
evidence, so there is a maximum posterior no pair can exceed — the app reports
that ceiling and suggests a workable cutoff.

## How to choose the posterior cutoff

0.999 is a convention, not a result. It doesn't have to be a guess: the method
gives every pair a probability, so the cost of any cutoff can be added up.

```
E[false merges] = sum over ACCEPTED pairs of (1 - p)     two animals joined
E[missed pairs] = sum over REJECTED pairs of p           one animal split
```

Both are in pairs, so the choice becomes a stated trade. The Individuals tab
plots both curves against the cutoff and tells you where you are on them.

On the AITRC panel:

| Cutoff | Individuals | E[false merges] | E[missed pairs] |
|---|---|---|---|
| 0.95 | 56 | 0.03 | 2.4 |
| 0.999 | 59 | 0.02 | 9.4 |
| 0.9999 | 66 | 0.01 | 37 |
| 0.99999 | **100** | 0.00 | **241** |

Every cutoff from 0.5 to 0.9995 costs under one expected false merge while the
count moves only 55 → 59, so within that range the cutoff is barely doing any
work. Push to 0.99999 and you buy a 0.02 reduction in false merges for 241
missed pairs — that is collapse, not caution.

**A defensible procedure:** check the ceiling first; set an error budget in
pairs rather than a probability; take the most permissive cutoff inside it;
check the answer sits on a plateau; report the range, not the point. *"55 to 59
animals across every defensible cutoff"* is a better result than a single
number.

The two errors aren't symmetric: a split invents an animal **and** manufactures
a spurious "caught once" record, which is what capture-recapture estimators are
most sensitive to. For abundance work, splitting is usually the more damaging
direction.

## Jargon

Tables show plain-English headings, not the raw column names. Hover any heading
for a definition. The full list — label, the column name it carries in
downloaded files, and what it means — is on the Methods tab under **Glossary of
every column**.

## How the app is laid out

Pick a method at the top of the sidebar. The sidebar then shows **only that
method's settings**, the Individuals tab reports **only that method's answer**,
with a plain-language description and a button through to its full derivation
and equations on the Methods tab.

Every method still runs behind the scenes, so the **Method comparison** tab
shows all five side by side whatever you picked. Switching method on the sidebar
re-reports the Individuals tab instantly, without recomputing.

Results are not computed until you press **Identify individuals** — the Data &
QC and Panel power tabs fill in on their own, but Individuals and Method
comparison wait for you. Both say so, and if a run fails they show the actual R
error rather than going blank.

## Deploying your own copy

Fork it, then enable **Settings → Pages → Source: GitHub Actions**. Every push
to `main` that touches `app/` rebuilds and redeploys the site. The workflow is
in `.github/workflows/deploy.yml`.

## Repository layout

```
app/
  app.R                  the Shiny app
  genoID_core.R          all methods and statistics; base R only, no dependencies
  methods_page.R         the Methods & mathematics tab
  demo/                  simulated example dataset
run_analysis.R           batch pipeline, writes result CSVs and figures
make_demo_data.R         regenerates the example (and checks it is recoverable)
compare_reps_vs_consensus.R  simulation experiment: replicates vs consensus
report/                  written report (not committed; see .gitignore)
```

No real genotype data is committed to this repository.

## References

Broquet & Petit (2004) *Mol Ecol* 13:3601 · Galpern et al. (2012) *Mol Ecol Res*
12:771 · Paetkau (2003) *Mol Ecol* 12:1375 · Peakall & Smouse (2006, 2012) ·
Sethi et al. (2016) *R Soc Open Sci* 3:160457 · Taberlet et al. (1996) *NAR*
24:3189 · Waits, Luikart & Taberlet (2001) *Mol Ecol* 10:249

## License

MIT
