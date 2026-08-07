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
| **Method comparison** | All six methods on the same data at once, adjusted Rand index between them, the samples they disagree about, and a sensitivity grid across 18 combinations of settings |
| **Methods** | Notation, a worked example carried through all six methods, the statistical model behind each one, and a glossary defining every column the app can show you |

Everything downloads as CSV, plus a zip of the whole run, the settings you used,
and an R script that reproduces it.

## The six methods

| Method | What it does | Use it for |
|---|---|---|
| **Exact match** | Identical at every co-typed locus | Baseline sanity check |
| **Mismatch threshold** | Allow up to *k* differing loci; sweeps *k* to show how sensitive the answer is | The field standard (Paetkau 2003) |
| **GenAlEx Matches** | Reproduces GenAlEx's Multilocus → Matches: match distribution, near-match list, P(ID) / P(ID)<sub>sib</sub> | Continuity with what your lab already runs |
| **allelematch** | Galpern et al. (2012), with its own automatic threshold selection | Independent published check |
| **Likelihood ratio** | Explicit dropout / false-allele error model, tested against a stated alternative (unrelated, half sib, full sib, parent–offspring), returning a posterior probability per pair | **Recommended** |
| **Sethi et al. (2016)** | Same likelihood, but divided by whichever relationship best explains the pair — so you don't pick one — and decided on evidence alone (Λ > 1) | **Recommended** |

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

## Try it without your own data

Click **"or load the example dataset"** in the sidebar. It's 55 simulated
samples from 28 known individuals across two species, generated with 3% allelic
dropout and 1% false alleles. The true answer is in the `TrueIndividual` column,
so you can score any method against it:

| Method | Individuals found (truth: 28) | Agreement with truth |
|---|---|---|
| Exact match | 47 | 0.22 |
| GenAlEx Matches | 47 | 0.22 |
| Mismatch threshold | 35 | 0.86 |
| allelematch | 31 | 0.95 |
| Likelihood ratio | 30 | 0.98 |
| **Sethi et al. (2016)** | **28** | **1.00** |

Sethi et al. recovers the truth exactly. That gap — 47 animals versus 28 — is
what a 3% dropout rate does to methods that only count differences.

## PCR replicates

If your file has several rows per sample — one per PCR replicate — the app
detects it and offers to **use every replicate observation directly** instead of
collapsing them to a consensus first. The two likelihood methods extend
naturally: for sample *i* at locus *ℓ* with replicate observations
*o*₁…*o*<sub>R</sub>, the genotype likelihood is

```
L_i(g) = prod_r P(o_r | g)
```

and it drops straight into the same two hypotheses. Missing observations
contribute a factor of 1, which makes both hypotheses collapse to the same
marginal — so missing data says nothing rather than needing a special rule.

**Does it help?** Across 180 simulations with a known answer (20–40 loci, 3–25%
dropout), the replicate route landed closer to the truth in 114 cases, the
consensus route in 3, and they tied in 63. The gain appears where you are
marginal:

| Situation | What replicates buy |
|---|---|
| Plenty of loci, low dropout | Nothing measurable — both already recover the truth |
| 30 loci, 25% dropout | Agreement with truth 0.73 → **0.97** |
| Too few loci for the job | Not enough; replicates cannot manufacture information the panel never had |

What improves unconditionally is the **margin** — the gap between the weakest
true recapture and the strongest coincidental match. It roughly doubled in every
condition tested. That does not change an answer that was already right, but it
means the answer survives a worse choice of threshold.

Reproduce the experiment with `Rscript compare_reps_vs_consensus.R`.

> **Watch the error rates.** With replicates you supply *per-reaction* dropout
> and false-allele rates (a few percent). With consensus genotypes you supply the
> much smaller residual rates that survive the multi-tube rule (a few tenths of a
> percent). Mixing them up is the easiest way to get this wrong, so the app says
> which one it wants.

Only the two likelihood methods can use replicates. Exact matching, the mismatch
threshold, GenAlEx and allelematch compare one genotype per sample, so they keep
running on a consensus and the comparison tab stays meaningful.

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

A Format B file can be analysed **either way**, and you can switch between them
in the sidebar to see whether it changes anything.

### Writing a genotype

| Thing | Write it like this | Notes |
|---|---|---|
| SNP genotype | `AG`, `CC`, `T/C` | Two alleles in one cell; separator optional |
| Microsatellite | `120/124`, `120\|124` | Separator **required**, so the two sizes can be told apart |
| Missing | `00`, `NA`, `-`, `?`, blank | All recognised. Missing loci are skipped, not counted against a pair |
| Allele order | `AG` = `GA` | Sorted before comparison, so reversed cells never become false differences |
| Quality flags | `AG*`, `CT?` | Flag stripped, genotype kept. If a flag means "do not trust this", set those cells to missing first |

> **One column per locus, not two.** If your file has `LOC01a` and `LOC01b` as
> separate columns — the layout GenAlEx and allelematch use internally — join
> each pair into a single column first. The app converts back internally when it
> hands data to allelematch.

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
2. **Set "analyse separately by"** to your species or population column, so
   samples from different groups are never compared and each gets its own
   frequencies.

## Running it locally

```bash
git clone https://github.com/taaltree/genoID.git
```

```bash
cd genoID && Rscript -e "shiny::runApp('app', port = 4599)"
```

Needs `shiny`, `bslib`, `DT`, `ggplot2`. `allelematch` is optional — without it
that method falls back to a built-in equivalent instead of failing.

## Using the R functions directly

`app/genoID_core.R` is **base R only** — no Shiny, no compiled packages, no
dependencies at all. Source it anywhere.

```r
source("app/genoID_core.R")

raw  <- gid_read("genotypes.csv")
loci <- setdiff(gid_detect_loci(raw, exclude = "SampleID"), c("SexMarker", "DIAG01"))
gt   <- gid_matrix(raw, "SampleID", loci)

res <- gid_by_group(gt, raw$Species, gid_method_lr,
                    dropout = 0.005, false_allele = 0.002,
                    kinship = "full_sib", post_cut = 0.999, min_loci = 15)

write.csv(res$assignment, "individuals.csv", row.names = FALSE)
```

```r
gid_estimate_error(gt, reps)            # picks the best method your data supports
gid_error_ml(rep_gt, sample_ids)        # dropout + false allele by maximum likelihood
gid_error_from_fis(gt)                  # rough dropout from a heterozygote deficit
gid_propagate_error(d, f, rule = "taberlet")   # per-reaction -> consensus residual
```

## Where do the error rates come from?

Guessing the dropout and false-allele rates is the weakest link in the analysis,
and you usually don't have to. Press **Estimate these from my data** in the
sidebar; the app picks the best method your data supports.

| What you have | Method | On the AITRC panel |
|---|---|---|
| **PCR replicates** | Maximum likelihood over the replicates, integrating over the unknown true genotype. No consensus involved, so nothing is circular. Returns confidence intervals. | dropout **2.9%** [2.5–3.4], false allele **1.5%** [1.3–1.7] per reaction |
| **No replicates** | Heterozygote deficit: dropout turns hets into homs, so *d* ≈ 1 − H<sub>o</sub>/H<sub>e</sub> = F<sub>IS</sub> | ≈ 6.0% — right order of magnitude, not a number to report |
| **Nothing** | Literature values, then check with the sensitivity panel | answer moves by 1 animal across a 100× change |

Two replicates on a *subset* of samples is enough — you don't have to replicate
everything, and the app uses whatever replicates exist even when you're
analysing consensus calls.

> **The rate depends on what you're analysing.** Raw replicates want the
> per-reaction rate (a few percent). Consensus calls want the much smaller
> residual rate that survives the multi-tube rule (a few tenths of a percent) —
> a factor of ~300 apart on this panel. The app measures the per-reaction rate
> and converts automatically when you're on consensus calls.

**A trap worth knowing.** The standard Broquet & Petit estimator scores each
replicate against a consensus. If your consensus calls a homozygote only when
all replicates agree, no replicate can ever be seen disagreeing with a
homozygous consensus, and the false-allele rate is forced to **exactly zero**.
On this dataset that estimator returns 0.00% where maximum likelihood, given the
same replicates, returns 1.5%.

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
shows all six side by side whatever you picked. Switching method on the sidebar
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
