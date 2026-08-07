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

## Your data format

One row per sample, one column per locus. Genotypes as two alleles in one cell
(`AG`) or separated (`120/124`). Missing data as `00`, `NA`, or blank. CSV, TSV,
or Excel. Loci are detected automatically and shown in the sidebar for you to
correct.

Two things to do before you trust the output:

1. **Remove species-diagnostic and sex markers** from the locus list. They are
   near-fixed within a species, so they add no power and distort the allele
   frequencies every probability depends on.
2. **Set "analyse separately by"** to your species or population column, so
   samples from different groups are never matched and each gets its own
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

If you have replicate genotypes, measure your error rates instead of guessing:

```r
gid_error_ml(rep_gt, sample_ids)        # dropout + false allele, by maximum likelihood
gid_propagate_error(d, f, rule = "taberlet")   # residual rate after a multi-tube consensus
```

Do **not** estimate the false-allele rate by scoring replicates against a
consensus that required unanimous replicates to call a homozygote — that forces
the estimate to exactly zero. The Methods page in the app explains why.

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
