## methods_page.R -------------------------------------------------------------
##
## The "Methods & mathematics" tab. Plain-language description of every method
## first, then the actual model. Equations are LaTeX rendered by KaTeX; if the
## CDN is unreachable the raw LaTeX stays legible rather than disappearing.
## ---------------------------------------------------------------------------

katex_head <- function() {
  tagList(
    tags$link(rel = "stylesheet",
              href = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css",
              crossorigin = "anonymous"),
    tags$script(defer = NA, src = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js",
                crossorigin = "anonymous"),
    tags$script(defer = NA,
                src = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js",
                crossorigin = "anonymous",
                onload = "renderMathInElement(document.body,{delimiters:[
                  {left:'\\\\[',right:'\\\\]',display:true},
                  {left:'\\\\(',right:'\\\\)',display:false}],throwOnError:false});"),
    tags$script(HTML(
      "document.addEventListener('shiny:value', function(){
         if (window.renderMathInElement) setTimeout(function(){
           renderMathInElement(document.body,{delimiters:[
             {left:'\\\\[',right:'\\\\]',display:true},
             {left:'\\\\(',right:'\\\\)',display:false}],throwOnError:false});
         }, 30);
       });")),
    tags$style(HTML("
      .m-body{max-width:47rem;}
      .m-body p{margin:.7rem 0;line-height:1.62;}
      .m-lead{font-size:1.02rem;color:#33474d;}
      .m-eq{margin:1rem 0;overflow-x:auto;padding:.2rem 0;}
      .m-h{font-weight:600;color:#12262a;margin:2.1rem 0 .5rem;font-size:1.12rem;}
      .m-sub{font-weight:600;color:#12262a;margin:1.4rem 0 .35rem;font-size:.97rem;}
      .m-tag{display:inline-block;font-size:.66rem;font-weight:700;letter-spacing:.08em;
        text-transform:uppercase;padding:.16em .55em;border-radius:2px;margin-left:.5rem;
        vertical-align:.15em;}
      .m-tag.rec{background:#f6eddc;color:#a8762a;}
      .m-tag.chk{background:#e3efef;color:#0f6b73;}
      .m-def{border-left:2px solid #dae4e2;padding-left:1rem;margin:1rem 0;font-size:.93rem;
        color:#4a6067;}
      .m-def b{color:#12262a;}
      table.m-tbl{border-collapse:collapse;margin:1rem 0;font-size:.88rem;}
      table.m-tbl th{text-align:left;font-size:.7rem;letter-spacing:.06em;text-transform:uppercase;
        color:#6b7a80;font-weight:600;padding:.2rem .9rem .35rem 0;border-bottom:1px solid #dae4e2;}
      table.m-tbl td{padding:.35rem .9rem .35rem 0;border-bottom:1px solid #eef3f2;
        vertical-align:top;}
      .m-note{background:#f6eddc;border-left:3px solid #a8762a;padding:.7rem 1rem;
        margin:1.1rem 0;font-size:.9rem;border-radius:0 3px 3px 0;}
      .m-toc a{display:block;padding:.22rem 0;font-size:.86rem;text-decoration:none;}
    "))
  )
}

.eq <- function(x) tags$div(class = "m-eq", paste0("\\[", x, "\\]"))
.h  <- function(x, id = NULL, tag = NULL, tagclass = "rec")
  tags$div(class = "m-h", id = id, x,
           if (!is.null(tag)) tags$span(class = paste("m-tag", tagclass), tag))
.s  <- function(x) tags$div(class = "m-sub", x)
.p  <- function(...) tags$p(...)


methods_panel <- function() {
  nav_panel(
    "Methods", icon = icon("book-open"),
    layout_columns(
      col_widths = c(3, 9),

      card(card_header("Contents"), class = "m-toc",
        tags$a(href = "#m-problem", "The problem"),
        tags$a(href = "#m-prep", "Preparing the data"),
        tags$a(href = "#m-exact", "1. Exact match"),
        tags$a(href = "#m-thresh", "2. Mismatch threshold"),
        tags$a(href = "#m-genalex", "3. GenAlEx Matches"),
        tags$a(href = "#m-am", "4. allelematch"),
        tags$a(href = "#m-lr", "5. Likelihood ratio"),
        tags$a(href = "#m-graph", "From pairs to individuals"),
        tags$a(href = "#m-err", "Measuring error rates"),
        tags$a(href = "#m-set", "Choosing settings"),
        tags$a(href = "#m-refs", "References")),

      card(card_header("Methods & mathematics"), tags$div(class = "m-body",

        # ------------------------------------------------------------ problem
        .h("The problem", id = "m-problem"),
        .p(class = "m-lead",
          "You have one genotype per sample. Several samples may have come from the ",
          "same animal, and you do not know which. If genotyping were perfect this ",
          "would be a lookup: identical genotype, same animal. Three things break that."),
        tags$div(class = "m-def",
          .p(tags$b("Allelic dropout."), " One allele of a heterozygote fails to amplify, ",
             "so a true ", tags$code("A/G"), " is scored ", tags$code("A/A"), ". Two samples ",
             "from one animal then differ at that locus. Dropout is the dominant error in ",
             "degraded material such as scat or hair."),
          .p(tags$b("False alleles."), " An allele appears that the animal does not carry, ",
             "usually from PCR error or contamination."),
          .p(tags$b("Relatives."), " Siblings share a lot of their genome. Two different ",
             "animals from the same litter can match at more loci than two samples from one ",
             "animal that suffered dropout. This is why a method that only counts ",
             "differences cannot get the hard cases right.")),
        .p("Every method below is a different answer to: how many differing loci are ",
           "still consistent with one animal?"),

        # --------------------------------------------------------------- prep
        .h("Preparing the data", id = "m-prep"),
        .s("Genotypes are order-normalised"),
        .p("Alleles are sorted within each genotype, so ", tags$code("TC"), " and ",
           tags$code("CT"), " are one genotype. Panels commonly contain a handful of ",
           "reversed cells, and compared as raw text those count as mismatches."),
        .s("Species-diagnostic and sex markers are excluded"),
        .p("A locus that is nearly fixed between two species carries almost no information ",
           "within either one, and its allele frequencies are meaningless when the two are ",
           "pooled. Use those loci to assign species, then drop them. The sex marker is kept ",
           "aside as an independent check: two samples from one animal cannot disagree about sex."),
        .s("Analyses are blocked by species or population"),
        .p("Samples in different blocks are never compared, and allele frequencies are ",
           "estimated within each block. Pooling populations inflates apparent heterozygosity ",
           "and makes every probability below wrong in the optimistic direction."),
        .s("Allele frequencies"),
        .p("For a locus with alleles \\(1 \\dots m\\), \\(p_i\\) is the frequency of allele ",
           "\\(i\\) counted over all sampled genotypes. Under Hardy-Weinberg the genotype ",
           "probabilities are \\(P(ii) = p_i^2\\) and \\(P(ij) = 2p_ip_j\\) for \\(i \\ne j\\)."),
        tags$div(class = "m-note",
          tags$b("Recaptures bias allele frequencies."), " An animal sampled sixteen times ",
          "contributes its genotype sixteen times. This app estimates frequencies once, ",
          "clusters, then re-estimates weighting each sample by the reciprocal of its ",
          "cluster size, so each animal counts once."),

        # -------------------------------------------------------------- exact
        .h("1. Exact match", id = "m-exact", tag = "baseline", tagclass = "chk"),
        .p("Two samples are the same individual if they are identical at every locus where ",
           "both were called. Loci missing in either sample are skipped, and a minimum number ",
           "of compared loci is required."),
        .p("With genuinely clean genotypes this is close to correct and it is a useful ",
           "check. Its weakness is that it cannot express doubt: a single dropout splits one ",
           "animal into two, and a pair matching on only 15 shared loci is treated with the ",
           "same confidence as a pair matching on 40."),

        # ---------------------------------------------------------- threshold
        .h("2. Mismatch threshold", id = "m-thresh", tag = "gut check", tagclass = "chk"),
        .p("Allow up to \\(k\\) differing loci. For a pair with \\(L\\) co-typed loci, count"),
        .eq("M = \\sum_{\\ell=1}^{L} \\mathbb{1}\\!\\left[G_{1\\ell} \\neq G_{2\\ell}\\right]"),
        .p("and call the pair a match when \\(M \\le k\\) and \\(L \\ge L_{\\min}\\)."),
        .p("The app also counts ", tags$b("hard mismatches"), " separately: loci where the two ",
           "genotypes share no allele at all. Dropout cannot produce those, since it only ever ",
           "turns a heterozygote into a homozygote for an allele the animal already carries. A ",
           "pair with even one hard mismatch is far more likely to be two animals."),
        .s("Choosing k"),
        .p("With a good panel the distribution of \\(M\\) over all pairs is bimodal: a pile near ",
           "zero (recaptures) and a pile far away (different animals), separated by a gap. Set ",
           "\\(k\\) inside the gap. The sweep on the Method comparison tab shows what happens ",
           "across values of \\(k\\); a flat stretch means the answer is a property of the data, ",
           "a steady slide means \\(k\\) is deciding your result."),

        # ------------------------------------------------------------ genalex
        .h("3. GenAlEx Matches", id = "m-genalex", tag = "familiar", tagclass = "chk"),
        .p("Reproduces the ", tags$b("Multilocus \u2192 Matches"), " routine in GenAlEx ",
           "(Peakall & Smouse 2006, 2012), which many labs already run. It does three things."),
        .p(tags$b("(a) Tabulates the match distribution"), " \u2014 how many pairs differ at ",
           "0, 1, 2, \u2026 loci. This is the table GenAlEx prints, and it is the same ",
           "distribution the threshold method's sweep is built on."),
        .p(tags$b("(b) Reports exact matches as the same individual"), ", and pairs one or two ",
           "loci short as ", tags$i("near matches"), " to inspect by hand rather than merging ",
           "automatically. Near matches are exactly where dropout and siblings live."),
        .p(tags$b("(c) Reports the probability of identity"), ", the justification that the ",
           "panel is powerful enough to trust a match. For one locus,"),
        .eq("P_{\\mathrm{ID}} = \\sum_{i} p_i^{4} \\;+\\; \\sum_{i}\\sum_{j>i} \\left(2p_ip_j\\right)^{2}"),
        .p("This is the chance two ", tags$i("unrelated"), " animals share a genotype at that ",
           "locus. For full siblings, who share alleles by descent,"),
        .eq("P_{\\mathrm{ID}\\text{-}\\mathrm{sib}} = 0.25 + 0.5\\sum_i p_i^{2} + 0.5\\left(\\sum_i p_i^{2}\\right)^{2} - 0.25\\sum_i p_i^{4}"),
        .p("Across loci, assuming they are independent, the values multiply:"),
        .eq("P_{\\mathrm{ID}}^{\\text{total}} = \\prod_{\\ell=1}^{L} P_{\\mathrm{ID},\\ell}"),
        .p("A panel is conventionally called adequate when \\(P_{\\mathrm{ID}\\text{-}\\mathrm{sib}}\\) ",
           "falls below 0.01, and comfortable below 0.001. Quote the sibling value, not the ",
           "unrelated one \u2014 for social or family-structured species it is the honest number, ",
           "and it is often several orders of magnitude larger."),
        tags$div(class = "m-note",
          tags$b("Two things to know about these numbers."), " First, they assume loci are ",
          "independent; two SNPs on one amplicon are not, and multiplying them overstates your ",
          "power. The Data & QC tab reports r\u00b2 between loci so you can drop one of any ",
          "correlated pair before quoting a figure. Second, GenAlEx can also report a ",
          "small-sample-corrected \u201cunbiased\u201d P(ID); this app reports the uncorrected ",
          "Waits et al. (2001) estimators above. At typical sample sizes the difference is far ",
          "smaller than the error introduced by population substructure, which neither ",
          "version accounts for."),

        # ------------------------------------------------------------------ am
        .h("4. allelematch", id = "m-am", tag = "gut check", tagclass = "chk"),
        .p("The published semi-automated method of Galpern et al. (2012). It scores ",
           "dissimilarity between two samples as the fraction of mismatching alleles among ",
           "those compared, clusters samples, then matches every sample back against the ",
           "cluster representatives."),
        .p("Its real contribution is the ", tags$b("unclassified"), " category: samples that ",
           "are genuinely ambiguous get flagged rather than forced into a cluster. It selects ",
           "its own threshold by sweeping the allowed number of mismatching alleles and taking ",
           "the value that leaves fewest samples unresolved."),
        .p("It was designed for microsatellites, where tolerating 2 or 3 mismatching alleles ",
           "out of 20 is natural. On a dense biallelic SNP panel with missing data it tends to ",
           "over-split, because its dissimilarity penalises missing data and an incomplete ",
           "sample cannot reach a perfect score against a complete genotype. Worth running as ",
           "a check; read a large individual count from it with that in mind."),

        # ------------------------------------------------------------------ lr
        .h("5. Likelihood ratio", id = "m-lr", tag = "recommended"),
        .p(class = "m-lead",
          "Instead of counting differences, this asks how much more probable the observed ",
          "pair of genotypes is if the two samples came from one animal than if they came ",
          "from two. It is the only method here that uses an explicit error model, states ",
          "what it is testing against, and returns a probability per pair."),

        .s("The two hypotheses"),
        tags$div(class = "m-def",
          .p(tags$b("\\(H_1\\):"), " one animal with true genotype \\(g\\), observed twice ",
             "through the error process."),
          .p(tags$b("\\(H_0\\):"), " two animals, with genotypes drawn from the population at ",
             "a stated degree of relatedness.")),

        .s("The error model"),
        .p("Let \\(d\\) be the allelic dropout rate and \\(f\\) the false allele rate. An ",
           "observation is generated from a true genotype in two stages:"),
        tags$div(class = "m-def",
          .p(tags$b("Dropout."), " With probability \\(d\\), a true heterozygote ",
             "\\(\\{a,b\\}\\) loses one allele at random and is seen as \\(\\{a,a\\}\\) or ",
             "\\(\\{b,b\\}\\), each with probability \\(d/2\\). Homozygotes are unaffected, ",
             "because losing one of two identical alleles is invisible."),
          .p(tags$b("False allele."), " With probability \\(f\\), one of the two surviving ",
             "alleles, chosen at random, is replaced by an allele drawn from the population ",
             "frequencies.")),
        .p("Writing \\(E[G \\mid g]\\) for the resulting probability of observing genotype ",
           "\\(G\\) given true genotype \\(g\\), the probability of the observed pair under ",
           "each hypothesis is"),
        .eq("P(G_1, G_2 \\mid H_1) \\;=\\; \\sum_{g} P(g)\\, E[G_1 \\mid g]\\, E[G_2 \\mid g]"),
        .eq("P(G_1, G_2 \\mid H_0) \\;=\\; \\sum_{g_1}\\sum_{g_2} P(g_1, g_2)\\, E[G_1 \\mid g_1]\\, E[G_2 \\mid g_2]"),
        .p("Note that \\(H_1\\) does not require the two observations to be identical. A pair ",
           "that differs at one locus is still perfectly consistent with one animal \u2014 it ",
           "just costs likelihood, and how much it costs depends on \\(d\\) and \\(f\\). ",
           "That is the whole difference from a threshold rule."),

        .s("The alternative hypothesis, and why it matters"),
        .p("Two animals related by descent share alleles. Decomposing by the number of alleles ",
           "shared identical-by-descent, with probabilities \\(k_0, k_1, k_2\\),"),
        .eq("P(g_1, g_2) = k_0\\,P(g_1)P(g_2) \\;+\\; k_1\\,P_1(g_1,g_2) \\;+\\; k_2\\,P(g_1)\\,\\mathbb{1}[g_1 = g_2]"),
        .p("where the one-allele-shared term sums over the shared allele \\(s\\) and the two ",
           "independently drawn alleles \\(a\\) and \\(b\\):"),
        .eq("P_1(g_1, g_2) = \\sum_{s}\\sum_{a}\\sum_{b} p_s\\,p_a\\,p_b\\; \\mathbb{1}\\!\\left[g_1 = \\{s,a\\}\\right]\\,\\mathbb{1}\\!\\left[g_2 = \\{s,b\\}\\right]"),
        tags$table(class = "m-tbl",
          tags$thead(tags$tr(tags$th("Alternative"), tags$th("\\(k_0\\)"), tags$th("\\(k_1\\)"),
                             tags$th("\\(k_2\\)"), tags$th("Use when"))),
          tags$tbody(
            tags$tr(tags$td("Unrelated"), tags$td("1"), tags$td("0"), tags$td("0"),
                    tags$td("solitary or well-mixed population")),
            tags$tr(tags$td("Half siblings"), tags$td("0.5"), tags$td("0.5"), tags$td("0"),
                    tags$td("moderate family structure")),
            tags$tr(tags$td("Full siblings"), tags$td("0.25"), tags$td("0.5"), tags$td("0.25"),
                    tags$td("packs, prides, litters, colonies")),
            tags$tr(tags$td("Parent\u2013offspring"), tags$td("0"), tags$td("1"), tags$td("0"),
                    tags$td("known pedigree sampling")))),
        .p("For a group-living species the samples competing to be \u201ca different animal\u201d ",
           "are usually packmates, not random members of the population. Testing against ",
           "unrelated animals then overstates the evidence for a match, because it is asking ",
           "the wrong question. Full siblings is the conservative and usually the correct choice."),

        .s("Combining loci"),
        .p("Loci are treated as independent, so evidence adds on the log scale over the loci ",
           "where both samples were called:"),
        .eq("\\log_{10}\\mathrm{LR} \\;=\\; \\sum_{\\ell \\,:\\, \\text{both typed}} \\log_{10} \\frac{P(G_{1\\ell}, G_{2\\ell} \\mid H_1)}{P(G_{1\\ell}, G_{2\\ell} \\mid H_0)}"),
        .p("Missing loci contribute nothing rather than counting against the pair, which is why ",
           "a minimum number of compared loci is enforced separately."),

        .s("From evidence to a decision"),
        .p("A likelihood ratio is not a probability. To decide, combine it with a prior ",
           "\\(\\pi\\) that a randomly chosen pair of samples is a recapture:"),
        .eq("\\frac{P(H_1 \\mid G_1, G_2)}{P(H_0 \\mid G_1, G_2)} \\;=\\; \\frac{\\pi}{1-\\pi} \\times \\mathrm{LR}"),
        .p("and accept the match when the posterior probability exceeds your cutoff, ",
           "0.999 by default."),
        .p("The prior is estimated from the data itself, self-consistently: start from one ",
           "recapture per sample, count how many pairs match at that prior, feed the count ",
           "back as \\(\\pi = (m+1)/\\binom{n}{2}\\), and repeat until it stops moving. It ",
           "converges in a few iterations. You can also set it by hand if you have a prior ",
           "expectation of the recapture rate."),
        tags$div(class = "m-note",
          tags$b("A panel carries a finite amount of evidence."), " The strongest possible ",
          "\\(\\log_{10}\\mathrm{LR}\\) is roughly \\(-\\log_{10} P_{\\mathrm{ID}\\text{-}\\mathrm{sib}}\\). ",
          "Set a posterior cutoff above what your loci can supply and nothing matches, so ",
          "every sample is reported as its own individual \u2014 which looks like a result and ",
          "is not one. The app checks your cutoff against the best evidence in your own data ",
          "and warns you on the Individuals tab."),

        # --------------------------------------------------------------- graph
        .h("From matched pairs to individuals", id = "m-graph"),
        .p("Every method produces a set of matching pairs. Turning that into a list of animals ",
           "means treating the pairs as edges in a graph and taking groups."),
        tags$div(class = "m-def",
          .p(tags$b("Single linkage"), " (connected components). If A matches B and B matches ",
             "C, all three are one animal even if A and C were never matched. This recovers ",
             "recaptures linked through an intermediate sample, but one wrong edge merges two ",
             "animals."),
          .p(tags$b("Complete linkage"), " (maximal cliques). Every member must match every ",
             "other member. Conservative; splits when a pair happens to fall below threshold ",
             "because of missing data.")),
        .p("Run both. If the answer moves, some clusters are held together by a single edge ",
           "and deserve a look. The app always reports ", tags$b("internally inconsistent ",
           "clusters"), " \u2014 groups where some members do not match each other \u2014 with ",
           "their completeness, the fraction of within-cluster pairs that actually matched. ",
           "A cluster below 1.0 is the single most useful diagnostic on the page."),

        # ----------------------------------------------------------------- err
        .h("Measuring genotyping error", id = "m-err"),
        .p("The likelihood ratio needs \\(d\\) and \\(f\\). If you have replicate genotypes, ",
           "measure them rather than guessing."),
        .s("Scoring replicates against a consensus"),
        .p("The standard approach (Broquet & Petit 2004) compares each replicate to a ",
           "reference genotype:"),
        .eq("\\hat{d} = \\frac{\\#\\{\\text{true het scored as hom}\\}}{\\#\\{\\text{true het}\\}}, \\qquad \\hat{f} = \\frac{\\#\\{\\text{true hom carrying a foreign allele}\\}}{\\#\\{\\text{true hom}\\}}"),
        tags$div(class = "m-note",
          tags$b("This breaks when the consensus rule is strict."), " If your consensus calls a ",
          "homozygote only when all replicates agree, then no replicate can ever be seen ",
          "disagreeing with a homozygous consensus, and \\(\\hat{f}\\) is forced to exactly ",
          "zero. The estimator is measuring its own reference."),
        .s("Maximum likelihood from the replicates alone"),
        .p("Better: never use a consensus. Treat the true genotype at each sample and locus as ",
           "unknown, integrate over it, and maximise over the two rates. For replicate ",
           "observations \\(o_1 \\dots o_R\\),"),
        .eq("\\mathcal{L}(d, f) = \\prod_{\\text{sample},\\,\\text{locus}} \\; \\sum_{g} P(g) \\prod_{r=1}^{R} E[o_r \\mid g, d, f]"),
        .p("This is what ", tags$code("gid_error_ml()"), " computes, with profile-likelihood ",
           "confidence intervals."),
        .s("Consensus calls are much cleaner than single reactions"),
        .p("A multi-tube rule \u2014 accept a heterozygote on two replicates, a homozygote only ",
           "on three \u2014 converts a few percent dropout per reaction into a fraction of a ",
           "percent in the consensus. If you are uploading consensus genotypes, the rates the ",
           "model wants are the residual consensus rates, not the per-reaction ones. ",
           tags$code("gid_propagate_error()"), " simulates a consensus rule to get them."),

        # ----------------------------------------------------------------- set
        .h("Choosing settings", id = "m-set"),
        tags$table(class = "m-tbl",
          tags$thead(tags$tr(tags$th("Setting"), tags$th("What it does"), tags$th("Guidance"))),
          tags$tbody(
            tags$tr(tags$td("Alternative hypothesis"),
                    tags$td("What \u201ctwo different animals\u201d means"),
                    tags$td("Full siblings for group-living species. This is usually the ",
                            "setting that moves the answer most.")),
            tags$tr(tags$td("Posterior cutoff"),
                    tags$td("Confidence required to accept a match"),
                    tags$td("0.999 is a reasonable default. Check the warning that your panel ",
                            "can actually reach it.")),
            tags$tr(tags$td("Dropout, false allele"),
                    tags$td("The error model"),
                    tags$td("Measure from replicates if you can. For multi-replicate consensus ",
                            "genotypes both are usually well under 0.005.")),
            tags$tr(tags$td("Minimum loci per pair"),
                    tags$td("Refuses to judge sparse comparisons"),
                    tags$td("Set so that the loci compared could reach your cutoff on their ",
                            "own \u2014 the per-sample power plot shows which samples cannot.")),
            tags$tr(tags$td("Cluster rule"), tags$td("Single vs complete linkage"),
                    tags$td("Run both; report the difference if there is one.")))),
        .p("Then use the sensitivity panel on the Method comparison tab. It re-runs the whole ",
           "analysis across a grid of these settings. If the individual count barely moves, it ",
           "is a property of your data; if it swings, it is a property of your choices, and ",
           "the range belongs in your results."),

        # ---------------------------------------------------------------- refs
        .h("References", id = "m-refs"),
        tags$div(class = "m-def", style = "font-size:.88rem;",
          .p("Broquet, T. & Petit, E. (2004). Quantifying genotyping errors in noninvasive ",
             "population genetics. ", tags$i("Molecular Ecology"), " 13, 3601\u20133608."),
          .p("Galpern, P., Manseau, M., Hettinga, P., Smith, K. & Wilson, P. (2012). ",
             "allelematch: an R package for identifying unique multilocus genotypes where ",
             "genotyping error and missing data may be present. ",
             tags$i("Molecular Ecology Resources"), " 12, 771\u2013778."),
          .p("Paetkau, D. (2003). An empirical exploration of data quality in DNA-based ",
             "population inventories. ", tags$i("Molecular Ecology"), " 12, 1375\u20131387."),
          .p("Peakall, R. & Smouse, P.E. (2006, 2012). GenAlEx 6: genetic analysis in Excel. ",
             tags$i("Molecular Ecology Notes"), " 6, 288\u2013295; ", tags$i("Bioinformatics"),
             " 28, 2537\u20132539."),
          .p("Taberlet, P. et al. (1996). Reliable genotyping of samples with very low DNA ",
             "quantities using PCR. ", tags$i("Nucleic Acids Research"), " 24, 3189\u20133194."),
          .p("Waits, L.P., Luikart, G. & Taberlet, P. (2001). Estimating the probability of ",
             "identity among genotypes in natural populations: cautions and guidelines. ",
             tags$i("Molecular Ecology"), " 10, 249\u2013256."))
      ))
    )
  )
}
