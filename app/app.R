## genoID -- individual identification from consensus genotypes ---------------
##
## Upload a genotype table, get unique individuals. Species-agnostic: works for
## any diploid panel (SNPs or microsatellites), any number of loci.
##
## Run locally:  shiny::runApp("app")
## Deployed to GitHub Pages as a WebAssembly build -- see
## .github/workflows/deploy.yml. Everything runs in the visitor's browser.
##
## Note: shinylive scans this source for package references, comments included.
## Never name a package the app does not actually use in double-colon form,
## even inside a comment -- every visitor's browser would download it on
## startup.
## ---------------------------------------------------------------------------

library(shiny); library(bslib); library(DT); library(ggplot2)

core <- c("genoID_core.R", "app/genoID_core.R", "../app/genoID_core.R")
source(core[file.exists(core)][1])
mp <- c("methods_page.R", "app/methods_page.R", "../app/methods_page.R")
source(mp[file.exists(mp)][1])
gl <- c("glossary.R", "app/glossary.R", "../app/glossary.R")
source(gl[file.exists(gl)][1])

options(shiny.maxRequestSize = 60 * 1024^2)

INK <- "#1d3557"; ACCENT <- "#c1502e"; MUTED <- "#6b7a8f"; SOFT <- "#e8edf2"
PAL <- c("#1d3557", "#c1502e", "#4a7c59", "#8b6f9e", "#c9a227")

theme_gid <- function() {
  theme_minimal(base_size = 12.5) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = SOFT, linewidth = 0.4),
          plot.title    = element_text(face = "bold", colour = INK, size = 13.5),
          plot.subtitle = element_text(colour = MUTED, size = 10.5, lineheight = 1.2),
          axis.title    = element_text(colour = MUTED, size = 10.5),
          plot.margin   = margin(6, 10, 6, 6),
          legend.position = "bottom")
}

hint <- function(...) tags$p(class = "gid-hint", ...)

## One place describing every method: drives the picker, the sidebar sections,
## and the explainer at the top of the Individuals tab.
GID_METHODS <- list(
  probabilistic = list(
    label  = "Likelihood ratio",
    anchor = "m-lr", tag = "recommended",
    blurb  = paste("Asks how much more probable the two observed genotypes are if they",
                   "came from one animal than from two, using an explicit dropout and",
                   "false-allele error model, then converts that evidence into a",
                   "posterior probability using a prior estimated from your data."),
    good   = "You want a probability for every pair and control over how strict the cutoff is.",
    watch  = "You must state what \u201ca different animal\u201d means. For group-living species that is full siblings, not unrelated animals."),
  sethi = list(
    label  = "Sethi et al. (2016)",
    anchor = "m-sethi", tag = "recommended",
    blurb  = paste("Same likelihood, but divides by whichever relationship best explains",
                   "the pair rather than one you nominate, and accepts a match on evidence",
                   "alone (\u039b > 1) with no prior."),
    good   = "You do not want to choose an alternative hypothesis, or to argue about a prior.",
    watch  = "\u039b > 1 takes no account of how many pairs you tested, so it merges more on large datasets."),
  threshold = list(
    label  = "Mismatch threshold",
    anchor = "m-thresh", tag = "gut check",
    blurb  = paste("Calls two samples the same animal when they differ at no more than k",
                   "loci. The sweep on the comparison tab shows how much the answer",
                   "depends on k."),
    good   = "A quick, transparent check that everyone understands.",
    watch  = "With no gap in the mismatch distribution, k decides your answer rather than the data."),
  genalex = list(
    label  = "GenAlEx Matches",
    anchor = "m-genalex", tag = "familiar",
    blurb  = paste("Reproduces GenAlEx's Multilocus \u2192 Matches routine: the match",
                   "distribution, exact matches as individuals, near matches flagged for",
                   "inspection, and P(ID) / P(ID)sib."),
    good   = "Continuity with what your lab already runs, and the near-match list is genuinely useful.",
    watch  = "Identification is exact matching, so a single dropout splits one animal in two."),
  allelematch = list(
    label  = "allelematch",
    anchor = "m-am", tag = "gut check",
    blurb  = paste("Galpern et al. (2012). Scores dissimilarity as the fraction of",
                   "mismatching alleles, clusters, and picks its own threshold by",
                   "minimising samples it cannot classify."),
    good   = "An independent published check that chooses its own settings.",
    watch  = "Built for microsatellites; on dense SNP panels with missing data it tends to over-split."),
  exact = list(
    label  = "Exact match",
    anchor = "m-exact", tag = "baseline",
    blurb  = "Same animal only if identical at every locus where both samples were called.",
    good   = "A baseline. If it agrees with everything else, your genotypes are clean.",
    watch  = "It cannot express doubt, and one dropout splits an animal in two.")
)
GID_METHOD_CHOICES <- setNames(names(GID_METHODS),
                               vapply(GID_METHODS, `[[`, "", "label"))

## value_box()'s first formal is `title`, so the tooltip cannot be passed as a
## `title` attribute -- it would bind to the formal and displace the label.
## Wrap the box instead and put the tooltip on the wrapper.
vbox <- function(label, id, ico, theme = "light", key = NULL) {
  tip <- if (!is.null(key) && !is.null(GID_GLOSSARY[[key]])) GID_GLOSSARY[[key]][[2]]
  box <- value_box(label, textOutput(id, inline = TRUE), showcase = icon(ico),
                   theme = theme, showcase_layout = showcase_left_center(width = "58px"),
                   height = "104px")
  if (is.null(tip)) box
  else tags$div(class = "gid-tip", title = tip,
                style = "height:100%", box)
}

## Tables show plain-English headers; the raw column name and its definition
## are on the header tooltip, so nothing is hidden from anyone who needs it.
dt <- function(x, ..., relabel = TRUE) {
  if (relabel) x <- gid_relabel(x)
  tips <- attr(x, "gid_tips")
  num  <- which(vapply(x, is.numeric, TRUE)) - 1L
  hdr  <- if (!is.null(tips) && any(nzchar(tips)))
    JS(sprintf("function(thead){var t=%s;$(thead).find('th').each(function(i){
         if(t[i]){$(this).attr('title',t[i]).addClass('gid-th-tip');}});}",
       jsonlite::toJSON(tips))) else NULL
  datatable(
    x, rownames = FALSE, extensions = "Buttons",
    options = c(list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
                     buttons = list(list(extend = "csv", text = "Download CSV"))),
                if (length(num)) list(columnDefs = list(
                  list(className = "dt-body-right", targets = num))),
                if (!is.null(hdr)) list(headerCallback = hdr)),
    ...)
}

# ============================================================================ UI
ui <- page_navbar(
  id = "nav",
  fillable = FALSE,   # let explicit plot/card heights stand instead of being flexed
  title = tags$span(tags$strong("genoID"),
                    tags$span(class = "gid-sub", "unique individuals from consensus genotypes")),
  # local = FALSE links the webfont rather than downloading it at startup, and
  # the collection falls back to system fonts if the machine is offline.
  theme = bs_theme(
    version = 5, primary = INK, secondary = MUTED,
    base_font = font_collection(font_google("Source Sans 3", local = FALSE),
                                "system-ui", "Segoe UI", "Helvetica Neue", "sans-serif"),
    heading_font = font_collection(font_google("Source Sans 3", local = FALSE),
                                   "system-ui", "sans-serif"),
    "navbar-bg" = "#ffffff", "body-bg" = "#fbfcfd"),
  header = tags$head(
    katex_head(),
    # Files are built in the browser and handed to a Blob, rather than served
    # from a Shiny download route. Under WebAssembly there is no server to
    # serve that route, and the service-worker fallback fails in Chrome
    # (Chromium bug 468227). A Blob works the same everywhere.
    tags$script(HTML("
      // Registered once only. Shiny appends to its handler-order list every
      // time a name is registered, so registering twice dispatches every
      // message twice and the browser downloads two copies of each file.
      if (!window.__gidDownloadReady) {
      window.__gidDownloadReady = true;
      Shiny.addCustomMessageHandler('gid_download', function(m) {
        var blob;
        if (m.b64) {
          var bin = atob(m.data), arr = new Uint8Array(bin.length);
          for (var i = 0; i < bin.length; i++) arr[i] = bin.charCodeAt(i);
          blob = new Blob([arr], {type: m.type});
        } else {
          var txt = (m.type.indexOf('csv') >= 0 ? '\\ufeff' : '') + m.data;
          blob = new Blob([txt], {type: m.type});
        }
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url; a.download = m.filename; a.style.display = 'none';
        document.body.appendChild(a); a.click();
        setTimeout(function () { URL.revokeObjectURL(url); a.remove(); }, 2000);
      });
      // Plots rendered while their tab is hidden keep the width they had when
      // hidden, so they come back the wrong size. Nudge Shiny to re-measure
      // whenever a tab is shown. It must be a jQuery-triggered resize --
      // Shiny binds through jQuery and ignores a plain native resize event.
      $(document).on('shown.bs.tab', function () {
        setTimeout(function () { $(window).trigger('resize'); }, 60);
      });
      // The Methods tab lays out after the nav switch, so poll briefly for the
      // target section. The page body does not scroll -- bslib scrolls an inner
      // container -- so find that container and move it, rather than relying on
      // scrollIntoView finding the right ancestor.
      Shiny.addCustomMessageHandler('gid_scroll', function(m) {
        var tries = 0;
        (function seek() {
          var el = document.getElementById(m.id);
          if (el && el.getBoundingClientRect().top !== 0) {
            var c = el.parentElement;
            while (c && !(c.scrollHeight > c.clientHeight + 5 &&
                          /(auto|scroll)/.test(getComputedStyle(c).overflowY)))
              c = c.parentElement;
            // Honour prefers-reduced-motion. Smooth scrolling is also silently
            // ignored in some environments, so this doubles as the reliable path.
            var reduce = window.matchMedia &&
                         window.matchMedia('(prefers-reduced-motion: reduce)').matches;
            var how = reduce ? 'auto' : 'smooth';
            if (c) {
              var top = el.getBoundingClientRect().top - c.getBoundingClientRect().top +
                        c.scrollTop - 16;
              c.scrollTo({top: top, behavior: how});
              // if smooth was ignored, land it anyway
              setTimeout(function () {
                if (Math.abs(c.scrollTop - top) > 40) c.scrollTop = top;
              }, 600);
            } else {
              el.scrollIntoView({behavior: how, block: 'start'});
            }
            return;
          }
          if (++tries < 25) setTimeout(seek, 120);
        })();
      });
      }")),
    tags$style(HTML(sprintf("
    .gid-sub{font-weight:400;color:%s;font-size:.78rem;margin-left:.6rem;letter-spacing:.02em}
    .gid-hint{color:%s;font-size:.82rem;margin:.15rem 0 .8rem;line-height:1.45}
    .card{box-shadow:0 1px 3px rgba(29,53,87,.07);border:1px solid #e3e9ef}
    .card-header{background:#fff;border-bottom:1px solid %s;font-weight:600;color:%s;font-size:.92rem}
    .bslib-value-box .value-box-value{font-size:1.85rem!important}
    .bslib-value-box .value-box-title{font-size:.74rem!important;letter-spacing:.05em;text-transform:uppercase}
    .form-label{font-weight:600;font-size:.83rem;color:%s;margin-bottom:.15rem}
    .irs-bar,.irs-handle>i:first-child{background:%s!important;border-color:%s!important}
    .gid-flag{background:#fdf3ee;border-left:3px solid %s;padding:.55rem .8rem;margin:.4rem 0;
              font-size:.85rem;border-radius:0 4px 4px 0}
    .gid-ok{background:#eef5f0;border-left:3px solid #4a7c59}
    .nav-link{font-size:.9rem}
    .value-box-showcase{max-width:64px!important;flex-basis:64px!important;padding:0!important}
    .value-box-showcase svg,.value-box-showcase .fa,.value-box-showcase i{
      height:1.6rem!important;width:1.6rem!important;font-size:1.6rem!important;opacity:.8}
    .bslib-value-box .value-box-area{padding:.35rem .2rem}
    .gid-empty{text-align:center;padding:3.2rem 1rem;color:%s}
    .gid-empty h4{color:%s;font-weight:600;font-size:1.05rem;margin-bottom:.4rem}
    .gid-empty p{font-size:.88rem;max-width:34rem;margin:0 auto .35rem}
    .dataTables_wrapper{font-size:.86rem}
    th.gid-th-tip{cursor:help;border-bottom:1px dotted %s!important}
    th.gid-th-tip:hover{color:%s!important}
    .gid-tip{cursor:help}
    .gid-mhead{background:#fff;border:1px solid #e3e9ef;border-left:3px solid %s;
      border-radius:0 4px 4px 0;padding:.9rem 1.1rem;margin-bottom:1rem;
      box-shadow:0 1px 3px rgba(29,53,87,.06)}
    .gid-mhead h5{font-size:1.02rem;font-weight:600;color:%s;margin:0 0 .35rem}
    .gid-mhead p{font-size:.88rem;color:#4a6067;margin:.3rem 0;line-height:1.5}
    .gid-mtag{display:inline-block;font-size:.62rem;font-weight:700;letter-spacing:.08em;
      text-transform:uppercase;padding:.14em .5em;border-radius:2px;margin-left:.55rem;
      vertical-align:.18em;background:#f6eddc;color:#a8762a}
    .gid-mtag.chk{background:#e3efef;color:#0f6b73}
    .gid-mgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(15rem,1fr));
      gap:.2rem 1.6rem;margin-top:.5rem}
  ", MUTED, MUTED, SOFT, INK, INK, INK, INK, ACCENT, MUTED, INK, ACCENT, INK,
     MUTED, INK)))),

  sidebar = sidebar(
    width = 330, class = "bg-white",

    fileInput("file", "Genotype file", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
    actionLink("load_demo", "or load the example dataset", class = "gid-hint"),
    uiOutput("demo_note"),
    conditionalPanel("output.is_excel == true",
                     selectInput("sheet", "Worksheet", choices = NULL)),

    selectInput("id_col", "Sample ID column", choices = NULL),
    selectizeInput("group_col", "Analyse separately by (optional)", choices = NULL,
                   options = list(placeholder = "e.g. species, study area")),
    hint("Blocking stops samples from different species or populations ",
         "from ever being called the same animal, and keeps allele frequencies separate."),

    conditionalPanel(
      "output.has_reps == true",
      tags$div(class = "gid-flag gid-ok", style = "font-size:.8rem",
        tags$b("This file has several rows per sample."), " They look like PCR ",
        "replicates. Using them directly is better than collapsing them to a ",
        "consensus first \u2014 see Methods."),
      selectInput("rep_col", "Replicate label column", choices = NULL),
      selectizeInput("rep_drop", "Row labels to exclude", choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "e.g. consensus")),
      hint("Exclude any summary rows, such as a pre-computed consensus. ",
           "Everything left is treated as one PCR replicate."),
      radioButtons("rep_mode", "Genotypes to analyse",
                   c("Use every replicate observation (recommended)" = "reps",
                     "Collapse to a consensus first" = "consensus"),
                   selected = "reps"),
      uiOutput("rep_note")),

    selectizeInput("use_loci", "Loci to use", choices = NULL, multiple = TRUE,
                   options = list(plugins = list("remove_button"), placeholder = "detected automatically")),
    uiOutput("locus_note"),
    hint("Detected automatically, then yours to correct. Remove species-diagnostic ",
         "and sex markers: they are near-fixed within a species, so they add no ",
         "power and distort the frequency model."),
    actionLink("reset_loci", "Reset to detected", class = "gid-hint"),

    tags$hr(style = "margin:.9rem 0 .8rem"),

    selectInput("method", "Method for the Individuals tab",
                choices = GID_METHOD_CHOICES, selected = "probabilistic"),
    uiOutput("method_hint"),

    ## Settings for the chosen method. The other methods keep their values while
    ## hidden, so the comparison tab can still run all six.
    conditionalPanel(
      "input.method == 'probabilistic'",
      selectInput("kinship", "Compare against", selected = "full_sib",
                  choices = c("Full siblings (conservative)" = "full_sib",
                              "Half siblings" = "half_sib",
                              "Unrelated individuals" = "unrelated")),
      hint("The alternative hypothesis. For pack- or group-living species the ",
           "samples competing to be \"a different animal\" are usually relatives, ",
           "so full siblings is the honest null."),
      sliderInput("post_cut", "Posterior probability to accept a match",
                  0.9, 0.9999, 0.999, step = 0.0001),
      numericInput("prior", "Prior P(a random pair is a recapture)", NA,
                   min = 0, max = 0.5, step = 0.005),
      hint("Leave blank to estimate it from the data.")),

    conditionalPanel(
      "input.method == 'sethi'",
      checkboxGroupInput("sethi_rel", "Competing relationships",
                         choices = c("Unrelated" = "unrelated",
                                     "Full siblings" = "full_sib",
                                     "Parent-offspring" = "parent_offspring",
                                     "Half siblings" = "half_sib"),
                         selected = c("unrelated", "full_sib", "parent_offspring")),
      hint("Divides by whichever of these best explains the pair, so you do not ",
           "have to pick one. A match must beat every relationship on the list."),
      numericInput("lambda_cut", "Accept a match when \u039b exceeds", 1, 0.01, 1e6, 1),
      hint("The paper uses \u039b > 1. Raise it on large datasets.")),

    conditionalPanel(
      "input.method == 'threshold'",
      numericInput("max_mismatch", "Loci allowed to differ", 1, 0, 20, 1),
      hint("Set this inside the gap in the mismatch distribution shown on the ",
           "Individuals tab.")),

    conditionalPanel(
      "input.method == 'genalex'",
      numericInput("near_match", "Flag near matches within N loci", 2, 1, 10, 1),
      hint("Near matches are listed on the comparison tab for you to inspect ",
           "rather than merged automatically.")),

    conditionalPanel(
      "input.method == 'allelematch'",
      hint("allelematch selects its own threshold from your data. Nothing to set.")),

    conditionalPanel(
      "input.method == 'exact'",
      hint("No tolerance to set \u2014 that is the point of this baseline.")),

    ## error model: used by both likelihood methods
    conditionalPanel(
      "input.method == 'probabilistic' || input.method == 'sethi'",
      sliderInput("dropout", "Allelic dropout rate", 0, 0.20, 0.005, step = 0.001),
      sliderInput("false_allele", "False allele rate", 0, 0.10, 0.002, step = 0.001),
      hint("Rates of the genotypes you are uploading. Multi-replicate consensus ",
           "calls are usually well under 0.005.")),

    accordion(
      open = FALSE,
      accordion_panel(
        "Settings shared by every method", value = "shared", icon = icon("sliders"),
        numericInput("min_loci", "Minimum loci compared per pair", 15, 1, 500, 1),
        numericInput("min_sample_call", "Minimum sample call rate", 0.5, 0, 1, 0.05),
        numericInput("min_locus_call", "Minimum locus call rate", 0.25, 0, 1, 0.05),
        radioButtons("linkage", "Cluster rule",
                     c("Single linkage (transitive)" = "single",
                       "Complete linkage (all pairs must match)" = "complete"),
                     selected = "single"),
        hint("Single linkage lets A-B and B-C merge A, B and C. Complete linkage ",
             "requires every pair inside a cluster to match. Compare both: if the ",
             "answer moves, some clusters are held together by one edge.")),
      accordion_panel(
        "Which methods to run", value = "which", icon = icon("list-check"),
        hint("The comparison tab runs everything ticked here. Untick the slow ones ",
             "if you only care about your chosen method."),
        checkboxInput("run_sethi", "Sethi et al. (2016)", TRUE),
        checkboxInput("run_genalex", "GenAlEx Matches", TRUE),
        checkboxInput("run_am", "allelematch", TRUE))
    ),

    actionButton("run", "Identify individuals", icon = icon("play"),
                 class = "btn-primary w-100"),
    tags$div(class = "gid-hint", style = "margin-top:.6rem",
             "Nothing leaves this machine unless you deploy the app to a server.")
  ),

  # ------------------------------------------------------------------ panels
  nav_panel(
    "Data & QC", icon = icon("table"),
    uiOutput("empty_state"),
    conditionalPanel("output.has_data == true",
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      vbox("Samples", "vb_samples", "vial", "primary", key = "n_samples"),
      vbox("Loci used", "vb_loci", "dna", "secondary", key = "n_loci"),
      vbox("Median call rate", "vb_call", "percent", "light", key = "call_rate"),
      vbox("Groups", "vb_groups", "layer-group", "light", key = "group")
    ),
    layout_columns(
      col_widths = c(12),
      card(card_header("Data checks"), uiOutput("flags"))
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("Per-locus summary"), DTOutput("tbl_locus")),
      card(card_header("Per-sample summary"), DTOutput("tbl_sample"))
    ),
    card(card_header("Linkage between loci (r², top 25 pairs)"),
         hint("P(ID) multiplies across loci and assumes they are independent. ",
              "Two SNPs from one amplicon are not. Drop one of any pair with high r² ",
              "before quoting panel power."),
         DTOutput("tbl_ld")))
  ),

  nav_panel(
    "Panel power", icon = icon("chart-line"),
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("Cumulative probability of identity"), plotOutput("plot_pid", height = 380)),
      card(card_header("Per-sample resolving power"), plotOutput("plot_power", height = 380))
    ),
    layout_columns(
      col_widths = c(12),
      card(card_header("P(ID) by locus"),
           hint("P(ID) is the chance two unrelated animals share a genotype; ",
                "P(ID)sib is the same for full siblings and is the number to quote. ",
                "A panel is normally called adequate below 0.01, and comfortable below 0.001."),
           DTOutput("tbl_pid"))
    )
  ),

  nav_panel(
    "Individuals", icon = icon("fingerprint"),
    uiOutput("method_header"),
    uiOutput("run_status_ind"),
    uiOutput("power_warning"),
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      vbox("Individuals", "vb_ind", "paw", "primary", key = "n_individuals"),
      vbox("Sampled once", "vb_single", "circle-dot", "light", key = "n_singletons"),
      vbox("Recapture rate", "vb_recap", "repeat", "light", key = "recapture_rate"),
      vbox("Largest cluster", "vb_max", "maximize", "light", key = "max_cluster")
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header(textOutput("evid_title", inline = TRUE)),
           plotOutput("plot_lr", height = 340),
           uiOutput("evid_note")),
      card(card_header("Locus mismatches between samples"), plotOutput("plot_mm", height = 340),
           hint("A clean panel gives two piles with a gap. Pairs sitting in the gap ",
                "are the ones worth looking at by hand."))
    ),
    card(card_header("Individual assignments"), DTOutput("tbl_final")),
    card(card_header("Clusters that are not internally consistent"),
         hint("Every sample in a cluster should match every other sample in it. ",
              "Where that fails, the cluster is being held together by a chain and ",
              "may be two animals. Inspect these before reporting a population size."),
         DTOutput("tbl_conflict"))
  ),

  nav_panel(
    "Method comparison", icon = icon("code-compare"),
    uiOutput("run_status_cmp"),
    card(card_header("How each method resolved the same data"), DTOutput("tbl_cmp")),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Agreement between methods (adjusted Rand index)"),
           hint("1 = identical partitions. Values below about 0.95 mean the choice ",
                "of method is materially changing your answer."),
           DTOutput("tbl_ari")),
      card(card_header("Sensitivity to the mismatch threshold"), plotOutput("plot_sweep", height = 300),
           hint("A flat stretch means the answer is robust. A steady slide means the ",
                "threshold, not the data, is deciding how many animals you have."))
    ),
    conditionalPanel("input.run_sethi == true",
      card(card_header("Sethi et al. \u2014 which relationship competed?"),
        hint("For every pair called a match, the alternative hypothesis that came closest. ",
             "If one relationship dominates, that is the null your other analyses should ",
             "be using too."),
        uiOutput("sethi_alt"))),
    conditionalPanel("input.run_genalex == true",
      card(card_header("GenAlEx Matches output"),
        hint("The same three things GenAlEx's Multilocus \u2192 Matches routine reports: ",
             "how many pairs differ at each number of loci, the near matches worth ",
             "inspecting by hand, and the probability of identity for the panel."),
        layout_columns(col_widths = c(5, 7),
          tagList(tags$div(class = "gid-hint", tags$b("Match distribution")),
                  DTOutput("tbl_gx_dist")),
          tagList(tags$div(class = "gid-hint", tags$b("Near matches \u2014 inspect these by hand")),
                  DTOutput("tbl_gx_near"))),
        uiOutput("gx_pid"))),
    card(card_header("Samples the methods disagree about"),
         hint("These are the samples worth a second look, and the ones to mention ",
              "in a methods section."),
         DTOutput("tbl_disagree")),
    card(card_header("How many individuals under other settings"),
         hint("The same data, re-run across the settings that matter most. If the ",
              "count barely moves, the number is a property of your data. If it ",
              "swings, it is a property of your choices, and the range belongs in ",
              "your results. Your current setting is marked."),
         plotOutput("plot_sens", height = 330),
         DTOutput("tbl_sens"))
  ),

  methods_panel(),

  nav_panel(
    "Download", icon = icon("download"),
    card(
      card_header("Results"),
      tags$p(class = "gid-hint", style = "font-size:.9rem",
             "Every table in the app also has its own Download CSV button."),
      layout_columns(
        col_widths = c(4, 4, 4),
        actionButton("dl_assign", "Individual assignments",
                     icon = icon("download"), class = "btn-primary w-100"),
        actionButton("dl_cons", "Consensus genotype per individual",
                     icon = icon("download"), class = "btn-outline-primary w-100"),
        actionButton("dl_pairs", "All pairwise comparisons",
                     icon = icon("download"), class = "btn-outline-primary w-100")
      ),
      tags$hr(),
      layout_columns(
        col_widths = c(4, 4, 4),
        actionButton("dl_all", "Everything (zip)",
                     icon = icon("file-zipper"), class = "btn-primary w-100"),
        actionButton("dl_params", "Settings used (for your methods section)",
                     icon = icon("download"), class = "btn-outline-primary w-100"),
        actionButton("dl_script", "Equivalent R script",
                     icon = icon("download"), class = "btn-outline-primary w-100")
      )
    ),
    card(card_header("Settings used"), verbatimTextOutput("params_txt"))
  ),

  nav_spacer(),
  nav_item(tags$a(href = "#", onclick = "return false;", class = "gid-sub",
                  "P(ID) after Waits et al. 2001 · allelematch after Galpern et al. 2012"))
)

# ======================================================================== SERVER
server <- function(input, output, session) {

  raw <- reactiveVal(NULL)

  ## GENOID_PRELOAD=/path/to/file.csv opens a file at startup. Handy when a lab
  ## always analyses the same export, and it is how the app is smoke-tested.
  observe({
    pre <- Sys.getenv("GENOID_PRELOAD", "")
    if (nzchar(pre) && file.exists(pre) && is.null(raw())) raw(gid_read(pre))
  })

  demo_loaded <- reactiveVal(FALSE)
  observeEvent(input$load_demo, {
    f <- c("demo/demo_genotypes.csv", "app/demo/demo_genotypes.csv")
    f <- f[file.exists(f)]
    if (!length(f)) return(showNotification("Example file not found.", type = "error"))
    raw(gid_read(f[1])); demo_loaded(TRUE)
    # the example was simulated at these rates; start the user at the truth so
    # the demo shows the method working, then let them break it on purpose
    updateSliderInput(session, "dropout", value = 0.03)
    updateSliderInput(session, "false_allele", value = 0.01)
    updateNumericInput(session, "min_loci", value = 12)
  })
  observeEvent(input$file, demo_loaded(FALSE))

  output$demo_note <- renderUI({
    if (!isTRUE(demo_loaded())) return(NULL)
    tags$div(class = "gid-flag", style = "margin:.5rem 0 .2rem;font-size:.8rem",
      tags$b("Example data loaded."), " 55 simulated samples from 28 known individuals in ",
      "two species. The true answer is in the ", tags$code("TrueIndividual"), " column, so ",
      "you can check any method against it. Simulated with 3% dropout and 1% false alleles \u2014 ",
      "the sliders have been set to match. Lower them and re-run to see what happens when you ",
      "understate your error rate.")
  })

  output$has_data <- reactive(!is.null(raw()))
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)

  output$empty_state <- renderUI({
    if (!is.null(raw())) return(NULL)
    tags$div(class = "gid-empty",
      icon("dna", class = "fa-2x", style = paste0("color:", SOFT, ";margin-bottom:.8rem")),
      tags$h4("Upload a genotype table to begin"),
      tags$p("One row per sample, one column per locus. Genotypes as two alleles ",
             "in a single cell (AG) or separated (120/124). Missing data as 00, ",
             "NA or blank. CSV, TSV or Excel."),
      tags$p("Loci are detected automatically and shown in the sidebar for you to correct."))
  })

  is_excel <- reactive(!is.null(input$file) &&
                         grepl("xlsx?$", tolower(input$file$name)))
  output$is_excel <- reactive(isTRUE(is_excel()))
  outputOptions(output, "is_excel", suspendWhenHidden = FALSE)

  observeEvent(input$file, {
    if (is_excel()) {
      sh <- readxl::excel_sheets(input$file$datapath)
      updateSelectInput(session, "sheet", choices = sh, selected = sh[1])
    } else raw(gid_read(input$file$datapath))
  })

  observeEvent(input$sheet, {
    req(is_excel(), nzchar(input$sheet))
    raw(gid_read(input$file$datapath, sheet = input$sheet))
  })

  ## A file with repeated sample IDs is either duplicated samples or PCR
  ## replicates. Offer the replicate route when it looks like the latter.
  has_reps <- reactive({
    df <- raw(); if (is.null(df) || is.null(input$id_col)) return(FALSE)
    anyDuplicated(as.character(df[[input$id_col]])) > 0
  })
  output$has_reps <- reactive(isTRUE(has_reps()))
  outputOptions(output, "has_reps", suspendWhenHidden = FALSE)

  observeEvent(list(raw(), input$id_col), {
    df <- raw(); req(df, input$id_col)
    if (!has_reps()) return()
    nonloci <- setdiff(names(df), gid_detect_loci(df, sep = gid_guess_sep(df)))
    guess <- grep("^(rep|replicate|pcr|run)", nonloci, ignore.case = TRUE, value = TRUE)
    sel <- if (length(guess)) guess[1] else nonloci[1]
    updateSelectInput(session, "rep_col", choices = nonloci, selected = sel)
  })

  observeEvent(input$rep_col, {
    df <- raw(); req(df, nzchar(input$rep_col %||% ""))
    v <- sort(unique(as.character(df[[input$rep_col]])))
    ## anything that looks like a summary row, not a reaction
    auto <- v[grepl("consensus|final|combined|summary", v, ignore.case = TRUE)]
    updateSelectizeInput(session, "rep_drop", choices = v, selected = auto)
  })

  output$rep_note <- renderUI({
    if (!isTRUE(has_reps())) return(NULL)
    if (identical(input$rep_mode, "reps"))
      tags$p(class = "gid-hint",
        tags$b("The dropout and false-allele rates below now mean per-reaction "),
        tags$b("rates"), ", not the residual rates of a consensus call. They are ",
        "typically ten to a hundred times larger \u2014 a few percent rather than ",
        "a few tenths of a percent.")
    else
      tags$p(class = "gid-hint",
        "Replicates are collapsed with the multi-tube rule: a heterozygote ",
        "accepted on two replicates, a homozygote only on three.")
  })

  detected <- reactive({
    df <- raw(); req(df)
    gid_detect_loci(df, exclude = input$id_col, sep = gid_guess_sep(df, exclude = input$id_col))
  })

  ## once data is in, guess the ID column and the loci
  observeEvent(raw(), {
    df <- raw(); req(df)
    d <- gid_detect_loci(df, sep = gid_guess_sep(df))
    nonloci <- setdiff(names(df), d)
    idg <- gid_guess_id_col(df, exclude = d)
    updateSelectInput(session, "id_col", choices = names(df), selected = idg)
    updateSelectizeInput(session, "group_col", choices = c("(none)" = "", nonloci), selected = "")
    # every column is offerable, so nothing the detector missed is unreachable
    updateSelectizeInput(session, "use_loci", choices = names(df), selected = d)
  })

  observeEvent(input$reset_loci,
               updateSelectizeInput(session, "use_loci", selected = detected()))

  output$locus_note <- renderUI({
    df <- raw(); req(df)
    d <- detected(); u <- input$use_loci %||% character(0)
    added <- setdiff(u, d); removed <- setdiff(d, u)
    txt <- sprintf("%d of %d detected loci in use", length(intersect(u, d)), length(d))
    if (length(removed)) txt <- paste0(txt, "; removed ", paste(removed, collapse = ", "))
    if (length(added))   txt <- paste0(txt, "; added ", paste(added, collapse = ", "))
    tags$p(class = "gid-hint", style = "margin-top:-.5rem", txt)
  })

  loci <- reactive({
    df <- raw(); req(df, input$id_col)
    setdiff(input$use_loci %||% detected(), input$id_col)
  })

  ## sample x locus matrix, filtered
  prep <- reactive({
    df <- raw(); req(df, input$id_col, length(loci()) > 0)
    ids <- as.character(df[[input$id_col]])
    sep <- gid_guess_sep(df, exclude = input$id_col)
    use_reps <- isTRUE(has_reps()) && identical(input$rep_mode, "reps") &&
                nzchar(input$rep_col %||% "")

    ## ---- replicate route: several rows per sample, one genotype per sample
    if (isTRUE(has_reps()) && nzchar(input$rep_col %||% "")) {
      lab  <- as.character(df[[input$rep_col]])
      keep <- !lab %in% (input$rep_drop %||% character(0))
      validate(need(sum(keep) > 0,
        "Every row was excluded. Clear some entries from 'Row labels to exclude'."))
      rgt  <- gid_matrix(df[keep, , drop = FALSE], input$id_col, loci(), sep = sep)
      rsm  <- ids[keep]
      gt   <- gid_consensus_from_reps(rgt, rsm)
      reps <- if (use_reps) list(gt = rgt, sample = rsm) else NULL
      grp  <- if (nzchar(input$group_col %||% ""))
        as.character(df[[input$group_col]])[keep][match(rownames(gt), rsm)]
        else rep("all", nrow(gt))
      names(grp) <- rownames(gt)
      grp[is.na(grp) | !nzchar(grp)] <- "unassigned"
      f <- gid_filter(gt,
             min_locus_call  = if (is.na(input$min_locus_call))  0 else input$min_locus_call,
             min_sample_call = if (is.na(input$min_sample_call)) 0 else input$min_sample_call)
      validate(need(nrow(f$gt) >= 2,
        "Fewer than two samples survive the call-rate filters. Lower them in the sidebar."),
        need(ncol(f$gt) >= 1, "No loci survive the call-rate filter."))
      if (!is.null(reps))
        reps$gt <- reps$gt[, colnames(f$gt), drop = FALSE]
      return(list(gt = f$gt, grp = grp[rownames(f$gt)], raw_gt = gt, reps = reps,
                  dropped_loci = f$dropped_loci, dropped_samples = f$dropped_samples,
                  df = df, n_reps = length(rsm) / nrow(gt)))
    }

    if (anyDuplicated(ids))
      ids <- ave(ids, ids, FUN = function(z)
        if (length(z) == 1) z else paste0(z, "#", seq_along(z)))
    df$.gid_key <- ids
    gt <- gid_matrix(df, ".gid_key", loci(), sep = sep)
    grp <- if (nzchar(input$group_col %||% "")) as.character(df[[input$group_col]]) else rep("all", nrow(df))
    names(grp) <- ids
    grp[is.na(grp) | !nzchar(grp)] <- "unassigned"

    f <- gid_filter(gt,
                    min_locus_call  = if (is.na(input$min_locus_call))  0 else input$min_locus_call,
                    min_sample_call = if (is.na(input$min_sample_call)) 0 else input$min_sample_call)
    validate(need(nrow(f$gt) >= 2,
                  "Fewer than two samples survive the call-rate filters. Lower them in the sidebar."),
             need(ncol(f$gt) >= 1,
                  "No loci survive the call-rate filter. Lower it in the sidebar."))
    list(gt = f$gt, grp = grp[rownames(f$gt)], raw_gt = gt, reps = NULL,
         dropped_loci = f$dropped_loci, dropped_samples = f$dropped_samples, df = df)
  })

  ## ---------------------------------------------------------------- QC panel
  output$vb_samples <- renderText({ p <- prep(); req(p); nrow(p$gt) })
  output$vb_loci    <- renderText({ p <- prep(); req(p); ncol(p$gt) })
  output$vb_call    <- renderText({ p <- prep(); req(p); sprintf("%.0f%%", 100 * median(rowMeans(!is.na(p$gt)))) })
  output$vb_groups  <- renderText({ p <- prep(); req(p); length(unique(p$grp)) })

  output$flags <- renderUI({
    p <- prep(); req(p)
    f <- list()
    ok <- function(x) tags$div(class = "gid-flag gid-ok", x)
    bad <- function(x) tags$div(class = "gid-flag", x)

    nre <- sum(gsub("/", "", p$raw_gt) !=
                 toupper(gsub("[*?!#]", "", as.matrix(p$df[, colnames(p$raw_gt)]))), na.rm = TRUE)
    if (nre > 0) f <- c(f, list(bad(sprintf(
      "%d genotype cells had their two alleles written in the other order (AG vs GA). They are now normalised; compared as raw text they would have counted as mismatches.", nre))))

    nflag <- sum(grepl("[*?!#]", as.matrix(p$df[, colnames(p$raw_gt)])))
    if (nflag > 0) f <- c(f, list(bad(sprintf(
      "%d cells carry a quality flag character. The flag was stripped and the genotype kept. If your pipeline uses the flag to mean \"low confidence\", consider setting those cells to missing before upload.", nflag))))

    dupn <- sum(duplicated(as.character(p$df[[input$id_col]])))
    if (dupn > 0) f <- c(f, list(bad(sprintf(
      "%d duplicated sample IDs. They were kept and suffixed with #1, #2 -- if these are the same extract run twice they are a useful positive control, since any correct method must put them together.", dupn))))

    if (length(p$dropped_loci)) f <- c(f, list(bad(paste0(
      "Loci dropped for low call rate or being monomorphic: ",
      paste(p$dropped_loci, collapse = ", ")))))
    if (length(p$dropped_samples)) f <- c(f, list(bad(paste0(
      "Samples dropped for low call rate: ", paste(p$dropped_samples, collapse = ", ")))))

    ls_ <- gid_locus_stats(p$gt)
    mono <- ls_$locus[ls_$maf < 0.05 & !is.na(ls_$maf)]
    if (length(mono)) f <- c(f, list(bad(paste0(
      "Nearly monomorphic loci (minor allele under 5%), contributing almost nothing: ",
      paste(mono, collapse = ", ")))))

    if (!length(f)) f <- list(ok("No structural problems found in the uploaded table."))
    tagList(f)
  })

  output$tbl_locus <- renderDT({ p <- prep(); req(p)
    x <- gid_locus_stats(p$gt)
    dt(transform(x, call_rate = round(call_rate, 3), maf = round(maf, 3),
                 Ho = round(Ho, 3), He = round(He, 3), Fis = round(Fis, 3))) })

  output$tbl_sample <- renderDT({ p <- prep(); req(p)
    x <- gid_sample_stats(p$gt); x$group <- p$grp[x$sample]
    dt(transform(x, call_rate = round(call_rate, 3))) })

  output$tbl_ld <- renderDT({ p <- prep(); req(p)
    x <- gid_ld(p$gt); req(!is.null(x))
    dt(head(transform(x, r2 = round(r2, 3)), 25)) })

  ## ------------------------------------------------------------ power panel
  freqs <- reactive({ p <- prep(); req(p); gid_allele_freq(p$gt) })
  pid   <- reactive({ x <- gid_pid(freqs()); x <- x[order(x$pid), ]
                      x$pid_cum <- cumprod(x$pid); x$pid_sib_cum <- cumprod(x$pid_sib); x })

  output$plot_pid <- renderPlot({
    x <- pid()
    d <- rbind(data.frame(k = seq_len(nrow(x)), v = x$pid_cum,     s = "P(ID)"),
               data.frame(k = seq_len(nrow(x)), v = x$pid_sib_cum, s = "P(ID) siblings"))
    ggplot(d, aes(k, v, colour = s)) +
      geom_hline(yintercept = c(0.01, 0.001), linetype = 2, colour = MUTED, linewidth = 0.35) +
      geom_line(linewidth = 1) + geom_point(size = 1.6) +
      scale_y_log10(labels = function(z) format(z, scientific = TRUE)) +
      scale_colour_manual(values = PAL[1:2]) +
      annotate("text", x = 1, y = 0.011, label = "0.01", hjust = 0, size = 3, colour = MUTED) +
      labs(x = "Loci (most informative first)", y = "Cumulative probability",
           colour = NULL, subtitle = "Two animals with this many loci in common would share a genotype this often") +
      theme_gid()
  })

  power <- reactive({ p <- prep(); req(p)
    x <- gid_sample_power(p$gt, pid()); x$call_rate <- rowMeans(!is.na(p$gt))[x$sample]; x })

  output$plot_power <- renderPlot({
    x <- power()
    ggplot(x, aes(call_rate, pid_sib)) +
      geom_hline(yintercept = 0.01, linetype = 2, colour = ACCENT) +
      geom_point(colour = INK, alpha = 0.7, size = 2.4) +
      scale_y_log10() +
      labs(x = "Call rate of this sample", y = "P(ID) siblings for this sample",
           subtitle = "Computed from the loci each sample actually has.\nPoints above the line cannot be told apart from a sibling.") +
      theme_gid()
  })

  output$tbl_pid <- renderDT(dt(transform(pid(), pid = signif(pid, 4), pid_sib = signif(pid_sib, 4),
                                          pid_cum = signif(pid_cum, 3), pid_sib_cum = signif(pid_sib_cum, 3))))

  ## ------------------------------------------------------------- run methods
  run_error <- reactiveVal(NULL)

  res <- eventReactive(input$run, {
    p <- prep(); req(p, nrow(p$gt) >= 2)
    run_error(NULL)
    out <- tryCatch(withProgress(message = "Identifying individuals", value = 0, {
      pr <- if (is.na(input$prior)) NULL else input$prior
      incProgress(0.15, detail = "probabilistic")
      lr <- gid_by_group(p$gt, p$grp, gid_method_lr, dropout = input$dropout,
                         false_allele = input$false_allele, kinship = input$kinship,
                         post_cut = input$post_cut, min_loci = input$min_loci,
                         linkage = input$linkage, prior_same = pr, reps = p$reps)
      incProgress(0.3, detail = "exact match")
      ex <- gid_by_group(p$gt, p$grp, gid_method_exact,
                         min_loci = input$min_loci, linkage = input$linkage)
      incProgress(0.2, detail = "mismatch threshold")
      th <- gid_by_group(p$gt, p$grp, gid_method_threshold,
                         max_mismatch = input$max_mismatch, min_loci = input$min_loci,
                         linkage = input$linkage)
      out <- list(exact = ex, threshold = th, probabilistic = lr)
      if (isTRUE(input$run_sethi) && length(input$sethi_rel)) {
        incProgress(0.15, detail = "Sethi et al. match calling")
        out$sethi <- gid_by_group(p$gt, p$grp, gid_method_sethi,
                                  dropout = input$dropout, false_allele = input$false_allele,
                                  relationships = input$sethi_rel,
                                  lambda_cut = input$lambda_cut,
                                  min_loci = input$min_loci, linkage = input$linkage,
                                  reps = p$reps)
      }
      if (isTRUE(input$run_genalex)) {
        incProgress(0.1, detail = "GenAlEx Matches")
        out$genalex <- gid_by_group(p$gt, p$grp, gid_method_genalex,
                                    min_loci = input$min_loci,
                                    near_match_loci = input$near_match,
                                    linkage = input$linkage)
      }
      if (isTRUE(input$run_am)) {
        incProgress(0.25, detail = "allelematch")
        out$allelematch <- try(gid_by_group(p$gt, p$grp, gid_method_allelematch,
                                            min_loci = input$min_loci), silent = TRUE)
        if (inherits(out$allelematch, "try-error")) out$allelematch <- NULL
      }
      incProgress(0.1, detail = "comparing")
      ord <- c("exact", "threshold", "genalex", "allelematch", "sethi", "probabilistic")
      list(methods = out[ord[ord %in% names(out)]],
           cmp = gid_compare_methods(out), gt = p$gt, grp = p$grp,
           sweep = gid_threshold_sweep(p$gt, max_k = 6, min_loci = input$min_loci,
                                       linkage = input$linkage))
    }), error = function(e) { run_error(conditionMessage(e)); NULL })
    out
  })

  ## the method the Individuals tab is showing
  best <- reactive({
    r <- res(); req(r)
    m <- r$methods[[input$method]]
    validate(need(!is.null(m), sprintf(
      "%s was not run. Tick it under \"Which methods to run\" in the sidebar, then press Identify individuals.",
      GID_METHODS[[input$method]]$label)))
    m
  })

  meta <- reactive(GID_METHODS[[input$method %||% "probabilistic"]])

  output$method_hint <- renderUI(
    tags$p(class = "gid-hint", style = "margin-top:-.4rem", meta()$blurb))

  output$method_header <- renderUI({
    m <- meta()
    tags$div(class = "gid-mhead",
      tags$h5(m$label, tags$span(
        class = paste("gid-mtag", if (m$tag != "recommended") "chk" else ""), m$tag)),
      tags$p(m$blurb),
      tags$div(class = "gid-mgrid",
        tags$p(tags$b("Good for: "), m$good),
        tags$p(tags$b("Watch out: "), m$watch)),
      tags$div(style = "margin-top:.7rem",
        actionButton("go_methods", "Read the full method and its equations",
                     icon = icon("book-open"), class = "btn-sm btn-outline-primary")))
  })

  ## jump to the Methods tab, scrolled to this method's section
  observeEvent(input$go_methods, {
    nav_select("nav", "Methods")
    session$sendCustomMessage("gid_scroll", list(id = meta()$anchor))
  })

  output$vb_ind <- renderText({ length(unique(best()$assignment$individual)) })
  output$vb_single <- renderText({ sum(table(best()$assignment$individual) == 1) })
  output$vb_recap <- renderText({
    s <- gid_summarise_assignment(best()$assignment); sprintf("%.0f%%", 100 * s$recapture_rate) })
  output$vb_max <- renderText({ max(table(best()$assignment$individual)) })

  run_status <- function() {
    if (!is.null(run_error()))
      return(tags$div(class = "gid-flag", style = "margin-bottom:1rem",
        tags$b("The analysis could not finish. "),
        "R reported: ", tags$code(run_error()),
        tags$div(style = "margin-top:.4rem",
          "Common causes: too few loci surviving the call-rate filters, or a ",
          "grouping column with a category containing a single sample. Adjust the ",
          "settings in the sidebar and run again.")))
    if (input$run == 0)
      return(tags$div(class = "gid-flag gid-ok", style = "margin-bottom:1rem",
        tags$b("No results yet. "),
        "Choose your sample ID column and loci in the sidebar, then press ",
        tags$b("Identify individuals"), " to compute. Nothing on this tab is ",
        "calculated until you do \u2014 the other tabs update on their own."))
    NULL
  }
  output$run_status_ind <- renderUI(run_status())
  output$run_status_cmp <- renderUI(run_status())

  output$power_warning <- renderUI({
    if (!identical(input$method, "probabilistic")) return(NULL)
    r <- res(); req(r)
    mx <- vapply(r$methods$probabilistic$by_group,
                 function(e) e$max_posterior_observed %||% NA_real_, 0)
    mxlr <- vapply(r$methods$probabilistic$by_group,
                   function(e) e$max_log10_LR_observed %||% NA_real_, 0)
    mx <- mx[!is.na(mx)]; mxlr <- mxlr[!is.na(mxlr)]
    if (!length(mx) || max(mx) >= input$post_cut) return(NULL)
    tags$div(class = "gid-flag", style = "margin-bottom:1rem",
      tags$strong("Your acceptance threshold is beyond what this panel can reach. "),
      sprintf(paste0("The most similar pair in this dataset only reaches a posterior of %s ",
                     "(log10 LR %.2f) against a %s alternative -- %d loci carry a finite ",
                     "amount of evidence. At a cutoff of %s nothing can match, so every ",
                     "sample is being reported as its own individual. Lower the cutoff, add ",
                     "loci, or compare against a less demanding alternative."),
              signif(max(mx), 5), max(mxlr), gsub("_", " ", input$kinship),
              ncol(r$gt), input$post_cut))
  })

  ## What counts as "evidence" depends on the method: a likelihood ratio for
  ## the two probabilistic ones, and mismatching loci for the rest.
  evidence <- reactive({
    r <- res(); req(r)
    b <- best(); pp <- b$pairs
    ## allelematch does its own clustering internally and never exposes a
    ## pairwise table, so fall back to the plain locus-mismatch counts, which
    ## are a property of the data rather than of any method.
    borrowed <- FALSE
    if (is.null(pp) || !nrow(pp)) {
      pp <- r$methods$exact$pairs
      borrowed <- TRUE
    }
    if (is.null(pp) || !nrow(pp)) return(NULL)
    pp <- pp[pp$n_compared >= input$min_loci, ]
    if (borrowed)
      return(list(v = pp$n_mismatch, cut = NA_real_, x = "Loci that differ",
                  kind = "mm",
                  note = paste("allelematch clusters internally and does not expose a",
                               "score per pair, so this shows how many loci differ",
                               "between each pair of samples \u2014 a property of your",
                               "data, not of the method.")))
    if (!is.null(pp$log10_LR))
      list(v = pp$log10_LR, cut = suppressWarnings(min(pp$log10_LR[pp$posterior_same >= input$post_cut])),
           x = "log10 likelihood ratio", kind = "lr",
           note = sprintf("Positive values favour one animal. Alternative hypothesis: %s.",
                          gsub("_", " ", input$kinship)))
    else if (!is.null(pp$log10_lambda))
      list(v = pp$log10_lambda, cut = log10(input$lambda_cut), x = "log10 lambda", kind = "lr",
           note = "Positive values favour one animal over the best competing relationship.")
    else
      list(v = pp$n_mismatch, cut = input$max_mismatch + 0.5, x = "Loci that differ",
           kind = "mm",
           note = "This method decides by counting differences, so the evidence is the mismatch count itself.")
  })

  output$evid_title <- renderText({
    e <- evidence(); if (is.null(e)) "Evidence for every pair of samples"
    else if (e$kind == "lr") "Evidence for every pair of samples"
    else "Loci that differ between samples" })

  output$evid_note <- renderUI({
    e <- evidence(); req(e)
    tags$p(class = "gid-hint", e$note,
           if (is.finite(e$cut)) " The dashed line is where this method accepts a match." else "") })

  output$plot_lr <- renderPlot({
    e <- evidence(); req(e)
    d <- data.frame(v = e$v)
    g <- ggplot(d, aes(v)) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
      labs(x = e$x, y = "Sample pairs", subtitle = "Counts are log-scaled.") +
      theme_gid()
    g <- g + if (e$kind == "mm")
      geom_histogram(binwidth = 1, fill = INK, colour = "white", linewidth = 0.3)
    else geom_histogram(bins = 60, fill = INK, colour = NA)
    if (is.finite(e$cut)) g <- g + geom_vline(xintercept = e$cut, colour = ACCENT, linetype = 2)
    g
  })

  output$plot_mm <- renderPlot({
    r <- res(); pp <- r$methods$threshold$pairs
    pp <- pp[pp$n_compared >= input$min_loci, ]
    ggplot(pp, aes(n_mismatch)) +
      geom_histogram(binwidth = 1, fill = INK, colour = "white", linewidth = 0.3) +
      geom_vline(xintercept = input$max_mismatch + 0.5, colour = ACCENT, linetype = 2) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
      labs(x = "Loci that differ between the two samples", y = "Sample pairs",
           subtitle = "Line = mismatch threshold. Counts are log-scaled.") +
      theme_gid()
  })

  final_tbl <- reactive({
    r <- res(); a <- best()$assignment
    a$group <- r$grp[a$sample]
    a$call_rate <- round(rowMeans(!is.na(r$gt))[a$sample], 3)
    pw <- power(); a$pid_sib <- signif(pw$pid_sib[match(a$sample, pw$sample)], 3)
    a$n_samples_for_individual <- as.integer(table(a$individual)[a$individual])
    a[order(-a$n_samples_for_individual, a$individual, -a$call_rate), ]
  })

  output$tbl_final <- renderDT(dt(final_tbl()))

  output$tbl_conflict <- renderDT({
    cf <- best()$conflicts
    if (is.null(cf)) return(dt(data.frame(
      message = "Every cluster is internally consistent: all members match all other members.")))
    dt(transform(cf, completeness = round(completeness, 3)))
  })

  output$tbl_cmp  <- renderDT({ x <- res()$cmp$summary
    dt(transform(x, recapture_rate = round(recapture_rate, 3))) })
  output$tbl_ari  <- renderDT({ m <- round(res()$cmp$ari, 3)
    dt(data.frame(method = rownames(m), m, check.names = FALSE)) })

  output$tbl_disagree <- renderDT({ x <- res()$cmp$table
    x <- x[x$disputed, setdiff(names(x), "disputed"), drop = FALSE]
    if (!nrow(x)) return(dt(data.frame(message = "All methods agree on every sample.")))
    dt(x) })

  gx <- reactive({ r <- res(); req(r); r$methods$genalex })

  output$tbl_gx_dist <- renderDT({
    g <- gx(); req(g)
    d <- do.call(rbind, lapply(names(g$by_group), function(k)
      cbind(group = k, g$by_group[[k]]$match_distribution)))
    names(d)[names(d) == "mismatching_loci"] <- "loci differing"
    dt(d)
  })

  output$tbl_gx_near <- renderDT({
    g <- gx(); req(g)
    n <- do.call(rbind, lapply(names(g$by_group), function(k) {
      x <- g$by_group[[k]]$near_matches
      if (is.null(x) || !nrow(x)) NULL else cbind(group = k,
        x[, c("id1", "id2", "n_compared", "n_mismatch")])
    }))
    if (is.null(n)) return(dt(data.frame(
      message = "No near matches: every pair either matches everywhere or differs clearly.")))
    dt(n[order(n$n_mismatch), ])
  })

  output$gx_pid <- renderUI({
    g <- gx(); req(g)
    tags$div(style = "margin-top:1rem", lapply(names(g$by_group), function(k) {
      v <- g$by_group[[k]]$pid_total
      if (is.null(v)) return(NULL)
      tags$div(class = "gid-hint", style = "font-size:.86rem",
        tags$b(k), ": P(ID) = ", tags$code(signif(v[["pid"]], 3)),
        "  \u00b7  P(ID)sib = ", tags$code(signif(v[["pid_sib"]], 3)),
        if (v[["pid_sib"]] > 0.01)
          tags$span(style = "color:#c1502e", "  \u2014 above the 0.01 benchmark; this panel ",
                    "cannot reliably separate siblings")
        else if (v[["pid_sib"]] > 0.001)
          tags$span("  \u2014 adequate, below the 0.01 benchmark")
        else tags$span("  \u2014 comfortable, below 0.001"))
    }))
  })

  output$sethi_alt <- renderUI({
    r <- res(); req(r); st <- r$methods$sethi; req(st)
    mp <- st$matched_pairs
    if (is.null(mp) || !nrow(mp))
      return(tags$p(class = "gid-hint", "No pairs met the threshold."))
    tb <- sort(table(mp$best_alternative), decreasing = TRUE)
    mx <- vapply(st$by_group, function(e) e$max_log10_lambda_observed %||% NA_real_, 0)
    tagList(
      tags$div(style = "display:flex;gap:1.4rem;flex-wrap:wrap;margin:.3rem 0 .7rem",
        lapply(names(tb), function(k) tags$div(
          tags$div(style = "font-family:var(--mono);font-size:1.35rem;color:#1d3557", tb[[k]]),
          tags$div(class = "gid-hint", style = "margin:0", gsub("_", " ", k))))),
      tags$p(class = "gid-hint",
        sprintf("%d pairs called as matches. Strongest evidence in the data: log10 \u039b = %.2f.",
                nrow(mp), max(mx, na.rm = TRUE)),
        if (length(tb) == 1) sprintf(
          " Every matched pair was hardest to distinguish from %s, so that is the alternative worth reporting.",
          gsub("_", " ", names(tb)[1])) else ""))
  })

  output$plot_sweep <- renderPlot({
    s <- res()$sweep
    ggplot(s, aes(max_mismatch, n_individuals)) +
      geom_line(colour = INK, linewidth = 1) +
      geom_point(aes(size = n_conflicting_clusters + 1), colour = INK) +
      geom_vline(xintercept = input$max_mismatch, colour = ACCENT, linetype = 2) +
      scale_size_continuous(range = c(2.2, 6), guide = "none") +
      labs(x = "Mismatching loci allowed", y = "Individuals inferred",
           subtitle = "Bigger points = more clusters whose members do not all match each other") +
      theme_gid()
  })

  ## ------------------------------------------------------- settings sensitivity
  sens <- eventReactive(input$run, {
    p <- prep(); req(p)
    grid <- expand.grid(
      post_cut = sort(unique(c(0.95, 0.99, 0.999, 0.9999, input$post_cut))),
      dropout  = sort(unique(c(0.001, 0.01, 0.05, input$dropout))),
      kinship  = c("unrelated", "full_sib"),
      stringsAsFactors = FALSE)
    withProgress(message = "Re-running across settings", value = 0, {
      out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
        incProgress(1 / nrow(grid))
        r <- gid_by_group(p$gt, p$grp, gid_method_lr,
                          dropout = grid$dropout[k],
                          false_allele = input$false_allele,
                          kinship = grid$kinship[k], post_cut = grid$post_cut[k],
                          min_loci = input$min_loci, linkage = input$linkage,
                          prior_same = if (is.na(input$prior)) NULL else input$prior)
        data.frame(grid[k, ], n_individuals = length(unique(r$assignment$individual)),
                   n_conflicting = r$n_conflict, stringsAsFactors = FALSE)
      }))
    })
    out$current <- out$post_cut == input$post_cut & out$dropout == input$dropout &
      out$kinship == input$kinship
    out
  })

  output$plot_sens <- renderPlot({
    x <- sens()
    x$kin <- ifelse(x$kinship == "full_sib", "vs full siblings", "vs unrelated")
    ggplot(x, aes(factor(post_cut), n_individuals,
                  colour = factor(dropout), group = factor(dropout))) +
      geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
      geom_point(data = x[x$current, ], size = 5, shape = 21, stroke = 1.1,
                 fill = NA, colour = ACCENT) +
      facet_wrap(~kin) +
      scale_colour_manual(values = PAL, name = "assumed dropout") +
      labs(x = "Posterior probability required to accept a match",
           y = "Individuals inferred",
           subtitle = "Circled point is your current setting") +
      theme_gid() + theme(strip.text = element_text(colour = INK, face = "bold"))
  })

  output$tbl_sens <- renderDT({
    x <- sens()
    x$setting <- ifelse(x$current, "<- current", "")
    dt(x[order(x$kinship, x$dropout, x$post_cut),
         c("kinship", "dropout", "post_cut", "n_individuals", "n_conflicting", "setting")])
  })

  ## ------------------------------------------------------------------ output
  params <- reactive({
    r <- res()
    sprintf(paste0(
      "genoID settings\n",
      "----------------------------------------\n",
      "file                 %s\n",
      "samples analysed     %d\n",
      "loci used            %d  (%s)\n",
      "loci dropped by you  %s\n",
      "grouped by           %s\n",
      "\nMethod shown on the Individuals tab: %s\n",
      "\nProbabilistic method\n",
      "  alternative        %s\n",
      "  posterior cutoff   %s\n",
      "  dropout rate       %s\n",
      "  false allele rate  %s\n",
      "  prior P(recapture) %s\n",
      "  min loci per pair  %d\n",
      "  cluster rule       %s\n",
      "\nGut-check methods\n",
      "  mismatch threshold %d loci\n",
      "  Sethi et al. 2016  %s\n",
      "  GenAlEx Matches    %s\n",
      "  allelematch        %s\n",
      "\nResult: %d samples -> %d individuals\n"),
      input$file$name %||% "-", nrow(r$gt), ncol(r$gt),
      paste(colnames(r$gt), collapse = " "),
      if (length(setdiff(detected(), loci()))) paste(setdiff(detected(), loci()), collapse = " ") else "none",
      if (nzchar(input$group_col %||% "")) input$group_col else "not grouped",
      GID_METHODS[[input$method]]$label,
      input$kinship, input$post_cut, input$dropout, input$false_allele,
      if (is.na(input$prior)) "estimated from data" else as.character(input$prior),
      input$min_loci, input$linkage, input$max_mismatch,
      if (isTRUE(input$run_sethi))
        sprintf("run, lambda > %s vs {%s}", input$lambda_cut,
                paste(input$sethi_rel, collapse = ", ")) else "not run",
      if (isTRUE(input$run_genalex))
        sprintf("run, near matches within %d loci", input$near_match) else "not run",
      if (isTRUE(input$run_am)) "run" else "not run",
      nrow(r$gt), length(unique(best()$assignment$individual)))
  })
  output$params_txt <- renderText(params())

  ## ---- downloads -----------------------------------------------------------
  ## Built here, handed to the browser as a Blob. No server route involved, so
  ## this behaves identically in a local Shiny session and in the WebAssembly
  ## build where there is no server at all.
  send_file <- function(filename, content, type = "text/csv;charset=utf-8",
                        b64 = FALSE) {
    session$sendCustomMessage("gid_download", list(
      filename = filename, data = content, b64 = b64, type = type))
  }
  as_csv <- function(df)
    paste(utils::capture.output(utils::write.csv(df, row.names = FALSE)),
          collapse = "\n")

  dl_button <- function(id, name, data_fn) {
    observeEvent(input[[id]], {
      d <- try(data_fn(), silent = TRUE)
      if (inherits(d, "try-error") || is.null(d))
        return(showNotification(
          "Nothing to download yet - press Identify individuals first.",
          type = "warning"))
      send_file(sprintf("genoID_%s_%s.csv", name, format(Sys.Date())), as_csv(d))
    })
  }

  dl_button("dl_assign", "assignments",           function() final_tbl())
  dl_button("dl_pairs",  "pairwise",              function() best()$pairs)
  dl_button("dl_cons",   "individual_genotypes",  function() {
    ic <- gid_individual_consensus(res()$gt, best()$assignment)
    data.frame(individual = rownames(ic$genotypes), ic$genotypes) })

  observeEvent(input$dl_params, {
    p <- try(params(), silent = TRUE)
    if (inherits(p, "try-error"))
      return(showNotification("Press Identify individuals first.", type = "warning"))
    send_file(sprintf("genoID_settings_%s.txt", format(Sys.Date())), p,
              type = "text/plain;charset=utf-8")
  })

  ## the exact gid_ call behind whatever the Individuals tab is showing
  method_call <- reactive({
    a <- function(...) paste(c(...), collapse = ",\n                     ")
    switch(input$method,
      probabilistic = sprintf("gid_method_lr,\n                     %s",
        a(sprintf("dropout = %s", input$dropout),
          sprintf("false_allele = %s", input$false_allele),
          sprintf('kinship = "%s"', input$kinship),
          sprintf("post_cut = %s", input$post_cut),
          sprintf("min_loci = %s", input$min_loci),
          sprintf('linkage = "%s"', input$linkage))),
      sethi = sprintf("gid_method_sethi,\n                     %s",
        a(sprintf("dropout = %s", input$dropout),
          sprintf("false_allele = %s", input$false_allele),
          sprintf('relationships = c("%s")', paste(input$sethi_rel, collapse = '", "')),
          sprintf("lambda_cut = %s", input$lambda_cut),
          sprintf("min_loci = %s", input$min_loci),
          sprintf('linkage = "%s"', input$linkage))),
      threshold = sprintf("gid_method_threshold,\n                     %s",
        a(sprintf("max_mismatch = %s", input$max_mismatch),
          sprintf("min_loci = %s", input$min_loci),
          sprintf('linkage = "%s"', input$linkage))),
      genalex = sprintf("gid_method_genalex,\n                     %s",
        a(sprintf("near_match_loci = %s", input$near_match),
          sprintf("min_loci = %s", input$min_loci),
          sprintf('linkage = "%s"', input$linkage))),
      allelematch = sprintf("gid_method_allelematch,\n                     min_loci = %s",
        input$min_loci),
      exact = sprintf("gid_method_exact,\n                     %s",
        a(sprintf("min_loci = %s", input$min_loci),
          sprintf('linkage = "%s"', input$linkage))))
  })

  observeEvent(input$dl_script, {
    req(input$id_col, length(loci()) > 0)
    send_file("genoID_reproduce.R", sprintf(paste0(
      '## Reproduces exactly what the genoID app just did.\n',
      'source("genoID_core.R")\n\n',
      'raw  <- gid_read("%s")\n',
      'loci <- c(%s)\n',
      'gt   <- gid_matrix(raw, "%s", loci)\n',
      'f    <- gid_filter(gt, min_locus_call = %s, min_sample_call = %s)\n',
      'grp  <- %s\n\n',
      'res  <- gid_by_group(f$gt, grp, %s)\n\n',
      'write.csv(res$assignment, "individual_assignments.csv", row.names = FALSE)\n'),
      input$file$name %||% "genotypes.csv",
      paste0('"', paste(loci(), collapse = '", "'), '"'),
      input$id_col, input$min_locus_call, input$min_sample_call,
      if (nzchar(input$group_col %||% ""))
        sprintf('raw[["%s"]][match(rownames(f$gt), raw[["%s"]])]',
                input$group_col, input$id_col)
      else 'rep("all", nrow(f$gt))',
      method_call()), type = "text/plain;charset=utf-8")
  })

  ## The bundle. Zipping needs a compiled package and a writable filesystem,
  ## neither guaranteed in the browser, so fall back to one concatenated text
  ## file rather than failing.
  observeEvent(input$dl_all, {
    r <- try(res(), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r))
      return(showNotification("Press Identify individuals first.", type = "warning"))
    ic <- gid_individual_consensus(r$gt, best()$assignment)
    parts <- list(
      individual_assignments = final_tbl(),
      pairwise_comparisons   = best()$pairs,
      method_comparison      = r$cmp$summary,
      assignments_all_methods = r$cmp$table,
      probability_of_identity = pid(),
      locus_qc               = gid_locus_stats(r$gt),
      sample_power           = power(),
      individual_genotypes   = data.frame(individual = rownames(ic$genotypes),
                                          ic$genotypes))
    if (!is.null(best()$conflicts)) parts$cluster_conflicts <- best()$conflicts

    zipped <- try({
      d <- file.path(tempdir(), paste0("genoID_", as.integer(Sys.time())))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      for (nm in names(parts))
        utils::write.csv(parts[[nm]], file.path(d, paste0(nm, ".csv")), row.names = FALSE)
      writeLines(params(), file.path(d, "settings.txt"))
      zf <- file.path(tempdir(), "genoID_results.zip")
      unlink(zf)
      zip::zip(zf, list.files(d), root = d)
      gsub("[\r\n]", "", jsonlite::base64_enc(readBin(zf, "raw", file.size(zf))))
    }, silent = TRUE)

    if (!inherits(zipped, "try-error") && nzchar(zipped)) {
      send_file(sprintf("genoID_results_%s.zip", format(Sys.Date())), zipped,
                type = "application/zip", b64 = TRUE)
    } else {
      txt <- unlist(lapply(names(parts), function(nm)
        c(paste0("##### ", nm, " #####"), as_csv(parts[[nm]]), "")))
      txt <- c(txt, "##### settings #####", params())
      send_file(sprintf("genoID_results_%s.txt", format(Sys.Date())),
                paste(txt, collapse = "\n"), type = "text/plain;charset=utf-8")
      showNotification(
        "Your browser could not build a zip, so all results were combined into one text file.",
        type = "message", duration = 8)
    }
  })
}

shinyApp(ui, server)
