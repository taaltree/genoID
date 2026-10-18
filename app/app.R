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
pj <- c("project.R", "app/project.R", "../app/project.R")
mt <- c("map_tab.R", "app/map_tab.R", "../app/map_tab.R")
## local = TRUE matters: Shiny evaluates app.R in its own environment, so a
## default source() would define these in globalenv -- the parent -- where they
## could not see GID_METHODS or dt().
source(mt[file.exists(mt)][1], local = TRUE)
source(pj[file.exists(pj)][1], local = TRUE)

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

## A slider silently clamps a value past its end, which would leave the analysis
## running on a different error rate from the one reported beside it. Widen the
## range instead. The upper ends here match the sliders' own.
GID_RATE_MAX <- c(dropout = 0.20, false_allele = 0.10)
## Where the error-rate sliders start. Nobody chose these, so a run on a file
## with replicates measures the rates instead of trusting them.
GID_RATE_DEFAULT <- c(dropout = 0.005, false_allele = 0.002)
## The multi-tube rule can leave a residual rate of one in a million on paper,
## but no consensus call is that trustworthy: a swapped tube or a contaminated
## extract does not care how many reactions agreed. Never go below this.
GID_CONSENSUS_FLOOR <- 1e-4
## rates span four orders of magnitude, so show significant figures
gid_fmt_rate <- function(x) format(signif(x, 2), scientific = FALSE, drop0trailing = TRUE)
gid_set_rate <- function(session, id, v) {
  v <- suppressWarnings(as.numeric(v))
  if (length(v) == 1 && is.finite(v) && v > GID_RATE_MAX[[id]])
    updateSliderInput(session, id, value = v, max = ceiling(v * 20) / 20)
  else updateSliderInput(session, id, value = v)
}

## One place describing every method: drives the picker, the sidebar sections,
## and the explainer at the top of the Individuals tab. Sethi comes first
## because it is the default: it needs no alternative hypothesis and no prior,
## the two choices a newcomer is least equipped to make.
GID_METHODS <- list(
  sethi = list(
    label  = "Sethi et al. (2016)",
    anchor = "m-sethi", tag = "recommended",
    blurb  = paste("Asks whether two samples' genotypes are better explained by one animal",
                   "than by any relationship that could stand in for it \u2014 unrelated,",
                   "full siblings, parent and offspring \u2014 using an explicit dropout and",
                   "false-allele error model. A match is accepted on the evidence alone",
                   "(\u039b > 1), with no prior."),
    good   = "You do not want to choose an alternative hypothesis, or to argue about a prior. The default.",
    watch  = "\u039b > 1 takes no account of how many pairs you tested, so it merges more on large datasets."),
  probabilistic = list(
    label  = "Likelihood ratio",
    anchor = "m-lr", tag = "recommended",
    blurb  = paste("The same likelihood as Sethi, set against one alternative you name",
                   "rather than the best of several, and turned into a posterior",
                   "probability using a prior estimated from your data."),
    good   = "You want a probability for every pair and control over how strict the cutoff is.",
    watch  = "You must state what \u201ca different animal\u201d means. For group-living species that is full siblings, not unrelated animals."),
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
  ## DT renders server-side, so the browser only ever holds the page currently
  ## on screen. DataTables' own CSV export can see nothing else, which is why it
  ## produced a 12-row file and forced people to page through and download each
  ## page separately. Stash the whole frame under this output's id and let the
  ## button ask R for it instead.
  sess <- shiny::getDefaultReactiveDomain()
  oid  <- tryCatch(shiny::getCurrentOutputInfo()$name, error = function(e) NULL)
  if (!is.null(sess) && !is.null(oid)) {
    if (is.null(sess$userData$gid_tables))
      sess$userData$gid_tables <- new.env(parent = emptyenv())
    assign(oid, x, envir = sess$userData$gid_tables)
  }
  dl_btn <- list(extend = "csv", text = "Download CSV", className = "gid-dt-dl",
                 action = JS(
    "function(e, dt, node, config) {",
    "  var id = $(dt.table().container()).closest('.html-widget').attr('id');",
    "  Shiny.setInputValue('gid_dl_table', {id: id, nonce: Math.random()},",
    "                      {priority: 'event'});",
    "}"))

  datatable(
    x, rownames = FALSE, extensions = "Buttons",
    options = c(list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
                     buttons = list(dl_btn)),
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
    .gid-advice{background:#fff;border:1px solid #e3e9ef;border-left:3px solid #1d3557;
      border-radius:0 4px 4px 0;padding:.9rem 1.1rem;margin-bottom:1rem;
      box-shadow:0 1px 3px rgba(29,53,87,.06)}
    .gid-advice-head{margin-bottom:.55rem;font-size:.95rem;cursor:pointer}
    .gid-advice-row{display:grid;grid-template-columns:7.2rem 1fr;gap:.9rem;
      padding:.65rem 0;border-top:1px solid #eef2f6;align-items:start}
    .gid-advice-title{font-weight:600;font-size:.92rem;color:#1d3557}
    .gid-advice-why{font-size:.85rem;color:#4a6067;line-height:1.5;margin-top:.15rem}
    .gid-advice-ev{font-size:.82rem;color:#6b7a8f;margin-top:.3rem;
      font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
    .gid-sev{display:inline-block;font-size:.62rem;font-weight:700;letter-spacing:.07em;
      text-transform:uppercase;padding:.2em .55em;border-radius:2px;white-space:nowrap}
    .gid-sev-high{background:#fdf3ee;color:#c1502e}
    .gid-sev-medium{background:#f6eddc;color:#a8762a}
    .gid-sev-low{background:#eef3f2;color:#6b7a8f}
    @media (max-width:640px){.gid-advice-row{grid-template-columns:1fr;gap:.25rem}}
    .gid-label{font-weight:700;font-size:.78rem;letter-spacing:.05em;text-transform:uppercase;
               color:#6b7a8f;margin:.2rem 0 .4rem}
    .gid-legend{max-height:330px;overflow:auto}
    .gid-legend-row{display:flex;align-items:center;gap:.45rem;font-size:.82rem;padding:.1rem 0}
    .gid-swatch{width:13px;height:13px;border-radius:3px;flex:none;border:1px solid rgba(0,0,0,.18)}
    .gid-kv td{padding:.22rem .4rem;font-size:.84rem;border-color:#eef2f6}
    .gid-kv td:first-child{white-space:nowrap;width:1%%}
    .gid-kv td:last-child{word-break:break-word}
    .gid-geno td{padding:.12rem .4rem;font-size:.78rem;border-color:#f2f5f8}
    .gid-geno code{font-size:.76rem}
    .leaflet-container{border-radius:4px}
    .nav-link{font-size:.9rem}
    .value-box-showcase{max-width:64px!important;flex-basis:64px!important;padding:0!important}
    .value-box-showcase svg,.value-box-showcase .fa,.value-box-showcase i{
      height:1.6rem!important;width:1.6rem!important;font-size:1.6rem!important;opacity:.8}
    .bslib-value-box .value-box-area{padding:.35rem .2rem}
    .gid-empty{text-align:center;padding:3.2rem 1rem;color:%s}
    .gid-empty h4{color:%s;font-weight:600;font-size:1.05rem;margin-bottom:.4rem}
    .gid-empty p{font-size:.88rem;max-width:34rem;margin:0 auto .35rem}
    .dataTables_wrapper{font-size:.86rem}
    .gid-fmt{display:inline-block;vertical-align:top;text-align:left;margin:1rem .9rem 0;
      background:#fff;border:1px solid #e3e9ef;border-radius:4px;padding:.8rem 1rem}
    .gid-fmt-h{font-size:.7rem;font-weight:600;letter-spacing:.07em;text-transform:uppercase;
      color:%s;margin-bottom:.45rem}
    table.gid-ex{border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
      font-size:.78rem}
    table.gid-ex th{font-family:system-ui,sans-serif;font-size:.64rem;letter-spacing:.04em;
      text-transform:uppercase;color:#8fa3a6;font-weight:600;padding:.15rem .6rem .3rem;
      border-bottom:1px solid #dde5e4;text-align:left}
    table.gid-ex td{padding:.22rem .6rem;border-bottom:1px solid #f0f4f3;color:#33474d}
    table.gid-ex .hi{background:#f6eddc;color:#a8762a}
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
     MUTED, INK, INK)))),

  sidebar = sidebar(
    width = 330, class = "bg-white",

    fileInput("file", "Genotype file", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
    hint("One column per locus. One row per sample, or one row per PCR ",
         "replicate. See ", tags$b("Your file format"), " on the Methods tab."),
    actionLink("load_demo", "or load the example dataset", class = "gid-hint"),
    uiOutput("demo_note"),
    tags$hr(style = "margin:.9rem 0 .8rem"),
    gid_project_ui(),
    conditionalPanel("output.is_excel == true",
                     selectInput("sheet", "Worksheet", choices = NULL)),

    selectInput("id_col", "Sample ID column", choices = NULL),
    selectizeInput("group_col", "Analyze separately by (optional)", choices = NULL,
                   options = list(placeholder = "e.g. species, study area")),
    hint("Blocking stops samples from different species or populations ",
         "from ever being called the same animal, and keeps allele frequencies separate."),

    conditionalPanel(
      "output.has_reps == true",
      tags$div(class = "gid-flag gid-ok", style = "font-size:.8rem",
        tags$b("This file has several rows per sample."), " They look like PCR ",
        "replicates, so genoID analyzes them both ways \u2014 every reaction ",
        "directly, and collapsed to a consensus first \u2014 and reports both. ",
        "See Methods."),
      selectInput("rep_col", "Replicate label column", choices = NULL),
      selectizeInput("rep_drop", "Row labels to exclude", choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "e.g. consensus")),
      hint("Exclude any summary rows, such as a pre-computed consensus. ",
           "Everything left is treated as one PCR replicate."),
      radioButtons("rep_mode", "Route the results tabs show",
                   c("Every replicate reaction" = "reps",
                     "A consensus genotype per sample" = "consensus"),
                   selected = "reps"),
      hint("Both routes run either way, and the Method comparison tab sets ",
           "them side by side. Using every reaction, with an error rate ",
           "measured for each sample, was the more accurate in simulations. ",
           "Where the two disagree, check those samples by hand."),
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
                choices = GID_METHOD_CHOICES, selected = "sethi"),
    uiOutput("method_hint"),

    ## Settings for the chosen method. The other methods keep their values while
    ## hidden, so the comparison tab can still run all five.
    conditionalPanel(
      "input.method == 'probabilistic'",
      conditionalPanel(
        "output.has_coords == true",
        checkboxInput("use_spatial", "Also use where samples were found", FALSE),
        hint("Two scats close together are more likely to be one animal. This ",
             "adds that as evidence, after Augustine et al. (2020). It helps ",
             "when home ranges are small next to the area you sampled, and can ",
             "mislead when they overlap heavily \u2014 check the diagnostics on ",
             "the Individuals tab before trusting it.")),
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
      "input.method == 'exact'",
      hint("No tolerance to set \u2014 that is the point of this baseline.")),

    ## error model: used by both likelihood methods
    conditionalPanel(
      "input.method == 'probabilistic' || input.method == 'sethi'",
      ## A residual consensus rate is often a few in ten thousand, so the step
      ## has to be finer than that or the slider shows a rate in force as 0.
      sliderInput("dropout", "Allelic dropout rate", 0, 0.20,
                  GID_RATE_DEFAULT[["dropout"]], step = 0.0001),
      sliderInput("false_allele", "False allele rate", 0, 0.10,
                  GID_RATE_DEFAULT[["false_allele"]], step = 0.0001),
      actionButton("est_error", "Estimate these from my data",
                   icon = icon("wand-magic-sparkles"),
                   class = "btn-sm btn-outline-primary w-100"),
      uiOutput("est_error_note"),
      hint("Guessing these is the weakest part of the analysis. If you uploaded ",
           "replicates and leave these alone, they are measured from the ",
           "replicates when you press Identify individuals; without replicates ",
           "the button can at least get the order of magnitude right. See ",
           tags$b("Measuring your own error rates"), " on the Methods tab.")),

    accordion(
      open = FALSE,
      accordion_panel(
        "Settings shared by every method", value = "shared", icon = icon("sliders"),
        numericInput("min_loci", "Minimum loci compared per pair", 15, 1, 500, 1),
        hint("How much overlap two samples need before they may be called the ",
             "same animal. Pairs sharing fewer called loci than this are left ",
             "unmatched rather than guessed at. Too low and samples match on ",
             "thin evidence; too high and real recaptures are never compared."),

        numericInput("min_sample_call", "Minimum sample call rate", 0.5, 0, 1, 0.05),
        hint("Drops a whole sample that was called at fewer than this fraction ",
             "of loci. This is the curation step, and it is why you can upload ",
             "raw data: a sample with almost no calls matches nothing, so, left ",
             "in, it becomes a phantom individual and inflates your count."),

        numericInput("min_locus_call", "Minimum locus call rate", 0.25, 0, 1, 0.05),
        hint("Drops a whole locus that worked in fewer than this fraction of ",
             "samples. Applied before the sample filter, so a marker that failed ",
             "across the board cannot drag otherwise good samples below the line."),
        checkboxInput("drop_flagged", "Treat quality-flagged calls as missing", FALSE),
        hint("Cells like ", tags$code("CT*"), " carry a flag from your genotyping ",
             "pipeline. By default the flag is stripped and the call kept."),

        checkboxInput("drop_mixtures", "Exclude samples that look like two animals", FALSE),
        numericInput("mix_cut", "Heterozygosity excess to call a mixture", 0.18, 0, 1, 0.01),
        hint("A scat holding DNA from two animals reads heterozygous wherever ",
             "they differ, so it carries too many heterozygous loci and gets ",
             "reported as a new animal. Excess is observed minus expected ",
             "heterozygosity; 0.18 is the threshold in Eriksson et al. (2026). ",
             "See the mixture card on Data & QC before switching this on."),
        radioButtons("linkage", "Cluster rule",
                     c("Single linkage (transitive)" = "single",
                       "Complete linkage (all pairs must match)" = "complete"),
                     selected = "single"),
        hint("Single linkage lets A-B and B-C merge A, B and C. Complete linkage ",
             "requires every pair inside a cluster to match. Compare both: if the ",
             "answer moves, some clusters are held together by one edge.")),
      accordion_panel(
        "Which methods to run", value = "which", icon = icon("list-check"),
        hint("Sethi, the likelihood ratio, exact matching and the mismatch ",
             "threshold always run, so the comparison tab can set them side by ",
             "side. GenAlEx is optional."),
        checkboxInput("run_genalex", "GenAlEx Matches", TRUE))
    ),

    ## Map controls. Keyed on the open tab rather than the method, because the
    ## map deliberately lets you color by a model other than the one the
    ## Individuals tab is showing.
    conditionalPanel(
      "input.nav == 'Scat map'",
      tags$hr(),
      tags$p(class = "gid-label", "Map"),
      selectInput("map_model", "Color individuals by", choices = NULL),
      hint("Switch models to see which samples change hands. The points stay ",
           "put; only the coloring changes."),
      checkboxInput("map_grey", "Gray out animals seen once", TRUE),
      uiOutput("map_coord_ui"),
      uiOutput("map_meta_ui"),

      tags$hr(),
      tags$p(class = "gid-label", "Show on the map"),
      hint("Leave a filter empty to show everything. They combine, so you can ",
           "ask for one year of females and get only those."),
      checkboxGroupInput("map_show_sex", "Sex", inline = TRUE,
                         choices = c("Female" = "F", "Male" = "M",
                                     "Not called" = "U")),
      selectizeInput("map_show_year", "Year", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All years")),
      selectizeInput("map_show_who", "Individuals", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All individuals")),
      actionButton("map_show_reset", "Show everything",
                   class = "btn-sm btn-outline-secondary"),

      tags$hr(),
      tags$p(class = "gid-label", "Link samples of the same animal"),
      radioButtons("map_link_style", NULL,
                   c("No links" = "none",
                     "Spider (lines to center)" = "spider",
                     "Polygon (convex hull)" = "polygon"),
                   selected = "none"),
      selectizeInput("map_link_who", "Which animals", choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "All animals")),
      actionButton("map_link_none", "Clear selection",
                   class = "btn-sm btn-outline-secondary"),
      hint("Leave the list empty to link every animal. Pick one or more to link ",
           "only those. Only animals with two or more mapped samples appear."),
      tags$hr(),
      tags$p(class = "gid-label", "Save the map"),
      tags$div(
        class = "d-flex gap-2",
        actionButton("dl_map_pdf", "PDF", icon = icon("file-pdf"),
                     class = "btn-sm btn-outline-secondary flex-fill"),
        actionButton("dl_map_jpg", "JPG", icon = icon("file-image"),
                     class = "btn-sm btn-outline-secondary flex-fill")),
      hint("PDF is vector, for figures you will scale or edit. JPG is raster. ",
           "Both draw the points and links as a clean figure without the ",
           "basemap tiles \u2014 tiles are copyrighted images and cannot go into ",
           "a vector file."),
      actionButton("dl_map_html", "Shareable interactive map",
                   icon = icon("share-nodes"),
                   class = "btn-sm btn-outline-primary w-100"),
      hint("One HTML file holding this map, its points and every popup. Send ",
           "it to a collaborator and they can pan, zoom and click it in any ",
           "browser \u2014 no genoID and no R needed. Put it on any web host or ",
           "shared drive and it becomes a link."),
      numericInput("map_fig_width", "Figure width (inches)", 9, 3, 20, 0.5),
      checkboxInput("map_fig_labels", "Name the linked animals on the figure", TRUE),
      hint("Turn the names off when many animals overlap.")),

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
    conditionalPanel("output.has_reps == true",
      card(
        card_header("How consistent each sample's replicates are"),
        hint("The share of loci where a sample's own reactions do not all agree. ",
             "One error rate for the whole dataset fits when these sit close ",
             "together. When a minority of samples disagrees with itself far ",
             "more than the rest, one rate under-describes them, and their ",
             "recaptures get split off as new animals \u2014 which is why the ",
             "replicate route measures an error rate for every sample."),
        uiOutput("rep_quality_note"),
        layout_columns(
          col_widths = c(7, 5),
          plotOutput("plot_rep_quality", height = 280),
          DTOutput("tbl_rep_quality")))),
    card(
      card_header("Samples that may hold two animals"),
      hint("A scat carrying DNA from two animals reads heterozygous wherever ",
           "they differ, so it shows more heterozygous loci than any real ",
           "individual should \u2014 and, being a genotype nobody has, it is ",
           "reported as a new animal. After Eriksson et al. (2026), ",
           tags$i("Mol Ecol Resour"), " 26:e70164."),
      uiOutput("mix_summary"),
      layout_columns(
        col_widths = c(7, 5),
        plotOutput("plot_mixture", height = 300),
        DTOutput("tbl_mixture"))
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
    uiOutput("route_summary"),
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      vbox("Individuals", "vb_ind", "paw", "primary", key = "n_individuals"),
      vbox("Sampled once", "vb_single", "circle-dot", "light", key = "n_singletons"),
      vbox("Samples per animal", "vb_persample", "layer-group", "light",
           key = "median_samples"),
      vbox("Largest cluster", "vb_max", "maximize", "light", key = "max_cluster")
    ),
    ## The answer comes before anything explaining it. These cards used to sit
    ## above the table and pushed it past 2,500px down the page, which read as
    ## the assignments not being there at all.
    card(card_header("Individual assignments"), DTOutput("tbl_final")),
    uiOutput("advice_card"),
    uiOutput("spatial_diag"),
    uiOutput("recap_detail"),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header(textOutput("evid_title", inline = TRUE)),
           plotOutput("plot_lr", height = 340),
           uiOutput("evid_note")),
      card(card_header("Locus mismatches between samples"), plotOutput("plot_mm", height = 340),
           hint("A clean panel gives two piles with a gap. Pairs sitting in the gap ",
                "are the ones worth looking at by hand."))
    ),
    conditionalPanel("input.method == 'probabilistic'",
      card(card_header("Choosing your posterior cutoff"),
        hint("The cutoff is the one setting with no data behind it \u2014 unless ",
             "you use the fact that the method already gives every pair a ",
             "probability. Adding those up says how many mistakes each cutoff ",
             "costs, in both directions, so the choice becomes a stated trade ",
             "instead of a convention."),
        tags$div(
          class = "d-flex align-items-center gap-2 flex-wrap",
          style = "margin-bottom:.6rem",
          tags$span(class = "gid-hint", style = "margin:0",
                    "How many wrongly merged pairs across the whole dataset are ",
                    "you willing to accept?"),
          tags$div(style = "width:7rem",
                   numericInput("err_budget", NULL, 1, 0.05, 100, 0.05))),
        uiOutput("calib_verdict"),
        plotOutput("plot_calib", height = 320),
        DTOutput("tbl_calib"))),
    card(card_header("Clusters that are not internally consistent"),
         hint("Every sample in a cluster should match every other sample in it. ",
              "Where that fails, the cluster is being held together by a chain and ",
              "may be two animals. Inspect these before reporting a population size."),
         DTOutput("tbl_conflict"))
  ),

  nav_panel(
    "Sample map", icon = icon("sitemap"),
    uiOutput("run_status_map"),
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      vbox("Individuals", "vb_map_ind", "paw", "primary", key = "n_individuals"),
      vbox("Confident samples", "vb_map_ok", "circle-check", "light", key = "status"),
      vbox("Worth checking", "vb_map_chk", "magnifying-glass", "light", key = "status"),
      vbox("Uncertain", "vb_map_bad", "triangle-exclamation", "light", key = "status")
    ),
    card(card_header("Samples to look at before you report anything"),
         hint("Every sample whose assignment would change under a slightly ",
              "different cutoff, ordered worst first. The rival column names the ",
              "sample it nearly matched instead \u2014 that is the comparison to ",
              "make by hand."),
         uiOutput("map_review_note"),
         DTOutput("tbl_review")),
    card(card_header("How securely each sample is placed"),
         plotOutput("plot_margin", height = 330),
         hint("Each point is a sample. Far right means the assignment survives ",
              "any reasonable cutoff; near zero means this sample is the reason ",
              "your answer depends on where you drew the line.")),
    card(card_header("Every sample and the animal it belongs to"),
         hint("The full map. Download it as your capture history."),
         DTOutput("tbl_map")),
    card(card_header("Roster: one row per animal"),
         hint("The same information the other way around \u2014 each animal and ",
              "the samples that make it up."),
         DTOutput("tbl_roster"))
  ),

  gid_map_tab_ui(),

  nav_panel(
    "Method comparison", icon = icon("code-compare"),
    uiOutput("run_status_cmp"),
    card(card_header("How each method resolved the same data"), DTOutput("tbl_cmp")),
    conditionalPanel("output.has_reps == true",
      card(
        card_header("Every replicate, or a consensus first?"),
        hint("With replicates in the file, both likelihood methods run twice: on ",
             "every PCR reaction, with each sample's own error rate, and on one ",
             "consensus genotype per sample. The routes assume different things ",
             "about error \u2014 that a sample's reactions err independently at ",
             "that sample's rate, or that every consensus call is equally ",
             "reliable \u2014 so where they disagree, one assumption is failing ",
             "for those samples, and they are the ones to check by hand. Where ",
             "they agree, the count does not depend on the choice."),
        uiOutput("route_note"),
        DTOutput("tbl_routes"),
        tags$div(class = "gid-hint", style = "margin-top:.8rem",
                 tags$b("Samples placed differently by the two routes")),
        DTOutput("tbl_route_moved"))),
    card(
      card_header("Why the two recommended methods disagree"),
      hint("Both use the same genotype likelihood, so a disagreement can only ",
           "come from one of two places: the alternative each divides by, or ",
           "the rule each uses to decide. Naming which, pair by pair, is the ",
           "difference between an argument about a threshold and a sign that ",
           "one model is asking the wrong question."),
      uiOutput("disagree_summary"),
      DTOutput("tbl_lr_sethi")),
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
    card(card_header("Sethi et al. \u2014 which relationship competed?"),
      hint("For every pair called a match, the alternative hypothesis that came closest. ",
           "If one relationship dominates, that is the null your other analyses should ",
           "be using too."),
      uiOutput("sethi_alt")),
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
         hint("The same data, re-run with the likelihood ratio on consensus ",
              "genotypes across the settings that matter most. If the count barely ",
              "moves, the number is a property of your data. If it swings, it is a ",
              "property of your choices, and the range belongs in your results. ",
              "Your current setting is marked."),
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
        actionButton("dl_map", "Sample map with confidence",
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
                  "P(ID) after Waits et al. 2001 · match calling after Sethi et al. 2016"))
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
  ## The analysis is driven by a counter rather than by the button directly, so
  ## reopening a project can recompute its results without faking a click.
  run_count <- reactiveVal(0L)
  observeEvent(input$run, run_count(run_count() + 1L))
  ## Where the map was looking when a project was saved. renderLeaflet consumes
  ## this once, then goes back to following the user's own panning.
  restored_view <- reactiveVal(NULL)
  ## Map filter selections waiting for their menus to exist. The map module
  ## claims these the moment it builds those menus, which is the only point at
  ## which selectize will accept them.
  pending_map <- reactiveVal(NULL)
  ## What the loaded data should be called. Tracked separately from
  ## input$file$name because a project restores data without a file upload.
  source_name <- reactiveVal("")
  observeEvent(input$file, source_name(input$file$name))
  observeEvent(input$load_demo, source_name("example dataset"))
  observeEvent(input$load_demo, {
    f <- c("demo/demo_genotypes.csv", "app/demo/demo_genotypes.csv")
    f <- f[file.exists(f)]
    if (!length(f)) return(showNotification("Example file not found.", type = "error"))
    raw(gid_read(f[1])); demo_loaded(TRUE)
    ## The sliders are deliberately NOT set to the truth. The example ships
    ## with PCR replicates, so the rates come back out of the data, either from
    ## "Estimate these from my data" or on the first run.
    updateNumericInput(session, "min_loci", value = 12)
  })
  observeEvent(input$file, demo_loaded(FALSE))

  output$demo_note <- renderUI({
    if (!isTRUE(demo_loaded())) return(NULL)
    tags$div(class = "gid-flag", style = "margin:.5rem 0 .2rem;font-size:.8rem",
      tags$b("Example data loaded."), " 55 samples from 28 known individuals in two ",
      "species, each amplified ", tags$b("three times"), " and followed by the lab's ",
      "consensus row \u2014 220 rows in all. It also carries sex, collection dates ",
      "across two field seasons, and coordinates.",
      tags$div(style = "margin-top:.3rem",
        tags$b("Try this first: "), "press ", tags$b("Estimate these from my data"),
        ", or just ", tags$b("Identify individuals"), ", which measures the rates ",
        "itself. The data were simulated at ", tags$code("dropout 0.015"), " and ",
        tags$code("false allele 0.002"), " per reaction. Dropout comes back close ",
        "to 0.015; the false-allele estimate reads a little high, around 0.005, ",
        "with an interval that still covers 0.002 \u2014 that is the family group ",
        "below, not a mistake."),
      tags$div(style = "margin-top:.3rem",
        "The true answer is in ", tags$code("TrueIndividual"), ", so you can score ",
        "any method against it. The second species is a family group: relatives ",
        "are the hard case, and pooling them inflates the false-allele estimate."))
  })

  output$has_data <- reactive(!is.null(raw()))
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)

  output$empty_state <- renderUI({
    if (!is.null(raw())) return(NULL)
    ex <- function(hdr, rows, hi) tags$table(class = "gid-ex",
      tags$thead(tags$tr(lapply(hdr, function(h)
        tags$th(class = if (h %in% hi) "hi" else "", h)))),
      tags$tbody(lapply(rows, function(r) tags$tr(Map(function(v, h)
        tags$td(class = if (h %in% hi) "hi" else "", v), r, hdr)))))

    tags$div(class = "gid-empty",
      icon("dna", class = "fa-2x", style = paste0("color:", SOFT, ";margin-bottom:.8rem")),
      tags$h4("Upload a genotype table to begin"),
      tags$p("One column per locus. Genotypes as two alleles in one cell (",
             tags$code("AG"), ") or separated (", tags$code("120/124"), "). ",
             "Missing as ", tags$code("00"), ", ", tags$code("NA"), " or blank. ",
             "CSV, TSV or Excel."),

      tags$div(class = "gid-fmt",
        tags$div(class = "gid-fmt-h", "One row per sample"),
        ex(c("SampleID", "Species", "LOC01", "LOC02", "LOC03"),
           list(c("WFS_001", "wolf", "AG", "CC", "00"),
                c("WFS_002", "wolf", "AA", "CT", "TT")), "SampleID"),
        tags$p(class = "gid-hint", style = "margin:.45rem 0 0",
               "Use this when replicates are already collapsed to a consensus.")),

      tags$div(class = "gid-fmt",
        tags$div(class = "gid-fmt-h", "One row per PCR replicate"),
        ex(c("SampleID", "Rep", "Species", "LOC01", "LOC02", "LOC03"),
           list(c("WFS_001", "a", "wolf", "AG", "CC", "00"),
                c("WFS_001", "b", "wolf", "AA", "CC", "CT"),
                c("WFS_001", "c", "wolf", "AG", "00", "00"),
                c("WFS_002", "a", "wolf", "AA", "CT", "TT")),
           c("SampleID", "Rep")),
        tags$p(class = "gid-hint", style = "margin:.45rem 0 0",
               "Repeat the sample ID once per reaction. The app spots this and ",
               "offers to use every observation directly, which beats collapsing ",
               "them first.")),

      tags$p(style = "margin-top:1.2rem",
        "Loci are detected automatically and listed in the sidebar for you to ",
        "correct. Full details, including microsatellites and quality flags, are ",
        "on the ", tags$b("Methods"), " tab under ", tags$b("Your file format"), "."))
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

  ## ---- estimating the error rates from the data itself --------------------
  err_est <- reactiveVal(NULL)
  ## A measurement belongs to the data it was made on.
  observeEvent(list(raw(), input$id_col, input$rep_col, input$rep_drop, input$use_loci),
               err_est(NULL), ignoreInit = TRUE)

  ## The rates the app put in the sliders itself. Sliders still holding these
  ## were never chosen by anyone, so a run on replicates measures the rates
  ## rather than trusting them. That matters most for the replicate route: left
  ## at the defaults, it would read a residual consensus rate as a per-reaction
  ## one, expect ten times too little error, and split recaptures off as new
  ## animals.
  rates_auto <- reactiveVal(GID_RATE_DEFAULT)
  set_rates <- function(d, f) {
    d <- round(d, 4); f <- round(f, 4)
    rates_auto(c(dropout = d, false_allele = f))
    gid_set_rate(session, "dropout", d)
    gid_set_rate(session, "false_allele", f)
  }
  ## within about half a slider step, since the slider rounds what it is given
  rates_untouched <- function() {
    v <- c(input$dropout %||% NA_real_, input$false_allele %||% NA_real_)
    all(is.finite(v)) && all(abs(v - rates_auto()) < 6e-5)
  }

  ## Per-reaction rates measured from the replicates, and the residual rates
  ## they leave on a consensus call. The replicate route needs the first, the
  ## consensus route the second.
  measure_error <- function(p, progress = NULL) {
    rp <- p$reps_all
    e <- gid_estimate_error(p$gt, rp, progress = progress)
    if (!identical(e$method, "replicates")) return(e)
    pe <- gid_propagate_error(e$dropout, e$false_allele, n_rep = p$n_rep,
                              rule = "taberlet", hom_n = p$hom_n, het_n = p$het_n,
                              per_rep_missing = mean(is.na(rp$gt)))
    c(e, list(consensus_dropout = max(pe[["dropout"]], GID_CONSENSUS_FLOOR),
              consensus_false   = max(pe[["false_allele"]], GID_CONSENSUS_FLOOR)))
  }
  ## the measured rates on the scale of one route
  measured_for <- function(e, route) {
    if (is.null(e) || !identical(e$method, "replicates")) return(NULL)
    if (route == "reps") c(dropout = e$dropout, false_allele = e$false_allele)
    else c(dropout = e$consensus_dropout, false_allele = e$consensus_false)
  }

  observeEvent(input$est_error, {
    p <- try(prep(), silent = TRUE)
    if (inherits(p, "try-error") || is.null(p))
      return(showNotification("Load a file and choose your loci first.", type = "warning"))
    ## The bar used to be pinned at 0.4 for the whole computation, so a slow
    ## estimate was indistinguishable from a hang. It now reports each stage.
    e <- try(withProgress(
      message = "Measuring genotyping error", value = 0.05,
      detail = "preparing",
      if (is.null(p$reps_all)) gid_estimate_error(p$gt)
      else measure_error(p, progress = function(v, d) setProgress(value = v, detail = d))),
      silent = TRUE)
    if (inherits(e, "try-error"))
      return(showNotification(paste("Could not estimate:", conditionMessage(attr(e, "condition"))),
                              type = "error"))
    err_est(e)
    ## The sliders hold the rate for the route the tabs show: per reaction on
    ## the replicate route, the residual rate of a consensus call otherwise.
    r <- measured_for(e, if (is.null(p$reps)) "consensus" else "reps")
    if (!is.null(r)) set_rates(r[["dropout"]], r[["false_allele"]])
    else set_rates(e$dropout, max(0.001, round(e$dropout / 3, 3)))
  })

  output$est_error_note <- renderUI({
    e <- err_est(); if (is.null(e)) return(NULL)
    if (identical(e$method, "replicates"))
      tags$div(class = "gid-flag gid-ok", style = "font-size:.78rem;margin:.5rem 0",
        tags$b(if (isTRUE(e$auto)) "Measured from your replicates when you ran the analysis."
               else "Measured from your replicates."),
        " Maximum likelihood over ", e$n_samples, " samples, integrating over ",
        "the unknown true genotype so no consensus is involved. Per reaction:",
        tags$div(style = "margin-top:.3rem;font-family:var(--bs-font-monospace)",
          sprintf("dropout %.3f [%.3f-%.3f]   false allele %.3f [%.3f-%.3f]",
                  e$dropout, e$dropout_ci[1], e$dropout_ci[2],
                  e$false_allele, e$false_ci[1], e$false_ci[2])),
        tags$div(style = "margin-top:.3rem",
          "These drive the replicate route, where a sample whose reactions err ",
          "more often is also given higher rates of its own. The consensus ",
          "route uses the much smaller rates left ",
          "after the multi-tube rule: ",
          tags$code(sprintf("dropout %s, false allele %s",
                            gid_fmt_rate(e$consensus_dropout), gid_fmt_rate(e$consensus_false))), ". ",
          "The sliders show the rates for the route the tabs show."))
    else
      tags$div(class = "gid-flag", style = "font-size:.78rem;margin:.5rem 0",
        tags$b("No replicates, so this is a rough estimate only."), " Dropout ",
        "turns heterozygotes into homozygotes, so it shows up as a heterozygote ",
        "deficit: ", tags$code(sprintf("d ~ %.3f", e$dropout)), " across ",
        e$n_samples, " unique individuals.",
        tags$div(style = "margin-top:.3rem",
          "Treat this as the right order of magnitude, not a number to report. ",
          "Inbreeding and population structure push it up; false alleles push it ",
          "down by manufacturing heterozygotes, so it is not a ceiling either. ",
          "The false-allele rate cannot be estimated this way at all and has been set ",
          "to a third of the dropout estimate."),
        tags$div(style = "margin-top:.3rem",
          tags$b("To measure it properly, genotype some samples twice"), " and ",
          "upload the replicates. Two reactions on a subset is enough."))
  })

  output$rep_note <- renderUI({
    if (!isTRUE(has_reps())) return(NULL)
    if (identical(input$rep_mode, "reps"))
      tags$p(class = "gid-hint",
        tags$b("The dropout and false-allele rates below mean per-reaction "),
        tags$b("rates"), ", not the residual rates of a consensus call. They are ",
        "typically ten to a hundred times larger \u2014 a few percent rather than ",
        "a few tenths of a percent. Each sample's own rate is measured on top.")
    else
      tags$p(class = "gid-hint",
        "Replicates are collapsed with the multi-tube rule: a heterozygote ",
        "accepted on two replicates, a homozygote only on three. The rates ",
        "below are the small residual rates of those consensus calls.")
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
    ## Rows with no usable identifier are removed here, before anything is keyed
    ## by them. Left in, they all collapse into one pseudo-sample whose
    ## "replicates" are unrelated animals, which both crashes the error
    ## estimator and, when it does run, wrecks the rates it reports.
    blank <- gid_blank_ids(df, input$id_col)
    if (any(blank)) df <- df[!blank, , drop = FALSE]
    validate(need(nrow(df) > 0, paste(
      "Every row has a blank", input$id_col, "\u2014 choose a different Sample ID column.")))
    ids <- as.character(df[[input$id_col]])
    sep <- gid_guess_sep(df, exclude = input$id_col)
    use_reps <- isTRUE(has_reps()) && identical(input$rep_mode, "reps") &&
                nzchar(input$rep_col %||% "")

    ## Two animals in one scat genotype as a blend of both, and because that
    ## blend is a genotype nobody has, it is reported as a new animal. Checked
    ## on the filtered matrix; dropped only if asked, since whether a flagged
    ## sample is a mixture or just a heterozygous animal is the user's call.
    mix_step <- function(g) {
      none <- list(gt = g, mixture = NULL, dropped = character(0))
      if (nrow(g) < 10) return(none)
      mxt <- try(gid_mixture_check(
        g, alpha = 0.05,
        excess_cut = if (is.na(input$mix_cut %||% NA)) NA_real_ else input$mix_cut,
        min_loci = max(10, input$min_loci %||% 10)), silent = TRUE)
      if (inherits(mxt, "try-error")) return(none)
      bad <- mxt$sample[mxt$flag == "possible mixture"]
      if (isTRUE(input$drop_mixtures) && length(bad) && nrow(g) - length(bad) >= 2)
        list(gt = g[setdiff(rownames(g), bad), , drop = FALSE],
             mixture = mxt, dropped = bad)
      else list(gt = g, mixture = mxt, dropped = character(0))
    }

    ## ---- replicate route: several rows per sample, one genotype per sample
    if (isTRUE(has_reps()) && nzchar(input$rep_col %||% "")) {
      lab  <- as.character(df[[input$rep_col]])
      keep <- !lab %in% (input$rep_drop %||% character(0))
      validate(need(sum(keep) > 0,
        "Every row was excluded. Clear some entries from 'Row labels to exclude'."))
      rgt  <- gid_matrix(df[keep, , drop = FALSE], input$id_col, loci(), sep = sep,
                         strip_flags = !isTRUE(input$drop_flagged))
      rsm  <- ids[keep]
      ## A homozygote needs three agreeing reactions under the multi-tube
      ## rule, which a lab that runs duplicates can never supply: every
      ## homozygous locus would come out missing. So the rule asks for no more
      ## agreeing reactions than a typical sample has -- with two, both must
      ## agree; with one, the call stands as it is.
      n_rep <- max(1L, as.integer(round(stats::median(table(rsm)))))
      hom_n <- min(3L, n_rep)
      het_n <- min(2L, n_rep)
      gt   <- gid_consensus_from_reps(rgt, rsm, hom_n = hom_n, het_n = het_n)
      ## A sample's label can sit on any of its rows -- some pipelines write
      ## the species only on the summary row excluded above -- so take the
      ## first non-blank one across all of them.
      grp  <- if (nzchar(input$group_col %||% "")) {
        gv <- trimws(as.character(df[[input$group_col]]))
        gv[!nzchar(gv)] <- NA
        first <- tapply(gv, ids, function(v) if (any(!is.na(v))) v[!is.na(v)][1] else NA_character_)
        unname(first[rownames(gt)])
      } else rep("all", nrow(gt))
      names(grp) <- rownames(gt)
      grp[is.na(grp) | !nzchar(grp)] <- "unassigned"
      f <- gid_filter(gt,
             min_locus_call  = if (is.na(input$min_locus_call))  0 else input$min_locus_call,
             min_sample_call = if (is.na(input$min_sample_call)) 0 else input$min_sample_call)
      validate(need(nrow(f$gt) >= 2,
        "Fewer than two samples survive the call-rate filters. Lower them in the sidebar."),
        need(ncol(f$gt) >= 1, "No loci survive the call-rate filter."))
      ## Every reaction is kept whichever route the tabs show, because both
      ## routes run: reps is the route on screen, reps_all is always there.
      reps_all <- list(gt = rgt[, colnames(f$gt), drop = FALSE], sample = rsm)
      mx <- mix_step(f$gt)
      return(list(gt = mx$gt, grp = grp[rownames(mx$gt)], raw_gt = gt,
                  reps = if (use_reps) reps_all else NULL, reps_all = reps_all,
                  n_rep = n_rep, hom_n = hom_n, het_n = het_n,
                  dropped_loci = f$dropped_loci, dropped_samples = f$dropped_samples,
                  mixture = mx$mixture, dropped_mixtures = mx$dropped,
                  df = df, n_reps = length(rsm) / nrow(gt)))
    }

    if (anyDuplicated(ids))
      ids <- ave(ids, ids, FUN = function(z)
        if (length(z) == 1) z else paste0(z, "#", seq_along(z)))
    df$.gid_key <- ids
    gt <- gid_matrix(df, ".gid_key", loci(), sep = sep,
                     strip_flags = !isTRUE(input$drop_flagged))
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
    mx <- mix_step(f$gt)
    list(gt = mx$gt, grp = grp[rownames(mx$gt)], raw_gt = gt, reps = NULL,
         reps_all = NULL,
         dropped_loci = f$dropped_loci, dropped_samples = f$dropped_samples,
         mixture = mx$mixture, dropped_mixtures = mx$dropped, df = df)
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

    ## p$raw_gt has one row per SAMPLE; p$df has one row per REACTION on a
    ## replicate file, so the two cannot be compared cell for cell. Normalize the
    ## uploaded rows themselves instead -- gid_norm_gt() is column-wise and keeps
    ## every row, so the shapes always match.
    lcols  <- colnames(p$raw_gt)
    sep_qc <- gid_guess_sep(p$df, exclude = input$id_col)
    norm_m <- vapply(lcols, function(L) gid_norm_gt(p$df[[L]], sep = sep_qc),
                     character(nrow(p$df)))
    if (is.null(dim(norm_m))) norm_m <- matrix(norm_m, nrow = nrow(p$df))
    nre <- sum(gsub("/", "", norm_m) !=
                 toupper(gsub("[*?!#]", "", as.matrix(p$df[, lcols]))), na.rm = TRUE)
    if (nre > 0) f <- c(f, list(bad(sprintf(
      "%d genotype cells had their two alleles written in the other order (AG vs GA). They are now normalized; compared as raw text they would have counted as mismatches.", nre))))

    flagged <- grepl("[*?!#]", as.matrix(p$df[, lcols]))
    nflag <- sum(flagged)
    ## A flag that only ever lands on heterozygotes usually marks an allele
    ## ratio outside the heterozygote range -- exactly the call a false allele
    ## produces -- so it is worth naming.
    het_only <- nflag >= 20 && all(vapply(which(flagged), function(k) {
      a <- strsplit(norm_m[k], "/", fixed = TRUE)[[1]]
      length(a) == 2 && a[1] != a[2] }, TRUE))
    if (nflag > 0) f <- c(f, list(bad(sprintf(
      "%d cells carry a quality flag character (%s).%s %s", nflag,
      if (isTRUE(input$drop_flagged)) "discarded as missing"
      else "the flag was stripped and the call kept",
      if (het_only) " Every one of them is a heterozygote, which usually means the flag marks an allele ratio outside the heterozygote range \u2014 the call a false allele makes." else "",
      if (isTRUE(input$drop_flagged))
        "Uncheck \"Treat quality-flagged calls as missing\" in the sidebar to keep them instead."
      else "If your pipeline uses the flag to mean \"low confidence\", check \"Treat quality-flagged calls as missing\" in the sidebar and re-run to see how much it moves the answer."))))

    nblank <- sum(gid_blank_ids(raw(), input$id_col))
    if (nblank > 0) {
      alt <- setdiff(names(raw()), c(input$id_col, colnames(p$raw_gt)))
      alt <- alt[vapply(alt, function(k)
        !any(gid_blank_ids(raw(), k)) && anyDuplicated(as.character(raw()[[k]])) > 0, TRUE)]
      f <- c(f, list(bad(tagList(
        sprintf("%d rows have a blank %s and were dropped. ", nblank, input$id_col),
        "They cannot be looked up by name, and left in they merge into a single ",
        "fake sample whose replicates are different animals \u2014 which corrupts ",
        "the estimated error rates.",
        if (length(alt)) tagList(" A column with no blanks that still repeats: ",
                                 tags$b(alt[1]),
                                 ". Switch to it if that is your sample identifier.")))))
    }

    ## On a replicate file the repeated ids ARE the replicates and nothing is
    ## suffixed, so the warning below would be both wrong and alarming.
    on_reps <- isTRUE(has_reps()) && nzchar(input$rep_col %||% "")
    dupn <- sum(duplicated(as.character(p$df[[input$id_col]])))
    if (dupn > 0 && !on_reps) f <- c(f, list(bad(sprintf(
      "%d duplicated sample IDs. They were kept and numbered #1, #2, and so on. If these are the same extract run twice, they are a useful positive control: any correct method must put them together.", dupn))))
    if (on_reps) f <- c(f, list(ok(sprintf(
      "%d rows for %d samples, read as replicate reactions from the %s column. %s",
      nrow(p$df), nrow(p$gt), input$rep_col,
      if (identical(input$rep_mode, "reps"))
        "The tabs show every reaction used directly; the consensus route runs alongside."
      else "The tabs show consensus genotypes; the route using every reaction runs alongside."))))

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

  ## ---- samples that may be two animals mixed -------------------------------
  output$mix_summary <- renderUI({
    p <- prep(); req(p)
    m <- p$mixture
    if (is.null(m)) return(tags$p(class = "gid-hint",
      "Too few samples to test. The check needs at least ten."))
    n_test <- sum(m$flag != "not tested")
    n_bad  <- sum(m$flag == "possible mixture")
    n_skip <- sum(m$flag == "not tested")
    skip_note <- if (n_skip) tags$div(class = "gid-hint", style = "margin-top:.4rem",
      sprintf(paste("%d samples were not tested: they are called at under 75%% of the",
                    "panel, and a patchy genotype has a noisy heterozygosity that clears",
                    "any fixed threshold by chance. Eriksson et al. required 77%% before",
                    "analyzing anything, for this reason."), n_skip))
    if (!n_bad) return(tags$div(class = "gid-flag gid-ok",
      tags$b("No sample looks like a mixture. "),
      sprintf("Across %d tested samples the largest heterozygosity excess is %.3f, under the %s threshold.",
              n_test, max(m$excess, na.rm = TRUE), input$mix_cut),
      skip_note))
    tags$div(
      class = "gid-flag",
      tags$b(sprintf("%d of %d samples carry more heterozygous loci than one animal should. ",
                     n_bad, n_test)),
      "Each is either a genuinely heterozygous animal or two animals in one ",
      "sample. The second kind is the dangerous one: it becomes an individual ",
      "that does not exist and inflates your count.",
      tags$div(style = "margin-top:.4rem",
        if (length(p$dropped_mixtures))
          sprintf("These %d are currently excluded from the analysis.", length(p$dropped_mixtures))
        else tagList("They are ", tags$b("included"), " right now. Check ",
                     tags$b("Exclude samples that look like two animals"),
                     " in the sidebar to drop them and see how much the count moves.")),
      tags$div(class = "gid-hint", style = "margin-top:.4rem",
        sprintf("%d flagged by the published 0.18 excess threshold, %d by a test that also accounts for how many loci the sample has. ",
                sum(m$by_excess, na.rm = TRUE), sum(m$by_test, na.rm = TRUE)),
        "The test is the stricter of the two. An even two-way mixture is the ",
        "case that slips through either way: on simulated mixtures this caught ",
        "17 of 20, and the misses were mostly lopsided ones that simply ",
        "genotype as the dominant animal and so introduce nobody new."),
      skip_note)
  })

  output$tbl_mixture <- renderDT({
    p <- prep(); req(p); m <- p$mixture
    if (is.null(m)) return(dt(data.frame(message = "Not enough samples to test.")))
    x <- m[m$flag != "not tested",
           c("sample", "n_loci", "Ho", "He", "excess", "q", "by_test", "by_excess", "flag")]
    x$Ho <- round(x$Ho, 3); x$He <- round(x$He, 3)
    x$excess <- round(x$excess, 3); x$q <- signif(x$q, 3)
    dt(x)
  })

  output$plot_mixture <- renderPlot({
    p <- prep(); req(p); m <- p$mixture; req(m)
    x <- m[!is.na(m$excess), ]
    cut <- if (is.na(input$mix_cut %||% NA)) NA_real_ else input$mix_cut
    ggplot(x, aes(excess, fill = flag)) +
      geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
      {if (!is.na(cut)) geom_vline(xintercept = cut, linetype = 2, colour = ACCENT)} +
      scale_fill_manual(values = c("single individual" = INK,
                                   "possible mixture" = ACCENT), name = NULL) +
      labs(x = "Heterozygosity excess (observed minus expected)", y = "Samples",
           subtitle = paste("One animal sits near zero. A blend of two sits to the right,",
                            "\nbecause the two differ at many loci and every difference reads heterozygous.")) +
      theme_gid()
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

  ## The residual rates a consensus call is left with, given per-reaction ones.
  to_consensus <- function(p, r) {
    pe <- gid_propagate_error(r[["dropout"]], r[["false_allele"]], n_rep = p$n_rep,
                              rule = "taberlet", hom_n = p$hom_n, het_n = p$het_n,
                              per_rep_missing = mean(is.na(p$reps_all$gt)))
    c(dropout = max(pe[["dropout"]], GID_CONSENSUS_FLOOR),
      false_allele = max(pe[["false_allele"]], GID_CONSENSUS_FLOOR))
  }

  ## The error rates each route runs on. The sliders hold the rate for the
  ## route the tabs show, and the other route's rate follows from it or from
  ## the replicates. Sliders nobody has touched are replaced by rates measured
  ## from the replicates, so someone who never opens the error panel still
  ## runs on their own data's rates rather than on defaults.
  route_rates <- function(p, progress = NULL) {
    slid <- c(dropout = input$dropout, false_allele = input$false_allele)
    if (is.null(p$reps_all))
      return(list(primary = "consensus", consensus = slid, reps = NULL,
                  source = c(consensus = "set in the sidebar")))
    primary <- if (is.null(p$reps)) "consensus" else "reps"
    manual  <- !rates_untouched()
    e <- err_est()
    if (!(manual && primary == "reps") && is.null(measured_for(e, "reps"))) {
      e <- tryCatch(measure_error(p, progress), error = function(err) NULL)
      if (!is.null(e)) err_est(c(e, list(auto = TRUE)))
    }
    R <- if (manual && primary == "reps") slid else measured_for(e, "reps")
    C <- if (manual && primary == "consensus") slid
         else if (!is.null(R)) to_consensus(p, R)
    ## too few replicated samples to fit: the sliders are all there is
    if (is.null(R) && primary == "reps") R <- slid
    if (is.null(C)) C <- slid
    if (!manual) {
      v <- if (primary == "reps") R else C
      set_rates(v[["dropout"]], v[["false_allele"]])
    }
    meas <- !is.null(measured_for(e, "reps"))
    list(primary = primary, reps = R, consensus = C,
         source = c(
           reps = if (manual && primary == "reps") "set in the sidebar"
                  else if (meas) "measured from the replicates"
                  else "the defaults, since too few samples were replicated to measure them",
           consensus = if (manual && primary == "consensus") "set in the sidebar"
                       else if (manual) "the sidebar's per-reaction rates, through the multi-tube rule"
                       else if (meas) "the measured rates, through the multi-tube rule"
                       else "the defaults"))
  }

  ## Not ignoreInit: Shiny skips an eventReactive's first evaluation, not its
  ## first event, so anything read for the first time after a run would sit
  ## blank until the next one. Waiting on the counter does the same job safely.
  res <- eventReactive(run_count(), {
    req(run_count() > 0)
    p <- prep(); req(p, nrow(p$gt) >= 2)
    run_error(NULL)
    out <- tryCatch(withProgress(message = "Identifying individuals", value = 0, {
      rates <- route_rates(p, progress = function(v, d)
        setProgress(value = 0.3 * v, detail = paste("measuring error:", d)))
      pr <- if (is.na(input$prior)) NULL else input$prior

      ## Sethi and the likelihood ratio on one route. Exact matching, the
      ## mismatch threshold and GenAlEx compare one genotype per sample, so
      ## they always run on the consensus.
      likelihood_methods <- function(route) {
        rr <- rates[[route]]
        rp <- if (route == "reps") p$reps_all else NULL
        list(
          sethi = if (length(input$sethi_rel))
            gid_by_group(p$gt, p$grp, gid_method_sethi,
                         dropout = rr[["dropout"]], false_allele = rr[["false_allele"]],
                         relationships = input$sethi_rel, lambda_cut = input$lambda_cut,
                         min_loci = input$min_loci, linkage = input$linkage, reps = rp),
          probabilistic = gid_by_group(p$gt, p$grp, gid_method_lr,
                         dropout = rr[["dropout"]], false_allele = rr[["false_allele"]],
                         kinship = input$kinship, post_cut = input$post_cut,
                         min_loci = input$min_loci, linkage = input$linkage,
                         prior_same = pr, reps = rp))
      }
      setProgress(0.32, detail = "Sethi and the likelihood ratio")
      routes <- list()
      routes[[rates$primary]] <- likelihood_methods(rates$primary)
      ## With replicates in the file, the route not on screen runs too.
      other <- setdiff(c("reps", "consensus"), rates$primary)
      if (!is.null(p$reps_all) && !is.null(rates[[other]])) {
        setProgress(0.5, detail = if (other == "reps") "every replicate reaction"
                                  else "consensus genotypes")
        routes[[other]] <- likelihood_methods(other)
      }
      lr <- routes[[rates$primary]]$probabilistic

      setProgress(0.68, detail = "exact match")
      ex <- gid_by_group(p$gt, p$grp, gid_method_exact,
                         min_loci = input$min_loci, linkage = input$linkage)
      setProgress(0.75, detail = "mismatch threshold")
      th <- gid_by_group(p$gt, p$grp, gid_method_threshold,
                         max_mismatch = input$max_mismatch, min_loci = input$min_loci,
                         linkage = input$linkage)
      ## Spatial evidence, if asked for and if the file carries coordinates.
      ## Applied only to the likelihood-ratio result: it is the one method with
      ## a per-pair posterior for the distance term to multiply into.
      spatial_note <- NULL
      if (isTRUE(input$use_spatial)) {
        cc <- try(map_api$coords(), silent = TRUE)
        if (!inherits(cc, "try-error") && !is.null(cc) && nrow(cc) > 1) {
          setProgress(0.8, detail = "spatial evidence")
          lr2 <- try(gid_apply_spatial(lr, cc, post_cut = input$post_cut,
                                       linkage = input$linkage), silent = TRUE)
          if (!inherits(lr2, "try-error")) lr <- lr2
          else spatial_note <- "Spatial evidence could not be calibrated; genetics only."
        } else spatial_note <- "No usable coordinates, so spatial evidence was skipped."
      }
      out <- list(exact = ex, threshold = th, probabilistic = lr,
                  sethi = routes[[rates$primary]]$sethi)
      if (isTRUE(input$run_genalex)) {
        setProgress(0.85, detail = "GenAlEx Matches")
        out$genalex <- gid_by_group(p$gt, p$grp, gid_method_genalex,
                                    min_loci = input$min_loci,
                                    near_match_loci = input$near_match,
                                    linkage = input$linkage)
      }
      out <- Filter(Negate(is.null), out)

      ## Where the two routes part ways, for each likelihood method. Compared
      ## on genetics alone, since spatial evidence is added to one route only.
      route_cmp <- NULL
      if (length(routes) == 2) {
        keys <- intersect(names(Filter(Negate(is.null), routes$reps)),
                          names(Filter(Negate(is.null), routes$consensus)))
        route_cmp <- setNames(lapply(keys, function(k)
          gid_compare_routes(routes$reps[[k]], routes$consensus[[k]])), keys)
      }
      ## how much more often each sample errs than typical, as the replicate
      ## route measured it
      own_error <- NULL
      src <- routes$reps$sethi %||% routes$reps$probabilistic
      if (!is.null(src))
        own_error <- unlist(unname(lapply(src$by_group, function(e) e$sample_error)))

      setProgress(0.95, detail = "comparing")
      ord <- c("exact", "threshold", "genalex", "sethi", "probabilistic")
      list(methods = out[ord[ord %in% names(out)]], spatial_note = spatial_note,
           cmp = gid_compare_methods(out), gt = p$gt, grp = p$grp,
           rates = rates, routes = routes, route_cmp = route_cmp,
           own_error = own_error,
           sweep = gid_threshold_sweep(p$gt, max_k = 6, min_loci = input$min_loci,
                                       linkage = input$linkage))
    }), error = function(e) { run_error(conditionMessage(e)); NULL })
    out
  })

  ## the method the Individuals tab is showing
  best <- reactive({
    r <- res(); req(r)
    m <- r$methods[[input$method]]
    validate(need(!is.null(m), switch(input$method,
      sethi = "Sethi needs at least one competing relationship. Choose one in the sidebar, then press Identify individuals.",
      sprintf("%s was not run. Check it under \"Which methods to run\" in the sidebar, then press Identify individuals.",
              GID_METHODS[[input$method]]$label))))
    m
  })

  meta <- reactive(GID_METHODS[[input$method %||% "sethi"]])

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
  output$vb_persample <- renderText({
    s <- gid_summarise_assignment(best()$assignment)
    sprintf("%.0f / %.1f", s$median_samples, s$mean_samples) })

  output$recap_detail <- renderUI({
    a <- best()$assignment; req(a)
    s <- gid_summarise_assignment(a)
    n <- as.integer(table(a$individual))
    tags$p(class = "gid-hint", style = "margin:-.4rem 0 1rem",
      tags$b("Samples per animal: "),
      sprintf("median %g, mean %.2f, range %d to %d. ",
              s$median_samples, s$mean_samples, min(n), max(n)),
      tags$b("Recaptures per animal "), "(samples after the first): ",
      sprintf("median %g, mean %.2f. ", s$median_recaptures, s$mean_recaptures),
      sprintf("Overall recapture rate %.0f%%: that share of your samples were repeats of an animal already seen. ",
              100 * s$recapture_rate),
      if (s$mean_samples > 1.5 * max(s$median_samples, 1))
        tags$span(style = "color:#c1502e",
          sprintf("The mean is well above the median, so a few animals dominate: the top animal alone accounts for %d of %d samples.",
                  max(n), s$n_samples)))
  })
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
    if (run_count() == 0)
      return(tags$div(class = "gid-flag gid-ok", style = "margin-bottom:1rem",
        tags$b("No results yet. "),
        "Choose your sample ID column and loci in the sidebar, then press ",
        tags$b("Identify individuals"), " to compute. Nothing on this tab is ",
        "calculated until you do \u2014 the other tabs update on their own."))
    NULL
  }
  ## ---- what to change before the next run ---------------------------------
  ## Enough samples typed more than once to measure error from: a file can
  ## repeat its sample ids yet hold one usable reaction each once rows are
  ## excluded.
  replicated <- reactive({
    p <- tryCatch(prep(), error = function(e) NULL)
    !is.null(p$reps_all) && sum(table(p$reps_all$sample) >= 2) >= 5
  })
  advice <- reactive({
    r <- res(); req(r)
    cal <- tryCatch({
      tb <- calib()
      g <- tb$group[which.max(tb$n_ids)]
      x <- tb[tb$group == g, ]
      attr(x, "max_posterior") <- x$max_post[1]
      x
    }, error = function(e) NULL)
    ## Rates in force and rates measured, both on the scale of the route on
    ## screen, or the comparison is between a per-reaction and a residual rate.
    ## Only a measurement from replicates is worth holding the sliders to; the
    ## heterozygote-deficit figure is an order of magnitude, not a target.
    inforce <- r$rates[[r$rates$primary]]
    e <- err_est()
    m <- measured_for(e, r$rates$primary)
    e <- if (is.null(m)) NULL
         else modifyList(e, list(dropout = m[["dropout"]], false_allele = m[["false_allele"]]))
    gid_advise(
      r,
      settings = list(post_cut = input$post_cut, dropout = inforce[["dropout"]],
                      false_allele = inforce[["false_allele"]], kinship = input$kinship,
                      min_loci = input$min_loci, linkage = input$linkage,
                      lambda_cut = input$lambda_cut,
                      budget = input$err_budget %||% 1),
      err_measured = e,
      calib = cal,
      pid = tryCatch(pid(), error = function(e) NULL),
      has_reps = replicated(),
      shown = input$method,
      routes = r$route_cmp[[input$method]])
  })

  output$advice_card <- renderUI({
    a <- advice()
    if (is.null(a) || !nrow(a))
      return(tags$div(class = "gid-flag gid-ok", style = "margin-bottom:1rem",
        tags$b("Nothing to change. "),
        "Every setting is consistent with what this data supports: ",
        if (identical(input$method, "probabilistic"))
          "the cutoff is reachable and on a plateau, ",
        if (replicated()) "the error rates match your replicates, ",
        "and no cluster is held together by a chain."))
    chip <- function(sev) tags$span(
      class = paste0("gid-sev gid-sev-", sev),
      switch(sev, high = "change this", medium = "worth changing", "consider"))
    ## Folded shut unless something wants changing urgently: open, this card is
    ## nearly 900px tall and buries the results it is commenting on.
    urgent <- any(a$severity == "high")
    tags$details(
      class = "gid-advice", open = if (urgent) NA else NULL,
      tags$summary(
        class = "gid-advice-head",
        tags$b(sprintf("%d thing%s to consider before your next run", nrow(a),
                       if (nrow(a) == 1) "" else "s")),
        tags$span(class = "gid-hint", style = "margin-left:.5rem",
                  paste(unique(as.character(a$severity)), collapse = ", ")),
        tags$span(class = "gid-hint", style = "margin-left:.5rem",
                  "\u2014 click to open")),
      tags$div(class = "gid-hint", style = "margin:.2rem 0 .5rem",
               "Each one is drawn from your own data. Applying a change does ",
               "not re-run anything \u2014 press Identify individuals when ready."),
      lapply(seq_len(nrow(a)), function(i) tags$div(
        class = "gid-advice-row",
        tags$div(chip(as.character(a$severity[i]))),
        tags$div(
          tags$div(class = "gid-advice-title", a$title[i]),
          tags$div(class = "gid-advice-why", a$why[i]),
          if (!is.na(a$evidence[i]))
            tags$div(class = "gid-advice-ev", a$evidence[i]),
          if (!is.na(a$label[i])) tags$div(
            style = "margin-top:.4rem",
            actionButton(paste0("adv_", a$id[i]), a$label[i],
                         class = "btn-sm btn-outline-primary"))))))
  })

  ## One handler per suggestion. They are a fixed vocabulary, so this stays
  ## explicit rather than dispatching on a string.
  apply_advice <- function(id) {
    a <- advice(); req(a)
    row <- a[a$id == id, ][1, ]
    switch(row$param,
      post_cut   = updateSliderInput(session, "post_cut", value = row$value),
      kinship    = updateSelectInput(session, "kinship", selected = "full_sib"),
      linkage    = updateRadioButtons(session, "linkage", selected = "complete"),
      min_loci   = updateNumericInput(session, "min_loci", value = row$value),
      lambda_cut = updateNumericInput(session, "lambda_cut", value = row$value),
      error_rates = {
        m <- measured_for(err_est(), if (identical(input$rep_mode, "consensus")) "consensus" else "reps")
        req(m)
        set_rates(m[["dropout"]], m[["false_allele"]])
      })
    showNotification("Setting changed. Press Identify individuals to re-run.",
                     type = "message", duration = 6)
  }
  for (.id in c("cutoff_unreachable", "cutoff_plateau", "error_rates",
                "kinship", "linkage", "min_loci", "sethi_gap"))
    local({
      id <- .id
      observeEvent(input[[paste0("adv_", id)]], apply_advice(id), ignoreInit = TRUE)
    })

  ## ---- why the two recommended methods disagree ---------------------------
  lr_sethi <- reactive({
    r <- res(); req(r)
    lr <- r$methods$probabilistic; se <- r$methods$sethi
    req(lr, se)
    if (!is.null(lr$by_group)) {
      do.call(rbind, lapply(names(lr$by_group), function(g) {
        x <- gid_explain_disagreement(lr$by_group[[g]], se$by_group[[g]],
                                      post_cut = input$post_cut,
                                      lambda_cut = input$lambda_cut,
                                      min_loci = input$min_loci)
        if (is.null(x)) NULL else cbind(group = g, x)
      }))
    } else gid_explain_disagreement(lr, se, post_cut = input$post_cut,
                                    lambda_cut = input$lambda_cut,
                                    min_loci = input$min_loci)
  })

  output$disagree_summary <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL); req(r)
    if (is.null(r$methods$sethi))
      return(tags$div(class = "gid-flag",
        "Sethi was not run, so there is nothing to compare. Choose at least one ",
        "competing relationship for it in the sidebar."))
    d <- tryCatch(lr_sethi(), error = function(e) NULL)
    n_lr <- length(unique(r$methods$probabilistic$assignment$individual))
    n_se <- length(unique(r$methods$sethi$assignment$individual))
    if (is.null(d) || !nrow(d))
      return(tags$div(class = "gid-flag gid-ok",
        tags$b(sprintf("The two agree on every pair: %d individuals each. ", n_lr)),
        "Where the recommended methods agree, the answer does not depend on ",
        "which decision rule you prefer, and that is worth saying in a paper."))
    rule <- sum(grepl("^Decision rule", d$cause))
    alt  <- nrow(d) - rule
    tags$div(class = "gid-flag",
      tags$b(sprintf("%d pair%s decided differently. Likelihood ratio: %d individuals. Sethi: %d.",
                     nrow(d), if (nrow(d) == 1) "" else "s", n_lr, n_se)),
      tags$ul(style = "margin:.45rem 0 0;padding-left:1.1rem",
        if (rule) tags$li(sprintf(
          "%d down to the decision rule: the same evidence, judged against a different bar. Argue about the threshold, not the model.",
          rule)),
        if (alt) tags$li(sprintf(
          "%d down to the alternative hypothesis: the two are dividing by different relationships, so one of them is asking the wrong question of these pairs.",
          alt))),
      if (rule && !alt) tags$div(style = "margin-top:.4rem", class = "gid-hint",
        "Every difference here is the threshold. The models are not in conflict."))
  })

  output$tbl_lr_sethi <- renderDT({
    d <- lr_sethi()
    if (is.null(d) || !nrow(d))
      return(dt(data.frame(message = "The two methods agree on every pair.")))
    dt(d)
  })

  ## ---- every replicate, or a consensus first ------------------------------
  route_label <- function(rt) if (rt == "reps") "every replicate reaction" else "consensus first"

  output$route_summary <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL); req(r)
    if (length(r$routes) < 2) return(NULL)
    rc <- r$route_cmp[[input$method]]
    if (is.null(rc))
      return(tags$p(class = "gid-hint", style = "margin:-.4rem 0 1rem",
        "This method compares one genotype per sample, so it runs on the ",
        "consensus whichever route you chose. Sethi and the likelihood ratio ",
        "run on both; the Method comparison tab sets them side by side."))
    n <- rc$n_individuals
    shown <- if (r$rates$primary == "reps")
      "every replicate reaction, with each sample's own error rate"
      else "one consensus genotype per sample"
    if (rc$n_moved == 0)
      return(tags$div(class = "gid-flag gid-ok", style = "margin-bottom:1rem",
        tags$b(sprintf("Both routes agree: %d individuals ", n[["replicates"]]),
               "whether every replicate is used or a consensus is taken first. "),
        sprintf("The results below use %s.", shown)))
    tags$div(class = "gid-flag", style = "margin-bottom:1rem",
      tags$b(sprintf("Every replicate gives %d individuals; a consensus first gives %d. ",
                     n[["replicates"]], n[["consensus"]])),
      sprintf("The results below use %s. ", shown),
      sprintf("%d sample%s placed differently by the two routes",
              rc$n_moved, if (rc$n_moved == 1) " is" else "s are"),
      " — they are listed on the ", tags$b("Method comparison"), " tab, and ",
      "they are the ones to check by hand before you report a count.",
      if (identical(input$method, "probabilistic") && isTRUE(input$use_spatial))
        " Both counts are from genetics alone.")
  })

  output$route_note <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL); req(r)
    rr <- r$rates
    if (length(r$routes) < 2)
      return(tags$div(class = "gid-flag",
        if (is.null(rr$reps))
          paste("The replicate route was not run: too few samples have two or",
                "more reactions to measure the per-reaction error rates it needs.")
        else "Only one route was run."))
    oe <- r$own_error
    n_hi <- sum(oe > 1, na.rm = TRUE)
    tags$p(class = "gid-hint",
      tags$b("Every replicate reaction: "),
      sprintf("dropout %s and false allele %s per reaction, %s. ",
              gid_fmt_rate(rr$reps[["dropout"]]), gid_fmt_rate(rr$reps[["false_allele"]]),
              rr$source[["reps"]]),
      if (length(oe)) sprintf("%d of %d samples err more often than that, by up to %s times, and %s given rates of %s own. ",
                              n_hi, length(oe), format(max(oe, na.rm = TRUE)),
                              if (n_hi == 1) "was" else "were", if (n_hi == 1) "its" else "their"),
      tags$b("Consensus first: "),
      sprintf("dropout %s and false allele %s per consensus call, from %s.",
              gid_fmt_rate(rr$consensus[["dropout"]]), gid_fmt_rate(rr$consensus[["false_allele"]]),
              rr$source[["consensus"]]))
  })

  output$tbl_routes <- renderDT({
    r <- res(); req(length(r$routes) == 2)
    x <- do.call(rbind, lapply(names(r$route_cmp), function(k) {
      rc <- r$route_cmp[[k]]
      do.call(rbind, lapply(c("reps", "consensus"), function(rt) {
        s <- gid_summarise_assignment(r$routes[[rt]][[k]]$assignment)
        data.frame(method = k, route = route_label(rt), shown = rt == r$rates$primary,
                   dropout = signif(r$rates[[rt]][["dropout"]], 3),
                   false_allele = signif(r$rates[[rt]][["false_allele"]], 3),
                   n_individuals = s$n_individuals, n_singletons = s$n_singletons,
                   agreement_between_routes = round(rc$ari, 3),
                   samples_placed_differently = rc$n_moved,
                   stringsAsFactors = FALSE)
      }))
    }))
    dt(x)
  })

  route_moved <- reactive({
    r <- res(); req(length(r$routes) == 2)
    x <- do.call(rbind, lapply(names(r$route_cmp), function(k) {
      mv <- r$route_cmp[[k]]$moved
      if (!nrow(mv)) return(NULL)
      if (!is.null(r$own_error)) mv$own_error <- unname(r$own_error[mv$sample])
      cbind(method = k, mv, stringsAsFactors = FALSE)
    }))
    if (is.null(x)) data.frame() else x
  })

  output$tbl_route_moved <- renderDT({
    x <- route_moved()
    if (!nrow(x))
      return(dt(data.frame(message = "The two routes group every sample the same way.")))
    dt(x)
  })

  ## ---- how consistent each sample's replicates are ------------------------
  rep_quality <- reactive({
    p <- prep(); req(p, p$reps_all)
    d <- gid_replicate_discordance(p$reps_all$gt, p$reps_all$sample)
    d <- d[d$sample %in% rownames(p$gt), , drop = FALSE]
    r <- tryCatch(res(), error = function(e) NULL)
    if (!is.null(r$own_error)) d$own_error <- unname(r$own_error[d$sample])
    d
  })

  output$rep_quality_note <- renderUI({
    d <- rep_quality()
    q <- gid_quality_spread(d)
    if (is.null(q))
      return(tags$p(class = "gid-hint",
        "Too few samples have loci typed twice or more to judge the spread."))
    r <- tryCatch(res(), error = function(e) NULL)
    own <- if (!is.null(d$own_error)) {
      n_hi <- sum(d$own_error > 1, na.rm = TRUE)
      sprintf(" On the last run, %d sample%s given error rates of %s own, above the dataset's.",
              n_hi, if (n_hi == 1) " was" else "s were", if (n_hi == 1) "its" else "their")
    }
    if (q$uneven)
      tags$div(class = "gid-flag", style = "margin-bottom:.6rem",
        tags$b(sprintf("%d of %d samples (%.0f%%) are much noisier than the rest. ",
                       q$n_noisy, q$n, 100 * q$share_noisy)),
        sprintf(paste("Their reactions disagree at more than %.0f%% of their replicated loci,",
                      "against %.0f%% for the typical sample. One error rate for the whole",
                      "dataset would under-describe them and split their recaptures off as",
                      "new animals, so the replicate route gives every sample its own",
                      "error rates, and the consensus route is reported alongside."),
                100 * q$cut, 100 * q$typical), own)
    else
      tags$div(class = "gid-flag gid-ok", style = "margin-bottom:.6rem",
        tags$b("The samples behave alike. "),
        sprintf(paste("The typical sample's reactions disagree at %.0f%% of its replicated",
                      "loci, and %d of %d exceed %.0f%%, so one error rate describes",
                      "this dataset well."),
                100 * q$typical, q$n_noisy, q$n, 100 * q$cut), own)
  })

  output$plot_rep_quality <- renderPlot({
    d <- rep_quality(); d <- d[!is.na(d$discordance), , drop = FALSE]; req(nrow(d))
    q <- gid_quality_spread(d)
    g <- ggplot(d, aes(discordance)) +
      geom_histogram(binwidth = 0.02, boundary = 0, fill = INK, colour = "white",
                     linewidth = 0.3) +
      labs(x = "Share of a sample's replicated loci whose reactions disagree",
           y = "Samples",
           subtitle = if (!is.null(q))
             "Dashed line: the point past which a sample counts as much noisier than typical"
           else NULL) +
      theme_gid()
    if (!is.null(q)) g <- g + geom_vline(xintercept = q$cut, colour = ACCENT, linetype = 2)
    g
  })

  output$tbl_rep_quality <- renderDT({
    d <- rep_quality()
    d$discordance <- round(d$discordance, 3)
    dt(d[order(-d$discordance, -d$loci_discordant), , drop = FALSE])
  })

  output$run_status_ind <- renderUI(run_status())

  ## What the spatial term actually did, so it is judged rather than trusted.
  output$spatial_diag <- renderUI({
    r <- tryCatch(res(), error = function(e) NULL); req(r)
    if (!is.null(r$spatial_note))
      return(tags$div(class = "gid-flag", tags$b("Spatial evidence: "), r$spatial_note))
    sp <- r$methods$probabilistic$spatial
    if (is.null(sp) || !identical(input$method, "probabilistic")) return(NULL)
    if (!is.list(sp) || is.null(names(sp))) sp <- list(all = sp)
    rows <- Filter(Negate(is.null), lapply(names(sp), function(g) {
      a <- sp[[g]]
      if (is.null(a)) return(NULL)
      if (!isTRUE(a$applied))
        return(tags$li(sprintf("%s: not applied \u2014 %s", g, a$note %||% "")))
      tags$li(sprintf(
        "%s: home-range scale %.1f km, from %d matched pairs. Matched pairs sit %.1f km apart on average, unrelated pairs %.1f km. %s",
        g, a$sigma / 1000, a$n_matched, a$median_matched_m / 1000,
        a$median_ruled_out_m / 1000,
        if (a$median_ruled_out_m > 3 * a$median_matched_m)
          "Well separated, so location is informative here."
        else
          "Not well separated \u2014 home ranges overlap too much for location to add much, and it may mislead."))
    }))
    if (!length(rows)) return(NULL)
    tags$div(class = "gid-flag gid-ok", style = "margin-bottom:1rem",
             tags$b("Spatial evidence is switched on. "),
             tags$ul(style = "margin:.4rem 0 0;padding-left:1.1rem", rows))
  })
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
      sprintf(paste0("The most similar pair in this dataset reaches a posterior of only %s ",
                     "(log10 LR %.2f) against a %s alternative, because %d loci carry a finite ",
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
    ## A method that clusters internally exposes no pairwise table. Fall back
    ## to plain locus-mismatch counts, which are a property of the data rather
    ## than of any method.
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
                  note = paste("This method clusters internally and does not expose a",
                               "score per pair, so this shows how many loci differ",
                               "between each pair of samples \u2014 a property of your",
                               "data, not of the method.")))
    if (!is.null(pp$log10_LR))
      list(v = pp$log10_LR, cut = suppressWarnings(min(pp$log10_LR[pp$posterior_same >= input$post_cut])),
           x = "log10 likelihood ratio", kind = "lr",
           note = sprintf("Positive values favor one animal. Alternative hypothesis: %s.",
                          gsub("_", " ", input$kinship)))
    else if (!is.null(pp$log10_lambda))
      list(v = pp$log10_lambda, cut = log10(input$lambda_cut), x = "log10 lambda", kind = "lr",
           note = "Positive values favor one animal over the best competing relationship.")
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

  ## ---- sample -> individual map with per-sample confidence ----------------
  conf <- reactive({
    r <- res(); req(r); req(best())
    gid_method_confidence(r, input$method, pid_tab = pid(),
                          post_cut = input$post_cut, min_loci = input$min_loci)
  })

  map_tbl <- reactive({
    cf <- conf(); req(cf)
    a <- best()$assignment
    cf$sex <- if (!is.null(res()$df) && nzchar(input$group_col %||% "")) NULL else NULL
    keep <- c("sample", "individual", "group", "n_in_individual", "held_by",
              "n_loci", "support", "rival", "rival_sample", "rival_individual",
              "margin", "status")
    cf[, intersect(keep, names(cf)), drop = FALSE]
  })

  output$run_status_map <- renderUI(run_status())
  output$vb_map_ind <- renderText({ length(unique(best()$assignment$individual)) })
  output$vb_map_ok  <- renderText({ sum(conf()$status == "Confident") })
  output$vb_map_chk <- renderText({ sum(conf()$status %in% c("Check", "Underpowered")) })
  output$vb_map_bad <- renderText({ sum(conf()$status == "Uncertain") })

  output$map_review_note <- renderUI({
    cf <- conf(); req(cf)
    n <- sum(cf$status != "Confident")
    if (!n) return(tags$div(class = "gid-flag gid-ok",
      tags$b("Nothing to review. "), "Every sample sits well clear of the cutoff, ",
      "so no assignment here turns on where you drew the line."))
    tags$div(class = "gid-flag",
      tags$b(sprintf("%d of %d samples are worth a look. ", n, nrow(cf))),
      "These are the ones whose assignment would change under a slightly ",
      "different cutoff. Nothing here is necessarily wrong \u2014 it is where ",
      "your answer is soft, and where a second line of evidence (collection date, ",
      "location, sex, a repeat extraction) earns its keep.")
  })

  output$tbl_review <- renderDT({
    x <- map_tbl(); x <- x[x$status != "Confident", , drop = FALSE]
    if (!nrow(x)) return(dt(data.frame(
      message = "No samples flagged: every assignment is clear of the cutoff.")))
    x$support <- signif(x$support, 4); x$rival <- signif(x$rival, 4)
    x$margin <- signif(x$margin, 3)
    dt(x)
  })

  output$plot_margin <- renderPlot({
    cf <- conf(); req(cf)
    d <- cf[!is.na(cf$margin), ]
    d$status <- factor(d$status, levels = c("Confident", "Check", "Underpowered", "Uncertain"))
    ggplot(d, aes(pmax(margin, 1e-6), status, colour = status)) +
      geom_jitter(height = 0.22, width = 0, size = 2.3, alpha = 0.8) +
      scale_x_log10(labels = function(z) format(z, scientific = FALSE, drop0trailing = TRUE)) +
      scale_colour_manual(values = c(Confident = "#4a7c59", Check = "#c9a227",
                                     Underpowered = MUTED, Uncertain = ACCENT),
                          guide = "none") +
      labs(x = "Margin: how far this sample's own cluster beats its nearest rival",
           y = NULL,
           subtitle = "Log scale. Points near the left are assignments that hinge on the cutoff.") +
      theme_gid()
  })

  output$tbl_map <- renderDT(dt(within(map_tbl(), {
    support <- signif(support, 4); rival <- signif(rival, 4); margin <- signif(margin, 3)
  })))

  output$tbl_roster <- renderDT({
    cf <- conf(); req(cf)
    sp <- split(cf, cf$individual)
    r <- do.call(rbind, lapply(names(sp), function(k) {
      x <- sp[[k]]
      data.frame(individual = k,
                 group = if (!is.null(x$group)) x$group[1] else NA_character_,
                 n_in_individual = nrow(x),
                 weakest_margin = min(x$margin, na.rm = TRUE),
                 status = as.character(x$status[which.min(x$margin)]),
                 members = paste(x$sample, collapse = ", "),
                 stringsAsFactors = FALSE)
    }))
    r$weakest_margin <- signif(r$weakest_margin, 3)
    dt(r[order(-r$n_in_individual, r$weakest_margin), ])
  })

  ## ---- posterior cutoff calibration ---------------------------------------
  calib <- reactive({
    r <- res(); req(r)
    b <- r$methods$probabilistic; req(b)
    do.call(rbind, lapply(names(b$by_group), function(g) {
      e <- b$by_group[[g]]
      if (is.null(e$pairs) || is.null(e$pairs$posterior_same)) return(NULL)
      ids <- e$assignment$sample
      tb <- gid_calibrate_threshold(e$pairs, ids, min_loci = input$min_loci,
                                    linkage = input$linkage)
      if (is.null(tb)) return(NULL)
      cbind(group = g, tb, max_post = attr(tb, "max_posterior"),
            n_ids = length(ids))
    }))
  })

  ## The recommendation, and the reasoning behind it. Users kept asking what to
  ## pick, which meant the table alone was not answering the question.
  recommendation <- reactive({
    tb <- calib(); req(tb)
    g <- tb$group[which.max(tb$n_ids)]
    x <- tb[tb$group == g, ]
    r <- gid_recommend_cutoff(x, budget = input$err_budget %||% 1,
                              max_post = x$max_post[1])
    if (is.null(r)) return(NULL)
    r$group <- g; r$max_post <- x$max_post[1]; r$table <- x
    r
  })

  output$calib_verdict <- renderUI({
    r <- recommendation(); req(r)
    x <- r$table; cur <- input$post_cut
    at_cur <- x[which.min(abs(x$cutoff - cur)), ]

    head <- if (is.na(r$cutoff)) {
      tags$div(class = "gid-flag",
        tags$b("No cutoff works on this data, and that is the finding. "),
        sprintf("The most similar pair reaches a posterior of %s. ", signif(r$max_post, 4)),
        "Every sample will come back as its own animal whatever you choose, ",
        "because nothing can clear any bar. Check the loci in use, the minimum ",
        "loci per pair, and whether the error rates are far from reality \u2014 ",
        "an overstated dropout rate flattens every posterior.")
    } else {
      tags$div(
        class = if (r$reason == "plateau") "gid-flag gid-ok" else "gid-flag",
        tags$div(style = "font-size:1.02rem;margin-bottom:.35rem",
          tags$b(sprintf("Recommended cutoff: %s", signif(r$cutoff, 6))),
          sprintf(" \u2014 %d individuals", r$n_individuals)),
        r$note,
        tags$div(style = "margin-top:.4rem",
          sprintf("At that cutoff: %.2f expected false merges and %.1f expected missed recapture pairs.",
                  r$exp_false_merges, r$exp_missed_pairs)),
        if (abs(cur - r$cutoff) > 1e-12) tags$div(
          style = "margin-top:.5rem",
          actionButton("apply_cutoff",
                       sprintf("Use %s", signif(r$cutoff, 6)),
                       class = "btn-sm btn-primary"),
          tags$span(class = "gid-hint", style = "margin-left:.6rem",
            sprintf("You are on %s: %.2f false merges, %.1f missed pairs, %d individuals.",
                    cur, at_cur$exp_false_merges, at_cur$exp_missed_pairs,
                    at_cur$n_individuals)))
        else tags$div(class = "gid-hint", style = "margin-top:.4rem",
                      "This is what you are already using."))
    }

    tagList(head,
      if (!is.na(r$cutoff)) tags$p(class = "gid-hint", style = "margin-top:.5rem",
        sprintf("Across every cutoff this panel can reach, the answer moves between %d and %d animals. ",
                r$range_individuals[1], r$range_individuals[2]),
        "A wide plateau means the cutoff is barely doing any work; a narrow one ",
        "means your answer turns on it, and the honest thing to report is the range. ",
        tags$a(href = "#m-cutoff", class = "gid-jump",
               onclick = "return false;", "How this is decided")))
  })

  observeEvent(input$apply_cutoff, {
    r <- recommendation(); req(r, !is.na(r$cutoff))
    updateSliderInput(session, "post_cut", value = r$cutoff)
    showNotification(sprintf("Cutoff set to %s. Press Identify individuals to re-run.",
                             signif(r$cutoff, 6)), type = "message", duration = 7)
  })

  output$plot_calib <- renderPlot({
    tb <- calib(); req(tb)
    g <- tb$group[which.max(tb$n_ids)]
    x <- tb[tb$group == g, ]
    d <- rbind(data.frame(cutoff = x$cutoff, v = x$exp_false_merges,
                          what = "Expected false merges (two animals joined)"),
               data.frame(cutoff = x$cutoff, v = x$exp_missed_pairs,
                          what = "Expected missed pairs (one animal split)"))
    ggplot(d, aes(factor(cutoff), v + 0.01, colour = what, group = what)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
      geom_vline(xintercept = which.min(abs(x$cutoff - input$post_cut)),
                 linetype = 2, colour = ACCENT) +
      scale_y_log10() +
      scale_colour_manual(values = c(INK, ACCENT)) +
      labs(x = "Posterior cutoff", y = "Expected number of errors (log scale)",
           colour = NULL,
           subtitle = paste("Dashed line = your current setting. Both curves are in the same units,",
                            "so\nwhere they cross is where the two kinds of mistake are equally likely.")) +
      theme_gid()
  })

  output$tbl_calib <- renderDT({
    tb <- calib(); req(tb)
    x <- tb[, c("group", "cutoff", "n_individuals", "n_accepted",
                "exp_false_merges", "exp_missed_pairs")]
    x$exp_false_merges <- round(x$exp_false_merges, 3)
    x$exp_missed_pairs <- round(x$exp_missed_pairs, 1)
    dt(x)
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
  ## Re-run on consensus genotypes, where each run takes a fraction of a
  ## second, at the residual rates the consensus route uses.
  sens <- eventReactive(run_count(), {
    req(run_count() > 0)
    p <- prep(); req(p)
    r <- res(); req(r)
    cr <- r$rates$consensus
    d0 <- signif(cr[["dropout"]], 3)
    grid <- expand.grid(
      post_cut = sort(unique(c(0.95, 0.99, 0.999, 0.9999, input$post_cut))),
      dropout  = sort(unique(c(0.001, 0.01, 0.05, d0))),
      kinship  = c("unrelated", "full_sib"),
      stringsAsFactors = FALSE)
    withProgress(message = "Re-running across settings", value = 0, {
      out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
        incProgress(1 / nrow(grid))
        x <- gid_by_group(p$gt, p$grp, gid_method_lr,
                          dropout = grid$dropout[k],
                          false_allele = cr[["false_allele"]],
                          kinship = grid$kinship[k], post_cut = grid$post_cut[k],
                          min_loci = input$min_loci, linkage = input$linkage,
                          prior_same = if (is.na(input$prior)) NULL else input$prior)
        data.frame(grid[k, ], n_individuals = length(unique(x$assignment$individual)),
                   n_conflicting = x$n_conflict, stringsAsFactors = FALSE)
      }))
    })
    out$current <- out$post_cut == input$post_cut & out$dropout == d0 &
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
    r <- res(); req(r)
    rr <- r$rates
    two <- length(r$routes) == 2
    num <- function(x) format(signif(x, 3), scientific = FALSE, drop0trailing = TRUE)
    rate_line <- function(lab, v, src)
      sprintf("  %-19s dropout %s, false allele %s  (%s)", lab, num(v[["dropout"]]),
              num(v[["false_allele"]]), src)
    oe <- r$own_error
    route_lines <- if (two) unlist(lapply(names(r$route_cmp), function(k) {
      rc <- r$route_cmp[[k]]
      sprintf("  %-20s every replicate %d, consensus first %d, %d sample%s placed differently",
              GID_METHODS[[k]]$label, rc$n_individuals[["replicates"]],
              rc$n_individuals[["consensus"]], rc$n_moved, if (rc$n_moved == 1) "" else "s")
    }))
    paste(c(
      "genoID settings",
      "----------------------------------------",
      sprintf("file                 %s", if (nzchar(source_name())) source_name() else "-"),
      sprintf("samples analyzed     %d", nrow(r$gt)),
      sprintf("loci used            %d  (%s)", ncol(r$gt), paste(colnames(r$gt), collapse = " ")),
      sprintf("loci dropped by you  %s",
              if (length(setdiff(detected(), loci()))) paste(setdiff(detected(), loci()), collapse = " ") else "none"),
      sprintf("grouped by           %s",
              if (nzchar(input$group_col %||% "")) input$group_col else "not grouped"),
      "",
      sprintf("Method shown on the Individuals tab: %s", GID_METHODS[[input$method]]$label),
      sprintf("Genotypes shown:                     %s",
              if (is.null(rr$reps)) "one genotype per sample"
              else if (rr$primary == "reps") "every replicate reaction" else "consensus genotypes"),
      "",
      "Error model",
      if (!is.null(rr$reps)) c(
        rate_line("per reaction", rr$reps, rr$source[["reps"]]),
        if (length(oe)) sprintf("  %-19s each sample's own error scale; %d of %d above 1, up to %s",
                                "per sample", sum(oe > 1, na.rm = TRUE), length(oe),
                                format(max(oe, na.rm = TRUE))),
        rate_line("consensus calls", rr$consensus, rr$source[["consensus"]]))
      else rate_line("per genotype", rr$consensus, rr$source[["consensus"]]),
      "",
      "Sethi et al. 2016",
      sprintf("  relationships      %s",
              if (length(input$sethi_rel)) paste(input$sethi_rel, collapse = ", ") else "none chosen, not run"),
      sprintf("  accept when        lambda > %s", input$lambda_cut),
      "",
      "Likelihood ratio",
      sprintf("  alternative        %s", input$kinship),
      sprintf("  posterior cutoff   %s", input$post_cut),
      sprintf("  prior P(recapture) %s",
              if (is.na(input$prior)) "estimated from data" else as.character(input$prior)),
      "",
      "Shared",
      sprintf("  min loci per pair  %d", input$min_loci),
      sprintf("  cluster rule       %s", input$linkage),
      sprintf("  mismatch threshold %d loci", input$max_mismatch),
      sprintf("  GenAlEx Matches    %s",
              if (isTRUE(input$run_genalex))
                sprintf("run, near matches within %d loci", input$near_match) else "not run"),
      "",
      sprintf("Result: %d samples -> %d individuals", nrow(r$gt),
              length(unique(best()$assignment$individual))),
      if (two) c("", "Every replicate or a consensus first", route_lines)),
      collapse = "\n")
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
          "Nothing to download yet. Press Identify individuals first.",
          type = "warning"))
      send_file(sprintf("genoID_%s_%s.csv", name, format(Sys.Date())), as_csv(d))
    })
  }

  gid_project_server(input, output, session, list(
    raw = raw, demo_loaded = demo_loaded, source_name = source_name,
    send_file = send_file, run_count = run_count, res = res,
    prep = prep, restored_view = restored_view, pending_map = pending_map))

  map_api <- gid_map_server(input, output, session, list(
    res = res, best = best, conf = conf, pid = pid, prep = prep, run_status = run_status,
    send_file = send_file, run_count = run_count, source_name = source_name,
    restored_view = restored_view, pending_map = pending_map))

  ## Every table's Download CSV button lands here. The table sends its output
  ## id and dt() stashed the complete frame under that id when it rendered, so
  ## what leaves is the whole table regardless of which page is on screen.
  observeEvent(input$gid_dl_table, {
    id  <- input$gid_dl_table$id
    env <- session$userData$gid_tables
    d <- if (!is.null(env) && !is.null(id) &&
             exists(id, envir = env, inherits = FALSE)) get(id, envir = env)
    if (is.null(d) || !NROW(d))
      return(showNotification("Nothing to download from this table yet.",
                              type = "warning"))
    send_file(sprintf("genoID_%s_%s.csv", sub("^tbl_", "", id), format(Sys.Date())),
              as_csv(d))
  })

  dl_button("dl_assign", "assignments",           function() final_tbl())
  dl_button("dl_pairs",  "pairwise",              function() best()$pairs)
  dl_button("dl_map",    "sample_map",            function() map_tbl())
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

  ## The gid_ call behind a method, at one route's error rates.
  method_call <- function(key, rate = NULL, reps = FALSE) {
    num  <- function(x) format(signif(x, 4), scientific = FALSE, drop0trailing = TRUE)
    rt   <- if (!is.null(rate)) c(sprintf("dropout = %s", num(rate[["dropout"]])),
                                  sprintf("false_allele = %s", num(rate[["false_allele"]])))
    tail <- c(sprintf("min_loci = %s", input$min_loci), sprintf('linkage = "%s"', input$linkage))
    body <- switch(key,
      sethi = c(rt, sprintf('relationships = c("%s")', paste(input$sethi_rel, collapse = '", "')),
                sprintf("lambda_cut = %s", input$lambda_cut), tail),
      probabilistic = c(rt, sprintf('kinship = "%s"', input$kinship),
                        sprintf("post_cut = %s", input$post_cut),
                        if (!is.na(input$prior)) sprintf("prior_same = %s", input$prior), tail),
      threshold = c(sprintf("max_mismatch = %s", input$max_mismatch), tail),
      genalex = c(sprintf("near_match_loci = %s", input$near_match), tail),
      exact = tail)
    fn <- c(sethi = "gid_method_sethi", probabilistic = "gid_method_lr",
            threshold = "gid_method_threshold", genalex = "gid_method_genalex",
            exact = "gid_method_exact")[[key]]
    sprintf("%s,\n                     %s", fn,
            paste(c(body, if (reps) "reps = reps"), collapse = ",\n                     "))
  }

  observeEvent(input$dl_script, {
    req(input$id_col, length(loci()) > 0)
    r <- tryCatch(res(), error = function(e) NULL)
    p <- tryCatch(prep(), error = function(e) NULL)
    key <- input$method
    lik <- key %in% c("sethi", "probabilistic")
    on_reps <- !is.null(p$reps_all)
    rr <- r$rates %||% list(primary = "consensus",
                            consensus = c(dropout = input$dropout, false_allele = input$false_allele))
    q  <- function(x) paste0('"', paste(x, collapse = '", "'), '"')
    id <- input$id_col
    flags <- if (isTRUE(input$drop_flagged)) ", strip_flags = FALSE" else ""
    grp <- if (nzchar(input$group_col %||% ""))
      sprintf('raw[["%s"]][match(rownames(f$gt), raw[["%s"]])]', input$group_col, id)
      else 'rep("all", nrow(f$gt))'
    other <- setdiff(c("reps", "consensus"), rr$primary)
    lines <- c(
      "## Reproduces what the genoID app just did.",
      'source("genoID_core.R")',
      "",
      sprintf('raw  <- gid_read("%s")', input$file$name %||% "genotypes.csv"),
      sprintf('raw  <- raw[!gid_blank_ids(raw, "%s"), ]', id),
      if (on_reps && length(input$rep_drop))
        sprintf('raw  <- raw[!raw[["%s"]] %%in%% c(%s), ]   # rows excluded in the app',
                input$rep_col, q(input$rep_drop)),
      sprintf("loci <- c(%s)", q(loci())),
      if (on_reps) c(
        "",
        "## one row per PCR reaction, then a consensus genotype per sample",
        sprintf('rgt  <- gid_matrix(raw, "%s", loci%s)', id, flags),
        sprintf('rsm  <- as.character(raw[["%s"]])', id),
        sprintf("gt   <- gid_consensus_from_reps(rgt, rsm, hom_n = %d, het_n = %d)",
                p$hom_n, p$het_n))
      else sprintf('gt   <- gid_matrix(raw, "%s", loci%s)', id, flags),
      sprintf("f    <- gid_filter(gt, min_locus_call = %s, min_sample_call = %s)",
              input$min_locus_call, input$min_sample_call),
      if (isTRUE(input$drop_mixtures)) c(
        sprintf("mx   <- gid_mixture_check(f$gt, alpha = 0.05, excess_cut = %s, min_loci = %d)",
                input$mix_cut, max(10L, as.integer(input$min_loci))),
        'f$gt <- f$gt[setdiff(rownames(f$gt), mx$sample[mx$flag == "possible mixture"]), , drop = FALSE]'),
      if (on_reps) "reps <- list(gt = rgt[, colnames(f$gt)], sample = rsm)",
      sprintf("grp  <- %s", grp),
      "")
    if (lik && on_reps) {
      cf <- function(route) method_call(key, rr[[route]], reps = route == "reps")
      lines <- c(lines,
        "## The error rates are the ones the app used: per reaction for the route that",
        "## uses every replicate, and the residual rate of a consensus call for the",
        "## other. gid_estimate_error(f$gt, reps) measures the per-reaction pair afresh,",
        "## and each sample's own error scale is measured inside the method.",
        sprintf("res  <- gid_by_group(f$gt, grp, %s)   # %s, as the app showed",
                cf(rr$primary), route_label(rr$primary)),
        if (!is.null(rr[[other]])) c(
          "",
          sprintf("## the other route: %s", route_label(other)),
          sprintf("res_other <- gid_by_group(f$gt, grp, %s)", cf(other)),
          sprintf("print(gid_compare_routes(%s)$moved)   # samples the routes place differently",
                  if (rr$primary == "reps") "res, res_other" else "res_other, res")),
        "")
    } else
      lines <- c(lines,
        sprintf("res  <- gid_by_group(f$gt, grp, %s)",
                method_call(key, if (lik) rr$consensus)), "")
    if (identical(key, "probabilistic") && isTRUE(input$use_spatial))
      lines <- c(lines, "## The app also added spatial evidence (gid_apply_spatial), not repeated here.")
    lines <- c(lines, 'write.csv(res$assignment, "individual_assignments.csv", row.names = FALSE)')
    send_file("genoID_reproduce.R", paste(c(lines, ""), collapse = "\n"),
              type = "text/plain;charset=utf-8")
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
      sample_map_confidence  = map_tbl(),
      pairwise_comparisons   = best()$pairs,
      method_comparison      = r$cmp$summary,
      assignments_all_methods = r$cmp$table,
      probability_of_identity = pid(),
      locus_qc               = gid_locus_stats(r$gt),
      sample_power           = power(),
      individual_genotypes   = data.frame(individual = rownames(ic$genotypes),
                                          ic$genotypes))
    if (!is.null(best()$conflicts)) parts$cluster_conflicts <- best()$conflicts
    ## both routes, side by side, when the file had replicates
    if (length(r$routes) == 2) {
      ar <- data.frame(sample = rownames(r$gt), stringsAsFactors = FALSE)
      for (k in names(r$route_cmp)) for (rt in c("reps", "consensus")) {
        a <- r$routes[[rt]][[k]]$assignment
        ar[[sprintf("%s_%s", k, if (rt == "reps") "every_reaction" else "consensus_first")]] <-
          a$individual[match(ar$sample, a$sample)]
      }
      if (!is.null(r$own_error)) ar$own_error <- unname(r$own_error[ar$sample])
      parts$assignments_by_route <- ar
      mv <- route_moved()
      if (nrow(mv)) parts$samples_placed_differently_by_route <- mv
    }
    rq <- tryCatch(rep_quality(), error = function(e) NULL)
    if (!is.null(rq)) parts$replicate_consistency <- rq
    ## the mapped samples, if coordinates were supplied
    geo <- map_api$table()
    if (!is.null(geo) && nrow(geo)) parts$scat_map <- geo

    zipped <- try({
      d <- file.path(tempdir(), paste0("genoID_", as.integer(Sys.time())))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      for (nm in names(parts))
        utils::write.csv(parts[[nm]], file.path(d, paste0(nm, ".csv")), row.names = FALSE)
      ## and the interactive map, so the zip is the whole shareable package
      mh <- map_api$map_html()
      if (!is.null(mh)) writeLines(mh, file.path(d, "interactive_map.html"))
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
