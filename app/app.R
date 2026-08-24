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
mt <- c("map_tab.R", "app/map_tab.R", "../app/map_tab.R")
## local = TRUE matters: Shiny evaluates app.R in its own environment, so a
## default source() would define these in globalenv -- the parent -- where they
## could not see GID_METHODS or dt().
source(mt[file.exists(mt)][1], local = TRUE)

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
    ## hidden, so the comparison tab can still run all five.
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
      "input.method == 'exact'",
      hint("No tolerance to set \u2014 that is the point of this baseline.")),

    ## error model: used by both likelihood methods
    conditionalPanel(
      "input.method == 'probabilistic' || input.method == 'sethi'",
      sliderInput("dropout", "Allelic dropout rate", 0, 0.20, 0.005, step = 0.001),
      sliderInput("false_allele", "False allele rate", 0, 0.10, 0.002, step = 0.001),
      actionButton("est_error", "Estimate these from my data",
                   icon = icon("wand-magic-sparkles"),
                   class = "btn-sm btn-outline-primary w-100"),
      uiOutput("est_error_note"),
      hint("Guessing these is the weakest part of the analysis. If you uploaded ",
           "replicates the app can measure them properly; if not it can at least ",
           "get the order of magnitude. See ", tags$b("Measuring your own error rates"),
           " on the Methods tab.")),

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
        checkboxInput("run_genalex", "GenAlEx Matches", TRUE))
    ),

    ## Map controls. Keyed on the open tab rather than the method, because the
    ## map deliberately lets you colour by a model other than the one the
    ## Individuals tab is showing.
    conditionalPanel(
      "input.nav == 'Scat map'",
      tags$hr(),
      tags$p(class = "gid-label", "Map"),
      selectInput("map_model", "Colour individuals by", choices = NULL),
      hint("Switch models to see which samples change hands. The points stay ",
           "put; only the colouring changes."),
      checkboxInput("map_grey", "Grey out animals seen once", TRUE),
      uiOutput("map_coord_ui"),
      tags$hr(),
      tags$p(class = "gid-label", "Link samples of the same animal"),
      radioButtons("map_link_style", NULL,
                   c("No links" = "none",
                     "Spider (lines to centre)" = "spider",
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
      vbox("Samples per animal", "vb_persample", "layer-group", "light",
           key = "median_samples"),
      vbox("Largest cluster", "vb_max", "maximize", "light", key = "max_cluster")
    ),
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
    card(card_header("Individual assignments"), DTOutput("tbl_final")),
    conditionalPanel("input.method == 'probabilistic'",
      card(card_header("Choosing your posterior cutoff"),
        hint("The cutoff is the one setting with no data behind it \u2014 unless ",
             "you use the fact that the method already gives every pair a ",
             "probability. Adding those up says how many mistakes each cutoff ",
             "costs, in both directions, so the choice becomes a stated trade ",
             "instead of a convention."),
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
         hint("The same information the other way round \u2014 each animal and ",
              "the samples that make it up."),
         DTOutput("tbl_roster"))
  ),

  gid_map_tab_ui(),

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

  ## Replicates the file contains, independent of whether the user chose to
  ## analyse with them. Error rates can always be measured from replicates that
  ## exist, even if the analysis is running on consensus calls.
  reps_available <- reactive({
    df <- raw()
    if (is.null(df) || !isTRUE(has_reps()) || !nzchar(input$rep_col %||% "")) return(NULL)
    if (!length(loci())) return(NULL)
    lab  <- as.character(df[[input$rep_col]])
    keep <- !lab %in% (input$rep_drop %||% character(0))
    if (!any(keep)) return(NULL)
    list(gt = gid_matrix(df[keep, , drop = FALSE], input$id_col, loci(),
                         sep = gid_guess_sep(df, exclude = input$id_col)),
         sample = as.character(df[[input$id_col]])[keep])
  })

  ## ---- estimating the error rates from the data itself --------------------
  err_est <- reactiveVal(NULL)
  observeEvent(input$est_error, {
    p <- try(prep(), silent = TRUE)
    if (inherits(p, "try-error") || is.null(p))
      return(showNotification("Load a file and choose your loci first.", type = "warning"))
    rp <- reps_available()
    e <- try(withProgress(message = "Measuring genotyping error", value = 0.4,
                          gid_estimate_error(p$gt, rp)), silent = TRUE)
    if (inherits(e, "try-error"))
      return(showNotification(paste("Could not estimate:", conditionMessage(attr(e, "condition"))),
                              type = "error"))
    err_est(e)
    ## Replicate rates apply to a single reaction. If the analysis is running on
    ## consensus calls, the model wants the much smaller rate that survives the
    ## consensus rule, so convert before filling the sliders in.
    d <- e$dropout; f <- e$false_allele
    converted <- FALSE
    if (identical(e$method, "replicates") && is.null(p$reps)) {
      n_per <- mean(table(rp$sample))
      pr <- gid_propagate_error(d, f, n_rep = max(2L, round(n_per)),
                                rule = "taberlet", hom_n = 3, het_n = 2,
                                per_rep_missing = mean(is.na(rp$gt)),
                                nsim = 20000)
      d <- max(pr[["dropout"]], 1e-4); f <- max(pr[["false_allele"]], 1e-4)
      converted <- TRUE
    }
    err_est(c(e, list(applied_dropout = d, applied_false = f, converted = converted)))
    if (is.na(f)) f <- max(0.001, round(d / 3, 3))
    updateSliderInput(session, "dropout", value = round(d, 3))
    updateSliderInput(session, "false_allele", value = round(f, 3))
  })

  output$est_error_note <- renderUI({
    e <- err_est(); if (is.null(e)) return(NULL)
    if (identical(e$method, "replicates"))
      tags$div(class = "gid-flag gid-ok", style = "font-size:.78rem;margin:.5rem 0",
        tags$b("Measured from your replicates."), " Maximum likelihood over ",
        e$n_samples, " samples, integrating over the unknown true genotype so no ",
        "consensus is involved.",
        tags$div(style = "margin-top:.3rem;font-family:var(--bs-font-monospace)",
          sprintf("dropout %.3f [%.3f-%.3f]   false allele %.3f [%.3f-%.3f]",
                  e$dropout, e$dropout_ci[1], e$dropout_ci[2],
                  e$false_allele, e$false_ci[1], e$false_ci[2])),
        tags$div(style = "margin-top:.3rem",
          if (isTRUE(e$converted)) tagList(
            "You are analysing consensus calls, so the sliders were set to the ",
            "much smaller residual rate that survives the multi-tube rule: ",
            tags$code(sprintf("dropout %.4f, false allele %.4f",
                              e$applied_dropout, e$applied_false)), ".")
          else "These are per-reaction rates, which is what the replicate route wants."))
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
          "The false-allele rate cannot be got this way at all and has been set ",
          "to a third of the dropout estimate."),
        tags$div(style = "margin-top:.3rem",
          tags$b("To measure it properly, genotype some samples twice"), " and ",
          "upload the replicates. Two reactions on a subset is enough."))
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
      incProgress(0.1, detail = "comparing")
      ord <- c("exact", "threshold", "genalex", "sethi", "probabilistic")
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
      sprintf("Overall recapture rate %.0f%% \u2014 that share of your samples were repeats of an animal already seen. ",
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

  ## ---- sample -> individual map with per-sample confidence ----------------
  conf <- reactive({
    r <- res(); req(r); b <- best(); req(b)
    ## A method that clusters internally exposes no per-pair score, so borrow
    ## the locus-mismatch counts, which are a property of the data.
    if (is.null(b$pairs) ||
        (!is.null(b$by_group) &&
         all(vapply(b$by_group, function(e) is.null(e$pairs), TRUE)))) {
      src <- r$methods$exact
      b <- if (!is.null(b$by_group)) {
        b$by_group <- Map(function(e, s) { e$pairs <- s$pairs; e },
                          b$by_group, src$by_group[names(b$by_group)]); b
      } else { b$pairs <- src$pairs; b }
    }
    gid_sample_confidence(b, gt = r$gt, pid_tab = pid(),
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

  output$calib_verdict <- renderUI({
    tb <- calib(); req(tb)
    ## report on the block with the most samples; a handful of samples in a
    ## minor block cannot reach a high posterior and would misdescribe the run
    g  <- tb$group[which.max(tb$n_ids)]
    x  <- tb[tb$group == g, ]
    mx <- x$max_post[1]
    cur <- input$post_cut
    usable <- x$cutoff <= mx
    rng <- if (any(usable)) range(x$n_individuals[usable]) else c(NA, NA)
    within <- x$cutoff[x$exp_false_merges <= 1 & usable]

    if (mx < cur) {
      sug <- if (any(usable)) max(x$cutoff[usable]) else NA
      return(tags$div(class = "gid-flag",
        tags$b("This is why every sample is coming back as its own animal. "),
        sprintf("The most similar pair in your data reaches a posterior of %s, and your cutoff is %s. ",
                signif(mx, 6), cur),
        "Nothing can clear the bar, so nothing matches. That is a statement about ",
        "the bar, not about your animals.",
        if (!is.na(sug)) tags$div(style = "margin-top:.4rem",
          tags$b(sprintf("Try %s instead.", sug)),
          sprintf(" At that cutoff you would expect %.2f false merges across the whole dataset.",
                  x$exp_false_merges[x$cutoff == sug]))
        else tags$div(style = "margin-top:.4rem",
          "No cutoff on the scale works here, which means the panel itself is ",
          "not separating your samples. Check the loci in use, the minimum loci ",
          "per pair, and whether the error rates are wildly off.")))
    }

    tags$div(class = "gid-flag gid-ok",
      tags$b(sprintf("Across every cutoff your panel can actually reach, the answer moves between %d and %d animals.",
                     rng[1], rng[2])),
      if (length(within) == sum(usable))
        sprintf(" Expected false merges stay under 1 throughout, so within this range the cutoff is barely doing any work \u2014 report the range and move on.")
      else sprintf(" Expected false merges stay under 1 for cutoffs of %s and above.",
                   min(within)),
      tags$div(style = "margin-top:.4rem",
        sprintf("At your current %s: expect %.2f false merges and %.1f missed recapture pairs. ",
                cur, x$exp_false_merges[which.min(abs(x$cutoff - cur))],
                x$exp_missed_pairs[which.min(abs(x$cutoff - cur))]),
        "Raising it further trades a fraction of a merge for many missed matches."))
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

  gid_map_server(input, output, session, list(
    res = res, best = best, conf = conf, prep = prep, run_status = run_status,
    send_file = send_file))

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
