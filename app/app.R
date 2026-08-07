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

vbox <- function(title, id, ico, theme = "light")
  value_box(title, textOutput(id, inline = TRUE), showcase = icon(ico),
            theme = theme, showcase_layout = showcase_left_center(width = "58px"),
            height = "104px")

dt <- function(x, ...) {
  num <- which(vapply(x, is.numeric, TRUE)) - 1L
  datatable(
    x, rownames = FALSE, extensions = "Buttons",
    options = c(list(pageLength = 12, scrollX = TRUE, dom = "Bfrtip",
                     buttons = list(list(extend = "csv", text = "Download CSV"))),
                if (length(num)) list(columnDefs = list(
                  list(className = "dt-body-right", targets = num)))),
    ...)
}

# ============================================================================ UI
ui <- page_navbar(
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
  header = tags$head(katex_head(), tags$style(HTML(sprintf("
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
  ", MUTED, MUTED, SOFT, INK, INK, INK, INK, ACCENT, MUTED, INK)))),

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

    selectizeInput("use_loci", "Loci to use", choices = NULL, multiple = TRUE,
                   options = list(plugins = list("remove_button"), placeholder = "detected automatically")),
    uiOutput("locus_note"),
    hint("Detected automatically, then yours to correct. Remove species-diagnostic ",
         "and sex markers: they are near-fixed within a species, so they add no ",
         "power and distort the frequency model."),
    actionLink("reset_loci", "Reset to detected", class = "gid-hint"),

    accordion(
      open = "probabilistic",

      accordion_panel(
        "Probabilistic settings", value = "probabilistic", icon = icon("wave-square"),
        selectInput("kinship", "Compare against", selected = "full_sib",
                    choices = c("Full siblings (conservative)" = "full_sib",
                                "Half siblings" = "half_sib",
                                "Unrelated individuals" = "unrelated")),
        hint("The alternative hypothesis. For pack- or group-living species the ",
             "samples competing to be \"a different animal\" are usually relatives, ",
             "so full siblings is the honest null."),
        sliderInput("post_cut", "Posterior probability to accept a match",
                    0.9, 0.9999, 0.999, step = 0.0001),
        sliderInput("dropout", "Allelic dropout rate", 0, 0.20, 0.005, step = 0.001),
        sliderInput("false_allele", "False allele rate", 0, 0.10, 0.002, step = 0.001),
        hint("These are the rates of the genotypes you are uploading. If those are ",
             "already multi-replicate consensus calls, both are far lower than the ",
             "per-PCR rates -- often under 0.005."),
        numericInput("prior", "Prior P(a random pair is a recapture)", NA,
                     min = 0, max = 0.5, step = 0.005),
        hint("Leave blank to estimate it from the data.")
      ),

      accordion_panel(
        "Sethi et al. (2016)", value = "sethi", icon = icon("scale-balanced"),
        checkboxInput("run_sethi", "Run error-tolerant match calling", TRUE),
        checkboxGroupInput("sethi_rel", "Competing relationships",
                           choices = c("Unrelated" = "unrelated",
                                       "Full siblings" = "full_sib",
                                       "Parent-offspring" = "parent_offspring",
                                       "Half siblings" = "half_sib"),
                           selected = c("unrelated", "full_sib", "parent_offspring")),
        hint("This method divides by whichever of these best explains the pair, so you ",
             "do not have to pick one. A match must beat every relationship on the list."),
        numericInput("lambda_cut", "Accept a match when \u039b exceeds", 1, 0.01, 1e6, 1),
        hint("The paper uses \u039b > 1: the match hypothesis simply has to be more likely ",
             "than the best alternative. That is assumption-free but takes no account of how ",
             "many pairs you are testing, so raise it on large datasets.")
      ),

      accordion_panel(
        "Shared settings", value = "shared", icon = icon("sliders"),
        numericInput("min_loci", "Minimum loci compared per pair", 15, 1, 500, 1),
        numericInput("min_sample_call", "Minimum sample call rate", 0.5, 0, 1, 0.05),
        numericInput("min_locus_call", "Minimum locus call rate", 0.25, 0, 1, 0.05),
        radioButtons("linkage", "Cluster rule",
                     c("Single linkage (transitive)" = "single",
                       "Complete linkage (all pairs must match)" = "complete"),
                     selected = "single"),
        hint("Single linkage lets A-B and B-C merge A, B and C. Complete linkage ",
             "requires every pair inside a cluster to match. Compare both: if the ",
             "answer moves, some clusters are being held together by one edge.")
      ),

      accordion_panel(
        "Gut-check methods", value = "gut", icon = icon("check-double"),
        numericInput("max_mismatch", "Mismatch threshold: loci allowed to differ", 1, 0, 20, 1),
        checkboxInput("run_genalex", "Run GenAlEx Matches (Peakall & Smouse)", TRUE),
        numericInput("near_match", "GenAlEx: flag near matches within N loci", 2, 1, 10, 1),
        checkboxInput("run_am", "Run allelematch (Galpern et al. 2012)", TRUE),
        hint("allelematch picks its own threshold. It was designed for ",
             "microsatellites and tends to over-split dense SNP panels with ",
             "missing data.")
      )
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
      vbox("Samples", "vb_samples", "vial", "primary"),
      vbox("Loci used", "vb_loci", "dna", "secondary"),
      vbox("Median call rate", "vb_call", "percent", "light"),
      vbox("Groups", "vb_groups", "layer-group", "light")
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
    uiOutput("run_status_ind"),
    uiOutput("power_warning"),
    layout_columns(
      fill = FALSE, col_widths = c(3, 3, 3, 3),
      vbox("Individuals", "vb_ind", "paw", "primary"),
      vbox("Sampled once", "vb_single", "circle-dot", "light"),
      vbox("Recapture rate", "vb_recap", "repeat", "light"),
      vbox("Largest cluster", "vb_max", "maximize", "light")
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Evidence for every pair of samples"), plotOutput("plot_lr", height = 340),
           hint("Each point is a pair. Pairs to the right are supported as the same ",
                "animal; the vertical line is your acceptance threshold.")),
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
        downloadButton("dl_assign", "Individual assignments", class = "btn-primary w-100"),
        downloadButton("dl_cons",   "Consensus genotype per individual", class = "btn-outline-primary w-100"),
        downloadButton("dl_pairs",  "All pairwise comparisons", class = "btn-outline-primary w-100")
      ),
      tags$hr(),
      layout_columns(
        col_widths = c(4, 4, 4),
        downloadButton("dl_all",    "Everything (zip)", class = "btn-primary w-100"),
        downloadButton("dl_params", "Settings used (for your methods section)",
                       class = "btn-outline-primary w-100"),
        downloadButton("dl_script", "Equivalent R script", class = "btn-outline-primary w-100")
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
    if (anyDuplicated(ids))
      ids <- ave(ids, ids, FUN = function(z)
        if (length(z) == 1) z else paste0(z, "#", seq_along(z)))
    df$.gid_key <- ids
    gt <- gid_matrix(df, ".gid_key", loci(), sep = gid_guess_sep(df, exclude = input$id_col))
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
    list(gt = f$gt, grp = grp[rownames(f$gt)], raw_gt = gt,
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
                         linkage = input$linkage, prior_same = pr)
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
                                  min_loci = input$min_loci, linkage = input$linkage)
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

  best <- reactive({ r <- res(); req(r); r$methods$probabilistic })

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

  output$plot_lr <- renderPlot({
    r <- res(); pp <- r$methods$probabilistic$pairs
    pp <- pp[pp$n_compared >= input$min_loci, ]
    cut_lr <- min(pp$log10_LR[pp$posterior_same >= input$post_cut], Inf)
    ggplot(pp, aes(log10_LR)) +
      geom_histogram(bins = 60, fill = INK, colour = NA) +
      { if (is.finite(cut_lr)) geom_vline(xintercept = cut_lr, colour = ACCENT, linetype = 2) } +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
      labs(x = expression(log[10] ~ "likelihood ratio"), y = "Sample pairs",
           subtitle = sprintf("Alternative hypothesis: %s. Counts are log-scaled.",
                              gsub("_", " ", input$kinship))) +
      theme_gid()
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
    a$n_samples <- as.integer(table(a$individual)[a$individual])
    a[order(-a$n_samples, a$individual, -a$call_rate), ]
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
      "\nProbabilistic method (reported result)\n",
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

  dlcsv <- function(name, data) downloadHandler(
    filename = function() sprintf("genoID_%s_%s.csv", name, format(Sys.Date())),
    content = function(f) write.csv(data(), f, row.names = FALSE))

  output$dl_assign <- dlcsv("assignments", final_tbl)
  output$dl_pairs  <- dlcsv("pairwise", reactive(best()$pairs))
  output$dl_cons   <- dlcsv("individual_genotypes", reactive({
    ic <- gid_individual_consensus(res()$gt, best()$assignment)
    data.frame(individual = rownames(ic$genotypes), ic$genotypes) }))
  output$dl_params <- downloadHandler(
    filename = function() sprintf("genoID_settings_%s.txt", format(Sys.Date())),
    content = function(f) writeLines(params(), f))

  output$dl_script <- downloadHandler(
    filename = "genoID_reproduce.R",
    content = function(f) writeLines(sprintf(paste0(
      '## Reproduces exactly what the genoID app just did.\n',
      'source("genoID_core.R")\n\n',
      'raw  <- gid_read("%s")\n',
      'loci <- c(%s)\n',
      'gt   <- gid_matrix(raw, "%s", loci)\n',
      'f    <- gid_filter(gt, min_locus_call = %s, min_sample_call = %s)\n',
      'grp  <- %s\n\n',
      'res  <- gid_by_group(f$gt, grp, gid_method_lr,\n',
      '                     dropout = %s, false_allele = %s,\n',
      '                     kinship = "%s", post_cut = %s,\n',
      '                     min_loci = %s, linkage = "%s")\n\n',
      'write.csv(res$assignment, "individual_assignments.csv", row.names = FALSE)\n'),
      input$file$name %||% "genotypes.csv",
      paste0('"', paste(loci(), collapse = '", "'), '"'),
      input$id_col, input$min_locus_call, input$min_sample_call,
      if (nzchar(input$group_col %||% "")) sprintf('raw[["%s"]][match(rownames(f$gt), raw[["%s"]])]',
                                                   input$group_col, input$id_col)
      else 'rep("all", nrow(f$gt))',
      input$dropout, input$false_allele, input$kinship, input$post_cut,
      input$min_loci, input$linkage), f))

  output$dl_all <- downloadHandler(
    filename = function() sprintf("genoID_results_%s.zip", format(Sys.Date())),
    content = function(f) {
      d <- tempfile(); dir.create(d)
      r <- res()
      write.csv(final_tbl(), file.path(d, "individual_assignments.csv"), row.names = FALSE)
      write.csv(best()$pairs, file.path(d, "pairwise_comparisons.csv"), row.names = FALSE)
      write.csv(r$cmp$summary, file.path(d, "method_comparison.csv"), row.names = FALSE)
      write.csv(r$cmp$table, file.path(d, "assignments_all_methods.csv"), row.names = FALSE)
      write.csv(pid(), file.path(d, "probability_of_identity.csv"), row.names = FALSE)
      write.csv(gid_locus_stats(r$gt), file.path(d, "locus_qc.csv"), row.names = FALSE)
      write.csv(power(), file.path(d, "sample_power.csv"), row.names = FALSE)
      ic <- gid_individual_consensus(r$gt, best()$assignment)
      write.csv(data.frame(individual = rownames(ic$genotypes), ic$genotypes),
                file.path(d, "individual_genotypes.csv"), row.names = FALSE)
      if (!is.null(best()$conflicts))
        write.csv(best()$conflicts, file.path(d, "cluster_conflicts.csv"), row.names = FALSE)
      writeLines(params(), file.path(d, "settings.txt"))
      if (requireNamespace("zip", quietly = TRUE)) zip::zip(f, list.files(d), root = d)
      else utils::zip(f, file.path(d, list.files(d)), flags = "-j9X")
    })
}

shinyApp(ui, server)
