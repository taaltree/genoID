## project.R -- save a whole session and pick it up again --------------------
##
## A project is one JSON file holding the data you uploaded plus every setting
## you chose. It is written and read entirely in the browser, so it works the
## same in a local Shiny session and in the WebAssembly build, and the genotypes
## never leave the machine.
##
## The data is carried inside the file rather than referenced by path: a path
## would break the moment the file moved, was renamed, or was opened on a
## colleague's laptop, which is exactly when you want the project to work.
## ---------------------------------------------------------------------------

GID_PROJECT_VERSION <- 1L

## Every input worth restoring, grouped by the updater it needs. Action buttons
## and the file input are deliberately absent: they are events, not state.
GID_SAVED_INPUTS <- list(
  select = c("id_col", "kinship", "method", "rep_col", "sheet",
             "map_date_col", "map_lat", "map_lon", "map_model", "map_sex_col",
             "map_utm_e", "map_utm_n"),
  selectize = c("group_col", "rep_drop", "use_loci",
                "map_link_who", "map_show_who", "map_show_year"),
  numeric = c("lambda_cut", "map_fig_width", "map_utm_zone", "max_mismatch",
              "min_loci", "min_locus_call", "min_sample_call", "near_match",
              "prior"),
  slider = c("dropout", "false_allele", "post_cut"),
  checkbox = c("map_fig_labels", "map_grey", "map_utm_south",
               "run_genalex", "run_sethi"),
  checkboxgroup = c("map_show_sex", "sethi_rel"),
  radio = c("linkage", "map_link_style", "rep_mode"))

#' Everything needed to reopen this session, as a JSON string.
gid_project_bundle <- function(input, data, source_name = "", demo = FALSE,
                               had_run = FALSE) {
  vals <- list()
  for (type in names(GID_SAVED_INPUTS))
    for (id in GID_SAVED_INPUTS[[type]]) {
      v <- input[[id]]
      ## NULL means the widget never existed this session; keeping it out lets
      ## the restore leave that control at its own default.
      if (!is.null(v)) vals[[id]] <- v
    }

  csv <- paste(utils::capture.output(utils::write.csv(data, row.names = FALSE)),
               collapse = "\n")
  ## as.character() matters: toJSON() returns a "json"-classed value, and Shiny
  ## would re-serialise that as a nested object rather than as the text of the
  ## file, so the browser would be handed [object Object] to save.
  as.character(jsonlite::toJSON(list(
    genoID_project = GID_PROJECT_VERSION,
    saved          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    source         = source_name,
    demo           = isTRUE(demo),
    had_run        = isTRUE(had_run),
    n_rows         = nrow(data),
    n_cols         = ncol(data),
    settings       = vals,
    ## where the map was looking, so it reopens on the same ground rather than
    ## snapping back to the full extent
    map_view       = if (!is.null(input$geo_map_center))
                       list(lng = input$geo_map_center$lng,
                            lat = input$geo_map_center$lat,
                            zoom = input$geo_map_zoom) else NULL,
    map_selected   = input$geo_map_marker_click$id,
    data_gz_b64    = jsonlite::base64_enc(
                       memCompress(charToRaw(csv), "gzip"))),
    auto_unbox = TRUE, null = "null"))
}

#' Read a project file back. Errors carry what is actually wrong with the file,
#' since "could not open project" tells nobody anything.
gid_project_parse <- function(txt) {
  p <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE),
                error = function(e)
                  stop("This is not a readable JSON file.", call. = FALSE))
  if (is.null(p$genoID_project))
    stop("This JSON file is not a genoID project.", call. = FALSE)
  if (as.integer(p$genoID_project) > GID_PROJECT_VERSION)
    stop(sprintf(paste("This project was saved by a newer version of genoID",
                       "(format %s, this app reads %s). Update the app."),
                 p$genoID_project, GID_PROJECT_VERSION), call. = FALSE)
  if (is.null(p$data_gz_b64))
    stop("This project file carries no data.", call. = FALSE)

  csv <- tryCatch(
    rawToChar(memDecompress(jsonlite::base64_dec(p$data_gz_b64), "gzip")),
    error = function(e)
      stop("The data inside this project is corrupt.", call. = FALSE))
  ## every column stays character: genotypes such as 0112 or 1/2 must not be
  ## guessed into numbers on the way back in
  p$data <- utils::read.csv(text = csv, check.names = FALSE,
                            colClasses = "character", stringsAsFactors = FALSE)
  p
}

# ------------------------------------------------------------------------ UI
gid_project_ui <- function() {
  tagList(
    tags$div(
      class = "d-flex gap-2 align-items-center",
      actionButton("proj_save", "Save project", icon = icon("floppy-disk"),
                   class = "btn-sm btn-outline-secondary flex-fill"),
      tags$label(
        class = "btn btn-sm btn-outline-secondary flex-fill mb-0",
        style = "cursor:pointer;text-align:center",
        icon("folder-open"), " Open project",
        tags$input(type = "file", id = "proj_open_native", accept = ".json",
                   style = "display:none"))),
    hint("A project holds your data and every setting in one file, so you can ",
         "close the tab and pick up exactly where you left off."),
    ## The file is read in the browser and handed over as text, which keeps the
    ## whole round trip working under WebAssembly where there is no upload route.
    tags$script(HTML("
      document.addEventListener('change', function (e) {
        if (e.target && e.target.id === 'proj_open_native' && e.target.files.length) {
          var fr = new FileReader();
          fr.onload = function () {
            Shiny.setInputValue('proj_open_text',
              {name: e.target.files[0].name, body: fr.result, nonce: Math.random()},
              {priority: 'event'});
            e.target.value = '';
          };
          fr.readAsText(e.target.files[0]);
        }
      });")))
}

# -------------------------------------------------------------------- SERVER
#' @param deps raw (reactiveVal), demo_loaded (reactiveVal), send_file, and
#'   source_name (a reactive giving the current file's name)
gid_project_server <- function(input, output, session, deps) {

  observeEvent(input$proj_save, {
    d <- deps$raw()
    if (is.null(d))
      return(showNotification("Load a genotype file before saving a project.",
                              type = "warning"))
    txt <- gid_project_bundle(input, d, deps$source_name(),
                              isTRUE(deps$demo_loaded()), (input$run %||% 0) > 0)
    deps$send_file(sprintf("genoID_project_%s.json", format(Sys.Date())),
                   txt, type = "application/json")
    showNotification("Project saved.", type = "message", duration = 4)
  })

  ## Restoring happens in two phases because the controls do not all exist at
  ## the same time. Phase one puts back everything the data alone supports, then
  ## starts the analysis. Phase two waits for results and puts back the map
  ## controls, whose menus -- the model picker, the individual and year filters
  ## -- are built from the results and so were empty a moment earlier.
  ##
  ## Within each phase the settings are applied several times over about a
  ## second, because loading data fires the observers that rebuild the column
  ## pickers from scratch and a single pass gets overwritten a beat later.
  pending <- reactiveVal(NULL)
  phase   <- reactiveVal(0L)
  passes  <- reactiveVal(0L)
  waited  <- reactiveVal(0L)   # ticks spent waiting for results, kept apart
                               # from passes so waiting cannot eat the passes

  apply_settings <- function(cfg) {
    for (type in names(GID_SAVED_INPUTS))
      for (id in GID_SAVED_INPUTS[[type]]) {
        if (is.null(cfg[[id]]) || id %in% GID_MAP_DEFERRED) next
        v <- cfg[[id]]
        switch(type,
          select        = updateSelectInput(session, id, selected = v),
          selectize     = updateSelectizeInput(session, id, selected = v),
          numeric       = updateNumericInput(session, id, value = v),
          slider        = updateSliderInput(session, id, value = v),
          checkbox      = updateCheckboxInput(session, id, value = as.logical(v)),
          checkboxgroup = updateCheckboxGroupInput(session, id, selected = v),
          radio         = updateRadioButtons(session, id, selected = v))
      }
  }

  ## phase 1 -- wait for the app to finish ingesting the data.
  ##
  ## This wait is the whole trick. Setting raw() kicks off the observers that
  ## detect the ID column and the loci and fill those pickers in. Pushing the
  ## saved choices in before those menus exist does not just get overwritten --
  ## selectize silently drops a selection it has no option for, so the server
  ## ends up holding "" while the widget later displays something sensible, and
  ## every downstream req() fails on a value the user can plainly see.
  ##
  ## Event-driven rather than polled: a poll on a short timer competes with the
  ## very flush cycle that carries the new input value back from the browser,
  ## and the value can take seconds to land. Waiting on the event itself has no
  ## such problem and needs no timeout.
  observeEvent(input$id_col, {
    if (is.null(pending()) || phase() != 1L) return()
    if (!nzchar(input$id_col %||% "")) return()
    phase(2L); passes(0L)
  }, ignoreInit = FALSE, ignoreNULL = FALSE)

  ## phase 2 -- now the menus exist, put the saved choices back. Repeated a few
  ## times because the detection observers settle over several flushes.
  observe({
    cfg <- pending()
    if (is.null(cfg) || phase() != 2L) return()
    n <- passes()
    if (n >= 4L) { phase(3L); passes(0L); return() }
    apply_settings(cfg$settings)
    passes(n + 1L)
    invalidateLater(260, session)
  })

  ## phase 3 -- start the analysis, but only once the data actually resolves.
  ## Launching early is not merely slow: res() is an eventReactive, so a run
  ## that fails its own req() caches the failure and never recomputes.
  observe({
    cfg <- pending()
    if (is.null(cfg) || phase() != 3L) return()
    if (!isTRUE(cfg$had_run)) { phase(4L); passes(0L); waited(0L); return() }
    ok <- tryCatch({ p <- deps$prep(); !is.null(p) && nrow(p$gt) >= 2 },
                   error = function(e) FALSE)
    n <- passes()
    if (!ok) {
      if (n >= 20L) {                       # ~6 s
        pending(NULL); phase(0L)
        showNotification(paste("Reopened the project, but its settings do not",
                               "produce a runnable dataset. Adjust the sidebar",
                               "and press Identify individuals."),
                         type = "warning", duration = 12)
        return()
      }
      passes(n + 1L); invalidateLater(300, session); return()
    }
    ## isolate: this observer must not depend on the counter it is setting
    deps$run_count(isolate(deps$run_count()) + 1L)
    phase(4L); passes(0L); waited(0L)
  })

  ## phase 4 -- results exist, so the map controls finally have menus to choose
  ## from: the model picker and the individual and year filters are all built
  ## from the results and were empty until now.
  observe({
    cfg <- pending()
    if (is.null(cfg) || phase() != 4L) return()
    if (isTRUE(cfg$had_run)) {
      r <- tryCatch(deps$res(), error = function(e) NULL)
      if (is.null(r)) {
        w <- waited()
        if (w >= 30L) {                     # ~10 s
          pending(NULL); phase(0L)
          showNotification(paste("Reopened the project and restored its",
                                 "settings, but the analysis did not finish.",
                                 "Press Identify individuals to retry."),
                           type = "warning", duration = 12)
          return()
        }
        waited(w + 1L); invalidateLater(330, session); return()
      }
    }
    n <- passes()
    if (n >= 4L) {
      if (!is.null(cfg$map_view)) deps$restored_view(cfg$map_view)
      pending(NULL); phase(0L)
      showNotification("Project fully restored.", type = "message", duration = 5)
      return()
    }
    apply_settings(cfg$settings)
    passes(n + 1L)
    invalidateLater(260, session)
  })

  ## These three are set by the map module instead, when it builds their menus.
  GID_MAP_DEFERRED <- c("map_show_year", "map_show_who", "map_link_who")

  observeEvent(input$proj_open_text, {
    p <- tryCatch(gid_project_parse(input$proj_open_text$body),
                  error = function(e) e)
    if (inherits(p, "error"))
      return(showNotification(conditionMessage(p), type = "error", duration = 12))

    deps$raw(p$data)
    deps$demo_loaded(isTRUE(p$demo))
    deps$source_name(p$source %||% input$proj_open_text$name)
    deps$pending_map(list(year     = p$settings$map_show_year,
                          who      = p$settings$map_show_who,
                          link_who = p$settings$map_link_who))
    passes(0L); phase(1L); pending(p)

    showNotification(
      sprintf("Opened %s - %d samples, %d columns, saved %s.",
              input$proj_open_text$name, p$n_rows, p$n_cols, p$saved),
      type = "message", duration = 8)
    if (isTRUE(p$had_run))
      showNotification("Recomputing the results this project was saved with...",
                       type = "default", duration = 6)
  })
}
