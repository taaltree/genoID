## map_tab.R -- spatial view of samples, coloured by individual ---------------
##
## Everything here is driven by whatever the chosen method decided, so switching
## the model in the picker recolours the same points rather than recomputing
## anything. Coordinates come from the uploaded file; nothing is geocoded and
## nothing is sent anywhere.
## ---------------------------------------------------------------------------

## Column-name patterns for coordinates. Deliberately strict about bare "x"/"y"
## so a locus called X or a read-count column cannot be mistaken for a position.
GID_LAT_PAT <- "^(lat|latitude|lat_dd|latdd|dec_?lat|decimal_?lat(itude)?|lat_wgs84|ycoord|y_coord)$"
GID_LON_PAT <- "^(lon|lng|long|longitude|lon_dd|londd|long_dd|dec_?long?|decimal_?long?(itude)?|lon_wgs84|xcoord|x_coord)$"
GID_UTME_PAT <- "^(easting|utm_?e(ast)?(ing)?|utmx|utm_x)$"
GID_UTMN_PAT <- "^(northing|utm_?n(orth)?(ing)?|utmy|utm_y)$"

#' Find candidate coordinate columns. Returns names, never guesses silently at
#' a column that does not parse as a plausible coordinate.
gid_find_coord_cols <- function(df) {
  nm <- names(df)
  numeric_enough <- function(k) {
    v <- suppressWarnings(as.numeric(as.character(df[[k]])))
    sum(!is.na(v)) >= max(2, 0.5 * sum(!is.na(df[[k]])))
  }
  pick <- function(pat, lo, hi) {
    hit <- nm[grepl(pat, tolower(nm))]
    hit <- hit[vapply(hit, numeric_enough, TRUE)]
    hit <- hit[vapply(hit, function(k) {
      v <- suppressWarnings(as.numeric(as.character(df[[k]])))
      v <- v[!is.na(v)]
      length(v) > 0 && all(v >= lo & v <= hi)
    }, TRUE)]
    if (length(hit)) hit[1] else ""
  }
  list(lat  = pick(GID_LAT_PAT, -90, 90),
       lon  = pick(GID_LON_PAT, -180, 180),
       utm_e = pick(GID_UTME_PAT, 1e5, 1e6),
       utm_n = pick(GID_UTMN_PAT, 0, 1e7))
}

#' Normalise however a lab writes sex into M / F / U.
#'
#' Deliberately conservative: numeric codings (1/2) are left unknown, because
#' 1 = female and 1 = male are both in common use and guessing wrong would
#' silently mislabel every animal in the dataset.
gid_norm_sex <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v <- gsub("[^A-Z]", "", v)                    # X/Y -> XY, "M " -> M
  out <- rep("U", length(v))
  out[v %in% c("XY", "M", "MALE")]   <- "M"
  out[v %in% c("XX", "F", "FEMALE")] <- "F"
  out[is.na(x) | !nzchar(v)] <- "U"
  factor(out, levels = c("F", "M", "U"))
}

#' Does this column look like a sex call? Judged on its values, not its name,
#' so an oddly-named marker column is still found.
gid_looks_like_sex <- function(v) {
  s <- gid_norm_sex(v)
  known <- sum(s != "U")
  nd <- length(unique(toupper(trimws(as.character(v[!is.na(v)])))))
  known >= 0.5 * sum(!is.na(v)) && nd <= 5
}

#' Collection year from whatever the date column holds: ISO or common text
#' formats, Excel serial numbers stored as text, or a bare four-digit year.
gid_parse_year <- function(x) {
  v <- trimws(as.character(x))
  out <- rep(NA_integer_, length(v))
  ok <- !is.na(v) & nzchar(v)

  bare <- ok & grepl("^(19|20)[0-9]{2}$", v)
  out[bare] <- as.integer(v[bare])

  ser <- ok & !bare & grepl("^[0-9]{5}$", v)
  if (any(ser)) out[ser] <- as.integer(format(
    as.Date(as.numeric(v[ser]), origin = "1899-12-30"), "%Y"))

  rest <- which(ok & is.na(out))
  for (k in rest) {
    y <- regmatches(v[k], regexpr("(19|20)[0-9]{2}", v[k]))
    if (length(y)) { out[k] <- as.integer(y); next }
    for (fmt in c("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%d/%m/%Y")) {
      dd <- suppressWarnings(as.Date(v[k], format = fmt))
      if (!is.na(dd)) { out[k] <- as.integer(format(dd, "%Y")); break }
    }
  }
  out
}

#' Name individuals M1, M2, F1, F2, U1 ... within each sex.
#'
#' Ordered by how often the animal was detected, then by first year, then by
#' cluster id, so the numbering means something and is stable: it is computed
#' from the whole dataset, so filtering the map never renumbers anything.
gid_name_by_sex <- function(individual, sex, year = NULL) {
  ind <- as.character(individual)
  sx  <- as.character(sex)
  tab <- table(ind)
  first <- if (is.null(year)) setNames(rep(NA_integer_, length(tab)), names(tab))
           else vapply(split(as.integer(year), ind), function(z)
                  if (all(is.na(z))) NA_integer_ else as.integer(min(z, na.rm = TRUE)), 1L)

  ## one sex per animal: the majority call among its typed samples
  per <- vapply(split(sx, ind), function(z) {
    z <- z[z != "U"]
    if (!length(z)) "U" else names(sort(table(z), decreasing = TRUE))[1]
  }, "")

  ids <- names(per)
  ord <- order(factor(per[ids], levels = c("F", "M", "U")),
               -as.integer(tab[ids]),
               ifelse(is.na(first[ids]), .Machine$integer.max, first[ids]),
               ids)
  ids <- ids[ord]
  lab <- character(length(ids)); names(lab) <- ids
  n <- c(F = 0L, M = 0L, U = 0L)
  for (i in ids) {
    p <- per[[i]]
    n[p] <- n[p] + 1L
    lab[i] <- paste0(p, n[p])
  }
  list(label = lab, sex = per)
}

#' Sort M1/F2/U3-style names properly: by sex, then numerically, so F10 comes
#' after F9 rather than after F1. Falls back to plain sorting for other labels.
gid_sort_animals <- function(x) {
  x <- unique(x)
  m <- regmatches(x, regexec("^([FMU])([0-9]+)$", x))
  ok <- vapply(m, length, 1L) == 3L
  if (!all(ok)) return(sort(x))
  sx <- vapply(m, `[`, "", 2)
  no <- as.integer(vapply(m, `[`, "", 3))
  x[order(factor(sx, levels = c("F", "M", "U")), no)]
}

#' Inverse UTM (WGS84) -- easting/northing to decimal degrees.
#'
#' Standard inverse transverse Mercator on the WGS84 ellipsoid, series form
#' (Snyder 1987, eqns 8-17 to 8-25). Accurate to well under a metre inside a
#' zone, which is far finer than any GPS fix on a scat.
gid_utm_to_ll <- function(easting, northing, zone, south = FALSE) {
  a <- 6378137; f <- 1 / 298.257223563
  e2 <- f * (2 - f); ep2 <- e2 / (1 - e2); k0 <- 0.9996
  x <- easting - 500000
  y <- if (south) northing - 1e7 else northing
  M  <- y / k0
  mu <- M / (a * (1 - e2/4 - 3*e2^2/64 - 5*e2^3/256))
  e1 <- (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2))
  phi1 <- mu + (3*e1/2 - 27*e1^3/32) * sin(2*mu) +
               (21*e1^2/16 - 55*e1^4/32) * sin(4*mu) +
               (151*e1^3/96) * sin(6*mu) + (1097*e1^4/512) * sin(8*mu)
  C1 <- ep2 * cos(phi1)^2
  T1 <- tan(phi1)^2
  N1 <- a / sqrt(1 - e2 * sin(phi1)^2)
  R1 <- a * (1 - e2) / (1 - e2 * sin(phi1)^2)^1.5
  D  <- x / (N1 * k0)
  lat <- phi1 - (N1 * tan(phi1) / R1) *
    (D^2/2 - (5 + 3*T1 + 10*C1 - 4*C1^2 - 9*ep2) * D^4/24 +
     (61 + 90*T1 + 298*C1 + 45*T1^2 - 252*ep2 - 3*C1^2) * D^6/720)
  lon <- (D - (1 + 2*T1 + C1) * D^3/6 +
          (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*ep2 + 24*T1^2) * D^5/120) / cos(phi1)
  list(lat = lat * 180/pi, lon = (zone * 6 - 183) + lon * 180/pi)
}

#' A colour per individual. Individuals seen more than once get saturated hues
#' spread around the wheel; the order is interleaved so neighbouring labels do
#' not land on neighbouring hues.
gid_ind_colours <- function(inds, n_samples = NULL, grey_singletons = TRUE) {
  inds <- unique(inds)
  multi <- if (is.null(n_samples) || !grey_singletons) inds
           else inds[n_samples[inds] > 1]
  multi <- multi[!is.na(multi)]
  cols <- setNames(rep("#9aa6b2", length(inds)), inds)
  n <- length(multi)
  if (n) {
    h <- seq(12, 372, length.out = n + 1)[seq_len(n)]
    ## interleave so W01 and W02 are not adjacent hues
    ord <- as.vector(t(matrix(c(seq_len(n), rep(NA, (3 - n %% 3) %% 3)), ncol = 3)))
    ord <- ord[!is.na(ord)]
    cols[multi] <- grDevices::hcl(h = h[order(ord)], c = 72, l = 55)[seq_len(n)]
  }
  cols
}

#' Convex hull of a set of points, returned closed. Two points give a segment.
gid_hull <- function(lon, lat) {
  k <- length(lon)
  if (k < 2) return(NULL)
  if (k == 2) return(data.frame(lon = lon, lat = lat))
  h <- grDevices::chull(lon, lat)
  data.frame(lon = lon[c(h, h[1])], lat = lat[c(h, h[1])])
}

# ---------------------------------------------------------------------------- UI
gid_map_tab_ui <- function() {
  nav_panel(
    "Scat map", icon = icon("map-location-dot"),
    uiOutput("run_status_geo"),
    uiOutput("map_no_coords"),
    conditionalPanel(
      "output.has_coords == true",
      uiOutput("map_empty"),
      layout_columns(
        col_widths = c(8, 4),
        card(
          card_header(
            "Where each animal was sampled",
            tags$span(class = "gid-hint", style = "font-weight:400;margin-left:.5rem",
                      "Click a scat for its genotype.")),
          leaflet::leafletOutput("geo_map", height = "560px")),
        tagList(
          card(card_header("Selected scat"), uiOutput("map_detail")),
          card(card_header("Colour key"), uiOutput("map_legend")))),
      card(
        card_header("Mapped samples"),
        DTOutput("tbl_geo"))))
}

#' A clean point map of the same thing the leaflet view is showing.
#'
#' Basemap tiles are deliberately absent: they are copyrighted raster images and
#' cannot be embedded in a vector file, so what gets exported is the figure --
#' points, links, scale bar -- which is what a paper wants anyway.
gid_map_figure <- function(d, cols, who = character(0), style = "none",
                           model = "", label_linked = TRUE) {
  lat0 <- mean(range(d$lat))
  ## one degree of longitude is cos(latitude) as long as one of latitude, so
  ## without this the map is stretched east-west
  asp <- 1 / cos(lat0 * pi / 180)

  p <- ggplot2::ggplot()

  if (length(who) && style != "none") {
    for (ind in who) {
      sset <- d[d$animal == ind, , drop = FALSE]
      if (nrow(sset) < 2) next
      col <- unname(cols[ind])
      if (style == "spider") {
        seg <- data.frame(x = mean(sset$lon), y = mean(sset$lat),
                          xe = sset$lon, ye = sset$lat)
        p <- p + ggplot2::geom_segment(
          data = seg, ggplot2::aes(x = x, y = y, xend = xe, yend = ye),
          colour = col, linewidth = 0.4, alpha = 0.85)
      } else {
        h <- gid_hull(sset$lon, sset$lat)
        if (is.null(h)) next
        p <- p + if (nrow(sset) == 2)
          ggplot2::geom_path(data = h, ggplot2::aes(lon, lat),
                             colour = col, linewidth = 0.5)
        else
          ggplot2::geom_polygon(data = h, ggplot2::aes(lon, lat),
                                fill = col, alpha = 0.18,
                                colour = col, linewidth = 0.45)
      }
    }
  }

  d$fill <- unname(cols[d$animal])
  p <- p +
    ggplot2::geom_point(data = d,
      ggplot2::aes(lon, lat, size = n_samples > 1),
      shape = 21, fill = d$fill, colour = "#33383d", stroke = 0.3) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 2.9), guide = "none")

  if (label_linked && length(who)) {
    cen <- do.call(rbind, lapply(who, function(i) {
      sset <- d[d$animal == i, , drop = FALSE]
      if (!nrow(sset)) NULL else
        data.frame(lon = mean(sset$lon), lat = max(sset$lat), lab = i,
                   stringsAsFactors = FALSE)
    }))
    if (!is.null(cen) && nrow(cen))
      p <- p + ggplot2::geom_text(data = cen, ggplot2::aes(lon, lat, label = lab),
                                  vjust = -0.9, size = 2.5, colour = "#33383d")
  }

  ## scale bar: a round number of km that spans about a fifth of the width
  km_per_deg <- 111.32 * cos(lat0 * pi / 180)
  span_km <- diff(range(d$lon)) * km_per_deg
  nice <- c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)
  bar_km <- nice[which.min(abs(nice - span_km / 5))]
  x0 <- min(d$lon); y0 <- min(d$lat) - diff(range(d$lat)) * 0.06
  bar <- data.frame(x = x0, xe = x0 + bar_km / km_per_deg, y = y0, ye = y0)
  p <- p +
    ggplot2::geom_segment(data = bar, ggplot2::aes(x = x, y = y, xend = xe, yend = ye),
                          linewidth = 0.8, colour = "#33383d") +
    ggplot2::annotate("text", x = x0 + bar_km / km_per_deg / 2, y = y0,
                      label = paste0(bar_km, " km"), vjust = -0.6, size = 2.6,
                      colour = "#33383d")

  n_ind <- length(unique(d$animal))
  p +
    ggplot2::coord_fixed(ratio = asp, clip = "off") +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Samples by individual",
      subtitle = sprintf("%d samples, %d individuals%s%s", nrow(d), n_ind,
                         if (nzchar(model)) paste0(" \u00b7 ", model) else "",
                         if (length(who) && style != "none")
                           sprintf(" \u00b7 %s linking %d animal%s", style,
                                   length(who), if (length(who) == 1) "" else "s")
                         else "")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e8edf2", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "#6b7a8f", size = 7),
      plot.title = ggplot2::element_text(face = "bold", colour = "#1d3557", size = 12),
      plot.subtitle = ggplot2::element_text(colour = "#6b7a8f", size = 8.5),
      plot.margin = ggplot2::margin(8, 12, 8, 8))
}

# ------------------------------------------------------------------------ SERVER
#' @param deps res, best, conf, prep, run_status, and send_file -- the last is
#'   defined inside app.R's server body, which is a child of the environment
#'   these functions live in, so it has to be handed over rather than inherited.
gid_map_server <- function(input, output, session, deps) {

  ## ---- coordinates ---------------------------------------------------------
  coord_cols <- reactive({
    p <- deps$prep(); req(p)
    gid_find_coord_cols(p$df)
  })

  ## Sex and collection date are found by looking at what the values are, not
  ## just what the column is called, so an oddly-named marker column is still
  ## picked up. Both are optional.
  meta_cols <- reactive({
    p <- deps$prep(); req(p)
    df <- p$df; nm <- names(df)
    cand <- setdiff(nm, colnames(p$gt))
    sex <- cand[vapply(cand, function(k)
      gid_looks_like_sex(df[[k]]), TRUE)]
    ## prefer a sex-sounding name when more than one column qualifies
    if (length(sex) > 1) {
      named <- sex[grepl("sex|gender|sry|zfx|omy", tolower(sex))]
      if (length(named)) sex <- named
    }
    yr <- cand[vapply(cand, function(k)
      sum(!is.na(gid_parse_year(df[[k]]))) >= 0.5 * sum(!is.na(df[[k]])) &&
      sum(!is.na(gid_parse_year(df[[k]]))) > 0, TRUE)]
    if (length(yr) > 1) {
      named <- yr[grepl("date|year|coll", tolower(yr))]
      if (length(named)) yr <- named
    }
    list(sex = if (length(sex)) sex[1] else "",
         date = if (length(yr)) yr[1] else "")
  })

  output$map_meta_ui <- renderUI({
    p <- deps$prep(); req(p)
    mc <- meta_cols(); nm <- c("(none)" = "", names(p$df))
    tagList(
      selectInput("map_sex_col", "Sex column", nm, selected = mc$sex),
      selectInput("map_date_col", "Collection date column", nm, selected = mc$date),
      hint("Leave a column as (none) to skip it. With a sex column set, animals ",
           "are named F1, F2, M1, M2 in order of how often each was detected."))
  })

  ## Which pair of columns to use: decimal degrees if present, otherwise UTM.
  ## The user can override both, because a file can carry more than one.
  output$map_coord_ui <- renderUI({
    p <- deps$prep(); req(p)
    cc <- coord_cols(); nm <- c("", names(p$df))
    tagList(
      selectInput("map_lat", "Latitude column",  nm, selected = cc$lat),
      selectInput("map_lon", "Longitude column", nm, selected = cc$lon),
      if (!nzchar(cc$lat) && nzchar(cc$utm_e)) tagList(
        selectInput("map_utm_e", "UTM easting",  nm, selected = cc$utm_e),
        selectInput("map_utm_n", "UTM northing", nm, selected = cc$utm_n),
        numericInput("map_utm_zone", "UTM zone", 8, 1, 60, 1),
        checkboxInput("map_utm_south", "Southern hemisphere", FALSE)))
  })

  ## sample -> lat/lon, restricted to samples that actually entered the analysis
  geo <- reactive({
    p <- deps$prep(); req(p)
    df <- p$df
    ids <- as.character(df[[input$id_col]])
    num <- function(k) suppressWarnings(as.numeric(as.character(df[[k]])))

    lat <- lon <- NULL
    if (nzchar(input$map_lat %||% "") && nzchar(input$map_lon %||% "")) {
      lat <- num(input$map_lat); lon <- num(input$map_lon)
    } else if (nzchar(input$map_utm_e %||% "") && nzchar(input$map_utm_n %||% "")) {
      ll <- gid_utm_to_ll(num(input$map_utm_e), num(input$map_utm_n),
                          input$map_utm_zone %||% 8, isTRUE(input$map_utm_south))
      lat <- ll$lat; lon <- ll$lon
    }
    if (is.null(lat) || all(is.na(lat)) || all(is.na(lon))) return(NULL)

    d <- data.frame(sample = ids, lat = lat, lon = lon, stringsAsFactors = FALSE)
    d <- d[!is.na(d$lat) & !is.na(d$lon), , drop = FALSE]
    ## a replicate file repeats each sample; one position per sample is enough
    d <- d[!duplicated(d$sample), , drop = FALSE]
    keep <- rownames(deps$res()$gt)
    d <- d[d$sample %in% keep, , drop = FALSE]
    if (!nrow(d)) NULL else d
  })

  output$has_coords <- reactive(!is.null(geo()))
  outputOptions(output, "has_coords", suspendWhenHidden = FALSE)
  ## These live in a conditionalPanel, so they would not render until the tab is
  ## first opened -- and a project restored before then would have no widgets to
  ## put its coordinate and sex choices into.
  outputOptions(output, "map_coord_ui", suspendWhenHidden = FALSE)
  outputOptions(output, "map_meta_ui", suspendWhenHidden = FALSE)

  output$run_status_geo <- renderUI(deps$run_status())

  output$map_no_coords <- renderUI({
    if (deps$run_count() == 0) return(NULL)
    if (!is.null(geo())) return(NULL)
    tags$div(
      class = "gid-flag", style = "margin:1rem 0",
      tags$b("No coordinates found in this file. "),
      "The map needs a position for each sample. Add two columns to the file you ",
      "upload and re-run: ", tags$code("Latitude"), " and ", tags$code("Longitude"),
      ", in decimal degrees (for example ", tags$code("55.4821"), " and ",
      tags$code("-132.8194"), "). ",
      "UTM is also read if the columns are named ", tags$code("Easting"), " and ",
      tags$code("Northing"), " — you will be asked for the zone. ",
      tags$br(), tags$br(),
      "One row per sample is enough; if your file has one row per PCR replicate, ",
      "put the same position on every row of a sample. Everything else on this ",
      "tab works the moment those columns are present.")
  })

  ## ---- which model colours the map ----------------------------------------
  observeEvent(deps$res(), {
    r <- deps$res(); req(r)
    ch <- setNames(names(r$methods),
                   vapply(names(r$methods),
                          function(k) GID_METHODS[[k]]$label %||% k, ""))
    updateSelectInput(session, "map_model", choices = ch,
                      selected = if (!is.null(input$map_model) &&
                                     input$map_model %in% ch) input$map_model
                                 else input$method)
  })

  ## assignment under the model chosen for the map (not necessarily the one the
  ## Individuals tab is showing -- that is the point of the picker)
  map_assign <- reactive({
    r <- deps$res(); req(r)
    m <- r$methods[[input$map_model %||% input$method]]
    req(m)
    m$assignment
  })

  pts <- reactive({
    g <- geo(); req(g)
    a <- map_assign()
    g$individual <- a$individual[match(g$sample, a$sample)]
    g <- g[!is.na(g$individual), , drop = FALSE]
    req(nrow(g) > 0)
    tab <- table(a$individual)
    g$n_samples <- as.integer(tab[g$individual])
    cf <- tryCatch(deps$conf(), error = function(e) NULL)
    g$status <- if (!is.null(cf)) as.character(cf$status[match(g$sample, cf$sample)]) else NA
    g$margin <- if (!is.null(cf)) cf$margin[match(g$sample, cf$sample)] else NA

    df  <- deps$prep()$df
    ids <- as.character(df[[input$id_col]])

    ## A replicate file holds several rows per sample and the first of them may
    ## be the reaction that failed, so summarise over a sample's rows rather
    ## than taking whichever row comes first.
    per_sample <- function(col, f) {
      if (!nzchar(col %||% "")) return(NULL)
      by <- split(df[[col]], ids)
      unname(vapply(g$sample, function(k)
        if (is.null(by[[k]])) NA_character_ else f(by[[k]]), ""))
    }
    sx <- per_sample(input$map_sex_col, function(z) {
      z <- as.character(gid_norm_sex(z)); z <- z[z != "U"]
      if (!length(z)) "U" else names(sort(table(z), decreasing = TRUE))[1]
    })
    g$sex <- if (is.null(sx)) factor(rep("U", nrow(g)), levels = c("F", "M", "U"))
             else factor(sx, levels = c("F", "M", "U"))
    yv <- per_sample(input$map_date_col, function(z) {
      y <- gid_parse_year(z); y <- y[!is.na(y)]
      if (!length(y)) NA_character_ else as.character(min(y))
    })
    g$year <- if (is.null(yv)) NA_integer_ else as.integer(yv)

    ## Names come from the whole dataset, before any filter, so hiding animals
    ## never renumbers the ones still on screen.
    if (nzchar(input$map_sex_col %||% "")) {
      nm <- gid_name_by_sex(g$individual, as.character(g$sex), g$year)
      g$animal <- unname(nm$label[g$individual])
      g$animal_sex <- unname(nm$sex[g$individual])
    } else {
      g$animal <- g$individual
      g$animal_sex <- "U"
    }
    g[order(-g$n_samples, g$animal, g$sample), ]
  })

  ## ---- what is actually on the map ----------------------------------------
  ## Filters are applied here and nowhere else, so the map, the legend, the
  ## table and the exported figure can never disagree about what is shown.
  shown <- reactive({
    d <- pts()
    keep <- rep(TRUE, nrow(d))
    if (length(input$map_show_sex))
      keep <- keep & as.character(d$sex) %in% input$map_show_sex
    if (length(input$map_show_year))
      keep <- keep & !is.na(d$year) & as.character(d$year) %in% input$map_show_year
    if (length(input$map_show_who))
      keep <- keep & d$animal %in% input$map_show_who
    d[keep, , drop = FALSE]
  })

  ## Built from the unfiltered set so an animal keeps its colour when others
  ## are hidden.
  pal <- reactive({
    d <- pts()
    n <- tapply(d$sample, d$animal, length)
    gid_ind_colours(d$animal, n, isTRUE(input$map_grey))
  })

  ## Only what is on screen can be linked, and only animals with two or more
  ## samples still visible after filtering.
  linkable <- reactive({
    d <- shown()
    if (!nrow(d)) return(character(0))
    n <- table(d$animal)
    gid_sort_animals(names(n)[n > 1])
  })

  ## Keep the filter menus in step with the data. Choices come from the
  ## unfiltered set, so a filter never removes its own options and strand the
  ## user with an empty menu they cannot undo.
  observe({
    d <- tryCatch(pts(), error = function(e) NULL)
    req(d)
    yrs <- as.character(sort(unique(d$year[!is.na(d$year)])))
    who <- gid_sort_animals(unique(d$animal))

    ## A reopened project's filter choices arrive before these menus have any
    ## options, and selectize discards a selection it has no option for. So the
    ## project leaves them here and they are claimed at the one moment they can
    ## be applied: as the options are created.
    pm <- deps$pending_map()
    want_year <- if (is.null(pm)) NULL else intersect(pm$year, yrs)
    want_who  <- if (is.null(pm)) NULL else intersect(pm$who,  who)
    sel_year <- if (!is.null(want_year)) want_year else intersect(input$map_show_year, yrs)
    sel_who  <- if (!is.null(want_who))  want_who  else intersect(input$map_show_who,  who)

    ## Held, not cleared on first use. This observer re-runs as the restored
    ## settings settle, and releasing the selection before the browser has
    ## echoed it back would let the next run overwrite it with an empty input.
    if (!is.null(pm) &&
        setequal(input$map_show_year %||% character(0), want_year) &&
        setequal(input$map_show_who  %||% character(0), want_who))
      deps$pending_map(NULL)

    updateSelectizeInput(session, "map_show_year", choices = yrs,
                         selected = sel_year, server = FALSE)
    updateSelectizeInput(session, "map_show_who", choices = who,
                         selected = sel_who, server = FALSE)
  })

  observeEvent(input$map_show_reset, {
    updateCheckboxGroupInput(session, "map_show_sex", selected = character(0))
    updateSelectizeInput(session, "map_show_year", selected = character(0))
    updateSelectizeInput(session, "map_show_who", selected = character(0))
  })

  ## Nothing left to draw is a filter result, not an error, so say so rather
  ## than leaving a blank card.
  output$map_empty <- renderUI({
    d <- tryCatch(shown(), error = function(e) NULL)
    if (is.null(d) || nrow(d)) return(NULL)
    tags$div(class = "gid-flag", style = "margin:.6rem 0",
      tags$b("No samples match these filters. "),
      "Widen them, or press ", tags$b("Show everything"), " in the sidebar.")
  })

  ## An empty selection means every animal, not none. Requiring a pick before
  ## the radio does anything made Spider and Polygon look broken: you chose one
  ## and the map did not change.
  link_targets <- reactive({
    who <- input$map_link_who %||% character(0)
    who <- intersect(who, linkable())
    if (length(who)) who else linkable()
  })

  observeEvent(linkable(), {
    lk <- linkable()
    pm <- deps$pending_map()
    sel <- if (!is.null(pm)) intersect(pm$link_who, lk) else
           intersect(input$map_link_who, lk)
    ## deliberately does not clear pm: the filter observer above owns that, and
    ## releases it only once every restored selection has come back from the
    ## browser.
    updateSelectizeInput(session, "map_link_who", choices = lk,
                         selected = sel, server = FALSE)
  })
  observeEvent(input$map_link_none, {
    updateSelectizeInput(session, "map_link_who", selected = character(0))
  })

  ## ---- the map -------------------------------------------------------------
  ## Everything is drawn inside the render rather than pushed through
  ## leafletProxy(). Outputs on a hidden tab are suspended, so proxy messages
  ## sent before the tab was first opened are dropped on the floor -- which left
  ## the map tiled but empty. Drawing here means it is always complete the
  ## moment it appears, at the cost of a redraw when the colouring changes.
  output$geo_map <- leaflet::renderLeaflet({
    d <- shown(); req(nrow(d) > 0)
    cols  <- pal()
    style <- input$map_link_style %||% "none"
    who   <- link_targets()

    m <- leaflet::leaflet() |>
      leaflet::addProviderTiles("Esri.WorldTopoMap", group = "Topographic") |>
      leaflet::addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
      leaflet::addProviderTiles("CartoDB.Positron",  group = "Plain") |>
      leaflet::addLayersControl(
        baseGroups = c("Topographic", "Satellite", "Plain"),
        options = leaflet::layersControlOptions(collapsed = TRUE))

    ## Hold the view the user has panned to, so switching model or linking does
    ## not throw them back to the full extent.
    ## isolate() is load-bearing: reading this reactively and then clearing it
    ## here would invalidate the render that just wrote it, forever.
    rv <- isolate(deps$restored_view())
    if (!is.null(rv)) {
      ctr <- list(lng = rv$lng, lat = rv$lat); zm <- rv$zoom
      deps$restored_view(NULL)          # once only; the user owns the view after
    } else {
      ctr <- isolate(input$geo_map_center); zm <- isolate(input$geo_map_zoom)
    }
    m <- if (!is.null(ctr) && !is.null(zm))
      leaflet::setView(m, ctr$lng, ctr$lat, zm)
    else
      leaflet::fitBounds(m, min(d$lon), min(d$lat), max(d$lon), max(d$lat))

    ## links first, so the scats sit on top of them
    if (length(who) && style != "none") {
      for (ind in who) {
        sset <- d[d$animal == ind, , drop = FALSE]
        if (nrow(sset) < 2) next
        col <- unname(cols[ind])
        if (style == "spider") {
          cx <- mean(sset$lon); cy <- mean(sset$lat)
          for (i in seq_len(nrow(sset)))
            m <- leaflet::addPolylines(m, lng = c(cx, sset$lon[i]),
                                       lat = c(cy, sset$lat[i]),
                                       color = col, weight = 2, opacity = 0.8)
          m <- leaflet::addCircleMarkers(m, lng = cx, lat = cy, radius = 3.5,
                                         color = col, fillColor = col, weight = 1,
                                         fillOpacity = 1,
                                         label = sprintf("%s centre", ind))
        } else {
          h <- gid_hull(sset$lon, sset$lat)
          if (is.null(h)) next
          m <- if (nrow(sset) == 2)
            leaflet::addPolylines(m, lng = h$lon, lat = h$lat, color = col,
                                  weight = 2.5, opacity = 0.85, label = ind)
          else
            leaflet::addPolygons(m, lng = h$lon, lat = h$lat, color = col,
                                 weight = 2, opacity = 0.9, fillColor = col,
                                 fillOpacity = 0.18, label = ind)
        }
      }
    }

    lab <- sprintf("<b>%s</b><br/>%s%s", d$sample, d$animal,
                   ifelse(d$n_samples > 1, sprintf(" (%d samples)", d$n_samples), ""))
    leaflet::addCircleMarkers(
      m, lng = d$lon, lat = d$lat, layerId = d$sample,
      radius = ifelse(d$n_samples > 1, 7, 5),
      color = "#33383d", weight = 1,
      fillColor = unname(cols[d$animal]), fillOpacity = 0.85,
      label = lapply(lab, htmltools::HTML))
  })

  ## If the map is already on screen when a project finishes restoring, no fresh
  ## render is coming, so pan it directly. When it is not yet on screen this does
  ## nothing and the next render picks the view up instead.
  observeEvent(deps$restored_view(), {
    rv <- deps$restored_view()
    req(rv, !is.null(input$geo_map_zoom))
    leaflet::leafletProxy("geo_map") |>
      leaflet::setView(rv$lng, rv$lat, rv$zoom)
    deps$restored_view(NULL)
  }, ignoreNULL = TRUE)

  ## ---- click a scat --------------------------------------------------------
  output$map_detail <- renderUI({
    id <- input$geo_map_marker_click$id
    if (is.null(id)) return(tags$p(class = "gid-hint",
      "Click any scat on the map to see its sample ID, the individual it was ",
      "assigned to, and its full genotype."))
    d <- shown(); row <- d[d$sample == id, , drop = FALSE]
    if (!nrow(row)) return(tags$p(class = "gid-hint", "Sample not on the map."))
    gt <- deps$res()$gt
    g  <- gt[id, ]
    mates <- setdiff(d$sample[d$animal == row$animal[1]], id)
    tagList(
      tags$table(class = "table table-sm gid-kv",
        tags$tr(tags$td(tags$b("Sample")), tags$td(tags$code(id))),
        tags$tr(tags$td(tags$b("Individual")), tags$td(tags$code(row$animal[1]),
                if (!identical(row$animal[1], row$individual[1]))
                  tags$span(class = "gid-hint",
                            sprintf(" (cluster %s)", row$individual[1])))),
        tags$tr(tags$td(tags$b("Sex")),
                tags$td(switch(as.character(row$sex[1]), F = "Female", M = "Male",
                               "not called"))),
        if (!is.na(row$year[1]))
          tags$tr(tags$td(tags$b("Year")), tags$td(row$year[1])),
        tags$tr(tags$td(tags$b("Model")),
                tags$td(GID_METHODS[[input$map_model %||% input$method]]$label)),
        tags$tr(tags$td(tags$b("Samples of this animal")), tags$td(row$n_samples[1])),
        if (!is.na(row$status[1]))
          tags$tr(tags$td(tags$b("Confidence")),
                  tags$td(sprintf("%s (margin %.3g)", row$status[1], row$margin[1]))),
        tags$tr(tags$td(tags$b("Position")),
                tags$td(sprintf("%.5f, %.5f", row$lat[1], row$lon[1]))),
        if (length(mates))
          tags$tr(tags$td(tags$b("Other samples")),
                  tags$td(paste(mates, collapse = ", ")))),
      tags$p(class = "gid-hint", style = "margin-top:.4rem",
             sprintf("Genotype: %d of %d loci called",
                     sum(!is.na(g)), length(g))),
      tags$div(style = "max-height:190px;overflow:auto",
        tags$table(class = "table table-sm gid-geno",
          tags$tbody(lapply(seq_along(g), function(i)
            tags$tr(tags$td(tags$code(names(g)[i])),
                    tags$td(if (is.na(g[i])) tags$span(class = "gid-hint", "—")
                            else tags$code(g[i]))))))))
  })

  output$map_legend <- renderUI({
    d <- shown(); cols <- pal()
    if (!nrow(d)) return(tags$p(class = "gid-hint", "Nothing shown."))
    n <- table(d$animal)
    multi <- names(n)[n > 1]
    multi <- multi[order(-as.integer(n[multi]), multi)]
    if (!length(multi)) return(tags$p(class = "gid-hint",
      "No animal was sampled more than once under this model."))
    tagList(
      tags$p(class = "gid-hint",
             sprintf("%d animals sampled more than once. %d seen once%s.",
                     length(multi), sum(d$n_samples == 1),
                     if (isTRUE(input$map_grey)) ", shown grey" else "")),
      tags$div(class = "gid-legend",
        lapply(multi, function(i) tags$div(
          class = "gid-legend-row",
          tags$span(class = "gid-swatch",
                    style = sprintf("background:%s", unname(cols[i]))),
          tags$span(i),
          tags$span(class = "gid-hint",
                    sprintf(" %d", sum(d$animal == i)))))))
  })

  ## ---- save the map --------------------------------------------------------
  ## Rendered to a temp file and handed over as base64 through the same Blob
  ## path the CSVs use, so it works with no server behind it.
  save_map <- function(fmt) {
    d <- tryCatch(shown(), error = function(e) NULL)
    if (is.null(d) || !nrow(d))
      return(showNotification("Nothing to save yet - run the analysis first.",
                              type = "warning"))
    w <- input$map_fig_width %||% 9
    if (!is.finite(w) || w < 3) w <- 9
    h <- w * 0.78
    fig <- gid_map_figure(d, pal(), link_targets(),
                          input$map_link_style %||% "none",
                          GID_METHODS[[input$map_model %||% input$method]]$label %||% "",
                          label_linked = isTRUE(input$map_fig_labels))

    ext <- fmt; dev_ok <- TRUE
    f <- tempfile(fileext = paste0(".", ext))
    if (fmt == "pdf") {
      grDevices::pdf(f, width = w, height = h, useDingbats = FALSE)
      print(fig); grDevices::dev.off()
    } else {
      ## Not every R build ships a JPEG device -- webR in particular -- so fall
      ## back to PNG rather than failing, and say so.
      if (isTRUE(unname(capabilities("jpeg")))) {
        grDevices::jpeg(f, width = w, height = h, units = "in", res = 300, quality = 95)
      } else if (isTRUE(unname(capabilities("png")))) {
        ext <- "png"; f <- tempfile(fileext = ".png"); dev_ok <- FALSE
        grDevices::png(f, width = w, height = h, units = "in", res = 300)
      } else {
        return(showNotification(
          "This build of R has no raster graphics device. Use the PDF button - it is vector and scales to any size.",
          type = "warning", duration = 10))
      }
      print(fig); grDevices::dev.off()
      if (!dev_ok)
        showNotification("Saved as PNG: this build of R has no JPEG device.",
                         type = "message", duration = 7)
    }

    raw <- readBin(f, "raw", file.info(f)$size)
    unlink(f)
    deps$send_file(sprintf("genoID_map_%s.%s", format(Sys.Date()), ext),
                   jsonlite::base64_enc(raw), b64 = TRUE,
                   type = switch(ext, pdf = "application/pdf",
                                 jpg = "image/jpeg", png = "image/png",
                                 "application/octet-stream"))
  }

  observeEvent(input$dl_map_pdf, save_map("pdf"))
  observeEvent(input$dl_map_jpg, save_map("jpg"))

  output$tbl_geo <- renderDT({
    d <- shown()
    d$margin <- signif(d$margin, 3)
    keep <- c("sample", "animal", "individual", "sex", "year", "n_samples",
              "lat", "lon", "status", "margin")
    dt(d[, intersect(keep, names(d))])
  })
}
