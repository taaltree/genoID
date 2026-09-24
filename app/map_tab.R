## map_tab.R -- spatial view of samples, colored by individual ---------------
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

#' Normalize however a lab writes sex into M / F / U.
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
#' (Snyder 1987, eqns 8-17 to 8-25). Accurate to well under a meter inside a
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

#' A color per individual. Individuals seen more than once get saturated hues
#' spread around the wheel; the order is interleaved so neighboring labels do
#' not land on neighboring hues.
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
        ## Stack below a laptop-width screen: beside the sidebar, a third of
        ## the remaining width is too narrow for the scat details to read.
        col_widths = breakpoints(sm = c(12, 12), lg = c(7, 5), xl = c(8, 4)),
        card(
          card_header(
            "Where each animal was sampled",
            tags$span(class = "gid-hint", style = "font-weight:400;margin-left:.5rem",
                      "Click a scat for its genotype.")),
          leaflet::leafletOutput("geo_map", height = "560px")),
        tagList(
          card(card_header("Selected scat"), uiOutput("map_detail")),
          card(card_header("Color key"), uiOutput("map_legend")))),
      card(
        card_header("Mapped samples"),
        DTOutput("tbl_geo"))))
}

#' Write an htmlwidget to one self-contained HTML file.
#'
#' htmlwidgets::saveWidget(selfcontained = TRUE) shells out to pandoc, which
#' does not exist in the WebAssembly build, so the inlining is done here: every
#' script and stylesheet is read off disk and embedded, and any url() a
#' stylesheet points at becomes a data URI. The result opens from a file:// path
#' or an email attachment with nothing else alongside it.
gid_widget_html <- function(widget, title = "genoID map") {
  tags <- htmltools::renderTags(widget)

  mime <- function(f) {
    switch(tolower(tools::file_ext(f)),
           png = "image/png", gif = "image/gif", jpg = , jpeg = "image/jpeg",
           svg = "image/svg+xml", woff = "font/woff", woff2 = "font/woff2",
           ttf = "font/ttf", "application/octet-stream")
  }
  data_uri <- function(path) {
    raw <- readBin(path, "raw", file.info(path)$size)
    sprintf("data:%s;base64,%s", mime(path), jsonlite::base64_enc(raw))
  }
  read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

  head_parts <- character(0)
  for (d in tags$dependencies) {
    base <- d$src$file
    if (is.null(base)) next                      # href-only dependency: skip
    if (!is.null(d$src$package))
      base <- system.file(base, package = d$src$package)

    for (css in d$stylesheet) {
      f <- file.path(base, css)
      if (!file.exists(f)) next
      txt <- read_text(f)
      ## pull in whatever the stylesheet points at, or the marker and layer
      ## icons would be broken links in the shared file
      for (u in unique(regmatches(txt, gregexpr("url\\(([^)]+)\\)", txt))[[1]])) {
        ref <- trimws(gsub("^url\\(|\\)$|[\"']", "", u))
        ref <- sub("[?#].*$", "", ref)
        if (!nzchar(ref) || grepl("^(data:|https?:|//)", ref)) next
        img <- file.path(dirname(f), ref)
        ## a bare url() can resolve to the stylesheet's own directory
        if (file.exists(img) && !dir.exists(img))
          txt <- gsub(u, sprintf("url(%s)", data_uri(img)), txt, fixed = TRUE)
      }
      head_parts <- c(head_parts, sprintf("<style>\n%s\n</style>", txt))
    }
    for (js in d$script) {
      f <- file.path(base, js)
      if (!file.exists(f)) next
      head_parts <- c(head_parts,
                      sprintf("<script>\n%s\n</script>", read_text(f)))
    }
  }

  paste0(
    "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"utf-8\"/>\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>\n",
    "<title>", htmltools::htmlEscape(title), "</title>\n",
    paste(head_parts, collapse = "\n"), "\n",
    "<style>html,body{margin:0;padding:0;height:100%;}",
    ".gid-wrap{height:100%;display:flex;flex-direction:column;",
    "font-family:system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;}",
    ".gid-head{padding:.6rem .9rem;border-bottom:1px solid #e3e9ef;background:#fff;}",
    ".gid-head h1{margin:0;font-size:1rem;color:#1d3557;}",
    ".gid-head p{margin:.15rem 0 0;font-size:.78rem;color:#6b7a8f;}",
    ".gid-map{flex:1;min-height:0;}",
    ".gid-map .leaflet,.gid-map .html-widget{height:100%!important;width:100%!important;}",
    "</style>\n",
    tags$head, "\n</head>\n<body>\n",
    "<div class=\"gid-wrap\">",
    "<div class=\"gid-head\"><h1>", htmltools::htmlEscape(title), "</h1>",
    "<p>Pan, zoom and click a scat for its details. Switch basemaps with the ",
    "control at the top right. Made with genoID.</p></div>",
    "<div class=\"gid-map\">", tags$html, "</div></div>\n",
    "</body>\n</html>\n")
}

#' A map marker as an inline SVG data URI.
#'
#' Shapes follow the pedigree convention every geneticist already reads without
#' a legend: circle female, square male, diamond unknown. Color still carries
#' the individual, so shape and color are independent channels and the map
#' stays readable in grayscale or to a color-blind reader.
gid_marker_svg <- function(shape, fill, size = 17, stroke = "#33383d") {
  h <- size / 2
  r <- size * 0.33
  body <- switch(
    shape,
    F = sprintf('<circle cx="%.2f" cy="%.2f" r="%.2f"', h, h, r),
    M = sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f"',
                h - r * 0.9, h - r * 0.9, r * 1.8, r * 1.8),
    sprintf('<polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f"',
            h, h - r * 1.15, h + r * 1.15, h, h, h + r * 1.15, h - r * 1.15, h))
  svg <- sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">%s fill="%s" stroke="%s" stroke-width="1.1"/></svg>',
    size, size, size, size, body, fill, stroke)
  sprintf("data:image/svg+xml;base64,%s",
          jsonlite::base64_enc(charToRaw(svg)))
}

#' Yardsticks for reading a sample against its closest other animal.
#'
#' "Differs at 4 loci" means nothing on its own: it depends on how many loci
#' there are and how noisy the genotyping was. So the pair is read against two
#' distributions from this same dataset -- how much two samples grouped into
#' ONE animal differ (that is genotyping error at work), and how much each
#' sample's closest OTHER animal differs.
#'
#' @param asg assignment (sample, individual) that defines "one animal" for the
#'   error baseline
#' @param gt genotype matrix
#' @param rivals data.frame(sample, rival_sample) under the map's model
#' @param pts mapped samples (sample, animal, lon, lat), for typical distances
#' @param min_loci pairs typed at fewer loci than this are left out
#' @param baseline how to name the model behind `asg`, in words
gid_rival_context <- function(asg, gt, rivals, pts, min_loci = 1,
                              baseline = "the model") {
  qq <- function(x, p)
    if (length(x)) unname(stats::quantile(x, p, na.rm = TRUE)) else NA_real_

  by <- split(as.character(asg$sample), asg$individual)
  by <- by[lengths(by) > 1]
  pr <- do.call(rbind, lapply(by, function(s) {
    k <- utils::combn(length(s), 2)
    data.frame(a = s[k[1, ]], b = s[k[2, ]], stringsAsFactors = FALSE)
  }))
  ## One huge cluster would dominate the baseline and slow the tab down. Even
  ## thinning keeps the shape without drawing random numbers.
  if (!is.null(pr) && nrow(pr) > 20000)
    pr <- pr[unique(round(seq(1, nrow(pr), length.out = 20000))), , drop = FALSE]
  w <- if (is.null(pr)) numeric(0) else {
    x <- gid_pair_mismatch(gt, pr$a, pr$b)
    x$prop_mismatch[x$n_compared >= min_loci]
  }

  rv <- rivals[!is.na(rivals$rival_sample), , drop = FALSE]
  r <- gid_pair_mismatch(gt, rv$sample, rv$rival_sample)
  r <- r$prop_mismatch[r$n_compared >= min_loci]

  sp <- split(pts[, c("lon", "lat")], pts$animal)
  dd <- unlist(lapply(sp, function(x) {
    if (nrow(x) < 2) return(NULL)
    k <- utils::combn(nrow(x), 2)
    gid_dist_m(x$lon[k[1, ]], x$lat[k[1, ]], x$lon[k[2, ]], x$lat[k[2, ]])
  }), use.names = FALSE)

  list(baseline = baseline,
       n_within = length(w), within_median = qq(w, 0.5), within_q95 = qq(w, 0.95),
       n_rival = length(r), rival_median = qq(r, 0.5),
       typical_m = qq(dd, 0.5), range_m = qq(dd, 0.95),
       dist_sorted = sort(dd[is.finite(dd)]))
}

#' The closest other animal: how different it is genetically, and how far away
#' it was found.
#'
#' A false split leaves a signature two ways at once: the sample's closest
#' genetic match sits in a DIFFERENT individual, and that sample was picked up
#' nearby. Either alone is weak -- animals have neighbors, and relatives look
#' alike -- but together they are the pattern a split produces, because both
#' halves came off the same animal in the same place.
#'
#' Both yardsticks come from this dataset rather than fixed numbers: genetic
#' difference is read against how much two samples of one animal differ here,
#' and distance against how far apart samples of one animal usually are, which
#' depends on home-range size and how the survey was walked.
#'
#' @param row the clicked sample: one row of the map's point table
#' @param all_pts every mapped sample, to find the rival's position
#' @param ctx gid_rival_context() output
#' @param mm gid_pair_mismatch() for the sample and its rival
#' @param ev the model's own score for the pair: list(scale, posterior,
#'   log10_lr, log10_lambda), any of which may be NULL or NA
#' @param kin_label what the likelihood ratio is weighed against, in words
gid_rival_block <- function(row, all_pts, ctx, mm, ev = list(),
                            post_cut = NA_real_, lambda_cut = NA_real_,
                            kin_label = "full siblings") {
  fmt_d <- function(m) if (!is.finite(m)) "unknown" else
    if (m < 1000) sprintf("%.0f m", m) else sprintf("%.1f km", m / 1000)
  pct <- function(p) if (!is.finite(p)) "?" else sprintf("%.0f%%", 100 * p)
  ## Rounded down near 1, so a posterior just short of the cutoff never
  ## prints as the cutoff itself.
  fmt_p <- function(p) if (!is.finite(p)) "not available" else
    if (p >= 0.99) sprintf("%.4f", floor(p * 1e4) / 1e4)
    else if (p >= 0.001) sprintf("%.3f", p)
    else if (p > 0) sprintf("%.1e", p) else "below 1e-300"
  chance <- function(p) sprintf(if (p >= 0.1) "%.0f%%" else "%.1f%%", 100 * p)
  num <- function(x) suppressWarnings(as.numeric(x %||% NA)[1])
  kv <- function(k, ...) tags$tr(tags$td(tags$b(k)), tags$td(...))
  hint <- function(...) tags$span(class = "gid-hint", ...)

  rv <- row$rival_sample[1]
  if (is.null(rv) || is.na(rv) || !nzchar(rv))
    return(tags$div(class = "gid-hint", style = "margin-top:.5rem",
      tags$b("No other animal to compare. "),
      "No sample assigned to another animal shares enough typed loci with this ",
      "one to be compared, so a split cannot be checked from here."))

  b <- all_pts[all_pts$sample == rv, , drop = FALSE]
  d_m <- if (nrow(b)) gid_dist_m(row$lon[1], row$lat[1], b$lon[1], b$lat[1]) else NA_real_
  typical <- num(ctx$typical_m); range_m <- num(ctx$range_m)

  k <- num(mm$n_mismatch); n <- num(mm$n_compared); p <- num(mm$prop_mismatch)
  hard <- num(mm$n_mismatch_2allele)
  q95 <- if (isTRUE(ctx$n_within >= 10)) num(ctx$within_q95) else NA_real_
  post <- num(ev$posterior); lr <- num(ev$log10_lr); lam <- num(ev$log10_lambda)
  scale <- as.character(ev$scale %||% NA)

  ## Genetically close: no more different than genotyping error makes two
  ## samples of one animal in this dataset, or -- with too few animals sampled
  ## twice to learn that -- two loci or fewer. A model that still gives the
  ## pair a real chance of being one animal counts too.
  by_loci  <- is.finite(p) && (if (is.finite(q95)) p <= q95 else k <= 2)
  by_model <- identical(scale, "posterior") && is.finite(post) && post > 0.01
  gen_close <- by_loci || by_model
  ## Where the distance falls among same-animal pairs: the share of them found
  ## at least this far apart. Only the farthest 5% counts against a split -- a
  ## median cut would call half of all true pairs "far" -- and with no animal
  ## mapped twice there is nothing to judge distance by, so it is not used.
  dd <- ctx$dist_sorted
  F_d <- if (is.finite(d_m) && length(dd)) mean(dd <= d_m) else NA_real_
  spa_close <- is.finite(d_m) && (d_m < 1000 || !is.finite(F_d) || F_d <= 0.95)

  diff_txt <- if (is.finite(p)) sprintf("%d of %d loci typed in both (%s)", k, n, pct(p))
              else "no locus typed in both"
  why_close <- if (by_loci && k == 0)
    sprintf("They match at all %d loci typed in both. ", n)
  else if (by_loci)
    sprintf("They differ at only %s, %s. ", diff_txt,
            if (is.finite(q95))
              sprintf("no more than genotyping error makes two samples of one animal differ here (up to %s)",
                      pct(q95))
            else "few enough for genotyping error alone to explain")
  else sprintf("The model still gives a %s chance that they are the same animal. ",
               chance(post))
  held_apart <- if (identical(scale, "posterior") && is.finite(post) &&
                    is.finite(post_cut) && post >= 0.5 && post < post_cut)
    sprintf("The model kept them apart only because %s falls short of the %s it needs to join them. ",
            fmt_p(post), post_cut)
  advice <- if (is.finite(k) && k == 0 && !is.null(held_apart))
      paste0("More loci would settle it. If many samples look like this, the cutoff ",
             "is beyond what this panel can reach; see \u201cChoosing your posterior ",
             "cutoff\u201d on the Individuals tab.")
    else if (is.finite(k) && k == 0) "More loci would settle it."
    else if (is.finite(k))
      sprintf("Re-genotype both at the %s where they differ before reporting them as separate animals.",
              if (k == 1) "locus" else sprintf("%d loci", k))
  where_txt <- if (!is.finite(F_d)) ""
    else if (F_d <= 0.5)
      sprintf(", closer than most same-animal pairs here (%s of them are farther apart)", pct(1 - F_d))
    else sprintf(", within the range one animal covers here (%s of same-animal pairs are farther apart)",
                 pct(1 - F_d))

  verdict <- if (gen_close && spa_close)
    tags$div(class = "gid-flag", style = "margin-top:.4rem;font-size:.82rem",
      tags$b("Worth a second look. "), why_close, held_apart,
      sprintf("They were found %s apart%s. ", fmt_d(d_m), where_txt),
      if (isTRUE(F_d <= 0.5) || d_m < 1000)
        "That is the pattern one animal split in two leaves. ",
      advice)
  else if (gen_close && !is.finite(d_m))
    tags$div(class = "gid-flag", style = "margin-top:.4rem;font-size:.82rem",
      tags$b("Worth a second look. "), why_close, held_apart,
      "The other sample has no coordinates, so distance cannot help decide. ", advice)
  else if (gen_close)
    tags$div(class = "gid-hint", style = "margin-top:.4rem",
      why_close, held_apart,
      sprintf("But they were found %s apart, farther than 95%% of same-animal pairs here (within %s). ",
              fmt_d(d_m), fmt_d(range_m)),
      "That makes a split less likely, though a dispersing animal can travel that far.")
  else if (!is.finite(p))
    tags$div(class = "gid-hint", style = "margin-top:.4rem",
      "These two could not be compared locus by locus.")
  else
    tags$div(class = "gid-hint", style = "margin-top:.4rem",
      tags$b("Different animals. "),
      sprintf("They differ at %s", diff_txt),
      if (is.finite(q95))
        sprintf(", more than genotyping error makes two samples of one animal differ here (up to %s)",
                pct(q95))
      else ", more than the two or fewer that would suggest a split",
      ", so the distance between them does not matter.",
      if ((isTRUE(F_d <= 0.5) || isTRUE(d_m < 1000)) && isTRUE(p < ctx$rival_median))
        paste0(" They are more alike than most samples are to their closest other ",
               "animal, and were found close together, which fits a relative ",
               "sharing the area."))

  scale_note <- if (isTRUE(ctx$n_within >= 10))
    tags$p(class = "gid-hint", style = "margin:.35rem 0 0",
      sprintf(paste0("For scale: same-animal pairs (two samples %s places in one animal) ",
                     "differ at a median of %s of loci, and 95%% of them at %s or less; that is ",
                     "what genotyping error does here. A sample's closest other animal differs ",
                     "at a median of %s."),
              ctx$baseline %||% "the model", pct(ctx$within_median), pct(ctx$within_q95),
              pct(ctx$rival_median)))
  else
    tags$p(class = "gid-hint", style = "margin:.35rem 0 0",
      "Too few animals were sampled more than once to measure how much genotyping ",
      "error separates two samples of one animal here, so a difference at two loci ",
      "or fewer is taken as the mark of a possible split.")

  who <- if (nrow(b)) b$animal[1] else row$rival_individual[1] %||% NA
  tagList(
    tags$div(class = "gid-label", style = "margin-top:.7rem", "Closest other animal"),
    tags$table(class = "table table-sm gid-kv",
      kv("Sample", tags$code(rv), if (!is.na(who)) hint(sprintf(" (%s)", who))),
      kv("Loci that differ", diff_txt,
         if (is.finite(hard) && hard > 0)
           hint(sprintf("; no allele shared at %d", hard))),
      if (is.finite(post))
        kv("Probability same animal", fmt_p(post),
           if (is.finite(post_cut)) hint(sprintf(" (match at %s)", post_cut))),
      if (is.finite(lr))
        kv("Likelihood ratio, log₁₀", sprintf("%.1f", lr),
           hint(sprintf(" (same animal vs. %s)", kin_label))),
      if (is.finite(lam))
        kv("Sethi Λ, log₁₀", sprintf("%.1f", lam),
           if (is.finite(lambda_cut) && lambda_cut > 0)
             hint(sprintf(" (match above %.1f)", log10(lambda_cut)))),
      kv("Distance apart", fmt_d(d_m),
         if (is.finite(F_d))
           hint(sprintf("; %s of same-animal pairs are farther apart", pct(1 - F_d)))),
      if (is.finite(typical))
        kv("Same-animal pairs", sprintf("typically %s apart", fmt_d(typical)),
           if (is.finite(range_m)) hint(sprintf("; 95%% within %s", fmt_d(range_m))))),
    scale_note,
    verdict,
    if (!nrow(b)) tags$p(class = "gid-hint",
      "That sample has no coordinates, so the two cannot be compared in space."))
}

#' Build the interactive map. Shared by the on-screen view and the exported
#' HTML file, so the two cannot drift apart.
#'
#' @param view NULL to fit the data, or list(lng, lat, zoom) to open somewhere
gid_leaflet_map <- function(d, cols, who = character(0), style = "none",
                            view = NULL) {
  m <- leaflet::leaflet() |>
    leaflet::addProviderTiles("Esri.WorldTopoMap", group = "Topographic") |>
    leaflet::addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
    leaflet::addProviderTiles("CartoDB.Positron",  group = "Plain") |>
    leaflet::addLayersControl(
      baseGroups = c("Topographic", "Satellite", "Plain"),
      options = leaflet::layersControlOptions(collapsed = TRUE))

  m <- if (!is.null(view))
    leaflet::setView(m, view$lng, view$lat, view$zoom)
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
                                       label = sprintf("%s center", ind))
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

  ## Shape carries sex, color carries the individual. Icons are built once per
  ## distinct shape/color/size combination rather than once per sample, so a
  ## few thousand scats do not become a few thousand data URIs.
  sx  <- as.character(d$sex); sx[is.na(sx)] <- "U"
  fil <- unname(cols[d$animal])
  siz <- ifelse(d$n_samples > 1, 19L, 14L)
  key <- paste(sx, fil, siz, sep = "|")
  uni <- !duplicated(key)
  lut <- setNames(vapply(which(uni), function(i)
    gid_marker_svg(sx[i], fil[i], siz[i]), ""), key[uni])

  leaflet::addMarkers(
    m, lng = d$lon, lat = d$lat, layerId = d$sample,
    icon = leaflet::icons(iconUrl = unname(lut[key]),
                          iconWidth = siz, iconHeight = siz,
                          iconAnchorX = siz / 2, iconAnchorY = siz / 2),
    label = lapply(lab, htmltools::HTML),
    popup = sprintf(
      "<b>%s</b><br/>Individual: <b>%s</b><br/>Samples of this animal: %d%s%s",
      d$sample, d$animal, d$n_samples,
      ifelse(is.na(d$year), "", sprintf("<br/>Year: %s", d$year)),
      ifelse(is.na(d$sex) | d$sex == "U", "",
             sprintf("<br/>Sex: %s", ifelse(d$sex == "F", "Female", "Male")))))
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
  ## 21 circle, 22 square, 23 diamond -- the same pedigree convention as the
  ## interactive map, and all three take a fill and a border
  sx <- as.character(d$sex); sx[is.na(sx)] <- "U"
  d$shape <- c(F = 21L, M = 22L, U = 23L)[sx]
  p <- p +
    ggplot2::geom_point(data = d,
      ggplot2::aes(lon, lat, size = n_samples > 1),
      shape = d$shape, fill = d$fill, colour = "#33383d", stroke = 0.3) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.9, `TRUE` = 3.1), guide = "none")

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
      caption = if (any(sx != "U"))
        "Circle female \u00b7 square male \u00b7 diamond sex not called" else NULL,
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
      plot.caption = ggplot2::element_text(colour = "#6b7a8f", size = 7.5, hjust = 0),
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
    ## prep(), not res(): the analysis can now ask geo() for coordinates, and
    ## depending on res() here would put the two reactives in a cycle.
    keep <- rownames(deps$prep()$gt)
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

  ## ---- which model colors the map ----------------------------------------
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

  ## Confidence under the model coloring the map. The Individuals tab's table is
  ## reused when both show the same model; otherwise it is worked out again,
  ## because a rival under one model can be a cluster-mate under another.
  map_conf <- reactive({
    r <- deps$res(); req(r)
    key <- input$map_model %||% input$method
    if (identical(key, input$method))
      return(tryCatch(deps$conf(), error = function(e) NULL))
    tryCatch(gid_method_confidence(
      r, key, pid_tab = tryCatch(deps$pid(), error = function(e) NULL),
      post_cut = input$post_cut, min_loci = input$min_loci),
      error = function(e) NULL)
  })

  pts <- reactive({
    g <- geo(); req(g)
    a <- map_assign()
    g$individual <- a$individual[match(g$sample, a$sample)]
    g <- g[!is.na(g$individual), , drop = FALSE]
    req(nrow(g) > 0)
    tab <- table(a$individual)
    g$n_samples <- as.integer(tab[g$individual])
    cf <- map_conf()
    g$status <- if (!is.null(cf)) as.character(cf$status[match(g$sample, cf$sample)]) else NA
    g$margin <- if (!is.null(cf)) cf$margin[match(g$sample, cf$sample)] else NA
    ## The strongest link to a sample assigned to a DIFFERENT animal. If that
    ## rival is also right next door, the two may be one animal split in two.
    k <- if (!is.null(cf)) match(g$sample, cf$sample) else NA
    g$rival        <- if (!is.null(cf)) cf$rival[k] else NA_real_
    g$rival_sample <- if (!is.null(cf)) as.character(cf$rival_sample[k]) else NA_character_
    g$rival_individual <- if (!is.null(cf)) as.character(cf$rival_individual[k]) else NA_character_
    attr(g, "scale") <- if (!is.null(cf)) attr(cf, "scale") else NA_character_

    df  <- deps$prep()$df
    ids <- as.character(df[[input$id_col]])

    ## A replicate file holds several rows per sample and the first of them may
    ## be the reaction that failed, so summarize over a sample's rows rather
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

  ## Yardsticks for the rival panel, worked out once per run and model rather
  ## than on every click. "One animal" for the genotyping-error baseline comes
  ## from the likelihood-ratio model whatever colors the map: it is the one
  ## model built to expect error, so its clusters show what error does, where
  ## an exact-match cluster differs at no locus by construction.
  rival_ctx <- reactive({
    r <- deps$res(); req(r)
    d <- pts(); cf <- map_conf()
    lr <- r$methods$probabilistic
    base <- if (!is.null(lr)) lr$assignment else map_assign()
    rivals <- if (is.null(cf)) data.frame(sample = character(0), rival_sample = character(0))
              else cf[, c("sample", "rival_sample")]
    gid_rival_context(base, r$gt, rivals, d, min_loci = input$min_loci %||% 1,
                      baseline = if (!is.null(lr)) "the likelihood-ratio model"
                                 else "this model")
  })

  ## Built from the unfiltered set so an animal keeps its color when others
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
  ## moment it appears, at the cost of a redraw when the coloring changes.
  output$geo_map <- leaflet::renderLeaflet({
    d <- shown(); req(nrow(d) > 0)

    ## Hold the view the user has panned to, so switching model or linking does
    ## not throw them back to the full extent.
    ## isolate() is load-bearing: reading this reactively and then clearing it
    ## here would invalidate the render that just wrote it, forever.
    rv <- isolate(deps$restored_view())
    if (!is.null(rv)) {
      view <- rv
      deps$restored_view(NULL)          # once only; the user owns the view after
    } else {
      ctr <- isolate(input$geo_map_center); zm <- isolate(input$geo_map_zoom)
      view <- if (!is.null(ctr) && !is.null(zm))
        list(lng = ctr$lng, lat = ctr$lat, zoom = zm) else NULL
    }

    gid_leaflet_map(d, pal(), link_targets(),
                    input$map_link_style %||% "none", view = view)
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

  ## A dashed line to the nearest rival, drawn on click. Only reachable while
  ## the map is on screen, so leafletProxy() is safe here.
  observeEvent(input$geo_map_marker_click, {
    id <- input$geo_map_marker_click$id
    p  <- leaflet::leafletProxy("geo_map") |> leaflet::clearGroup("selection")
    d  <- tryCatch(pts(), error = function(e) NULL)
    if (is.null(d) || is.null(id)) return(invisible(p))
    a <- d[d$sample == id, , drop = FALSE]
    if (!nrow(a)) return(invisible(p))
    p <- leaflet::addCircleMarkers(
      p, lng = a$lon[1], lat = a$lat[1], radius = 13, group = "selection",
      color = "#1d3557", weight = 2, fill = FALSE)
    rv <- a$rival_sample[1]
    if (!is.na(rv) && rv %in% d$sample) {
      b <- d[d$sample == rv, , drop = FALSE]
      p <- p |>
        leaflet::addPolylines(lng = c(a$lon[1], b$lon[1]), lat = c(a$lat[1], b$lat[1]),
                              color = "#c1502e", weight = 2, opacity = 0.9,
                              dashArray = "6,6", group = "selection") |>
        leaflet::addCircleMarkers(lng = b$lon[1], lat = b$lat[1], radius = 10,
                                  group = "selection", color = "#c1502e",
                                  weight = 2, fill = FALSE,
                                  label = sprintf("%s: closest other animal", rv))
    }
    invisible(p)
  })

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
      local({
        rv <- row$rival_sample[1]
        if (is.na(rv %||% NA)) return(gid_rival_block(row, pts(), list(), NULL))
        prs <- deps$res()$methods[[input$map_model %||% input$method]]$pairs
        hit <- if (is.null(prs) || !nrow(prs)) NA_integer_ else
          which((prs$id1 == id & prs$id2 == rv) | (prs$id1 == rv & prs$id2 == id))[1]
        pick <- function(col) if (is.na(hit) || is.null(prs[[col]])) NA_real_ else prs[[col]][hit]
        gid_rival_block(
          row, pts(), rival_ctx(), gid_pair_mismatch(gt, id, rv),
          ev = list(scale = attr(map_conf(), "scale"),
                    posterior = pick("posterior_same"), log10_lr = pick("log10_LR"),
                    log10_lambda = pick("log10_lambda")),
          post_cut = input$post_cut %||% NA_real_,
          lambda_cut = input$lambda_cut %||% NA_real_,
          kin_label = switch(input$kinship %||% "full_sib", half_sib = "half siblings",
                             unrelated = "unrelated animals", "full siblings"))
      }),
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
    shapes <- sort(unique(as.character(d$sex[!is.na(d$sex)])))
    tagList(
      tags$p(class = "gid-hint",
             sprintf("%d animals sampled more than once. %d seen once%s.",
                     length(multi), sum(d$n_samples == 1),
                     if (isTRUE(input$map_grey)) ", shown gray" else "")),
      if (length(shapes) > 1) tags$div(
        class = "gid-legend-row", style = "margin-bottom:.5rem;gap:.8rem",
        lapply(shapes, function(k) tags$span(
          class = "gid-legend-row", style = "gap:.3rem",
          tags$img(src = gid_marker_svg(k, "#8fa3a6", 15), width = 15, height = 15),
          tags$span(class = "gid-hint",
                    switch(k, F = "female", M = "male", "sex not called"))))),
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
      return(showNotification("Nothing to save yet. Run the analysis first.",
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
          "This build of R has no raster graphics device. Use the PDF button instead: it is a vector file and scales to any size.",
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

  ## The interactive map as one shareable file. Collaborators open it in any
  ## browser and pan, zoom and click it exactly as here -- no genoID, no R, no
  ## server. Basemap tiles still come from the internet; the samples, the links
  ## and every popup are inside the file.
  map_html <- function() {
    d <- shown()
    if (!nrow(d)) return(NULL)
    ttl <- sprintf("%s \u00b7 %d samples, %d individuals",
                   deps$source_name() %||% "genoID map",
                   nrow(d), length(unique(d$animal)))
    view <- if (!is.null(input$geo_map_center) && !is.null(input$geo_map_zoom))
      list(lng = input$geo_map_center$lng, lat = input$geo_map_center$lat,
           zoom = input$geo_map_zoom) else NULL
    gid_widget_html(
      gid_leaflet_map(d, pal(), link_targets(),
                      input$map_link_style %||% "none", view = view), ttl)
  }

  observeEvent(input$dl_map_html, {
    h <- tryCatch(map_html(), error = function(e) NULL)
    if (is.null(h))
      return(showNotification("Nothing to share yet. Run the analysis first.",
                              type = "warning"))
    deps$send_file(sprintf("genoID_map_%s.html", format(Sys.Date())), h,
                   type = "text/html;charset=utf-8")
    showNotification(paste("Saved an interactive map. Send that one file to",
                           "anyone; it opens in any browser."),
                     type = "message", duration = 8)
  })

  output$tbl_geo <- renderDT({
    d <- shown()
    d$margin <- signif(d$margin, 3)
    keep <- c("sample", "animal", "individual", "sex", "year", "n_samples",
              "lat", "lon", "status", "margin")
    dt(d[, intersect(keep, names(d))])
  })

  ## Handed back so the Download tab can bundle the mapped table and the
  ## shareable map without reaching into this module's internals.
  invisible(list(
    table     = function() tryCatch(shown(),   error = function(e) NULL),
    coords    = function() tryCatch(geo(),     error = function(e) NULL),
    map_html  = function() tryCatch(map_html(), error = function(e) NULL)))
}
