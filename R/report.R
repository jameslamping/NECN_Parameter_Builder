# report.R -- self-contained HTML provenance report (no pandoc needed), and a PDF of it printed by
# headless Chrome / Chromium / Edge (no R packages needed; set NECN_CHROME to the browser if it is not found)

.png64 <- function(p, w = 11, h = 9) {
  f <- tempfile(fileext = ".png")
  ggplot2::ggsave(f, p, width = w, height = h, dpi = 130)
  paste0("data:image/png;base64,", base64enc::base64encode(f))
}
# multi-species figures split into page-sized images of `per` species (3 rows x 4 columns)
.fig_pages <- function(fig, species, per = 12, ...) {
  if (!length(species)) return(NULL)
  n_pages <- ceiling(length(species) / per)
  chunks <- split(species, ceiling(seq_along(species) / ceiling(length(species) / n_pages)))   # balanced pages
  lapply(chunks, function(sp) htmltools::tags$img(class = "page-fig",
    src = .png64(fig(species = sp, ...), 11, 1.2 + 2.3 * ceiling(length(sp) / 4))))
}
.tab_html <- function(d, digits = 2) {
  d <- as.data.frame(d)
  for (j in seq_along(d)) {
    if (is.numeric(d[[j]])) {
      d[[j]][!is.na(d[[j]]) & d[[j]] == 0] <- 0   # no "-0.000"
      x <- formatC(d[[j]], format = "f", digits = if (all(abs(d[[j]] - round(d[[j]])) < 1e-9, na.rm = TRUE)) 0 else digits, big.mark = ",")
      x[is.na(d[[j]])] <- ""
      d[[j]] <- x
    } else d[[j]][is.na(d[[j]])] <- ""
  }
  htmltools::tags$table(class = "t",
    htmltools::tags$thead(htmltools::tags$tr(lapply(names(d), htmltools::tags$th))),
    htmltools::tags$tbody(lapply(seq_len(nrow(d)), function(i) htmltools::tags$tr(lapply(d[i, ], function(v) htmltools::tags$td(as.character(v)))))))
}

#' Write the report
#' @param s list with settings, species_map, sla_tab, plot_lai, fit, final, necn_out, mode, f
write_report <- function(s, path) {
  css <- "body{font:14px/1.45 -apple-system,Helvetica,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#111}
  h1{font-size:24px}h2{font-size:18px;border-top:1px solid #ddd;padding-top:10px;margin-top:26px}
  .t{border-collapse:collapse;font-size:12px;margin:6px 0 14px}.t th,.t td{border-bottom:1px solid #e3e2de;padding:3px 7px;text-align:right}
  .t th:first-child,.t td:first-child{text-align:left}.t th{color:#555}.meta{color:#666;font-size:12px}img{max-width:100%}
  .box{background:#f6f5f2;border-left:4px solid #111;padding:8px 14px;margin:10px 0}pre{white-space:pre-wrap;font-size:11px}
  @page{size:letter landscape;margin:0.45in}
  @media print{body{max-width:none;margin:0;padding:0;font-size:11px}h1{font-size:20px}h2{font-size:15px;break-after:avoid}
    .t{font-size:9px}.t th,.t td{padding:2px 5px}tr,img,.box{break-inside:avoid}img{max-height:7.3in;object-fit:contain}.page-fig{display:block;break-before:page}
    .box{-webkit-print-color-adjust:exact;print-color-adjust:exact}}"
  st <- s$settings
  has_fit <- !is.null(s$fit)
  title <- if (has_fit) "NECN parameters from FIA" else if (!is.null(s$ic)) "NECN initial communities from TreeMap" else "LANDIS run check against FIA"
  doc <- htmltools::tags$html(htmltools::tags$head(htmltools::tags$meta(charset = "utf-8"), htmltools::tags$title(title),
                                                   htmltools::tags$style(css)),
    htmltools::tags$body(
      htmltools::h1(title),
      htmltools::div(class = "meta", sprintf("Built %s with NECN Parameter Builder. FIA folder: %s. States: %s.",
                                             format(Sys.time(), "%Y-%m-%d %H:%M"), st$fia_dir, paste(st$states, collapse = ", "))),
      if (has_fit) .lai_sections(s, st),
      if (!is.null(s$lfit)) .light_section(s),
      if (!is.null(s$ic)) .ic_section(s$ic),
      if (!is.null(s$chk)) .chk_section(s$chk),
      if (!is.null(s$species_map)) htmltools::tagList(htmltools::h2("Species matching"), .tab_html(s$species_map)),
      if (!is.null(s$sla_tab)) htmltools::tagList(htmltools::h2("Specific leaf area used"), .tab_html(s$sla_tab[SPCD %in% expand_species_map(s$species_map)$SPCD])),
      htmltools::h2("Settings"), htmltools::tags$pre(yaml::as.yaml(st))
    ))
  htmltools::save_html(doc, path)
  invisible(path)
}

.lai_sections <- function(s, st) {
  htmltools::tagList(
      htmltools::div(class = "box", htmltools::HTML(paste0(
        "<b>Method.</b> Observed plot LAI = live-tree foliage biomass (FIA DRYBIO_FOLIAGE) &times; species SLA. NECN's curve LAI = MaximumLAI &times; C / (KLAI + C), C = wood biomass &times; 0.47 (g C/m&sup2;), ",
        "is fitted three ways: <i>monospecific</i> (plots where the species holds &ge; ", st$dom_min * 100, "% of live wood; plot LAI), ",
        "<i>mixed</i> (plots where it holds &ge; ", st$share_min * 100, "%; the species' own LAI divided by its wood share, weighted by share &mdash; biased high for minority hardwoods and subcanopy trees), ",
        "and <i>group-pooled</i> (a functional-group curve, with each species' monospecific estimate shrunk toward it by n/(n+", st$n0, ")). ",
        "In each, MaximumLAI = ", st$plateau_q * 100, "th percentile of LAI in the largest ", st$plateau_top * 100, "% of stands (capped at ", st$max_cap, ") and KLAI is a least-squares fit. ",
        "KLAI does not depend on SLA; MaximumLAI does.<br>",
        "<b>Parameterized for:</b> ", if (identical(s$mode, "current_code")) sprintf("NECN 8.3 and earlier (site LAI = sum of every cohort's own LAI). The export applies the growth-neutral rescale f = %.2f: MaximumLAI &times; f, GrowthLAI &divide; f, MinLAI &times; f, CompetitionIndex &divide; f, LightLAIScale and LightLAIAdjust &times; f, LightLAILocation &divide; f.", s$f)
                                else "NECN with stand-level site LAI (each cohort adds MaximumLAI &times; C<sub>cohort</sub> / (KLAI + C<sub>site</sub>)). Values are exported as fitted; GrowthLAI, MinLAI and CompetitionIndex are unchanged. Recheck CompetitionIndex and retune MaximumANPP with single-cell runs."))),
      htmltools::h2("Plot selection"),
      htmltools::p(sprintf("%s plots kept: sampled annual-inventory plots, %s inventory year, at least %d%% forested, one visit per plot location%s%s.",
                           format(nrow(s$plot_lai), big.mark = ","), paste0(st$min_year, "+"), round(st$forest_min * 100),
                           if (isTRUE(st$exclude_disturbed)) ", no recorded disturbance or treatment" else "",
                           if (!is.null(st$region_file)) sprintf(", inside %s (+%s km buffer)", basename(st$region_file), st$buffer_km) else "")),
      htmltools::h2("Parameters"),
      .tab_html(s$final[, .(Species = LandisName, Group = group, Mono_plots = n_mono, Mixed_plots = n_mixed, Confidence = confidence,
                            Method = method_used, MaxLAI_current = MaxLAI_current, KLAI_current = KLAI_current,
                            KLAI_mono = round(KLAI_mono), KLAI_mixed = round(KLAI_mixed), KLAI_pooled = round(KLAI_pooled),
                            MaximumLAI_used = round(MaximumLAI, 2), KLAI_used = round(KLAI), Source = value_source, Note = note)]),
      htmltools::h2("LAI against stand wood"),
      htmltools::p(class = "meta", "Grey = monospecific plots (plot LAI); green = mixed plots (species LAI / wood share). Curves: blue monospecific, green dashed mixed, purple dotted group-pooled. Black vertical = KLAI used; orange dashed = current KLAI."),
      .fig_pages(function(species) fig_species_fit(s$fit, s$final, species),
                 sort(s$fit$params[!is.na(KLAI_mono) | !is.na(KLAI_mixed) | !is.na(KLAI_pooled), LandisName])),
      htmltools::h2("Group curves"), .tab_html(s$fit$groups),
      if (!is.null(s$check)) htmltools::tagList(
        htmltools::h2("Check: site LAI on FIA plots with these parameters"),
        htmltools::p(sprintf("Site LAI that NECN with stand-level LAI computes on %s FIA plots (LANDIS species >= 90%% of live wood) with the final MaximumLAI / KLAI, against observed plot LAI: R2 %.2f, median NECN / observed %.2f, mean difference %+.2f.",
                             format(s$check$stats$plots, big.mark = ","), s$check$stats$r2, s$check$stats$median_ratio, s$check$stats$bias)),
        .tab_html(s$check$bins[, .(Wood_class_Mgha = as.character(wood_class), Plots = plots, Observed_median = observed, Observed_q25 = obs_q25,
                                   Observed_q75 = obs_q75, NECN_stand_level_median = necn)]),
        htmltools::tags$img(src = .png64(fig_lai_check(s$check), 9, 5))))
}

.ic_section <- function(ic) {
  htmltools::tagList(
    htmltools::h2("Initial communities"),
    htmltools::div(class = "box", htmltools::HTML(paste0(
      "TreeMap 2022 resampled onto the ecoregion map's grid; each TreeMap ID's FIA plot supplies live trees (STATUSCD 1) of the LANDIS species, ",
      "grouped by species and 10-year age class. Leaf = DRYBIO_FOLIAGE, wood = DRYBIO_AG (stem, bark, branches), both &times; TPA_UNADJ, in g/m&sup2; of dry biomass. ",
      "Active cells with no TreeMap ID get map code ", ic$empty_code, ", which has no rows in the CSV (they start empty)."))),
    .tab_html(ic$summary[, .(Item = item, Value = value)]),
    htmltools::tags$img(class = "page-fig", src = .png64(fig_ic_biomass(ic), 10, 7.2)),
    htmltools::h2("Initial communities by species"),
    .tab_html(ic$species[, .(Species = SpeciesName, Cells = cells, Landscape_Tg = round(landscape_Mg / 1e6, 3), Cohorts = cohorts,
                             Mean_age = round(mean_age), Max_age = max_age)]),
    htmltools::tags$img(class = "page-fig", src = .png64(fig_ic_check(ic$crosswalk), 10, 7)))
}

#' Path to Chrome, Chromium or Edge (or NECN_CHROME), NA if none
find_chrome <- function() {
  cand <- c(Sys.getenv("NECN_CHROME"), Sys.getenv("CHROME_PATH"),
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/Applications/Chromium.app/Contents/MacOS/Chromium",
            "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
            unname(Sys.which(c("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "microsoft-edge"))),
            "C:/Program Files/Google/Chrome/Application/chrome.exe",
            "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
            "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[1] else NA_character_
}

#' Print an HTML file to PDF with headless Chrome
html_to_pdf <- function(html, pdf, browser = find_chrome(), timeout = 180) {
  if (is.na(browser)) stop("PDF needs Google Chrome, Chromium or Microsoft Edge. Install one, or set NECN_CHROME to its path.")
  pdf <- normalizePath(pdf, winslash = "/", mustWork = FALSE)
  url <- paste0("file://", if (.Platform$OS.type == "windows") "/" else "", normalizePath(html, winslash = "/"))
  prof <- tempfile("chrome-profile"); on.exit(unlink(prof, recursive = TRUE), add = TRUE)
  unlink(pdf)
  system2(browser, c("--headless", "--disable-gpu", "--no-first-run", "--no-pdf-header-footer", "--print-to-pdf-no-header",
                     shQuote(paste0("--user-data-dir=", prof)), shQuote(paste0("--print-to-pdf=", pdf)), shQuote(url)),
          stdout = FALSE, stderr = FALSE, timeout = timeout)
  if (!file.exists(pdf) || file.size(pdf) == 0) stop("Chrome did not produce the PDF.")
  invisible(pdf)
}

#' Write the report as HTML and PDF next to each other
write_report_pdf <- function(s, pdf_path, html_path = sub("\\.pdf$", ".html", pdf_path)) {
  write_report(s, html_path)
  html_to_pdf(html_path, pdf_path)
}

.light_section <- function(s) {
  ls <- s$lfit$settings
  seed <- switch(ls$seed_how, plot = sprintf("a live tree of the species (DBH &ge; %s in) on the plot", ls$seed_dia),
                 radius = sprintf("the species (DBH &ge; %s in) on any selected plot within %s km", ls$seed_dia, ls$seed_km), none = "no condition (all plots)")
  htmltools::tagList(
    htmltools::h2("Establishment light curves"),
    htmltools::div(class = "box", htmltools::HTML(paste0(
      "NECN: P(light) = min(1, Adjust &times; (WeibullPDF(site LAI; Shape, Scale) + Location)). Fitted by maximum likelihood to presence of FIA seedlings of the species ",
      "on plots with a seed source (", seed, "). LAI scale: ", LIGHT_SCALE_LABEL[[ls$scale]], ". ",
      if (isTRUE(ls$normalize)) "Adjust rescaled so each curve peaks at 1 over the observed LAI range (light only limits establishment; climate and seed supply are modelled separately). " else "Adjust kept at the fitted FIA presence rate. ",
      if (identical(ls$bare, "lai05")) paste0("Bare ground (site LAI 0, a cell with no cohorts): where the fitted Weibull curve falls toward 0 at LAI 0 (Shape &gt; 1), the curve was refitted so its value there is at least its value at LAI ", ls$bare_lai, " (about 80% of full light). ") else "Bare-ground values are as fitted. ",
      "Species with fewer than ", ls$min_pres, " plots with seedlings keep their current values unless overridden. Values below are before any export rescale."))),
    htmltools::p("P_bare is the light probability on bare ground (site LAI 0, a cell with no cohorts); it controls how fast disturbed cells regenerate. P_LAI10 is the probability under a closed canopy. Plots_LAI_under_0.5 shows how many fitting plots support the low-LAI end."),
    .tab_html(s$lfinal[, .(Species = LandisName, Plots = n_plots, With_seedlings = n_present, Status = status, Dev_explained = round(dev_expl, 3),
                           Plots_LAI_under_0.5 = if ("n_low_lai" %in% names(s$lfinal)) n_low_lai else NA_integer_,
                           P_bare = if ("p0_used" %in% names(s$lfinal)) round(p0_used, 3) else NA_real_,
                           P_LAI10 = if ("p10_used" %in% names(s$lfinal)) round(p10_used, 3) else NA_real_,
                           Shape = Shape_used, Scale = Scale_used, Location = Location_used, Adjust = Adjust_used, Source = light_source,
                           Shape_now = LightLAIShape, Scale_now = LightLAIScale, Location_now = LightLAILocation, Adjust_now = LightLAIAdjust, Note = note)], 3),
    htmltools::p(class = "meta", "Figures: share of seed-source plots with seedlings in equal-count LAI bins (points, 90% intervals). Blue: fitted presence; black dashed: exported light probability; orange: current table."),
    .fig_pages(function(species) fig_light_fit(s$lfit, s$lfinal, species, ls$scale), sort(s$lfit$params[n_plots > 0, LandisName])))
}

.chk_section <- function(k) {
  cmp <- k$cmp; w <- k$window
  sp <- landscape_species(k$lr, cmp); hh <- chk_fig_height(length(sp), 2.1, 1.4)
  htmltools::tagList(
    htmltools::h2("Check: regeneration in a LANDIS run against FIA"),
    htmltools::div(class = "box", htmltools::HTML(sprintf(paste0(
      "Run: %s. Regeneration = cohorts (LANDIS) or live trees (FIA) up to %d years old. LANDIS: %s cells forested at year 0, years %d&ndash;%d. ",
      "FIA: %s plots within %s km of the landscape, aged with the same source as the initial communities. ",
      "Young biomass (Mg/ha) does not depend on plot or cell size; frequencies do (a 90 m cell samples far more ground than FIA microplots)."),
      htmltools::htmlEscape(k$run$dir), k$lr$regen_age, format(k$lr$n_forest, big.mark = ","), w[1], w[2],
      format(k$fr$n_plots, big.mark = ","), k$fr$buffer_km))),
    if (!is.null(n <- fia_age_note(k$fr))) htmltools::div(class = "box", htmltools::tags$b("Caution. "), n),
    .tab_html(cmp[, .(Species = LandisName, Assessment = assessment, FIA_Mgha = round(fia_Mgha, 3), FIA_lo = round(fia_lo, 3), FIA_hi = round(fia_hi, 3),
                      LANDIS_Mgha = round(landis_Mgha, 3), Ratio = round(ratio, 2), FIA_pct_plots = round(fia_pct_plots, 1), LANDIS_pct_cells = round(landis_pct_cells, 1),
                      FIA_where_present = round(fia_Mgha_where_present, 2), LANDIS_where_present = round(landis_Mgha_where_present, 2))], 3),
    htmltools::tags$img(class = "page-fig", src = .png64(fig_regen_compare(cmp), 11, min(14, 2 + 0.3 * nrow(cmp)))),
    htmltools::tags$img(class = "page-fig", src = .png64(fig_regen_time(k$lr, cmp, sp), 11, hh)),
    if (!is.null(k$cohorts)) htmltools::tags$img(class = "page-fig", src = .png64(fig_cohort_biomass(k$cohorts, sp) +
      ggplot2::labs(title = "Cohort biomass by cohort age", subtitle = "Each species' growth curve. LANDIS above FIA at young ages points to MaximumANPP / KLAI; departures at old ages to mortality and longevity."), 11, hh + 0.5)),
    if (!is.null(k$fr$plot_lai)) htmltools::tagList(
      htmltools::tags$img(class = "page-fig", src = .png64(fig_regen_lai(k$lr, k$fr, w, sp), 11, hh)),
      htmltools::tags$img(src = .png64(fig_landscape_lai(k$lr, k$fr), 10, 3.6))))
}
