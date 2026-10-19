# pipeline.R -- the whole build without the app, driven by a config.yml (see config_example.yml)

default_settings <- function() list(
  fia_dir = "", states = character(), cache_dir = "cache",
  forest_min = 0.9, annual_only = TRUE, min_year = 2000, exclude_disturbed = FALSE,
  region_file = NULL, buffer_km = 0,
  species_file = "", necn_table = "", species_overrides = NULL, sla_overrides = NULL,
  dom_min = 0.8, share_min = 0.3, n0 = 30, plateau_top = 0.2, plateau_q = 0.9, max_cap = 20, min_n = 12, nboot = 200, seed = 1,
  sla_softwood = 5, sla_hardwood = 12, download = FALSE, download_method = "rFIA",
  choices_file = NULL, mode = "current_code", rescale_f = 0.47, out_dir = "output",
  light_fit = FALSE, light_scale = NULL, age_dir = "", seed_source = "plot", seed_min_dia = 5, seed_radius_km = 5,
  light_min_pres = 20, light_normalize = TRUE, light_bare = "lai05", light_choices_file = NULL,
  check_run_dir = NULL, check_mode = NULL, check_regen_age = 20, check_every = 5, check_buffer_km = 25, check_window = NULL,
  ic_build = FALSE, ic_template = "", treemap_raster = "", treemap_table = NULL, ic_species_file = "",
  ic_resample = "near", ic_codes = "tm_id", ic_empty_code = 1, ic_min_biomass = 1, ic_only = FALSE, report_pdf = TRUE
)

# HTML report, plus a PDF of it when report_pdf is true and Chrome / Chromium / Edge is available
.write_reports <- function(r, s, stamp, say) {
  h <- file.path(s$out_dir, sprintf("NECN_builder_report_%s.html", stamp))
  write_report(r, h)
  if (isTRUE(s$report_pdf)) tryCatch(html_to_pdf(h, sub("\\.html$", ".pdf", h)),
                                     error = function(e) say("PDF report skipped: ", conditionMessage(e)))
}

run_builder <- function(config, verbose = TRUE) {
  s <- utils::modifyList(default_settings(), if (is.character(config)) yaml::read_yaml(config) else config)
  if (!s$mode %in% c("stand_level", "current_code")) stop("mode must be stand_level (NECN with stand-level site LAI) or current_code (NECN 8.3 and earlier)")
  say <- function(...) if (verbose) message(...)
  ref <- load_ref_species(); sref <- load_sla_reference()

  if (isTRUE(s$download)) { say("Downloading FIA ..."); fia_download(s$states, s$fia_dir, s$download_method) }
  say("Reading FIA ...")
  fia <- fia_read(s$fia_dir, s$states, s$cache_dir)
  region <- if (!is.null(s$region_file)) sf::st_read(s$region_file, quiet = TRUE) else NULL
  plots <- select_plots(fia, s$forest_min, s$annual_only, s$min_year, s$exclude_disturbed, region, s$buffer_km)
  say(nrow(plots), " plots selected")

  if (isTRUE(s$ic_only)) s$ic_build <- TRUE
  spp <- read_landis_species(s$species_file)
  map <- match_species(spp, ref)
  if (!is.null(s$species_overrides)) {        # csv: LandisName, SPCD and/or group
    ov <- data.table::fread(s$species_overrides, colClasses = list(character = "SPCD"))
    if ("SPCD" %in% names(ov)) map[ov[!is.na(SPCD) & nzchar(SPCD)], on = "LandisName", `:=`(SPCD = i.SPCD, status = "user override")]
    map <- add_default_groups(map, ref)
    if ("group" %in% names(ov)) map[ov[!is.na(group) & nzchar(group)], on = "LandisName", group := i.group]
  }
  map <- add_default_groups(map, ref)

  ic <- NULL
  if (isTRUE(s$ic_build)) {
    say("Initial communities from TreeMap ...")
    ic <- build_ic(s$ic_template, s$treemap_raster, fia, map, fit_age_models(fia$TREE[STATUSCD == 1], ref), read_lu_ages(s$age_dir, s$states),
                   treemap_table = s$treemap_table, longevity = read_longevity(s$ic_species_file), resample = s$ic_resample,
                   codes = s$ic_codes, empty_code = as.integer(s$ic_empty_code), min_biomass = s$ic_min_biomass)
    f <- write_ic(ic, s$out_dir, paste0("IC_", format(Sys.Date(), "%Y%m%d")))
    if (verbose) print(ic$summary)
    say("Wrote ", paste(basename(f), collapse = ", "))
    if (isTRUE(s$ic_only)) {
      .write_reports(list(settings = s, species_map = map, ic = ic, mode = s$mode, f = s$rescale_f), s, format(Sys.Date(), "%Y%m%d"), say)
      return(invisible(list(ic = ic, species_map = map)))
    }
  }
  ovs <- if (!is.null(s$sla_overrides)) data.table::fread(s$sla_overrides) else NULL
  sla_tab <- sla_lookup(unique(fia$TREE[PLT_CN %in% plots$PLT_CN, SPCD]), ref, sref, ovs, c(S = s$sla_softwood, H = s$sla_hardwood))

  say("Plot LAI ...")
  lai <- compute_plot_lai(fia, plots, sla_tab)
  sp_tab <- species_plot_table(lai, map)
  say("Fitting MaximumLAI and KLAI ...")
  fit <- fit_lai_params(sp_tab, species_groups(map), s$dom_min, s$share_min, s$plateau_top, s$plateau_q,
                        s$max_cap, s$min_n, s$nboot, s$n0, s$seed)

  necn <- data.table::fread(s$necn_table)
  choices <- if (!is.null(s$choices_file)) data.table::fread(s$choices_file) else NULL
  final <- final_values(fit$params, choices, necn)

  lfit <- lfinal <- NULL
  if (isTRUE(s$light_fit)) {
    scale <- if (is.null(s$light_scale)) s$mode else s$light_scale
    say("Establishment light curves (", scale, ") ...")
    xt <- lai$plots[, .(PLT_CN, x = LAI)]
    if (scale == "stand_level") {
      xt <- merge(xt[, .(PLT_CN)], stand_site_lai(sp_tab, final[, .(LandisName, MaximumLAI, KLAI)])[, .(PLT_CN, x = LAI_necn)], by = "PLT_CN", all.x = TRUE)
      xt[is.na(x), x := 0]
    } else if (scale == "current_code") {
      co <- plot_cohorts(fia, lai$plots$PLT_CN, map, fit_age_models(fia$TREE[STATUSCD == 1], ref), read_lu_ages(s$age_dir, s$states))
      minlai <- if ("MinLAI" %in% names(necn)) as.numeric(necn$MinLAI) else rep(0.1, nrow(necn))
      pars <- merge(final[, .(LandisName, MaximumLAI, KLAI)], data.table::data.table(LandisName = necn$SpeciesCode, MinLAI = minlai), by = "LandisName")
      xt <- merge(xt[, .(PLT_CN)], necn_site_lai(co, pars, scale)[, .(PLT_CN, x = LAI_necn)], by = "PLT_CN", all.x = TRUE)
      xt[is.na(x), x := 0]
    }
    treeless <- setdiff(plots$PLT_CN, fia$TREE[STATUSCD == 1 & TPA_UNADJ > 0, unique(PLT_CN)])
    xt <- rbind(xt, data.table::data.table(PLT_CN = treeless, x = 0))
    ld <- light_data(xt, seedling_presence(fia, plots$PLT_CN, map),
                     seed_source(fia, plots, map, s$seed_source, s$seed_min_dia, s$seed_radius_km))
    lfit <- fit_light_all(ld, NULL, s$light_min_pres, normalize = s$light_normalize, bare = s$light_bare,
                          c_max = if (scale == "current_code") min(1, s$rescale_f) else 1)
    lfit$settings <- c(lfit$settings, list(scale = scale, seed_how = s$seed_source, seed_dia = s$seed_min_dia, seed_km = s$seed_radius_km, age_dir = s$age_dir))
    lch <- if (!is.null(s$light_choices_file)) data.table::fread(s$light_choices_file) else NULL
    lfinal <- light_final(lfit, lch, necn)
  }
  check <- lai_check_fia(lai, sp_tab, final[, .(LandisName, MaximumLAI, KLAI)])
  say(sprintf("Check vs FIA (stand-level site LAI): %d plots, R2 %.2f, median NECN/observed %.2f", check$stats$plots, check$stats$r2, check$stats$median_ratio))
  out <- build_necn_table(necn, final, s$mode, s$rescale_f, light = lfinal)

  chk <- NULL
  if (!is.null(s$check_run_dir) && nzchar(s$check_run_dir)) {
    say("Checking LANDIS run ", s$check_run_dir, " ...")
    run <- read_landis_run(s$check_run_dir)
    lr <- landis_regen(run, unique(c(0, seq(0, max(run$community_years), by = s$check_every), max(run$community_years))), s$check_regen_age)
    tb <- data.table::fread(run$species_table)
    rpars <- tb[, .(LandisName = SpeciesCode, MaximumLAI = as.numeric(MaximumLAI), KLAI = as.numeric(KLAI), MinLAI = if ("MinLAI" %in% names(tb)) as.numeric(MinLAI) else 0.1)]
    fr <- fia_regen(fia, plots, landscape_outline(run), s$check_buffer_km, map, fit_age_models(fia$TREE[STATUSCD == 1], ref),
                    read_lu_ages(s$age_dir, s$states), s$check_regen_age, sp_tab, rpars, if (is.null(s$check_mode)) s$mode else s$check_mode)
    win <- if (is.null(s$check_window)) c(s$check_regen_age + 1, max(run$community_years)) else unlist(s$check_window)
    chk <- list(run = run, lr = lr, fr = fr, cmp = compare_regen(lr, fr, win), cohorts = cohort_biomass_by_age(lr, fr, win), window = win)
    if (!is.null(n <- fia_age_note(fr))) warning(n, call. = FALSE)
    say(sprintf("  %d forested cells, %d FIA plots within %s km", lr$n_forest, fr$n_plots, s$check_buffer_km))
  }

  dir.create(s$out_dir, showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.Date(), "%Y%m%d")
  data.table::fwrite(out, file.path(s$out_dir, sprintf("NECN_species_%s_%s.csv", s$mode, stamp)))
  data.table::fwrite(final, file.path(s$out_dir, sprintf("LAI_parameters_%s.csv", stamp)))
  if (!is.null(lfinal)) data.table::fwrite(lfinal, file.path(s$out_dir, sprintf("light_curves_%s.csv", stamp)))
  if (!is.null(chk)) data.table::fwrite(chk$cmp, file.path(s$out_dir, sprintf("LANDIS_regeneration_check_%s.csv", stamp)))
  data.table::fwrite(lai$plots, file.path(s$out_dir, sprintf("plot_LAI_%s.csv", stamp)))
  yaml::write_yaml(s, file.path(s$out_dir, sprintf("config_used_%s.yml", stamp)))
  .write_reports(list(settings = s, species_map = map, sla_tab = sla_tab, plot_lai = lai$plots, fit = fit, final = final,
                      mode = s$mode, f = s$rescale_f, lfit = lfit, lfinal = lfinal, ic = ic, check = check, chk = chk), s, stamp, say)
  say("Done: ", normalizePath(s$out_dir))
  invisible(list(fit = fit, final = final, light = lfit, light_final = lfinal, ic = ic, necn = out, plots = lai$plots, species_map = map, sla = sla_tab))
}
