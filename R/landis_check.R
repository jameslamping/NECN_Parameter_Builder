# landis_check.R -- compare regeneration in a LANDIS-II NECN run with FIA plots around the landscape
#
# Regeneration = cohorts (LANDIS) or live trees (FIA) up to `regen_age` years old, by species.
#   LANDIS: community-input-file-{t}.csv (Output Biomass Community) + output-community-0.tif (site map codes).
#           Only cells forested at year 0 count, and only years after the initial communities' young
#           cohorts have aged out (t > regen_age) are compared, so all regeneration counted was simulated.
#   FIA:    live trees on the selected plots within `buffer_km` of the landscape, aged with Lu et al. (2025)
#           or the TOTAGE ~ DIA + HT models (as the initial communities are built), plus SEEDLING presence.
# Main metric: young biomass density (leaf + wood, Mg/ha), which does not depend on plot or cell size.
# Presence frequency is reported too, but a 90 m cell (8,100 m2) samples far more ground than FIA's
# microplots (54 m2) or subplots (672 m2), so cell frequencies run higher for any species.
# Site LAI at regenerating cells / plots checks the light curves: LANDIS LAI-{t}.tif maps (whole numbers)
# against the site LAI the run's species table gives on the FIA plots.

#' Paths and settings of a LANDIS run folder
read_landis_run <- function(run_dir) {
  if (!dir.exists(run_dir)) stop("Folder not found: ", run_dir)
  sc <- file.path(run_dir, "Scenario.txt")
  if (!file.exists(sc)) stop("No Scenario.txt in ", run_dir)
  rd <- function(f) { l <- readLines(f, warn = FALSE); l <- sub(">>.*", "", l); trimws(l[nzchar(trimws(l))]) }
  kv <- function(lines, key) { m <- grep(sprintf("^%s\\s", key), lines, value = TRUE); if (length(m)) trimws(sub(sprintf("^%s\\s+", key), "", m[1])) else NA_character_ }
  s <- rd(sc)
  necn_line <- grep('^"NECN Succession"', s, value = TRUE)
  necn_file <- if (length(necn_line)) file.path(run_dir, trimws(sub('^"NECN Succession"\\s+', "", necn_line[1]))) else NA_character_
  n <- if (!is.na(necn_file) && file.exists(necn_file)) rd(necn_file) else character()
  yrs <- as.integer(sub("community-input-file-(\\d+)\\.csv", "\\1", list.files(run_dir, "^community-input-file-\\d+\\.csv$")))
  lai_yrs <- as.integer(sub("LAI-(\\d+)\\.tif", "\\1", list.files(file.path(run_dir, "NECN"), "^LAI-\\d+\\.tif$")))
  sp_tab <- kv(n, "SpeciesParameters")
  list(dir = run_dir, duration = as.integer(kv(s, "Duration")), ecoregion_map = file.path(run_dir, kv(s, "EcoregionsMap")),
       necn_file = necn_file, species_table = if (is.na(sp_tab)) NA_character_ else file.path(run_dir, sp_tab),
       community_years = sort(yrs), lai_years = sort(lai_yrs),
       code_map = file.path(run_dir, "output-community-0.tif"), reproduction_log = file.path(run_dir, "NECN-reproduction-log.csv"))
}

#' Young-cohort summaries from a LANDIS run
#' @return list(regen = year x species, lai_young = year x species x LAI, lai_all = year x LAI, n_forest, years)
landis_regen <- function(run, years, regen_age = 20, progress = function(msg, value) NULL) {
  years <- sort(intersect(years, run$community_years))
  if (!length(years)) stop("None of the chosen years have a community-input-file in the run folder.")
  codes <- terra::values(terra::rast(run$code_map), mat = FALSE)          # site map code per cell (LANDIS output orientation)
  cells <- data.table::data.table(cell = seq_along(codes), MapCode = codes)[!is.na(MapCode) & MapCode > 0]
  rd <- function(t) data.table::fread(file.path(run$dir, sprintf("community-input-file-%d.csv", t)),
                                      select = c("MapCode", "SpeciesName", "CohortAge", "CohortBiomass"), fill = TRUE, showProgress = FALSE)
  progress("Year 0 (forested cells)", 0.02)
  forest <- unique(rd(0)$MapCode)
  n_forest <- length(forest)
  out <- list(); lai_y <- list(); lai_a <- list()
  for (i in seq_along(years)) {
    t <- years[i]; progress(sprintf("Year %d", t), 0.05 + 0.9 * i / length(years))
    d <- rd(t)[CohortAge <= regen_age & MapCode %in% forest]
    out[[i]] <- d[, .(cells_young = data.table::uniqueN(MapCode), young_gm2_sum = sum(as.numeric(CohortBiomass))), by = SpeciesName][, year := t]
    lf <- file.path(run$dir, "NECN", sprintf("LAI-%d.tif", t))
    if (file.exists(lf)) {
      lv <- terra::values(terra::rast(lf), mat = FALSE)
      lc <- cells[, .(MapCode, lai = lv[cell])][MapCode %in% forest]
      lai_a[[i]] <- lc[, .N, by = lai][, year := t]
      yl <- merge(unique(d[, .(MapCode, SpeciesName)]), lc, by = "MapCode")
      lai_y[[i]] <- yl[, .N, by = .(SpeciesName, lai)][, year := t]
    }
  }
  regen <- data.table::rbindlist(out)
  regen[, `:=`(young_Mgha = young_gm2_sum / n_forest / 100, pct_cells = 100 * cells_young / n_forest)]
  list(regen = regen[], lai_young = data.table::rbindlist(lai_y), lai_all = data.table::rbindlist(lai_a),
       n_forest = n_forest, years = years, regen_age = regen_age)
}

#' Landscape outline (active ecoregion cells) as an sf polygon in EPSG:5070
landscape_outline <- function(run) {
  e <- terra::rast(run$ecoregion_map)
  if (!nzchar(terra::crs(e))) stop("The ecoregion map has no coordinate reference system, so FIA plots around it cannot be located.")
  f <- max(1, floor(max(dim(e)[1:2]) / 300))
  a <- if (f > 1) terra::aggregate(e > 0, f, fun = "max", na.rm = TRUE) else e > 0
  a <- terra::ifel(a > 0, 1, NA)
  p <- terra::as.polygons(a, dissolve = TRUE)
  sf::st_transform(sf::st_make_valid(sf::st_as_sf(p)), 5070)
}

#' FIA plots near the landscape and their young-tree biomass, seedling presence and site LAI by species
#' @param pars data.table(LandisName, MaximumLAI, KLAI, MinLAI) from the run's species table
#' @param lai_mode "stand_level" or "current_code" (how the run's NECN computes site LAI)
fia_regen <- function(fia, plots, outline, buffer_km, species_map, models, lu = NULL, regen_age = 20,
                      sp_tab = NULL, pars = NULL, lai_mode = "stand_level", progress = function(msg, value) NULL) {
  progress("FIA plots near the landscape", 0.1)
  pts <- sf::st_transform(sf::st_as_sf(plots[!is.na(LAT) & !is.na(LON)], coords = c("LON", "LAT"), crs = 4269, remove = FALSE), 5070)
  reg <- sf::st_union(if (buffer_km > 0) sf::st_buffer(outline, buffer_km * 1000) else outline)
  near <- pts$PLT_CN[lengths(sf::st_intersects(pts, reg)) > 0]
  if (!length(near)) stop("No selected FIA plots within ", buffer_km, " km of the landscape. Load the right states in step 1 or widen the distance.")
  m <- expand_species_map(species_map)
  progress("Tree ages", 0.3)
  tr <- fia$TREE[PLT_CN %in% near & STATUSCD == 1 & TPA_UNADJ > 0 & SPCD %in% m$SPCD]
  tr <- merge(assign_tree_ages(tr, models, lu), m, by = "SPCD")
  young <- tr[age <= regen_age, .(young_gm2 = sum((DRYBIO_AG + DRYBIO_FOLIAGE) * TPA_UNADJ * LBAC_TO_GM2, na.rm = TRUE)), by = .(PLT_CN, LandisName)]
  grid <- data.table::CJ(PLT_CN = near, LandisName = unique(m$LandisName))
  d <- merge(grid, young, by = c("PLT_CN", "LandisName"), all.x = TRUE)[is.na(young_gm2), young_gm2 := 0]
  progress("Seedlings", 0.6)
  if (!is.null(fia$SEEDLING)) {
    sd <- seedling_presence(fia, near, species_map)[, seedlings := TRUE]
    d <- merge(d, sd[, .(PLT_CN, LandisName, seedlings)], by = c("PLT_CN", "LandisName"), all.x = TRUE)[is.na(seedlings), seedlings := FALSE]
  } else d[, seedlings := NA]
  lai <- NULL
  if (!is.null(pars)) {
    progress("Site LAI with the run's MaximumLAI / KLAI", 0.8)
    lai <- if (lai_mode == "stand_level" && !is.null(sp_tab)) stand_site_lai(sp_tab[PLT_CN %in% near], pars)
           else necn_site_lai(plot_cohorts(fia, near, species_map, models, lu), pars, "current_code")[, .(PLT_CN, LAI_necn)]
    lai <- merge(data.table::data.table(PLT_CN = near), lai, by = "PLT_CN", all.x = TRUE)[is.na(LAI_necn), LAI_necn := 0]
  }
  list(plots = d[], n_plots = length(near), plot_lai = lai, buffer_km = buffer_km, regen_age = regen_age)
}

#' Species table: LANDIS (window of years) against FIA, with a bootstrap 90% interval for the FIA mean
compare_regen <- function(lr, fr, window, nboot = 500, seed = 1) {
  set.seed(seed)
  w <- lr$regen[year >= window[1] & year <= window[2]]
  ls <- w[, .(landis_Mgha = mean(young_Mgha), landis_min = min(young_Mgha), landis_max = max(young_Mgha), landis_pct_cells = mean(pct_cells)), by = .(LandisName = SpeciesName)]
  n <- fr$n_plots
  fs <- fr$plots[, {
    v <- young_gm2 / 100
    b <- if (nboot > 0) replicate(nboot, mean(v[sample.int(length(v), replace = TRUE)])) else NA_real_
    .(fia_Mgha = mean(v), fia_lo = stats::quantile(b, .05, na.rm = TRUE), fia_hi = stats::quantile(b, .95, na.rm = TRUE),
      fia_pct_plots = 100 * mean(young_gm2 > 0), fia_pct_seedlings = if (all(is.na(seedlings))) NA_real_ else 100 * mean(seedlings))
  }, by = LandisName]
  x <- merge(fs, ls, by = "LandisName", all = TRUE)
  for (cc in c("landis_Mgha", "landis_min", "landis_max", "landis_pct_cells")) data.table::set(x, which(is.na(x[[cc]])), cc, 0)
  x[, ratio := landis_Mgha / fia_Mgha]
  # density = how often x how much where present: separates establishment from young-cohort growth
  x[, `:=`(fia_Mgha_where_present = data.table::fifelse(fia_pct_plots > 0, fia_Mgha / (fia_pct_plots / 100), NA_real_),
           landis_Mgha_where_present = data.table::fifelse(landis_pct_cells > 0, landis_Mgha / (landis_pct_cells / 100), NA_real_))]
  x[, assessment := data.table::fcase(
    fia_Mgha < 0.01 & landis_Mgha < 0.01, "little regeneration in either",
    fia_Mgha < 0.01, "LANDIS regenerates; FIA has almost none",
    landis_Mgha < fia_lo / 2, "too little in LANDIS",
    landis_Mgha > fia_hi * 2, "too much in LANDIS",
    landis_Mgha < fia_lo, "somewhat low",
    landis_Mgha > fia_hi, "somewhat high",
    default = "within the FIA range")]
  x[, `:=`(fia_share = fia_Mgha / sum(fia_Mgha), landis_share = landis_Mgha / sum(landis_Mgha))]
  x[order(-pmax(fia_Mgha, landis_Mgha))]
}

#' New cohorts per year from NECN-reproduction-log.csv, per 1,000 forested cells
read_reproduction_log <- function(run, n_forest) {
  if (!file.exists(run$reproduction_log)) return(NULL)
  d <- data.table::fread(run$reproduction_log, fill = TRUE)
  d[, SpeciesName := trimws(SpeciesName)]
  keep <- intersect(c("NumCohortsPlanting", "NumCohortsSerotiny", "NumCohortsResprout", "NumCohortsSeed"), names(d))
  d <- data.table::melt(d[, c("Time", "SpeciesName", keep), with = FALSE], id.vars = c("Time", "SpeciesName"), variable.name = "type", value.name = "n")
  d[, type := sub("NumCohorts", "", type)][, per_1000_cells := 1000 * n / n_forest]
  d[n > 0]
}
