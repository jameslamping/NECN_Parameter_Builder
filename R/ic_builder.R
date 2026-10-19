# ic_builder.R -- NECN initial communities (map + CSV) from TreeMap 2022 and FIA
#
# The ecoregion map is the template: the IC map gets its extent, resolution, CRS and active cells.
#   1. TreeMap (30 m, TM_ID per pixel) is cropped to the template and resampled onto its grid
#      ("near" = nearest 30 m pixel, the default: unbiased area shares; "mode" = majority class, with
#      treeless pixels counted as a class).
#   2. Active cells (ecoregion > 0) with no TreeMap ID (no trees) get the empty-site code, which has no
#      rows in the CSV, so LANDIS starts them empty. Inactive cells get 0.
#   3. TM_ID -> FIA PLT_CN from the raster's attribute table (or the TreeMap tree table). Live trees
#      (STATUSCD 1) of the LANDIS species on that plot become cohorts: species x 10-year age class,
#      leaf = DRYBIO_FOLIAGE, wood = DRYBIO_AG (stem, bark, branches), both x TPA_UNADJ, in g/m2.
#      Ages: Lu et al. (2025) where available, else per-species TOTAGE ~ DIA + HT (see cohorts.R).
#      Ages are capped below species longevity when a LANDIS species file is given.

TONAC_TO_GM2 <- 224.1702   # short tons/acre -> g/m2

#' TM_ID -> PLT_CN lookup from the TreeMap raster attribute table, else from the tree table CSV
treemap_lookup <- function(treemap, table_path = NULL) {
  ct <- terra::cats(treemap)[[1]]
  if (!is.null(ct) && all(c("TM_ID", "PLT_CN") %in% names(ct))) {
    lk <- data.table::as.data.table(ct)[, .(TM_ID = as.integer(TM_ID), PLT_CN = format(as.numeric(PLT_CN), scientific = FALSE, trim = TRUE),
                                            DRYBIO_L = if ("DRYBIO_L" %in% names(ct)) as.numeric(DRYBIO_L) else NA_real_,
                                            FORTYPCD = if ("FORTYPCD" %in% names(ct)) FORTYPCD else NA)]
    return(lk[!duplicated(TM_ID)])
  }
  if (is.null(table_path) || !file.exists(table_path)) stop("The TreeMap raster has no attribute table with TM_ID and PLT_CN; give the TreeMap tree table CSV.")
  lk <- data.table::fread(table_path, select = c("TM_ID", "TM_PLT_CN"), colClasses = list(character = "TM_PLT_CN"))
  unique(lk[, .(TM_ID = as.integer(TM_ID), PLT_CN = TM_PLT_CN, DRYBIO_L = NA_real_, FORTYPCD = NA)], by = "TM_ID")
}

#' TreeMap IDs on the template grid
#' @param template ecoregion SpatRaster (or path); cells > 0 are active
#' @param resample "mode" or "near"
treemap_on_template <- function(treemap_path, template, resample = c("near", "mode"), threads = TRUE) {
  resample <- match.arg(resample)
  tm <- terra::rast(treemap_path)
  eco <- if (is.character(template)) terra::rast(template) else template
  # crop TreeMap to the template footprint (plus 2 cells) in TreeMap's CRS before projecting
  box <- terra::project(terra::as.polygons(terra::ext(eco), crs = terra::crs(eco)), terra::crs(tm))
  box <- terra::buffer(box, 2 * max(terra::res(eco)))
  tmc <- terra::crop(tm, box)
  levels(tmc) <- NULL                                    # resample the IDs, not their labels
  # "no trees" must count as a class, or mode would give a cell a plot whenever any 30 m pixel has one
  tmc <- terra::classify(tmc, cbind(NA, 0))
  out <- terra::project(tmc, eco, method = resample, threads = threads)
  out <- terra::classify(out, cbind(0, NA))
  names(out) <- "TM_ID"
  out
}

#' Read Longevity from a LANDIS core species file (.txt)
read_longevity <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  l <- readLines(path, warn = FALSE)
  l <- trimws(sub(">>.*", "", l))
  l <- l[nzchar(l) & !grepl("^LandisData", l)]
  p <- strsplit(l, "\\s+")
  p <- p[lengths(p) >= 2 & !is.na(suppressWarnings(as.numeric(vapply(p, `[`, "", 2))))]
  data.table::data.table(LandisName = vapply(p, `[`, "", 1), Longevity = as.numeric(vapply(p, `[`, "", 2)))
}

#' Build the IC map and CSV
#' @param fia list from fia_read() (TREE must include the TreeMap plots' states)
#' @param species_map from step 2; models / lu from fit_age_models() / read_lu_ages()
#' @param codes "tm_id" (map code = TreeMap ID) or "renumber" (2..N; crosswalk written alongside)
#' @return list(map, csv, crosswalk, summary, plots_missing, check)
build_ic <- function(template, treemap_path, fia, species_map, models, lu = NULL, treemap_table = NULL,
                     longevity = NULL, resample = "near", codes = c("tm_id", "renumber"), empty_code = 1L,
                     min_biomass = 1, progress = function(msg, value) NULL) {
  codes <- match.arg(codes)
  eco <- if (is.character(template)) terra::rast(template) else template
  progress("Resampling TreeMap onto the ecoregion grid", 0.1)
  tm <- treemap_on_template(treemap_path, eco, resample)
  lk <- treemap_lookup(terra::rast(treemap_path), treemap_table)

  ev <- terra::values(eco, mat = FALSE); tv <- terra::values(tm, mat = FALSE)
  active <- !is.na(ev) & ev > 0
  ids <- sort(unique(tv[active & !is.na(tv)]))

  progress("Cohorts from FIA trees", 0.4)
  pl <- lk[TM_ID %in% ids]
  m <- expand_species_map(species_map)
  tr_all <- fia$TREE[PLT_CN %in% pl$PLT_CN]
  found <- unique(tr_all$PLT_CN)
  tr <- tr_all[STATUSCD == 1 & TPA_UNADJ > 0 & SPCD %in% m$SPCD]
  tr <- merge(assign_tree_ages(tr, models, lu), m, by = "SPCD")
  tr[, age := cohort_class(age)]
  if (!is.null(longevity)) {
    tr <- merge(tr, longevity, by = "LandisName", all.x = TRUE)
    tr[!is.na(Longevity) & age >= Longevity, `:=`(age = pmax(10, floor((Longevity - 1) / 10) * 10), capped = TRUE)]
  }
  co <- tr[, .(LeafBiomass = sum(DRYBIO_FOLIAGE * TPA_UNADJ * LBAC_TO_GM2, na.rm = TRUE),
               WoodBiomass = sum(DRYBIO_AG * TPA_UNADJ * LBAC_TO_GM2, na.rm = TRUE)), by = .(PLT_CN, LandisName, CohortAge = age)]
  co <- merge(co, pl[, .(TM_ID, PLT_CN)], by = "PLT_CN", allow.cartesian = TRUE)
  co[, `:=`(LeafBiomass = round(LeafBiomass), WoodBiomass = round(WoodBiomass))]
  co[, CohortBiomass := LeafBiomass + WoodBiomass]
  co <- co[CohortBiomass >= min_biomass]

  progress("Map codes", 0.75)
  if (codes == "tm_id") {
    if (empty_code %in% ids) {
      new <- max(lk$TM_ID, ids) + 1L
      warning(sprintf("Empty-site code %d is a TreeMap ID on this landscape; using %d instead.", empty_code, new))
      empty_code <- new
    }
    xw <- data.table::data.table(TM_ID = ids, MapCode = ids)
  } else {
    empty_code <- 1L
    xw <- data.table::data.table(TM_ID = ids, MapCode = seq_along(ids) + 1L)
  }
  xw <- merge(xw, pl[, .(TM_ID, PLT_CN, FORTYPCD, DRYBIO_L)], by = "TM_ID", all.x = TRUE)
  xw[, in_fia := PLT_CN %in% found]
  mc <- rep(NA_integer_, length(tv))
  mc[active] <- empty_code
  hit <- active & !is.na(tv)
  mc[hit] <- xw$MapCode[match(tv[hit], xw$TM_ID)]
  mc[!active] <- 0L
  icmap <- terra::rast(eco); terra::values(icmap) <- mc; names(icmap) <- "MapCode"

  co <- merge(co, xw[, .(TM_ID, MapCode)], by = "TM_ID")
  csv <- co[order(MapCode, LandisName, CohortAge), .(MapCode, SpeciesName = LandisName, CohortAge, CohortBiomass, LeafBiomass, WoodBiomass)]

  progress("Summaries", 0.9)
  cells <- data.table::data.table(MapCode = mc[active])[, .N, by = MapCode]
  xw <- merge(xw, cells[, .(MapCode, cells = N)], by = "MapCode", all.x = TRUE)
  tot <- csv[, .(wood_gm2 = sum(WoodBiomass), leaf_gm2 = sum(LeafBiomass), n_cohorts = .N, n_species = data.table::uniqueN(SpeciesName)), by = MapCode]
  xw <- merge(xw, tot, by = "MapCode", all.x = TRUE)
  xw[, treemap_live_gm2 := DRYBIO_L * TONAC_TO_GM2]
  sp <- merge(csv, cells, by = "MapCode")[, .(cells = sum(N), landscape_Mg = sum((CohortBiomass) * N) * prod(terra::res(eco)) / 1e6,
                                               cohorts = .N, mean_age = stats::weighted.mean(CohortAge, CohortBiomass), max_age = max(CohortAge)),
                                           by = SpeciesName][order(-landscape_Mg)]
  n_act <- sum(active); n_trees <- sum(hit)
  lost <- xw[in_fia == FALSE, sum(cells, na.rm = TRUE)]
  summary <- data.table::data.table(
    item = c("Active cells (ecoregion > 0)", "Active cells with a TreeMap ID", "Active cells with no TreeMap ID (empty-site code)",
             "TreeMap IDs (plots) on the landscape", "Plots found in the loaded FIA TREE tables", "Cells whose plot is missing from FIA (start empty)",
             "Map codes with at least one cohort", "Cohorts in the CSV", "Empty-site map code",
             "Trees with ages capped at longevity"),
    value = c(n_act, n_trees, n_act - n_trees, length(ids), sum(xw$in_fia), lost, data.table::uniqueN(csv$MapCode), nrow(csv), empty_code,
              if ("capped" %in% names(tr)) sum(tr$capped, na.rm = TRUE) else NA))
  progress("done", 1)
  list(map = icmap, treemap = tm, csv = csv, crosswalk = xw[order(MapCode)], species = sp, summary = summary,
       age_sources = tr[, .N, by = age_src], empty_code = empty_code)
}

#' Write the IC map (INT4S GeoTIFF) and CSV
write_ic <- function(ic, dir, stem = paste0("IC_", format(Sys.Date(), "%Y%m%d"))) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- c(map = file.path(dir, paste0(stem, "_map.tif")), csv = file.path(dir, paste0(stem, ".csv")),
         crosswalk = file.path(dir, paste0(stem, "_crosswalk.csv")))
  terra::writeRaster(ic$map, f[["map"]], datatype = "INT4S", overwrite = TRUE, NAflag = -9999)
  data.table::fwrite(ic$csv, f[["csv"]])
  data.table::fwrite(ic$crosswalk, f[["crosswalk"]])
  invisible(f)
}
