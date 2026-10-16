# plots.R -- choosing FIA plots for calibration

#' Select plots: sampled, annual inventory, mostly forested, one visit per plot location,
#' optionally inside a study-area polygon (plus buffer)
#' @param region optional sf polygon(s); plots are kept if their (fuzzed) location falls inside
#' @param buffer_km buffer applied to `region`
select_plots <- function(fia, forest_min = 0.9, annual_only = TRUE, min_year = 2000,
                         exclude_disturbed = FALSE, region = NULL, buffer_km = 0) {
  p <- data.table::copy(fia$PLOT)
  p <- p[PLOT_STATUS_CD == 1 & INVYR >= min_year & INVYR < 9999]
  if (annual_only) p <- p[KINDCD %in% 1:3]

  cnd <- fia$COND[COND_STATUS_CD == 1, .(forest_prop = sum(CONDPROP_UNADJ, na.rm = TRUE),
                                         STDAGE = suppressWarnings(max(as.numeric(STDAGE), na.rm = TRUE)),
                                         disturbed = any(!is.na(DSTRBCD1) & DSTRBCD1 > 0) |
                                                     any(!is.na(TRTCD1) & TRTCD1 > 0)), by = PLT_CN]
  p <- merge(p, cnd, by.x = "CN", by.y = "PLT_CN")
  p <- p[forest_prop >= forest_min]
  if (exclude_disturbed) p <- p[disturbed == FALSE]
  p[!is.finite(STDAGE), STDAGE := NA]

  if (!is.null(region)) {
    pts <- sf::st_as_sf(p[!is.na(LAT) & !is.na(LON)], coords = c("LON", "LAT"), crs = 4269, remove = FALSE)
    reg <- sf::st_transform(sf::st_make_valid(region), 5070)
    if (buffer_km > 0) reg <- sf::st_buffer(reg, buffer_km * 1000)
    inside <- lengths(sf::st_intersects(sf::st_transform(pts, 5070), sf::st_union(reg))) > 0
    p <- p[CN %in% pts$CN[inside]]
  }

  # one visit per plot location: the most recent
  p[, loc := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "_")]
  p <- p[order(-INVYR)][!duplicated(loc)]
  data.table::setnames(p, "CN", "PLT_CN")
  p[]
}
