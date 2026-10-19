# figures.R -- plots used in the app and the report

app_dir <- function() getOption("necnbuilder.dir", ".")

COL <- c(mono = "#2a78d6", mixed = "#1baf7a", pooled = "#4a3aa7", current = "#eb6834")

theme_builder <- function(base = 10) ggplot2::theme_classic(base_size = base) +
  ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"))

#' LAI vs stand wood per species: monospecific plots (grey), mixed plots as pure-stand equivalents
#' (green), and the three fitted curves; vertical lines = chosen and current KLAI
fig_species_fit <- function(fit, final = NULL, species = NULL, ncol = 4, base = 9) {
  p <- fit$params
  sp <- if (!is.null(species)) species else p[!is.na(KLAI_mono) | !is.na(KLAI_mixed) | !is.na(KLAI_pooled), LandisName]
  if (!length(sp)) return(ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = "No species with enough plots") + ggplot2::theme_void())
  mono <- fit$mono[LandisName %in% sp]; mix <- fit$mixed[LandisName %in% sp]
  xmax <- rbind(mono[, .(LandisName, C)], mix[, .(LandisName, C)])[, .(xmax = if (.N >= 5) stats::quantile(C / C_FRAC, .98) else 50000), by = LandisName]
  xmax <- merge(data.table::data.table(LandisName = sp), xmax, by = "LandisName", all.x = TRUE)[is.na(xmax), xmax := 50000]
  mono <- merge(mono, xmax, by = "LandisName")[C / C_FRAC <= xmax]
  mix <- merge(mix, xmax, by = "LandisName")[C / C_FRAC <= xmax]
  pp <- merge(p[LandisName %in% sp], xmax, by = "LandisName")
  cv <- pp[, data.table::rbindlist(list(
    if (!is.na(KLAI_mono)) curve_xy(MaxLAI_mono, KLAI_mono, xmax)[, est := "mono"],
    if (!is.na(KLAI_mixed)) curve_xy(MaxLAI_mixed, KLAI_mixed, xmax)[, est := "mixed"],
    if (!is.na(KLAI_pooled)) curve_xy(MaxLAI_pooled, KLAI_pooled, xmax)[, est := "pooled"])), by = LandisName]
  lab <- pp[, .(LandisName, lab = sprintf("mono n=%d | mixed n=%d", n_mono, n_mixed))]
  one <- length(sp) == 1
  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = mix, ggplot2::aes(wood_gm2 / 100, y), colour = COL[["mixed"]], alpha = if (one) 0.35 else 0.12, size = if (one) 1.3 else 0.5) +
    ggplot2::geom_point(data = mono, ggplot2::aes(wood_gm2 / 100, y), colour = "grey25", alpha = if (one) 0.55 else 0.25, size = if (one) 1.6 else 0.6) +
    ggplot2::geom_line(data = cv, ggplot2::aes(wood_Mgha, LAI, colour = est, linetype = est), linewidth = if (one) 1.2 else 0.8) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.4, size = base / 3.4, colour = "grey30")
  if (!is.null(final)) {
    vl <- final[LandisName %in% sp, .(LandisName, chosen = KLAI / C_FRAC / 100, current = KLAI_current / C_FRAC / 100)]
    vl <- merge(vl, xmax, by = "LandisName")
    g <- g + ggplot2::geom_vline(data = vl[chosen <= xmax / 100], ggplot2::aes(xintercept = chosen), colour = "black", linewidth = 0.6) +
      ggplot2::geom_vline(data = vl[current <= xmax / 100], ggplot2::aes(xintercept = current), colour = COL[["current"]], linetype = "22", linewidth = 0.6)
  }
  g + ggplot2::scale_colour_manual(NULL, values = COL[c("mono", "mixed", "pooled")],
                                   labels = c(mono = "monospecific fit", mixed = "mixed-plot fit", pooled = "group-pooled")) +
    ggplot2::scale_linetype_manual(NULL, values = c(mono = "solid", mixed = "22", pooled = "11"),
                                   labels = c(mono = "monospecific fit", mixed = "mixed-plot fit", pooled = "group-pooled")) +
    ggplot2::facet_wrap(~LandisName, ncol = ncol, scales = "free") + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "Live wood in the stand (Mg/ha)", y = "LAI",
                  caption = "Points: grey = monospecific plots (plot LAI); green = mixed plots (species LAI / wood share).\nBlack vertical = KLAI used; orange dashed = current KLAI.") +
    theme_builder(base) + ggplot2::theme(legend.position = "top", axis.title.y = ggplot2::element_text(size = base - 1))
}

#' LAI vs stand age for monospecific plots (intuition only)
fig_species_age <- function(fit, species = NULL, ncol = 4, base = 9) {
  d <- fit$mono[!is.na(STDAGE) & STDAGE > 0 & STDAGE <= 500]
  sp <- if (!is.null(species)) species else fit$params[n_mono >= 12, LandisName]
  d <- d[LandisName %in% sp]
  big <- fit$params[n_mono >= 50, LandisName]
  ggplot2::ggplot(d, ggplot2::aes(STDAGE, y)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.6, colour = "grey45") +
    ggplot2::geom_smooth(data = d[LandisName %in% big], method = "loess", span = 0.9, se = FALSE,
                         colour = COL[["pooled"]], linewidth = 0.8, formula = y ~ x) +
    ggplot2::facet_wrap(~LandisName, ncol = ncol, scales = "free") + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "FIA stand age (years)", y = "Plot LAI") + theme_builder(base)
}

#' All selected plots: LAI vs wood (2D bins)
fig_all_plots <- function(plot_lai) {
  ggplot2::ggplot(plot_lai, ggplot2::aes(wood_gm2 / 100, LAI)) +
    ggplot2::geom_bin2d(bins = 60) + ggplot2::scale_fill_viridis_c(trans = "log10", name = "Plots") +
    ggplot2::coord_cartesian(xlim = c(0, stats::quantile(plot_lai$wood_gm2 / 100, .99))) +
    ggplot2::labs(x = "Live wood (Mg/ha)", y = "Plot LAI") + theme_builder(11)
}

fig_plot_map <- function(plot_lai) {
  ggplot2::ggplot(plot_lai, ggplot2::aes(LON, LAT, colour = LAI)) + ggplot2::geom_point(size = 0.4, alpha = 0.6) +
    ggplot2::scale_colour_viridis_c(limits = c(0, stats::quantile(plot_lai$LAI, .98)), oob = scales::squish) +
    ggplot2::coord_quickmap() + ggplot2::labs(x = NULL, y = NULL, colour = "LAI") + theme_builder(11)
}

LIGHT_SCALE_LABEL <- c(current_code = "Site LAI as NECN 8.3 and earlier compute it (sum of cohort LAI, before rescale)",
                       stand_level = "Stand-level site LAI (NECN with the stand-level change)",
                       observed = "Observed plot LAI (foliage x SLA)")

#' Seedling presence vs site LAI per species: binned presence (points, 90% intervals), the fitted
#' presence curve (blue), the curve exported to NECN (black dashed) and the current table curve (orange)
fig_light_fit <- function(lfit, lfinal, species = NULL, scale = "current_code", ncol = 4, base = 9) {
  sp <- if (!is.null(species)) species else lfit$params[n_plots > 0, LandisName]
  b <- lfit$bins[LandisName %in% sp]
  xmax <- lfit$data[LandisName %in% sp, .(xmax = max(1, stats::quantile(x, .99))), by = LandisName]
  grid <- xmax[, .(x = seq(0.02, xmax, length.out = 200)), by = LandisName]
  p <- merge(lfit$params[LandisName %in% sp, .(LandisName, Shape, Scale, Location, Adjust_fit)],
             lfinal[, .(LandisName, Shape_used, Scale_used, Location_used, Adjust_used, LightLAIShape, LightLAIScale, LightLAILocation, LightLAIAdjust)],
             by = "LandisName", all.x = TRUE)
  cv <- merge(grid, p, by = "LandisName")
  cv <- rbind(
    cv[!is.na(Shape), .(LandisName, x, p = light_prob(x, Shape, Scale, Location, Adjust_fit), curve = "fit")],
    cv[!is.na(Shape_used), .(LandisName, x, p = light_prob(x, Shape_used, Scale_used, Location_used, Adjust_used), curve = "used")],
    cv[!is.na(LightLAIShape), .(LandisName, x, p = light_prob(x, LightLAIShape, LightLAIScale, LightLAILocation, LightLAIAdjust), curve = "current")])
  lab <- lfit$params[LandisName %in% sp, .(LandisName, lab = sprintf("%d plots, %d with seedlings", n_plots, n_present))]
  one <- length(sp) == 1
  labs_c <- c(fit = "fitted presence (FIA)", used = "exported light probability", current = "current table")
  ggplot2::ggplot() +
    ggplot2::geom_linerange(data = b, ggplot2::aes(x, ymin = lo, ymax = hi), colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_point(data = b, ggplot2::aes(x, p, size = n), colour = "grey25", alpha = 0.8) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x, p, colour = curve, linetype = curve), linewidth = if (one) 1.1 else 0.7) +
    { if ("p0_used" %in% names(lfinal)) ggplot2::geom_point(data = lfinal[LandisName %in% sp & !is.na(p0_used), .(LandisName, x = 0, p = p0_used)],
                                                             ggplot2::aes(x, p), shape = 23, fill = "white", colour = "black", size = if (one) 3 else 1.8) } +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.4, size = base / 3.4, colour = "grey30") +
    ggplot2::scale_colour_manual(NULL, values = c(fit = COL[["mono"]], used = "black", current = COL[["current"]]), labels = labs_c) +
    ggplot2::scale_linetype_manual(NULL, values = c(fit = "solid", used = "22", current = "solid"), labels = labs_c) +
    ggplot2::scale_size_area(max_size = if (one) 4 else 2.2, guide = "none") +
    ggplot2::facet_wrap(~LandisName, ncol = ncol, scales = "free_x") + ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(x = LIGHT_SCALE_LABEL[[scale]], y = "Probability",
                  caption = "Points: share of plots (with a seed source) that have seedlings of the species, in equal-count LAI bins; bars = 90% intervals.\nDiamond at LAI 0: exported light probability on bare ground (a cell with no cohorts).") +
    theme_builder(base) + ggplot2::theme(legend.position = "top")
}

#' Observed LAI vs the LAI NECN computes on the same plots
fig_lai_scales <- function(xt) {
  keep <- intersect(c("stand", "coded"), names(xt))
  lab <- c(stand = "NECN with stand-level site LAI", coded = "NECN 8.3 and earlier (cohort sum)")
  d <- data.table::melt(xt[, c("PLT_CN", "obs", keep), with = FALSE], id.vars = c("PLT_CN", "obs"), variable.name = "scale", value.name = "necn")
  d[, scale := factor(lab[as.character(scale)], levels = lab)]
  lim <- stats::quantile(c(d$obs, d$necn), .99, na.rm = TRUE)
  ggplot2::ggplot(d[!is.na(necn)], ggplot2::aes(obs, necn)) + ggplot2::geom_bin2d(bins = 60) +
    ggplot2::geom_abline(linetype = "22") + ggplot2::scale_fill_viridis_c(trans = "log10", name = "Plots") +
    ggplot2::facet_wrap(~scale) + ggplot2::coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
    ggplot2::labs(x = "Observed plot LAI (foliage x SLA)", y = "Site LAI NECN computes from the plot's cohorts") + theme_builder(11)
}

#' IC map shown as live aboveground biomass per cell (Mg/ha); grey = empty site, white = inactive
fig_ic_biomass <- function(ic, max_cols = 700) {
  tot <- ic$csv[, .(bio = sum(CohortBiomass) / 100), by = MapCode]   # g/m2 -> Mg/ha
  r <- ic$map
  f <- max(1, ceiling(terra::ncol(r) / max_cols))
  if (f > 1) r <- terra::aggregate(r, f, fun = "modal", na.rm = TRUE)
  b <- terra::subst(r, tot$MapCode, tot$bio, others = NA)
  empty <- terra::ifel(r > 0 & is.na(b), 1, NA)
  d <- data.table::as.data.table(terra::as.data.frame(b, xy = TRUE, na.rm = TRUE)); data.table::setnames(d, 3, "bio")
  e <- data.table::as.data.table(terra::as.data.frame(empty, xy = TRUE, na.rm = TRUE))
  ggplot2::ggplot() + ggplot2::geom_raster(data = e, ggplot2::aes(x, y), fill = "grey80") +
    ggplot2::geom_raster(data = d, ggplot2::aes(x, y, fill = bio)) +
    ggplot2::scale_fill_viridis_c("Live biomass\n(Mg/ha)", limits = c(0, stats::quantile(d$bio, .98)), oob = scales::squish) +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL, caption = "Grey: active cells with no trees (empty-site code).") + theme_builder(11)
}

#' IC live biomass per map code against TreeMap's own live biomass (DRYBIO_L) for the same plot
fig_ic_check <- function(xw) {
  d <- xw[!is.na(treemap_live_gm2) & !is.na(wood_gm2) & cells > 0]
  lim <- stats::quantile(c(d$treemap_live_gm2, d$wood_gm2 + d$leaf_gm2) / 100, .99, na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(treemap_live_gm2 / 100, (wood_gm2 + leaf_gm2) / 100, size = cells)) +
    ggplot2::geom_point(alpha = 0.35, colour = "#1f3b57") + ggplot2::geom_abline(linetype = "22") +
    ggplot2::scale_size_area(max_size = 5, name = "Cells") + ggplot2::coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
    ggplot2::labs(x = "TreeMap DRYBIO_L, all live trees (Mg/ha)", y = "IC: leaf + wood, LANDIS species (Mg/ha)",
                  caption = "One point per map code (FIA plot). Below the line: trees of species not in the LANDIS table, or foliage/branch definitions differ.") + theme_builder(11)
}

#' Observed plot LAI against the site LAI stand-level NECN computes with the final MaximumLAI / KLAI, by stand wood
fig_lai_check <- function(chk) {
  b <- chk$bins[plots >= 10]
  st <- chk$stats
  ggplot2::ggplot(b, ggplot2::aes(wood)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = obs_q25, ymax = obs_q75), fill = "grey88") +
    ggplot2::geom_line(ggplot2::aes(y = observed, colour = "FIA observed (median, 25th-75th percentile)"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = necn, colour = "NECN stand-level site LAI with these parameters (median)"), linewidth = 1, linetype = "22") +
    ggplot2::geom_point(ggplot2::aes(y = observed, colour = "FIA observed (median, 25th-75th percentile)"), size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(y = necn, colour = "NECN stand-level site LAI with these parameters (median)"), size = 1.8) +
    ggplot2::scale_colour_manual(NULL, values = c("grey20", COL[["mono"]])) + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "Live wood in the stand (Mg/ha; median of each class)", y = "Site LAI",
                  caption = sprintf("%s FIA plots where the LANDIS species hold >= 90%% of live wood. R2 = %.2f, median NECN/observed = %.2f.",
                                    format(st$plots, big.mark = ","), st$r2, st$median_ratio)) +
    theme_builder(11) + ggplot2::theme(legend.position = "top", legend.direction = "vertical")
}

# ---- LANDIS run check -----------------------------------------------------------------------------------------
REGEN_COL <- c(FIA = "grey15", LANDIS = "#2a78d6")

#' Young-cohort regeneration by species: density (scale-free), how often, and how much where present
fig_regen_compare <- function(cmp, top = Inf) {
  x <- cmp[seq_len(min(.N, top))]
  lev <- rev(x$LandisName)
  f <- function(v) pmax(v, 0.001)
  d <- rbind(
    x[, .(LandisName, panel = "a. Young biomass (Mg/ha)", src = "FIA", v = f(fia_Mgha), lo = f(fia_lo), hi = f(fia_hi))],
    x[, .(LandisName, panel = "a. Young biomass (Mg/ha)", src = "LANDIS", v = f(landis_Mgha), lo = f(landis_min), hi = f(landis_max))],
    x[, .(LandisName, panel = "b. Plots / cells with young trees (%)", src = "FIA", v = fia_pct_plots, lo = NA_real_, hi = NA_real_)],
    x[, .(LandisName, panel = "b. Plots / cells with young trees (%)", src = "LANDIS", v = landis_pct_cells, lo = NA_real_, hi = NA_real_)],
    if (!all(is.na(x$fia_pct_seedlings))) x[, .(LandisName, panel = "b. Plots / cells with young trees (%)", src = "FIA seedlings", v = fia_pct_seedlings, lo = NA_real_, hi = NA_real_)],
    x[, .(LandisName, panel = "c. Young biomass where present (Mg/ha)", src = "FIA", v = f(fia_Mgha_where_present), lo = NA_real_, hi = NA_real_)],
    x[, .(LandisName, panel = "c. Young biomass where present (Mg/ha)", src = "LANDIS", v = f(landis_Mgha_where_present), lo = NA_real_, hi = NA_real_)])
  d[, LandisName := factor(LandisName, levels = lev)]
  d[, y := as.numeric(LandisName) + data.table::fifelse(src == "LANDIS", -0.15, 0.15)]
  ggplot2::ggplot(d, ggplot2::aes(v, y, colour = src)) +
    ggplot2::geom_segment(data = d[!is.na(lo)], ggplot2::aes(x = lo, xend = hi, yend = y), linewidth = 0.8, alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = src), size = 2.2) +
    ggplot2::scale_colour_manual(NULL, values = c(REGEN_COL, `FIA seedlings` = REGEN_COL[["FIA"]]),
                                 labels = c(FIA = "FIA plots (90% interval)", `FIA seedlings` = "FIA plots with seedlings", LANDIS = "LANDIS cells (range across years)")) +
    ggplot2::scale_shape_manual(NULL, values = c(FIA = 16, `FIA seedlings` = 1, LANDIS = 16),
                                labels = c(FIA = "FIA plots (90% interval)", `FIA seedlings` = "FIA plots with seedlings", LANDIS = "LANDIS cells (range across years)")) +
    ggplot2::scale_y_continuous(breaks = seq_along(lev), labels = lev) +
    ggplot2::facet_wrap(~panel, nrow = 1, scales = "free_x") +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.01, base = 10), labels = function(b) vapply(b, format, "", scientific = FALSE, drop0trailing = TRUE),
                                breaks = function(l) { b <- c(0, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300); b[b >= l[1] & b <= l[2]] }) +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = "a = b x c (filled points). (a) does not depend on plot or cell size and is the main comparison. (b) a 90 m cell samples far more ground than FIA microplots,\nso cell frequencies run higher for any species; open circles = FIA plots with seedlings (< 1 in DBH, no biomass in FIA).\nA gap in (c) with similar (b) points to young-cohort growth (MaximumANPP, KLAI) rather than establishment.") +
    theme_builder(10) + ggplot2::theme(legend.position = "top", panel.grid.major.y = ggplot2::element_line(colour = "grey92"))
}

#' LANDIS young biomass through time with the FIA mean and interval
fig_regen_time <- function(lr, cmp, top = 9) {
  sp <- cmp[seq_len(min(.N, top)), LandisName]
  d <- lr$regen[SpeciesName %in% sp][, LandisName := factor(SpeciesName, levels = sp)]
  f <- cmp[LandisName %in% sp][, LandisName := factor(LandisName, levels = sp)]
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = f, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = fia_lo, ymax = fia_hi), fill = "grey88") +
    ggplot2::geom_hline(data = f, ggplot2::aes(yintercept = fia_Mgha), colour = REGEN_COL[["FIA"]], linetype = "22") +
    ggplot2::geom_vline(xintercept = lr$regen_age, colour = "grey55", linetype = "13") +
    ggplot2::geom_line(data = d, ggplot2::aes(year, young_Mgha), colour = REGEN_COL[["LANDIS"]], linewidth = 0.9) +
    ggplot2::geom_point(data = d, ggplot2::aes(year, young_Mgha), colour = REGEN_COL[["LANDIS"]], size = 1) +
    ggplot2::facet_wrap(~LandisName, scales = "free_y") + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "Simulation year", y = "Young biomass (Mg/ha)",
                  caption = sprintf("Blue: LANDIS (cohorts <= %d yr, cells forested at year 0). Grey band: FIA mean (dashed) and 90%% interval. Dotted line: year %d,\nafter which all young cohorts were established in the simulation rather than in the initial communities.", lr$regen_age, lr$regen_age)) +
    theme_builder(10)
}

#' Species shares of young biomass
fig_regen_composition <- function(cmp, top = 10) {
  x <- data.table::copy(cmp)[order(-pmax(fia_share, landis_share))]
  keep <- x[seq_len(min(.N, top)), LandisName]
  x[, sp := data.table::fifelse(LandisName %in% keep, LandisName, "Other")]
  d <- rbind(x[, .(share = sum(fia_share)), by = sp][, src := "FIA plots"], x[, .(share = sum(landis_share)), by = sp][, src := "LANDIS cells"])
  d[, sp := factor(sp, levels = c(keep, "Other"))]
  ggplot2::ggplot(d, ggplot2::aes(src, 100 * share, fill = sp)) + ggplot2::geom_col(width = 0.6) +
    ggplot2::scale_fill_manual(NULL, values = c(scales::hue_pal(l = 55, c = 70)(length(keep)), "grey80")) +
    ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Share of young biomass (%)") + theme_builder(11) + ggplot2::theme(legend.position = "right")
}

#' Regeneration rate by site LAI: share of LANDIS cells / FIA plots in each LAI class with young cohorts or trees
#' (and FIA seedlings). Rates, not distributions, so a landscape that is denser overall does not shift the curves.
fig_regen_lai <- function(lr, fr, window, sp, cap = 16, min_plots = 15) {
  if (is.null(fr$plot_lai) || !nrow(lr$lai_young)) return(ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = "No LAI data (run step 3 for FIA plot LAI)") + ggplot2::theme_void())
  all_c <- lr$lai_all[year >= window[1] & year <= window[2], .(n_all = sum(N)), by = .(lai = pmin(lai, cap))]
  w <- lr$lai_young[year >= window[1] & year <= window[2] & SpeciesName %in% sp, .(n = sum(N)), by = .(LandisName = SpeciesName, lai = pmin(lai, cap))]
  w <- merge(w, all_c, by = "lai")[, .(LandisName, lai, rate = n / n_all, src = "LANDIS cells with young cohorts")]
  pl <- merge(fr$plots[LandisName %in% sp], fr$plot_lai, by = "PLT_CN")[, lai := pmin(floor(LAI_necn), cap)]
  f <- pl[, .(n = .N, young = mean(young_gm2 > 0), seed = mean(seedlings %in% TRUE)), by = .(LandisName, lai)][n >= min_plots]
  d <- rbind(w, f[, .(LandisName, lai, rate = young, src = "FIA plots with young trees")],
             f[, .(LandisName, lai, rate = seed, src = "FIA plots with seedlings")])
  d[, LandisName := factor(LandisName, levels = sp)]
  ggplot2::ggplot(d, ggplot2::aes(lai, 100 * rate, colour = src)) + ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_colour_manual(NULL, values = c("LANDIS cells with young cohorts" = REGEN_COL[["LANDIS"]], "FIA plots with young trees" = "grey15", "FIA plots with seedlings" = "#1baf7a")) +
    ggplot2::facet_wrap(~LandisName, scales = "free_y") + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = sprintf("Site LAI (NECN, whole numbers; %d+ grouped)", cap), y = "% of cells / plots in the LAI class",
                  caption = sprintf("LANDIS: the run's LAI maps, years %d-%d. FIA: site LAI the run's MaximumLAI / KLAI give on each plot; classes with < %d plots dropped.\nCompare the shapes: where LANDIS rises under dense canopy and FIA falls, establishment there is too easy (light curves).", window[1], window[2], min_plots)) +
    theme_builder(10) + ggplot2::theme(legend.position = "top")
}

#' Site LAI across the landscape through time against FIA plots near it
fig_landscape_lai <- function(lr, fr) {
  if (is.null(fr$plot_lai) || !nrow(lr$lai_all)) return(ggplot2::ggplot() + ggplot2::theme_void())
  q <- function(v, w, p) { o <- order(v); cw <- cumsum(w[o]) / sum(w); v[o][which(cw >= p)[1]] }
  a <- lr$lai_all[, .(q25 = q(lai, N, .25), q50 = q(lai, N, .5), q75 = q(lai, N, .75)), by = year]
  f <- floor(stats::quantile(fr$plot_lai$LAI_necn, c(.25, .5, .75)))
  ggplot2::ggplot(a, ggplot2::aes(year)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = f[1], ymax = f[3], fill = "grey88") +
    ggplot2::geom_hline(yintercept = f[2], linetype = "22", colour = "grey20") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), fill = REGEN_COL[["LANDIS"]], alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = q50), colour = REGEN_COL[["LANDIS"]], linewidth = 1) +
    ggplot2::labs(x = "Simulation year", y = "Site LAI (whole numbers)",
                  caption = "Blue: forested LANDIS cells, median and 25th-75th percentile. Grey: FIA plots near the landscape with the run's MaximumLAI / KLAI.") + theme_builder(11)
}

#' New cohorts per year from the reproduction log
fig_establish_log <- function(rl, sp) {
  d <- rl[SpeciesName %in% sp][, SpeciesName := factor(SpeciesName, levels = sp)]
  ggplot2::ggplot(d, ggplot2::aes(Time, per_1000_cells, colour = type)) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(NULL, values = c(Seed = "#2a78d6", Resprout = "#eb6834", Serotiny = "#4a3aa7", Planting = "grey40")) +
    ggplot2::facet_wrap(~SpeciesName, scales = "free_y") +
    ggplot2::labs(x = "Simulation year", y = "New cohorts per 1,000 forested cells", caption = "NECN-reproduction-log.csv (whole landscape).") +
    theme_builder(10) + ggplot2::theme(legend.position = "top")
}

#' Cohort biomass by cohort age per species: median and 25th-75th percentile, LANDIS against FIA
fig_cohort_biomass <- function(cb, sp, max_age = 400) {
  d <- cb$table[LandisName %in% sp & cls <= max_age][, LandisName := factor(LandisName, levels = sp)]
  ggplot2::ggplot(d, ggplot2::aes(cls, q50, colour = source, fill = source)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.18, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(NULL, values = REGEN_COL, labels = c(FIA = "FIA cohorts near the landscape", LANDIS = "LANDIS cohorts")) +
    ggplot2::scale_fill_manual(NULL, values = REGEN_COL, labels = c(FIA = "FIA cohorts near the landscape", LANDIS = "LANDIS cohorts")) +
    ggplot2::facet_wrap(~LandisName, scales = "free_y") + ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::labs(x = "Cohort age (10-year classes)", y = "Cohort biomass (Mg/ha, leaf + wood)",
                  caption = sprintf(paste0("Cohort = one species' trees in one 10-year age class on a cell (LANDIS) or plot (FIA), as the initial communities are built. ",
                                           "Lines: median; bands: 25th-75th percentile of cohorts that exist.\nLANDIS: years %d-%d pooled. FIA: tree ages from Lu et al. or the age models; classes with < %d FIA cohorts dropped. ",
                                           "LANDIS above FIA at young ages points to MaximumANPP / KLAI (young-cohort growth)."), cb$window[1], cb$window[2], cb$min_fia)) +
    theme_builder(10) + ggplot2::theme(legend.position = "top")
}

#' Cohort age structure per species: share of the species' cohorts in each age class
fig_cohort_ages <- function(cb, sp, max_age = 400) {
  d <- cb$table[LandisName %in% sp & cls <= max_age][, share := n / sum(n), by = .(LandisName, source)][, LandisName := factor(LandisName, levels = sp)]
  ggplot2::ggplot(d, ggplot2::aes(cls, 100 * share, colour = source)) + ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(NULL, values = REGEN_COL, labels = c(FIA = "FIA cohorts", LANDIS = "LANDIS cohorts")) +
    ggplot2::facet_wrap(~LandisName, scales = "free_y") +
    ggplot2::labs(x = "Cohort age (10-year classes)", y = "% of the species' cohorts",
                  caption = "Age structure of each species' cohorts. FIA classes with too few cohorts are dropped, so FIA shares are of the classes shown.") +
    theme_builder(10) + ggplot2::theme(legend.position = "top")
}
