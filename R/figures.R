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
