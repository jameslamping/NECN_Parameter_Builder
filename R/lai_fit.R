# lai_fit.R -- observed plot LAI and species MaximumLAI / KLAI
#
# NECN cohort LAI = MaximumLAI * C / (KLAI + C), C = wood carbon (g C/m2) = wood biomass x 0.47.
# KLAI is the wood carbon at which LAI reaches half of MaximumLAI, so it is estimated against
# stand wood (not age). Three estimates per species:
#   mono   : plots where the species holds >= dom_min of live wood; y = plot LAI
#   mixed  : plots where it holds >= share_min; y = the species' own LAI / its wood share
#            ("pure-stand equivalent"), weighted by share. More data, but biased high for
#            minority hardwoods and subcanopy trees (more foliage per unit wood than canopy dominants).
#   pooled : group value = geometric mean of the monospecific estimates of the group's species
#            (species with >= 20 monospecific plots, else all with an estimate), so each species
#            counts once; each species' mono estimate is shrunk toward it with weight n / (n + n0)
#            on the log scale (species without one get the group value).

LBAC_TO_GM2 <- 0.112085   # lb/acre -> g/m2
C_FRAC <- 0.47            # NECN carbon fraction

#' Plot LAI from live-tree foliage (DRYBIO_FOLIAGE) x SLA, and live wood (DRYBIO_AG) by species
#' @return list(plots = one row per plot, plot_sp = one row per plot x SPCD)
compute_plot_lai <- function(fia, plots, sla_tab, max_missing = 0.1) {
  tr <- fia$TREE[PLT_CN %in% plots$PLT_CN & STATUSCD == 1 & TPA_UNADJ > 0]
  tr <- merge(tr, sla_tab[, .(SPCD, sla)], by = "SPCD", all.x = TRUE)
  tr[, `:=`(wood = DRYBIO_AG * TPA_UNADJ * LBAC_TO_GM2, fol = DRYBIO_FOLIAGE * TPA_UNADJ * LBAC_TO_GM2)]
  ps <- tr[, .(wood_gm2 = sum(wood, na.rm = TRUE), foliage_gm2 = sum(fol, na.rm = TRUE),
               lai = sum(fol / 1000 * sla, na.rm = TRUE), n_trees = .N,
               miss = sum(is.na(fol) | is.na(sla)) / .N), by = .(PLT_CN, SPCD)]
  pl <- ps[, .(LAI = sum(lai), wood_gm2 = sum(wood_gm2), foliage_gm2 = sum(foliage_gm2),
               missing_frac = sum(miss * n_trees) / sum(n_trees)), by = PLT_CN]
  pl <- merge(plots[, .(PLT_CN, STATECD, INVYR, LAT, LON, ELEV, STDAGE)], pl, by = "PLT_CN")
  pl <- pl[missing_frac <= max_missing & wood_gm2 > 0]
  list(plots = pl[], plot_sp = ps[PLT_CN %in% pl$PLT_CN][])
}

#' One row per plot x LANDIS species: the species' own LAI and wood, its share, plot totals
species_plot_table <- function(plot_lai, species_map) {
  m <- expand_species_map(species_map)
  s <- merge(plot_lai$plot_sp, m, by = "SPCD")[, .(lai_s = sum(lai), wood_s = sum(wood_gm2)), by = .(PLT_CN, LandisName)]
  s <- merge(s, plot_lai$plots[, .(PLT_CN, LAI, wood_gm2, STDAGE)], by = "PLT_CN")
  s[, `:=`(share = wood_s / wood_gm2, C = wood_gm2 * C_FRAC)]
  s[share > 0][]
}

#' kept for the earlier interface: share of each plot's live wood by LANDIS species
species_shares <- function(plot_sp, species_map) {
  m <- expand_species_map(species_map)
  tot <- plot_sp[, .(tot = sum(wood_gm2)), by = PLT_CN]
  s <- merge(plot_sp, m, by = "SPCD")[, .(w = sum(wood_gm2)), by = .(PLT_CN, LandisName)]
  merge(s, tot, by = "PLT_CN")[, share := w / tot][, .(PLT_CN, LandisName, share)]
}

# core fitter: MaximumLAI from the plateau, KLAI by (weighted) least squares
.fit_curve <- function(x, plateau_top, plateau_q, max_cap, w = NULL) {
  M <- min(max_cap, unname(stats::quantile(x[C >= stats::quantile(C, 1 - plateau_top), y], plateau_q)))
  if (is.null(w)) w <- rep(1, nrow(x))
  f <- tryCatch(minpack.lm::nlsLM(y ~ M * C / (K + C), data = x, weights = w, start = list(K = unname(stats::median(x$C))),
                                  lower = 10, upper = 2e5), error = function(e) NULL)
  c(M = M, K = if (is.null(f)) NA_real_ else unname(stats::coef(f)[1]))
}
.fit_with_boot <- function(x, nboot, weighted, ...) {
  a <- list(...)
  fit_i <- function(i) do.call(.fit_curve, c(list(x[i], w = if (weighted) x$share[i] else NULL), a))
  f <- fit_i(seq_len(nrow(x)))
  b <- if (nboot > 0) t(vapply(seq_len(nboot), function(k) fit_i(sample(nrow(x), replace = TRUE)), c(M = 0, K = 0))) else NULL
  list(M = f[["M"]], K = f[["K"]],
       lo = if (is.null(b)) NA_real_ else stats::quantile(b[, "K"], .1, na.rm = TRUE),
       hi = if (is.null(b)) NA_real_ else stats::quantile(b[, "K"], .9, na.rm = TRUE))
}

#' Fit all three estimates and pick a recommended one per species
#' @param groups data.table(LandisName, group, softwood) -- functional group per species
fit_lai_params <- function(sp_tab, groups, dom_min = 0.8, share_min = 0.3, plateau_top = 0.2, plateau_q = 0.9,
                           max_cap = 20, min_n = 12, nboot = 200, n0 = 30, seed = 1,
                           progress = function(msg, value) NULL) {
  set.seed(seed)
  mono_d <- sp_tab[share >= dom_min][, y := LAI]
  mix_d  <- sp_tab[share >= share_min][, y := lai_s / share]
  spp <- sort(unique(c(groups$LandisName)))
  args <- list(plateau_top = plateau_top, plateau_q = plateau_q, max_cap = max_cap)
  res <- data.table::rbindlist(lapply(seq_along(spp), function(i) {
    s <- spp[i]; progress(paste("Fitting", s), i / (length(spp) + 2))
    a <- mono_d[LandisName == s]; b <- mix_d[LandisName == s]
    fa <- if (nrow(a) >= min_n) do.call(.fit_with_boot, c(list(a, nboot, FALSE), args)) else list(M = NA_real_, K = NA_real_, lo = NA_real_, hi = NA_real_)
    fb <- if (nrow(b) >= min_n) do.call(.fit_with_boot, c(list(b, max(0, nboot %/% 2), TRUE), args)) else list(M = NA_real_, K = NA_real_, lo = NA_real_, hi = NA_real_)
    data.table::data.table(LandisName = s, n_any = sp_tab[LandisName == s, .N], n_mono = nrow(a), n_mixed = nrow(b),
                           MaxLAI_mono = fa$M, KLAI_mono = fa$K, KLAI_mono_lo = fa$lo, KLAI_mono_hi = fa$hi,
                           MaxLAI_mixed = fb$M, KLAI_mixed = fb$K, KLAI_mixed_lo = fb$lo, KLAI_mixed_hi = fb$hi)
  }))
  res <- merge(res, groups[, .(LandisName, group, softwood)], by = "LandisName", all.x = TRUE)

  # group values from the species-level monospecific estimates (each species counts once)
  progress("Group values", (length(spp) + 1) / (length(spp) + 2))
  gfit <- res[!is.na(group), {
    ok <- !is.na(KLAI_mono) & n_mono >= 20
    if (!any(ok)) ok <- !is.na(KLAI_mono)
    .(n_species = sum(ok), species = paste(LandisName[ok], collapse = ", "),
      MaxLAI_group = if (any(ok)) exp(mean(log(MaxLAI_mono[ok]))) else NA_real_,
      KLAI_group = if (any(ok)) exp(mean(log(KLAI_mono[ok]))) else NA_real_)
  }, by = group]
  res <- merge(res, gfit[, .(group, MaxLAI_group, KLAI_group)], by = "group", all.x = TRUE)
  res[, w_own := n_mono / (n_mono + n0)]
  res[, `:=`(MaxLAI_pooled = data.table::fifelse(is.na(MaxLAI_mono), MaxLAI_group, exp(w_own * log(MaxLAI_mono) + (1 - w_own) * log(MaxLAI_group))),
             KLAI_pooled   = data.table::fifelse(is.na(KLAI_mono),   KLAI_group,   exp(w_own * log(KLAI_mono)   + (1 - w_own) * log(KLAI_group))))]

  res[, confidence := data.table::fifelse(n_mono >= 50, "good", data.table::fifelse(n_mono >= 20, "fair", data.table::fifelse(n_mono >= min_n, "low", "none")))]
  res[, disagree := !is.na(KLAI_mono) & !is.na(KLAI_mixed) & abs(log(KLAI_mixed / KLAI_mono)) > log(1.3)]
  # recommendation: mono with >= 20 plots; mono shrunk to the group (pooled) with min_n-19 plots;
  # with fewer, the mixed-plot fit for softwoods with >= 20 mixed plots, else the group value
  res[, recommended := data.table::fcase(
    n_mono >= 20 & !is.na(KLAI_mono), "mono",
    !is.na(KLAI_mono) & !is.na(KLAI_pooled), "pooled",
    softwood %in% TRUE & n_mixed >= 20 & !is.na(KLAI_mixed), "mixed",
    !is.na(KLAI_pooled), "pooled",
    !is.na(KLAI_mixed), "mixed",
    default = "none")]
  res[, why := data.table::fcase(
    recommended == "mono", sprintf("%d monospecific plots", n_mono),
    recommended == "pooled" & !is.na(KLAI_mono), sprintf("%d monospecific plots; shrunk toward the %s group (%.0f%% own data)", n_mono, group, 100 * w_own),
    recommended == "mixed", sprintf("%d monospecific plots; softwood with %d plots at >= %d%% share", n_mono, n_mixed, round(share_min * 100)),
    recommended == "pooled", sprintf("%d monospecific plots; %s group value", n_mono, group),
    default = "no usable data; keep current value or borrow")]
  pick <- function(m, a, b, c) data.table::fcase(m == "mono", a, m == "mixed", b, m == "pooled", c, default = NA_real_)
  res[, check := !is.na(KLAI_group) & recommended != "none" &
        abs(log(pick(recommended, KLAI_mono, KLAI_mixed, KLAI_pooled) / KLAI_group)) > log(1.5)]
  progress("done", 1)
  data.table::setorder(res, -n_mono)
  list(params = res[], mono = mono_d[], mixed = mix_d[], groups = gfit[],
       settings = list(dom_min = dom_min, share_min = share_min, n0 = n0))
}

#' NECN curve points for plotting
curve_xy <- function(M, K, xmax_gm2, n = 150) {
  w <- seq(0, xmax_gm2, length.out = n)
  data.table::data.table(wood_Mgha = w / 100, LAI = M * (w * C_FRAC) / (K + w * C_FRAC))
}
