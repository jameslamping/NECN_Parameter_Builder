# light_fit.R -- establishment light curves (LightLAIShape, Scale, Location, Adjust) from FIA seedlings
#
# NECN v8 (PlugIn.SufficientLight): each year a species can establish on a site with probability
#   P(light) = min(1, Adjust * (WeibullPDF(siteLAI; Shape, Scale) + Location))
# multiplied by seed arrival and the climate establishment probability, which NECN models separately.
#
# Data: one row per FIA plot x LANDIS species where the species has a seed source (see seed_source()).
#   y = seedlings of the species present on the plot's microplots (FIA SEEDLING table)
#   x = site LAI on the plot, on the scale NECN will evaluate the curve at:
#       "current_code": cohort-sum LAI NECN v8 computes from the plot's cohorts with your MaximumLAI/KLAI
#                       (before the growth-neutral rescale; the export rescales the curve with everything else)
#       "stand_level" : stand-level LAI from the proposed code change
#       "observed"    : foliage x SLA (no NECN parameters involved)
# The curve is fitted by maximum likelihood (binomial), not least squares on binned proportions, and kept
# only when it beats a constant presence rate (likelihood-ratio test, 3 df). Otherwise the species gets a
# flat curve (Shape 1, Scale 1000: light never limits establishment).
# FIA presence also reflects climate, seed supply and the years a seedling persists, so by default the
# curve is normalized to peak at 1 over the observed LAI range: light then only scales establishment down
# where the site is too dark (or too open) for the species.

light_prob <- function(x, a, b, c, adj) pmin(1, adj * (stats::dweibull(x, a, b) + c))

#' Seedling presence per plot x LANDIS species
seedling_presence <- function(fia, plt_cn, species_map) {
  if (is.null(fia$SEEDLING)) stop("No SEEDLING table loaded. Add <ST>_SEEDLING.csv files to the FIA folder (the download step gets them) and reload.")
  m <- expand_species_map(species_map)
  s <- fia$SEEDLING[PLT_CN %in% plt_cn & SPCD %in% m$SPCD]
  s[, n := data.table::fifelse(!is.na(TREECOUNT_CALC), as.numeric(TREECOUNT_CALC), as.numeric(TREECOUNT))]
  s <- merge(s, m, by = "SPCD")[, .(seedlings = sum(n, na.rm = TRUE)), by = .(PLT_CN, LandisName)]
  s[seedlings > 0]
}

#' Plots with a seed source for each LANDIS species
#' @param how "plot" (a live tree of the species >= min_dia inches on the plot), "radius" (on any selected
#'   plot within radius_km; FIA coordinates are fuzzed by up to ~1.6 km), or "none" (every plot)
#' @return data.table(PLT_CN, LandisName)
seed_source <- function(fia, plots, species_map, how = c("plot", "radius", "none"), min_dia = 5, radius_km = 5) {
  how <- match.arg(how)
  m <- expand_species_map(species_map)
  if (how == "none") return(data.table::CJ(PLT_CN = plots$PLT_CN, LandisName = unique(m$LandisName)))
  src <- unique(merge(fia$TREE[PLT_CN %in% plots$PLT_CN & STATUSCD == 1 & !is.na(DIA) & DIA >= min_dia, .(PLT_CN, SPCD)], m, by = "SPCD")[, .(PLT_CN, LandisName)])
  if (how == "plot") return(src)
  p <- plots[!is.na(LAT) & !is.na(LON)]
  pts <- sf::st_transform(sf::st_as_sf(p, coords = c("LON", "LAT"), crs = 4269), 5070)
  nb <- sf::st_is_within_distance(pts, pts, dist = radius_km * 1000)
  nb_dt <- data.table::data.table(PLT_CN = rep(p$PLT_CN, lengths(nb)), nbr = p$PLT_CN[unlist(nb)])
  unique(merge(nb_dt, src, by.x = "nbr", by.y = "PLT_CN", allow.cartesian = TRUE)[, .(PLT_CN, LandisName)])
}

#' Assemble the fitting table
#' @param x_tab data.table(PLT_CN, x) -- site LAI per plot (0 for plots with no live trees)
light_data <- function(x_tab, presence, source) {
  d <- merge(source, x_tab, by = "PLT_CN")
  d <- merge(d, presence, by = c("PLT_CN", "LandisName"), all.x = TRUE)
  d[, y := as.integer(!is.na(seedlings) & seedlings > 0)]
  d[!is.na(x)][]
}

#' Maximum-likelihood fit of NECN's light curve to presence / absence
#' @param c_max upper bound for Location (NECN limit 1; keep <= f so the rescaled Location stays <= 1)
fit_light_curve <- function(x, y, c_max = 1) {
  x <- pmax(x, 1e-3)
  nll <- function(p) {
    pr <- light_prob(x, p[1], p[2], p[3], p[4])
    pr <- pmin(pmax(pr, 1e-9), 1 - 1e-9)
    -sum(y * log(pr) + (1 - y) * log(1 - pr))
  }
  q <- stats::quantile(x, c(.25, .5, .75, .95))
  starts <- expand.grid(a = c(0.8, 1.5, 3, 5), b = unique(pmax(0.2, q)), c = c(0, 0.01))
  lo <- c(0.3, 0.05, 0, 1e-4); hi <- c(6, 1000, c_max, 100)
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    s <- unlist(starts[i, ])
    adj0 <- min(50, max(1e-3, mean(y) / mean(stats::dweibull(x, s[1], s[2]) + s[3])))
    o <- tryCatch(stats::optim(c(s, adj0), nll, method = "L-BFGS-B", lower = lo, upper = hi), error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(list(a = NA_real_, b = NA_real_, c = NA_real_, adj = NA_real_, nll = NA_real_, converged = FALSE))
  p <- pmin(pmax(best$par, lo), hi)
  # deviance explained vs. a constant presence rate
  m0 <- mean(y); nll0 <- -sum(y * log(m0) + (1 - y) * log(1 - m0))
  list(a = p[1], b = p[2], c = p[3], adj = p[4], nll = best$value, dev_expl = 1 - best$value / nll0,
       p_lai = stats::pchisq(2 * max(0, nll0 - best$value), df = 3, lower.tail = FALSE), converged = best$convergence == 0)
}

#' Fit with a bare-ground condition on P(site LAI 0).
#' NECN's curve is Adjust x (WeibullPDF + Location). With Shape > 1 the PDF is 0 at LAI 0 (and far out on the LAI axis),
#' so the bare-ground value is Location's floor, which also applies under a closed canopy. With Shape <= 1 the curve
#' is highest at LAI 0 and falls with LAI.
#'   type "min":   P(0) >= p0. The fit chooses between a curve falling from bare ground (Shape <= 1) and a rising
#'                 curve with a floor high enough (Shape > 1).
#'   type "exact": P(0) = p0 with Shape > 1 (Location = p0 x maxPDF / (1 - p0) when normalized, else p0 / Adjust).
#' With normalization (peak 1 over (0, xmax]) P(0) = Location / (maxPDF + Location) for Shape > 1.
#' Location is capped at c_max (the achieved P(0) is then lower).
fit_light_curve_p0 <- function(x, y, p0, xmax, c_max = 1, normalize = TRUE, type = c("min", "exact")) {
  type <- match.arg(type)
  x <- pmax(x, 1e-3)
  g <- seq(max(0.02, xmax / 400), xmax, length.out = 400)
  need <- function(a, b, adj) if (normalize) p0 * max(stats::dweibull(g, a, b)) / (1 - p0) else p0 / adj
  loc_of <- function(a, b, cf, adj) {
    if (type == "exact") return(min(c_max, need(a, b, adj)))
    if (a <= 1) return(cf)
    min(c_max, max(cf, need(a, b, adj)))
  }
  nll <- function(p) {
    cf <- if (type == "exact") 0 else p[4]; adj <- if (type == "exact") p[3] else p[3]
    loc <- loc_of(p[1], p[2], cf, adj)
    pr <- pmin(pmax(light_prob(x, p[1], p[2], loc, adj), 1e-9), 1 - 1e-9)
    -sum(y * log(pr) + (1 - y) * log(1 - pr))
  }
  q <- stats::quantile(x, c(.25, .5, .75, .95))
  if (type == "exact") {
    starts <- expand.grid(a = c(1.2, 1.5, 3, 5), b = unique(pmax(0.2, q)), c = 0)
    lo <- c(1.01, 0.05, 1e-4); hi <- c(6, 1000, 100)
  } else {
    starts <- expand.grid(a = c(0.6, 0.9, 1.5, 3, 5), b = unique(pmax(0.2, q)), c = c(0, 0.02))
    lo <- c(0.3, 0.05, 1e-4, 0); hi <- c(6, 1000, 100, c_max)
  }
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    s <- unlist(starts[i, ])
    adj0 <- min(50, max(1e-3, mean(y) / mean(stats::dweibull(x, s[1], s[2]) + s[3] + 0.01)))
    p_start <- if (type == "exact") c(s[1], s[2], adj0) else c(s[1], s[2], adj0, s[3])
    o <- tryCatch(stats::optim(p_start, nll, method = "L-BFGS-B", lower = lo, upper = hi), error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(NULL)
  p <- pmin(pmax(best$par, lo), hi)
  cf <- if (type == "exact") 0 else p[4]
  m0 <- mean(y); nll0 <- -sum(y * log(m0) + (1 - y) * log(1 - m0))
  list(a = p[1], b = p[2], c = loc_of(p[1], p[2], cf, p[3]), adj = p[3], nll = best$value, dev_expl = 1 - best$value / nll0)
}

#' Rescale Adjust so the curve peaks at 1 over [0, xmax]
normalize_light <- function(a, b, c, adj, xmax) {
  g <- seq(max(0.02, xmax / 400), xmax, length.out = 400)
  pk <- max(adj * (stats::dweibull(g, a, b) + c))
  if (!is.finite(pk) || pk <= 0) return(adj)
  adj / pk
}

#' Binned presence for display: equal-count bins of x with Wilson 90% intervals
bin_presence <- function(d, nbins = 12) {
  d <- data.table::copy(d)
  br <- unique(stats::quantile(d$x, seq(0, 1, length.out = nbins + 1)))
  if (length(br) < 3) br <- unique(c(min(d$x), max(d$x)))
  d[, bin := cut(x, br, include.lowest = TRUE)]
  b <- d[, .(x = mean(x), n = .N, k = sum(y)), by = bin][order(x)]
  z <- 1.645
  b[, p := k / n]
  b[, `:=`(lo = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
           hi = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n))]
  b[]
}

#' Flat curve in NECN's form: Shape 1, Scale 1000 (Weibull PDF ~0.001 everywhere) plus Location
.flat_curve <- function(level, c_max) {
  loc <- min(c_max, 1); base <- stats::dweibull(1, 1, 1000) + loc
  list(a = 1, b = 1000, c = loc, adj = level / base)
}

#' Fit every species
#' @param ld fitting table from light_data()
#' @param current data.table(LandisName, LightLAIShape, LightLAIScale, LightLAILocation, LightLAIAdjust) or NULL
#' @return list(params, data, bins, settings)
#' @param bare "lai05": where the fitted curve drops toward 0 at bare ground (site LAI 0, a cell with no cohorts),
#'   refit with the bare-ground probability fixed at the curve's value at LAI `bare_lai` (nearly full light);
#'   "fit": keep the unconstrained fit (with Shape > 1 it is 0 at LAI 0)
fit_light_all <- function(ld, current = NULL, min_pres = 20, min_plots = 50, normalize = TRUE, c_max = 1, alpha = 0.05,
                          bare = c("lai05", "fit"), bare_lai = 0.5, progress = function(msg, value) NULL) {
  bare <- match.arg(bare)
  spp <- sort(unique(ld$LandisName))
  res <- data.table::rbindlist(lapply(seq_along(spp), function(i) {
    s <- spp[i]; progress(paste("Fitting", s), i / length(spp))
    d <- ld[LandisName == s]
    out <- data.table::data.table(LandisName = s, n_plots = nrow(d), n_present = sum(d$y), pct_present = 100 * mean(d$y),
                                  n_low_lai = sum(d$x < 0.5), x_p95 = unname(stats::quantile(d$x, .95)))
    if (sum(d$y) < min_pres || nrow(d) < min_plots || sum(d$y) == nrow(d))
      return(out[, `:=`(Shape = NA_real_, Scale = NA_real_, Location = NA_real_, Adjust_fit = NA_real_, Adjust = NA_real_,
                        dev_expl = NA_real_, p_lai = NA_real_, status = "too few presences")])
    f <- fit_light_curve(d$x, d$y, c_max)
    st <- "fitted"
    if (!is.finite(f$p_lai) || f$p_lai > alpha) {
      fl <- .flat_curve(mean(d$y), c_max); fl$dev_expl <- f$dev_expl; fl$p_lai <- f$p_lai; f <- fl
      st <- "flat (no LAI effect)"
    }
    adj <- if (normalize) normalize_light(f$a, f$b, f$c, f$adj, out$x_p95) else f$adj
    out[, `:=`(Shape = f$a, Scale = f$b, Location = f$c, Adjust_fit = f$adj, Adjust = min(100, adj), dev_expl = f$dev_expl,
               p_lai = f$p_lai, status = st, bare_rule = "as fitted", p_lai0_unconstrained = light_p0(f$a, f$b, f$c, min(100, adj)))]
    # bare ground: the Weibull PDF is 0 at LAI 0 when Shape > 1, so an empty cell would get no establishment
    if (bare == "lai05" && st == "fitted" && f$a > 1) {
      target <- light_prob(bare_lai, f$a, f$b, f$c, min(100, adj))
      if (out$p_lai0_unconstrained < target - 0.02) {
        g2 <- fit_light_curve_p0(d$x, d$y, min(target, 0.99), out$x_p95, c_max, normalize, type = "min")
        if (!is.null(g2)) {
          adj2 <- if (normalize) normalize_light(g2$a, g2$b, g2$c, g2$adj, out$x_p95) else g2$adj
          out[, `:=`(Shape = g2$a, Scale = g2$b, Location = g2$c, Adjust_fit = g2$adj, Adjust = min(100, adj2), dev_expl = g2$dev_expl,
                     bare_rule = sprintf("bare ground >= value at LAI %s", bare_lai))]
        }
      }
    }
    out
  }), fill = TRUE)
  # light probability on bare ground (site LAI exactly 0: a cell with no cohorts). For Shape > 1 the Weibull
  # term is 0 there, so it is Adjust x Location. Few fitting plots below LAI 0.5 means it is extrapolated.
  res[, p_lai0 := light_p0(Shape, Scale, Location, Adjust)]
  if (!"bare_rule" %in% names(res)) res[, bare_rule := NA_character_]
  if (!is.null(current)) res <- merge(res, current, by = "LandisName", all.x = TRUE)
  bins <- ld[, bin_presence(.SD), by = LandisName]
  list(params = res[order(-n_present)], data = ld, bins = bins,
       settings = list(min_pres = min_pres, min_plots = min_plots, normalize = normalize, c_max = c_max, bare = bare, bare_lai = bare_lai))
}

#' Light probability at site LAI exactly 0 (bare ground)
light_p0 <- function(a, b, c, adj) {
  w0 <- data.table::fifelse(is.na(a), NA_real_, data.table::fifelse(a > 1, 0, data.table::fifelse(a == 1, 1 / b, Inf)))
  pmin(1, adj * (w0 + c))
}

#' Final LightLAI values per species. Order: user choice of curve ("fit", "current", or "flat" = light never
#' limits), then a borrowed species' fit, then an optional bare-ground target (p0_use: LightLAILocation is set
#' so the probability at site LAI 0 equals it, keeping Shape, Scale and Adjust).
#' @param choices data.table(LandisName, use, borrow_from, p0_use, note) or NULL
light_final <- function(lfit, choices, necn_tab) {
  cols <- c("LightLAIShape", "LightLAIScale", "LightLAILocation", "LightLAIAdjust")
  cur <- necn_tab[, c("SpeciesCode", intersect(cols, names(necn_tab))), with = FALSE]
  data.table::setnames(cur, "SpeciesCode", "LandisName")
  for (cc in setdiff(cols, names(cur))) cur[, (cc) := NA_real_]
  pk <- intersect(c("LandisName", "n_plots", "n_present", "pct_present", "n_low_lai", "Shape", "Scale", "Location", "Adjust",
                    "dev_expl", "p_lai", "p_lai0", "p_lai0_unconstrained", "bare_rule", "status", "x_p95"), names(lfit$params))
  p <- merge(cur, lfit$params[, pk, with = FALSE], by = "LandisName", all.x = TRUE)
  empty <- data.table::data.table(LandisName = character(), use = character(), borrow_from = character(), p0_use = numeric(), note = character())
  if (is.null(choices) || !nrow(choices)) choices <- empty
  choices <- data.table::copy(choices)
  for (cc in setdiff(names(empty), names(choices))) choices[, (cc) := empty[[cc]][NA_integer_]]
  choices[, p0_use := suppressWarnings(as.numeric(p0_use))]
  p <- merge(p, choices[, names(empty), with = FALSE], by = "LandisName", all.x = TRUE)
  p[is.na(use) | !nzchar(use), use := data.table::fifelse(!is.na(Shape), "fit", "current")]
  c_max <- if (is.null(lfit$settings$c_max)) 1 else lfit$settings$c_max
  p[, `:=`(Shape_used = LightLAIShape, Scale_used = LightLAIScale, Location_used = LightLAILocation, Adjust_used = LightLAIAdjust, light_source = "current table")]
  p[use == "fit" & !is.na(Shape), `:=`(Shape_used = Shape, Scale_used = Scale, Location_used = Location, Adjust_used = Adjust, light_source = "FIA seedlings")]
  for (i in which(!is.na(p$borrow_from) & nzchar(p$borrow_from))) {
    j <- match(p$borrow_from[i], p$LandisName)
    if (!is.na(j) && !is.na(p$Shape[j]))
      p[i, `:=`(Shape_used = p$Shape[j], Scale_used = p$Scale[j], Location_used = p$Location[j], Adjust_used = p$Adjust[j],
                light_source = paste("borrowed from", p$borrow_from[i]))]
  }
  fl <- .flat_curve(1, c_max)
  p[use == "flat", `:=`(Shape_used = fl$a, Scale_used = fl$b, Location_used = fl$c, Adjust_used = fl$adj, light_source = "flat (light never limits)")]
  # bare-ground override: refit the species' curve with its probability at site LAI 0 fixed at p0_use
  # (only for the species' own fit; a uniform Location shift would raise the whole curve)
  p[, p0_note := ""]
  norm <- !isFALSE(lfit$settings$normalize)
  for (i in which(!is.na(p$p0_use))) {
    sp <- p$LandisName[i]; d <- lfit$data[LandisName == sp]
    if (!nrow(d) || is.na(p$x_p95[i])) { p$p0_note[i] <- "no FIA data to refit"; next }
    if (p$p0_use[i] >= 0.99) { fl1 <- .flat_curve(1, c_max)
      p[i, `:=`(Shape_used = fl1$a, Scale_used = fl1$b, Location_used = fl1$c, Adjust_used = fl1$adj, light_source = "flat (bare ground set to 1)")]; next }
    g2 <- fit_light_curve_p0(d$x, d$y, p$p0_use[i], p$x_p95[i], c_max, norm, type = "exact")
    if (is.null(g2)) { p$p0_note[i] <- "refit failed"; next }
    adj2 <- if (norm) normalize_light(g2$a, g2$b, g2$c, g2$adj, p$x_p95[i]) else g2$adj
    p[i, `:=`(Shape_used = g2$a, Scale_used = g2$b, Location_used = g2$c, Adjust_used = min(100, adj2),
              light_source = sprintf("FIA seedlings, refitted with bare ground = %s", p0_use))]
  }
  p[, p0_used := light_p0(Shape_used, Scale_used, Location_used, Adjust_used)]
  p[, p10_used := light_prob(10, Shape_used, Scale_used, Location_used, Adjust_used)]
  p[]
}
