# cohorts.R -- tree ages and NECN-style cohorts (LANDIS species x 10-year age class) from FIA trees,
# and the site LAI NECN would compute for those cohorts. Shared by the establishment (LightLAI) step
# and the initial-communities builder.
#
# Ages: Lu et al. (2025, Forest Ecology and Management 122603) tree ages by tree CN where a folder of
# <ST>_TREEAGE.csv files is given; otherwise a per-species linear model TOTAGE ~ DIA + HT fitted on FIA
# trees with a recorded total age (DIA only when HT is missing; a softwood/hardwood model for species
# with too few aged trees).

#' Lu et al. tree ages from a folder of <ST>_TREEAGE.csv files (columns CN, Tree_Age)
read_lu_ages <- function(dir, states = NULL) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  f <- list.files(dir, "_TREEAGE\\.csv$", full.names = TRUE)
  if (!is.null(states)) f <- f[sub("_TREEAGE\\.csv$", "", basename(f)) %in% states]
  if (!length(f)) return(NULL)
  a <- data.table::rbindlist(lapply(f, data.table::fread, select = c("CN", "Tree_Age"), colClasses = list(character = "CN")))
  a[!is.na(Tree_Age) & Tree_Age > 0][!duplicated(CN)]
}

#' Per-species age models: TOTAGE ~ DIA + HT and TOTAGE ~ DIA; species with fewer than `min_n`
#' aged trees use the softwood or hardwood model
fit_age_models <- function(tree, ref, min_n = 30) {
  d <- tree[!is.na(TOTAGE) & TOTAGE > 0 & !is.na(DIA) & DIA > 0, .(SPCD, TOTAGE, DIA, HT)]
  if (!nrow(d)) stop("No FIA trees with a recorded total age (TOTAGE); give a folder of Lu et al. tree ages instead.")
  d[, sw := data.table::fifelse(ref$SFTWD_HRDWD[match(SPCD, ref$SPCD)] %in% "H", "H", "S")]
  fit2 <- function(x) {
    h <- x[!is.na(HT) & HT > 0]
    list(dh = if (nrow(h) >= 10) stats::coef(stats::lm(TOTAGE ~ DIA + HT, h)) else NULL,
         d = stats::coef(stats::lm(TOTAGE ~ DIA, x)), n = nrow(x))
  }
  sp <- d[, .N, by = SPCD][N >= min_n, SPCD]
  list(species = stats::setNames(lapply(sp, function(s) fit2(d[SPCD == s])), sp),
       group = list(S = fit2(d[sw == "S"]), H = if (d[sw == "H", .N] >= 10) fit2(d[sw == "H"]) else fit2(d[sw == "S"])),
       ref = ref[, .(SPCD, SFTWD_HRDWD)])
}

.predict_age <- function(m, DIA, HT) {
  out <- rep(NA_real_, length(DIA))
  if (!is.null(m$dh)) { i <- !is.na(HT) & HT > 0; out[i] <- m$dh[1] + m$dh[2] * DIA[i] + m$dh[3] * HT[i] }
  j <- is.na(out); out[j] <- m$d[1] + m$d[2] * DIA[j]
  out
}

#' Age for every tree: Lu et al. age where available, else the age model
#' @return tree with columns age (years, >= 1) and age_src ("Lu et al.", "species model", "group model")
assign_tree_ages <- function(tree, models, lu = NULL) {
  t <- data.table::copy(tree)
  t[, `:=`(age = NA_real_, age_src = NA_character_)]
  if (!is.null(lu)) { t[lu, on = "CN", `:=`(age = as.numeric(i.Tree_Age), age_src = "Lu et al.")] }
  todo <- t[is.na(age), unique(SPCD)]
  for (s in todo) {
    i <- which(t$SPCD == s & is.na(t$age))
    m <- models$species[[as.character(s)]]; src <- "species model"
    if (is.null(m)) { m <- models$group[[if (models$ref$SFTWD_HRDWD[match(s, models$ref$SPCD)] %in% "H") "H" else "S"]]; src <- "group model" }
    data.table::set(t, i, "age", .predict_age(m, t$DIA[i], t$HT[i]))
    data.table::set(t, i, "age_src", src)
  }
  t[!is.finite(age) | age < 1, age := 1]
  t[]
}

#' NECN-style cohorts on FIA plots: live trees of the LANDIS species grouped by species and
#' 10-year age class; wood = DRYBIO_AG, leaf = DRYBIO_FOLIAGE (g/m2)
plot_cohorts <- function(fia, plt_cn, species_map, models, lu = NULL) {
  m <- expand_species_map(species_map)
  tr <- fia$TREE[PLT_CN %in% plt_cn & STATUSCD == 1 & TPA_UNADJ > 0 & SPCD %in% m$SPCD]
  tr <- assign_tree_ages(tr, models, lu)
  tr <- merge(tr, m, by = "SPCD")
  tr[, age_class := pmax(10, round(age / 10) * 10)]
  co <- tr[, .(wood_gm2 = sum(DRYBIO_AG * TPA_UNADJ * LBAC_TO_GM2, na.rm = TRUE),
               leaf_gm2 = sum(DRYBIO_FOLIAGE * TPA_UNADJ * LBAC_TO_GM2, na.rm = TRUE), n_trees = .N),
           by = .(PLT_CN, LandisName, age_class)]
  attr(co, "age_sources") <- tr[, .N, by = age_src]
  co[wood_gm2 > 0]
}

#' Site LAI NECN would compute for these cohorts, before any rescale.
#'   current_code: sum over cohorts of max(MinLAI, MaximumLAI * C / (KLAI + C))  (NECN v8 as coded)
#'   stand_level : sum over cohorts of MaximumLAI * C / (KLAI + C_site)         (proposed code change)
#' @param params data.table(LandisName, MaximumLAI, KLAI, MinLAI)
necn_site_lai <- function(cohorts, params, mode = c("current_code", "stand_level")) {
  mode <- match.arg(mode)
  co <- merge(cohorts, params[, .(LandisName, M = MaximumLAI, K = KLAI, MinLAI)], by = "LandisName")
  co[, C := wood_gm2 * C_FRAC]
  if (mode == "current_code") co[, lai := pmax(MinLAI, M * C / (K + C))]
  else co[, lai := M * C / (K + sum(C)), by = PLT_CN]
  co[, .(LAI_necn = sum(lai), n_cohorts = .N), by = PLT_CN]
}

#' Site LAI under stand-level NECN computed directly from each plot's wood by species (no ages needed):
#'   sum over species of MaximumLAI_s * C_s / (KLAI_s + C_site), C = wood carbon of the LANDIS species.
#' Identical to necn_site_lai(..., "stand_level") because the stand-level formula does not depend on how
#' a species' wood is split into cohorts.
#' @param sp_tab species_plot_table() output (PLT_CN, LandisName, wood_s)
stand_site_lai <- function(sp_tab, params) {
  d <- merge(sp_tab[, .(PLT_CN, LandisName, wood_s)], params[, .(LandisName, M = MaximumLAI, K = KLAI)], by = "LandisName")
  d[, C := wood_s * C_FRAC][, C_site := sum(C), by = PLT_CN]
  d[, .(LAI_necn = sum(M * C / (K + C_site))), by = PLT_CN]
}

#' How well a parameter set reproduces observed plot LAI under stand-level NECN, on plots where the LANDIS
#' species hold at least `min_share` of live wood (observed LAI also counts other species' foliage)
#' @return list(plots = per-plot table, stats = one-row summary, bins = medians by stand wood)
lai_check_fia <- function(plot_lai, sp_tab, params, min_share = 0.9) {
  sh <- sp_tab[, .(wood_landis = sum(wood_s)), by = PLT_CN]
  d <- merge(plot_lai$plots[, .(PLT_CN, obs = LAI, wood_gm2)], sh, by = "PLT_CN")
  d <- d[wood_landis / wood_gm2 >= min_share]
  d <- merge(d, stand_site_lai(sp_tab, params), by = "PLT_CN", all.x = TRUE)[is.na(LAI_necn), LAI_necn := 0]
  r2 <- 1 - sum((d$obs - d$LAI_necn)^2) / sum((d$obs - mean(d$obs))^2)
  d[, wood_Mgha := wood_gm2 / 100]
  br <- c(0, 50, 100, 200, 300, 400, 600, Inf)
  bins <- d[, .(plots = .N, observed = stats::median(obs), obs_q25 = stats::quantile(obs, .25), obs_q75 = stats::quantile(obs, .75),
                necn = stats::median(LAI_necn), wood = stats::median(wood_Mgha)), keyby = .(wood_class = cut(wood_Mgha, br, dig.lab = 4))]
  list(plots = d[], bins = bins[],
       stats = data.table::data.table(plots = nrow(d), r2 = r2, median_ratio = stats::median(d$LAI_necn / pmax(d$obs, 0.1)),
                                      bias = mean(d$LAI_necn - d$obs), rmse = sqrt(mean((d$LAI_necn - d$obs)^2))))
}
