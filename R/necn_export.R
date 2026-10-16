# necn_export.R -- turning fitted values into an NECN species table

METHOD_LABEL <- c(mono = "monospecific plots", mixed = "mixed plots (pure-stand equivalent)", pooled = "group-pooled")

#' Final per-species values.
#' Order: user value override > borrowed species > chosen method (user choice, else recommended) > current table.
#' @param choices data.table(LandisName, method, MaxLAI_use, KLAI_use, borrow_from, note); may be NULL
final_values <- function(params, choices, necn_tab) {
  cur <- necn_tab[, .(LandisName = SpeciesCode, MaxLAI_current = as.numeric(MaximumLAI), KLAI_current = as.numeric(KLAI))]
  keep <- c("LandisName", "group", "n_any", "n_mono", "n_mixed", "confidence", "disagree", "check", "recommended", "why",
            "MaxLAI_mono", "KLAI_mono", "MaxLAI_mixed", "KLAI_mixed", "MaxLAI_pooled", "KLAI_pooled")
  f <- merge(cur, params[, intersect(keep, names(params)), with = FALSE], by = "LandisName", all.x = TRUE)
  empty <- data.table::data.table(LandisName = character(), method = character(), MaxLAI_use = numeric(),
                                  KLAI_use = numeric(), borrow_from = character(), note = character())
  if (is.null(choices) || !nrow(choices)) choices <- empty
  for (col in setdiff(names(empty), names(choices))) choices[, (col) := empty[[col]][NA_integer_]]
  f <- merge(f, choices[, names(empty), with = FALSE], by = "LandisName", all.x = TRUE)
  f[is.na(method) | !nzchar(method) | method == "auto", method := NA_character_]
  f[, method_used := data.table::fifelse(is.na(method), recommended, method)]
  pick <- function(m, a, b, c) data.table::fcase(m == "mono", a, m == "mixed", b, m == "pooled", c, default = NA_real_)
  f[, `:=`(MaxLAI_chosen = pick(method_used, MaxLAI_mono, MaxLAI_mixed, MaxLAI_pooled),
           KLAI_chosen   = pick(method_used, KLAI_mono, KLAI_mixed, KLAI_pooled))]
  f[, `:=`(MaximumLAI = MaxLAI_current, KLAI = KLAI_current, value_source = "current table")]
  f[!is.na(MaxLAI_chosen) & !is.na(KLAI_chosen), `:=`(MaximumLAI = MaxLAI_chosen, KLAI = KLAI_chosen,
                                                       value_source = paste0("FIA: ", METHOD_LABEL[method_used]))]
  for (i in which(!is.na(f$borrow_from) & nzchar(f$borrow_from))) {
    j <- match(f$borrow_from[i], f$LandisName)
    if (!is.na(j) && !is.na(f$KLAI_chosen[j]))
      f[i, `:=`(MaximumLAI = f$MaxLAI_chosen[j], KLAI = f$KLAI_chosen[j], value_source = paste("borrowed from", f$borrow_from[i]))]
  }
  f[!is.na(MaxLAI_use), `:=`(MaximumLAI = MaxLAI_use, value_source = "user override")]
  f[!is.na(KLAI_use), `:=`(KLAI = KLAI_use, value_source = "user override")]
  f[]
}

#' Growth-neutral LAI rescale for NECN v8 as coded (site LAI = sum of cohort LAI).
#' Multiplies every LAI by f and compensates every place NECN uses LAI, so growth,
#' competition and establishment are unchanged. Requires GrowthLAI / f <= 1.
necn_rescale <- function(tab, f = 0.47) {
  t <- data.table::copy(tab)
  g0 <- if ("GrowthLAI" %in% names(t)) t$GrowthLAI else rep(0.47, nrow(t))
  m0 <- if ("MinLAI" %in% names(t)) t$MinLAI else rep(0.1, nrow(t))
  if (any(g0 / f > 1 + 1e-9)) stop(sprintf("Rescale factor %.3f is below the smallest allowed (%.3f): NECN caps GrowthLAI at 1.", f, max(g0)))
  t[, MaximumLAI := round(MaximumLAI * f, 4)]
  t[, GrowthLAI := round(g0 / f, 4)]
  t[, MinLAI := round(m0 * f, 5)]
  if ("CompetitionIndex" %in% names(t)) t[, CompetitionIndex := round(CompetitionIndex / f, 5)]
  if ("LightLAIScale" %in% names(t)) t[, LightLAIScale := round(LightLAIScale * f, 4)]
  if ("LightLAIAdjust" %in% names(t)) t[, LightLAIAdjust := round(LightLAIAdjust * f, 4)]
  if ("LightLAILocation" %in% names(t)) {
    if (any(t$LightLAILocation / f > 1)) warning("LightLAILocation / f exceeds NECN's limit of 1 for some species; capped at 1.")
    t[, LightLAILocation := round(pmin(1, LightLAILocation / f), 5)]
  }
  t[]
}

#' Build the NECN species table to write
#' @param mode "current_code" (apply growth-neutral rescale with factor f), "stand_level" (values as fitted)
#' @param light optional output of light_final(); its *_used values replace the LightLAI columns
#'   (before the rescale, which then moves the curves onto the rescaled LAI scale with everything else)
build_necn_table <- function(necn_tab, final, mode = c("current_code", "stand_level"), f = 0.47, update_lai = TRUE, light = NULL) {
  mode <- match.arg(mode)
  t <- data.table::copy(necn_tab)
  for (col in intersect(c("MaximumLAI", "KLAI", "CompetitionIndex", "LightLAIScale", "LightLAIAdjust", "LightLAILocation"), names(t)))
    data.table::set(t, j = col, value = as.numeric(t[[col]]))
  if (update_lai) {
    i <- match(t$SpeciesCode, final$LandisName)
    t[, MaximumLAI := round(final$MaximumLAI[i], 3)]
    t[, KLAI := round(final$KLAI[i], 1)]
  }
  if (!is.null(light)) {
    i <- match(t$SpeciesCode, light$LandisName)
    t[, `:=`(LightLAIShape = round(light$Shape_used[i], 4), LightLAIScale = round(light$Scale_used[i], 4),
             LightLAILocation = round(light$Location_used[i], 5), LightLAIAdjust = round(light$Adjust_used[i], 4))]
  }
  if (any(t$MaximumLAI > 20, na.rm = TRUE)) warning("MaximumLAI above NECN's limit of 20 for some species")
  if (mode == "current_code") t <- necn_rescale(t, f)
  t[]
}
