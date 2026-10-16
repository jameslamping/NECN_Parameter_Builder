# sla.R -- specific leaf area (m2 hemi-surface leaf area per kg dry mass) for every FIA species on the plots
# Lookup order: user override > species value in the reference table > genus mean > generic by
# softwood/hardwood. SLA scales observed LAI (and so MaximumLAI) but does not affect KLAI.

GENERIC_SLA <- c(S = 5, H = 12)

load_sla_reference <- function(path = file.path(app_dir(), "data", "sla_reference.csv")) {
  s <- data.table::fread(path)
  s[, key := tolower(sub("^(\\S+ \\S+).*", "\\1", scientific_name))]
  s[, genus := sub(" .*", "", key)]
  s[]
}

#' SLA for a set of SPCDs
#' @param spcd integer vector of FIA species codes present
#' @param overrides optional data.table(SPCD, sla, source)
sla_lookup <- function(spcd, ref_species, sla_ref, overrides = NULL, generic = GENERIC_SLA) {
  r <- ref_species[SPCD %in% spcd, .(SPCD, SCIENTIFIC_NAME, GENUS, SFTWD_HRDWD, key)]
  sp <- sla_ref[!is.na(sla_m2_per_kg), .(sla = mean(sla_m2_per_kg), n = sum(as.numeric(n), na.rm = TRUE),
                                          source = paste(unique(source), collapse = "; ")), by = key]
  gn <- sla_ref[!is.na(sla_m2_per_kg), .(sla_g = mean(sla_m2_per_kg),
                                          src_g = paste0("genus mean of ", data.table::uniqueN(key), " species (", paste(unique(source), collapse = "; "), ")")), by = genus]
  r <- merge(r, sp, by = "key", all.x = TRUE)
  r[, genus := tolower(GENUS)]
  r <- merge(r, gn, by = "genus", all.x = TRUE)
  r[, basis := data.table::fifelse(!is.na(sla), "species", data.table::fifelse(!is.na(sla_g), "genus", "generic"))]
  r[basis == "genus", `:=`(sla = sla_g, source = src_g)]
  r[basis == "generic", `:=`(sla = unname(generic[data.table::fifelse(SFTWD_HRDWD == "S", "S", "H")]),
                             source = paste0("generic ", data.table::fifelse(SFTWD_HRDWD == "S", "softwood", "hardwood")))]
  if (!is.null(overrides) && nrow(overrides)) {
    ov <- overrides[!is.na(sla)]
    r[ov, on = "SPCD", `:=`(sla = i.sla, source = i.source, basis = "override")]
  }
  r[, .(SPCD, SCIENTIFIC_NAME, sla, basis, source, n)][order(SPCD)]
}
