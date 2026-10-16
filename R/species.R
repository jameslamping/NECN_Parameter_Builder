# species.R -- matching LANDIS species names to FIA species codes (SPCD)

load_ref_species <- function(path = file.path(app_dir(), "data", "ref_species.csv")) {
  r <- data.table::fread(path)
  r[, key := tolower(paste(GENUS, SPECIES))]
  r[]
}

#' "PseudotsugaMenziesii" -> "Pseudotsuga menziesii"; names that already contain a space are kept
landis_to_scientific <- function(x) {
  out <- sub("^([A-Z][a-z]+)([A-Z][a-z]+.*)$", "\\1 \\2", x)
  split <- grepl(" ", out) & !grepl(" ", x)
  out[split] <- paste(sub(" .*", "", out[split]), tolower(sub("^\\S+ ", "", out[split])))
  out
}

#' Read species names from a NECN species CSV (SpeciesCode column) or a LANDIS species .txt
read_landis_species <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    d <- data.table::fread(path)
    col <- intersect(c("SpeciesCode", "SpeciesName", "Name"), names(d))[1]
    if (is.na(col)) stop("CSV needs a SpeciesCode (or SpeciesName) column")
    return(unique(as.character(d[[col]])))
  }
  l <- readLines(path, warn = FALSE)
  l <- l[!grepl("^\\s*(>>|LandisData|$)", l)]
  unique(sub("\\s.*", "", trimws(l[grepl("^\\s*[A-Za-z]", l)])))
}

#' Match LANDIS names to SPCD. Returns one row per LANDIS species; SPCD may hold several codes
#' separated by commas (to lump FIA taxa into one LANDIS species).
match_species <- function(landis_names, ref) {
  data.table::rbindlist(lapply(landis_names, function(nm) {
    sci <- landis_to_scientific(nm)
    k <- tolower(sci)
    cand <- ref[key == k]
    status <- "exact"
    if (nrow(cand) == 0) {
      g <- sub(" .*", "", k); sp <- sub("^[^ ]+ ", "", k)
      cand <- ref[tolower(GENUS) == g & startsWith(tolower(SPECIES), substr(sp, 1, 4))]
      status <- if (nrow(cand)) "fuzzy - check" else "no match"
    }
    if (nrow(cand) > 1) {
      plain <- cand[is.na(VARIETY) | VARIETY == ""][is.na(SUBSPECIES) | SUBSPECIES == ""]
      others <- setdiff(cand$SPCD, plain$SPCD[1])
      pick <- if (nrow(plain)) plain[1] else cand[1]
      note <- paste(sprintf("%d %s", cand[SPCD %in% others]$SPCD, cand[SPCD %in% others]$SCIENTIFIC_NAME), collapse = "; ")
      return(data.table::data.table(LandisName = nm, SPCD = as.character(pick$SPCD), FIA_name = pick$SCIENTIFIC_NAME,
                                    common = pick$COMMON_NAME, status = "ambiguous - check", other_codes = note))
    }
    if (nrow(cand) == 0)
      return(data.table::data.table(LandisName = nm, SPCD = NA_character_, FIA_name = NA_character_, common = NA_character_,
                                    status = status, other_codes = ""))
    data.table::data.table(LandisName = nm, SPCD = as.character(cand$SPCD), FIA_name = cand$SCIENTIFIC_NAME,
                           common = cand$COMMON_NAME, status = status, other_codes = "")
  }))
}

#' Default functional group per LANDIS species from FIA's softwood/hardwood code (editable in the app)
add_default_groups <- function(map, ref) {
  m <- data.table::copy(map)
  first <- suppressWarnings(as.integer(sub("[,; ].*", "", m$SPCD)))
  sw <- ref$SFTWD_HRDWD[match(first, ref$SPCD)]
  m[, softwood := sw == "S"]
  if (!"group" %in% names(m)) m[, group := data.table::fifelse(is.na(sw), NA_character_, data.table::fifelse(sw == "S", "conifers", "hardwoods"))]
  m[]
}
species_groups <- function(map) map[, .(LandisName, group, softwood)]

#' Expand a species map (SPCD column may hold "741,747") to one row per SPCD
expand_species_map <- function(map) {
  m <- map[!is.na(SPCD) & nzchar(SPCD)]
  m[, .(SPCD = as.integer(trimws(strsplit(SPCD, "[,; ]+")[[1]]))), by = LandisName][!is.na(SPCD)]
}
