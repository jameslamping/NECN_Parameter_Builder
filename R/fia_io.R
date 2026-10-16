# fia_io.R -- reading FIADB CSV tables (one file per state and table, e.g. WA_TREE.csv)
# Only the columns the builder needs are read. Tables are cached as parquet (if arrow is
# installed) so later loads of the same files are fast.

FIA_COLS <- list(
  TREE = c("CN", "PLT_CN", "INVYR", "SPCD", "STATUSCD", "DIA", "HT", "TOTAGE", "TPA_UNADJ", "DRYBIO_AG", "DRYBIO_FOLIAGE"),
  PLOT = c("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR", "PLOT_STATUS_CD", "KINDCD",
           "LAT", "LON", "ELEV", "MEASYEAR"),
  COND = c("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "STDAGE", "FORTYPCD", "DSTRBCD1", "TRTCD1"),
  SEEDLING = c("PLT_CN", "INVYR", "SUBP", "CONDID", "SPCD", "TREECOUNT", "TREECOUNT_CALC")
)
CACHE_SCHEMA <- "s2"         # bump when FIA_COLS changes so old parquet caches are not reused
FIA_CHAR <- c("CN", "PLT_CN")

#' States that have a TREE, PLOT and COND file in a folder
fia_available_states <- function(dir) {
  f <- list.files(dir, pattern = "^[A-Z]{2}_(TREE|PLOT|COND)\\.csv$")
  st <- sub("_.*", "", f)
  tab <- table(st)
  sort(names(tab)[tab == 3])
}

.fia_read_one <- function(path, table, cache_dir) {
  cols_want <- FIA_COLS[[table]]
  key <- NULL
  if (!is.null(cache_dir) && requireNamespace("arrow", quietly = TRUE)) {
    info <- file.info(path)
    key <- file.path(cache_dir, sprintf("%s_%s_%s_%s.parquet", sub("\\.csv$", "", basename(path)),
                                        as.integer(info$mtime), info$size, CACHE_SCHEMA))
    if (file.exists(key)) return(data.table::as.data.table(arrow::read_parquet(key)))
  }
  hdr <- names(data.table::fread(path, nrows = 0))
  missing <- setdiff(cols_want, hdr)
  if (table == "TREE" && any(c("DRYBIO_FOLIAGE", "DRYBIO_AG") %in% missing))
    stop(sprintf("%s has no DRYBIO_FOLIAGE/DRYBIO_AG columns. The builder needs FIADB 9.1 or later (NSVB biomass).",
                 basename(path)))
  cols <- intersect(cols_want, hdr)
  cc <- intersect(FIA_CHAR, cols)
  d <- data.table::fread(path, select = cols, colClasses = list(character = cc), showProgress = FALSE)
  for (m in missing) d[, (m) := NA]
  if (!is.null(key)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    arrow::write_parquet(d, key)
  }
  d
}

#' Read TREE, PLOT and COND (and SEEDLING where present) tables for a set of states
#' @return list(TREE, PLOT, COND, SEEDLING) of data.tables; SEEDLING is NULL when no state has it,
#'   and attr(out, "seedling_states") lists the states it came from
fia_read <- function(dir, states, cache_dir = NULL, progress = function(msg, value) NULL) {
  stopifnot(length(states) > 0)
  out <- list()
  seed_st <- states[file.exists(file.path(dir, sprintf("%s_SEEDLING.csv", states)))]
  k <- 0; n <- length(states) * 3 + length(seed_st)
  for (tab in c("PLOT", "COND", "TREE", "SEEDLING")) {
    st <- if (tab == "SEEDLING") seed_st else states
    out[tab] <- list(if (length(st)) data.table::rbindlist(lapply(st, function(s) {
      k <<- k + 1
      progress(sprintf("Reading %s_%s", s, tab), k / n)
      .fia_read_one(file.path(dir, sprintf("%s_%s.csv", s, tab)), tab, cache_dir)
    }), use.names = TRUE, fill = TRUE))
  }
  attr(out, "seedling_states") <- seed_st
  out
}

#' Download TREE, PLOT, COND and SEEDLING CSVs for states from the FIA DataMart into `dir`
#' @param method "rFIA" (rFIA::getFIA) or "direct" (DataMart zip files); rFIA falls back to direct on error
fia_download <- function(states, dir, method = c("rFIA", "direct"), tables = c("TREE", "PLOT", "COND", "SEEDLING"),
                         progress = function(msg, value) NULL, overwrite = FALSE) {
  method <- match.arg(method)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
  todo <- data.table::CJ(state = toupper(states), table = tables)
  if (!overwrite) todo <- todo[!file.exists(file.path(dir, sprintf("%s_%s.csv", state, table)))]
  if (!nrow(todo)) return(invisible(character()))
  done <- character()
  for (st in unique(todo$state)) {
    tabs <- todo[state == st, table]
    progress(sprintf("Downloading %s (%s)", st, paste(tabs, collapse = ", ")), match(st, unique(todo$state)) / length(unique(todo$state)))
    ok <- FALSE
    if (method == "rFIA" && requireNamespace("rFIA", quietly = TRUE)) {
      ok <- tryCatch({ rFIA::getFIA(states = st, dir = dir, common = FALSE, tables = tabs, load = FALSE); TRUE },
                     error = function(e) { message("rFIA failed for ", st, ": ", conditionMessage(e), " -- trying direct download"); FALSE })
    }
    if (!ok) for (tb in tabs) {
      z <- tempfile(fileext = ".zip")
      utils::download.file(sprintf("https://apps.fs.usda.gov/fia/datamart/CSV/%s_%s.zip", st, tb), z, mode = "wb", quiet = TRUE)
      utils::unzip(z, exdir = dir); unlink(z)
    }
    done <- c(done, file.path(dir, sprintf("%s_%s.csv", st, tabs)))
  }
  missing <- done[!file.exists(done)]
  if (length(missing)) stop("Download finished but these files are missing: ", paste(basename(missing), collapse = ", "))
  invisible(done)
}
