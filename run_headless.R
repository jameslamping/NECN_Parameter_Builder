# Run the whole build without the app:  Rscript run_headless.R config.yml
args <- commandArgs(trailingOnly = TRUE)
dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
options(necnbuilder.dir = dir)
suppressMessages({library(data.table)})
for (f in list.files(file.path(dir, "R"), full.names = TRUE)) source(f)
run_builder(if (length(args)) args[1] else file.path(dir, "config_example.yml"))
