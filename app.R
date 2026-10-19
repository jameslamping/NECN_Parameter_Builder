# NECN Parameter Builder -- estimate NECN LAI parameters (MaximumLAI, KLAI) from FIA data
# Run:  shiny::runApp("path/to/NECN_Parameter_Builder")

suppressPackageStartupMessages({
  library(shiny); library(bslib); library(DT); library(data.table); library(ggplot2)
})
options(necnbuilder.dir = normalizePath("."), shiny.maxRequestSize = 200 * 1024^2)
for (f in list.files("R", full.names = TRUE)) source(f)

REF <- load_ref_species()
SLA_REF <- load_sla_reference()
US_STATES <- c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID","IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS",
               "MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")
METHODS <- c("auto", "mono", "mixed", "pooled")

step_card <- function(title, ...) card(card_header(title), ..., full_screen = TRUE)
help_text <- function(...) div(class = "text-muted small mb-2", ...)

ui <- page_navbar(
  title = "NECN Parameter Builder",
  theme = bs_theme(version = 5, primary = "#2a78d6", "navbar-bg" = "#1f3b57"),
  fillable = FALSE,

  nav_panel("Start",
    layout_columns(col_widths = c(7, 5),
      step_card("What this does",
        p("Builds NECN LAI parameters (", code("MaximumLAI"), ", ", code("KLAI"), ") and establishment light curves (", code("LightLAI*"), ") for each LANDIS species from FIA plots, for any US region."),
        tags$ol(
          tags$li(strong("FIA data:"), " download state tables from the FIA DataMart, or point to a folder that has them; set plot filters."),
          tags$li(strong("Species:"), " load your NECN species table, check FIA species-code matches and functional groups."),
          tags$li(strong("Plot LAI:"), " observed LAI on every plot = live foliage biomass x specific leaf area (SLA)."),
          tags$li(strong("MaximumLAI & KLAI:"), " three estimates per species (monospecific plots, mixed plots, group-pooled) with a recommendation; review, override, and check the site LAI they give on FIA plots."),
          tags$li(strong("Establishment:"), " LightLAI curves (Shape, Scale, Location, Adjust) fitted to FIA seedling presence against site LAI on the scale NECN evaluates it."),
          tags$li(strong("Export:"), " NECN species table for the NECN you run (stand-level site LAI, or 8.3 and earlier with the growth-neutral rescale), parameter table, report, config."),
          tags$li(strong("Initial communities:"), " IC map and CSV from TreeMap 2022 on your ecoregion map's grid, with cohorts from live FIA trees (needs steps 1 and 2 only)."),
          tags$li(strong("Check a LANDIS run:"), " compare species regeneration in a finished (historical) run with FIA plots around the landscape: young biomass, frequency, composition, and regeneration along site LAI.")),
        help_text("KLAI is the stand wood carbon (g C/m2) at which LAI reaches half of MaximumLAI. It is read from LAI against wood, not age, and it does not depend on SLA.")),
      div(
      step_card("Which NECN will you run?",
        radioButtons("mode", NULL, selected = character(0), choices = c(
          "NECN with stand-level site LAI (the October 2026 change: each cohort adds its share of a stand-level curve)" = "stand_level",
          "NECN 8.3 and earlier (site LAI = sum of every cohort's own LAI)" = "current_code")),
        help_text(HTML("This sets how MaximumLAI / KLAI are checked, which LAI scale the light curves are fitted on, and what is exported. <b>Stand-level</b>: values are exported as fitted. <b>8.3 and earlier</b>: the export adds the growth-neutral rescale, so the same stand-level values work with the cohort sum. If you are unsure, check whether your NECN's <code>CohortBiomass.cs</code> has <code>SiteLAIContribution</code>."))),
      step_card("Load a saved configuration (optional)",
        fileInput("cfg_in", NULL, accept = c(".yml", ".yaml")),
        help_text("A config saved from the Export step restores every setting. Data still has to be loaded."),
        verbatimTextOutput("cfg_status"))))),

  nav_panel("1 · FIA data",
    layout_columns(col_widths = c(4, 8),
      div(
        step_card("FIA source",
          radioButtons("src", NULL, inline = TRUE, choices = c("Use a folder I already have" = "folder", "Download from FIA DataMart" = "download")),
          conditionalPanel("input.src == 'download'",
            selectizeInput("dl_states", "States to download", choices = US_STATES, multiple = TRUE, options = list(placeholder = "e.g. WA, OR, ID")),
            textInput("dl_dir", "Save to folder", file.path(path.expand("~"), "FIA_data")),
            radioButtons("dl_method", "Download with", inline = TRUE, choices = c("rFIA::getFIA" = "rFIA", "Direct (DataMart zip files)" = "direct")),
            help_text("Downloads TREE, PLOT, COND and SEEDLING for each state (TREE is 100–500 MB for large western states). Files already in the folder are skipped. rFIA falls back to direct download if it fails."),
            actionButton("dl_go", "Download", class = "btn-outline-primary w-100 mb-3")),
          textInput("fia_dir", "FIADB CSV folder (files like WA_TREE.csv)", ""),
          shinyFiles::shinyDirButton("fia_browse", "Browse…", "Select FIADB CSV folder", class = "btn-sm btn-outline-secondary mb-2"),
          uiOutput("states_ui")),
        step_card("Plot filters",
          sliderInput("forest_min", "Minimum forested share of plot", 0.5, 1, 0.9, 0.05),
          numericInput("min_year", "Earliest inventory year", 2000, 1990, 2030),
          checkboxInput("annual_only", "Annual inventory plots only", TRUE),
          checkboxInput("exclude_disturbed", "Drop plots with recorded disturbance or treatment", FALSE),
          fileInput("region_file", "Study-area polygon (optional: .gpkg, .geojson, or zipped shapefile)", accept = c(".gpkg", ".geojson", ".json", ".zip")),
          numericInput("buffer_km", "Buffer around study area (km)", 0, 0, 500),
          checkboxInput("use_cache", "Cache tables as parquet for faster reloads", TRUE),
          actionButton("load_fia", "Load FIA and select plots", class = "btn-primary w-100"))),
      div(
        layout_columns(fill = FALSE,
          value_box("States", textOutput("vb_states"), theme = "secondary"),
          value_box("Live trees read", textOutput("vb_trees"), theme = "secondary"),
          value_box("Plots selected", textOutput("vb_plots"), theme = "primary")),
        step_card("Selected plots", plotOutput("plot_map", height = "420px"))))),

  nav_panel("2 · Species",
    layout_columns(col_widths = c(4, 8),
      step_card("Your NECN species table",
        fileInput("necn_file", "NECN species parameter CSV (SpeciesCode column)", accept = ".csv"),
        textInput("necn_path", "…or a path to it on this computer", ""),
        actionButton("necn_path_load", "Load from path", class = "btn-sm btn-outline-secondary mb-3"),
        help_text("Species names like PseudotsugaMenziesii are matched to FIA species codes (SPCD). Check rows marked 'ambiguous' or 'no match'. To lump several FIA taxa, enter codes separated by commas, e.g. 741,747."),
        help_text("Group is used for group-pooled estimates. It defaults to conifers / hardwoods; finer groups (e.g. pines, shade-tolerant conifers, subalpine) work better when each has several well-sampled species."),
        actionButton("rematch", "Re-run automatic matching", class = "btn-sm btn-outline-secondary")),
      step_card("Matches (double-click SPCD or group to edit)", DTOutput("species_tbl")))),

  nav_panel("3 · Plot LAI",
    layout_columns(col_widths = c(4, 8),
      step_card("Specific leaf area and plot LAI",
        help_text("SLA (m2 hemi-surface area per kg dry mass) is looked up for every species on the plots: species value from the reference table, else genus mean, else a generic softwood/hardwood value. Overrides take priority."),
        help_text(HTML("The built-in SLA table is mostly <b>Pacific Northwest</b> (Oregon and northern California, Berner &amp; Law 2016). Elsewhere, check the 'SLA used' tab: species marked genus or generic need a value. Please add your species' values to <code>data/sla_reference.csv</code> on <a href='https://github.com/jameslamping/NECN_Parameter_Builder/blob/main/data/README.md' target='_blank'>GitHub</a> so everyone benefits.")),
        numericInput("sla_sw", "Generic SLA, softwoods", 5, 1, 30, 0.1),
        numericInput("sla_hw", "Generic SLA, hardwoods", 12, 1, 40, 0.1),
        actionButton("compute_lai", "Compute plot LAI", class = "btn-primary w-100")),
      navset_card_tab(full_screen = TRUE,
        nav_panel("SLA used", help_text("SLA for your species. Double-click the override column to edit, then click 'Compute plot LAI' again."), DTOutput("sla_tbl")),
        nav_panel("LAI vs wood", plotOutput("all_plots_fig", height = "460px")),
        nav_panel("Map", plotOutput("lai_map", height = "460px")),
        nav_panel("Coverage", help_text("Plots per species: anywhere, as a mixed-plot member (share ≥ step-4 mixed threshold) and monospecific (share ≥ step-4 dominance)."), DTOutput("coverage_tbl"))))),

  nav_panel("4 · MaximumLAI & KLAI",
    layout_columns(col_widths = c(3, 9),
      step_card("Fit settings",
        uiOutput("mode_badge4"),
        sliderInput("dom_min", "Monospecific: species share of live wood ≥", 0.6, 1, 0.8, 0.05),
        sliderInput("share_min", "Mixed plots: species share of live wood ≥", 0.1, 0.8, 0.3, 0.05),
        numericInput("n0", "Group pooling strength (plots)", 30, 0, 500, 5),
        help_text("A species' weight on its own monospecific estimate is n / (n + this)."),
        sliderInput("plateau_top", "Plateau = largest fraction of stands", 0.05, 0.5, 0.2, 0.05),
        sliderInput("plateau_q", "MaximumLAI = this LAI quantile of the plateau stands", 0.5, 0.99, 0.9, 0.01),
        numericInput("max_cap", "Cap MaximumLAI at (NECN limit 20)", 20, 1, 20),
        numericInput("min_n", "Minimum plots for a fit", 12, 3, 200),
        numericInput("nboot", "Bootstrap resamples", 200, 0, 2000, 50),
        actionButton("fit", "Fit", class = "btn-primary w-100")),
      div(
        step_card("Parameters",
          help_text(HTML("Three KLAI estimates per species: <b>mono</b> (monospecific plots), <b>mixed</b> (all plots where the species has at least the mixed share; its own LAI ÷ its wood share — runs high for minority hardwoods and subcanopy trees), <b>pooled</b> (mono shrunk toward the group value). <b>recommended</b>: mono with ≥ 20 plots; pooled with fewer; mixed for softwoods with no mono fit. <b>disagree</b>: mono and mixed differ by > 30%. <b>check</b>: the recommended value is > 1.5× from the group value. Double-click <b>method</b> (auto/mono/mixed/pooled), <b>MaxLAI_use</b>, <b>KLAI_use</b>, <b>borrow_from</b> or <b>note</b> to override.")),
          DTOutput("param_tbl")),
        navset_card_tab(full_screen = TRUE,
          nav_panel("Selected species", selectInput("sp_pick", NULL, choices = NULL, width = "100%"), plotOutput("sp_fig", height = "460px")),
          nav_panel("All species", uiOutput("all_fig_ui")),
          nav_panel("Groups", help_text("Group values are geometric means of the monospecific estimates of species with ≥ 20 monospecific plots (else all with an estimate)."), DTOutput("group_tbl")),
          nav_panel("LAI vs stand age", help_text("For intuition only: the age at which LAI rises reflects how fast wood accumulates. KLAI comes from the wood plot."),
                    uiOutput("age_fig_ui")),
          nav_panel("Check vs FIA", help_text("Site LAI that NECN with stand-level LAI computes on each FIA plot with the final MaximumLAI / KLAI below, against observed plot LAI. Plots where your LANDIS species hold at least 90% of live wood. For NECN 8.3 and earlier the exported values are rescaled, and site LAI there also depends on the number of cohorts (see step 5, LAI scales)."),
                    tableOutput("check_stats"), plotOutput("check_fig", height = "400px")))))),

  nav_panel("5 · Establishment",
    layout_columns(col_widths = c(3, 9),
      div(
        step_card("Establishment light curves",
          help_text(HTML("NECN lets a species establish with probability <code>min(1, Adjust &times; (WeibullPDF(site LAI; Shape, Scale) + Location))</code>, times seed arrival and climate. Here the curve is fitted to FIA seedling presence on plots where the species has a seed source.")),
          uiOutput("mode_badge5"),
          radioButtons("lc_scale", "Fit against site LAI", choices = c(
            "stand-level site LAI (NECN with the stand-level change)" = "stand_level",
            "sum of cohort LAI (NECN 8.3 and earlier)" = "current_code",
            "observed (foliage x SLA)" = "observed")),
          help_text("Follows your NECN choice on the Start tab: the curve must be fitted on the LAI scale your NECN evaluates it on. Stand-level LAI comes straight from each plot's wood by species with your step-4 MaximumLAI and KLAI. The 8.3 cohort sum rebuilds each plot's cohorts (species x 10-year age class); the export's rescale then moves the curves with everything else."),
          conditionalPanel("input.lc_scale == 'current_code'",
            textInput("age_dir", "Lu et al. (2025) tree ages folder (optional, <ST>_TREEAGE.csv)", ""),
            help_text("Without it, tree ages come from per-species models TOTAGE ~ DIA + HT fitted to FIA trees with a recorded age.")),
          radioButtons("seed_how", "Seed source", choices = c(
            "a live tree of the species on the plot" = "plot",
            "the species on any plot within a distance" = "radius",
            "none (all plots)" = "none")),
          conditionalPanel("input.seed_how != 'none'", numericInput("seed_dia", "Minimum seed-tree diameter (in)", 5, 1, 30, 1)),
          conditionalPanel("input.seed_how == 'radius'", numericInput("seed_km", "Distance (km)", 5, 1, 50, 1)),
          numericInput("lc_min_pres", "Minimum plots with seedlings for a fit", 20, 5, 500, 5),
          checkboxInput("lc_norm", "Scale the curve to peak at 1 (light only limits establishment)", TRUE),
          radioButtons("lc_bare", "Bare ground (site LAI 0, a cell with no cohorts)", choices = c(
            "at least the curve's value at LAI 0.5 (recommended)" = "lai05",
            "as fitted" = "fit")),
          help_text("NECN's Weibull curve is 0 at LAI exactly 0 whenever Shape > 1, so cleared cells would get no establishment. The recommended rule refits those species so bare ground is at least as good as LAI 0.5 (about 80% of full light): the fit chooses a curve that falls from bare ground, or a rising curve with a floor."),
          help_text("FIA presence also reflects climate, seed supply and how long seedlings persist, which NECN handles elsewhere, so by default only the shape is kept."),
          actionButton("fit_light", "Fit light curves", class = "btn-primary w-100")),
        uiOutput("light_status")),
      div(
        step_card("Light-curve parameters",
          help_text(HTML("<b>Shape, Scale, Location, Adjust</b>: fitted values on the chosen LAI scale (before any rescale). <b>dev_expl</b>: share of deviance explained by LAI versus a constant presence rate; <b>p_lai</b>: likelihood-ratio test against that constant. Species where LAI adds nothing (p &gt; 0.05) get a flat curve: light never limits them. <b>P_bare</b>: light probability on bare ground (site LAI 0, a cell with no cohorts) = Adjust &times; Location when Shape &gt; 1; it decides how fast burns and harvest units regenerate. <b>n_low_lai</b>: fitting plots with site LAI &lt; 0.5; under about 20, P_bare is an extrapolation (flagged). Double-click <b>use</b> (fit / current / flat), <b>borrow_from</b>, <b>P_bare_use</b> (set the bare-ground probability; Location is adjusted) or <b>note</b> to override.")),
          DTOutput("light_tbl")),
        navset_card_tab(full_screen = TRUE,
          nav_panel("Selected species", selectInput("lc_pick", NULL, choices = NULL, width = "100%"), plotOutput("lc_fig", height = "460px")),
          nav_panel("All species", uiOutput("lc_all_ui")),
          nav_panel("LAI scales", help_text("Observed LAI against the site LAI NECN computes from the same plots' cohorts with the step-4 values. The cohort sum runs well above observed LAI, which is why curves fitted on observed LAI do not transfer to NECN as coded."),
                    plotOutput("lai_scales_fig", height = "420px"), tableOutput("age_src_tbl")))))),

  nav_panel("6 · Export",
    layout_columns(col_widths = c(4, 8),
      step_card("Export",
        uiOutput("mode_badge6"),
        conditionalPanel("input.mode == 'current_code'",
          numericInput("rescale_f", "Rescale factor f", 0.47, 0.01, 1, 0.01),
          help_text("NECN 8.3 and earlier: the rescale multiplies every LAI by f and compensates GrowthLAI, MinLAI, CompetitionIndex and the LightLAI curves, so growth, competition and establishment are unchanged. f cannot go below the table's GrowthLAI (default 0.47) because NECN caps GrowthLAI at 1.")),
        conditionalPanel("input.mode == 'stand_level'",
          help_text("Stand-level NECN: MaximumLAI, KLAI and the LightLAI curves are written as fitted; GrowthLAI, MinLAI and CompetitionIndex are left as they are. Do not apply the old rescale. Next: recheck CompetitionIndex and retune MaximumANPP with single-cell runs, because site LAI and young-cohort growth change.")),
        checkboxInput("update_lai", "Replace MaximumLAI / KLAI with the values from step 4", TRUE),
        checkboxInput("update_light", "Replace LightLAI parameters with the values from step 5", TRUE),
        uiOutput("export_warn"),
        hr(),
        downloadButton("dl_necn", "NECN species table (.csv)", class = "w-100 mb-2"),
        downloadButton("dl_params", "Parameter table (.csv)", class = "w-100 mb-2"),
        downloadButton("dl_report", "Report (.html)", class = "w-100 mb-2"),
        downloadButton("dl_report_pdf", "Report (.pdf)", class = "w-100 mb-2"),
        downloadButton("dl_cfg", "Configuration (.yml)", class = "w-100")),
      step_card("Preview: changed columns", DTOutput("preview_tbl")))),

  nav_panel("7 · Initial communities",
    layout_columns(col_widths = c(3, 9),
      div(
        step_card("Inputs",
          help_text("Uses the FIA tables from step 1 (all inventory years, so TreeMap's plots are found) and the species matches from step 2."),
          textInput("ic_template", "Ecoregion map (.tif): sets extent, resolution, CRS and active cells (> 0)", ""),
          shinyFiles::shinyFilesButton("ic_template_browse", "Browse…", "Select ecoregion map", multiple = FALSE, class = "btn-sm btn-outline-secondary mb-2"),
          textInput("ic_treemap", "TreeMap 2022 raster (TreeMap_CONUS_2022.tif)", ""),
          shinyFiles::shinyFilesButton("ic_treemap_browse", "Browse…", "Select TreeMap raster", multiple = FALSE, class = "btn-sm btn-outline-secondary mb-2"),
          help_text("TreeMap 2022 is published by the USDA Forest Service (Research Data Archive). Its attribute table links each TM_ID to an FIA plot; if your copy lacks it, give the tree table below."),
          textInput("ic_tm_table", "TreeMap tree table CSV (only if the raster has no attribute table)", ""),
          textInput("ic_species_file", "LANDIS species file (.txt, optional: caps ages below Longevity)", ""),
          textInput("ic_age_dir", "Lu et al. (2025) tree ages folder (optional, same as step 5)", "")),
        step_card("Options",
          radioButtons("ic_resample", "TreeMap to your grid", choices = c("nearest 30 m pixel" = "near", "most common class in the cell (mode)" = "mode")),
          help_text("'Nearest' samples one 30 m pixel per cell, so each plot (and treeless ground) keeps its share of the landscape. 'Mode' takes the majority (treeless counts as a class), which smooths the map but drops scattered and minority stands."),
          radioButtons("ic_codes", "Map codes", choices = c("TreeMap ID (TM_ID)" = "tm_id", "renumber 2…N (1 = empty site)" = "renumber")),
          conditionalPanel("input.ic_codes == 'tm_id'", numericInput("ic_empty", "Map code for active cells with no trees", 1, 0, 1e6, 1),
                           help_text("Has no rows in the CSV, so those cells start empty. If it is also a TreeMap ID on your landscape, an unused code is chosen.")),
          numericInput("ic_min_bio", "Drop cohorts below (g/m2)", 1, 0, 100, 1),
          actionButton("ic_build", "Build initial communities", class = "btn-primary w-100"))),
      div(
        layout_columns(fill = FALSE,
          value_box("Active cells", textOutput("vb_ic_active"), theme = "secondary"),
          value_box("With trees", textOutput("vb_ic_trees"), theme = "secondary"),
          value_box("Plots found in FIA", textOutput("vb_ic_plots"), theme = "primary")),
        navset_card_tab(full_screen = TRUE,
          nav_panel("Map", plotOutput("ic_map_fig", height = "520px")),
          nav_panel("Summary", tableOutput("ic_summary"), tableOutput("ic_age_src")),
          nav_panel("Species", help_text("Landscape totals (biomass x cells x cell area) and biomass-weighted mean cohort age."), DTOutput("ic_species_tbl")),
          nav_panel("Check vs TreeMap", plotOutput("ic_check_fig", height = "460px")),
          nav_panel("Map codes", help_text("One row per map code: TreeMap ID, FIA plot, whether the plot was found in the loaded FIA tables, cells, and IC totals (g/m2). Plots not found start empty: load the state that holds them."), DTOutput("ic_xw_tbl"))),
        step_card("Save",
          layout_columns(
            downloadButton("dl_ic_map", "IC map (.tif)", class = "w-100"),
            downloadButton("dl_ic_csv", "IC CSV", class = "w-100"),
            downloadButton("dl_ic_xw", "Map-code crosswalk (.csv)", class = "w-100")),
          help_text("The report (HTML or PDF, step 6) includes an initial-communities section once this step has run."),
          layout_columns(col_widths = c(8, 4),
            textInput("ic_out", "…or write all three to a folder", ""),
            div(style = "margin-top: 32px", actionButton("ic_write", "Write", class = "btn-outline-primary w-100"))))))),

  nav_panel("8 · Check a LANDIS run",
    layout_columns(col_widths = c(3, 9),
      div(
        step_card("LANDIS run",
          help_text("Compares regeneration in a finished LANDIS-II NECN run with FIA plots around the landscape. Needs steps 1–2 (FIA tables with SEEDLING, species matches); step 3 adds the site-LAI views."),
          textInput("chk_dir", "Run folder (holds Scenario.txt, community-input-file-*.csv, output-community-0.tif)", ""),
          shinyFiles::shinyDirButton("chk_browse", "Browse…", "Select LANDIS run folder", class = "btn-sm btn-outline-secondary mb-2"),
          uiOutput("chk_run_info"),
          radioButtons("chk_mode", "How this run's NECN computes site LAI", choices = c(
            "stand-level (NECN with the stand-level change)" = "stand_level", "sum of cohort LAI (NECN 8.3 and earlier)" = "current_code")),
          numericInput("chk_age", "Regeneration = cohorts / trees up to this age (years)", 20, 5, 60, 5),
          help_text("FIA has biomass only for trees >= 1 in DBH, which are rarely under 10-15 years old, so use 20 or more. Set the Lu et al. tree ages folder (step 5) for realistic ages of small trees."),
          numericInput("chk_every", "Read every N years of the run", 5, 1, 20, 1),
          numericInput("chk_buffer", "FIA plots within this distance of the landscape (km)", 25, 0, 200, 5),
          help_text("Requires Output Biomass Community output (community-input-file-{t}.csv) and output-community-0.tif. Only cells forested at year 0 count. FIA tree ages come from the same age source as step 5/7 (Lu et al. folder or age models)."),
          actionButton("chk_go", "Read run and compare", class = "btn-primary w-100")),
        step_card("Comparison window",
          sliderInput("chk_window", "LANDIS years to compare", 0, 100, c(25, 100), step = 1),
          help_text("Default: years after the initial communities' young cohorts have aged out (> the regeneration age), so all regeneration counted was simulated."),
          downloadButton("dl_chk", "Comparison table (.csv)", class = "w-100 mb-2"),
          downloadButton("dl_chk_figs", "All figures (.zip of PNG)", class = "w-100 mb-2"),
          downloadButton("dl_chk_pdf", "Report (.pdf)", class = "w-100 mb-2"),
          downloadButton("dl_chk_html", "Report (.html)", class = "w-100"),
          help_text("The report includes every figure here plus anything from steps 4, 5 and 7."))),
      div(
        layout_columns(fill = FALSE,
          value_box("Forested cells (year 0)", textOutput("vb_chk_cells"), theme = "secondary"),
          value_box("FIA plots near the landscape", textOutput("vb_chk_plots"), theme = "secondary"),
          value_box("Years read", textOutput("vb_chk_years"), theme = "primary")),
        navset_card_tab(full_screen = TRUE,
          nav_panel("By species", uiOutput("chk_age_note"), uiOutput("chk_species_ui")),
          nav_panel("Table", help_text(HTML("<b>Young biomass</b> (leaf + wood of cohorts / trees up to the regeneration age, Mg/ha) is the main comparison: it does not depend on plot or cell size. <b>% plots / cells</b> with young trees: cells sample far more ground than FIA microplots, so cell values run higher for any species. <b>Where present</b> = young biomass ÷ share present: a gap here with similar frequencies points to young-cohort growth (MaximumANPP, KLAI) rather than establishment.")),
                    DTOutput("chk_tbl")),
          nav_panel("Through time", uiOutput("chk_time_ui")),
          nav_panel("Cohort biomass by age", help_text("Biomass of each species' cohorts by cohort age (species x 10-year age class on a cell or plot, built like the initial communities), LANDIS against FIA cohorts near the landscape. This is each species' growth curve: it tests MaximumANPP and KLAI (young-cohort growth) and, at older ages, mortality and longevity."),
                    uiOutput("chk_coh_ui"), uiOutput("chk_cohage_ui")),
          nav_panel("Composition", uiOutput("chk_comp_ui")),
          nav_panel("By site LAI", help_text("Regeneration rate in each site-LAI class: LANDIS cells with young cohorts against FIA plots with young trees and with seedlings. Compare shapes (cells are larger than FIA plots, so levels differ). Below: site LAI across the landscape through time against FIA plots."),
                    uiOutput("chk_lai_ui"), plotOutput("chk_fig_lsl", height = "300px")),
          nav_panel("Establishment log", help_text("New cohorts per year by reproduction type from NECN-reproduction-log.csv: seed, resprout, serotiny."), uiOutput("chk_log_ui")))))),

  nav_panel("Method",
    card(card_body(
      h4("Observed LAI"), p("Each FIA plot's LAI is the sum over live trees of foliage biomass (FIADB DRYBIO_FOLIAGE, NSVB) x trees per acre x SLA. Wood is DRYBIO_AG (stem, bark and branches). FIADB 9.1 or later is required."),
      h4("MaximumLAI and KLAI"), p("NECN computes cohort LAI as MaximumLAI x C / (KLAI + C), with C the cohort's wood carbon (wood biomass x 0.47). In a stand of one species this curve describes the whole stand. MaximumLAI is set from the upper plateau (a high quantile of LAI among the largest stands) and KLAI is fitted. KLAI is the wood carbon at which LAI reaches half of MaximumLAI; dividing by 0.47 and then 100 gives Mg/ha of wood."),
      h4("Three estimates"), tags$ul(
        tags$li(strong("Monospecific: "), "plots where the species holds at least the dominance share of live wood; y = plot LAI. Most direct, but some species rarely form pure stands."),
        tags$li(strong("Mixed plots: "), "every plot where the species holds at least the mixed share; y = the species' own LAI divided by its wood share (a pure-stand equivalent), weighted by share. Uses far more plots. Agrees with the monospecific estimate within about 10–30% for most conifers, but runs high for minority hardwoods and subcanopy trees, which carry more foliage per unit wood than canopy dominants."),
        tags$li(strong("Group-pooled: "), "the group value is the geometric mean of its species' monospecific estimates; each species' estimate is shrunk toward it with weight n / (n + n0). Species with no monospecific fit get the group value.")),
      h4("Which NECN"), p("NECN 8.3 and earlier sum every cohort's own LAI curve over all cohorts on a site. With 10-year age cohorts a site holds about 20 cohorts, so site LAI overshoots and depends on how a stand is split into cohorts. The stand-level change (October 2026) makes each cohort add MaximumLAI x C_cohort / (KLAI + C_site), its share of the species' curve at the stand's wood, so site LAI no longer depends on the cohort split. The builder's MaximumLAI and KLAI are stand-level values in both cases: for stand-level NECN they are exported as fitted; for 8.3 and earlier the export adds the growth-neutral rescale. A cohort's own growth limit (GrowthLAI) uses its own curve in both versions; re-tune MaximumANPP with single-cell runs, because KLAI also sets early cohort growth."),
      h4("Establishment light curves"), p("NECN v8 lets a species establish on a site with probability min(1, Adjust x (WeibullPDF(site LAI; Shape, Scale) + Location)), multiplied by seed arrival and the climate establishment probability. The builder fits this curve by maximum likelihood to presence/absence of FIA seedlings (SEEDLING table, microplots) of the species, using only plots where the species has a seed source. The LAI on the x-axis must be the LAI NECN will evaluate the curve at: for stand-level NECN, the stand-level site LAI from each plot's wood by species with your MaximumLAI and KLAI; for NECN 8.3 and earlier, the sum of cohort LAI, rebuilt from each plot's trees (age from Lu et al. 2025 or per-species TOTAGE ~ DIA + HT models, 10-year age classes). Check each species' light probability on bare ground (site LAI 0): it controls how fast disturbed cells regenerate, and with a seed-tree rule few fitting plots are near LAI 0. By default Adjust is rescaled so the curve peaks at 1 over the observed LAI range, because FIA presence also reflects climate, seed supply and seedling persistence."),
      h4("Initial communities"), p("The ecoregion map is the template. TreeMap 2022 (30 m; each pixel holds the ID of the FIA plot imputed to it) is cropped and resampled onto the template grid. Active cells with no TreeMap ID get an empty-site code with no rows in the CSV; inactive cells get 0. Each TreeMap ID's FIA plot supplies live trees (STATUSCD 1) of the LANDIS species, grouped into cohorts by species and 10-year age class: leaf biomass = DRYBIO_FOLIAGE and wood biomass = DRYBIO_AG (stem, bark, branches), both x TPA_UNADJ, in g/m2. Ages come from Lu et al. (2025) where available, else per-species TOTAGE ~ DIA + HT models, and are capped below species longevity when the LANDIS species file is given."),
      h4("SLA sources"), p("Built-in values: Berner & Law (2016) Scientific Data 3:160002, Table 2 (Oregon and northern California, converted from cm2 hemi-surface area per g C using measured leaf C); Bond-Lamberty et al. (2002) Tree Physiology 22:993 for aspen. Species without a value use the genus mean, then a generic value. SLA affects MaximumLAI, not KLAI."))))
)

server <- function(input, output, session) {
  rv <- reactiveValues(fia = NULL, plots = NULL, region = NULL, necn = NULL, map = NULL, sla = NULL,
                       sla_over = data.table(SPCD = integer(), sla = numeric(), source = character()),
                       lai = NULL, sp_tab = NULL, fit = NULL, choices = NULL, states_refresh = 0,
                       ic = NULL, cohorts = NULL, chk = NULL, lxt = NULL, lfit = NULL, lchoices = NULL, lfit_key = NULL)
  vols <- c(Home = path.expand("~"), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "fia_browse", roots = vols, session = session)
  observeEvent(input$fia_browse, {
    p <- shinyFiles::parseDirPath(vols, input$fia_browse)
    if (length(p)) updateTextInput(session, "fia_dir", value = as.character(p))
  })
  states_avail <- reactive({ rv$states_refresh; req(dir.exists(input$fia_dir)); fia_available_states(input$fia_dir) })
  output$states_ui <- renderUI({
    s <- tryCatch(states_avail(), error = function(e) character())
    if (!length(s)) return(help_text("Enter a folder that holds <STATE>_TREE.csv, _PLOT.csv and _COND.csv files."))
    checkboxGroupInput("states", sprintf("States to use (%d found)", length(s)), choices = s, selected = s, inline = TRUE)
  })

  # ---- download ----
  observeEvent(input$dl_go, {
    req(length(input$dl_states) > 0, nzchar(input$dl_dir))
    withProgress(message = "Downloading FIA", value = 0, {
      res <- tryCatch(fia_download(input$dl_states, input$dl_dir, input$dl_method, progress = function(m, v) setProgress(v * 0.95, detail = m)),
                      error = function(e) { showNotification(paste("Download failed:", conditionMessage(e)), type = "error", duration = NULL); NULL })
    })
    if (!is.null(res)) {
      updateTextInput(session, "fia_dir", value = input$dl_dir)
      rv$states_refresh <- rv$states_refresh + 1
      showNotification(sprintf("%d file(s) downloaded to %s", length(res), input$dl_dir), type = "message")
    }
  })

  # ---- config in ----
  output$cfg_status <- renderText("No configuration loaded.")
  observeEvent(input$cfg_in, {
    cfg <- yaml::read_yaml(input$cfg_in$datapath)
    if (!is.null(cfg$fia_dir)) updateTextInput(session, "fia_dir", value = cfg$fia_dir)
    for (id in c("forest_min", "dom_min", "share_min", "plateau_top", "plateau_q")) if (!is.null(cfg[[id]])) updateSliderInput(session, id, value = cfg[[id]])
    if (!is.null(cfg$light_scale)) updateRadioButtons(session, "lc_scale", selected = cfg$light_scale)
    if (!is.null(cfg$seed_source)) updateRadioButtons(session, "seed_how", selected = cfg$seed_source)
    if (!is.null(cfg$age_dir)) updateTextInput(session, "age_dir", value = cfg$age_dir)
    for (p in list(c("seed_min_dia", "seed_dia"), c("seed_radius_km", "seed_km"), c("light_min_pres", "lc_min_pres"))) if (!is.null(cfg[[p[1]]])) updateNumericInput(session, p[2], value = cfg[[p[1]]])
    for (p in list(c("ic_template", "ic_template"), c("treemap_raster", "ic_treemap"), c("treemap_table", "ic_tm_table"), c("ic_species_file", "ic_species_file")))
      if (!is.null(cfg[[p[1]]])) updateTextInput(session, p[2], value = cfg[[p[1]]])
    if (!is.null(cfg$ic_resample)) updateRadioButtons(session, "ic_resample", selected = cfg$ic_resample)
    if (!is.null(cfg$ic_codes)) updateRadioButtons(session, "ic_codes", selected = cfg$ic_codes)
    if (!is.null(cfg$ic_empty_code)) updateNumericInput(session, "ic_empty", value = cfg$ic_empty_code)
    if (!is.null(cfg$light_bare)) updateRadioButtons(session, "lc_bare", selected = cfg$light_bare)
    if (!is.null(cfg$light_normalize)) updateCheckboxInput(session, "lc_norm", value = cfg$light_normalize)
    for (id in c("min_year", "buffer_km", "max_cap", "min_n", "nboot", "n0", "rescale_f")) if (!is.null(cfg[[id]])) updateNumericInput(session, id, value = cfg[[id]])
    if (!is.null(cfg$sla_softwood)) updateNumericInput(session, "sla_sw", value = cfg$sla_softwood)
    if (!is.null(cfg$sla_hardwood)) updateNumericInput(session, "sla_hw", value = cfg$sla_hardwood)
    for (id in c("annual_only", "exclude_disturbed")) if (!is.null(cfg[[id]])) updateCheckboxInput(session, id, value = cfg[[id]])
    if (!is.null(cfg$mode)) updateRadioButtons(session, "mode", selected = cfg$mode)
    if (!is.null(cfg$states)) session$onFlushed(function() updateCheckboxGroupInput(session, "states", selected = cfg$states), once = TRUE)
    if (!is.null(cfg$necn_table) && file.exists(cfg$necn_table)) updateTextInput(session, "necn_path", value = cfg$necn_table)
    output$cfg_status <- renderText(paste("Loaded", input$cfg_in$name))
  })

  # ---- step 1 ----
  observeEvent(input$region_file, {
    f <- input$region_file$datapath
    if (grepl("\\.zip$", f)) { d <- tempfile(); utils::unzip(f, exdir = d); f <- list.files(d, "\\.shp$", recursive = TRUE, full.names = TRUE)[1] }
    rv$region <- tryCatch(sf::st_read(f, quiet = TRUE), error = function(e) { showNotification(e$message, type = "error"); NULL })
  })
  observeEvent(input$load_fia, {
    req(input$fia_dir, length(input$states) > 0)
    withProgress(message = "Loading FIA", value = 0, {
      rv$fia <- tryCatch(fia_read(input$fia_dir, input$states, if (input$use_cache) file.path(app_dir(), "cache") else NULL,
                                  progress = function(m, v) setProgress(v * 0.9, detail = m)),
                         error = function(e) { showNotification(e$message, type = "error", duration = NULL); NULL })
      req(rv$fia)
      setProgress(0.95, detail = "Selecting plots")
      rv$plots <- select_plots(rv$fia, input$forest_min, input$annual_only, input$min_year, input$exclude_disturbed, rv$region, input$buffer_km)
      rv$lai <- NULL; rv$fit <- NULL; rv$cohorts <- NULL; rv$ages <- NULL; rv$lfit <- NULL
      if (is.null(rv$fia$SEEDLING)) showNotification("No SEEDLING tables in this folder: step 5 (establishment curves) will need them.", type = "warning", duration = 8)
    })
  })
  output$vb_states <- renderText(if (is.null(rv$fia)) "–" else paste(input$states, collapse = " "))
  output$vb_trees <- renderText(if (is.null(rv$fia)) "–" else format(sum(rv$fia$TREE$STATUSCD == 1, na.rm = TRUE), big.mark = ","))
  output$vb_plots <- renderText(if (is.null(rv$plots)) "–" else format(nrow(rv$plots), big.mark = ","))
  output$plot_map <- renderPlot({
    req(rv$plots)
    g <- ggplot(rv$plots, aes(LON, LAT)) + geom_point(size = 0.3, alpha = 0.4, colour = "#1f3b57") + coord_quickmap() + theme_builder(11) + labs(x = NULL, y = NULL)
    if (!is.null(rv$region)) g <- g + geom_sf(data = sf::st_transform(rv$region, 4269), fill = NA, colour = "#eb6834", inherit.aes = FALSE)
    g
  })

  # ---- step 2 ----
  do_match <- function() { req(rv$necn); rv$map <- add_default_groups(match_species(rv$necn$SpeciesCode, REF), REF) }
  necn_name <- reactiveVal("")
  observeEvent(input$necn_file, { rv$necn <- fread(input$necn_file$datapath); necn_name(input$necn_file$name); do_match() })
  observeEvent(input$necn_path_load, {
    req(file.exists(input$necn_path)); rv$necn <- fread(input$necn_path); necn_name(input$necn_path); do_match()
  })
  observeEvent(input$rematch, do_match())
  species_view <- reactive({ req(rv$map); rv$map[, .(LandisName, SPCD, FIA_name, common, status, other_codes, softwood, group)] })
  output$species_tbl <- renderDT({
    v <- species_view()
    datatable(v, rownames = FALSE, selection = "none", options = list(pageLength = 50, dom = "ft", scrollX = TRUE),
              editable = list(target = "cell", disable = list(columns = c(0, 2, 3, 4, 5, 6)))) |>
      formatStyle("status", backgroundColor = styleEqual(c("exact", "ambiguous - check", "fuzzy - check", "no match", "user override"),
                                                         c("#e8f5e9", "#fff3cd", "#fff3cd", "#f8d7da", "#e3f2fd")))
  })
  observeEvent(input$species_tbl_cell_edit, {
    e <- input$species_tbl_cell_edit; m <- copy(rv$map); col <- names(species_view())[e$col + 1]
    if (col == "SPCD") {
      m[e$row, SPCD := as.character(e$value)]
      codes <- suppressWarnings(as.integer(strsplit(as.character(e$value), "[,; ]+")[[1]]))
      m[e$row, `:=`(FIA_name = paste(REF[SPCD %in% codes, SCIENTIFIC_NAME], collapse = " + "),
                    common = paste(REF[SPCD %in% codes, COMMON_NAME], collapse = " + "), status = "user override")]
      m[e$row, softwood := REF$SFTWD_HRDWD[match(codes[1], REF$SPCD)] == "S"]
      if (is.na(m$group[e$row])) m[e$row, group := data.table::fifelse(isTRUE(softwood), "conifers", "hardwoods")]
      rv$sp_tab <- NULL; rv$cohorts <- NULL
    } else if (col == "group") m[e$row, group := trimws(as.character(e$value))]
    rv$map <- m; rv$fit <- NULL
  })

  # ---- step 3 ----
  observeEvent(input$compute_lai, {
    req(rv$fia, rv$plots, rv$map)
    withProgress(message = "Computing plot LAI", value = 0.2, {
      spcd <- unique(rv$fia$TREE[PLT_CN %in% rv$plots$PLT_CN & STATUSCD == 1, SPCD])
      rv$sla <- sla_lookup(spcd, REF, SLA_REF, rv$sla_over, c(S = input$sla_sw, H = input$sla_hw))
      setProgress(0.5)
      rv$lai <- compute_plot_lai(rv$fia, rv$plots, rv$sla)
      rv$sp_tab <- species_plot_table(rv$lai, rv$map)
      rv$fit <- NULL
    })
  })
  sla_view <- reactive({
    req(rv$sla, rv$map)
    merge(expand_species_map(rv$map), rv$sla, by = "SPCD")[, .(LandisName, SPCD, SCIENTIFIC_NAME, sla = round(sla, 2), basis, source, override = NA_real_)]
  })
  output$sla_tbl <- renderDT(datatable(sla_view(), rownames = FALSE, selection = "none", options = list(pageLength = 30, dom = "t", scrollX = TRUE),
                                       editable = list(target = "cell", disable = list(columns = 0:5))))
  observeEvent(input$sla_tbl_cell_edit, {
    e <- input$sla_tbl_cell_edit; v <- sla_view()[e$row]
    val <- suppressWarnings(as.numeric(e$value))
    o <- rv$sla_over[SPCD != v$SPCD]
    if (!is.na(val)) o <- rbind(o, data.table(SPCD = v$SPCD, sla = val, source = "user override"))
    rv$sla_over <- o
    showNotification("Override stored. Click 'Compute plot LAI' to apply.", type = "message")
  })
  output$all_plots_fig <- renderPlot({ req(rv$lai); fig_all_plots(rv$lai$plots) })
  output$lai_map <- renderPlot({ req(rv$lai); fig_plot_map(rv$lai$plots) })
  output$coverage_tbl <- renderDT({
    req(rv$sp_tab)
    cv <- rv$sp_tab[, .(plots_any = .N, mixed = sum(share >= input$share_min), monospecific = sum(share >= input$dom_min)), by = LandisName]
    cv <- merge(rv$map[, .(LandisName, group)], cv, by = "LandisName", all.x = TRUE)
    datatable(cv[order(-monospecific)], rownames = FALSE, options = list(pageLength = 50, dom = "t"))
  })

  # ---- NECN version (Start tab) ----
  mode_label <- c(stand_level = "NECN with stand-level site LAI", current_code = "NECN 8.3 and earlier (cohort-sum LAI)")
  mode_badge <- function() renderUI({
    if (is.null(input$mode) || !nzchar(input$mode))
      return(div(class = "alert alert-warning small py-1 px-2", "Choose your NECN version on the Start tab."))
    div(class = "alert alert-info small py-1 px-2", strong("Parameterizing for: "), mode_label[[input$mode]])
  })
  output$mode_badge4 <- mode_badge(); output$mode_badge5 <- mode_badge(); output$mode_badge6 <- mode_badge()
  observeEvent(input$mode, if (nzchar(input$mode)) updateRadioButtons(session, "lc_scale", selected = input$mode), ignoreInit = FALSE)
  need_mode <- function() {
    if (is.null(input$mode) || !nzchar(input$mode)) { showNotification("Choose your NECN version on the Start tab first.", type = "warning"); return(FALSE) }
    TRUE
  }

  # ---- step 4 ----
  observeEvent(input$fit, {
    req(rv$sp_tab, rv$map)
    withProgress(message = "Fitting", value = 0, {
      rv$fit <- fit_lai_params(rv$sp_tab, species_groups(rv$map), input$dom_min, input$share_min, input$plateau_top, input$plateau_q,
                               input$max_cap, input$min_n, input$nboot, input$n0, progress = function(m, v) setProgress(v, detail = m))
    })
    rv$choices <- data.table(LandisName = rv$necn$SpeciesCode, method = "auto", MaxLAI_use = NA_real_, KLAI_use = NA_real_, borrow_from = "", note = "")
    ok <- rv$fit$params[recommended != "none"][order(-n_mono), LandisName]
    updateSelectInput(session, "sp_pick", choices = ok, selected = ok[1])
  })
  final <- reactive({ req(rv$fit, rv$necn); final_values(rv$fit$params, rv$choices, rv$necn) })
  param_view <- reactive({
    f <- final(); p <- rv$fit$params[, .(LandisName, KLAI_mono_lo, KLAI_mono_hi)]
    v <- merge(f, p, by = "LandisName", all.x = TRUE)
    m <- if (is.null(rv$choices)) rep("auto", nrow(v)) else rv$choices$method[match(v$LandisName, rv$choices$LandisName)]
    v[, .(LandisName, group, n_mono, n_mixed, confidence, recommended, why, disagree, check,
          MaxLAI_current, KLAI_current,
          KLAI_mono = round(KLAI_mono), KLAI_mono_80pct = ifelse(is.na(KLAI_mono_lo), "", sprintf("%.0f–%.0f", KLAI_mono_lo, KLAI_mono_hi)),
          KLAI_mixed = round(KLAI_mixed), KLAI_pooled = round(KLAI_pooled),
          method = ifelse(is.na(m) | !nzchar(m), "auto", m), MaxLAI_use, KLAI_use, borrow_from, note,
          MaximumLAI_final = round(MaximumLAI, 2), KLAI_final = round(KLAI), KLAI_final_Mgha_wood = round(KLAI / C_FRAC / 100), value_source,
          flag = fifelse(value_source == "current table", "current-table value (tuned for the cohort sum); borrow or override", ""))]
  })
  EDITABLE <- c("method", "MaxLAI_use", "KLAI_use", "borrow_from", "note")
  output$param_tbl <- renderDT({
    v <- param_view(); ed <- which(names(v) %in% EDITABLE) - 1
    datatable(v, rownames = FALSE, selection = "none", class = "compact nowrap stripe", options = list(pageLength = 50, dom = "ft", scrollX = TRUE),
              editable = list(target = "cell", disable = list(columns = setdiff(seq_along(v) - 1, ed)))) |>
      formatStyle("confidence", backgroundColor = styleEqual(c("good", "fair", "low", "none"), c("#e8f5e9", "#fff3cd", "#fde2cf", "#f8d7da"))) |>
      formatStyle(c("disagree", "check"), backgroundColor = styleEqual(c(TRUE, FALSE), c("#fff3cd", "transparent"))) |>
      formatStyle("method", backgroundColor = "#eef4fc")
  })
  observeEvent(input$param_tbl_cell_edit, {
    e <- input$param_tbl_cell_edit; v <- param_view(); col <- names(v)[e$col + 1]; sp <- v$LandisName[e$row]
    ch <- copy(rv$choices)
    if (!sp %in% ch$LandisName) ch <- rbind(ch, data.table(LandisName = sp, method = "auto", MaxLAI_use = NA_real_, KLAI_use = NA_real_, borrow_from = "", note = ""))
    val <- if (col %in% c("MaxLAI_use", "KLAI_use")) suppressWarnings(as.numeric(e$value)) else trimws(as.character(e$value))
    if (col == "method" && !val %in% METHODS) { showNotification("method must be auto, mono, mixed or pooled", type = "warning"); return() }
    if (col == "borrow_from" && nzchar(val) && !val %in% v$LandisName) showNotification(sprintf("'%s' is not a species in the table.", val), type = "warning")
    set(ch, which(ch$LandisName == sp), col, val)
    rv$choices <- ch
  })
  lai_check <- reactive({ req(rv$fit, rv$lai, rv$sp_tab); lai_check_fia(rv$lai, rv$sp_tab, final()[, .(LandisName, MaximumLAI, KLAI)]) })
  output$check_stats <- renderTable({
    st <- lai_check()$stats
    data.table(`FIA plots` = format(st$plots, big.mark = ","), R2 = sprintf("%.2f", st$r2), `Median NECN / observed` = sprintf("%.2f", st$median_ratio),
               `Mean difference` = sprintf("%+.2f", st$bias), RMSE = sprintf("%.2f", st$rmse))
  })
  output$check_fig <- renderPlot(fig_lai_check(lai_check()))
  output$group_tbl <- renderDT({
    req(rv$fit)
    datatable(rv$fit$groups[, .(group, n_species, species, MaxLAI_group = round(MaxLAI_group, 2), KLAI_group = round(KLAI_group))],
              rownames = FALSE, options = list(dom = "t"))
  })
  output$sp_fig <- renderPlot({ req(rv$fit, input$sp_pick); fig_species_fit(rv$fit, final(), input$sp_pick, ncol = 1, base = 13) })
  n_fit <- reactive({ req(rv$fit); sum(rv$fit$params$recommended != "none") })
  output$all_fig_ui <- renderUI(plotOutput("all_fig", height = paste0(max(320, 240 * ceiling(n_fit() / 4) + 60), "px")))
  output$all_fig <- renderPlot(fig_species_fit(rv$fit, final()))
  output$age_fig_ui <- renderUI(plotOutput("age_fig", height = paste0(max(320, 230 * ceiling(n_fit() / 4)), "px")))
  output$age_fig <- renderPlot(fig_species_age(rv$fit))

  # ---- tree ages (steps 5 and 7) ----
  age_models <- reactive({ req(rv$fia); fit_age_models(rv$fia$TREE[STATUSCD == 1], REF) })
  lu_ages <- reactive({ req(rv$fia); read_lu_ages(input$age_dir, input$states) })
  observeEvent(input$age_dir, if (!identical(input$age_dir, input$ic_age_dir)) updateTextInput(session, "ic_age_dir", value = input$age_dir))
  observeEvent(input$ic_age_dir, if (!identical(input$age_dir, input$ic_age_dir)) updateTextInput(session, "age_dir", value = input$ic_age_dir))

  # ---- step 5: establishment light curves ----
  lai_key <- reactive({ f <- final(); paste(round(f$MaximumLAI, 3), round(f$KLAI), collapse = ",") })
  observeEvent(input$fit_light, {
    req(rv$fia, rv$plots, rv$lai, rv$map, rv$necn)
    if (!need_mode()) return()
    if (input$lc_scale != "observed" && is.null(rv$fit)) { showNotification("Run step 4 first: the NECN LAI scales use your MaximumLAI and KLAI.", type = "warning"); return() }
    withProgress(message = "Establishment curves", value = 0, {
      ok <- tryCatch({
        pres <- seedling_presence(rv$fia, rv$plots$PLT_CN, rv$map)
        pl <- rv$lai$plots$PLT_CN
        treeless <- setdiff(rv$plots$PLT_CN, rv$fia$TREE[STATUSCD == 1 & TPA_UNADJ > 0, unique(PLT_CN)])
        xt <- rv$lai$plots[, .(PLT_CN, obs = LAI)]
        if (input$lc_scale != "observed") {
          f <- final()
          xt <- merge(xt, stand_site_lai(rv$sp_tab, f[, .(LandisName, MaximumLAI, KLAI)])[, .(PLT_CN, stand = LAI_necn)], by = "PLT_CN", all.x = TRUE)
          xt[is.na(stand), stand := 0]
        }
        if (input$lc_scale == "current_code") {
          if (is.null(rv$cohorts) || !identical(attr(rv$cohorts, "age_dir"), input$age_dir)) {
            setProgress(0.1, detail = "Tree ages")
            models <- age_models(); lu <- lu_ages()
            if (nzchar(input$age_dir) && is.null(lu)) showNotification("No <ST>_TREEAGE.csv files found for these states; using age models.", type = "warning")
            setProgress(0.25, detail = "Cohorts")
            co <- plot_cohorts(rv$fia, pl, rv$map, models, lu)
            attr(co, "age_dir") <- input$age_dir
            rv$cohorts <- co
          }
          minlai <- if ("MinLAI" %in% names(rv$necn)) as.numeric(rv$necn$MinLAI) else rep(0.1, nrow(rv$necn))
          pars <- merge(f[, .(LandisName, MaximumLAI, KLAI)], data.table(LandisName = rv$necn$SpeciesCode, MinLAI = minlai), by = "LandisName")
          xt <- merge(xt, necn_site_lai(rv$cohorts, pars, "current_code")[, .(PLT_CN, coded = LAI_necn, n_cohorts)], by = "PLT_CN", all.x = TRUE)
          xt[is.na(coded), coded := 0]
        }
        rv$lxt <- xt
        xcol <- c(current_code = "coded", stand_level = "stand", observed = "obs")[[input$lc_scale]]
        x_tab <- rbind(xt[, .(PLT_CN, x = get(xcol))], data.table(PLT_CN = treeless, x = 0))
        setProgress(0.4, detail = "Seed sources")
        src <- seed_source(rv$fia, rv$plots, rv$map, input$seed_how, input$seed_dia, input$seed_km)
        ld <- light_data(x_tab, pres, src)
        cols <- c("LightLAIShape", "LightLAIScale", "LightLAILocation", "LightLAIAdjust")
        cur <- if (all(cols %in% names(rv$necn))) rv$necn[, c("SpeciesCode", cols), with = FALSE][, LandisName := SpeciesCode][, SpeciesCode := NULL][] else NULL
        cmax <- if (input$lc_scale == "current_code") min(1, input$rescale_f) else 1
        lf <- fit_light_all(ld, NULL, input$lc_min_pres, normalize = input$lc_norm, c_max = cmax, bare = input$lc_bare,
                            progress = function(m, v) setProgress(0.45 + 0.5 * v, detail = m))
        lf$settings <- c(lf$settings, list(scale = input$lc_scale, seed_how = input$seed_how, seed_dia = input$seed_dia,
                                           seed_km = input$seed_km, age_dir = input$age_dir))
        rv$lfit <- lf; rv$lfit_key <- if (input$lc_scale == "observed") NULL else lai_key()
        rv$lchoices <- data.table(LandisName = character(), use = character(), borrow_from = character(), p0_use = numeric(), note = character())
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
    })
    if (ok) {
      sp <- rv$lfit$params[!is.na(Shape), LandisName]
      updateSelectInput(session, "lc_pick", choices = c(sp, setdiff(rv$lfit$params$LandisName, sp)), selected = sp[1])
    }
  })
  lfinal <- reactive({ req(rv$lfit, rv$necn); light_final(rv$lfit, rv$lchoices, rv$necn) })
  output$light_status <- renderUI({
    req(rv$lfit)
    msgs <- list()
    if (!is.null(rv$lfit_key) && !is.null(rv$fit) && !identical(rv$lfit_key, lai_key()))
      msgs <- c(msgs, list(div(class = "alert alert-warning small", "MaximumLAI or KLAI changed since these curves were fitted. Fit them again so the LAI scale matches.")))
    if (rv$lfit$settings$scale != "observed" && rv$lfit$settings$scale != input$mode)
      msgs <- c(msgs, list(div(class = "alert alert-warning small", sprintf("Curves were fitted on the '%s' scale but the export (step 6) is set to '%s'.", rv$lfit$settings$scale, input$mode))))
    tagList(msgs)
  })
  output$export_warn <- renderUI({
    if (is.null(rv$lfit) || !isTRUE(input$update_light)) return(NULL)
    sc <- rv$lfit$settings$scale
    if (sc == input$mode) return(NULL)
    div(class = "alert alert-warning small", if (sc == "observed")
      "Light curves were fitted on observed LAI. NECN evaluates them on its own site LAI, so they will be shifted; refit on the scale that matches this export."
      else sprintf("Light curves were fitted on the '%s' LAI scale; this export is '%s'. Refit in step 5.", sc, input$mode))
  })
  light_view <- reactive({
    v <- lfinal(); sc <- rv$lfit$settings$scale
    v[, flag := ""]
    v[status == "fitted" & !is.na(n_low_lai) & n_low_lai < 20 & is.na(p0_use), flag := sprintf("P_bare extrapolated (%d plots < LAI 0.5)", n_low_lai)]
    v[!is.na(p0_used) & p0_used < 0.05, flag := "P_bare ~ 0: no establishment on cleared cells; set P_bare_use or use the bare-ground rule"]
    if (sc != "current_code") v[light_source == "current table", flag := "current-table curve is on the old cohort-sum LAI scale: borrow, flat, or refit"]
    v[, .(LandisName, n_plots, n_present, pct_present = round(pct_present, 1), status, dev_expl = round(dev_expl, 3), p_lai = signif(p_lai, 2),
          n_low_lai, bare_rule, P_bare_fit = round(p_lai0, 3),
          Shape = round(Shape, 3), Scale = round(Scale, 3), Location = round(Location, 4), Adjust = round(Adjust, 4),
          use, borrow_from, P_bare_use = p0_use, note, flag, P_bare = round(p0_used, 3), P_LAI10 = round(p10_used, 3), p0_note,
          Shape_used = round(Shape_used, 3), Scale_used = round(Scale_used, 3),
          Location_used = round(Location_used, 4), Adjust_used = round(Adjust_used, 4), light_source,
          Shape_now = LightLAIShape, Scale_now = LightLAIScale, Location_now = LightLAILocation, Adjust_now = LightLAIAdjust)][order(-n_present)]
  })
  LEDIT <- c("use", "borrow_from", "P_bare_use", "note")
  output$light_tbl <- renderDT({
    v <- light_view(); ed <- which(names(v) %in% LEDIT) - 1
    datatable(v, rownames = FALSE, selection = "none", class = "compact nowrap stripe", options = list(pageLength = 50, dom = "ft", scrollX = TRUE),
              editable = list(target = "cell", disable = list(columns = setdiff(seq_along(v) - 1, ed)))) |>
      formatStyle("status", backgroundColor = styleEqual(c("fitted", "flat (no LAI effect)", "too few presences"), c("#e8f5e9", "#eef4fc", "#f8d7da"))) |>
      formatStyle(c("use", "borrow_from", "P_bare_use", "note"), backgroundColor = "#eef4fc") |>
      formatStyle("flag", color = "#9a3412")
  })
  observeEvent(input$light_tbl_cell_edit, {
    e <- input$light_tbl_cell_edit; v <- light_view(); col <- names(v)[e$col + 1]; sp <- v$LandisName[e$row]
    val <- trimws(as.character(e$value))
    if (col == "use" && !val %in% c("fit", "current", "flat", "")) { showNotification("use must be fit, current or flat", type = "warning"); return() }
    if (col == "P_bare_use") {
      col <- "p0_use"; val <- suppressWarnings(as.numeric(val))
      if (!is.na(val) && (val < 0 || val > 1)) { showNotification("P_bare_use must be between 0 and 1 (blank to clear)", type = "warning"); return() }
    }
    if (col == "borrow_from" && nzchar(val) && !val %in% v$LandisName) showNotification(sprintf("'%s' is not a species in the table.", val), type = "warning")
    ch <- copy(rv$lchoices)
    if (!sp %in% ch$LandisName) ch <- rbind(ch, data.table(LandisName = sp, use = NA_character_, borrow_from = "", p0_use = NA_real_, note = ""), fill = TRUE)
    set(ch, which(ch$LandisName == sp), col, if (col == "use" && !nzchar(val)) NA_character_ else val)
    rv$lchoices <- ch
  })
  output$lc_fig <- renderPlot({ req(rv$lfit, input$lc_pick); fig_light_fit(rv$lfit, lfinal(), input$lc_pick, rv$lfit$settings$scale, ncol = 1, base = 13) })
  n_lc <- reactive({ req(rv$lfit); nrow(rv$lfit$params) })
  output$lc_all_ui <- renderUI(plotOutput("lc_all", height = paste0(max(320, 230 * ceiling(n_lc() / 4) + 60), "px")))
  output$lc_all <- renderPlot(fig_light_fit(rv$lfit, lfinal(), NULL, rv$lfit$settings$scale))
  output$lai_scales_fig <- renderPlot({ req(rv$lxt, any(c("coded", "stand") %in% names(rv$lxt))); fig_lai_scales(rv$lxt) })
  output$age_src_tbl <- renderTable({
    req(rv$cohorts); a <- attr(rv$cohorts, "age_sources"); req(a)
    a[, .(`Tree age source` = age_src, Trees = format(N, big.mark = ","), Share = sprintf("%.0f%%", 100 * N / sum(N)))]
  })

  # ---- step 7: initial communities ----
  for (id in c("ic_template", "ic_treemap")) local({
    i <- id
    shinyFiles::shinyFileChoose(input, paste0(i, "_browse"), roots = vols, session = session, filetypes = c("tif", "tiff"))
    observeEvent(input[[paste0(i, "_browse")]], {
      p <- shinyFiles::parseFilePaths(vols, input[[paste0(i, "_browse")]])
      if (nrow(p)) updateTextInput(session, i, value = as.character(p$datapath))
    })
  })
  observeEvent(input$ic_build, {
    req(rv$fia, rv$map)
    if (!file.exists(input$ic_template)) { showNotification("Ecoregion map not found.", type = "error"); return() }
    if (!file.exists(input$ic_treemap)) { showNotification("TreeMap raster not found.", type = "error"); return() }
    withProgress(message = "Initial communities", value = 0, {
      rv$ic <- tryCatch(withCallingHandlers(
        build_ic(input$ic_template, input$ic_treemap, rv$fia, rv$map, age_models(), lu_ages(),
                 treemap_table = if (nzchar(input$ic_tm_table)) input$ic_tm_table else NULL,
                 longevity = read_longevity(input$ic_species_file), resample = input$ic_resample, codes = input$ic_codes,
                 empty_code = as.integer(input$ic_empty), min_biomass = input$ic_min_bio,
                 progress = function(m, v) setProgress(v, detail = m)),
        warning = function(w) { showNotification(conditionMessage(w), type = "warning", duration = 10); invokeRestart("muffleWarning") }),
        error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
    })
  })
  ic_val <- function(label) { req(rv$ic); s <- rv$ic$summary; v <- s$value[startsWith(s$item, label)][1]; format(v, big.mark = ",") }
  output$vb_ic_active <- renderText(if (is.null(rv$ic)) "–" else ic_val("Active cells (ecoregion"))
  output$vb_ic_trees <- renderText(if (is.null(rv$ic)) "–" else sprintf("%s cells", ic_val("Active cells with a TreeMap")))
  output$vb_ic_plots <- renderText(if (is.null(rv$ic)) "–" else sprintf("%s of %s", ic_val("Plots found"), ic_val("TreeMap IDs")))
  output$ic_map_fig <- renderPlot({ req(rv$ic); fig_ic_biomass(rv$ic) })
  output$ic_check_fig <- renderPlot({ req(rv$ic); fig_ic_check(rv$ic$crosswalk) })
  output$ic_summary <- renderTable({ req(rv$ic); rv$ic$summary[, .(Item = item, Value = format(value, big.mark = ","))] })
  output$ic_age_src <- renderTable({ req(rv$ic); a <- rv$ic$age_sources; a[, .(`Tree age source` = age_src, Trees = format(N, big.mark = ","), Share = sprintf("%.0f%%", 100 * N / sum(N)))] })
  output$ic_species_tbl <- renderDT({
    req(rv$ic)
    datatable(rv$ic$species[, .(SpeciesName, cells, landscape_Tg = round(landscape_Mg / 1e6, 3), cohorts, mean_age = round(mean_age), max_age)],
              rownames = FALSE, options = list(pageLength = 50, dom = "t"))
  })
  output$ic_xw_tbl <- renderDT({
    req(rv$ic)
    datatable(rv$ic$crosswalk[, .(MapCode, TM_ID, PLT_CN, FORTYPCD, in_fia, cells, n_species, n_cohorts, wood_gm2, leaf_gm2, treemap_live_gm2 = round(treemap_live_gm2))],
              rownames = FALSE, filter = "top", options = list(pageLength = 25, scrollX = TRUE)) |>
      formatStyle("in_fia", backgroundColor = styleEqual(c(TRUE, FALSE), c("transparent", "#f8d7da")))
  })
  ic_stem <- function() paste0("IC_", stamp())
  output$dl_ic_map <- downloadHandler(function() paste0(ic_stem(), "_map.tif"), function(file) {
    req(rv$ic); f <- tempfile(fileext = ".tif"); terra::writeRaster(rv$ic$map, f, datatype = "INT4S", NAflag = -9999); file.copy(f, file) })
  output$dl_ic_csv <- downloadHandler(function() paste0(ic_stem(), ".csv"), function(file) { req(rv$ic); fwrite(rv$ic$csv, file) })
  output$dl_ic_xw <- downloadHandler(function() paste0(ic_stem(), "_crosswalk.csv"), function(file) { req(rv$ic); fwrite(rv$ic$crosswalk, file) })
  observeEvent(input$ic_write, {
    req(rv$ic, nzchar(input$ic_out))
    f <- tryCatch(write_ic(rv$ic, input$ic_out, ic_stem()), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    if (!is.null(f)) showNotification(paste("Wrote", paste(basename(f), collapse = ", ")), type = "message", duration = 8)
  })

  # ---- step 8: check a LANDIS run ----
  shinyFiles::shinyDirChoose(input, "chk_browse", roots = vols, session = session)
  observeEvent(input$chk_browse, { p <- shinyFiles::parseDirPath(vols, input$chk_browse); if (length(p)) updateTextInput(session, "chk_dir", value = as.character(p)) })
  observeEvent(input$mode, if (nzchar(input$mode)) updateRadioButtons(session, "chk_mode", selected = input$mode))
  chk_run <- reactive({ req(nzchar(input$chk_dir)); tryCatch(read_landis_run(input$chk_dir), error = function(e) NULL) })
  output$chk_run_info <- renderUI({
    r <- chk_run(); if (is.null(r)) return(help_text("No LANDIS run found in this folder (needs Scenario.txt)."))
    help_text(sprintf("%d years; %d community files; %d LAI maps; species table: %s", r$duration, length(r$community_years), length(r$lai_years),
                      if (is.na(r$species_table)) "not found" else basename(r$species_table)))
  })
  observeEvent(list(chk_run(), input$chk_age), {
    r <- chk_run(); req(r); mx <- max(r$community_years)
    updateSliderInput(session, "chk_window", min = 0, max = mx, value = c(min(mx, input$chk_age + 1), mx))
  })
  observeEvent(input$chk_go, {
    req(rv$fia, rv$plots, rv$map); r <- chk_run()
    if (is.null(r)) { showNotification("Point to a LANDIS run folder first.", type = "warning"); return() }
    withProgress(message = "Checking LANDIS run", value = 0, {
      ok <- tryCatch({
        yrs <- unique(c(0, seq(0, max(r$community_years), by = input$chk_every), max(r$community_years)))
        lr <- landis_regen(r, yrs, input$chk_age, progress = function(m, v) setProgress(0.6 * v, detail = m))
        setProgress(0.65, detail = "FIA plots near the landscape")
        pars <- NULL
        if (!is.na(r$species_table) && file.exists(r$species_table)) {
          tb <- fread(r$species_table)
          pars <- tb[, .(LandisName = SpeciesCode, MaximumLAI = as.numeric(MaximumLAI), KLAI = as.numeric(KLAI),
                         MinLAI = if ("MinLAI" %in% names(tb)) as.numeric(MinLAI) else 0.1)]
          miss <- setdiff(unique(lr$regen$SpeciesName), rv$map$LandisName)
          if (length(miss)) showNotification(paste("Run species not in step 2's table (no FIA match):", paste(miss, collapse = ", ")), type = "warning", duration = 12)
        }
        if (input$chk_mode == "stand_level" && is.null(rv$sp_tab)) { pars <- NULL; showNotification("Run step 3 to add the site-LAI views.", type = "message") }
        fr <- fia_regen(rv$fia, rv$plots, landscape_outline(r), input$chk_buffer, rv$map, age_models(), lu_ages(), input$chk_age,
                        rv$sp_tab, pars, input$chk_mode, progress = function(m, v) setProgress(0.65 + 0.3 * v, detail = m))
        rv$chk <- list(run = r, lr = lr, fr = fr, rl = read_reproduction_log(r, lr$n_forest))
        if (!is.null(n <- fia_age_note(fr))) showNotification(n, type = "warning", duration = 30)
        TRUE
      }, error = function(e) { showNotification(conditionMessage(e), type = "error", duration = NULL); FALSE })
    })
  })
  chk_cmp <- reactive({ req(rv$chk); compare_regen(rv$chk$lr, rv$chk$fr, input$chk_window) })
  chk_species <- reactive({ req(rv$chk); landscape_species(rv$chk$lr, chk_cmp()) })
  chk_h <- function(per_row = 210) sprintf("%dpx", chk_fig_height(length(chk_species()), per_row))
  output$vb_chk_cells <- renderText(if (is.null(rv$chk)) "–" else format(rv$chk$lr$n_forest, big.mark = ","))
  output$vb_chk_plots <- renderText(if (is.null(rv$chk)) "–" else sprintf("%s (within %s km)", format(rv$chk$fr$n_plots, big.mark = ","), rv$chk$fr$buffer_km))
  output$vb_chk_years <- renderText(if (is.null(rv$chk)) "–" else sprintf("%d (%d–%d)", length(rv$chk$lr$years), min(rv$chk$lr$years), max(rv$chk$lr$years)))
  output$chk_species_ui <- renderUI(plotOutput("chk_fig_species", height = sprintf("%dpx", 180 + 34 * nrow(chk_cmp()))))
  output$chk_fig_species <- renderPlot(fig_regen_compare(chk_cmp()))
  output$chk_age_note <- renderUI({ req(rv$chk); n <- fia_age_note(rv$chk$fr); if (!is.null(n)) div(class = "alert alert-warning", n) })
  output$chk_tbl <- renderDT({
    x <- chk_cmp()[, .(Species = LandisName, assessment, FIA_Mgha = round(fia_Mgha, 3), FIA_90pct = sprintf("%.3f–%.3f", fia_lo, fia_hi),
                       LANDIS_Mgha = round(landis_Mgha, 3), LANDIS_range = sprintf("%.3f–%.3f", landis_min, landis_max), ratio = round(ratio, 2),
                       FIA_pct_plots = round(fia_pct_plots, 1), FIA_pct_seedlings = round(fia_pct_seedlings, 1), LANDIS_pct_cells = round(landis_pct_cells, 1),
                       FIA_where_present = round(fia_Mgha_where_present, 2), LANDIS_where_present = round(landis_Mgha_where_present, 2),
                       FIA_share = round(100 * fia_share, 1), LANDIS_share = round(100 * landis_share, 1))]
    datatable(x, rownames = FALSE, class = "compact nowrap stripe", options = list(pageLength = 50, dom = "ft", scrollX = TRUE)) |>
      formatStyle("assessment", backgroundColor = styleEqual(
        c("within the FIA range", "somewhat low", "somewhat high", "too little in LANDIS", "too much in LANDIS", "LANDIS regenerates; FIA has almost none", "little regeneration in either"),
        c("#e8f5e9", "#fff3cd", "#fff3cd", "#f8d7da", "#f8d7da", "#fde2cf", "transparent")))
  })
  output$chk_time_ui <- renderUI(plotOutput("chk_fig_time", height = chk_h()))
  output$chk_fig_time <- renderPlot(fig_regen_time(rv$chk$lr, chk_cmp(), chk_species()))
  output$chk_fig_comp <- renderPlot(fig_regen_composition(chk_cmp()))
  output$chk_comp_ui <- renderUI(plotOutput("chk_fig_comp", height = sprintf("%dpx", 140 + 30 * nrow(chk_cmp()))))
  chk_cohorts <- reactive({ req(rv$chk); cohort_biomass_by_age(rv$chk$lr, rv$chk$fr, input$chk_window) })
  output$chk_coh_ui <- renderUI(plotOutput("chk_fig_coh", height = chk_h()))
  output$chk_fig_coh <- renderPlot(fig_cohort_biomass(chk_cohorts(), chk_species()))
  output$chk_cohage_ui <- renderUI(plotOutput("chk_fig_cohage", height = chk_h(190)))
  output$chk_fig_cohage <- renderPlot(fig_cohort_ages(chk_cohorts(), chk_species()))
  output$chk_lai_ui <- renderUI(plotOutput("chk_fig_lai", height = chk_h()))
  output$chk_fig_lai <- renderPlot({ req(rv$chk); fig_regen_lai(rv$chk$lr, rv$chk$fr, input$chk_window, chk_species()) })
  output$chk_fig_lsl <- renderPlot({ req(rv$chk); fig_landscape_lai(rv$chk$lr, rv$chk$fr) })
  output$chk_log_ui <- renderUI(plotOutput("chk_fig_log", height = chk_h(190)))
  output$chk_fig_log <- renderPlot({ req(rv$chk$rl); fig_establish_log(rv$chk$rl, chk_species()) })
  output$dl_chk <- downloadHandler(function() sprintf("LANDIS_regeneration_check_%s.csv", stamp()), function(file) {
    x <- chk_cmp(); x[, `:=`(run = rv$chk$run$dir, window = paste(input$chk_window, collapse = "-"), regen_age = input$chk_age, fia_buffer_km = input$chk_buffer)]
    fwrite(x, file)
  })

  # ---- step 6: export ----
  observe({
    req(rv$necn)
    g0 <- if ("GrowthLAI" %in% names(rv$necn)) max(rv$necn$GrowthLAI) else 0.47
    updateNumericInput(session, "rescale_f", min = g0)
  })
  out_tab <- reactive({
    req(rv$necn, rv$fit, input$mode)
    tryCatch(build_necn_table(rv$necn, final(), input$mode, input$rescale_f, input$update_lai,
                              light = if (!is.null(rv$lfit) && isTRUE(input$update_light)) lfinal() else NULL),
             error = function(e) { showNotification(e$message, type = "error"); NULL })
  })
  output$preview_tbl <- renderDT({
    o <- out_tab(); req(o)
    cols <- intersect(c("MaximumLAI", "KLAI", "GrowthLAI", "MinLAI", "CompetitionIndex", "LightLAIShape", "LightLAIScale", "LightLAILocation", "LightLAIAdjust"), names(o))
    before <- rv$necn[, c("SpeciesCode", intersect(cols, names(rv$necn))), with = FALSE]
    setnames(before, setdiff(names(before), "SpeciesCode"), paste0(setdiff(names(before), "SpeciesCode"), "_now"))
    v <- merge(before, o[, c("SpeciesCode", cols), with = FALSE], by = "SpeciesCode")
    datatable(v, rownames = FALSE, options = list(pageLength = 50, dom = "t", scrollX = TRUE))
  })
  settings <- reactive(list(
    fia_dir = input$fia_dir, states = input$states, cache_dir = "cache", forest_min = input$forest_min,
    annual_only = input$annual_only, min_year = input$min_year, exclude_disturbed = input$exclude_disturbed,
    region_file = if (!is.null(input$region_file)) input$region_file$name else NULL, buffer_km = input$buffer_km,
    necn_table = necn_name(), species_file = necn_name(),
    sla_softwood = input$sla_sw, sla_hardwood = input$sla_hw, dom_min = input$dom_min, share_min = input$share_min, n0 = input$n0,
    plateau_top = input$plateau_top, plateau_q = input$plateau_q, max_cap = input$max_cap, min_n = input$min_n, nboot = input$nboot,
    seed = 1, mode = input$mode, rescale_f = input$rescale_f,
    light_fit = !is.null(rv$lfit), light_scale = input$lc_scale, age_dir = input$age_dir, seed_source = input$seed_how,
    seed_min_dia = input$seed_dia, seed_radius_km = input$seed_km, light_min_pres = input$lc_min_pres, light_normalize = input$lc_norm, light_bare = input$lc_bare,
    ic_build = !is.null(rv$ic), ic_template = input$ic_template, treemap_raster = input$ic_treemap, treemap_table = input$ic_tm_table,
    ic_species_file = input$ic_species_file, ic_resample = input$ic_resample, ic_codes = input$ic_codes, ic_empty_code = input$ic_empty,
    ic_min_biomass = input$ic_min_bio))
  stamp <- function() format(Sys.Date(), "%Y%m%d")
  output$dl_necn <- downloadHandler(function() sprintf("NECN_species_%s_%s.csv", input$mode, stamp()), function(file) fwrite(out_tab(), file))
  output$dl_params <- downloadHandler(function() sprintf("LAI_parameters_%s.csv", stamp()), function(file) {
    f <- final()
    if (!is.null(rv$lfit)) f <- merge(f, lfinal()[, .(LandisName, light_status = status, light_n_plots = n_plots, light_n_present = n_present,
                                                      LightLAIShape = Shape_used, LightLAIScale = Scale_used, LightLAILocation = Location_used,
                                                      LightLAIAdjust = Adjust_used, light_source, light_note = note)], by = "LandisName", all.x = TRUE)
    fwrite(f, file)
  })
  output$dl_cfg <- downloadHandler(function() sprintf("necn_builder_config_%s.yml", stamp()), function(file) yaml::write_yaml(settings(), file))
  report_input <- function() {
    if (is.null(rv$fit) && is.null(rv$ic) && is.null(rv$chk)) { showNotification("Nothing to report yet: run step 4, 7 or 8 first.", type = "warning"); req(FALSE) }
    if (!need_mode()) req(FALSE)
    list(settings = settings(), species_map = rv$map, sla_tab = rv$sla, plot_lai = rv$lai$plots, fit = rv$fit,
         check = if (!is.null(rv$fit)) lai_check() else NULL,
         final = if (!is.null(rv$fit)) final() else NULL, mode = input$mode, f = input$rescale_f,
         lfit = rv$lfit, lfinal = if (!is.null(rv$lfit)) lfinal() else NULL, ic = rv$ic,
         chk = if (!is.null(rv$chk)) c(rv$chk, list(cmp = chk_cmp(), cohorts = chk_cohorts(), window = input$chk_window)) else NULL)
  }
  html_report <- function(file) withProgress(message = "Writing report", value = 0.5, write_report(report_input(), file))
  pdf_report <- function(file) withProgress(message = "Writing PDF report", value = 0.3, {
    h <- tempfile(fileext = ".html"); write_report(report_input(), h)
    setProgress(0.7, detail = "Printing with Chrome")
    tryCatch(html_to_pdf(h, file), error = function(e) showNotification(conditionMessage(e), type = "error", duration = NULL))
  })
  output$dl_report <- downloadHandler(function() sprintf("NECN_builder_report_%s.html", stamp()), html_report)
  output$dl_report_pdf <- downloadHandler(function() sprintf("NECN_builder_report_%s.pdf", stamp()), pdf_report)
  output$dl_chk_html <- downloadHandler(function() sprintf("LANDIS_run_check_report_%s.html", stamp()), html_report)
  output$dl_chk_pdf <- downloadHandler(function() sprintf("LANDIS_run_check_report_%s.pdf", stamp()), pdf_report)
  output$dl_chk_figs <- downloadHandler(function() sprintf("LANDIS_run_check_figures_%s.zip", stamp()), function(file) {
    req(rv$chk)
    withProgress(message = "Writing figures", value = 0.5,
                 chk_figures_zip(c(rv$chk, list(cmp = chk_cmp(), cohorts = chk_cohorts(), window = input$chk_window)), file))
  })
}

shinyApp(ui, server)
