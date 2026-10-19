# NECN Parameter Builder

A Shiny app (plus a script-only runner) that estimates LANDIS-II NECN Succession LAI parameters
(`MaximumLAI`, `KLAI`) and establishment light curves (`LightLAIShape`, `LightLAIScale`, `LightLAILocation`,
`LightLAIAdjust`) for each species from FIA plots, and builds NECN initial communities (IC map + CSV) from TreeMap 2022. It works for any US region with FIADB data.

## Which NECN?

On the Start tab (or `mode:` in the config), choose the NECN you run. The choice sets the FIA check in step 4, the LAI scale the light curves are fitted on in step 5, and what is exported in step 6.

| Choice | NECN | Site LAI | Export |
|---|---|---|---|
| `stand_level` | NECN with the October 2026 stand-level change | each cohort adds `MaximumLAI × C_cohort / (KLAI + C_site)` | values as fitted |
| `current_code` | NECN 8.3 and earlier | sum of every cohort's own LAI curve | adds the growth-neutral rescale |

Your NECN has the change if its `CohortBiomass.cs` contains `SiteLAIContribution`.

## What it does

1. **FIA data.** Either downloads TREE, PLOT, COND and SEEDLING tables for chosen states from the FIA DataMart (via `rFIA::getFIA`, with direct DataMart download as a fallback), or uses a folder that already has them. Reads FIADB CSV files (`<ST>_TREE.csv`, `<ST>_PLOT.csv`, `<ST>_COND.csv`) for chosen states and selects plots with these filters:
   - sampled, annual inventory, from a minimum year onward;
   - at least 90% forested (adjustable);
   - one visit per plot location;
   - optionally, inside a study-area polygon plus a buffer.
2. **Species.** Matches LANDIS species names (for example `PseudotsugaMenziesii`) to FIA species codes and flags ambiguous matches. You can edit codes, or lump taxa with entries like `741,747`.
3. **Plot LAI.** Computes observed LAI per plot as live-tree foliage biomass (`DRYBIO_FOLIAGE`) × specific leaf area (SLA). SLA comes from the built-in reference table, then a genus mean, then a generic softwood/hardwood value; your overrides take priority.
4. **`MaximumLAI` and `KLAI`.** Three estimates per species. In each, `MaximumLAI` is the upper plateau (the 90th percentile of LAI in the largest 20% of stands, capped at 20), and `KLAI` is fitted on NECN's curve, `LAI = MaximumLAI × C / (KLAI + C)` with C = wood biomass × 0.47.
   - **mono:** plots where the species holds at least 80% of live wood; y = plot LAI.
   - **mixed:** every plot where it holds at least 30%; y = the species' own LAI divided by its wood share (a pure-stand equivalent), weighted by share. It uses many more plots and agrees with mono within about 10–30% for most conifers, but runs high for minority hardwoods and subcanopy trees.
   - **pooled:** the group value is the geometric mean of the group's species-level mono estimates. Each species is shrunk toward it with weight n / (n + n0). Groups default to conifers and hardwoods and can be edited in step 2.
   - **Recommended method:** mono with 20+ plots; pooled when mono has fewer; mixed for softwoods with no mono fit; otherwise the group value.
   - **Flags:** `disagree` when mono and mixed differ by more than 30%; `check` when the recommended value is more than 1.5× from the group value.
   - Per-species overrides: method, values, borrow from another species, notes.
   - **Check vs FIA:** the site LAI that stand-level NECN computes on each FIA plot with your final values, against observed plot LAI (R², median ratio, by stand wood). For NOCA species: R² 0.77, median ratio 1.01.
5. **Establishment light curves.** NECN lets a species establish with probability
   `min(1, LightLAIAdjust × (WeibullPDF(site LAI; Shape, Scale) + Location))`, times seed arrival and climate.
   - y = seedlings of the species on the plot's microplots (FIA SEEDLING table); only plots with a seed source count: a live tree of the species (≥ 5 in DBH by default) on the plot, or within a chosen distance, or no condition.
   - x = site LAI on the scale NECN will evaluate the curve on. For NECN v8 as coded, each plot's cohorts are rebuilt (species × 10-year age class; tree ages from Lu et al. 2025 where a folder of `<ST>_TREEAGE.csv` files is given, else per-species `TOTAGE ~ DIA + HT` models) and cohort LAI is summed with your MaximumLAI and KLAI. The cohort sum runs 1.5–3× above observed LAI, so curves fitted on observed LAI don't transfer.
   - Fitted by maximum likelihood (binomial). A fit is kept only if it beats a constant presence rate (likelihood-ratio test); otherwise the species gets a flat curve (light never limits it). Species with fewer than 20 plots with seedlings keep their current values unless you override them.
   - By default `Adjust` is scaled so the curve peaks at 1: FIA presence also reflects climate, seed supply and seedling persistence, which NECN models separately.
   - The export's growth-neutral rescale moves the curves with everything else.
   - **Bare ground.** A cell with no cohorts has site LAI exactly 0, where NECN's Weibull term is 0 whenever Shape > 1, so cleared cells would get no establishment.
     - The default rule (`light_bare: lai05`) refits those species so bare ground is at least as good as LAI 0.5.
     - The fit then chooses a curve that falls from bare ground, or a rising curve with a floor.
     - The table shows each species' bare-ground value (`P_bare`) and closed-canopy value (`P_LAI10`), and how many fitting plots sit below LAI 0.5.
     - `P_bare_use` refits a species with an exact bare-ground value.
     - `use = flat` makes light never limit a species.
   - **Flags:** bare-ground values that are extrapolated, values near 0, and species still on current-table curves (which were tuned for the old cohort-sum LAI).
6. **Export.**
   - The NECN species table, either for NECN v8 as coded (with the growth-neutral rescale applied) or for stand-level LAI.
   - A parameter table, a report (HTML, and PDF for sharing), and a `config.yml` that reproduces the run.
   - The PDF is printed by headless Google Chrome, Chromium or Microsoft Edge, so one of them must be installed. If it isn't found automatically, set the `NECN_CHROME` environment variable to its path. The report covers whichever steps were run, including initial communities.

7. **Initial communities.** Needs only steps 1 and 2 (FIA tables and species matches).
   - Your ecoregion map is the template: the IC map gets its extent, resolution, CRS and active cells (> 0).
   - TreeMap 2022 (30 m; each pixel holds the ID of an FIA plot) is cropped and resampled onto that grid. "Nearest pixel" is the default because it keeps each plot's share of the landscape; "mode" takes the majority, with treeless pixels counted as a class.
   - Active cells with no TreeMap ID get an empty-site code that has no CSV rows, so they start empty. The code defaults to 1, and if 1 is a real TreeMap ID on your landscape an unused code is picked instead. Inactive cells get 0.
   - TM_ID → FIA plot from the raster's attribute table (or the TreeMap tree table). Live trees (STATUSCD 1) of your LANDIS species become cohorts by species and 10-year age class. Leaf = `DRYBIO_FOLIAGE` and wood = `DRYBIO_AG`, both × `TPA_UNADJ`, in g/m².
   - Ages come from Lu et al. (2025) where given, else per-species `TOTAGE ~ DIA + HT` models. They are capped below species longevity when you give the LANDIS species file.
   - Outputs: IC map (INT4S GeoTIFF), IC CSV (`MapCode, SpeciesName, CohortAge, CohortBiomass, LeafBiomass, WoodBiomass`), and a map-code crosswalk (TM_ID, PLT_CN, whether the plot was found in FIA, cells, totals, TreeMap's own live biomass for checking).
   - TreeMap plots missing from the loaded FIA tables start empty and are listed. Load the states that hold them.

8. **Check a LANDIS run.** Compares species regeneration in a finished (e.g. historical) NECN run with FIA plots around the landscape.
   - **Regeneration** = cohorts (LANDIS) or live trees (FIA) up to a set age (default 20 years).
     - LANDIS: cells forested at year 0, in years after the initial communities' young cohorts have aged out.
     - FIA: plots within a set distance of the landscape. Trees are aged the same way the initial communities are, and FIA seedlings are reported too.
   - **Young biomass (Mg/ha)** is the main comparison, because it does not depend on plot or cell size. It is split into how often the species regenerates and how much young biomass it carries where present. A gap in the second, with similar frequencies, points to young-cohort growth (`MaximumANPP`, `KLAI`) rather than establishment.
   - **Other views:**
     - young biomass through time against the FIA range;
     - **cohort biomass by cohort age** per species: each species' growth curve, with cohorts built like the initial communities (species × 10-year class per cell or plot). It tests `MaximumANPP`/`KLAI` and, at older ages, mortality;
     - species composition;
     - regeneration rate by site LAI: the run's LAI maps against the site LAI the run's species table gives on FIA plots, which tests the light curves;
     - site LAI across the landscape through time;
     - the reproduction log by type (seed, resprout, serotiny).
   - **Needs** Output Biomass Community files (`community-input-file-{t}.csv`) and `output-community-0.tif`. In a script run, set `check_run_dir` in the config.

## Run

```r
# in R, from this folder
shiny::runApp(".")
```

Script only, without the app. Copy `config_example.yml` to `config_local.yml` (git ignores it), fill in your paths, then run it. Set `download: true` to fetch missing state tables first.

```bash
Rscript run_headless.R config_local.yml
```

## Requirements

- R ≥ 4.3 with these packages: shiny, bslib, DT, shinyFiles, data.table, ggplot2, minpack.lm, sf, yaml, htmltools, base64enc, scales. arrow is optional and is used for the parquet cache.
- FIADB 9.1 or later, which has the NSVB biomass columns `DRYBIO_FOLIAGE` and `DRYBIO_AG`.

## Notes on the method

- **KLAI** is the stand wood carbon (g C/m²) at which LAI reaches half of `MaximumLAI`. Divide by 0.47 and then by 100 to get Mg/ha of wood. It's estimated against wood, not stand age, and SLA doesn't affect it.
- **NECN 8.3 and earlier sum LAI over every cohort on a site**, so stand-level values used per cohort overshoot site LAI. For those versions the export applies a growth-neutral rescale:
  - multiply `MaximumLAI` by f and divide `GrowthLAI` by f;
  - multiply `MinLAI` by f and divide `CompetitionIndex` by f;
  - multiply `LightLAIScale` and `LightLAIAdjust` by f, and divide `LightLAILocation` by f.

  Growth, competition and establishment stay the same. f can't go below the table's `GrowthLAI`, because NECN caps `GrowthLAI` at 1.
- **`KLAI` also sets how fast young cohorts grow.** After changing it, re-tune `MaximumANPP` with single-cell runs.
- **Specific leaf area (SLA): the built-in table is Pacific Northwest focused.**
  - Sources: 32 species from Berner & Law (2016), *Scientific Data* 3:160002, Table 2 (Oregon and northern California), plus aspen from Bond-Lamberty et al. (2002), *Tree Physiology* 22:993.
  - Species elsewhere fall back to a genus mean or a generic value; step 3 shows which. SLA directly scales `MaximumLAI`, so a good value matters.
  - **Please contribute:** add your species' values to [`data/sla_reference.csv`](data/sla_reference.csv) with a pull request (or open an issue with the value and citation), so the table fills in for everyone. [`data/README.md`](data/README.md) explains the columns, units and conversions.

## Files

- `app.R`: the Shiny app.
- `R/`: core functions, shared by the app and `run_headless.R`:
  - `fia_io.R`, `plots.R`, `species.R`, `sla.R`, `lai_fit.R`, `cohorts.R` (tree ages, cohorts, NECN site LAI), `light_fit.R`, `ic_builder.R`, `necn_export.R`, `figures.R`, `report.R`, `pipeline.R`.
- `data/ref_species.csv`: FIA `REF_SPECIES` subset.
- `data/sla_reference.csv`: SLA reference table with sources (Pacific Northwest focused; contributions welcome).
- `data/README.md`: how the SLA table is used, its units, and how to add values.
- `config_example.yml`: template settings for `run_headless.R`.

## Planned (phase 2 and 3)

- Single-cell LANDIS setup generator, with Docker runs and comparison against FIA trajectories.
- Wiring FIA reference lines into the LANDIS output viewer.

## Contributing

- **SLA values:** the most useful contribution. See [`data/README.md`](data/README.md).
- **Bugs and suggestions:** open an issue with your NECN version, the states you used, and what you saw.
