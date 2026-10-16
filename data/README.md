# Specific leaf area table (`sla_reference.csv`)

Observed plot LAI in the builder is live-tree foliage biomass (FIA `DRYBIO_FOLIAGE`) × specific leaf area (SLA). SLA sets the height of the LAI curves, so it directly scales `MaximumLAI`. It does not affect `KLAI`.

## Coverage: mostly Pacific Northwest

The table started as the values for a Pacific Northwest project, so most of it covers **Oregon and northern California species**:

- 32 species from Berner & Law (2016), *Scientific Data* 3:160002, Table 2.
- Aspen from Bond-Lamberty, Wang & Gower (2002), *Tree Physiology* 22:993 (boreal Manitoba).

For other regions, many species will fall back to a genus mean, or to the generic softwood (5 m²/kg) or hardwood (12 m²/kg) value. The app's step 3 "SLA used" tab shows which basis each species got.

**Please add values for your species.** Every row makes the builder more accurate for the next person working in that region.

## How the table is used

For each FIA species on the plots, the builder uses the first of these that exists:

1. **Your override:** in the app (step 3, "SLA used" tab), or `sla_overrides` in the config.
2. **Species value:** the mean of every row in this table for that species. Several rows per species are fine, for example from different regions.
3. **Genus mean:** the mean of this table's species in the same genus.
4. **Generic value:** for softwood or hardwood (adjustable in step 3).

## Columns

| Column | What to enter |
|---|---|
| `scientific_name` | Genus and species as in FIA `REF_SPECIES` (`SCIENTIFIC_NAME`), e.g. `Pinus contorta`. Varieties and subspecies after the second word are ignored for matching. |
| `sla_m2_per_kg` | SLA in **m² of hemi-surface (half of total) leaf area per kg of dry leaf mass**. |
| `n` | Number of samples or trees behind the value (blank if unknown). |
| `source` | Full citation: authors, year, journal, volume, pages or table. |
| `notes` | Region, and any conversion you made (see below). |

## Units: converting a published value

- **cm²/g → m²/kg:** divide by 10.
- **Per gram of carbon (cm² / g C) → per gram of dry mass:** multiply by the leaf carbon fraction (e.g. 0.48), then divide by 10. This is how the Berner & Law values were converted; `notes` gives each conversion.
- **Area basis:**
  - For broad, flat leaves, one-sided projected area equals hemi-surface area.
  - For conifer needles, the two differ. Use hemi-surface area if the source gives it, and say which basis you used in `notes`.
  - Total (all-sided) area: divide by 2.
- **Leaf age:** prefer values for the whole canopy (all needle ages) over current-year foliage only. Say which in `notes`.

## Contributing on GitHub

- **Edit the table directly:** open `data/sla_reference.csv` on GitHub, click the pencil icon, add your rows, and open a pull request. Keep one species per row and fill in `source` and `notes`.
- **Or open an issue:** post the species, value, units and citation, and someone will add it.

Values from peer-reviewed studies or agency data are preferred. If a species already has a row from another region, add a new row rather than replacing it; the builder averages them.
