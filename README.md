# PutativeQTLmapper

**An interactive R Shiny app for visualizing putative QTLs from SSR-based GWAS results on a real rice chromosome map.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![R](https://img.shields.io/badge/R-%E2%89%A54.1-276DC3?logo=r)
![Shiny](https://img.shields.io/badge/built%20with-Shiny-0b8a8f)
![Status](https://img.shields.io/badge/status-active-brightgreen)

Developed by **Dr. Zafir Ahmad Naik**, PhD Quantitative Genetics — [DataHarvest Labs](https://www.dataharvestlabs.com/)

---

## Table of Contents

- [What this tool is (and isn't)](#what-this-tool-is-and-isnt)
- [Screenshots](#screenshots)
- [Quick Start](#quick-start)
- [Input Data Schema](#input-data-schema)
- [Reference Genome & Provenance](#reference-genome--provenance)
- [Methodology](#methodology)
- [Limitations](#limitations)
- [Project Structure](#project-structure)
- [Roadmap / Changelog](#roadmap--changelog)
- [How to Cite](#how-to-cite)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## What this tool is (and isn't)

You've already run GWAS on SSR marker data in **TASSEL** (GLM/MLM), applied
your own significance threshold, and pulled out the markers associated with
your trait(s). **PutativeQTLmapper** takes that finished marker table and
places it onto a real rice chromosome ideogram — with genuine centromere
positions, not estimated ones — so you can see at a glance where your
putative QTLs sit relative to each other and to the centromere, and export a
publication-ready figure.

**It will:**
- Turn a spreadsheet of `Marker / Chr / Trait / Model / p_value / R² / Mb`
  into an interactive and static chromosome map.
- Size markers by `-log10(p)` and color by trait.
- Show real centromere positions per chromosome (see [Reference Genome &
  Provenance](#reference-genome--provenance)).
- Export a publication-quality PNG with full provenance in the caption.

**It will not:**
- Run GWAS or recompute significance — that's TASSEL's job, upstream of this
  tool.
- Generate a genome-wide Manhattan plot — that needs *all* markers tested,
  not just the significant subset this app expects; TASSEL already produces
  the correct one from your full dataset.
- Plot a genetic (cM) linkage map — cM positions require a linkage-mapping
  population and a genetic map (TASSEL/JoinMap/MapQTL/R-qtl output), which is
  a different kind of experiment from an association panel.
- Identify the causal gene, or tell you how precise your marker-to-QTL
  distance really is. See [Limitations](#limitations).

---

## Screenshots

*(Run the app locally — see [Quick Start](#quick-start) — to see the five
tabs: How This Works, Your Data, QTL Map, Summary Table, and Methodology &
Limitations. Screenshots will be added here in a future update.)*

---

## Quick Start

### Requirements

- R ≥ 4.1
- The following packages:

```r
install.packages(c(
  "shiny", "bslib", "DT", "ggplot2", "ggrepel",
  "dplyr", "scales", "readxl", "writexl",
  "plotly", "RColorBrewer"
))
```

### Run locally

```r
# From the repository root
shiny::runApp("app.R")
```

Or, if you use RStudio, just open `app.R` and click **Run App**.

### Try it with the bundled example data

The app loads a small, clearly-labeled demo dataset on launch (also available
at `example_data/example_qtl_data.csv`). Use the **Load Demo Data** button on
the "Your Data" tab any time you want to see it again, and **Clear table** to
start fresh before pasting in your own results.

### Bring your own data

1. Get physical (Mb) positions for your significant SSR markers — via
   Gramene, or by BLASTing each marker's flanking primer sequence against the
   reference genome assembly for a more precise coordinate (see
   [Methodology](#methodology)).
2. Format your marker table to match the [schema below](#input-data-schema).
3. Upload it on the "Your Data" tab (CSV or XLSX), or paste values directly
   into the editable table.
4. View and download your map on the "QTL Map" tab.

---

## Input Data Schema

| Column    | Type      | Description                                                        |
|-----------|-----------|--------------------------------------------------------------------|
| `Marker`  | text      | SSR / marker name (e.g. `RM285`)                                   |
| `Chr`     | integer   | Chromosome number, `1`–`12` (plain integer, not `"Chr09"`)          |
| `Trait`   | text      | Phenotype associated with the marker                               |
| `Model`   | text      | Statistical model used in TASSEL (e.g. `GLM`, `MLM`)                |
| `p_value` | numeric   | Association p-value from TASSEL                                    |
| `R2`      | numeric   | Proportion of phenotypic variance explained (R²), from TASSEL       |
| `Mb`      | numeric   | Physical position in megabases (Gramene lookup or BLAST result)     |

Only markers you've already determined to be significant belong in this
table — the app does not re-run multiple-testing correction.

---

## Reference Genome & Provenance

Chromosome lengths and centromere positions are **built into the app** — you
never need to supply or upload them.

- **Genome build:** Os-Nipponbare-Reference-IRGSP-1.0 pseudomolecules,
  MSU / Rice Genome Annotation Project (RGAP) Release 7
  (release date 31 Oct 2011; data files updated 1 Sep 2024).
  <https://rice.uga.edu/download_osa1r7.shtml>
- **Centromere positions (Mb):** identified via the CentO centromeric repeat
  (GenBank [AY101510](https://www.ncbi.nlm.nih.gov/nuccore/AY101510)) plus
  sequencing/FISH evidence.
  <https://rice.uga.edu/annotation_pseudo_centromeres.shtml>
- **Citation for the assembly/annotation:**
  Kawahara Y, de la Bastide M, Hamilton JP, et al. Improvement of the *Oryza
  sativa* Nipponbare reference genome using next generation sequence and
  optical map data. *Rice* 6, 4 (2013).
  [doi:10.1186/1939-8433-6-4](https://doi.org/10.1186/1939-8433-6-4)

If your marker positions originally came from an older MSU/RGAP build,
reconcile them first using the [Pseudomolecule Version
Converter](https://rice.uga.edu/analyses_search_converter.shtml) rather than
assuming coordinates carry over unchanged between releases.

---

## Methodology

1. **Genotyping** — SSR markers scored across your mapping/association panel.
2. **Association testing (TASSEL)** — GLM and/or MLM models test each marker
   against phenotype(s), producing a p-value and R² per marker–trait pair.
3. **Significance filtering** (done by you, in TASSEL, *before* using this
   app) — markers surviving your chosen threshold (Bonferroni, FDR, or a
   suggestive cutoff) are called putatively associated.
4. **Physical positioning** — each significant SSR is located on the
   reference genome via its cataloged Gramene position, ideally cross-checked
   or refined by BLASTing its primer sequence against the actual assembly,
   since Gramene positions can lag behind the current genome build.
5. **Visualization** (this app) — markers are placed on a chromosome
   ideogram at their Mb position, alongside the real centromere position for
   that chromosome, sized by `-log10(p)` and colored by trait.

---

## Limitations

- **This is a putative map, not a fine-mapped one.** SSR marker density is
  typically low, so the true causal gene could sit well away from a marker's
  exact position — this plot shows the marker, not a confidence interval
  around the QTL.
- **No support interval is drawn.** A single point per marker overstates
  precision; if you have LOD support intervals or LD-decay-based interval
  estimates, those are more honest to report alongside this figure.
- **Coordinate provenance matters.** Gramene-cataloged positions may
  reference an older genome build than the one used here. Reconcile with the
  Pseudomolecule Version Converter when in doubt.
- **Centromere positions are Nipponbare-specific.** They come from
  BAC-tiling/CentO-repeat evidence in the Nipponbare reference and may not
  exactly match centromere boundaries in the cultivars used in your own
  panel, due to structural variation between genotypes.
- **This app trusts your significance calls.** It does not re-run
  multiple-testing correction — unfiltered or non-significant markers pasted
  in will be plotted as if they were putative QTLs.
- **Co-localization can be coincidental**, especially with modest population
  size and sparse SSR coverage. Two traits mapping near the same marker is
  suggestive, not proof, of pleiotropy or linkage.
- **Physical distance ≠ recombination distance.** Mb position alone doesn't
  indicate how much recombination separates a marker from a nearby causal
  locus — that depends on local recombination rate, often suppressed near
  centromeres.

*In short: use this tool to communicate and explore where your significant
markers fall genome-wide — not as standalone evidence of a validated QTL.*

---

## Project Structure

```
PutativeQTLmapper/
├── app.R                          # Shiny app (single-file, standard entry point)
├── example_data/
│   └── example_qtl_data.csv       # Schema-matching example marker table
├── CHANGELOG.md                   # Version history and fixes
├── LICENSE                        # MIT License
└── README.md
```

---

## Roadmap / Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full version history, including
bug fixes (fabricated centromere placeholder replaced with real RGAP data,
broken download handlers, non-persisting table edits, a plotly/ggrepel
rendering conflict) and scope changes (genetic map and Manhattan plot
removed once the app's role was narrowed to complement, not duplicate,
TASSEL's own output).

---

## How to Cite

If this tool contributed to a publication or thesis, please cite both the
tool and the underlying genome resource it draws on:

> Naik, Z.A. *PutativeQTLmapper*: an interactive Shiny app for visualizing
> putative QTLs from SSR-based GWAS results. DataHarvest Labs.
> https://github.com/zaff198/PutativeQTLmapper

> Kawahara Y, de la Bastide M, Hamilton JP, et al. Improvement of the *Oryza
> sativa* Nipponbare reference genome using next generation sequence and
> optical map data. *Rice* 6, 4 (2013). doi:10.1186/1939-8433-6-4

---

## Contributing

Issues and pull requests are welcome — in particular, bug reports with a
minimal reproducible marker table are the fastest way to get something
fixed. Please open an issue before submitting a large structural change.

---

## License

Released under the [MIT License](LICENSE).

---

## Contact

**Dr. Zafir Ahmad Naik**
PhD Quantitative Genetics
DataHarvest Labs — [dataharvestlabs.com](https://www.dataharvestlabs.com/)
