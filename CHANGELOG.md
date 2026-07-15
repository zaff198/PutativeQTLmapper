# Changelog

All notable changes to PutativeQTLmapper are documented here.

## v3.1 (current)

**Fixed**
- Interactive plot threw `subscript out of bounds` when rendered. Root cause: the
  interactive view was built by running `ggplotly()` on a `ggplot2` object that
  included a `geom_text_repel()` layer. Plotly's conversion code does not
  reliably translate ggrepel's internal segments/labels, and mismatched legend
  entries between layers is a known trigger for this exact error. Fixed by
  building the interactive view directly with plotly's native API
  (`plot_ly()` / `add_markers()` / `add_segments()`) instead of converting a
  ggplot object, which removes the incompatible conversion step entirely.
- A second, unrelated `subscript out of bounds` occurred when editing the data
  table. `datatable()` was rendered with its default `rownames = TRUE` (an
  extra row-number column shown in the UI) while `editData()` was called with
  `rownames = FALSE`, shifting column indices by one. Fixed by rendering with
  `rownames = FALSE` so both sides agree, plus added a bounds check on the
  edit handler so a bad edit shows a warning instead of crashing.

**Added**
- Explicit "Load Demo Data" and "Clear table" buttons, so the illustrative
  dataset can always be restored on demand rather than only appearing on
  first launch.
- Static preview panel (plain ggplot2, no plotly) shown alongside the
  interactive plot, so users always have a guaranteed-correct view that also
  matches the PNG export exactly.
- Developer attribution footer on every page and in the plot caption.

## v3.0

**Changed (scope reduction, by design)**
- Removed the genetic (cM) map entirely. cM positions require a linkage
  mapping population and a genetic map (TASSEL/JoinMap/MapQTL/R-qtl); a
  TASSEL GWAS association result is not that, so cM was never a value this
  tool could honestly plot.
- Removed the Manhattan plot. TASSEL already produces a genome-wide Manhattan
  plot from the *full* marker set tested. Reconstructing one here from only
  the significant hits would be visually similar but statistically
  misleading (it would look like a full scan when it isn't one).
- As a result, the app now requires a single input table (marker, chromosome,
  trait, model, p-value, R², Mb) instead of two separate uploads.

**Added**
- Self-explaining UI: a "How this works" landing tab, inline help text next
  to every input, and a dedicated "Methodology & Limitations" tab that
  documents the full pipeline (genotyping → TASSEL association test →
  significance filtering → physical positioning → visualization) and its
  honest limitations (marker density vs. resolution, no support interval,
  genome-build provenance, centromere positions being Nipponbare-specific,
  no re-validation of the user's significance calls, physical vs.
  recombination distance).
- Clear "DEMO DATA" badge shown whenever the illustrative dataset is active,
  so it can't be mistaken for real results.

## v2.0

**Fixed**
- Download handlers (`downloadHandler`) were defined *inside*
  `observeEvent(input$generate, ...)`, so the download buttons did nothing
  until after the Generate button had been clicked once. Moved to top-level
  reactives guarded with `req()`.
- Manual-entry tables were `editable = TRUE` but had no observer capturing
  cell edits back into the underlying reactive data, so edits were silently
  discarded. Added `editData()` handlers.
- No schema validation on upload; a mistyped column name produced a raw,
  unreadable ggplot2 error. Added explicit column-presence and type checks
  with user-facing messages.
- `p_value == 0` or `NA` produced `-log10(p) = Inf`, silently distorting
  point sizing. Values are now floored and flagged.
- `geom_segment(size = ...)` is deprecated in ggplot2 ≥ 3.4 (superseded by
  `linewidth`); updated accordingly.
- **Centromere positions were fabricated** as `Length_Mb * 0.4` on every
  chromosome and mislabeled "centromere" in the plot caption. Replaced with
  real per-chromosome centromere coordinates from the MSU / Rice Genome
  Annotation Project (RGAP), sourced from
  <https://rice.uga.edu/annotation_pseudo_centromeres.shtml>.
- No genome build, version, or citation traveled with the output. Every plot
  now carries a provenance caption (genome build, citation, centromere
  source, generation timestamp).

**Added**
- Genome-wide Manhattan plot with a configurable significance threshold
  (Bonferroni or custom -log10(p) cutoff). *(Later removed in v3.0 — see
  above — once the scope was narrowed to complement, not duplicate, TASSEL.)*

## v1.0

Initial script (`PutativeQTLmapper.R`): a Shiny app pairing user-supplied
TASSEL-like GWAS output (marker, chromosome, trait, model, p-value, R²,
genetic and physical positions) with a chromosome backbone plot, using
placeholder centromere markers and dummy demonstration data.
