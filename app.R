# =============================================================================
# PutativeQTLmapper v3.1
# Visualizes putative QTLs from your OWN TASSEL GWAS results (SSR markers)
# onto a physical (Mb) chromosome map with real centromere positions.
#
# Developed by Dr. Zafir Ahmad Naik, PhD Quantitative Genetics
# DataHarvest Labs -- https://www.dataharvestlabs.com/
#
# SCOPE:
#   - No genetic (cM) map. cM needs a linkage-mapping population; TASSEL-GWAS
#     results are association results, not linkage-map results.
#   - No Manhattan plot. TASSEL already produces this from your full marker
#     set; recreating it here from only your significant hits would be
#     misleading (a real Manhattan plot needs ALL tested markers).
#   - What this app does: take your list of significant SSR markers (already
#     filtered/thresholded by you in TASSEL) and lay them onto a real
#     chromosome ideogram with real centromeres, using the physical (Mb)
#     position you looked up for each marker (Gramene / BLAST).
#
# REFERENCE GENOME BUILD (Mb lengths + centromeres only -- built in):
#   Os-Nipponbare-Reference-IRGSP-1.0 pseudomolecules,
#   MSU / Rice Genome Annotation Project (RGAP) Release 7
#   (release date 31 Oct 2011; data files updated 1 Sep 2024)
#   https://rice.uga.edu/download_osa1r7.shtml
#
# CENTROMERE POSITIONS (Mb):
#   https://rice.uga.edu/annotation_pseudo_centromeres.shtml
#   Identified via the CentO centromeric repeat (GenBank AY101510;
#   Cheng et al. 2002, Plant Cell) plus sequencing/FISH evidence.
#
# CITATION for the assembly/annotation:
#   Kawahara Y, de la Bastide M, Hamilton JP, et al. Improvement of the
#   Oryza sativa Nipponbare reference genome using next generation
#   sequence and optical map data. Rice 6, 4 (2013).
#   https://doi.org/10.1186/1939-8433-6-4
#
# FIX LOG (v3 -> v3.1):
#   - "subscript out of bounds" on edit: datatable() was rendered with its
#     default rownames = TRUE (a row-number column shown in the UI), while
#     editData() was called with rownames = FALSE. That mismatch shifts the
#     column index by one, so an edited cell could reference a column past
#     the last real column. Fixed by rendering with rownames = FALSE so both
#     sides agree on column indexing.
#   - Added a dedicated "Load Demo Data" button (separate, explicit action)
#     so the illustrative dataset can always be brought back on demand.
#   - Added defensive bounds/NA checks around the edit handler so a bad edit
#     shows a clear message instead of a raw R error.
#   - Added developer attribution footer.
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(readxl)
library(writexl)
library(plotly)
library(RColorBrewer)

# ----------------------------- Reference data -------------------------------

GENOME_BUILD_LABEL <- "Os-Nipponbare-Reference-IRGSP-1.0 / MSU-RGAP Release 7 (2011; files updated 2024)"
GENOME_CITATION     <- "Kawahara et al. 2013, Rice 6:4, doi:10.1186/1939-8433-6-4"
CENTROMERE_SOURCE   <- "rice.uga.edu/annotation_pseudo_centromeres.shtml (CentO repeat; Cheng et al. 2002)"

rice_genome_ref <- data.frame(
  Chr           = 1:12,
  Length_Mb     = c(43.27, 35.94, 36.41, 35.50, 29.96, 31.25,
                     29.70, 28.44, 23.01, 23.21, 29.02, 27.53),
  Centromere_Mb = c(16.7, 13.6, 19.4, 9.7, 12.4, 15.3,
                     12.1, 12.9, 2.8, 8.2, 12.0, 11.9),
  stringsAsFactors = FALSE
)

# --------------------------- Demo / dummy data ------------------------------
# THIS IS NOT REAL DATA. It exists only so the app has something to show
# before you bring your own TASSEL output. Every value here is illustrative.
DEMO_QTL <- data.frame(
  Marker  = c("RM285", "RM27462", "RM240", "RM531", "RM101", "RM317"),
  Chr     = c(9, 12, 2, 6, 12, 3),
  Trait   = c("1000-seed weight", "Fertility (50F)", "Node/tiller no.",
              "Panicle length", "Plant height", "Grain length"),
  Model   = c("GLM", "GLM", "GLM", "MLM", "MLM", "MLM"),
  p_value = c(0.0035, 0.0050, 0.0083, 0.0021, 0.0087, 0.0044),
  R2      = c(0.056, 0.0518, 0.0433, 0.061, 0.0441, 0.0512),
  Mb      = c(6.8, 15.0, 29.6, 18.2, 9.0, 22.1),
  stringsAsFactors = FALSE
)

REQUIRED_COLS <- c("Marker", "Chr", "Trait", "Model", "p_value", "R2", "Mb")

# ----------------------------- Validation helpers ----------------------------

validate_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("No data loaded.")
  missing_cols <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing_cols) > 0) {
    return(paste0(
      "Your data is missing required column(s): ", paste(missing_cols, collapse = ", "),
      ". Expected columns: ", paste(REQUIRED_COLS, collapse = ", "), "."
    ))
  }
  bad_chr <- suppressWarnings(as.integer(df$Chr))
  if (any(is.na(bad_chr)) || any(!bad_chr %in% 1:12)) {
    return("Chr column must contain plain chromosome numbers 1-12 (e.g. 9, not 'Chr09' or 'chr9').")
  }
  bad_p <- suppressWarnings(as.numeric(df$p_value))
  if (any(is.na(bad_p))) {
    return("p_value column must be numeric for every row.")
  }
  bad_mb <- suppressWarnings(as.numeric(df$Mb))
  if (any(is.na(bad_mb))) {
    return("Mb column must be numeric for every row.")
  }
  NULL
}

clean_pvalues <- function(df) {
  df$p_value <- suppressWarnings(as.numeric(df$p_value))
  n_bad <- sum(is.na(df$p_value) | df$p_value <= 0)
  df$p_value <- pmax(df$p_value, 1e-300, na.rm = FALSE)
  df$p_value[is.na(df$p_value)] <- 1e-300
  attr(df, "n_invalid_pvalues") <- n_bad
  df
}

# ----------------------------- Plot builder ---------------------------------

theme_qtl <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      plot.title          = element_text(face = "bold", size = 17),
      plot.subtitle       = element_text(color = "gray40", size = 11),
      plot.caption        = element_text(color = "gray55", size = 8, hjust = 0),
      axis.title          = element_text(face = "bold"),
      legend.position     = "right"
    )
}

build_provenance_caption <- function() {
  paste0(
    "Genome build: ", GENOME_BUILD_LABEL, "\nCitation: ", GENOME_CITATION,
    "\nCentromere source: ", CENTROMERE_SOURCE,
    " | Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"),
    "\nApp by Dr. Zafir Ahmad Naik, DataHarvest Labs (dataharvestlabs.com)"
  )
}

create_qtl_map <- function(qtl_df, is_demo) {
  qtl_df <- qtl_df %>%
    mutate(
      Chr    = as.integer(as.character(Chr)),
      Mb     = as.numeric(Mb),
      R2     = suppressWarnings(as.numeric(R2)),
      logP   = -log10(p_value)
    )

  chr_levels <- sort(unique(rice_genome_ref$Chr))

  # Defensive: drop any rows whose Chr fell outside the reference set instead
  # of letting a stray value cause an out-of-range factor/index downstream.
  qtl_df <- qtl_df %>% filter(Chr %in% chr_levels)
  validate(need(nrow(qtl_df) > 0, "No valid rows to plot (check Chr values are 1-12)."))

  qtl_df$Chr <- factor(qtl_df$Chr, levels = chr_levels)

  chr_backbone <- rice_genome_ref %>%
    mutate(Chr = factor(Chr, levels = chr_levels))

  qtl_df <- qtl_df %>%
    mutate(Label = ifelse(
      duplicated(paste(Marker, Chr)) | duplicated(paste(Marker, Chr), fromLast = TRUE),
      paste0(Marker, " | ", Trait, " (", Model, ")"),
      paste0(Marker, " | ", Trait)
    ))

  max_pos <- max(chr_backbone$Length_Mb) * 1.08

  p <- ggplot() +
    geom_segment(data = chr_backbone,
                 aes(x = Chr, xend = Chr, y = 0, yend = Length_Mb),
                 linewidth = 6, color = "gray85", lineend = "round") +
    geom_point(data = chr_backbone, aes(x = Chr, y = Centromere_Mb),
               shape = 15, size = 4, color = "gray40") +
    geom_point(data = qtl_df,
               aes(x = Chr, y = Mb, size = logP, color = Trait,
                   text = paste0("Marker: ", Marker,
                                  "<br>Trait: ", Trait,
                                  "<br>Model: ", Model,
                                  "<br>p = ", signif(p_value, 3),
                                  "<br>R\u00b2 = ", R2,
                                  "<br>Position: ", Mb, " Mb")),
               alpha = 0.85, stroke = 1.1) +
    geom_text_repel(data = qtl_df,
                     aes(x = Chr, y = Mb, label = Label, color = Trait),
                     box.padding = 0.7, point.padding = 0.3,
                     segment.color = "gray50", size = 3.2, fontface = "bold",
                     max.overlaps = 30, show.legend = FALSE) +
    scale_x_discrete(name = "Chromosome", labels = paste0("Chr", chr_levels), drop = FALSE) +
    scale_y_continuous(name = "Physical Position (Mb)",
                        limits = c(-max_pos * 0.03, max_pos), expand = c(0, 0)) +
    scale_size_continuous(name = "-log10(p)") +
    scale_color_brewer(name = "Trait", palette = "Set2") +
    labs(
      title    = if (is_demo) "Putative QTL Map -- DEMO DATA" else "Putative QTL Map",
      subtitle = paste0(GENOME_BUILD_LABEL, "  |  Gray squares = real centromere positions"),
      caption  = build_provenance_caption()
    ) +
    theme_qtl()

  p
}

# NOTE ON THE INTERACTIVE VIEW:
# An earlier version built the interactive plot by running ggplotly() on the
# ggplot2/ggrepel object above. That is a known-fragile combination -- plotly's
# conversion code does not reliably handle the extra segments/labels that
# ggrepel adds, and mismatched legend entries between the point layer and the
# text layer is a common cause of exactly the "subscript out of bounds" error
# seen here. Rather than patch around that, the interactive tab below is built
# directly with plotly's own API (no ggplotly conversion at all), which avoids
# the incompatibility entirely. The ggplot2 version above is kept as-is and
# used only for the static preview and the PNG download, where ggrepel works
# fine because there is no plotly conversion involved.
create_qtl_map_plotly <- function(qtl_df, is_demo) {
  qtl_df <- qtl_df %>%
    mutate(
      Chr  = as.integer(as.character(Chr)),
      Mb   = as.numeric(Mb),
      R2   = suppressWarnings(as.numeric(R2)),
      logP = -log10(p_value)
    ) %>%
    filter(Chr %in% rice_genome_ref$Chr)

  validate(need(nrow(qtl_df) > 0, "No valid rows to plot (check Chr values are 1-12)."))

  chr_levels <- sort(unique(rice_genome_ref$Chr))
  max_pos <- max(rice_genome_ref$Length_Mb) * 1.08

  traits <- sort(unique(qtl_df$Trait))
  base_palette <- brewer.pal(max(3, min(8, length(traits))), "Set2")
  palette <- if (length(traits) > length(base_palette)) {
    colorRampPalette(base_palette)(length(traits))
  } else {
    base_palette[seq_along(traits)]
  }
  trait_colors <- setNames(palette, traits)

  p <- plot_ly()

  # Chromosome backbone, drawn once per chromosome.
  for (i in seq_len(nrow(rice_genome_ref))) {
    p <- p %>% add_segments(
      x = rice_genome_ref$Chr[i], xend = rice_genome_ref$Chr[i],
      y = 0, yend = rice_genome_ref$Length_Mb[i],
      line = list(color = "lightgray", width = 14),
      showlegend = FALSE, hoverinfo = "skip"
    )
  }

  # Real centromere positions.
  p <- p %>% add_markers(
    data = rice_genome_ref, x = ~Chr, y = ~Centromere_Mb,
    marker = list(symbol = "square", size = 10, color = "gray40"),
    name = "Centromere", hoverinfo = "text",
    text = ~paste0("Chr", Chr, " centromere: ", Centromere_Mb, " Mb"),
    showlegend = TRUE
  )

  # QTL markers, one trace per trait so the legend is genuinely toggleable.
  for (tr in traits) {
    sub <- qtl_df %>% filter(Trait == tr)
    p <- p %>% add_markers(
      data = sub, x = ~Chr, y = ~Mb,
      marker = list(size = ~(logP * 6 + 8), color = trait_colors[[tr]],
                    line = list(width = 1, color = "white")),
      name = tr, hoverinfo = "text",
      text = ~paste0("Marker: ", Marker,
                      "<br>Trait: ", Trait,
                      "<br>Model: ", Model,
                      "<br>p = ", signif(p_value, 3),
                      "<br>R2 = ", R2,
                      "<br>Position: ", Mb, " Mb")
    )
  }

  p %>% layout(
    title  = list(text = if (is_demo) "Putative QTL Map -- DEMO DATA" else "Putative QTL Map"),
    xaxis  = list(title = "Chromosome", tickvals = chr_levels,
                  ticktext = paste0("Chr", chr_levels), range = c(0.3, 12.7)),
    yaxis  = list(title = "Physical Position (Mb)", range = c(-max_pos * 0.03, max_pos)),
    legend = list(title = list(text = "Trait / Feature"))
  )
}

# ================================ UI ========================================

app_footer <- tags$footer(
  class = "text-center text-muted py-3 mt-4",
  style = "font-size:0.85em; border-top:1px solid #ddd;",
  "Developed by ", tags$b("Dr. Zafir Ahmad Naik"), ", PhD Quantitative Genetics, ",
  tags$a(href = "https://www.dataharvestlabs.com/", target = "_blank", "DataHarvest Labs"),
  " (dataharvestlabs.com)"
)

ui <- page_navbar(
  title = "PutativeQTLmapper",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2c5f2d"),
  fillable = FALSE,
  footer = app_footer,

  nav_panel(
    "1. How this works",
    layout_column_wrap(
      width = 1/2, heights_equal = "row",
      card(
        card_header("What this app does"),
        card_body(
          p("You already ran GWAS on your SSR marker data in TASSEL (GLM/MLM) ",
            "and pulled out the markers that were significant for your trait(s). ",
            "This app does not repeat that analysis -- it takes your finished, ",
            "significant-marker table and draws it onto a real rice chromosome ",
            "map, so you can see where your putative QTLs sit relative to each ",
            "other and to the centromere."),
          p(strong("It will NOT:"), " recompute significance, run GWAS, generate ",
            "a genome-wide Manhattan plot (TASSEL already gives you that from ",
            "your full marker set), or tell you the true causal gene."),
          p(strong("It WILL:"), " turn a spreadsheet of marker/Chr/p-value/R\u00b2/Mb ",
            "into a publication-style ideogram with correct centromere positions.")
        )
      ),
      card(
        card_header("Three steps"),
        card_body(
          tags$ol(
            tags$li(strong("Get physical (Mb) positions"), " for your significant SSR markers: ",
                    "look up each marker in Gramene, or BLAST its flanking primer sequence ",
                    "against the reference genome assembly for a more precise coordinate."),
            tags$li(strong("Fill in the table"), " on the 'Your Data' tab -- either edit the ",
                    "built-in example directly, or upload your own CSV/XLSX with the same columns."),
            tags$li(strong("View and download"), " your putative QTL map on the 'QTL Map' tab.")
          ),
          p(em("Tip: "), "the demo data currently loaded is fake -- replace it before drawing ",
            "any conclusions or sharing a figure. Use ", strong("Load Demo Data"), " any time ",
            "you want to see the example again.")
        )
      )
    )
  ),

  nav_panel(
    "2. Your Data",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Upload your TASSEL-derived marker table"),
        fileInput("qtl_file", NULL, accept = c(".csv", ".xlsx")),
        helpText("Required columns:"),
        tags$ul(
          tags$li(tags$b("Marker"), " -- SSR/marker name"),
          tags$li(tags$b("Chr"), " -- chromosome number, 1-12"),
          tags$li(tags$b("Trait"), " -- phenotype associated"),
          tags$li(tags$b("Model"), " -- GLM or MLM (or your model name)"),
          tags$li(tags$b("p_value"), " -- from TASSEL"),
          tags$li(tags$b("R2"), " -- variance explained, from TASSEL"),
          tags$li(tags$b("Mb"), " -- physical position you looked up (Gramene/BLAST)")
        ),
        uiOutput("data_status"),
        hr(),
        actionButton("load_demo", "Load Demo Data", class = "btn-warning btn-sm w-100"),
        br(), br(),
        actionButton("clear_data", "Clear table", class = "btn-outline-secondary btn-sm w-100")
      ),
      card(
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
              span("Marker table (double-click a cell to edit)"),
              uiOutput("demo_badge"))
        ),
        card_body(DTOutput("qtl_table"))
      )
    )
  ),

  nav_panel(
    "3. QTL Map",
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        downloadButton("download_plot", "Download PNG", class = "btn-primary w-100"),
        br(), br(),
        downloadButton("download_guide", "Download example data (Excel)", class = "btn-outline-secondary w-100 btn-sm"),
        hr(),
        helpText("Hover over points in the interactive view for marker details. ",
                 "Point size = -log10(p), i.e. bigger point = stronger association.")
      ),
      card(
        card_header("Interactive view"),
        card_body(plotlyOutput("plot_interactive", height = "520px"))
      ),
      card(
        card_header("Static preview (this is exactly what the PNG download produces)"),
        card_body(plotOutput("plot_static", height = "520px"))
      )
    )
  ),

  nav_panel(
    "4. Summary Table",
    card(card_body(DTOutput("summary_table")))
  ),

  nav_panel(
    "5. Methodology & Limitations",
    layout_column_wrap(
      width = 1,
      card(
        card_header("Method, step by step"),
        card_body(
          tags$ol(
            tags$li(strong("Genotyping: "), "SSR markers scored across your mapping/association panel."),
            tags$li(strong("Association testing (TASSEL): "), "GLM and/or MLM models test each ",
                    "marker against phenotype(s), producing a p-value and R\u00b2 per marker-trait pair."),
            tags$li(strong("Significance filtering (done by you, in TASSEL, before this app): "),
                    "markers surviving your chosen threshold (Bonferroni, FDR, or a suggestive ",
                    "cutoff) are called putatively associated."),
            tags$li(strong("Physical positioning: "), "each significant SSR is located on the ",
                    "reference genome -- via its cataloged Gramene position, ideally cross-checked ",
                    "or refined by BLASTing its primer sequence against the actual assembly, since ",
                    "Gramene positions can lag behind the current genome build."),
            tags$li(strong("Visualization (this app): "), "markers are placed on a chromosome ",
                    "ideogram at their Mb position, alongside the real centromere position for ",
                    "that chromosome, sized by -log10(p) and colored by trait.")
          )
        )
      ),
      card(
        card_header("Limitations you should keep in mind"),
        card_body(
          tags$ul(
            tags$li(strong("This is a putative map, not a fine-mapped one. "),
                    "SSR marker density is typically low (tens of cM apart), so the true ",
                    "causal gene could sit well away from the marker's exact position -- this ",
                    "plot shows the marker, not a confidence interval around the QTL."),
            tags$li(strong("No support interval is drawn. "),
                    "A single point per marker overstates precision; if you have LOD support ",
                    "intervals or LD-decay-based interval estimates, those are more honest to report."),
            tags$li(strong("Coordinate provenance matters. "),
                    "Positions pulled from Gramene may reference an older genome build than ",
                    "the one used here (IRGSP-1.0 / RGAP Release 7). Always sanity-check with ",
                    "the MSU/RGAP Pseudomolecule Version Converter if your markers came from ",
                    "an older source."),
            tags$li(strong("Centromere positions are Nipponbare-specific. "),
                    "They come from BAC-tiling/CentO-repeat evidence in the Nipponbare reference ",
                    "and may not exactly match centromere boundaries in the specific cultivars ",
                    "used in your own panel, due to structural variation between genotypes."),
            tags$li(strong("This app trusts your significance calls. "),
                    "It does not re-run multiple-testing correction -- if unfiltered/non-significant ",
                    "markers are pasted in, they will be plotted as if they were putative QTLs."),
            tags$li(strong("Co-localization can be coincidental, "),
                    "especially with modest population size and sparse SSR coverage -- two ",
                    "different traits mapping near the same marker is suggestive, not proof, of ",
                    "pleiotropy or linked genes."),
            tags$li(strong("Physical distance is not the same as recombination distance. "),
                    "Mb position alone doesn't tell you how much recombination separates a marker ",
                    "from a nearby causal locus -- that depends on local recombination rate, which ",
                    "varies substantially across the genome (often suppressed near centromeres).")
          ),
          p(class = "text-muted mt-2",
            em("In short: use this tool to communicate and explore where your significant markers ",
               "fall genome-wide -- not as evidence of a validated QTL on its own."))
        )
      )
    )
  )
)

# ============================== Server ======================================

server <- function(input, output, session) {

  qtl_data <- reactiveVal(DEMO_QTL)
  is_demo  <- reactiveVal(TRUE)

  observeEvent(input$qtl_file, {
    file <- input$qtl_file
    ext <- tools::file_ext(file$name)
    df <- tryCatch(
      if (ext == "csv") read.csv(file$datapath, stringsAsFactors = FALSE) else as.data.frame(read_excel(file$datapath)),
      error = function(e) NULL
    )
    if (is.null(df)) {
      showNotification("Could not read that file. Check it's a valid CSV/XLSX.", type = "error")
      return()
    }
    err <- validate_schema(df)
    if (!is.null(err)) {
      showNotification(err, type = "error", duration = 12)
      return()
    }
    df <- clean_pvalues(df)
    if (attr(df, "n_invalid_pvalues") > 0) {
      showNotification(paste0(attr(df, "n_invalid_pvalues"),
                               " row(s) had p_value <= 0 or NA and were floored so the plot doesn't break."),
                        type = "warning", duration = 10)
    }
    qtl_data(df)
    is_demo(FALSE)
    showNotification("Your data is loaded. Check the QTL Map tab.", type = "message")
  })

  observeEvent(input$load_demo, {
    qtl_data(DEMO_QTL)
    is_demo(TRUE)
    showNotification("Demo data loaded.", type = "message")
  })

  observeEvent(input$clear_data, {
    empty_df <- DEMO_QTL[0, ]
    qtl_data(empty_df)
    is_demo(FALSE)
    showNotification("Table cleared. Upload your data or click 'Load Demo Data'.", type = "message")
  })

  output$demo_badge <- renderUI({
    if (is_demo()) {
      span(class = "badge bg-warning text-dark", "DEMO DATA -- not real results")
    } else {
      span(class = "badge bg-success", "Your data")
    }
  })

  # rownames = FALSE here MUST match the rownames = FALSE used in editData()
  # below -- this was the source of the earlier "subscript out of bounds" bug.
  output$qtl_table <- renderDT({
    datatable(qtl_data(), editable = TRUE, rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE))
  }, server = FALSE)

  observeEvent(input$qtl_table_cell_edit, {
    edit_info <- input$qtl_table_cell_edit
    current <- qtl_data()

    # Defensive bounds check: ignore any edit event that doesn't map onto a
    # real row/column of the current data, instead of letting it propagate
    # into an out-of-range index later in the plotting code.
    if (is.null(edit_info) ||
        edit_info$row < 1 || edit_info$row > nrow(current) ||
        edit_info$col < 0 || edit_info$col > (ncol(current) - 1)) {
      showNotification("Edit ignored: cell reference was out of range. Try again.", type = "warning")
      return()
    }

    df <- tryCatch(
      editData(current, edit_info, "qtl_table", rownames = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) {
      showNotification("That edit could not be applied.", type = "error")
      return()
    }

    err <- validate_schema(df)
    if (!is.null(err)) {
      showNotification(err, type = "error", duration = 10)
      return()
    }
    df <- clean_pvalues(df)
    qtl_data(df)
    is_demo(FALSE)
  })

  output$data_status <- renderUI({
    err <- validate_schema(qtl_data())
    if (!is.null(err)) {
      div(style = "color:#b02a37; font-size: 0.85em;", p(err))
    } else {
      div(style = "color:#2c5f2d; font-size: 0.85em;",
          p(paste0(nrow(qtl_data()), " marker row(s) loaded and valid.")))
    }
  })

  qtl_plot_reactive <- reactive({
    req(nrow(qtl_data()) > 0, is.null(validate_schema(qtl_data())))
    create_qtl_map(clean_pvalues(qtl_data()), is_demo())
  })

  output$plot_interactive <- renderPlotly({
    req(nrow(qtl_data()) > 0, is.null(validate_schema(qtl_data())))
    create_qtl_map_plotly(clean_pvalues(qtl_data()), is_demo())
  })

  output$plot_static <- renderPlot({
    qtl_plot_reactive()
  })

  output$summary_table <- renderDT({
    req(nrow(qtl_data()) > 0, is.null(validate_schema(qtl_data())))
    df <- qtl_data() %>%
      mutate(logP = round(-log10(pmax(as.numeric(p_value), 1e-300)), 3)) %>%
      select(Marker, Chr, Trait, Model, p_value, R2, Mb, logP) %>%
      arrange(Chr, Mb)
    datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$download_plot <- downloadHandler(
    filename = function() if (is_demo()) "PutativeQTL_DEMO.png" else "PutativeQTL_Map.png",
    content = function(file) {
      req(qtl_plot_reactive())
      ggsave(file, qtl_plot_reactive(), width = 12, height = 8, dpi = 300, bg = "white")
    }
  )

  output$download_guide <- downloadHandler(
    filename = function() "PutativeQTLmapper_Example_Data.xlsx",
    content = function(file) {
      write_xlsx(list(Demo_QTL_Data = DEMO_QTL, Reference_Genome_Mb = rice_genome_ref), path = file)
    }
  )
}

shinyApp(ui, server)
