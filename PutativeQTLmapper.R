# QTL Mapper Shiny App
# Author: Zafir Ahmad Naik (PhD Genetics and Plant Breeding)
# Purpose: Automated integrated GWAS-QTL visualization and demonstration with dummy data

library(shiny)
library(DT)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(gridExtra)
library(readxl)
library(writexl)

#------------------ Example (Dummy) Data ------------------
# This dummy dataset mimics the format of TASSEL GWAS output and chromosome map files.

dummy_chr <- data.frame(
  Chr = 1:12,
  Length_cM = c(158.7,132.5,152.8,123.4,114.5,128.6,108.4,111.2,99.8,97.6,107.2,88.3),
  Length_Mb = c(43.3,35.9,36.4,35.5,30.0,31.2,29.7,28.4,23.0,23.2,29.0,27.5)
)

dummy_qtl <- data.frame(
  Marker = c("RM285", "RM27462", "RM240", "RM27462", "RM101"),
  Chr = c(9, 12, 2, 12, 12),
  Trait = c("1000SW", "50F", "NT", "50F", "PH"),
  Model = c("GLM", "GLM", "GLM", "MLM", "MLM"),
  p_value = c(0.0035, 0.005, 0.0083, 0.0089, 0.0087),
  R2 = c(0.056, 0.0518, 0.0433, 0.045, 0.0441),
  cM = c(1.0, 45.0, 85.0, 45.0, 48.2),
  Mb = c(6.8, 15.0, 29.6, 15.0, 9.0)
)

#------------------ QTL Plot Function ------------------
create_qtl_map <- function(qtl_data, chr_data, position_type = "cM", title_suffix = "") {
  if (position_type == "cM") {
    qtl_data$position <- qtl_data$cM
    chr_lengths <- chr_data$Length_cM
    y_label <- "Genetic Position (cM)"
    max_pos <- max(chr_lengths) + 10
  } else {
    qtl_data$position <- qtl_data$Mb
    chr_lengths <- chr_data$Length_Mb
    y_label <- "Physical Position (Mb)"
    max_pos <- max(chr_lengths) + 5
  }
  
  qtl_data <- qtl_data %>%
    mutate(logP = -log10(p_value),
           Chr = factor(Chr, levels = unique(chr_data$Chr)),
           Label = ifelse(duplicated(paste(Marker, Chr)) |
                            duplicated(paste(Marker, Chr), fromLast = TRUE),
                          paste0(Marker, "\n", Trait, " (", Model, ")"),
                          paste0(Marker, "\n", Trait)))
  
  chr_backbone <- data.frame(
    Chr = factor(chr_data$Chr, levels = chr_data$Chr),
    Length = chr_lengths,
    Start = 0
  )
  
  ggplot() +
    geom_segment(data = chr_backbone, aes(x = Chr, xend = Chr, y = Start, yend = Length),
                 size = 8, color = "lightgray", alpha = 0.7) +
    geom_point(data = chr_backbone, aes(x = Chr, y = Length * 0.4),
               size = 4, shape = 15, color = "darkgray") +
    geom_point(data = qtl_data, aes(x = Chr, y = position, size = -log10(p_value), color = Trait),
               alpha = 0.8, stroke = 1.2) +
    geom_text_repel(data = qtl_data,
                    aes(x = Chr, y = position, label = Label, color = Trait),
                    box.padding = 0.8, point.padding = 0.3,
                    segment.color = "gray50", size = 3.5, fontface = "bold", max.overlaps = 20) +
    scale_x_discrete(name = "Chromosomes", labels = paste0("Chr", chr_data$Chr)) +
    scale_y_continuous(name = y_label, limits = c(-max_pos * 0.05, max_pos), expand = c(0.02, 0)) +
    scale_color_brewer(name = "Trait", palette = "Set2") +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("QTL Map ", title_suffix),
      subtitle = paste0("Chromosomal distribution of associations (", position_type, " positions)"),
      caption = "Marker size = -log10(p-value); Centromeres shown as gray squares"
    )
}

#------------------ UI ------------------
ui <- fluidPage(
  titlePanel("📊 Integrated GWAS–QTL Mapper"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Data Input Options"),
      tabsetPanel(
        tabPanel("Upload Excel/CSV",
                 fileInput("chr_file", "Upload Chromosome Data", accept = c(".csv", ".xlsx")),
                 fileInput("qtl_file", "Upload QTL/GWAS Data", accept = c(".csv", ".xlsx")),
                 helpText("Expected format: 
                         • Chromosome file: Chr, Length_cM, Length_Mb
                         • QTL file: Marker, Chr, Trait, Model, p_value, R2, cM, Mb")
        ),
        tabPanel("Manual Entry / Dummy Data",
                 p("This example demonstrates the structure and visual output expected 
                   when using TASSEL or similar GWAS outputs."),
                 DTOutput("chr_manual"),
                 br(),
                 DTOutput("qtl_manual"),
                 br(),
                 h5("📘 Dummy Data Description:"),
                 helpText(HTML("
                 <ul>
                   <li><b>Marker:</b> SSR or SNP identifier, as listed in TASSEL results.</li>
                   <li><b>Chr:</b> Chromosome number where the marker is located.</li>
                   <li><b>Trait:</b> Phenotypic trait under association study (e.g., 1000SW = 1000-seed weight).</li>
                   <li><b>Model:</b> Statistical model used in GWAS (e.g., GLM, MLM).</li>
                   <li><b>p_value:</b> Significance level of marker-trait association.</li>
                   <li><b>R2:</b> Percentage of phenotypic variance explained by the marker (R² value).</li>
                   <li><b>cM:</b> Genetic position of marker on the chromosome (centiMorgans).</li>
                   <li><b>Mb:</b> Physical position of marker in megabases (from reference genome data).</li>
                 </ul>
                 <p>This format replicates the structure of TASSEL output tables, where
                 p-values and R² values are directly obtained from GWAS results, while
                 chromosomal locations (cM and Mb) are aligned using published marker
                 position databases.</p>
                 "))
        )
      ),
      hr(),
      actionButton("generate", "Generate QTL Maps", class = "btn-primary"),
      hr(),
      downloadButton("download_plots", "📥 Download Combined Plot (PNG)"),
      br(), br(),
      downloadButton("download_guide", "📄 Download Dummy Example Data (Excel)")
    ),
    
    mainPanel(
      h3("QTL Map Outputs"),
      plotOutput("plot_cM", height = "500px"),
      plotOutput("plot_Mb", height = "500px"),
      hr(),
      h4("Summary of Associations"),
      DTOutput("summary_table")
    )
  )
)

#------------------ Server ------------------
server <- function(input, output, session) {
  
  # Render editable dummy tables
  output$chr_manual <- renderDT({
    datatable(dummy_chr, editable = TRUE)
  })
  output$qtl_manual <- renderDT({
    datatable(dummy_qtl, editable = TRUE)
  })
  
  # Reactive values
  chr_data <- reactiveVal(dummy_chr)
  qtl_data <- reactiveVal(dummy_qtl)
  
  # Read uploaded files
  observeEvent(input$chr_file, {
    file <- input$chr_file
    if (is.null(file)) return()
    ext <- tools::file_ext(file$name)
    df <- if (ext == "csv") read.csv(file$datapath) else read_excel(file$datapath)
    chr_data(df)
  })
  
  observeEvent(input$qtl_file, {
    file <- input$qtl_file
    if (is.null(file)) return()
    ext <- tools::file_ext(file$name)
    df <- if (ext == "csv") read.csv(file$datapath) else read_excel(file$datapath)
    qtl_data(df)
  })
  
  # Generate plots
  observeEvent(input$generate, {
    req(chr_data(), qtl_data())
    
    plot_cM <- create_qtl_map(qtl_data(), chr_data(), "cM", " - Genetic Map")
    plot_Mb <- create_qtl_map(qtl_data(), chr_data(), "Mb", " - Physical Map")
    
    output$plot_cM <- renderPlot({ plot_cM })
    output$plot_Mb <- renderPlot({ plot_Mb })
    
    output$summary_table <- renderDT({
      qtl_data() %>%
        select(Marker, Chr, Trait, Model, p_value, R2, cM, Mb) %>%
        arrange(Chr, cM)
    })
    
    output$download_plots <- downloadHandler(
      filename = function() "Combined_QTL_Map.png",
      content = function(file) {
        ggsave(file, grid.arrange(plot_cM, plot_Mb, ncol = 2),
               width = 18, height = 8, dpi = 300, bg = "white")
      }
    )
  })
  
  # Download dummy example files
  output$download_guide <- downloadHandler(
    filename = function() "Dummy_QTL_Guide.xlsx",
    content = function(file) {
      write_xlsx(list(Chromosome_Data = dummy_chr, QTL_Data = dummy_qtl), path = file)
    }
  )
}

shinyApp(ui, server)
