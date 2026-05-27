library(shiny)
library(bslib)
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(patchwork)
library(httr)
library(jsonlite)
library(future)
library(promises)

# ─── Async setup ──────────────────────────────────────────────────────────────
plan(multisession)

# ─── Increase upload size limit ───────────────────────────────────────────────
options(shiny.maxRequestSize = 10 * 1024^3)  # 10 GB

# ─── Load API key from config.json ────────────────────────────────────────────
if (file.exists("config.json")) {
  config  <- fromJSON("config.json")
  api_key <- config$ANTHROPIC_API_KEY
} else if (Sys.getenv("ANTHROPIC_API_KEY") != "") {
  api_key <- Sys.getenv("ANTHROPIC_API_KEY")
} else {
  warning("No API key found — set ANTHROPIC_API_KEY env var or create config.json")
  api_key <- ""
}

# ─── UI ───────────────────────────────────────────────────────────────────────
ui <- page_sidebar(
  title = "scAdvisorAI — AI-powered scRNA-seq QC",
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2C7BB6",
    base_font  = font_google("DM Sans")
  ),

  sidebar = sidebar(
    width = 300,

    fileInput("h5_files", "Upload 10X .h5 file(s)",
              multiple = TRUE,
              accept   = ".h5"),

    hr(),
    h6("QC Thresholds", style = "font-weight:bold;"),

    sliderInput("min_umi",     "Min UMI per cell",
                min = 0, max = 5000, value = 500, step = 50),
    sliderInput("min_genes",   "Min genes per cell",
                min = 0, max = 3000, value = 300, step = 50),
    sliderInput("max_mito",    "Max mitochondrial ratio",
                min = 0, max = 1,    value = 0.2, step = 0.01),
    sliderInput("min_novelty", "Min novelty score",
                min = 0, max = 1,    value = 0.8, step = 0.01),

    hr(),
    actionButton("run_qc", "Run QC",
                 class = "btn-primary w-100"),

    br(), br(),
    actionButton("run_ai", " Get AI Recommendations",
                 class = "btn-warning w-100"),

    br(), br(),
    downloadButton("download_rds", "Download Filtered Seurat Object",
                   class = "btn-success w-100")
  ),

  navset_tab(

    nav_panel(" Cell Counts",
      plotOutput("plot_cell_counts", height = "450px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("Cell counts per sample:"),
          " Total number of cells detected per sample after initial
           CellRanger filtering. Blue = kept, grey = filtered out.
           Large differences between samples may indicate technical
           variation in library preparation or sequencing depth."
        )
      )
    ),

    nav_panel(" UMI Distribution",
      plotOutput("plot_umi", height = "450px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("UMI counts per cell:"),
          " UMIs (Unique Molecular Identifiers) represent the number of
           transcripts captured per cell. Blue = cells kept, grey = cells
           filtered out by current threshold. Very low UMI counts suggest
           empty droplets or dead cells."
        )
      )
    ),

    nav_panel(" Genes per Cell",
      plotOutput("plot_genes", height = "600px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("Number of genes detected per cell:"),
          " Higher gene counts per cell generally indicate better quality
           cells. Blue = cells kept, grey = cells filtered out. Very low
           gene counts may indicate empty droplets or dead cells, while
           extremely high counts may suggest doublets."
        )
      )
    ),

    nav_panel(" UMIs vs Genes",
      plotOutput("plot_umi_genes", height = "450px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("UMIs vs genes detected:"),
          " Each point represents a cell. Blue = kept, grey = filtered out
           by current thresholds. A linear relationship is expected —
           cells deviating significantly may indicate low quality or doublets."
        )
      )
    ),

    nav_panel(" Mitochondrial Ratio",
      plotOutput("plot_mito", height = "450px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("Mitochondrial ratio:"),
          " Percentage of reads from mitochondrial genes. Blue = kept,
           grey = filtered out. High mitochondrial content (>20%) typically
           indicates damaged or dying cells. Drag the slider to see how
           the threshold affects cell retention."
        )
      )
    ),

    nav_panel(" Complexity",
      plotOutput("plot_novelty", height = "450px"),
      br(),
      div(style = "padding: 0 1rem;",
        p(style = "color:#555; font-size:13px; background:#f8f9fa;
                   border-left:3px solid #2C7BB6; padding:0.75rem;
                   border-radius:4px;",
          " ", strong("Transcriptome complexity (novelty score):"),
          " Calculated as log10(genes) / log10(UMIs). Blue = kept,
           grey = filtered out. A low novelty score indicates few unique
           genes relative to UMI count, suggesting a low-complexity or
           potentially low-quality cell. Threshold is typically 0.8."
        )
      )
    ),

    nav_panel(" Filtering Summary",
      fluidRow(
        column(6,
          h6("Before Filtering", style = "font-weight:bold; color:#2C7BB6;"),
          tableOutput("summary_before")
        ),
        column(6,
          h6("After Filtering", style = "font-weight:bold; color:#D7191C;"),
          tableOutput("summary_after")
        )
      ),
      plotOutput("plot_filter_comparison", height = "350px")
    ),

    nav_panel(" AI Advisor",
      br(),
      uiOutput("ai_status"),
      br(),
      uiOutput("ai_recommendations")
    )
  )
)

# ─── SERVER ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive values ───────────────────────────────────────────────────────
  ai_text       <- reactiveVal(NULL)
  ai_running    <- reactiveVal(FALSE)
  ai_status_msg <- reactiveVal(NULL)

  # ── Load and merge Seurat objects ─────────────────────────────────────────
  scrna_raw <- eventReactive(input$run_qc, {
    req(input$h5_files)

    withProgress(message = "Loading 10X data...", {

      scrna_list <- list()

      for (i in seq_len(nrow(input$h5_files))) {
        setProgress(i / nrow(input$h5_files),
                    detail = paste("Processing", input$h5_files$name[i]))
        sample_name <- tools::file_path_sans_ext(input$h5_files$name[i])
        counts      <- Read10X_h5(input$h5_files$datapath[i],
                                  use.names = TRUE, unique.features = TRUE)
        obj <- CreateSeuratObject(
          counts       = counts,
          min.cells    = 10,
          min.features = 100,
          project      = sample_name
        )
        obj[["Group"]]  <- sample_name
        obj <- RenameCells(obj, add.cell.id = sample_name)
        scrna_list[[i]] <- obj
      }

      setProgress(0.7, message = "Merging samples...")
      if (length(scrna_list) > 1) {
        scrna <- Reduce(merge, scrna_list)
      } else {
        scrna <- scrna_list[[1]]
      }

      setProgress(0.85, message = "Computing QC metrics...")
      scrna$log10GenesPerUMI <- log10(scrna$nFeature_RNA) /
                                log10(scrna$nCount_RNA)
      scrna$mitoPercent      <- PercentageFeatureSet(scrna, pattern = "^MT-")
      scrna$mitoRatio        <- scrna$mitoPercent / 100

      scrna@meta.data <- scrna@meta.data %>%
        dplyr::rename(nUMI  = nCount_RNA,
                      nGene = nFeature_RNA)
      scrna@meta.data$sample <- scrna@meta.data$orig.ident

      setProgress(1, message = "Done!")
      scrna
    })
  })

  # ── Reactive meta with keep flag — live with sliders ──────────────────────
  meta_flagged <- reactive({
    req(scrna_raw())
    meta       <- scrna_raw()@meta.data
    meta$kept  <- meta$nUMI             >= input$min_umi    &
                  meta$nGene            >= input$min_genes  &
                  meta$mitoRatio        <= input$max_mito   &
                  meta$log10GenesPerUMI >= input$min_novelty
    meta
  })

  # ── Filtered Seurat object ────────────────────────────────────────────────
  scrna_filtered <- reactive({
    req(scrna_raw())
    subset(scrna_raw(),
           subset = nUMI             >= input$min_umi    &
                    nGene            >= input$min_genes  &
                    mitoRatio        <= input$max_mito   &
                    log10GenesPerUMI >= input$min_novelty)
  })

  # ── Shared theme ──────────────────────────────────────────────────────────
  qc_theme <- function() {
    theme_bw() +
    theme(
      axis.text        = element_text(size = 11, face = "bold"),
      axis.title       = element_text(size = 13, face = "bold"),
      plot.title       = element_text(size = 14, face = "bold",
                                      color = "darkblue", hjust = 0.5),
      strip.text       = element_text(size = 11, face = "bold"),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black")
    )
  }

  kept_colors <- c("TRUE" = "#2C7BB6", "FALSE" = "grey80")
  kept_labels <- c("TRUE" = "Kept",    "FALSE" = "Filtered")

  # ── Tab 1: Cell counts ────────────────────────────────────────────────────
  output$plot_cell_counts <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()
    ggplot(meta, aes(x = sample, fill = kept)) +
      geom_bar() +
      scale_fill_manual(values = kept_colors, labels = kept_labels) +
      qc_theme() +
      theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "right") +
      labs(title = "Cell Counts per Sample (kept vs filtered)",
           x = "Sample", y = "Number of Cells", fill = "")
  })

  # ── Tab 2: UMI distribution ───────────────────────────────────────────────
  output$plot_umi <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()
    ggplot(meta, aes(x = nUMI, fill = kept, color = kept)) +
      geom_density(alpha = 0.35) +
      scale_x_log10() +
      scale_fill_manual(values  = kept_colors, labels = kept_labels) +
      scale_color_manual(values = kept_colors, labels = kept_labels) +
      geom_vline(xintercept = input$min_umi,
                 color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~sample) +
      qc_theme() +
      theme(legend.position = "right") +
      labs(title = "UMI Counts per Cell",
           x = "UMI Count (log10)", y = "Cell Density",
           fill = "", color = "")
  })

  # ── Tab 3: Genes per cell ─────────────────────────────────────────────────
  output$plot_genes <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()

    p1 <- ggplot(meta, aes(x = nGene, fill = kept, color = kept)) +
      geom_density(alpha = 0.35) +
      scale_x_log10() +
      scale_fill_manual(values  = kept_colors, labels = kept_labels) +
      scale_color_manual(values = kept_colors, labels = kept_labels) +
      geom_vline(xintercept = input$min_genes,
                 color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~sample) +
      qc_theme() +
      theme(legend.position = "right") +
      labs(title = "Genes Detected per Cell",
           x = "Gene Count (log10)", y = "Cell Density",
           fill = "", color = "")

    p2 <- ggplot(meta, aes(x = sample, y = log10(nGene), fill = kept)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
      scale_fill_manual(values = kept_colors, labels = kept_labels) +
      qc_theme() +
      theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "right",
            plot.margin     = margin(t = 20, r = 10, b = 10, l = 10)) +
      labs(title = "Genes per Cell by Sample",
           x = "Sample", y = "log10(nGene)", fill = "")

    p1 / plot_spacer() / p2 + plot_layout(heights = c(5, 0.3, 5))
  })

  # ── Tab 4: UMIs vs genes ──────────────────────────────────────────────────
  output$plot_umi_genes <- renderPlot({
    req(scrna_raw())
    meta         <- meta_flagged()
    meta_ordered <- meta[order(meta$kept), ]
    ggplot(meta_ordered, aes(x = nUMI, y = nGene, color = kept)) +
      geom_point(alpha = 0.4, size = 0.6) +
      scale_color_manual(values = kept_colors, labels = kept_labels) +
      scale_x_log10() +
      scale_y_log10() +
      geom_vline(xintercept = input$min_umi,
                 color = "red", linetype = "dashed", linewidth = 0.8) +
      geom_hline(yintercept = input$min_genes,
                 color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_wrap(~sample) +
      qc_theme() +
      theme(legend.position = "right") +
      labs(title = "UMIs vs Genes Detected",
           x = "UMI Count (log10)", y = "Gene Count (log10)",
           color = "")
  })

  # ── Tab 5: Mitochondrial ratio ────────────────────────────────────────────
  output$plot_mito <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()
    ggplot(meta, aes(x = mitoRatio, fill = kept, color = kept)) +
      geom_density(alpha = 0.35) +
      scale_x_log10() +
      scale_fill_manual(values  = kept_colors, labels = kept_labels) +
      scale_color_manual(values = kept_colors, labels = kept_labels) +
      geom_vline(xintercept = input$max_mito,
                 color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~sample) +
      qc_theme() +
      theme(legend.position = "right") +
      labs(title = "Mitochondrial Read Ratio per Cell",
           x = "Mitochondrial Ratio (log10)", y = "Cell Density",
           fill = "", color = "")
  })

  # ── Tab 6: Complexity ─────────────────────────────────────────────────────
  output$plot_novelty <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()
    ggplot(meta, aes(x = log10GenesPerUMI, fill = kept, color = kept)) +
      geom_density(alpha = 0.35) +
      scale_fill_manual(values  = kept_colors, labels = kept_labels) +
      scale_color_manual(values = kept_colors, labels = kept_labels) +
      geom_vline(xintercept = input$min_novelty,
                 color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~sample) +
      qc_theme() +
      theme(legend.position = "right") +
      labs(title = "Transcriptome Complexity (Novelty Score)",
           x = "log10(Genes per UMI)", y = "Cell Density",
           fill = "", color = "")
  })

  # ── Tab 7: Filtering summary ──────────────────────────────────────────────
  output$summary_before <- renderTable({
    req(scrna_raw())
    scrna_raw()@meta.data %>%
      group_by(sample) %>%
      summarise(
        Cells        = n(),
        Median_UMI   = round(median(nUMI)),
        Median_Genes = round(median(nGene)),
        Median_Mito  = round(median(mitoRatio), 3),
        .groups = "drop"
      ) %>%
      rename(Sample = sample) %>%
      mutate(Stage = "Before Filtering") %>%
      select(Stage, everything())
  })

  output$summary_after <- renderTable({
    req(scrna_raw())
    meta_flagged() %>%
      filter(kept) %>%
      group_by(sample) %>%
      summarise(
        Cells        = n(),
        Median_UMI   = round(median(nUMI)),
        Median_Genes = round(median(nGene)),
        Median_Mito  = round(median(mitoRatio), 3),
        .groups = "drop"
      ) %>%
      rename(Sample = sample) %>%
      mutate(Stage = "After Filtering") %>%
      select(Stage, everything())
  })

  output$plot_filter_comparison <- renderPlot({
    req(scrna_raw())
    meta <- meta_flagged()

    before <- meta %>%
      group_by(sample) %>%
      summarise(Cells = n(), .groups = "drop") %>%
      mutate(Stage = "Before")

    after <- meta %>%
      filter(kept) %>%
      group_by(sample) %>%
      summarise(Cells = n(), .groups = "drop") %>%
      mutate(Stage = "After")

    combined       <- bind_rows(before, after)
    combined$Stage <- factor(combined$Stage, levels = c("Before", "After"))

    ggplot(combined, aes(x = sample, y = Cells, fill = Stage)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Before" = "#2C7BB6",
                                   "After"  = "#D7191C")) +
      qc_theme() +
      theme(legend.position = "right",
            axis.text.x     = element_text(angle = 45,
                                           vjust = 1, hjust = 1)) +
      labs(title = "Cells Before vs After Filtering",
           x = "Sample", y = "Number of Cells", fill = "Stage")
  })

  # ── Tab 8: AI Advisor status ──────────────────────────────────────────────
  output$ai_status <- renderUI({
    if (isTRUE(ai_running())) {
      div(
        style = paste(
          "background:#fff3cd; border-left:4px solid #ffc107;",
          "padding:1rem; border-radius:4px; margin-bottom:1rem;"
        ),
        tags$span(
          style = "color:#856404; font-weight:500;",
          tags$span(class = "spinner-border spinner-border-sm",
                    role  = "status",
                    style = "margin-right:8px;"),
          "Consulting AI advisor — analyzing your QC metrics...
           This may take 10–20 seconds."
        )
      )
    } else if (!is.null(ai_status_msg())) {
      div(
        style = paste(
          "background:#d1e7dd; border-left:4px solid #198754;",
          "padding:1rem; border-radius:4px; margin-bottom:1rem;"
        ),
        tags$span(style = "color:#0f5132; font-weight:500;",
                  "✅ ", ai_status_msg())
      )
    }
  })

  # ── Tab 8: AI Advisor — async API call ───────────────────────────────────
  observeEvent(input$run_ai, {
    req(scrna_raw())

    ai_running(TRUE)
    ai_status_msg(NULL)
    ai_text(NULL)

    meta <- scrna_raw()@meta.data

    qc_stats <- meta %>%
      group_by(sample) %>%
      summarise(
        n_cells        = n(),
        median_umi     = round(median(nUMI)),
        q1_umi         = round(quantile(nUMI, 0.25)),
        q3_umi         = round(quantile(nUMI, 0.75)),
        median_genes   = round(median(nGene)),
        pct_high_mito  = round(mean(mitoRatio > 0.2) * 100, 1),
        median_mito    = round(median(mitoRatio), 3),
        median_novelty = round(median(log10GenesPerUMI), 3),
        .groups = "drop"
      )

    stats_text <- apply(qc_stats, 1, function(row) {
      paste0(
        "Sample: ", row["sample"], "\n",
        "  Total cells: ", row["n_cells"], "\n",
        "  Median UMI: ", row["median_umi"],
        " (Q1=", row["q1_umi"], ", Q3=", row["q3_umi"], ")\n",
        "  Median genes per cell: ", row["median_genes"], "\n",
        "  % cells with >20% mito reads: ", row["pct_high_mito"], "%\n",
        "  Median mitochondrial ratio: ", row["median_mito"], "\n",
        "  Median novelty score: ", row["median_novelty"], "\n"
      )
    })

    prompt <- paste0(
      "You are an expert single-cell RNA-seq bioinformatician. ",
      "Analyze the following QC metrics from a 10X Genomics scRNA-seq ",
      "experiment and provide filtering recommendations.\n\n",
      paste(stats_text, collapse = "\n"),
      "\nBased on these metrics, please provide:\n",
      "1. Recommended thresholds: min_nUMI, min_nGene, ",
      "max_mitoRatio, min_novelty_score\n",
      "2. Reasoning for each threshold based on the data\n",
      "3. Any quality concerns or flags about the data\n",
      "4. A 2-3 sentence Methods section describing the QC ",
      "filtering performed\n\n",
      "Format your response with clear headers for each section."
    )

    key <- api_key

    future_promise({
      httr::POST(
        url  = "https://api.anthropic.com/v1/messages",
        httr::add_headers(
          "Content-Type"      = "application/json",
          "x-api-key"         = key,
          "anthropic-version" = "2023-06-01"
        ),
        body = jsonlite::toJSON(list(
          model      = "claude-sonnet-4-20250514",
          max_tokens = 1000L,
          messages   = list(
            list(role = "user", content = prompt)
          )
        ), auto_unbox = TRUE),
        encode = "json",
        httr::timeout(60)
      )
    }) %...>% (function(response) {
      result <- httr::content(response, as = "parsed")
      text   <- result$content[[1]]$text
      ai_text(text)
      ai_running(FALSE)
      ai_status_msg("AI recommendations ready — sliders updated automatically.")
    }) %...!% (function(e) {
      ai_text(paste0("Error contacting AI advisor: ", e$message,
                     "\n\nCheck your API key in config.json and try again."))
      ai_running(FALSE)
      ai_status_msg(NULL)
    })
  })

  # ── Update sliders when AI result arrives ─────────────────────────────────
  observeEvent(ai_text(), {
    req(ai_text())
    text <- ai_text()

    extract_num <- function(pattern, text, default) {
      m <- regmatches(text, regexpr(pattern, text, perl = TRUE))
      if (length(m) > 0) as.numeric(gsub("[^0-9.]", "", m)) else default
    }

    updateSliderInput(session, "min_umi",
      value = min(extract_num("min.nUMI[^0-9]*([0-9]+)",    text, 500), 5000))
    updateSliderInput(session, "min_genes",
      value = min(extract_num("min.nGene[^0-9]*([0-9]+)",   text, 300), 3000))
    updateSliderInput(session, "max_mito",
      value = min(extract_num("max.mito[^0-9]*([0-9.]+)",   text, 0.2), 1))
    updateSliderInput(session, "min_novelty",
      value = min(extract_num("min.novelty[^0-9]*([0-9.]+)",text, 0.8), 1))
  })

  # ── Render AI response ────────────────────────────────────────────────────
  output$ai_recommendations <- renderUI({
    req(ai_text())
    tagList(
      div(
        style = paste(
          "background:#f8f9fa; border-left:4px solid #2C7BB6;",
          "padding:1.5rem; border-radius:4px;",
          "font-family:'DM Sans', sans-serif; line-height:1.7;"
        ),
        h4(" AI QC Recommendations",
           style = "color:#2C7BB6; margin-bottom:1rem;"),
        p(style = "color:#555; font-size:13px; margin-bottom:1rem;",
          "Thresholds have been applied to the sliders automatically. ",
          "Review and adjust as needed before downloading."),
        HTML(gsub("\n", "<br>", ai_text()))
      )
    )
  })

  # ── Download filtered Seurat object ───────────────────────────────────────
  output$download_rds <- downloadHandler(
    filename = function() {
      paste0("filtered_seurat_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(scrna_filtered())
      withProgress(message = "Saving filtered Seurat object...", {
        saveRDS(scrna_filtered(), file)
      })
    }
  )
}

# ─── RUN ──────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
