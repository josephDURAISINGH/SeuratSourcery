options(shiny.maxRequestSize = 4096 * 1024^2)  # 4 GB

# Make sure these come from your package
# source("R/summonData.R")
# source("R/getSourceryReport.R")

# helper to decompress the commonly found zip files
gunzip_file <- function(src, dst) {
  con_in  <- gzfile(src, "rb")
  con_out <- file(dst, "wb")
  writeBin(readBin(con_in, what="raw", n = 1e8), con_out)
  close(con_in); close(con_out)
}


# Helper for building the local directory
prepare_upload_dir <- function(uploaded) {

  temp_dir <- tempfile("sourcery_upload_")
  dir.create(temp_dir, showWarnings = FALSE)

  # Loop through uploaded files
  for (i in seq_len(nrow(uploaded))) {

    src <- uploaded$datapath[i]
    orig_name <- uploaded$name[i]
    dst <- file.path(temp_dir, orig_name)

    # Detect folders manually (Shiny does not natively handle folders well)
    if (dir.exists(src)) {
      # Copy full directory (macOS only supports some)
      file.copy(src, dst, recursive = TRUE)
      next
    }

    # ----------- Case 1: GZ FILES -----------
    if (grepl("\\.gz$", orig_name, ignore.case = TRUE)) {
      message("Decompressing: ", orig_name)

      unzipped_name <- sub("\\.gz$", "", orig_name)
      dst_unzipped <- file.path(temp_dir, unzipped_name)

      gunzip_file(src, dst_unzipped)
      next
    }

    # ----------- Case 2: NORMAL FILES -----------
    file.copy(src, dst)
  }

  # ----------- Auto-group 10X bundles into folders -----------
  # Detect presence of matrix.mtx, barcodes, features
  files <- list.files(temp_dir, full.names = TRUE)

  tenx_patterns <- c("matrix\\.mtx$", "barcodes", "features|genes")

  tenx_files <- files[
    sapply(files, function(f)
      any(grepl(tenx_patterns, basename(f), ignore.case = TRUE))
    )
  ]

  if (length(tenx_files) > 0) {
    tenx_dir <- file.path(temp_dir, "10x_bundle")
    dir.create(tenx_dir, showWarnings = FALSE)

    for (f in tenx_files)
      file.rename(f, file.path(tenx_dir, basename(f)))
  }

  return(temp_dir)
}



# ============================================================
# UI
# ============================================================

ui <- shiny::fluidPage(

  titlePanel("SeuratSourcery — Dataset Inspection Dashboard"),

  sidebarLayout(
    sidebarPanel(
      h3("1. Load Datasets"),
      fileInput("files", "Upload files or folders",
                multiple = TRUE,
                accept = c(
                  ".rds", ".rds.gz",
                  ".h5ad", ".h5ad.gz",
                  ".csv", ".csv.gz",
                  ".tsv", ".tsv.gz",
                  ".mtx", ".mtx.gz",
                  ".h5"
                )),

      actionButton("load_data", "Summon Data", class = "btn-primary"),

      hr(),

      h3("2. Activate Runes (Cleaning)"),

      checkboxInput("fix_names", "Fix gene names (uppercase, strip versions)", TRUE),
      checkboxInput("drop_zero", "Drop zero-expression genes", TRUE),
      checkboxInput("drop_dups", "Drop duplicated genes", TRUE),
      checkboxInput("normalize", "Normalize (log + HVGs)", FALSE),

      actionButton("activate", "Activate Rune", class = "btn-warning"),

      hr(),

      h3("3. Generate Report"),

      checkboxGroupInput("plots_to_show", "Select Plots:",
                         choices = c(
                           "Gene Count" = "genes",
                           "Cell Count" = "cells",
                           "Nonzero Features" = "features",
                           "Mean UMI per Cell" = "umi",
                           "Mito %" = "mito",
                           "Conserved Genes" = "conserved"
                         ),
                         selected = c("genes")),


      actionButton("run_report", "Update Sourcery Report", class = "btn-success")

    ),

    mainPanel(
      h2("Sourcery Report"),
      uiOutput("plot_area"),
      hr(),
      h3("Dataset Summary"),
      tableOutput("summary_table")
    )
  )
)

server <- function(input, output, session) {

  # STORAGE --------------------------------------------------------

  raw_datasets <- reactiveVal(NULL)
  activated_datasets <- reactiveVal(NULL)
  report_cache <- reactiveVal(NULL)


  # 1. SUMMON DATA --------------------------------------------------

  observeEvent(input$load_data, {
    req(input$files)

    temp_dir <- prepare_upload_dir(input$files)

    ds <- tryCatch(
      summonData(temp_dir),
      error = function(e) {
        showNotification(
          paste("Error loading datasets:", e$message),
          type = "error"
        )
        return(NULL)
      }
    )

    raw_datasets(ds)
    activated_datasets(NULL)
    report_cache(NULL)

    showNotification("Datasets loaded successfully!", type = "message")
  })


  # 2. ACTIVATE RUNE ------------------------------------------------

  observeEvent(input$activate, {
    req(raw_datasets())

    cleaned <- activateRune(
      raw_datasets(),
      normalize = input$normalize,
      fix_gene_names = input$fix_names,
      drop_zero_genes = input$drop_zero,
      drop_duplicates = input$drop_dups
    )

    activated_datasets(cleaned)
    report_cache(NULL)

    showNotification("Runes activated!", type = "warning")
  })


  # 3. GENERATE REPORT ---------------------------------------------

  observeEvent(input$run_report, {
    req(activated_datasets())

    report <- runeInspection(activated_datasets())
    report_cache(report)

    showNotification("Report generated!", type = "message")
  })


  # OUTPUTS --------------------------------------------------------

  output$plot_area <- renderUI({
    req(report_cache())

    plots <- getSourceryReport(report_cache(), show_plots = FALSE)$plots
    selected <- input$plots_to_show

    plot_output_list <- lapply(selected, function(pname) {
      plotname <- paste0("plot_", pname)
      output[[plotname]] <- renderPlot({ plots[[pname]] })
      plotOutput(plotname, height = "300px")
    })

    do.call(tagList, plot_output_list)
  })

  output$summary_table <- renderTable({
    req(report_cache())
    report_cache()$summary
  })

}

shiny::shinyApp(ui, server)
