

# Install packages (if not already present) and Load libraries ---------------------------------------------------------------
list.of.packages <- c("shiny","DT","readr","shinyjs","shinycssloaders")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

options(shiny.maxRequestSize = 100 * 1024^2)

source("functions.R")


# Define UI ---------
ui <- fluidPage(
  
  useShinyjs(),  # Enable shinyjs in UI
  titlePanel("GSEA and Volcano Plot Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      #fileInput("file", "Upload File", accept = c(".tsv", ".txt")),
      #actionButton("use_example", "Use Example Dataset"),  # New button
      fluidRow(
        column(8, fileInput("file", "Upload File", accept = c(".tsv", ".txt"))),
        column(4, tags$span("...or use an example dataset"), 
               actionButton("use_example", "Example Dataset"))),
      selectInput("id_col", "Select ID Column", choices = NULL),
      selectInput("keytype", "Select Keytype", choices = c("Ensembl", "Uniprot", "Entrez", "Refseq", "Symbol")),
      selectInput("sort_col", "Select Sorting Column", choices = NULL),
      actionButton("run_gsea", "Run GSEA"),
      br(), br(),  # Adds space
      selectInput("logFC_col", "Select LogFC Column", choices = NULL),
      selectInput("pval_col", "Select P-Value Column", choices = NULL),   
      #selectInput("stat_col", "Select Stat Column", choices = NULL),
      actionButton("run_volcano", "Run Volcano Plot"),
      uiOutput("volcano_inputs"),
      hr(),  
      downloadButton("download_data", "Download Dataset"),
      downloadButton("download_gsea", "Download GSEA Results")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview",
                 DTOutput("table", width = "80%", height = "300px")
        ),
        tabPanel("GSEA Results", 
                 plotOutput("gsea_plot", height = "500px") %>% withSpinner(color="#0dc5c1"), 
                 DTOutput("gsea_table", width = "80%", height = "300px") 
        ),
        tabPanel("Volcano Plot", 
                 plotOutput("volcano_plot", height = "500px")
        )
      )
    )),
    
    tags$head(
      tags$style(HTML("
    table.dataTable.compact tbody td {
      padding: 2px 6px !important;
      font-size: 12px !important
    }
  "))
 )
)


# Helper function to find the first matching column -------
select_best_match <- function(preferred, col_names, fallback = NULL) {
  match <- preferred[preferred %in% col_names]
  if (length(match) > 0) {
    return(match[1])  # Use the first match
  } else if (!is.null(fallback) && fallback %in% col_names) {
    return(fallback)  # Use fallback if available
  } else {
    return(col_names[1])  # Default to the first column
  }
}

# Server ---------
server <- function(input, output, session) {
  example_url <- "https://raw.githubusercontent.com/opalasca/inflammatome_package_sandbox/main/data/test_datasets/02_GSE138614_MS_CTL.tsv"
  # Reactive value to track whether example data is being used
  example_mode <- reactiveVal(FALSE)
  # Reactive block to read the dataset
  dataset <- reactive({
    if (example_mode()) {
      # Load example dataset from GitHub
      showNotification("Loading example dataset...", type = "message", duration = 3)
      data <- read.delim(example_url, header = TRUE, sep = "\t")
    } else {
    req(input$file)  
    data <- read.delim(input$file$datapath, header = TRUE, sep = "\t")
    data$ID <- data[[1]]  # Set first column as ID
    col_names <- colnames(data)
    # Automatically select relevant columns
    default_sort_col <- select_best_match(c("stat", "t"), col_names, fallback = col_names[2])
    default_logFC_col <- select_best_match(c("log2FoldChange", "logFC"), col_names)
    default_pval_col <- select_best_match(c("pvalue", "P.Value", "pval"), col_names)
    # Update selectInputs
    updateSelectInput(session, "id_col", choices = col_names, selected = col_names[1])
    updateSelectInput(session, "sort_col", choices = col_names, selected = default_sort_col)
    updateSelectInput(session, "logFC_col", choices = col_names, selected = default_logFC_col)
    updateSelectInput(session, "pval_col", choices = col_names, selected = default_pval_col)
    return(data)
  }})
  
  observeEvent(input$use_example, {
    example_mode(TRUE)  # Switch to example mode
    data <- dataset()   # Load dataset
    
    # Ensure dataset is valid
    req(data)
    
    # Update column choices
    col_names <- colnames(data)
    
    # Automatically select relevant columns for the example dataset
    default_sort_col <- select_best_match(c("stat", "t"), col_names, fallback = col_names[2])
    default_logFC_col <- select_best_match(c("log2FoldChange", "logFC"), col_names)
    default_pval_col <- select_best_match(c("pvalue", "P.Value", "pval"), col_names)
    
    # Update UI with detected columns
    updateSelectInput(session, "id_col", choices = col_names, selected = col_names[1])
    updateSelectInput(session, "sort_col", choices = col_names, selected = default_sort_col)
    updateSelectInput(session, "logFC_col", choices = col_names, selected = default_logFC_col)
    updateSelectInput(session, "pval_col", choices = col_names, selected = default_pval_col)
    
    showNotification("Example dataset loaded!", type = "message", duration = 3)
  })
  
  # When a file is uploaded, switch back to user mode
  observeEvent(input$file, {
    example_mode(FALSE)  # Switch to user mode
    dataset()  # Reload dataset
  })
  
  # Processed Data
  processed_data <- reactive({
    req(dataset(), input$id_col, input$keytype)  # Ensure required inputs are available
    
    data <- dataset()
    id_col <- input$id_col
    keytype <- tolower(input$keytype)
    
    # Validate ID column
    if (!(id_col %in% colnames(data))) {
      showNotification("⚠️ Selected ID column does not exist!", type = "error", duration = 5)
      return(NULL)
    }
    
    # Try processing the data, handle keytype errors gracefully
    result <- tryCatch({
      processed <- process_input_data(data, id = id_col, keytype = keytype)
      
      if (is.null(processed)) {
        showNotification("⚠️ Identifiers could not be mapped. Please check your keytype selection.", 
                         type = "warning", duration = 5)
        return(NULL)  # Allow re-selection instead of stopping
      }
      
      return(processed)
    }, error = function(e) {
      showNotification("⚠️ Invalid keytype selected. Please try again.", 
                       type = "error", duration = 5)
      return(NULL)  # Keep app running instead of stopping
    })
    
    return(result)
  })
  
  
  # GSEA Analysis
  gsea_results <- eventReactive(input$run_gsea, {
    req(processed_data(), input$sort_col)
    
    showNotification("Running GSEA analysis...", type = "message", duration = 5)
    
    data <- processed_data()
    id_col <- isolate(input$id_col)  # Ensure the input is only read when the button is clicked
    sort_col <- isolate(input$sort_col)
    keytype <- isolate(input$keytype)
    
    
    if (is.null(data)) {
      showNotification("⚠️ Processed data is NULL. Check keytype selection.", type = "error", duration = 5)
      return(NULL)
    }

    gene_sets <- get_gene_sets(keytype)
    
    result <- gsea_analysis(data, gene_sets, sorting_value_col_name = sort_col, name = "data")
    result <- tryCatch({
      gsea_analysis(data, gene_sets, sorting_value_col_name = sort_col, name = "data")
      }, error = function(e) {
      showNotification("⚠️ GSEA failed: Check keytype or sorting column.", type = "error", duration = 5)
      return(NULL)
    })
    
    showNotification("GSEA analysis complete!", type = "message", duration = 3)
    
    return(result)
  })
  
  # Show a notification and disable the button when running GSEA
  observeEvent(input$run_gsea, {
    shinyjs::disable("run_gsea")  # Disable button
    showNotification("Running GSEA analysis...", type = "message", duration = 5)
    
    # Trigger the GSEA analysis by calling gsea_results
    gsea_results()  
    
    # Show completion message and re-enable button after GSEA finishes
    showNotification("GSEA analysis complete!", type = "message", duration = 3)
    shinyjs::enable("run_gsea")  # Re-enable button after completion
  })
  
  # Enable the runGSEA button when keytype or other relevant inputs change
  observeEvent(input$keytype, {
    shinyjs::enable("run_gsea")
    data <- processed_data()  # Re-trigger data processing
  })
  
  observeEvent(input$sort_col, {
    shinyjs::enable("run_gsea")
  })
  
  observeEvent(input$id_col, {
    shinyjs::enable("run_gsea")
  })
  
  # Volcano Plot
  volcano_results <- eventReactive(input$run_volcano, {
    req(processed_data(), input$logFC_col, input$pval_col)
    data <- processed_data()
    plot_volcano(data, keytype=input$keytype, logFC_col_name = input$logFC_col, pval_col_name = input$pval_col)
  })
  
  # Show a notification and disable the button when running volcano
  observeEvent(input$run_volcano, {
    shinyjs::disable("run_volcano")  # Disable button
    showNotification("Generating volcano plots...", type = "message", duration = 5)
    
    # Trigger the Volcano plot generation by calling gsea_results
    #volcano_results()  
    
    # Show completion message and re-enable button after GSEA finishes
    showNotification("Volcano plots complete!", type = "message", duration = 3)
    shinyjs::enable("run_volcano")  # Re-enable button after completion
  })
  
  
  # Render GSEA Plot
  output$gsea_plot <- renderPlot({
    req(gsea_results())
    plot_gsea(gsea_results(), "data")
  }, height = 500, width = 800)
  
  # Render GSEA Table
  output$gsea_table <- renderDT({
    req(gsea_results())
    datatable(gsea_results(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Render Data Table
  output$table <- renderDT({
    req(dataset())
    datatable(dataset(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Render Volcano Plot
  output$volcano_plot <- renderPlot({
    req(volcano_results())
    volcano_results()
  })
  
  # Download Dataset
  output$download_data <- downloadHandler(
    filename = function() { paste("dataset", Sys.Date(), ".csv", sep="") },
    content = function(file) {
      write.csv(dataset(), file, row.names = FALSE)
    }
  )
  
  # Download GSEA Results
  output$download_gsea <- downloadHandler(
    filename = function() { paste("GSEA_results", Sys.Date(), ".csv", sep="") },
    content = function(file) {
      write.csv(gsea_results(), file, row.names = FALSE)
    }
  )
  # Observe file input change and reset selections when file is uploaded
  observeEvent(input$file, {
    req(input$file)  # Wait for file input
    
    # Reset the selections after file is uploaded, setting them to default values
    col_names <- colnames(dataset())
    
    # Automatically select relevant columns for new file
    default_sort_col <- select_best_match(c("stat", "t"), col_names, fallback = col_names[2])
    default_logFC_col <- select_best_match(c("log2FoldChange", "logFC"), col_names)
    default_pval_col <- select_best_match(c("pvalue", "P.Value", "pval"), col_names)
    
    # Update selectInputs with default values
    updateSelectInput(session, "id_col", choices = col_names, selected = col_names[1])
    updateSelectInput(session, "sort_col", choices = col_names, selected = default_sort_col)
    updateSelectInput(session, "logFC_col", choices = col_names, selected = default_logFC_col)
    updateSelectInput(session, "pval_col", choices = col_names, selected = default_pval_col)
  })
}

# Run App
shinyApp(ui = ui, server = server)






























