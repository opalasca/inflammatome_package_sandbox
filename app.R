

# Install packages (if not already present) and Load libraries ---------------------------------------------------------------
list.of.packages <- c("shiny","DT","readr","shinyjs")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

options(shiny.maxRequestSize = 100 * 1024^2)

source("functions.R")



# Define UI
ui <- fluidPage(
  titlePanel("GSEA and Volcano Plot Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload File", accept = c(".tsv", ".txt")),
      selectInput("id_col", "Select ID Column", choices = NULL),
      selectInput("keytype", "Select Keytype", choices = c("ensembl", "uniprot", "entrez", "refseq", "symbol")),
      selectInput("sort_col", "Select Sorting Column", choices = NULL),
      actionButton("run_gsea", "Run GSEA"),
      actionButton("run_volcano", "Run Volcano Plot"),
      uiOutput("volcano_inputs")
    ),
    
    mainPanel(
      plotOutput("gsea_plot"),
      DTOutput("table"),
      plotOutput("volcano_plot")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive block to read the dataset (without rownames)
  dataset <- reactive({
    req(input$file)  # Ensure file is uploaded
    data <- read.delim(input$file$datapath, header = TRUE, sep = "\t")
    data$ID <- data[[1]]  # Set the first column as ID column
    updateSelectInput(session, "id_col", choices = colnames(data), selected = colnames(data)[1])
    updateSelectInput(session, "sort_col", choices = colnames(data), selected = colnames(data)[2])
    return(data)
  })
  
  # Validate the ID column exists in the dataset
  validate_dataset <- function(data, id_col) {
    if (!(id_col %in% colnames(data))) {
      showNotification("Selected ID column does not exist in the dataset!", type = "error")
      return(FALSE)  # Invalid dataset
    }
    return(TRUE)  # Valid dataset
  }
  
  # Reactive block to process the data once
  processed_data <- reactive({
    req(dataset(), input$id_col, input$keytype)  # Ensure the dataset and selections are available
    data <- dataset()  # Get the raw dataset
    id_col <- input$id_col  # Get the selected ID column
    
    # Validate the ID column exists
    if (!validate_dataset(data, id_col)) return(NULL)
    
    # Process the data for GSEA and volcano
    processed_data <- process_input_data(data, id = id_col, keytype = tolower(input$keytype))  # Convert keytype to lowercase
    
    return(processed_data)  # Return processed data
  })
  
  # Run GSEA analysis
  gsea_results <- eventReactive(input$run_gsea, {
    req(processed_data(), input$sort_col)
    
    data <- processed_data()  # Use the processed data
    id_col <- input$id_col  # Get the selected ID column
    
    # Ensure valid ID keys for the keytype
    valid_keys <- keys(org.Hs.eg.db, keytype = tolower(input$keytype))  # Convert to lowercase keytype
    if (!any(data[[id_col]] %in% valid_keys)) {
      showNotification("Invalid ID keys for the selected keytype!", type = "error")
      return(NULL)
    }
    
    gene_sets <- get_gene_sets(input$keytype)  # Fetch gene sets based on the selected keytype
    gsea_analysis(data, gene_sets, sorting_value_col_name = input$sort_col, name = "data")
  })
  
  # Run Volcano plot
  volcano_results <- eventReactive(input$run_volcano, {
    req(processed_data(), input$logFC_col, input$pval_col)
    
    data <- processed_data()  # Use the processed data
    id_col <- input$id_col  # Get the selected ID column
    
    plot_volcano(data, logFC_col_name = input$logFC_col, pval_col_name = input$pval_col, name = "volcano")
  })
  
  # Render GSEA plot
  output$gsea_plot <- renderPlot({
    req(gsea_results())  # Ensure GSEA results are available
    gsea_results()
  })
  
  # Render volcano plot
  output$volcano_plot <- renderPlot({
    req(volcano_results())  # Ensure volcano results are available
    volcano_results()
  })
  
  # Render data table
  output$table <- renderDT({
    req(dataset())  # Ensure dataset is available
    datatable(dataset())
  })
  
  # Dynamically render volcano plot inputs
  output$volcano_inputs <- renderUI({
    req(dataset())
    list(
      selectInput("logFC_col", "Select LogFC Column", choices = colnames(dataset())),
      selectInput("pval_col", "Select P-Value Column", choices = colnames(dataset()))
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

