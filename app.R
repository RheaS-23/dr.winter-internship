#Himani' shiny 

# Load packages ----------------------------------------------------------------
library(shiny)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(shinythemes)
# Load data --------------------------------------------------------------------

training_merged_obj <- readRDS("~/Downloads/training_merged_obj.rds")

input_genes <- list()

# Clustering -------------------------------------------------------------------

#UMAP by Annotation
DimPlot(training_merged_obj, reduction = "umap", label = TRUE)
# Original
DimPlot(training_merged_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

# Define UI --------------------------------------------------------------------

ui <- fluidPage(theme = shinytheme("flatly"),
                titlePanel("Visualizing Seurat Object Using Shiny"),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      inputId = "reduction",
                      label = "Reduction Choice:",
                      choices = c(
                        "UMAP" = "umap_",
                        "t-SNE" = "tsne_"
                      ),
                      selected = "umap_"
                    ),
                    selectInput(
                      inputId = "plot",
                      label = "Group By:",
                      choices = c(
                        "Clusters" = "seurat_clusters",
                        "Annotation" = "annotation.cluster",
                        "Location/Sample of Origin" = "location",
                        "Condition/Disease" = "condition"
                      ),
                      selected = "annotation.cluster"
                    ),
                    #allowing the user to subset the amount of cells
                    sliderInput(
                      inputId = "subset",
                      label = "Select an amount of cells to subset:",
                      min = 0,
                      max = 33000,
                      value = 5000
                    ),
                    selectInput(
                      inputId ="Binary_Variable",
                      label = "Select a Binary Variable to Split:",
                      choices = c("None", "nFeature_RNA", "nCount_RNA", "percent.mt"),
                      selected = "None"
                    ),
                    numericInput(
                      inputId = "cut_off",
                      label = "Enter a number to split the variable:",
                      value = 5
                    ),
                    #add option to select multiple textboxes
                    numericInput(
                      inputId = "numTextBoxes",
                      label = "Number of Genes to View Expression:",
                      value = 1,
                      min = 1),
                    actionButton("addButton", "Add Text Boxes"),
                    uiOutput("dynamicTextBoxes"),
                    textInput(
                      inputId = "single_gene_summary",
                      label = "Single Gene Summary Choice:",
                      value = "HLA-DRA"
                    ), 
                    selectInput(
                      inputId = "dataset",
                      label = "Dataset:",
                      choices = c(
                        "Training" = "training_",
                        "Testing" = "testing_",
                        "Validation" = "validation_"
                      ),
                      selected = "validation_"
                    ),
                    # Model Visualization - Logistic Regression
                    selectInput("data_logreg",
                                label = "Choose a data set to test model (logistic regression)",
                                choices = c("30% Held-out Data" = "test_set",
                                            "Independent Samples" = "ind_set"),
                                selected = "30% Held-Out Data")
                  ),
                  
                  mainPanel(
                    tabsetPanel(
                      type = "tabs",
                      tabPanel("Data Summary", imageOutput("selected_plot"), plotOutput("gene_expression_group.by")
                               ,plotOutput("single_gene_expression"), textOutput("gene_summary")),
                      tabPanel("ML Model Visualization", imageOutput("selected_visual"), imageOutput("selected_visual2")),
                      tabPanel("ML Model LogReg", imageOutput(outputId = "logreg_ml_model"))
                    )
                  )
                )
)



# Define server ----------------------------------------------------------------

server <- function(input, output) {
  output$selected_plot <- renderImage({
    filename <- normalizePath(file.path('./images',
                                        paste(input$reduction, input$plot, '.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename,
         width = 650,
         height = 400,
         alt = "plot")
  }, deleteFile = FALSE)
  
  #viewing the data by group.by violin plot
  output$gene_expression_group.by <- renderPlot({
    if (input$Binary_Variable == "None"){
      new_cells <- subset(training_merged_obj, downsample = input$subset)
      VlnPlot(new_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = input$plot, pt.size = 0)
    }
    else{
      training_merged_obj@meta.data[["binary"]] <- ifelse(training_merged_obj@meta.data[[input$Binary_Variable]] > input$cut_off, "Above Threshold", "Below Threshold")
      new_cells <- subset(training_merged_obj, downsample = input$subset)
      VlnPlot(new_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "binary", pt.size = 0)
    }
  })
  #making the amount of textboxes the user wants to appear
  observeEvent(input$addButton, {
    output$dynamicTextBoxes <- renderUI({
      numTextBoxes <- input$numTextBoxes
      textBoxes <- lapply(1:numTextBoxes, function(i) {
        textInput(inputId = paste0("Gene", i), label = paste0("Gene ", i), value = "HLA-DRA")
      })
      do.call(tagList, textBoxes)
    })
  })
  
  observeEvent(input$addButton, {
    numTextBoxes <- input$numTextBoxes
    
    output$single_gene_expression <- renderPlot({
      input_genes <- lapply(1:numTextBoxes, function(i) {
        inputId <- paste0("Gene", i)
        input[[inputId]]
      })
      print(input_genes)
      
      vec <- rep(0, times = 33000)
      training_merged_obj@meta.data[["all.cells"]] <- vec
      VlnPlot(training_merged_obj, features = input_genes, group.by = "all.cells", pt.size = 0)
      # You may need to adjust the remaining code to work with your specific data and plot functions
    })
  })
  output$gene_summary <- renderText({
    exp_val <- FetchData(training_merged_obj, vars = input$single_gene_summary)
    exp_val <- exp_val %>%
      rename("gene" = input$single_gene_summary)
    paste("Summary of Expression of Selected Gene: ", 
          "   Minimum: ", min(exp_val$gene),
          ",   Maximum: ", max(exp_val$gene),
          ",   Median: ", median(exp_val$gene),
          ",   Mean: ", mean(exp_val$gene), sep = "")
  })
  output$selected_visual <- renderImage({
    filename <- normalizePath(file.path('./images',
                                        paste(input$dataset,'conf.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename,
         width = 450,
         height = 275,
         alt = "plot")
  }, deleteFile = FALSE)
  output$selected_visual2 <- renderImage({
    filename <- normalizePath(file.path('./images',
                                        paste(input$dataset,'bar.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename,
         width = 450,
         height = 275,
         alt = "plot")
  }, deleteFile = FALSE)
  
  # ML Model - LogReg
  output$logreg_ml_model_ <- renderImage({
    filename <- normalizePath(file.path('~/Documents/GitHub/dr.winter-internship/images',
                                        paste(input$data_logreg,'_allcells.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename,
         width = 450,
         height = 275,
         alt = "plot")
  }, deleteFile = FALSE)
  
}

# Create the Shiny app object --------------------------------------------------

shinyApp(ui = ui, server = server)