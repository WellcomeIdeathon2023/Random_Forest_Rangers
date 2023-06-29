library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

ui <- navbarPage(
  title = "Numerical Boxplot App",
  
  tabPanel(
    "Upload Data",
    sidebarLayout(
      sidebarPanel(
        fileInput("files", "Choose CSV files", multiple = TRUE, accept = ".csv")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("File Previews", DTOutput("filePreviews")),
          tabPanel("Boxplots", selectInput("column", "Select a numeric column:", choices = NULL), plotOutput("boxplot"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    req(input$files)
    
    file_data <- lapply(input$files$datapath, read.csv)
    
    # Combine the data frames into a single data frame
    combined_data <- bind_rows(file_data)
    
    # Return the combined data frame
    combined_data
  })
  
  output$filePreviews <- renderDT({
    data_preview <- data()
    
    if (!is.null(data_preview)) {
      datatable(data_preview, options = list(pageLength = 5)) %>%
        formatStyle(names(data_preview), color = "white") %>%
        formatStyle(names(data_preview), color = "white", selector = "th")
    }
  })
  
  observeEvent(data(), {
    dataset <- data()
    numeric_columns <- sapply(dataset, is.numeric)
    updateSelectInput(session, "column", choices = names(dataset)[numeric_columns])
  })
  
  output$boxplot <- renderPlot({
    dataset <- data()
    if (!is.null(input$column)) {
      ggplot(dataset, aes_string(y = input$column)) +
        geom_boxplot() +
        labs(y = input$column, title = "Boxplot")
    }
  })
}

shinyApp(ui, server)
