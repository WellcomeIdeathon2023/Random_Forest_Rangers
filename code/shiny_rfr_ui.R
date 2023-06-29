library(shiny)
library(shinythemes)
library(vroom)
library(tidyverse)
library(DT)



data_upload <- tabPanel(
  title = "Data Import",
  sidebarLayout(
    sidebarPanel( 
      title = "Upload Options",
      fileInput("files","Choose CSV Files",multiple = TRUE,accept = ".csv")
    ),
        mainPanel(
        tabsetPanel(
          tabPanel("File Previews", style = "background-color: white;", DTOutput("filePreviews"))
        )
      )
    )
  )




box_page <- tabPanel(
  title = "Boxplot Anlaysis",
  sidebarLayout(
    sidebarPanel(
      title = "Box Options"
      tabPanel(selectInput("column", "Select a numeric column:", choices = NULL),
      tabPanel(selectInput("column", "Select a numeric column:", choices = NULL),
               
               
               plotOutput("boxplot")
    )
  )
)

umap_page <- tabPanel((
  title = "UMAP Analysis",
  titlePanel = "UMAP Options",
  sidebarLayout(
    sidebarPanel(
      title = "Inputs"
    )
  )
))


box_page <- tabPanel((
  title = "Boxplot Analysis",
  titlePanel = "BOX Options",
  sidebarLayout(
    sidebarPanel(
      title = "Inputs"
    )
  )
))


clust_page <- 

heatmap_page <- tabPanel((
  title = "Heatmap Analysis",
  titlePanel = "Map Options",
  sidebarLayout(
    sidebarPanel(
      title = "Inputs"
    )
  )
))

scatter_page <- tabPanel((
  title = "Scatterplot Analysis",
  titlePanel = "Plot Options",
  sidebarLayout(
    sidebarPanel(
      title = "Inputs"
    )
  )
))

#Define UI ---
ui <- shiny::navbarPage(
      title = "RFR: Data Analyser",
      theme = shinytheme('superhero'),
      data_upload,
      box_page,
      umap_page,
      heatmap_page,
      scatter_page
)
  




  
  fluidPage(
  
  #Application Title
  titlePanel("RFR:"),
  
  #Application Tabs
  tabsetPanel(
    tabPanel("About"),
    
    tabPanel("Import Data",
             fileInput("file","Data", buttonLabel = "Upload..")),
    tabPanel("Cytodelics Analysis"),
    tabPanel("PCA Explorer"),
    tabPanel("BoxPloter")
             

 )
)



shinyApp(ui = ui, server = server)
