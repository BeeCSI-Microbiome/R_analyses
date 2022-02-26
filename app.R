library(shiny)
library(tidyverse)
library(modules)
library(glue)

widen_results <- use("scripts/widenResults.R")

ui <- fluidPage(
  textInput(inputId = "dataset_name",
            label = "For output naming purposes, what is the name of your dataset?",
            placeholder = "cor_2020"),
  
  # fileInput(inputId = "count_file",
  #           label = "Select the count file produced by Pavian",
  #           accept = c(".csv", ".tsv")),
  # 
  # tableOutput(outputId = "count_table"),
  
  fileInput(inputId = "kraken_files",
            label = "Select the kraken report files",
            accept = c(".txt"),
            multiple = TRUE),
  
  actionButton("aggregate_button", "Aggregate reports", class = "btn-primary"),
  
  tableOutput(outputId = "k_analytic"),
  
  textOutput("widen_results_message")
)

server <- function(input, output, session) {

  k_analytic <- eventReactive(input$aggregate_button, {
    req(input$dataset_name,
        input$kraken_files)
    
    kraken_matrix_dir <- glue("results/{input$dataset_name}/aggregated_kraken_reports")
    
    print(glue("Saving analytic matrices to {kraken_matrix_dir}"))
    
    widen_results$widen_results_function(input$kraken_files$datapath,
                                         input$kraken_files$name,
                                         kraken_matrix_dir)
  })
  
  output$k_analytic <- renderTable({
    req(k_analytic)
    k_analytic()
  })
  
  # ct <- reactive({
  #   req(input$count_file)
  #   
  #   ext <- tools::file_ext(input$count_file$name)
  #   switch(ext,
  #          csv = vroom::vroom(input$count_file$datapath, delim = ","),
  #          tsv = vroom::vroom(input$count_file$datapath, delim = "\t"),
  #          validate("Invalid file; Please upload a .csv or .tsv"))
  # })
  # 
  # output$count_table <- renderTable({
  #   req(ct())
  #   head(ct(), n = 5)
  # })
  
  # output$report_files <- renderTable({
  #   req(krakenReportPaths())
  #   krakenReportPaths()
  # })
}

shinyApp(ui, server)