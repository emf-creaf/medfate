.forest_plot_choices<-function(){
  c("LeafAreaDensity", "RootDistribution", "FuelBulkDensity",
    "PARExtinction", "SWRExtinction", "WindExtinction")
}
plot.forest<-function(x, SpParams, type = "LeafAreaDensity", byCohorts = FALSE, bySpecies = FALSE, ...) {
  type = match.arg(type, .forest_plot_choices())
  if(type=="LeafAreaDensity") return(vprofile_leafAreaDensity(x, SpParams = SpParams, draw = TRUE, 
                                                              byCohorts = byCohorts, bySpecies = bySpecies, ...))
  else if(type=="RootDistribution") return(vprofile_rootDistribution(x, SpParams = SpParams,  draw = TRUE, 
                                                                     bySpecies = bySpecies, ...))
  else if(type=="FuelBulkDensity") return(vprofile_fuelBulkDensity(x, SpParams = SpParams,  draw = TRUE, ...))
  else if(type=="PARExtinction") return(vprofile_PARExtinction(x, SpParams = SpParams,  draw = TRUE, ...))
  else if(type=="SWRExtinction") return(vprofile_SWRExtinction(x, SpParams = SpParams,  draw = TRUE, ...))
  else if(type=="WindExtinction") return(vprofile_windExtinction(x, SpParams = SpParams,  draw = TRUE, ...))
}


shinyplot.forest<-function(x, SpParams, ...){
  plot_main_choices = .forest_plot_choices()

  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = "plot_main_type",
          label = "Category",
          choices = plot_main_choices,
          selected = plot_main_choices[1]
        ),
        checkboxInput(
          inputId = "bycohort_check",
          label = "Distinguish cohorts",
          value = FALSE
        ),
        checkboxInput(
          inputId = "byspecies_check",
          label = "Aggregation by species",
          value = FALSE
        )
      ),
      mainPanel(
        plotOutput("forest_plot")
      )
    ),
    h4("Summary"),
    verbatimTextOutput("summary")
  )
  server <- function(input, output, session) {
    output$forest_plot <- renderPlot({
      plot(x, SpParams = SpParams, 
           type = input$plot_main_type, 
           byCohorts = input$bycohort_check,
           bySpecies = input$byspecies_check)
    })
    output$summary <- renderPrint({
      summary(x, SpParams)
    })
  }
  
  shinyApp(ui = ui, server = server)
}
