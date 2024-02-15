.forest_plot_choices<-function(){
  c("LeafAreaDensity", "RootDistribution", "FuelBulkDensity",
    "PARExtinction", "SWRExtinction", "WindExtinction")
}

#' Plot forest attributes
#'
#' Convenient wrappers for vertical forest profiles (see \code{\link{vprofile_leafAreaDensity}}).
#'
#' @param x An object of class \code{\link{forest}}.
#' @param SpParams A data frame with species parameters (see \code{\link{SpParamsMED}}).
#' @param type A string of the plot type: "LeafAreaDensity", "RootDistribution", "FuelBulkDensity", "PARExtinction", "SWRExtinction" or "WindExtinction". 
#' @param byCohorts A logical flag to separate profiles for each cohort.
#' @param bySpecies A logical flag to aggregate results by species.
#' @param includeHerbs A logical flag to include herbaceous layer in the profile.
#' @param \dots Additional parameters to vertical profiles
#' 
#' @return A ggplot or a shiny application, depending on the function.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{vprofile_leafAreaDensity}}
#' 
#' @examples 
#' data(exampleforest)
#' data(SpParamsMED)
#' plot(exampleforest, SpParamsMED)
#' 
#' @name plot.forest
plot.forest<-function(x, SpParams, type = "LeafAreaDensity", byCohorts = FALSE, bySpecies = FALSE, includeHerbs = FALSE, ...) {
  type = match.arg(type, .forest_plot_choices())
  if(type=="LeafAreaDensity") return(vprofile_leafAreaDensity(x, SpParams = SpParams, draw = TRUE, 
                                                              byCohorts = byCohorts, bySpecies = bySpecies,
                                                              includeHerbs = includeHerbs, ...))
  else if(type=="RootDistribution") return(vprofile_rootDistribution(x, SpParams = SpParams,  draw = TRUE, 
                                                                     bySpecies = bySpecies, ...))
  else if(type=="FuelBulkDensity") return(vprofile_fuelBulkDensity(x, SpParams = SpParams,  draw = TRUE, ...))
  else if(type=="PARExtinction") return(vprofile_PARExtinction(x, SpParams = SpParams,  draw = TRUE,
                                                               includeHerbs = includeHerbs, ...))
  else if(type=="SWRExtinction") return(vprofile_SWRExtinction(x, SpParams = SpParams,  draw = TRUE,
                                                               includeHerbs = includeHerbs, ...))
  else if(type=="WindExtinction") return(vprofile_windExtinction(x, SpParams = SpParams,  draw = TRUE,
                                                                 includeHerbs = includeHerbs, ...))
}


#' @rdname plot.forest
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
        ),
        checkboxInput(
          inputId = "includeherbs_check",
          label = "Include herbaceous layer",
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
           bySpecies = input$byspecies_check,
           includeHerbs = input$includeherbs_check)
    })
    output$summary <- renderPrint({
      summary(x, SpParams)
    })
  }
  
  shinyApp(ui = ui, server = server)
}
