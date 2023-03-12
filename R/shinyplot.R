.shinyplot_day<-function(out){
  if(inherits(out, c("spwb_day", "pwb_day"))) {
    transpirationMode = out$spwbInput$control$transpirationMode
    cohorts_out = row.names(out$spwbInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$spwbInput$cohorts), 
                            " (",out$spwbInput$cohorts$Name, ")")
  } else {
    transpirationMode = out$growthInput$control$transpirationMode
    cohorts_out = row.names(out$growthInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$growthInput$cohorts), 
                            " (",out$growthInput$cohorts$Name, ")")
  }
  plot_main_choices = c("Soil", "Plants", "Sunlit/Shade", "Energy balance")
  if(inherits(out, c("growth_day"))) {
    plot_main_choices = c(plot_main_choices, 
                          "Labile carbon balance")
  }
  subdaily_soil_plot_choices = .getSubdailySoilPlotTypes()
  subdaily_plant_plot_choices = .getSubdailyPlantPlotTypes()
  subdaily_sunlitshade_plot_choices = .getSubdailySunlitShadePlotTypes()
  subdaily_energy_plot_choices = .getSubdailyEnergyBalancePlotTypes()
  subdaily_labile_plot_choices = .getSubdailyLabilePlotTypes()
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          inputId = "plot_main_type",
          label = "Category",
          choices = plot_main_choices,
          selected = plot_main_choices[1]
        ),
        selectInput(
          inputId = "plot_var",
          label = "Variable", 
          choices = subdaily_soil_plot_choices,
          selected = subdaily_soil_plot_choices[1]
        ),
        checkboxInput(
          inputId = "byspecies_check",
          label = "Aggregation by species",
          value = FALSE
        )
      ),
      mainPanel(
        plotOutput("spatial_plot")
      )
    )
  )
  server <- function(input, output, session) {
    observe({
      main_plot <- input$plot_main_type
      if(main_plot=="Soil") sub_choices = subdaily_soil_plot_choices
      else if(main_plot=="Plants") sub_choices = subdaily_plant_plot_choices
      else if(main_plot=="Sunlit/Shade") sub_choices = subdaily_sunlitshade_plot_choices
      else if(main_plot=="Energy balance") sub_choices = subdaily_energy_plot_choices
      else if(main_plot=="Labile carbon balance") sub_choices = subdaily_labile_plot_choices
      updateSelectInput(session, "plot_var",
                        choices = sub_choices)
    })
    output$spatial_plot <- renderPlot({
      plot(out, type = input$plot_var, bySpecies = input$byspecies_check)
    })
  }
  
  shinyApp(ui = ui, server = server)
}

.shinyplot_sim<-function(out, measuredData = NULL) {
  if(inherits(out, c("spwb_day", "pwb_day", "growth_day"))) return(.shinyplot_day(out))
  if(!inherits(out, c("growth", "pwb" , "spwb", "fordyn"))) stop("Wrong class for 'out'. Should either be 'spwb', 'growth' or 'fordyn'.")
  type_out = class(out)[1] #growth, spwb, fordyn
  if(type_out=="spwb") {
    transpirationMode = out$spwbInput$control$transpirationMode
    subdaily_out = out$spwbInput$control$subdailyResults
    cohorts_out = row.names(out$spwbInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$spwbInput$cohorts), 
                            " (",out$spwbInput$cohorts$Name, ")")
    dates_out = as.Date(row.names(out$WaterBalance))
  } else if(type_out=="growth") {
    transpirationMode = out$growthInput$control$transpirationMode
    cohorts_out = row.names(out$growthInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$growthInput$cohorts), 
                            " (",out$growthInput$cohorts$Name, ")")
    subdaily_out = out$growthInput$control$subdailyResults
    dates_out = as.Date(row.names(out$WaterBalance))
  } else { # fordyn
    out_1 = out$GrowthResults[[1]]
    transpirationMode = out_1$growthInput$control$transpirationMode
    cohorts_out = row.names(out_1$growthInput$cohorts)
    sp_out = out_1$growthInput$cohorts$Name
    if(length(out$GrowthResults)>1) {
      for(i in 2:length(out$GrowthResults)) {
        out_i = out$GrowthResults[[i]]
        cn = row.names(out_i$growthInput$cohorts)
        for(j in 1:length(cn)) {
          if(!(cn[j] %in% cohorts_out)) {
            cohorts_out = c(cohorts_out, cn[j])
            sp_out = c(sp_out, out_i$growthInput$cohorts$Name[j])
          }
        }
      }
    }
    cohorts_sp_out = paste0(cohorts_out, 
                            " (",sp_out, ")")
    subdaily_out = FALSE
    dates_out = NULL
    for(i in 1:length(out$GrowthResults)) {
      if(is.null(dates_out)) dates_out = as.Date(row.names(out$GrowthResults[[i]]$WaterBalance))
      else dates_out = c(dates_out, as.Date(row.names(out$GrowthResults[[i]]$WaterBalance)))
    }
  }
  cohort_choices = cohorts_out
  names(cohort_choices) = cohorts_sp_out
  
  wb_plot_choices = .getWaterBalancePlotTypes()
  stand_plot_choices = .getStandPlotTypes(type_out)
  soil_plot_choices = .getSoilPlotTypes(type_out, transpirationMode) 
  plant_plot_choices = .getPlantPlotTypes(transpirationMode)
  sunlitshade_plot_choices = .getSunlitShadePlotTypes(transpirationMode)
  energy_plot_choices = .getEnergyPlotTypes(transpirationMode)
  labile_plot_choices = .getLabileGROWTHPlotTypes(transpirationMode)
  plant_balance_plot_choices = .getCohortBiomassBalanceGROWTHPlotTypes(transpirationMode)
  plant_structure_plot_choices = .getStructuralGROWTHPlotTypes(transpirationMode)
  plant_growthmortality_plot_choices = .getGrowthMortalityGROWTHPlotTypes(transpirationMode)
  forest_dynamics_plot_choices = .getUniqueFORDYNPlotTypes(transpirationMode)
  
  plot_main_choices = c("Water balance", "Soil", "Stand", "Plants")
  if(transpirationMode %in% c("Sperry","Cochard")) {
    plot_main_choices = c(plot_main_choices,"Sunlit/Shade", "Energy balance")
  }
  if(type_out %in% c("growth", "fordyn")) {
    plot_main_choices = c(plot_main_choices, 
                          "Labile carbon balance",
                          "Biomass balance",
                          "Plant structure",
                          "Growth & mortality")
  }
  if(type_out=="fordyn") {
    plot_main_choices = c(plot_main_choices, "Forest structure & composition")
  }
  
  subdaily_soil_plot_choices = .getSubdailySoilPlotTypes()
  subdaily_plant_plot_choices = .getSubdailyPlantPlotTypes()
  subdaily_sunlitshade_plot_choices = .getSubdailySunlitShadePlotTypes()
  subdaily_labile_plot_choices = .getSubdailyLabilePlotTypes()
  subdaily_energy_plot_choices = .getSubdailyEnergyBalancePlotTypes()
  
  
  # Define UI for application that draws a histogram
  results <- tabPanel("Results",
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            tabPanel("Plot selection",
                                     checkboxInput(
                                       inputId = "subdaily_check",
                                       label = "Subdaily plots",
                                       value = FALSE
                                     ),                          
                                     selectInput(
                                       inputId = "plot_main_type",
                                       label = "Plot category", 
                                       choices = plot_main_choices,
                                       selected = "Water balance"
                                     ),
                                     selectInput(
                                       inputId = "plot_type",
                                       label = "Plot type", 
                                       choices = wb_plot_choices,
                                       selected = wb_plot_choices[1]
                                     ),
                                     sliderInput(
                                       inputId = "date_range",
                                       label = "Date range",
                                       value = c(dates_out[1],dates_out[length(dates_out)]),
                                       min = dates_out[1],
                                       max = dates_out[length(dates_out)]
                                     ),
                                     selectInput(
                                       inputId = "summary_type",
                                       label = "Temporal aggregation", 
                                       choices = c("None" = "day", 
                                                   "By weeks" = "week", 
                                                   "By months" = "month", 
                                                   "By seasons" = "quarter", 
                                                   "By years" = "year"),
                                       selected = ifelse(type_out=="fordyn","year","None")
                                     )
                            ),
                            tabPanel("Plants",
                                     radioButtons(inputId = "plant_group_selection",
                                                  label = "Plant group",
                                                  choices = c("all", "trees", "shrubs"),
                                                  selected = "all",
                                                  inline = TRUE),
                                     selectInput(
                                       inputId = "cohort_selection",
                                       label = "Plant cohorts",
                                       choices = cohort_choices,
                                       selected = cohort_choices,
                                       multiple = TRUE,
                                       selectize = FALSE
                                     ),
                                     checkboxInput(
                                       inputId = "byspecies_check",
                                       label = "Aggregation by species",
                                       value = FALSE
                                     )
                            )
                          )
                        ),
                        
                        # Show a plot of the generated distribution
                        mainPanel(
                          plotOutput("results_plot")
                        )
                      )
  )
  eval_choices = c()
  if(!is.null(measuredData)) {
    if("SWC" %in% names(measuredData)) {
      eval_choices[["Soil water content"]] = "SWC"
      eval_choices[["Relative soil water content"]] = "REW"
    } 
    if("ETR" %in% names(measuredData)) {
      eval_choices[["Total evapotranspiration"]] = "ETR"
      eval_choices[["Total ET against modelled SE+TR"]] = "SE+TR"
    }
    if(any(paste0("E_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Transpiration per leaf area"]] = "E"
    }
    if(any(paste0("FMC_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Live fuel moisture content"]] = "LFMC"
    }
    if(any(c(paste0("PD_", cohorts_out) %in% names(measuredData),
             paste0("MD_", cohorts_out) %in% names(measuredData)))) {
      eval_choices[["Leaf water potential"]] = "WP"
    }
    if(type_out=="growth" && any(paste0("BAI_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Basal area increment"]] = "BAI"
    }
    if(type_out=="growth" && any(paste0("DI_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Diameter increment"]] = "DI"
    }
    if(type_out=="growth" && any(paste0("DBH_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Diameter at breast height"]] = "DBH"
    }
    if(type_out=="growth" && any(paste0("Height_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Plant height"]] = "Height"
    }
  }
  evaluation <- tabPanel("Evaluation",
                         sidebarLayout(
                           sidebarPanel(
                             selectInput(
                               inputId = "eval_type",
                               label = "Evaluation type", 
                               choices = eval_choices,
                               selected = eval_choices[1]
                             ),
                             selectInput(
                               inputId = "eval_cohort",
                               label = "Plant cohort",
                               choices = c()
                             ),
                             selectInput(
                               inputId = "eval_summary_type",
                               label = "Temporal aggregation", 
                               choices = c("By days" = "day", 
                                           "By weeks" = "week", 
                                           "By months" = "month", 
                                           "By years" = "year"),
                               selected = c("None")
                             )
                           ),
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Dynamic",
                                        plotOutput("dynamic_eval_plot")
                               ),
                               tabPanel("Scatter",
                                        plotOutput("scatter_eval_plot")
                               ), 
                               tabPanel("Stats",
                                        verbatimTextOutput("eval_stats")
                               )
                             )
                           )
                         ),
  )
  if(is.null(measuredData)) {
    ui <- navbarPage("Interactive plots", results)
  } else {
    ui <- navbarPage("Interactive plots", results, evaluation)
  }
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    observe({
      if(input$subdaily_check && subdaily_out)  {
        sel <- plot_main_choices %in% c("Plants", "Sunlit/Shade","Labile carbon balance", "Soil", "Energy balance")
        updateSelectInput(session, "plot_main_type",
                          choices = plot_main_choices[sel])
      } else {
        updateSelectInput(session, "plot_main_type",
                          choices = plot_main_choices)
      }
    })
    observe({
      main_plot <- input$plot_main_type
      if(main_plot=="Water balance") sub_choices = wb_plot_choices
      else if(main_plot=="Stand") sub_choices = stand_plot_choices
      else if(main_plot=="Plants") sub_choices = plant_plot_choices
      else if(main_plot=="Sunlit/Shade") sub_choices = sunlitshade_plot_choices
      else if(main_plot=="Labile carbon balance") sub_choices = labile_plot_choices
      else if(main_plot=="Biomass balance") sub_choices = plant_balance_plot_choices
      else if(main_plot=="Energy balance") sub_choices = energy_plot_choices
      else if(main_plot=="Plant structure") sub_choices = plant_structure_plot_choices
      else if(main_plot=="Growth & mortality") sub_choices = plant_growthmortality_plot_choices
      else if(main_plot=="Forest structure & composition") sub_choices = forest_dynamics_plot_choices
      else sub_choices = soil_plot_choices
      
      if(input$subdaily_check && subdaily_out) {
        if(main_plot=="Plants") sub_choices = subdaily_plant_plot_choices
        else if(main_plot=="Sunlit/Shade") sub_choices = subdaily_sunlitshade_plot_choices
        else if(main_plot=="Labile carbon balance")  sub_choices = subdaily_labile_plot_choices
        else if(main_plot=="Energy balance") sub_choices = subdaily_energy_plot_choices
        else sub_choices = subdaily_soil_plot_choices
      }
      updateSelectInput(session, "plot_type",
                        choices = sub_choices)
    })
    observe({
      plant_group = input$plant_group_selection
      if(plant_group=="trees") {
        sel = (substr(cohorts_out,1,1)=="T")
      } else if(plant_group=="shrubs") {
        sel = (substr(cohorts_out,1,1)=="S")
      } else {
        sel = rep(TRUE, length(cohorts_out))
      }
      cohort_choices = cohorts_out[sel]
      names(cohort_choices) = cohorts_sp_out[sel]
      updateSelectInput(session, "cohort_selection",
                        choices = cohort_choices,
                        selected = cohort_choices)
    })
    output$results_plot <- renderPlot({
      date_lim = input$date_range
      date_range = dates_out[dates_out >= date_lim[1] & dates_out <= date_lim[2]]
      plot(out, type = input$plot_type, dates = date_range,
           summary.freq = input$summary_type, 
           subdaily = input$subdaily_check,
           cohorts = cohort_choices[cohort_choices %in% input$cohort_selection],
           bySpecies = input$byspecies_check)
    })
    if(!is.null(measuredData)) {
      observe({
        eval_plot <- input$eval_type
        sub_choices <- NULL
        if(eval_plot %in% c("E", "BAI", "DI", "DBH", "Height", "FMC")) {
          sel = paste0(eval_plot,"_", cohorts_out) %in% names(measuredData)
          sub_choices = cohorts_out[sel]
          names(sub_choices) = cohorts_sp_out[sel]
        }
        else if(eval_plot == "WP") {
          sel = paste0("MD_", cohorts_out) %in% names(measuredData) | paste0("PD_", cohorts_out) %in% names(measuredData)
          sub_choices = cohorts_out[sel]
          names(sub_choices) = cohorts_sp_out[sel]
        }
        updateSelectInput(session, "eval_cohort",
                          choices = sub_choices)
      })
      output$dynamic_eval_plot <- renderPlot(
        evaluation_plot(out, measuredData, 
                        type = input$eval_type,
                        temporalResolution = input$eval_summary_type,
                        cohort = input$eval_cohort,
                        plotType = "dynamic")
      )
      output$scatter_eval_plot <- renderPlot(
        evaluation_plot(out, measuredData, 
                        type = input$eval_type,
                        temporalResolution = input$eval_summary_type,
                        cohort = input$eval_cohort,
                        plotType = "scatter")
      )
      output$eval_stats <- renderPrint({
        evaluation_stats(out, measuredData, 
                         type = input$eval_type,
                         temporalResolution = input$eval_summary_type,
                         cohort = input$eval_cohort)
      })
    }
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}


#' Shiny app with interactive plots
#' 
#' Creates a shiny app with interactive plots for simulation results and evaluation
#' 
#' @param x An object of the right class.
#' @param measuredData A data frame with observed/measured values (see \code{\link{evaluation_plot}}).
#' @param ... Additional parameters.
#' 
#' @details Only run this function in interactive mode. When \code{measuredData} is not \code{NULL}, an additional panel is shown for evaluation plots.
#' 
#' @return An object that represents the shiny app
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso \code{\link{plot.spwb}}, \code{\link{evaluation_plot}}
#' 
#' @name shinyplot
shinyplot<-function(x, ...) {
  UseMethod("shinyplot", x)
}

#' @rdname shinyplot 
shinyplot.growth<-function(x, measuredData = NULL, ...) {
  .shinyplot_sim(x, measuredData = measuredData)
}

#' @rdname shinyplot
shinyplot.spwb<-function(x, measuredData = NULL, ...) {
  .shinyplot_sim(x, measuredData = measuredData)
}

#' @rdname shinyplot
shinyplot.pwb<-function(x, measuredData = NULL, ...) {
  .shinyplot_sim(x, measuredData = measuredData)
}

#' @rdname shinyplot
shinyplot.fordyn<-function(x, measuredData = NULL, ...) {
  .shinyplot_sim(x, measuredData = measuredData)
}

#' @rdname shinyplot
shinyplot.growth_day<-function(x, ...) {
  .shinyplot_day(x)
}

#' @rdname shinyplot
shinyplot.spwb_day<-function(x, ...) {
  .shinyplot_day(x)
}

#' @rdname shinyplot
shinyplot.pwb_day<-function(x, ...) {
  .shinyplot_day(x)
}



