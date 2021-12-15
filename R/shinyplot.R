
shinyplot<-function(out, measuredData = NULL, SpParams = NULL) {
  type_out = class(out)[1] #growth or spwb
  dates_out = as.Date(names(out$subdaily))
  if(type_out=="spwb") {
    transpirationMode = out$spwbInput$control$transpirationMode
    subdaily_out = out$spwbInput$control$subdailyResults
    cohorts_out = row.names(out$spwbInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$spwbInput$cohorts), 
                            " (",out$spwbInput$cohorts$Name, ")")
  } else {
    transpirationMode = out$growthInput$control$transpirationMode
    cohorts_out = row.names(out$growthInput$cohorts)
    cohorts_sp_out = paste0(row.names(out$growthInput$cohorts), 
                            " (",out$growthInput$cohorts$Name, ")")
    subdaily_out = out$growthInput$control$subdailyResults
  }
  cohort_choices = cohorts_out
  names(cohort_choices) = cohorts_sp_out
  
  plot_main_choices = c("Water balance", "Soil", "Plants")
  if(transpirationMode=="Sperry") {
    plot_main_choices = c(plot_main_choices, "Energy balance")
  }
  if(type_out=="growth") {
    plot_main_choices = c(plot_main_choices, "Carbon balance", "Growth")
  }
  wb_plot_choices = c("PET & Precipitation" = "PET_Precipitation",
                      "PET and Net rain" = "PET_NetRain",
                      "Snow" = "Snow",
                      "Water exported" = "Export",
                      "Evapotranspiration" = "Evapotranspiration")
  soil_plot_choices = c(
    "Soil water potential" = "SoilPsi",
    "Soil relative water content" = "SoilRWC",
    "Soil moisture (m3/m3) content" = "SoilTheta",
    "Soil volume (mm) content" = "SoilVol",
    "Water table depth" = "WTD",
    "Plant extraction from soil"= "PlantExtraction")
  if(transpirationMode=="Sperry") {
    soil_plot_choices <-c(soil_plot_choices,
                          "Hydraulic redistribution" = "HydraulicRedistribution"
                          )
  }
  subdaily_soil_plot_choices = c("Plant extraction from soil" = "PlantExtraction")
  
  plant_plot_choices = c("LAI" = "PlantLAI",
                         "Stress" = "PlantStress",
                         "Transpiration" = "PlantTranspiration",
                         "Transpiration per leaf" = "TranspirationPerLeaf",
                         "Gross photosynthesis" = "PlantGrossPhotosynthesis",
                         "Gross photosynthesis per leaf" = "GrossPhotosynthesisPerLeaf")
  if(transpirationMode=="Sperry") {
    plant_plot_choices <-c(plant_plot_choices,
                           "Water use efficiency" = "PlantWUE",
                           "Soil-plant conductance" = "SoilPlantConductance",
                           "Minimum leaf water potential" = "LeafPsiMin",
                           "Maximum leaf water potential" = "LeafPsiMax",
                           "Leaf water potential range" = "LeafPsiRange",
                           "Midday upper stem water potential" = "StemPsi",
                           "Midday root crown water potential" = "RootPsi",
                           "Stem PLC" = "StemPLC",
                           "Stem relative water content" = "StemRWC",
                           "Leaf relative water content" = "LeafRWC",
                           "Plant water balance" = "PlantWaterBalance",
                           "Absorbed short-wave radiation" = "PlantAbsorbedSWR",
                           "Absorbed short-wave radiation per leaf" = "AbsorbedSWRPerLeaf",
                           "Net long-wave radiation" = "PlantNetLWR",
                           "Net long-wave radiation per leaf" = "NetLWRPerLeaf"
                           )
  }
  subdaily_plant_plot_choices = c("Plant transpiration" = "PlantTranspiration",
                                  "Plant water balance" = "PlantWaterBalance",
                                  "Plant gross photosynthesis" = "PlantGrossPhotosynthesis",
                                  "Plant net photosynthesis" = "PlantNetPhotosynthesis",
                                  "Soil-plant conductance" = "SoilPlantConductance",
                                  "Leaf water potential" = "LeafPsi",
                                  "Root crown water potential" = "RootPsi",
                                  "Upper stem water potential" = "StemPsi",
                                  "Stem relative water content" = "StemRWC",
                                  "Leaf relative water content" = "LeafRWC",
                                  "Absorbed short-wave radiation" = "PlantAbsorbedSWR",
                                  "Leaf absorbed SWR" = "LeafAbsorbedSWR",
                                  "Leaf net LWR" = "LeafNetLWR",
                                  "Leaf transpiration" = "LeafTranspiration",
                                  "Leaf gross photosynthesis" = "LeafGrossPhotosynthesis",
                                  "Leaf net photosynthesis" = "LeafNetPhotosynthesis"
                                  )
  carbon_plot_choices = c("Gross photosynthesis per dry" = "GrossPhotosynthesis",
                          "Maintenance respiration per dry" = "MaintenanceRespiration",
                          "Carbon balance per dry" = "CarbonBalance",
                          "Leaf sugar concentration" = "SugarLeaf",
                          "Leaf starch concentration" = "StarchLeaf",
                          "Sapwood sugar concentration" = "SugarSapwood",
                          "Sapwood starch concentration" = "StarchSapwood",
                          "Sugar transport" = "SugarTransport",
                          "Root exudation" = "RootExudation",
                          "Leaf osmotic potential at full turgor" = "LeafPI0",
                          "Stem osmotic potential at full turgor" = "StemPI0")
  
  growth_plot_choices = c("Leaf area" = "LeafArea",
                          "Sapwood area" = "SapwoodArea",
                          "Fine root area" = "FineRootArea",
                          "Leaf biomass per individual" = "LeafBiomass",
                          "Sapwood biomass per individual" = "SapwoodBiomass",
                          "Fine root biomass per individual" = "FineRootBiomass",
                          "Labile carbon biomass per individual" = "LabileBiomass",
                          "Total dry biomass per individual" = "TotalLivingBiomass",
                          "Sapwood area growth" = "SAgrowth",
                          "Leaf area growth" = "LAgrowth",
                          "Fine root area growth" = "FRAgrowth",
                          "Sapwood area / Leaf area" = "HuberValue",
                          "Fine root area / Leaf area" = "RootAreaLeafArea")
  subdaily_carbon_plot_choices = c("Gross photosynthesis per dry" = "GrossPhotosynthesis")
  energy_plot_choices = c("Above-canopy air temperature" = "AirTemperature",
                          "Within-canopy air temperature" = "CanopyTemperature",
                          "Soil surface temperature" = "SoilTemperature",
                          "Canopy energy balance components" = "CanopyEnergyBalance",
                          "Soil energy balance components" = "SoilEnergyBalance")
  subdaily_energy_plot_choices = c("Air/canopy/soil temperature" ="Temperature")
  # Define UI for application that draws a histogram
  results <- tabPanel("Results",
                      # Sidebar with a slider input for number of bins 
                      sidebarLayout(
                        sidebarPanel(
                          tabsetPanel(
                            tabPanel("Plot type & period",
                                checkboxInput(
                                    inputId = "subdaily_check",
                                    label = "Subdaily plots",
                                    value = FALSE
                                ),                          
                                selectInput(
                                    inputId = "plot_main_type",
                                    label = "Plot category", 
                                    choices = plot_main_choices,
                                    selected = c("Water balance")
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
                                )
                            ),
                            tabPanel("Options",
                                     selectInput(
                                       inputId = "cohort_selection",
                                       label = "Plant cohorts",
                                       choices = cohort_choices,
                                       selected = cohort_choices,
                                       multiple = TRUE,
                                       selectize = FALSE
                                     ),
                                     selectInput(
                                       inputId = "summary_type",
                                       label = "Temporal aggregation", 
                                       choices = c("None" = "day", 
                                                   "By weeks" = "week", 
                                                   "By months" = "month", 
                                                   "By seasons" = "quarter", 
                                                   "By years" = "year"),
                                       selected = c("None")
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
      eval_choices[["Fuel moisture content"]] = "FMC"
    }
    if(any(c(paste0("PD_", cohorts_out) %in% names(measuredData),
             paste0("MD_", cohorts_out) %in% names(measuredData)))) {
      eval_choices[["Leaf water potential"]] = "WP"
    }
    if(type_out=="growth" && any(paste0("BAI_", cohorts_out) %in% names(measuredData))) {
      eval_choices[["Basal area increment"]] = "BAI"
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
        sel <- plot_main_choices %in% c("Plants", "Carbon balance", "Soil", "Energy balance")
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
      else if(main_plot=="Plants") sub_choices = plant_plot_choices
      else if(main_plot=="Carbon balance") sub_choices = carbon_plot_choices
      else if(main_plot=="Energy balance") sub_choices = energy_plot_choices
      else if(main_plot=="Growth") sub_choices = growth_plot_choices
      else sub_choices = soil_plot_choices

      if(input$subdaily_check && subdaily_out) {
        if(main_plot=="Plants") sub_choices = subdaily_plant_plot_choices
        else if(main_plot=="Carbon balance")  sub_choices = subdaily_carbon_plot_choices
        else if(main_plot=="Energy balance") sub_choices = subdaily_energy_plot_choices
        else sub_choices = subdaily_soil_plot_choices
      }
      updateSelectInput(session, "plot_type",
                        choices = sub_choices)
    })
    output$results_plot <- renderPlot({
      date_lim = input$date_range
      date_range = dates_out[dates_out >= date_lim[1] & dates_out <= date_lim[2]]
      plot(out, type = input$plot_type, dates = date_range,
           summary.freq = input$summary_type, 
           subdaily = input$subdaily_check,
           cohorts = (cohort_choices %in% input$cohort_selection),
           bySpecies = input$byspecies_check)
    })
    if(!is.null(measuredData)) {
      observe({
        eval_plot <- input$eval_type
        sub_choices <- NULL
        if(eval_plot %in% c("E", "BAI", "FMC")) {
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
                        SpParams = SpParams,
                        plotType = "dynamic")
      )
      output$scatter_eval_plot <- renderPlot(
        evaluation_plot(out, measuredData, 
                        type = input$eval_type,
                        temporalResolution = input$eval_summary_type,
                        cohort = input$eval_cohort,
                        SpParams = SpParams,
                        plotType = "scatter")
      )
      output$eval_stats <- renderPrint({
        evaluation_stats(out, measuredData, 
                         type = input$eval_type,
                         temporalResolution = input$eval_summary_type,
                         cohort = input$eval_cohort,
                         SpParams = SpParams)
      })
    }
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}