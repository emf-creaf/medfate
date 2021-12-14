
shinyplot<-function(out) {
  type_out = class(out)[1] #growth or spwb
  dates_out = as.Date(names(out$subdaily))
  if(type_out=="spwb") {
    transpirationMode = out$spwbInput$control$transpirationMode
  } else {
    transpirationMode = out$growthInput$control$transpirationMode
  }
  
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
  plant_plot_choices = c("LAI" = "PlantLAI",
                         "Stress" = "PlantStress",
                         "Transpiration" = "PlantTranspiration",
                         "Transpiration per leaf" = "TranspirationPerLeaf",
                         "Gross photosynthesis" = "PlantGrossPhotosynthesis",
                         "Gross photosynthesis per leaf" = "GrossPhotosynthesisPerLeaf")
  if(transpirationMode=="Sperry") {
    plant_plot_choices <-c(plant_plot_choices,
                           "Soil-plant conductance" = "SoilPlantConductance",
                           "Midday leaf water potential" = "LeafPsi",
                           "Midday upper stem water potential" = "StemPsi",
                           "Midday root crown water potential" = "RootPsi",
                           "Water use efficiency" = "PlantWUE",
                           "Absorbed short-wave radiation" = "PlantAbsorbedSWR",
                           "Absorbed short-wave radiation per leaf" = "AbsorbedSWRPerLeaf",
                           "Net long-wave radiation" = "PlantNetLWR",
                           "Net long-wave radiation per leaf" = "NetLWRPerLeaf"
                           )
  }
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
  
  energy_plot_choices = c("Above-canopy air temperature" = "AirTemperature",
                          "Within-canopy air temperature" = "CanopyTemperature",
                          "Soil surface temperature" = "SoilTemperature",
                          "Canopy energy balance components" = "CanopyEnergyBalance",
                          "Soil energy balance components" = "SoilEnergyBalance")
  
  # Define UI for application that draws a histogram
  ui <- fluidPage(
    
    # Application title
    titlePanel("Interactive result plots"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
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
        # checkboxInput(
        #   inputId = "subdaily_check",
        #   label = "Subdaily",
        #   value = FALSE
        # ),
        checkboxInput(
          inputId = "byspecies_check",
          label = "By species",
          value = FALSE
        )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("medfate_plot")
      )
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    
    observe({
      main_plot <- input$plot_main_type
      if(main_plot=="Water balance") sub_choices = wb_plot_choices
      else if(main_plot=="Plants") sub_choices = plant_plot_choices
      else if(main_plot=="Carbon balance") sub_choices = carbon_plot_choices
      else if(main_plot=="Energy balance") sub_choices = energy_plot_choices
      else if(main_plot=="Growth") sub_choices = growth_plot_choices
      else sub_choices = soil_plot_choices
      updateSelectInput(session, "plot_type",
                        choices = sub_choices)
    })
    output$medfate_plot <- renderPlot({
      date_lim = input$date_range
      date_range = dates_out[dates_out >= date_lim[1] & dates_out <= date_lim[2]]
      plot(out, type = input$plot_type, dates = date_range,
           summary.freq = input$summary_type, 
           # subdaily = input$subdaily_check,
           bySpecies = input$byspecies_check)
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}