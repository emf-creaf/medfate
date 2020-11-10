defaultControl<-function() {
  return(list(
    #For all
    verbose = TRUE,
    subdailyResults = FALSE,
    defaultWindSpeed = 2.5, #m/s
    soilFunctions = "SX",
    
    # For water balance
    snowpack = TRUE,
    rockyLayerDrainage = TRUE,
    unlimitedSoilWater = FALSE,
    plantWaterPools = FALSE,
    leafPhenology = TRUE,
    unfoldingDD = 300,
    transpirationMode = "Granier",
    cavitationRefill = "total",
    refillMaximumRate = 0.05,
    verticalLayerSize = 100,
    
    #spwb with sperry
    gainModifier = 1, 
    costModifier = 1, 
    costWater = "dEdP",
    cochard = FALSE,
    capacitance = FALSE,
    klatstem = 0.01, # stem symplastic-apoplastic lateral conductance
    klatleaf = 0.01, # leaf symplastic-apoplastic lateral conductance
    taper = TRUE,
    numericParams=list(maxNsteps = 400, ntrial = 200, psiTol = 0.0001, ETol = 0.0000001),
    fracLeafResistance = NA,
    fracRootResistance = 0.4,
    averageFracRhizosphereResistance = 0.15,
    Catm = 386,
    ndailysteps = 24,
    thermalCapacityLAI = 1000000,
    
    # growth
    allocationStrategy = "Plant_kmax",
    nonStomatalPhotosynthesisLimitation = TRUE,
    phloemConductanceFactor = 0.2, # phloem conductance per leaf area basis (l*m-2*MPa-1*s-1)
    nonSugarConcentration = 0.25, # mol · l-1
    equilibriumOsmoticConcentration = c(leaf = 0.8, sapwood = 0.6),  # (Paljakka et al. 2017)
    minimumSugarForGrowth = c(leaf = 0.25, sapwood=0.25, fineroot = 0.25),
    # Ogle and Pacala 2010, Tree Physiology 29, 587–605
    respirationRates = c(leaf = 0.00260274, sapwood = 6.849315e-05, fineroot = 0.002054795), # g gluc · g dw -1 · day -1
    turnoverRates = c(sapwood = 0.0001261398, # day-1 Equivalent to annual 4.5% 1-(1-0.045)^(1.0/365)
                      fineroot = 0.001897231), #day-1 Equivalent to annual 50% 1-(1-0.5)^(1.0/365)
    constructionCosts = c(leaf = 1.5, sapwood = 1.47, fineroot = 1.30), #  g gluc · g dw -1
    maximumRelativeGrowthRates = c(leaf = 0.01, # m2 leaf ·cm-2 sapwood· day-1
                                   sapwood = 0.001, # cm2 sapwood ·cm-2 sapwood· day-1
                                   fineroot = 0.1) # g dw · g dw -1 · day -1
    #For forest dynamics
#     freqZopt = 20,
#     sap2tree=TRUE,
#     mergeCohorts=TRUE,
#     removeDeadCohorts=TRUE,
#     initRadius = 5000,
#     setDefaults = TRUE,
#     diamClasses = c(0,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,7.5,10.0, seq(15,200, by=5)),
#     heightClasses = c(0,25,50,75,100,125,150,175, 200,seq(250,6000, by=50)),    
    #For fire
#     sparkLambda = 0, #Average number of sparks/day (0 means no fires)
#     dailySparkingTime = 960.0, #in minutes (960 min = 60 min/h * 16h)
#     sparkingProbabilityVector = NULL, #Relative values for unequal ignition probability
#     brigadeDensity =  0.001, #in number of fire brigades/square km2
#     brigadeDispatchTime = 1.0, #in minutes
#     brigadeFireArrivalTime = 50.0, #in minutes
#     extinctionCapacity = 5000, #in kW/(m*min)
#     fireBrandTimeStep = 1/60, #in minutes
#     numBrandsPerBurningTime = 1, # dimensionless
    #For fire behaviour
   # liveFMCmode = "swb",
   # useModelForLive = FALSE
    # initialMoisture = c(35,35,35,35,35),
    # surfaceToVolumeRatios = c(5600,358,98,6200,8000), #dimensionless
    # heatContent = c(18622, 18622, 18622, 19500, 20000), #kJ/kg
    # deadFuelMoistureExtinction = 30 #Percent
    
    #For dispersal
    # maxDispersal = 500,
    #For management
    # managementUnits = NULL,
    # managementPlan = NULL
))
}