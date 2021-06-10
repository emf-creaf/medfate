defaultControl<-function(transpirationMode = "Granier") {
  return(list(
    #For all functions
    modifyInput = TRUE,
    verbose = TRUE,
    subdailyResults = FALSE,
    
    # For water balance
    transpirationMode = transpirationMode,
    soilFunctions = ifelse(transpirationMode=="Sperry", "VG", "SX"),
    defaultWindSpeed = 2.5, #m/s
    snowpack = TRUE,
    leafPhenology = TRUE,
    rockyLayerDrainage = TRUE,
    unlimitedSoilWater = FALSE,
    plantWaterPools = FALSE,
    unfoldingDD = 300,
    verticalLayerSize = 100,
    windMeasurementHeight = 200,
    cavitationRefill = "total",
    
    #spwb with sperry
    ndailysteps = 24,
    nsubsteps = 3600,
    cochard = FALSE,
    capacitance = FALSE,
    taper = TRUE,
    multiLayerBalance = FALSE,
    gainModifier = 1, 
    costModifier = 1, 
    costWater = "dEdP",
    klatstem = 0.01, # stem symplastic-apoplastic lateral conductance
    klatleaf = 0.01, # leaf symplastic-apoplastic lateral conductance
    numericParams=list(maxNsteps = 400, ntrial = 200, psiTol = 0.0001, ETol = 0.0000001),
    fracLeafResistance = NA,
    fracRootResistance = 0.4,
    averageFracRhizosphereResistance = 0.15,
    Catm = 386,
    thermalCapacityLAI = 1000000,
    boundaryLayerSize = 2000,
    refillMaximumRate = 0.05,
    
    # growth/mortality
    allowDessication = TRUE,
    allowStarvation = TRUE,
    allowDefoliation = TRUE,
    sinkLimitation = TRUE,
    shrubDynamics = FALSE,
    allocationStrategy = "Plant_kmax",
    nonStomatalPhotosynthesisLimitation = TRUE,
    phloemConductanceFactor = 0.2, # phloem conductance per leaf area basis (l*m-2*MPa-1*s-1)
    nonSugarConcentration = 0.25, # mol · l-1
    equilibriumOsmoticConcentration = list(leaf = 0.8, sapwood = 0.6),  # (Paljakka et al. 2017)
    minimumRelativeSugarForGrowth = 0.5,
    # Ogle and Pacala 2010, Tree Physiology 29, 587–605
    respirationRates = list(leaf = 0.00260274, 
                            sapwood = 6.849315e-05, 
                            fineroot = 0.002054795), # g gluc · g dw -1 · day -1
    turnoverRates = list(sapwood = 0.0001261398, # day-1 Equivalent to annual 4.5% 1-(1-0.045)^(1.0/365)
                      fineroot = 0.001897231), #day-1 Equivalent to annual 50% 1-(1-0.5)^(1.0/365)
    constructionCosts = list(leaf = 1.5, 
                             sapwood = 1.47, 
                             fineroot = 1.30), #  g gluc · g dw -1
    maximumRelativeGrowthRates = list(leaf = 0.01, # m2 leaf ·cm-2 sapwood· day-1
                                      sapwood = 0.002, # cm2 sapwood ·cm-2 sapwood· day-1
                                      fineroot = 0.1), # g dw · g dw -1 · day -1
    mortalityMode = "density/deterministic",
    mortalityBaselineRate = 0.01,
    mortalityRelativeSugarThreshold = 0.3,
    mortalityRWCThreshold = 0.3,
    
    #dynamics
    recruitmentMode = "deterministic",
    removeDeadCohorts=TRUE,
    minimumCohortDensity = 1,
    seedRain = NULL,
    seedProductionTreeHeight = 300,
    seedProductionShrubHeight = 30,
    minTempRecr	= 0,
    minMoistureRecr	= 0.3,
    minFPARRecr = 10,
    recrTreeDBH = 1,
    recrTreeDensity = 100,
    recrTreeHeight = 100,
    recrShrubCover = 1,
    recrShrubHeight = 10,
    recrTreeZ50 = 100,
    recrShrubZ50 = 50,
    recrTreeZ95 = 1000,
    recrShrubZ95 = 500
    
    
    #For forest dynamics
#     freqZopt = 20,
#     sap2tree=TRUE,
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