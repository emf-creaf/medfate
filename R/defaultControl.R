defaultControl<-function() {
  return(list(
    #For all
    verbose = TRUE,
    # swb
    soilFunctions = "SX",
    snowpack = TRUE,
    drainage = TRUE,
    transpirationMode = "Simple", #vs. "Complex"
    hydraulicCostFunction = 1, #vs. 2
    verticalLayerSize = 100,
    nStemSegments = 1,
    capacitance = FALSE,
    cavitationRefill = TRUE,
    ksymver = 0, # whole-stem symplastic vertical conductance
    klat = 0, # symplastic-apoplastic lateral conductance
    taper = TRUE,
    numericParams=list(maxNsteps = 400, psiStep = -0.001, psiMax = -10.0, ntrial = 200, psiTol = 0.0001, ETol = 0.0001),
    averageFracRhizosphereResistance = 0.15,
    Catm = 386,
    ndailysteps = 24,
    thermalCapacityLAI = 1000000,
    defaultWindSpeed = 5,
    # growth
    storagePool = "none"
    
#     #For water balance
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