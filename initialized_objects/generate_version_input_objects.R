library(medfate)

version <- as.character(packageVersion("medfate"))

#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

#Initialize control parameters
control <- defaultControl("Granier")
control$verbose <- FALSE

#Initialize input
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
saveRDS(x1, file = paste0("initialized_objects/spwbInput_", version,".rds"))

x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)
saveRDS(x1, file = paste0("initialized_objects/growthInput_", version,".rds"))
