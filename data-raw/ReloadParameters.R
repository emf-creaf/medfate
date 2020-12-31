## Script to reload parameters from 'SpParams.xlsx'
SpParamsMED <-as.data.frame(readxl::read_xlsx("data-raw/SpParams.xlsx",
                                              sheet="SpParamsMED", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsMED, overwrite = T)
rm(SpParamsMED)
SpParamsUS <-as.data.frame(readxl::read_xlsx("data-raw/SpParams.xlsx",
                                              sheet="SpParamsUS", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsUS, overwrite = T)
rm(SpParamsUS)

## Builds a fake observed data set from simulation results
library(medfate)
data(examplemeteo)
data(exampleforestMED)
data(SpParamsMED)
examplesoil1 = soil(defaultSoilParams(4))
control = defaultControl()
x1 = forest2spwbInput(exampleforestMED,examplesoil1, SpParamsMED, control)
S1<-spwb(x1, examplesoil1, examplemeteo, latitude = 41.82592, elevation = 100)

exampleobs = data.frame(SWC = S1$Soil$W.1*(soil_thetaFC(examplesoil1)[1]), 
             ETR  = S1$WaterBalance$Evapotranspiration, 
             E_T1_54 = S1$Plants$Transpiration[,"T1_54"]/x1$above["T1_54","LAI_expanded"],
             E_T2_68 = S1$Plants$Transpiration[,"T2_68"]/x1$above["T2_68","LAI_expanded"])
#Add normal error
exampleobs$SWC = exampleobs$SWC + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$SWC)/4)
exampleobs$ETR = exampleobs$ETR + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$ETR)/4)
exampleobs$E_T1_54 = exampleobs$E_T1_54 + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_T1_54)/4)
exampleobs$E_T2_68 = exampleobs$E_T2_68 + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_T2_68)/4)
usethis::use_data(exampleobs, overwrite = T)

#Rebuild!!!