## Reload MEGAN parameters
MEGANParams<-read.csv("data-raw/MEGANParams.csv",skip=0)
usethis::use_data(MEGANParams, overwrite = T)
rm(MEGANParams)

## Reload parameters from 'SpParams.xlsx'
SpParamsDefinition <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsDefinition.xlsx",
                                              sheet="Definition", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsDefinition, overwrite = T)
SpParamsMED <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsMED.xlsx",
                                              sheet="SpParamsMED", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsMED, overwrite = T)
rm(SpParamsMED)
SpParamsUS <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsUS.xlsx",
                                              sheet="SpParamsUS", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsUS, overwrite = T)
rm(SpParamsUS)

## Trait family means
trait_family_means = read.csv2("data-raw/trait_family_means.csv", dec=".")
usethis::use_data(trait_family_means, internal=TRUE, overwrite=TRUE)

## Modify exampleforestMED (after rebuilding)
data(exampleforestMED)
data("SpParamsMED")
exampleforestMED$treeData$Species[1] = SpParamsMED$SpIndex[SpParamsMED$Name=="Pinus halepensis"]
exampleforestMED$treeData$Species[2] = SpParamsMED$SpIndex[SpParamsMED$Name=="Quercus ilex"]
exampleforestMED$shrubData$Species[1] = SpParamsMED$SpIndex[SpParamsMED$Name=="Quercus coccifera"]
usethis::use_data(exampleforestMED, overwrite = T)
##Rebuild!

## Builds a fake observed data set from simulation results
library(medfate)
data(examplemeteo)
data(exampleforestMED)
PH_cohName = paste0("T1_",exampleforestMED$treeData$Species[1])
QI_cohName = paste0("T2_",exampleforestMED$treeData$Species[2])
QC_cohName = paste0("S1_",exampleforestMED$shrubData$Species[1])
data(SpParamsMED)
examplesoil1 = soil(defaultSoilParams(4))
control = defaultControl("Granier")
x1 = forest2growthInput(exampleforestMED,examplesoil1, SpParamsMED, control)
S1<-growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
fmc<-moisture_cohortFMC(S1, SpParamsMED)
exampleobs = data.frame(SWC = S1$Soil$W.1*(soil_thetaFC(examplesoil1)[1]), 
             ETR  = S1$WaterBalance$Evapotranspiration, 
             E_PH = S1$Plants$Transpiration[,PH_cohName]/x1$above[PH_cohName,"LAI_expanded"],
             E_QI = S1$Plants$Transpiration[,QI_cohName]/x1$above[QI_cohName,"LAI_expanded"],
             FMC_PH= fmc[,PH_cohName],
             FMC_QI = fmc[,QI_cohName],
             BAI_PH = S1$PlantStructure$SapwoodArea[,PH_cohName]*S1$GrowthMortality$SAgrowth[,PH_cohName],
             BAI_QI = S1$PlantStructure$SapwoodArea[,QI_cohName]*S1$GrowthMortality$SAgrowth[,QI_cohName])
#Add normal error
exampleobs$SWC = exampleobs$SWC + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$SWC)/4)
exampleobs$ETR = exampleobs$ETR + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$ETR)/4)
exampleobs$E_PH = exampleobs$E_PH + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_PH)/4)
exampleobs$E_QI = exampleobs$E_QI + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_QI)/4)
exampleobs$FMC_PH = exampleobs$FMC_PH + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$FMC_PH)/4)
exampleobs$FMC_QI = exampleobs$FMC_QI + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$FMC_QI)/4)
exampleobs$BAI_PH = exampleobs$BAI_PH + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$BAI_PH)/4)
exampleobs$BAI_QI = exampleobs$BAI_QI + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$BAI_QI)/4)

names(exampleobs)[3:8] = c(paste0("E_",PH_cohName), paste0("E_",QI_cohName),
                           paste0("FMC_",PH_cohName),paste0("FMC_",QI_cohName),
                           paste0("BAI_",PH_cohName),paste0("BAI_",QI_cohName))

usethis::use_data(exampleobs, overwrite = T)

#Rebuild!!!