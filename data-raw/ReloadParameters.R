## Reload MEGAN parameters
# MEGANParams<-read.csv("data-raw/MEGANParams.csv",skip=0)
# usethis::use_data(MEGANParams, overwrite = T)
# rm(MEGANParams)

# Poblet tree data
poblet_trees <- openxlsx::read.xlsx("data-raw/PobletData.xlsx", sheet="TreeData")
usethis::use_data(poblet_trees, overwrite = T)


# SpParamsDefinition ----------------------------------------------------------------
SpParamsDefinition <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsDefinition.xlsx",
                                              sheet="Definition", na = "NA"), stringsAsFactors=FALSE)
SpParamsDefinition$Definition = stringi::stri_enc_toascii(SpParamsDefinition$Definition)
SpParamsDefinition$Units = stringi::stri_enc_toascii(SpParamsDefinition$Units)
SpParamsDefinition$Strict = as.logical(SpParamsDefinition$Strict)
usethis::use_data(SpParamsDefinition, overwrite = T)



# SpParamsMED [MODE A] -------------------------------------------------------------
SpParamsMED <-as.data.frame(readxl::read_xlsx("data-raw/InitialSpParamsMED.xlsx",
                                              sheet="InitialSpParamsMED", na = "NA"), stringsAsFactors=FALSE)
MFWdir = "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/"
harmonized_allometry_path <- "~/OneDrive/EMF_datasets/AllometryDatabases/Products/harmonized"
harmonized_trait_path <- "~/OneDrive/EMF_datasets/PlantTraitDatabases/Products/harmonized"

# Revised hydraulic/photosynthesis parameters
source(paste0(MFWdir, "Metamodelling_TR_WUE/R/utils.R"))
customParamsSpecies <- get_custom_params(paste0(MFWdir,"Metamodelling_TR_WUE/data-raw"))
SpParamsMED <- medfate::modifySpParams(SpParamsMED, customParamsSpecies, subsetSpecies = FALSE)
# Results of meta-modelling exercise
metamodellingParamsSpecies <- readRDS(paste0(MFWdir,"Metamodelling_TR_WUE/data/SpParamsMED/metamodelling_params.rds"))
metamodellingParamsSpecies$SpIndex <- NULL
SpParamsMED <- medfate::modifySpParams(SpParamsMED, metamodellingParamsSpecies, subsetSpecies = FALSE)
# Load growth calibration results
RGRcambiummaxTrees <- readRDS(paste0(MFWdir,"GrowthCalibration/data/SpParamsMED/output/RGRcambiummax_trees.rds")) |>
  dplyr::rename(Species = Name)
RGRcambiummaxTrees$SpIndex <- NULL
SpParamsMED <- medfate::modifySpParams(SpParamsMED, RGRcambiummaxTrees, subsetSpecies = FALSE)
# Load ingrowth calibration results
## SHOULD BE RECALIBRATED: THEY REFER TO INGROWTH (~7.5 cm) 
recruitmentParamsSpecies <- readRDS(paste0(MFWdir,"MortalityRegenerationCalibration/data/final_recruitment_params.rds"))
recruitmentParamsSpecies$SpIndex <- NULL
recruitmentParamsSpecies$IFNcodes <- NULL
recruitmentParamsSpecies$RecrTreeHeight <- recruitmentParamsSpecies$RecrTreeHeight/10
recruitmentParamsSpecies$IngrowthTreeDensity <- recruitmentParamsSpecies$RecrTreeDensity
recruitmentParamsSpecies$RecrTreeDensity <- NULL
SpParamsMED <- medfate::modifySpParams(SpParamsMED, recruitmentParamsSpecies, subsetSpecies = FALSE)
# Load Baseline mortality calibration results
mortalityParamsSpecies <- readRDS(paste0(MFWdir,"MortalityRegenerationCalibration/data/mort_rates.rds"))
mortalityParamsSpecies$SpIndex <- NULL
SpParamsMED <- medfate::modifySpParams(SpParamsMED, mortalityParamsSpecies, subsetSpecies = FALSE)
# Load SurvivalModel calibration results
survivalParamsSpecies <- readRDS(paste0(MFWdir,"MortalityRegenerationCalibration/data/survival_models.rds"))
SpParamsMED <- medfate::modifySpParams(SpParamsMED, survivalParamsSpecies, subsetSpecies = FALSE)
# Load SurvivalModel calibration results
resproutingParamsSpecies <- readxl::read_xlsx(paste0(MFWdir,"MortalityRegenerationCalibration/data-raw/ResproutingMED.xlsx"))
SpParamsMED <- medfate::modifySpParams(SpParamsMED, resproutingParamsSpecies, subsetSpecies = FALSE)
# Add bark thickness parameters from harmonized allometries
SpParamsMED <- traits4models::fill_medfate_allometries(SpParamsMED, 
                                                       harmonized_allometry_path = harmonized_allometry_path,
                                                       responses = "BarkThickness")
# Add maturation diameter parameters from harmonized traits
SpParamsMED <- traits4models::fill_medfate_traits(SpParamsMED, harmonized_trait_path = harmonized_trait_path,
                                                  parameters = "SeedProductionDiameter")
  
# Manual tuning
tree_all_cols = 31:43
names(SpParamsMED)[tree_all_cols] # CHECK!
#Use allometries of A. alba for P. abies
SpParamsMED[147,tree_all_cols] = SpParamsMED[1,tree_all_cols]
#Use allometries of A. alba for A. pinsapo
SpParamsMED[2,tree_all_cols] = SpParamsMED[1,tree_all_cols]
#Use allometries of A. alba for P. menziesii
SpParamsMED[163,tree_all_cols] = SpParamsMED[1,tree_all_cols]
#Use allometries of A. alba for Cedrus
SpParamsMED[34,tree_all_cols] = SpParamsMED[1,tree_all_cols]
#Use allometries of A. alba for T.baccata
SpParamsMED[201,tree_all_cols] = SpParamsMED[1,tree_all_cols]
#Use allometries of J. communis for J. phoenicea
SpParamsMED[118,tree_all_cols] = SpParamsMED[116,tree_all_cols]
#Use allometries of J. communis for J. phoenicea
SpParamsMED[121,tree_all_cols] = SpParamsMED[116,tree_all_cols]
#Use allometries of J. communis for Cupressus
SpParamsMED[62,tree_all_cols] = SpParamsMED[116,tree_all_cols]
pines = c("Pinus halepensis", "Pinus nigra", "Pinus pinea","Pinus sylvestris", "Pinus uncinata", "Pinus radiata", "Pinus pinaster")
SpParamsMED$fHDmin[SpParamsMED$Name %in% pines] = 80
SpParamsMED$fHDmax[SpParamsMED$Name %in% pines] = 160
#Save data
usethis::use_data(SpParamsMED, overwrite = T)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "SpParamsMED")
openxlsx::writeDataTable(wb, "SpParamsMED", SpParamsMED)
openxlsx::saveWorkbook(wb,"data-raw/SpParamsMED.xlsx", overwrite=TRUE)

rm(SpParamsMED)


# SpParamsMED [MODE B] -------------------------------------------------------------
# SpParamsMED <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsMED.xlsx",
#                                               sheet="SpParamsMED", na = "NA"), stringsAsFactors=FALSE)
# usethis::use_data(SpParamsMED, overwrite = T)


## Trait family means
trait_family_means = read.csv2("data-raw/trait_family_means.csv", dec=".")
usethis::use_data(trait_family_means, internal=TRUE, overwrite=TRUE)

## Modify exampleforest (after rebuilding)
data(exampleforest)
data("SpParamsMED")
exampleforest$treeData$Species[1] = "Pinus halepensis"
exampleforest$treeData$Species[2] = "Quercus ilex"
exampleforest$shrubData$Species[1] = "Quercus coccifera"
usethis::use_data(exampleforest, overwrite = T)
##Rebuild!


# Builds a fake observed data set from simulation results -------------------------------
library(medfate)
data(examplemeteo)
data(exampleforest)
PH_cohName = paste0("T1_",SpParamsMED$SpIndex[SpParamsMED$Name=="Pinus halepensis"])
QI_cohName = paste0("T2_",SpParamsMED$SpIndex[SpParamsMED$Name=="Quercus ilex"])
QC_cohName = paste0("S1_",SpParamsMED$SpIndex[SpParamsMED$Name=="Quercus coccifera"])
data(SpParamsMED)
examplesoil = defaultSoilParams(4)
control = defaultControl("Granier")
x1 = growthInput(exampleforest,examplesoil, SpParamsMED, control)
DBH_ini_PH = x1$above[PH_cohName, "DBH"]
DBH_ini_QI = x1$above[QI_cohName, "DBH"]
S1<-growth(x1, examplemeteo, latitude = 41.82592, elevation = 100)
fmc<-S1$Plants$LFMC
DI_PH = S1$PlantStructure$DBH[,PH_cohName] - c(DBH_ini_PH, S1$PlantStructure$DBH[-nrow(examplemeteo),PH_cohName])
DI_QI = S1$PlantStructure$DBH[,QI_cohName] - c(DBH_ini_QI, S1$PlantStructure$DBH[-nrow(examplemeteo),QI_cohName])
dates <- row.names(S1$Soil$SWC)
exampleobs = data.frame(dates = as.Date(dates),
                        SWC = S1$Soil$SWC[,1], 
                        ETR  = S1$WaterBalance$Evapotranspiration, 
                        E_PH = S1$Plants$Transpiration[,PH_cohName]/x1$above[PH_cohName,"LAI_expanded"],
                        E_QI = S1$Plants$Transpiration[,QI_cohName]/x1$above[QI_cohName,"LAI_expanded"],
                        FMC_PH= fmc[,PH_cohName],
                        FMC_QI = fmc[,QI_cohName],
                        BAI_PH = S1$GrowthMortality$SAgrowth[,PH_cohName],
                        BAI_QI = S1$GrowthMortality$SAgrowth[,QI_cohName],
                        DI_PH = DI_PH,
                        DI_QI = DI_QI)
#Add normal error
exampleobs$SWC = exampleobs$SWC + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$SWC)/4)
exampleobs$ETR = exampleobs$ETR + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$ETR)/4)
exampleobs$E_PH = exampleobs$E_PH + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_PH)/4)
exampleobs$E_QI = exampleobs$E_QI + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$E_QI)/4)
exampleobs$FMC_PH = exampleobs$FMC_PH + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$FMC_PH)/4)
exampleobs$FMC_QI = exampleobs$FMC_QI + rnorm(nrow(exampleobs), mean = 0, sd = sd(exampleobs$FMC_QI)/4)
exampleobs$BAI_PH = exampleobs$BAI_PH*exp(rnorm(nrow(exampleobs), mean = 0, sd = 0.5))
exampleobs$BAI_QI = exampleobs$BAI_QI*exp(rnorm(nrow(exampleobs), mean = 0, sd = 0.5))
exampleobs$DI_PH = exampleobs$DI_PH*exp(rnorm(nrow(exampleobs), mean = 0, sd = 0.5))
exampleobs$DI_QI = exampleobs$DI_QI*exp(rnorm(nrow(exampleobs), mean = 0, sd = 0.5))
names(exampleobs)[4:11] = c(paste0("E_",PH_cohName), paste0("E_",QI_cohName),
                            paste0("FMC_",PH_cohName),paste0("FMC_",QI_cohName),
                            paste0("BAI_",PH_cohName),paste0("BAI_",QI_cohName),
                            paste0("DI_",PH_cohName),paste0("DI_",QI_cohName))
row.names(exampleobs) <- NULL
usethis::use_data(exampleobs, overwrite = T)

#Rebuild!!!


# Check missing
apply(SpParamsMED,2, function(x) round(100*sum(is.na(x))/length(x),1))
apply(SpParamsMED[SpParamsMED$GrowthForm %in% c("Tree","Tree/Shrub"), ],2, function(x) round(100*sum(is.na(x))/length(x),1))
apply(SpParamsMED[SpParamsMED$GrowthForm %in% c("Shrub","Tree/Shrub"), ],2, function(x) round(100*sum(is.na(x))/length(x),1))
