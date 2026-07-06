## Reload MEGAN parameters
# MEGANParams<-read.csv("data-raw/MEGANParams.csv",skip=0)
# usethis::use_data(MEGANParams, overwrite = T)
# rm(MEGANParams)

# Build example forest
exampleforest <- emptyforest(ntree = 2, nshrub = 1)
exampleforest$treeData$Species <- c("Pinus halepensis", "Quercus ilex")
exampleforest$treeData$N <- c(168, 384)
exampleforest$treeData$DBH <- c(37.55, 14.60)
exampleforest$treeData$Height <- c(800, 660)
exampleforest$treeData$Z50 <- c(100, 300)
exampleforest$treeData$Z95 <- c(300, 1000)
exampleforest$shrubData$Species <- c("Quercus coccifera")
exampleforest$shrubData$Cover <- 3.75
exampleforest$shrubData$Height <- 80
exampleforest$shrubData$Z50 <- 200
exampleforest$shrubData$Z95 <- 1000
usethis::use_data(exampleforest, overwrite = T)

exampleforest2 <- exampleforest
exampleforest2$treeData$N <- c(NA, NA)
exampleforest2$treeData$DBH <- c(NA, NA)
exampleforest2$treeData$LAI <- c(0.8, 0.5)
exampleforest2$treeData$CrownRatio <- c(0.66, 0.60)
exampleforest2$shrubData$Cover <- NA
exampleforest2$shrubData$LAI <- 0.03
exampleforest2$shrubData$CrownRatio <- 0.8
exampleforest$herbHeight <- 20
exampleforest2$herbLAI <- 0.25
usethis::use_data(exampleforest2, overwrite = T)



# Poblet tree data
poblet_trees <- openxlsx::read.xlsx("data-raw/PobletData.xlsx", sheet="TreeData")
usethis::use_data(poblet_trees, overwrite = T)


# SpParamsDefinition ----------------------------------------------------------------
SpParamsDefinition <-as.data.frame(readxl::read_xlsx("data-raw/SpParamsDefinition.xlsx",
                                              sheet="Definition", na = "NA"), stringsAsFactors=FALSE)
SpParamsDefinition$Definition = stringi::stri_enc_toascii(SpParamsDefinition$Definition)
SpParamsDefinition$Units = stringi::stri_enc_toascii(SpParamsDefinition$Units)
SpParamsDefinition$Strict = as.logical(SpParamsDefinition$Strict)
if(any(is.na(SpParamsDefinition$ParameterName))) stop("Missing values in ParameterName")
if(any(is.na(SpParamsDefinition$Definition))) stop("Missing values in Definition")
if(any(is.na(SpParamsDefinition$Type))) stop("Missing values in Type")
if(any(is.na(SpParamsDefinition$Strict))) stop("Missing values in Strict")
usethis::use_data(SpParamsDefinition, overwrite = T)



# SpParamsMED -------------------------------------------------------------
SpParamsTaxonomy <-as.data.frame(readxl::read_xlsx("data-raw/InitialTaxonomySpParamsMED.xlsx", na = "NA"), stringsAsFactors=FALSE)

SpParamsMED <- traits4models::init_medfate_params(SpParamsTaxonomy, complete_rows = FALSE, sort = FALSE)
SpParamsMED$IFNcodes <- SpParamsTaxonomy$IFNcodes
SpParamsMED <- SpParamsMED|>
  dplyr::relocate(IFNcodes, .before=AcceptedName)

MFWdir = "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/"
harmonized_allometry_path <- "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/traits_and_models/data/harmonized_allometry_sources/"
harmonized_trait_path <- "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/traits_and_models/data/harmonized_trait_sources/"

# Filling structural parameters from inventory data 
cli::cli_h2("SpParamsMED filling parameters from IFN")
sf_IFN3 <- readRDS("~/OneDrive/mcaceres_work/model_initialisation/medfate_initialisation/IFN2medfate/data/SpParamsES/IFN3/soilmod/IFN3_spain_soilmod_WGS84.rds")
SpParamsMED<- traits4models::fill_medfate_inventory_traits(SpParamsMED, sf_IFN3,
                                                           progress = TRUE)

cli::cli_h2("SpParamsMED filling parameters from harmonized allometries")
SpParamsMED <- traits4models::fill_medfate_allometries(SpParamsMED, 
                                                       harmonized_allometry_path = harmonized_allometry_path)
cli::cli_h2("SpParamsMED filling parameters from harmonized traits")
SpParamsMED <- traits4models::fill_medfate_traits(SpParamsMED, harmonized_trait_path = harmonized_trait_path)

# Revised hydraulic/photosynthesis parameters
cli::cli_h2("Calibration/metamodelling parameters")
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

# Complete strict (for taxa) 
cli::cli_h2("SpParamsES completing strict")
SpParamsMED <- traits4models::complete_medfate_strict(SpParamsMED)

# Complete strict for non-taxa or delete them 
cli::cli_h2("Cleaning and checking")
traits4models::check_medfate_params(SpParamsMED, check_consistency = TRUE)

#Save data
usethis::use_data(SpParamsMED, overwrite = T)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "SpParamsMED")
openxlsx::writeDataTable(wb, "SpParamsMED", SpParamsMED)
openxlsx::saveWorkbook(wb,"data-raw/SpParamsMED.xlsx", overwrite=TRUE)

rm(SpParamsMED)


# Trait family means ------------------------------------------------------
harmonized_trait_path <- "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/traits_and_models/data/harmonized_trait_sources"
trait_family_means <- traits4models::taxon_trait_summary(harmonized_trait_path, taxonomic_level = "family", 
                                                         traits =c("WoodDensity", "LeafDensity", "WoodC", 
                                                                    "Ptlp", "LeafPI0", "LeafEPS", "LeafAF",
                                                                    "Gswmin", "Gswmax",
                                                                    "Nleaf", "Nsapwood", "Nfineroot",
                                                                    "Ks", "VCstem_P50", 
                                                                    "Al2As", "conduit2sapwood")) |>
  dplyr::filter(!is.na(family)) |>
  dplyr::rename("Kmax_stemxylem" = "Ks")
row.names(trait_family_means) <- trait_family_means$family
trait_family_means <- trait_family_means |> dplyr::select(-family)
write.table(trait_family_means, "data-raw/trait_family_means.csv", sep=";")
usethis::use_data(trait_family_means, internal=TRUE, overwrite=TRUE)


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

