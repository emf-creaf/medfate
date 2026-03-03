library(medfate)

#' Transforms Rcentury output for soil temperature and soil moisture from month to day resolution
env_daily<-function(out, yearIni = 2000) {
  yearMax <- yearIni + max(out$time) - 1
  dateIni <- as.Date(paste0(yearIni, "-01-01"))
  dateFin <- as.Date(paste0(yearMax, "-12-31"))
  days <- seq(dateIni, dateFin, by="day")
  months <- cut(days, breaks="months")
  unique_months <- unique(months)
  # Trim out to the complete years
  out <- out[1:length(unique_months), , drop = FALSE]
  df <- data.frame(date = days, Step = 1:length(days), SoilMoisture = NA, SoilTemperature = NA,
                   AirTemperature = NA, AirRelativeHumidity = NA)
  rel_mois <- out$asmos.1./ max(out$asmos.1.)
  for(i in 1:length(unique_months)) {
    sel_i <- months == unique_months[i]
    df$SoilMoisture[sel_i] <- rel_mois[i]
    df$SoilTemperature[sel_i] <- out$stemp[i]
  }
  return(df)
}

prd_daily<-function(out, yearIni = 2000, species = "") {
  yearMax <- yearIni + max(out$time) - 1
  dateIni <- as.Date(paste0(yearIni, "-01-01"))
  dateFin <- as.Date(paste0(yearMax, "-12-31"))
  days <- seq(dateIni, dateFin, by="day")
  months <- cut(days, breaks="months")
  unique_months <- unique(months)
  # Trim out to the complete years
  out <- out[1:length(unique_months), , drop = FALSE]
  # fbrprd -> production of fine branches (gC/m2/yr)
  # crtprd -> production of coarse roots (gC/m2/yr)
  # rlwprd -> large wood component C production (gC/m2/yr)
  # rlvprd -> leaves component C production (gC/m2/yr)
  df <- data.frame(date = days, Step = 1:length(days), Species = species, Leaves = NA, Twigs = NA, SmallBranches = NA, LargeWood = NA, CoarseRoots = NA, FineRoots = NA)
  for(i in 1:length(unique_months)) {
    sel_i <- months == unique_months[i]
    df$Leaves[sel_i] <- out$rlvprd[i]/365
    df$Twigs[sel_i] <- 0
    df$SmallBranches[sel_i] <- out$fbrprd[i]/365
    df$LargeWood[sel_i] <- out$rlwprd[i]/365
    df$CoarseRoots[sel_i] <- out$crtprd[i]/365
    df$FineRoots[sel_i] <- out$frtprd[i]/365
  }
  return(df)
}

sp_params <- function(rcentury_tree) {
  # % lignin leaves tree$`1`$df$WDLIG(1)
  # % lignin fine roots tree$`1`$df$WDLIG(2)
  # % lignin branches tree$`1`$df$WDLIG(3)
  # % lignin large wood tree$`1`$df$WDLIG(4)
  # % lignin coarse roots tree$`1`$df$WDLIG(5)
  # average 3-5 for sapwood lignin
  
  # read lignin/N ratios from
  # tree$cerfor(3,1,1) C/
  wdlig1 <- rcentury_tree$df[rcentury_tree$df[,2]=="WDLIG(1)",1]
  wdlig2 <- rcentury_tree$df[rcentury_tree$df[,2]=="WDLIG(2)",1]
  wdlig3 <- rcentury_tree$df[rcentury_tree$df[,2]=="WDLIG(3)",1]
  wdlig4 <- rcentury_tree$df[rcentury_tree$df[,2]=="WDLIG(4)",1]
  wdlig5 <- rcentury_tree$df[rcentury_tree$df[,2]=="WDLIG(5)",1]
  cn_leaves <- rcentury_tree$df[rcentury_tree$df[,2]=="CERFOR(3,1,1)",1]
  cn_fineroots <- rcentury_tree$df[rcentury_tree$df[,2]=="CERFOR(3,2,1)",1]
  cn_branches <- rcentury_tree$df[rcentury_tree$df[,2]=="CERFOR(3,3,1)",1]
  cn_largewood <- rcentury_tree$df[rcentury_tree$df[,2]=="CERFOR(3,4,1)",1]
  cn_coarseroots <- rcentury_tree$df[rcentury_tree$df[,2]=="CERFOR(3,5,1)",1]
  cn_wood <- mean(c(cn_branches, cn_largewood, cn_coarseroots))
  leaf_c <- 0.5
  wood_c <- 0.5
  fineroot_c <- 0.5
  n_leaves <- (1/cn_leaves)*leaf_c
  n_fineroot <- (1/cn_fineroots)*fineroot_c
  n_wood <- (1/cn_wood)*wood_c
  df <- data.frame(Species  = rcentury_tree$label,
                   LeafLignin = wdlig1*100, # %
                   WoodLignin = mean(wdlig3, wdlig4, wdlig5)*100,# %
                   FineRootLignin = wdlig2*100,# %
                   Nleaf = n_leaves*1000, # mg / g
                   Nsapwood = n_wood*1000, # mg / g
                   Nfineroot = n_fineroot*1000) # mg / g
  return(df)
}

test_that("DAYCENT can be run", {
  skip_on_ci()
  skip_on_cran()
  
  load("~/OneDrive/mcaceres_work/model_development/medfate_evaluation/RcenturyBenchmark/century.RData")
  out_ForestC <- readRDS("~/OneDrive/mcaceres_work/model_development/medfate_evaluation/RcenturyBenchmark/out_ForestC.rds")
  
  environmentalConditions <- env_daily(out_WaterTemp)
  litterProduction <- prd_daily(out_ForestC)
  
  paramsLitterDecomposition <- sp_params(tree$`1`)
  litterProduction$Species <- paramsLitterDecomposition$Species[1]
  
  nlitter <- 1
  litterData <- data.frame(Species = as.character(rep(paramsLitterDecomposition$Species[1], nlitter)),
                           Leaves = as.numeric(rep(0, nlitter)),
                           Twigs = as.numeric(rep(0, nlitter)),
                           SmallBranches = as.numeric(rep(0, nlitter)),
                           LargeWood = as.numeric(rep(0, nlitter)),
                           CoarseRoots = as.numeric(rep(0, nlitter)),
                           FineRoots = as.numeric(rep(0, nlitter)))
  
  nsnag <- 0
  snagData<- data.frame(Species = as.character(rep(NA, nsnag)),
                        DBH = as.numeric(rep(NA, nsnag)),
                        Height = as.numeric(rep(NA, nsnag)),
                        DeathAge = as.numeric(rep(NA, nsnag)),
                        SmallBranches = as.numeric(rep(NA, nsnag)),
                        LargeWood = as.numeric(rep(NA, nsnag)))
  
  SOCData <- c(SurfaceMetabolic = 0, SoilMetabolic = 0, SurfaceActive = 0, SoilActive = 0,
               SurfaceSlow = 0, SoilSlow = 0, SoilPassive = 0)
  
  control <- defaultControl()
  baseAnnualRates <- control$decompositionAnnualBaseRates
  annualTurnoverRate <- control$decompositionAnnualTurnoverRate
  
  l <- decomposition_DAYCENT(snagData, litterData, SOCData,
                        paramsLitterDecomposition,
                        baseAnnualRates,
                        annualTurnoverRate,
                        environmentalConditions,
                        litterProduction,
                        sand = 30,
                        clay = 20,
                        soilPH = 7, balanceCheck = FALSE)  
  expect_type(l, "list")
})
