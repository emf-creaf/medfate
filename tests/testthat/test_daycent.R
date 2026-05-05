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

prd_daily_new<-function(out, rcentury_tree, yearIni = 2000) {
  decid <- rcentury_tree$df[rcentury_tree$df[,2]=="DECID",1]
  species <- rcentury_tree$label
  yearMax <- yearIni + max(out$time) - 1
  dateIni <- as.Date(paste0(yearIni, "-01-01"))
  dateFin <- as.Date(paste0(yearMax, "-12-31"))
  days <- seq(dateIni, dateFin, by="day")
  months <- cut(days, breaks="months")
  unique_months <- unique(months)
  leafdr <- c(rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(1)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(2)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(3)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(4)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(5)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(6)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(7)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(8)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(9)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(10)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(11)",1],
              rcentury_tree$df[rcentury_tree$df[,2]=="LEAFDR(12)",1])
  wooddr1 <- rcentury_tree$df[rcentury_tree$df[,2]=="WOODDR(1)",1] #fraction of forest which is deciduous; the fraction of leaves which fall during senescence month or at the end of the growing season, used when decid = 1 or 2
  wooddr2 <- rcentury_tree$df[rcentury_tree$df[,2]=="WOODDR(2)",1] #monthly death rate fraction for fine root component
  wooddr3 <- rcentury_tree$df[rcentury_tree$df[,2]=="WOODDR(3)",1] #monthly death rate fraction for fine branch component
  wooddr4 <- rcentury_tree$df[rcentury_tree$df[,2]=="WOODDR(4)",1] #monthly death rate fraction for large wood component
  wooddr5 <- rcentury_tree$df[rcentury_tree$df[,2]=="WOODDR(5)",1] #monthly death rate fraction for coarse root component
  # Trim out to the complete years
  out <- out[1:length(unique_months), , drop = FALSE]
  # fbrchc -> C in forest system fine branch component (g m-2).
  # crootc -> C in forest system coarse root component (g m-2).
  # rlwodc -> C in forest system large wood component (g m-2).
  # frootc -> C in forest system fine root component (g m-2).
  # rleavc -> C in forest system leaf component (g m-2).
  # print(wooddr1)
  df <- data.frame(date = days, Step = 1:length(days), Species = species, Leaves = NA, Twigs = NA, SmallBranches = NA, LargeWood = NA, CoarseRoots = NA, FineRoots = NA)
  for(i in 1:length(unique_months)) {
    sel_i <- months == unique_months[i]
    nmonth_days <- sum(sel_i)
    month_number <- as.numeric(format(days[sel_i], "%m")[1])
    if(decid == 0) {
      df$Leaves[sel_i] <- out$rleavc[i]*leafdr[month_number]/nmonth_days
    } else {
      if((i %% 12)==10) df$Leaves[sel_i] <- out$rleavc[i]*wooddr1/nmonth_days
      else df$Leaves[sel_i] <- 0
    }
    df$Twigs[sel_i] <- 0
    df$SmallBranches[sel_i] <- out$fbrchc[i]*wooddr3/nmonth_days
    df$LargeWood[sel_i] <- out$rlwodc[i]*wooddr4/nmonth_days
    df$CoarseRoots[sel_i] <- out$crootc[i]*wooddr5/nmonth_days
    df$FineRoots[sel_i] <- out$frootc[i]*wooddr2/nmonth_days
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
  
  load("/home/miquel/OneDrive/mcaceres_work/model_development/medfate_evaluation/RcenturyBenchmark/Rdata/RC2m_Process.RData")
  duke_tree <- tree$`3`
  
  environmentalConditions <- env_daily(out_duke)
  litterProduction <- prd_daily_new(out_duke, duke_tree)
  
  paramsLitterDecomposition <- sp_params(duke_tree)
  paramsAnatomy <- paramsLitterDecomposition[,c("Species"), drop = FALSE]
  paramsAnatomy$WoodDensity <- 0.5
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
  
  SOCData <- c(SurfaceMetabolic = 0, BelowgroundMetabolic = 0, SurfaceActive = 0, BelowgroundActive = 0,
               SurfaceSlow = 0, BelowgroundSlow = 0, BelowgroundPassive = 0)
  
  control <- defaultControl()
  baseAnnualRates <- control$decompositionAnnualBaseRates
  annualTurnoverRate <- control$decompositionAnnualTurnoverRate
  
  l <- decomposition_DAYCENT(snagData, litterData, SOCData,
                             paramsLitterDecomposition,
                             paramsAnatomy,
                             baseAnnualRates,
                             annualTurnoverRate,
                             environmentalConditions,
                             litterProduction,
                             sand = 30,
                             clay = 20,
                             soilPH = 7, balanceCheck = FALSE)  
  expect_type(l, "list")
})
