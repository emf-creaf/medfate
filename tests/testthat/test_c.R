library(medfate)

data(exampleforest)

forest_herbs <- exampleforest
forest_herbs$herbData <- forest_herbs$shrubData

data(SpParamsMED)
data(examplemeteo)
examplemeteo2 <- examplemeteo
row.names(examplemeteo2) <- as.character(examplemeteo2$dates)
examplemeteo2$dates <- NULL
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

control_agriculture <- defaultControl()
control_agriculture$verbose <- FALSE
control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE
control_sperry <- defaultControl("Sperry")
control_sperry$subdailyResults <- TRUE
control_sperry$verbose <- FALSE
control_sureau <- defaultControl("Sureau")
control_sureau$verbose <- FALSE
control_sureau$subdailyResults <- TRUE

#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("aspwb_day and aspwb_day_old return the same result", {
  xa <- aspwbInput(0.75, control_agriculture, examplesoil)
  sda <- aspwb_day(xa, date, meteovec, latitude = 41.82592, elevation = 100, modifyInput = FALSE) 
  sda_c <- medfate:::.aspwb_day_old(xa, date, meteovec, latitude = 41.82592, elevation = 100, modifyInput = FALSE) 
  expect_equal(sda, sda_c) # Check for same output
})

test_that("spwb_day and spwb_day_old return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate::spwb_day_old(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1, sd1_c) # Check for same output

  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1_c <- medfate::spwb_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1, sd1_c) # Check for same output
})

test_that("spwb_day and spwb_day_old return the same result with sperry",{
  control_sperry$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1_c <- medfate::spwb_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$spwbOutput, sd1_c$spwbOutput) # Check for same output
  expect_equal(x1$soil, x1_c$soil) # Check for same modified input
  expect_equal(x1$canopy, x1_c$canopy) # Check for same modified input
  expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output

  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1_c <- medfate::spwb_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output
})


test_that("spwb_day and spwb_day_old return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate::spwb_day_old(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output

  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1_c <- medfate::spwb_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output
})
test_that("aspwb and aspwb_old return the same result", {
  xa <- aspwbInput(0.75, control_agriculture, examplesoil)
  S1 <- aspwb(xa, examplemeteo[1:10,], latitude = 41.82592, elevation = 100) 
  S1_c <- medfate:::.aspwb_old(xa, examplemeteo[1:10,], latitude = 41.82592, elevation = 100) 
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
})

test_that("spwb and spwb_old return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$spwbOutput$internalWater, S1_c$spwbOutput$internalWater) # Check for same outputs
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
})

test_that("spwb and spwb_old return the same result with sperry",{
  control_sperry$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  # expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  # expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
})


test_that("spwb and spwb_old return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.spwb_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
})


test_that("growth_day and growth_day_old return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate:::.growth_day_old(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output

  control_granier$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1_c <- medfate:::.growth_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output
  expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
})

test_that("growth_day and growth_day_old return the same result with sperry",{
  control_sperry$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate:::.growth_day_old(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output

  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  sd1_c <- medfate:::.growth_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output
})

test_that("growth_day and growth_day_old return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate:::.growth_day_old(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output

  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1_c <- medfate:::.growth_day_old(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1$CarbonBalance, sd1_c$CarbonBalance) # Check for same output
  expect_equal(sd1$LabileCarbonBalance, sd1_c$LabileCarbonBalance) # Check for same output
  expect_equal(sd1$PlantStructure, sd1_c$PlantStructure) # Check for same output
  expect_equal(sd1$GrowthMortality, sd1_c$GrowthMortality) # Check for same output
})

test_that("growth and .growth_old return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
})

test_that("growth and .growth_old return the same result with sperry",{
  control_sperry$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
})

test_that("growth and .growth_old return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::growth(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate:::.growth_old(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$growthOutput$internalCarbon, S1_c$growthOutput$internalCarbon)
  expect_equal(S1$Plants, S1_c$Plants) # Check for same outputs
  expect_equal(S1$CarbonBalance, S1_c$CarbonBalance) # Check for same outputs
  # expect_equal(S1$LabileCarbonBalance, S1_c$LabileCarbonBalance)
  expect_equal(S1$PlantStructure, S1_c$PlantStructure)
  expect_equal(S1$GrowthMortality, S1_c$GrowthMortality)
  expect_equal(S1$DecompositionPools, S1_c$DecompositionPools)
})
