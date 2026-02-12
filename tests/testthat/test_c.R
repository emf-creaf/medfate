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

test_that("spwb_day can be run after reorganizing code",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
})

test_that("spwb_day and spwb_day_c return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(sd1, sd1_c) # Check for same output

  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  sd1_c <- medfate::spwb_day_c(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(sd1, sd1_c) # Check for same output
})

test_that("spwb_day and spwb_day_c return the same result with sperry",{
  control_sperry$rhizosphereOverlap <- "total"
  # x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  # sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  # sd1_c <- medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  # expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  # expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  # expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  # expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  # expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  # expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  # expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  # expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  # expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  # expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  # expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  # expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output
  # 
  # control_sperry$rhizosphereOverlap <- "partial"
  # x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  # sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  # x1_c <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  # sd1_c <- medfate::spwb_day_c(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  # expect_equal(sd1$WaterBalance, sd1_c$WaterBalance) # Check for same output
  # expect_equal(sd1$Plants, sd1_c$Plants) # Check for same output
  # expect_equal(sd1$Sunlit, sd1_c$Sunlit) # Check for same output
  # expect_equal(sd1$Shade, sd1_c$Shade) # Check for same output
  # expect_equal(sd1$PlantsInst, sd1_c$PlantsInst) # Check for same output
  # expect_equal(sd1$SunlitInst, sd1_c$SunlitInst) # Check for same output
  # expect_equal(sd1$ShadeInst, sd1_c$ShadeInst) # Check for same output
  # expect_equal(sd1$Soil, sd1_c$Soil) # Check for same output
  # expect_equal(sd1$EnergyBalance, sd1_c$EnergyBalance) # Check for same output
  # expect_equal(sd1$LWRExtinction, sd1_c$LWRExtinction) # Check for same output
  # expect_equal(sd1$RhizoPsi, sd1_c$RhizoPsi) # Check for same output
  # expect_equal(sd1$CanopyTurbulence, sd1_c$CanopyTurbulence) # Check for same output
})


test_that("spwb_day and spwb_day_c return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  sd1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  sd1_c <- medfate::spwb_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
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
  sd1_c <- medfate::spwb_day_c(x1_c, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
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
test_that("spwb and spwb_c return the same result with granier",{
  control_granier$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate::spwb_c(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- medfate::spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate::spwb_c(x1, examplemeteo, latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
})

test_that("spwb and spwb_c return the same result with sperry",{
  # control_sperry$rhizosphereOverlap <- "total"
  # x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  # S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  # S1_c <- medfate::spwb_c(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  # expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  # expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  # expect_equal(S1$Soil, S1_c$Soil)
  # expect_equal(S1$Plants, S1_c$Plants)
  # control_sperry$rhizosphereOverlap <- "partial"
  # x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  # S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  # S1_c <- medfate::spwb_c(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  # expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  # expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  # expect_equal(S1$Soil, S1_c$Soil)
  # expect_equal(S1$Plants, S1_c$Plants)
})


test_that("spwb and spwb_c return the same result with sureau",{
  control_sureau$rhizosphereOverlap <- "total"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate::spwb_c(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  S1 <- medfate::spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  S1_c <- medfate::spwb_c(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100, slope=0, aspect=0)
  expect_equal(S1$WaterBalance, S1_c$WaterBalance) # Check for same outputs
  expect_equal(S1$Stand, S1_c$Stand) # Check for same outputs
  expect_equal(S1$Soil, S1_c$Soil)
  expect_equal(S1$Plants, S1_c$Plants)
})


test_that("growth_day can be run after reorganizing code",{
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

  control_granier$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

  control_sperry$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  expect_equal(x1$soil$W, W_ini) # Check that W has not changed
  expect_s3_class(medfate::growth_day_c(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "growth_day")
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed

})
