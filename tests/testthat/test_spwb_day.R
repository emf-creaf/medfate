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
control_sperry$verbose <- FALSE
control_sureau <- defaultControl("Sureau")
control_sureau$verbose <- FALSE


#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("spwb_day can be run with missing meteo data",{
  
  meteovec_mis <- meteovec
  meteovec_mis[4:7] <- NA
  
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["MinTemperature"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                    date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["MaxTemperature"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                 date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["Precipitation"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                 date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
})


test_that("spwb_day gives same result with inner and direct calls",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  ic <- medfate:::instance_communication_structures(x1, "spwb")
  s_dir <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  medfate::spwb_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x1, "spwb"), s_dir)
  # Second call (after modifying ic)
  medfate::spwb_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x1, "spwb"), s_dir)
  
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  ic <- medfate:::instance_communication_structures(x2, "spwb")
  s_dir <- medfate::spwb_day(x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  medfate::spwb_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x2, "spwb"), s_dir)
  # Second call (after modifying ic)
  medfate::spwb_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x2, "spwb"), s_dir)
  
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  ic <- medfate:::instance_communication_structures(x3, "spwb")
  s_dir <- medfate::spwb_day(x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  medfate::spwb_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x3, "spwb"), s_dir)
  # Second call (after modifying ic)
  medfate::spwb_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x3, "spwb"), s_dir)
})

test_that("spwb_day gives same result with inner and direct calls with general communication",{
  x1i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x1d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_equal(x1i, x1d)
  numCohorts <- nrow(x1i$cohorts)
  nlayers <- nrow(examplesoil)
  ncanlayers <- nrow(x1i$canopy)
  ntimesteps <- control_granier$ndailysteps
  ic <- medfate::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "spwb")
  medfate::spwb_day_inner(ic, x1i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x1d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x1i, x1d)
  expect_equal(medfate::copy_model_output(ic, x1i, "spwb"), s_dir)

  x2i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  x2d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_equal(x2i, x2d)
  ic <- medfate::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "spwb")
  medfate::spwb_day_inner(ic, x2i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x2d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x2i, x2d)
  expect_equal(medfate::copy_model_output(ic, x2i, "spwb"), s_dir)

  x3i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  x3d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_equal(x3i, x3d)
  ic <- medfate::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "spwb")
  medfate::spwb_day_inner(ic, x3i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x3d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x3i, x3d)
  expect_equal(medfate::copy_model_output(ic, x3i, "spwb"), s_dir)

})

test_that("spwb_day can be run with truncated root distribution",{
  exampleforest$shrubData$Z50 <- 100
  exampleforest$shrubData$Z95 <- 200
  control_granier$truncateRootDistribution<- TRUE
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_s3_class(spwb_day(x1, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sperry$truncateRootDistribution<- TRUE
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(spwb_day(x2, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sureau$truncateRootDistribution<- TRUE
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_s3_class(spwb_day(x3, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
})

test_that("spwb_day can be run with truncated root distribution and rhizosphere overlap",{
  exampleforest$shrubData$Z50 <- 100
  exampleforest$shrubData$Z95 <- 200
  control_granier$truncateRootDistribution<- TRUE
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_s3_class(spwb_day(x1, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sperry$truncateRootDistribution<- TRUE
  control_sperry$rhizosphereOverlap <- "partial"
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(spwb_day(x2, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sureau$truncateRootDistribution<- TRUE
  control_sureau$rhizosphereOverlap <- "partial"
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_s3_class(spwb_day(x3, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
})


test_that("spwb_day can be run after reorganizing code",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  Psi_ini <- rlang::duplicate(x1$internalWater$PlantPsi)
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$internalWater$PlantPsi, Psi_ini) # Check that psi has not changed
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$internalWater$PlantPsi == Psi_ini)) # Check that psi has changed
  
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  Psi_ini <- rlang::duplicate(x1$internalWater$PlantPsi)
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$internalWater$PlantPsi, Psi_ini) # Check that psi has not changed
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$internalWater$PlantPsi == Psi_ini)) # Check that psi has changed
  
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  Psi_ini <- rlang::duplicate(x1$internalWater$PlantPsi)
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$internalWater$PlantPsi, Psi_ini) # Check that psi has not changed
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$internalWater$PlantPsi == Psi_ini)) # Check that psi has changed
  
  control_sureau$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  Psi_ini <- rlang::duplicate(x1$internalWater$PlantPsi)
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE), "spwb_day")
  expect_equal(x1$internalWater$PlantPsi, Psi_ini) # Check that psi has not changed
  expect_s3_class(medfate::spwb_day_test(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE), "spwb_day")
  expect_false(all(x1$internalWater$PlantPsi == Psi_ini)) # Check that psi has changed
  
})
  