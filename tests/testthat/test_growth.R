library(medfate)

data(exampleforest)
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

test_that("growth can be run in example and empty forests",{
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(emptyforest(), examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(emptyforest(), examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(emptyforest(), examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
})

test_that("growth can be run using dates as columns",{
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
})

test_that("growth can be run using species codes",{
  f <- exampleforest
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  expect_s3_class(growth(growthInput(f, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(f, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  expect_s3_class(growth(growthInput(f, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "growth")
  
})

test_that("growth can be run disabling output",{
  control <- control_granier
  control$soilResults <- FALSE
  control$snowResults <- FALSE
  control$standResults <- FALSE
  control$temperatureResults <- FALSE
  control$leafResults <- FALSE
  control$plantResults <- FALSE
  control$plantLabileCarbonBalanceResults  <- FALSE
  control$plantStructureResults  <- FALSE
  control$growthMortalityResults  <- FALSE
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control), 
                         examplemeteo[1:10,],
                         latitude = 41.82592, elevation = 100), "growth")
  
  expect_s3_class(medfate::growth_day(growthInput(exampleforest, examplesoil, SpParamsMED, control), 
                                      date, meteovec, latitude = 41.82592, elevation = 100, 
                                      slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  
  control <- control_sperry
  control$soilResults <- FALSE
  control$snowResults <- FALSE
  control$standResults <- FALSE
  control$temperatureResults <- FALSE
  control$leafResults <- FALSE
  control$plantResults <- FALSE
  control$plantLabileCarbonBalanceResults  <- FALSE
  control$plantStructureResults  <- FALSE
  control$growthMortalityResults  <- FALSE
  expect_s3_class(growth(growthInput(exampleforest, examplesoil, SpParamsMED, control), 
                         examplemeteo[1:10,],
                         latitude = 41.82592, elevation = 100), "growth")
  
  expect_s3_class(medfate::growth_day(growthInput(exampleforest, examplesoil, SpParamsMED, control), 
                                      date, meteovec, latitude = 41.82592, elevation = 100, 
                                      slope=0, aspect=0, modifyInput = FALSE), "growth_day")
  
})

test_that("growth_day gives same result with inner and direct calls",{
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  ic <- medfate:::instance_communication_structures(x1, "growth")
  s_dir <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  medfate::growth_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x1, "growth"), s_dir)

  x2 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  ic <- medfate:::instance_communication_structures(x2, "growth")
  medfate::growth_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::growth_day(x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x2, "growth"), s_dir)
  # Second call (after modifying ic)
  medfate::growth_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x2, "growth"), s_dir)
  
  x3 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  ic <- medfate:::instance_communication_structures(x3, "growth")
  s_dir <- medfate::growth_day(x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  medfate::growth_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x3, "growth"), s_dir)
  # Second call (after modifying ic)
  medfate::growth_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(medfate::copy_model_output(ic, x3, "growth"), s_dir)
  
})

test_that("growth_day gives same result with inner and direct calls with general communication",{
  x1i <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  numCohorts <- nrow(x1i$cohorts)
  nlayers <- nrow(examplesoil)
  ncanlayers <- nrow(x1i$canopy)
  ntimesteps <- control_granier$ndailysteps
  x1d <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_equal(x1i, x1d)
  ic <- medfate::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "growth")
  s_dir <- medfate::growth_day(x1d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  medfate::growth_day_inner(ic, x1i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x1i, x1d)
  expect_equal(medfate::copy_model_output(ic, x1i, "growth"), s_dir)

  x2i <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  x2d <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_equal(x2i, x2d)
  ic <- medfate:::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "growth")
  s_dir <- medfate::growth_day(x2d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  medfate::growth_day_inner(ic, x2i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x2i, x2d)
  expect_equal(medfate::copy_model_output(ic, x2i, "growth"), s_dir)

  x3i <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  x3d <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_equal(x3i, x3d)
  ic <- medfate:::general_communication_structures(numCohorts, nlayers, ncanlayers, ntimesteps, "growth")
  s_dir <- medfate::growth_day(x3d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  medfate::growth_day_inner(ic, x3i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x3i, x3d)
  expect_equal(medfate::copy_model_output(ic, x3i, "growth"), s_dir)

})