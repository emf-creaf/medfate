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
})

test_that("growth_day gives same result with inner and direct calls",{
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  ic <- medfate:::.instanceCommunicationStructures(x1)
  s_inner <- medfate:::.growth_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.growth_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  
  x2 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  ic <- medfate:::.instanceCommunicationStructures(x2)
  s_inner <- medfate:::.growth_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::growth_day(x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.growth_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  
  x3 <- growthInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  ic <- medfate:::.instanceCommunicationStructures(x3)
  s_inner <- medfate:::.growth_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::growth_day(x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.growth_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  
})

