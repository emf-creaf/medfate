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

test_that("spwb can be run in example, herbs and empty forests",{
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest2, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest2, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest2, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(emptyforest(), examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(emptyforest(), examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(emptyforest(), examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(forest_herbs, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(forest_herbs, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(forest_herbs, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("pwb can be run in example and empty forests",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100)
  expect_s3_class(pwb(x1, examplemeteo[1:10,], W = as.matrix(S1$Soil$SWC[,1:4]), 
                       latitude = 41.82592, elevation = 100), "pwb")
})

test_that("spwb can be run with missing values columns",{
  examplemeteo_missing <- examplemeteo[1:2, ]
  examplemeteo_missing$Radiation <- NA
  expect_s3_class(suppressWarnings(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                        examplemeteo_missing,
                                        latitude = 41.82592, elevation = 100)), "spwb")
  examplemeteo_missing$MinRelativeHumidity <- NA
  examplemeteo_missing$MaxRelativeHumidity <- NA
  expect_s3_class(suppressWarnings(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                        examplemeteo_missing,
                                        latitude = 41.82592, elevation = 100)), "spwb")
  examplemeteo_missing$MaxTemperature <- NA
  expect_error(suppressWarnings(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                        examplemeteo_missing,
                                        latitude = 41.82592, elevation = 100)))
})

test_that("spwb can be run using dates as columns",{
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run using species codes",{
  f <- exampleforest
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  expect_s3_class(spwb(spwbInput(f, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(f, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(f, examplesoil, SpParamsMED, control_sureau), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run using single soilDomains",{
  control_granier_single <- control_granier
  control_granier_single$soilDomains <- "single"
  control_sperry_single <- control_sperry
  control_sperry_single$soilDomains <- "single"
  control_sureau_single <- control_sureau
  control_sureau_single$soilDomains <- "single"
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier_single), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry_single), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau_single), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run using dual soilDomains",{
  control_granier_dual <- control_granier
  control_granier_dual$soilDomains <- "dual"
  control_sperry_dual <- control_sperry
  control_sperry_dual$soilDomains <- "dual"
  control_sureau_dual <- control_sureau
  control_sureau_dual$soilDomains <- "dual"
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier_dual), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry_dual), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau_dual), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})
test_that("spwb can be run using partial rhizosphere overlap",{
  control_granier_partial_overlap <- control_granier
  control_granier_partial_overlap$rhizosphereOverlap <- "partial"
  control_sperry_partial_overlap <- control_sperry
  control_sperry_partial_overlap$rhizosphereOverlap <- "partial"
  control_sureau_partial_overlap <- control_sureau
  control_sureau_partial_overlap$rhizosphereOverlap <- "partial"
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier_partial_overlap), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry_partial_overlap), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau_partial_overlap), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run with fire hazard results",{
  control_fh <- control_granier
  control_fh$fireHazardResults <- TRUE
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_fh), 
                       examplemeteo[1:2,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run with less output",{
  control_less <- control_granier
  control_less$standResults <- FALSE
  control_less$soilResults <- FALSE
  control_less$plantResults <- FALSE
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_less), 
                       examplemeteo[1:2,],
                       latitude = 41.82592, elevation = 100), "spwb")
})
