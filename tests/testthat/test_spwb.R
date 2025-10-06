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

test_that("spwb can be run in example and empty forests",{
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
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
})

test_that("pwb can be run in example and empty forests",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  S1 <- spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100)
  expect_s3_class(pwb(x1, examplemeteo[1:10,], W = as.matrix(S1$Soil$SWC[,1:4]), 
                       latitude = 41.82592, elevation = 100), "pwb")
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
