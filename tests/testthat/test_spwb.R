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

test_that("spwb can be run using truncated root systems",{
  f <- exampleforest
  f$treeData$Z100 <- c(1500, 1900)
  f$shrubDataZ100 <- 800
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

test_that("spwb_day gives same result with inner and direct calls",{
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  ic <- medfate:::.instanceCommunicationStructures(x1)
  s_inner <- medfate:::.spwb_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.spwb_day_inner(ic, x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  ic <- medfate:::.instanceCommunicationStructures(x2)
  s_inner <- medfate:::.spwb_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::spwb_day(x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.spwb_day_inner(ic, x2, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)

  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  ic <- medfate:::.instanceCommunicationStructures(x3)
  s_inner <- medfate:::.spwb_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  s_dir <- medfate::spwb_day(x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  # Second call (after modifying ic)
  s_inner <- medfate:::.spwb_day_inner(ic, x3, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_equal(s_inner, s_dir)
  
})

test_that("spwb_day gives same result with inner and direct calls with general communication",{
  x1i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x1d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_equal(x1i, x1d)
  ic <- medfate:::.generalCommunicationStructures()
  s_inner <- medfate:::.spwb_day_inner(ic, x1i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x1d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x1i, x1d)
  expect_equal(s_inner, s_dir)

  x2i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  x2d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_equal(x2i, x2d)
  ic <- medfate:::.generalCommunicationStructures()
  s_inner <- medfate:::.spwb_day_inner(ic, x2i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x2d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x2i, x2d)
  expect_equal(s_inner, s_dir)

  x3i <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  x3d <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_equal(x3i, x3d)
  ic <- medfate:::.generalCommunicationStructures()
  s_inner <- medfate:::.spwb_day_inner(ic, x3i, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  s_dir <- medfate::spwb_day(x3d, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = TRUE)
  expect_equal(x3i, x3d)
  expect_equal(s_inner, s_dir)

})