library(medfate)

data(exampleforestMED)
data(SpParamsMED)
data(examplemeteo)
examplemeteo2 <- examplemeteo
examplemeteo2$dates <- as.Date(row.names(examplemeteo2))
row.names(examplemeteo2) <- NULL

control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE
control_sperry <- defaultControl("Sperry")
control_sperry$verbose <- FALSE
control_cochard <- defaultControl("Cochard")
control_cochard$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- soil(defaultSoilParams(4))

test_that("spwb can be run in example and empty forests",{
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_cochard), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(emptyforest(), examplesoil, SpParamsMED, control_granier), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(emptyforest(), examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(emptyforest(), examplesoil, SpParamsMED, control_cochard), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("pwb can be run in example and empty forests",{
  x1 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_granier)
  S1 <- spwb(x1, examplemeteo[1:10,], latitude = 41.82592, elevation = 100)
  expect_s3_class(pwb(x1, examplemeteo[1:10,], W = as.matrix(S1$Soil[,1:4]), 
                       latitude = 41.82592, elevation = 100), "pwb")
})
test_that("spwb can be run using dates as columns",{
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_cochard), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
})

test_that("spwb can be run using species codes",{
  f <- exampleforestMED
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  expect_s3_class(spwb(forest2spwbInput(f, examplesoil, SpParamsMED, control_granier), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(f, examplesoil, SpParamsMED, control_sperry), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  expect_s3_class(spwb(forest2spwbInput(f, examplesoil, SpParamsMED, control_cochard), 
                       examplemeteo2[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  
})
