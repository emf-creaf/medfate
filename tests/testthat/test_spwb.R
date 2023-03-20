library(medfate)

data(exampleforestMED)
data(SpParamsMED)
data(examplemeteo)


#Initialize soil with default soil params (4 layers)
examplesoil <- soil(defaultSoilParams(4))

test_that("spwb can be run in example and empty forests",{
  control <- defaultControl("Granier")
  control$verbose <- FALSE
  x1 = forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)
  expect_s3_class(spwb(x1, examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  control <- defaultControl("Sperry")
  control$verbose <- FALSE
  x2 = forest2spwbInput(emptyforest(), examplesoil, SpParamsMED, control)
  expect_s3_class(spwb(x2, examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "spwb")
  # control <- defaultControl("Cochard")
  # control$verbose <- FALSE
  # x2 = forest2spwbInput(emptyforest(),examplesoil, SpParamsMED, control)
  # expect_s3_class(spwb(x2, examplemeteo[1:10,],
  #                      latitude = 41.82592, elevation = 100), "spwb")
})

test_that("fordyn can be run using species codes",{
  control <- defaultControl("Granier")
  control$verbose <- FALSE
  f <- exampleforestMED
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  x1 = forest2spwbInput(f,examplesoil, SpParamsMED, control)
  expect_s3_class(spwb(x1, examplemeteo,
                       latitude = 41.82592, elevation = 100), "spwb")
})
