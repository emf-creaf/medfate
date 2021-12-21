library(medfate)

#Load example daily meteorological data
data(examplemeteo)

data(exampleforestMED)
data(SpParamsMED)
examplesoil = soil(defaultSoilParams(2))
control = defaultControl("Granier")
control$verbose = FALSE
x = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
S1<-spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)

examplesoil2 = soil(defaultSoilParams(4))
control = defaultControl("Sperry")
control$verbose = FALSE
x2 = forest2spwbInput(exampleforestMED,examplesoil2, SpParamsMED, control)
d = 100:110
S2<-spwb(x2, examplemeteo[d,], latitude = 41.82592, elevation = 100)


test_that("Can produce all basic plots",{
  expect_s3_class(plot(S1, "PET_Precipitation"), "ggplot")
  expect_s3_class(plot(S1, "PET_NetRain"), "ggplot")
  expect_s3_class(plot(S1, "Snow"), "ggplot")
  expect_s3_class(plot(S1, "Export"), "ggplot")
  expect_s3_class(plot(S1, "Evapotranspiration"), "ggplot")
  expect_s3_class(plot(S1, "SoilPsi"), "ggplot")
  expect_s3_class(plot(S1, "SoilRWC"), "ggplot")
  expect_s3_class(plot(S1, "SoilTheta"), "ggplot")
  expect_s3_class(plot(S1, "SoilVol"), "ggplot")
  expect_s3_class(plot(S1, "PlantExtraction"), "ggplot")
  expect_s3_class(plot(S1, "LAI"), "ggplot")
  expect_s3_class(plot(S1, "PlantLAI"), "ggplot")
  expect_s3_class(plot(S1, "PlantStress"), "ggplot")
  expect_s3_class(plot(S1, "PlantPsi"), "ggplot")
  expect_s3_class(plot(S1, "PlantTranspiration"), "ggplot")
  expect_s3_class(plot(S1, "TranspirationPerLeaf"), "ggplot")
})

test_that("Can produce all advanced plots",{
  expect_s3_class(plot(S2, "HydraulicRedistribution"), "ggplot")
  expect_s3_class(plot(S2, "SoilPlantConductance"), "ggplot")
  expect_s3_class(plot(S2, "LeafPsiMin"), "ggplot")
  expect_s3_class(plot(S2, "StemPsi"), "ggplot")
  expect_s3_class(plot(S2, "RootPsi"), "ggplot")
})
