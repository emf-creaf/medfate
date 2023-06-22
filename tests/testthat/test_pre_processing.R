library(medfate)

data(exampleforestMED)
data(SpParamsMED)

control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE
control_sperry <- defaultControl("Sperry")
control_sperry$verbose <- FALSE
control_cochard <- defaultControl("Cochard")
control_cochard$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- soil(defaultSoilParams(4))

test_that("pressure-volume curves can be shown",{
  x1 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_granier)
  expect_s3_class(moisture_pressureVolumeCurvePlot(x1, "stem"), "ggplot")
  x2 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(moisture_pressureVolumeCurvePlot(x2, "stem"), "ggplot")
  expect_s3_class(moisture_pressureVolumeCurvePlot(x2, "leaf"), "ggplot")
  x3 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_cochard)
  expect_s3_class(moisture_pressureVolumeCurvePlot(x3, "stem"), "ggplot")
  expect_s3_class(moisture_pressureVolumeCurvePlot(x3, "leaf"), "ggplot")
})

test_that("hydraulic vulnerability curves can be shown",{
  x2 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x2, type = "leaf"), "ggplot")
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x2, type = "stem"), "ggplot")
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x2, type = "root"), "ggplot")
  expect_null(hydraulics_vulnerabilityCurvePlot(x2, soil = examplesoil, type = "rhizosphere"))
  
  x3 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_cochard)
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x3, type = "leaf"), "ggplot")
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x3, type = "stem"), "ggplot")
  expect_s3_class(hydraulics_vulnerabilityCurvePlot(x3, type = "root"), "ggplot")
  expect_null(hydraulics_vulnerabilityCurvePlot(x3, soil = examplesoil, type = "rhizosphere"))
})

test_that("hydraulic supply functions can be shown",{
  x2 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(hydraulics_supplyFunctionPlot(x2, type = "E"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x2, type = "ERhizo"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x2, type = "StemPsi"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x2, type = "dEdP"), "ggplot")
  
  x3 <- forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control_cochard)
  expect_s3_class(hydraulics_supplyFunctionPlot(x3, type = "E"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x3, type = "ERhizo"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x3, type = "StemPsi"), "ggplot")
  expect_s3_class(hydraulics_supplyFunctionPlot(x3, type = "dEdP"), "ggplot")
})
