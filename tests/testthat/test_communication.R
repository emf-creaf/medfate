#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)


test_that("Control parameter Rcpp structures can be copied to C++ structures",{
  expect_type(medfate:::.testControlListToStructure(defaultControl()), "integer")
})
test_that("Model input Rcpp structures can be copied to C++ structures",{
  controlGranier <- defaultControl("Granier")
  controlSperry <- defaultControl("Sperry")
  controlSureau <- defaultControl("Sureau")
  expect_type(medfate:::.testModelInputToStructure(spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)), "double")
  expect_type(medfate:::.testModelInputToStructure(spwbInput(exampleforest, examplesoil, SpParamsMED, controlSperry)), "double")
  expect_type(medfate:::.testModelInputToStructure(spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)), "double")
  expect_type(medfate:::.testModelInputToStructure(growthInput(exampleforest, examplesoil, SpParamsMED, controlGranier)), "double")
  expect_type(medfate:::.testModelInputToStructure(growthInput(exampleforest, examplesoil, SpParamsMED, controlSperry)), "double")
  expect_type(medfate:::.testModelInputToStructure(growthInput(exampleforest, examplesoil, SpParamsMED, controlSureau)), "double")
})