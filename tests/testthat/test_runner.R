#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("MedfateRunner initializes correctly", {
  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  mfr <- new(MedfateRunner, x1)
  expect_s4_class(mfr, "Rcpp_MedfateRunner")
})
