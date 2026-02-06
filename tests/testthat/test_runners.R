#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("SPWBrunner initializes correctly", {
  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  mfr <- new(SPWBrunner, x1, 41.82592, 100, 0, 0)
  expect_s4_class(mfr, "Rcpp_SPWBrunner")
  d <- 100
  meteovec <- unlist(examplemeteo[d,-1])
  date <- as.character(examplemeteo$dates[d])
  mfr$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfr$get_output(), "spwb_day")
})
