#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

test_that("FCCS properties can be calculated", {
  expect_s3_class(fuel_FCCS(exampleforest, SpParamsMED), "data.frame")
})
  
test_that("FCCS behaviour can be calculated", {
  fccs_props <- fuel_FCCS(exampleforest, SpParamsMED)
  expect_type(fire_FCCS(fccs_props), "list")
})
