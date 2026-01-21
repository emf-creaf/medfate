# Default parameters
df_soil <- defaultSoilParams()

# Initializes soil
s <- soil(df_soil)

test_that("rainfall intensity routine", {
  control <- defaultControl()
  expect_type(hydrology_rainfallIntensity(10, 5, control$defaultRainfallIntensityPerMonth), "double")
})
test_that("infiltration routines", {
  expect_type(hydrology_infiltrationRepartition(100, s$widths, s$macro), "double")
  expect_type(hydrology_infiltrationBoughton(100, 30), "double")
})

test_that("soil evaporation routine", {
  expect_type(hydrology_soilEvaporationAmount(5, 10, 2), "double")
})
test_that("snow melting routine", {
  expect_type(hydrology_snowMelt(5, 2, 100,100), "double")
  expect_error(hydrology_snowMelt(5, NA, 100,100))
  expect_error(hydrology_snowMelt(5, 2, 100,NA))
})