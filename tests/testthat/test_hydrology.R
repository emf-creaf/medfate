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
  expect_true(sum(hydrology_infiltrationRepartition(100, s$widths, s$macro)) == 100)
  expect_type(hydrology_infiltrationBoughton(100, 30), "double")
})

test_that("soil evaporation routine", {
  W_ini_1 <- s$W[1]
  expect_type(hydrology_soilEvaporation(s, 0.0, "VG", 10, 50, FALSE), "double")
  expect_true(s$W[1] == W_ini_1) # Check that soil water content has not changed
  expect_type(hydrology_soilEvaporation(s, 0.0, "VG", 10, 50, TRUE), "double")
  expect_false(s$W[1] == W_ini_1) # Check that soil water content has changed
  expect_type(hydrology_soilEvaporationAmount(5, 10, 2), "double")
})

test_that("herbaceous transpiration routine", {
  W_ini <- rlang::duplicate(s$W)
  expect_type(hydrology_herbaceousTranspiration(10, 100.0, 1, s, "VG", FALSE), "double")
  expect_equal(s$W, W_ini) # Check that soil water content has not changed
  expect_type(hydrology_herbaceousTranspiration(10, 100.0, 1, s, "VG", TRUE), "double")
  expect_false(s$W[1] == W_ini[1]) # Check that soil water content has changed
})

test_that("snow melting routine", {
  expect_type(hydrology_snowMelt(5, 2, 100,100), "double")
  expect_error(hydrology_snowMelt(5, NA, 100,100))
  expect_error(hydrology_snowMelt(5, 2, 100,NA))
})