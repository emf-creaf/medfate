# Default parameters
df_soil <- defaultSoilParams()

# Initializes soil
s <- soil(df_soil)

#Initialize input
control <- defaultControl("Granier")
x1 <- spwbInput(exampleforest,s, SpParamsMED, control)

test_that("rainfall intensity routine", {
  control <- defaultControl()
  expect_type(hydrology_rainfallIntensity(10, 5, control$defaultRainfallIntensityPerMonth), "double")
})
test_that("infiltration routines", {
  expect_type(hydrology_infiltrationAmount(20, 4, s, "VG", "GreenAmpt1911"), "double")
  expect_type(hydrology_infiltrationAmount(20, 4, s, "VG", "Boughton1989"), "double")
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

test_that("water inputs routine", {
  x1$snowpack <- 5
  snowpack_ini <- x1$snowpack
  expect_type(hydrology_waterInputs(x1, 10, 4, 10, 15, 30, 100, 2, 100, 100, FALSE), "double")
  expect_equal(x1$snowpack, snowpack_ini) # Check that snowpack has not changed
  expect_type(hydrology_waterInputs(x1, 10, 4, 10, 15, 30, 100, 2, 100, 100, TRUE), "double")
  expect_false(x1$snowpack == snowpack_ini) # Check that snowpack has changed
})

test_that("soil water balance routine (buckets)", {
  s <- soil(df_soil)
  W_ini <- rlang::duplicate(s$W)
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "buckets", modifySoil = FALSE), "double")
  expect_equal(s$W, W_ini) # Check that soil water content has not changed
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "buckets", modifySoil = TRUE), "double")
  expect_false(all(s$W == W_ini)) # Check that soil water content has changed
})

test_that("soil water balance routine (single)", {
  s <- soil(df_soil)
  W_ini <- rlang::duplicate(s$W)
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "single", modifySoil = FALSE), "double")
  expect_equal(s$W, W_ini) # Check that soil water content has not changed
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "single", modifySoil = TRUE), "double")
  expect_false(all(s$W == W_ini)) # Check that soil water content has changed
})

test_that("soil water balance routine (dual)", {
  s <- soil(df_soil)
  W_ini <- rlang::duplicate(s$W)
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "dual", modifySoil = FALSE), "double")
  expect_equal(s$W, W_ini) # Check that soil water content has not changed
  expect_type(hydrology_soilWaterBalance(s, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
                                         soilDomains = "dual", modifySoil = TRUE), "double")
  expect_false(all(s$W == W_ini)) # Check that soil water content has changed
})