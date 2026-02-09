#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("SPWB_runner initializes and runs correctly", {
  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  mfr <- new(SPWB_runner, x1, 41.82592, 100, 0, 0)
  expect_s4_class(mfr, "Rcpp_SPWB_runner")
  d <- 100
  meteovec <- unlist(examplemeteo[d,-1])
  date <- as.character(examplemeteo$dates[d])
  mfr$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfr$get_output(), "spwb_day")
})

test_that("SPWB_multiple_runner initializes and runs correctly", {
  controlSureau <- defaultControl("Sureau")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  n = 100;
  x_vec <- vector("list", n)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
  latitude_vec <- rep(41, n)
  elevation_vec <- rep(100, n)
  slope_vec <- rep(0, n)
  aspect_vec <- rep(0, n)
  mfmr <- new(SPWB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_SPWB_multiple_runner")
  d <- 100
  meteovec1 <- unlist(examplemeteo[d,-1])
  meteo_vec <- vector("list", n)
  for(i in 1:n) meteo_vec[[i]] = rlang::duplicate(meteovec1)
  date <- as.character(examplemeteo$dates[d])
  system.time(mfmr$run_day(date, meteo_vec, FALSE))
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
  Sys.setenv(ARMA_NUM_THREADS = "1")
  system.time(mfmr$run_day(date, meteo_vec, TRUE))
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
})
