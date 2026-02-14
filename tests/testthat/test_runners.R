#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("WB_runner initializes and runs correctly", {
  d <- 100
  meteovec <- unlist(examplemeteo[d,-1])
  date <- as.character(examplemeteo$dates[d])

  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  mfra <- new(WB_runner, xa, 41.82592, 100, 0, 0)
  expect_s4_class(mfra, "Rcpp_WB_runner")
  mfra$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfra$get_output(), "aspwb_day")
  W_ini <- rlang::duplicate(xa$soil$W)
  mfra$update_input(xa)
  expect_false(all(xa$soil$W == W_ini)) # Check that W has changed

  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  mfr1 <- new(WB_runner, x1, 41.82592, 100, 0, 0)
  expect_s4_class(mfr1, "Rcpp_WB_runner")
  mfr1$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfr1$get_output(), "spwb_day")
  W_ini <- rlang::duplicate(x1$soil$W)
  mfr1$update_input(x1)
  expect_false(all(x1$soil$W == W_ini)) # Check that W has changed
  
  controlSperry <- defaultControl("Sperry")
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSperry)
  mfr2 <- new(WB_runner, x2, 41.82592, 100, 0, 0)
  expect_s4_class(mfr2, "Rcpp_WB_runner")
  mfr2$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfr2$get_output(), "spwb_day")
  W_ini <- rlang::duplicate(x2$soil$W)
  mfr2$update_input(x2)
  expect_false(all(x2$soil$W == W_ini)) # Check that W has changed
  
  controlSureau <- defaultControl("Sureau")
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  mfr3 <- new(WB_runner, x3, 41.82592, 100, 0, 0)
  expect_s4_class(mfr3, "Rcpp_WB_runner")
  mfr3$run_day(date, meteovec, 0, NULL, NA)
  expect_s3_class(mfr3$get_output(), "spwb_day")
  W_ini <- rlang::duplicate(x3$soil$W)
  mfr3$update_input(x3)
  expect_false(all(x3$soil$W == W_ini)) # Check that W has changed
})

test_that("WB_multiple_runner initializes and runs correctly", {
  n <- 100
  d <- 100
  meteovec1 <- unlist(examplemeteo[d,-1])
  meteo_vec <- vector("list", n)
  for(i in 1:n) meteo_vec[[i]] = rlang::duplicate(meteovec1)
  x_vec <- vector("list", n)
  latitude_vec <- rep(41, n)
  elevation_vec <- rep(100, n)
  slope_vec <- rep(0, n)
  aspect_vec <- rep(0, n)

  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(xa)
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  date <- as.character(examplemeteo$dates[d])
  system.time(mfmr$run_day(date, meteo_vec, FALSE))
  expect_s3_class(mfmr$get_output_at(1), "aspwb_day")
  expect_s3_class(mfmr$get_output_at(2), "aspwb_day")
  Sys.setenv(ARMA_NUM_THREADS = "1")
  system.time(mfmr$run_day(date, meteo_vec, TRUE))
  expect_s3_class(mfmr$get_output_at(1), "aspwb_day")
  expect_s3_class(mfmr$get_output_at(2), "aspwb_day")

  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  date <- as.character(examplemeteo$dates[d])
  system.time(mfmr$run_day(date, meteo_vec, FALSE))
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
  Sys.setenv(ARMA_NUM_THREADS = "1")
  system.time(mfmr$run_day(date, meteo_vec, TRUE))
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
})

test_that("GROWTH_runner initializes and runs correctly", {
  d <- 100
  meteovec <- unlist(examplemeteo[d,-1])
  date <- as.character(examplemeteo$dates[d])
  
  # controlGranier <- defaultControl("Granier")
  # x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  # mfr1 <- new(GROWTH_runner, x1, 41.82592, 100, 0, 0)
  # expect_s4_class(mfr1, "Rcpp_GROWTH_runner")
  # mfr1$run_day(date, meteovec, 0, NULL, NA)
  # expect_s3_class(mfr1$get_output(), "growth_day")
  # 
  # controlSperry <- defaultControl("Sperry")
  # x2 <- growthInput(exampleforest, examplesoil, SpParamsMED, controlSperry)
  # mfr2 <- new(GROWTH_runner, x2, 41.82592, 100, 0, 0)
  # expect_s4_class(mfr2, "Rcpp_GROWTH_runner")
  # mfr2$run_day(date, meteovec, 0, NULL, NA)
  # expect_s3_class(mfr2$get_output(), "growth_day")
  # 
  # controlSureau <- defaultControl("Sureau")
  # x3 <- growthInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  # mfr3 <- new(GROWTH_runner, x3, 41.82592, 100, 0, 0)
  # expect_s4_class(mfr3, "Rcpp_GROWTH_runner")
  # mfr3$run_day(date, meteovec, 0, NULL, NA)
  # expect_s3_class(mfr3$get_output(), "growth_day")
})
