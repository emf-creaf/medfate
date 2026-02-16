#Load example plot plant data
data(exampleforest)

# Example of aboveground parameters taken from a forest
# described using LAI and crown ratio
data(exampleforest2)

#Default species parameterization
data(SpParamsMED)


# Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

test_that("WB_runner initializes and runs correctly for aspwb", {
  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  W_ini <- rlang::duplicate(xa$soil$W)
  # builds runner
  mfra <- new(WB_runner, xa, 41.82592, 100, 0, 0)
  expect_s4_class(mfra, "Rcpp_WB_runner")
  # runs and gets output
  mfra$run_day(date, meteovec, 0, NULL, NA)
  sa <- mfra$get_output()
  expect_s3_class(sa, "aspwb_day")
  mfra$update_input(xa)
  # Check that W has changed
  expect_false(all(xa$soil$W == W_ini)) 
  # Check that result is equal to a stand-alone run
  xa2 <- aspwbInput(0.75, defaultControl(), examplesoil)
  sa2 <- aspwb_day(xa2, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(xa2, xa)
  expect_equal(sa2, sa)
})
test_that("WB_runner initializes and runs correctly for spwb with granier", {
  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  W_ini <- rlang::duplicate(x1$soil$W)
  # builds runner
  mfr1 <- new(WB_runner, x1, 41.82592, 100, 0, 0)
  expect_s4_class(mfr1, "Rcpp_WB_runner")
  # runs and gets output
  mfr1$run_day(date, meteovec, 0, NULL, NA)
  s1 <- mfr1$get_output()
  expect_s3_class(s1, "spwb_day")
  # Check that W has changed
  mfr1$update_input(x1)
  expect_false(all(x1$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x12 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  s12 <- spwb_day(x12, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x12, x1)
  expect_equal(s12, s1)
})

test_that("WB_runner initializes and runs correctly for spwb with sperry", {
  controlSperry <- defaultControl("Sperry")
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSperry)
  W_ini <- rlang::duplicate(x2$soil$W)
  # builds runner
  mfr2 <- new(WB_runner, x2, 41.82592, 100, 0, 0)
  expect_s4_class(mfr2, "Rcpp_WB_runner")
  # runs and gets output
  mfr2$run_day(date, meteovec, 0, NULL, NA)
  s2 <- mfr2$get_output()
  expect_s3_class(s2, "spwb_day")
  # Check that W has changed
  mfr2$update_input(x2)
  expect_false(all(x2$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x22 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSperry)
  s22 <- spwb_day(x22, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x22, x2)
  expect_equal(s22, s2)
})
test_that("WB_runner initializes and runs correctly for spwb with sureau", {
  controlSureau <- defaultControl("Sureau")
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  W_ini <- rlang::duplicate(x3$soil$W)
  # builds runner
  mfr3 <- new(WB_runner, x3, 41.82592, 100, 0, 0)
  expect_s4_class(mfr3, "Rcpp_WB_runner")
  # runs and gets output
  mfr3$run_day(date, meteovec, 0, NULL, NA)
  s3<-mfr3$get_output()
  expect_s3_class(s3, "spwb_day")
  # Check that W has changed
  mfr3$update_input(x3)
  expect_false(all(x3$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x32 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  s32 <- spwb_day(x32, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x32, x3)
  expect_equal(s32, s3)
})

n <- 100
d <- 100
meteovec1 <- unlist(examplemeteo[d,-1])
meteo_vec <- vector("list", n)
date <- as.character(examplemeteo$dates[d])
for(i in 1:n) meteo_vec[[i]] = rlang::duplicate(meteovec1)
x_vec <- vector("list", n)
latitude_vec <- rep(41, n)
elevation_vec <- rep(100, n)
slope_vec <- rep(0, n)
aspect_vec <- rep(0, n)

test_that("WB_multiple_runner initializes and runs correctly with agriculture", {
  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  W_ini <- rlang::duplicate(xa$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(xa)
  # Build runner
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  # Run with parallelization
  Sys.setenv(ARMA_NUM_THREADS = "1")
  mfmr$run_day(date, meteo_vec, TRUE)
  expect_s3_class(mfmr$get_output_at(1), "aspwb_day")
  expect_s3_class(mfmr$get_output_at(2), "aspwb_day")
  sa <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  xa2 <- rlang::duplicate(xa)
  sa2 <- aspwb_day(xa2, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(xa2, x_vec[[2]])
  expect_equal(sa2, sa)
})
test_that("WB_multiple_runner initializes and runs correctly with granier", {
  controlGranier <- defaultControl("Granier")
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlGranier)
  W_ini <- rlang::duplicate(x1$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
  # Build runner
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  # Run with parallelization
  Sys.setenv(ARMA_NUM_THREADS = "1")
  mfmr$run_day(date, meteo_vec, TRUE)
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
  s1 <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x12 <- rlang::duplicate(x1)
  s12 <- spwb_day(x12, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x12, x_vec[[2]])
  expect_equal(s12, s1)
})
test_that("WB_multiple_runner initializes and runs correctly with sperry", {
  controlSperry <- defaultControl("Sperry")
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSperry)
  W_ini <- rlang::duplicate(x2$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(x2)
  # Build runner
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  # Run with parallelization
  Sys.setenv(ARMA_NUM_THREADS = "1")
  mfmr$run_day(date, meteo_vec, TRUE)
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
  s2 <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x22 <- rlang::duplicate(x2)
  s22 <- spwb_day(x22, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x22, x_vec[[2]])
  expect_equal(s22, s2)
})
test_that("WB_multiple_runner initializes and runs correctly with sureau", {
  controlSureau <- defaultControl("Sureau")
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, controlSureau)
  W_ini <- rlang::duplicate(x3$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(x3)
  # Build runner
  mfmr <- new(WB_multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_WB_multiple_runner")
  # Run with parallelization
  Sys.setenv(ARMA_NUM_THREADS = "1")
  mfmr$run_day(date, meteo_vec, TRUE)
  expect_s3_class(mfmr$get_output_at(1), "spwb_day")
  expect_s3_class(mfmr$get_output_at(2), "spwb_day")
  s3 <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  x32 <- rlang::duplicate(x3)
  expect_equal(x32, x3)
  s32 <- spwb_day(x32, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(x32, x_vec[[2]])
  expect_equal(s32, s2)
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
