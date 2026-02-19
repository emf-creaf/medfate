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

test_that("single_runner initializes and runs correctly for aspwb", {
  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  W_ini <- rlang::duplicate(xa$soil$W)
  # builds runner
  mfra <- new(runners$single_runner, xa, 41.82592, 100, 0, 0)
  expect_s4_class(mfra, "Rcpp_single_runner")
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

test_that("single_runner initializes and runs correctly for spwb", {
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    for(soilDomains in c("buckets", "single", "dual")) {
      ctl <- defaultControl(transpirationMode, soilDomains)
      x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, ctl)
      W_ini <- rlang::duplicate(x1$soil$W)
      # builds runner
      mfr1 <- new(runners$single_runner, x1, 41.82592, 100, 0, 0)
      expect_s4_class(mfr1, "Rcpp_single_runner")
      # runs and gets output
      mfr1$run_day(date, meteovec, 0, NULL, NA)
      s1 <- mfr1$get_output()
      expect_s3_class(s1, "spwb_day")
      # Check that W has changed
      mfr1$update_input(x1)
      expect_false(all(x1$soil$W == W_ini))
      # Check that result is equal to a stand-alone run
      x12 <- spwbInput(exampleforest, examplesoil, SpParamsMED, ctl)
      s12 <- spwb_day(x12, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
      expect_equal(x12, x1)
      expect_equal(s12, s1)
    }
  }
})

test_that("single_runner initializes and runs correctly for growth", {
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    for(soilDomains in c("buckets", "single", "dual")) {
      ctl <- defaultControl(transpirationMode)
      x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, ctl)
      W_ini <- rlang::duplicate(x1$soil$W)
      # builds runner
      mfr1 <- new(runners$single_runner, x1, 41.82592, 100, 0, 0)
      expect_s4_class(mfr1, "Rcpp_single_runner")
      # runs and gets output
      mfr1$run_day(date, meteovec, 0, NULL, NA)
      s1 <- mfr1$get_output()
      expect_s3_class(s1, "growth_day")
      # # Check that W has changed
      mfr1$update_input(x1)
      expect_false(all(x1$soil$W == W_ini))
      # # Check that result is equal to a stand-alone run
      x12 <- growthInput(exampleforest, examplesoil, SpParamsMED, ctl)
      s12 <- growth_day(x12, date, meteovec, 41.82592, 100, 0, 0, modifyInput = TRUE)
      expect_equal(x12, x1)
      expect_equal(s12, s1)
    }
  }
})


n <- 50
d <- 100
meteovec1 <- unlist(examplemeteo[d,-1])
meteo_vec <- vector("list", n)
date <- as.character(examplemeteo$dates[d])
for(i in 1:n) meteo_vec[[i]] = rlang::duplicate(meteovec1)
x_vec <- vector("list", n)
latitude_vec <- rep(41.82592, n)
elevation_vec <- rep(100, n)
slope_vec <- rep(0, n)
aspect_vec <- rep(0, n)

test_that("multiple_runner initializes and runs correctly with agriculture", {
  parallelize = FALSE
  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  W_ini <- rlang::duplicate(xa$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(xa)
  # Build runner
  mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_multiple_runner")
  # Run without parallelization
  mfmr$run_day(date, meteo_vec, parallelize)
  expect_s3_class(mfmr$get_output_at(1), "aspwb_day")
  expect_s3_class(mfmr$get_output_at(2), "aspwb_day")
  sa <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  xa2 <- rlang::duplicate(xa)
  sa2 <- aspwb_day(xa2, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(xa2, x_vec[[2]])
  expect_equal(sa2, sa)
})

test_that("multiple_runner initializes and runs correctly for spwb", {
  parallelize = FALSE
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    ctl <- defaultControl(transpirationMode)
    x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, ctl)
    W_ini <- rlang::duplicate(x1$soil$W)
    for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
    # Build runner
    mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
    expect_s4_class(mfmr, "Rcpp_multiple_runner")
    # Run without parallelization
    mfmr$run_day(date, meteo_vec, parallelize)
    expect_s3_class(mfmr$get_output_at(1), "spwb_day")
    expect_s3_class(mfmr$get_output_at(2), "spwb_day")
    s1 <- mfmr$get_output_at(2)
    # Check that W has changed
    mfmr$update_input_at(2, x_vec[[2]])
    expect_false(all(x_vec[[2]]$soil$W == W_ini))
    # Check that result is equal to a stand-alone run
    x12 <- rlang::duplicate(x1)
    s12 <- spwb_day(x12, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
    expect_equal(x12, x_vec[[2]])
    expect_equal(s12, s1)
  }
})

test_that("multiple_runner initializes and runs correctly for growth", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  parallelize = FALSE
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    ctl <- defaultControl(transpirationMode)
    x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, ctl)
    W_ini <- rlang::duplicate(x1$soil$W)
    for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
    # Build runner
    mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
    expect_s4_class(mfmr, "Rcpp_multiple_runner")
    # Run without parallelization
    mfmr$run_day(date, meteo_vec, parallelize)
    expect_s3_class(mfmr$get_output_at(1), "growth_day")
    expect_s3_class(mfmr$get_output_at(2), "growth_day")
    s1 <- mfmr$get_output_at(2)
    # Check that W has changed
    mfmr$update_input_at(2, x_vec[[2]])
    expect_false(all(x_vec[[2]]$soil$W == W_ini))
    # Check that result is equal to a stand-alone run
    x12 <- rlang::duplicate(x1)
    s12 <- growth_day(x12, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
    expect_equal(x12, x_vec[[2]])
    expect_equal(s12, s1)
    
  }
})

test_that("multiple_runner initializes and runs correctly with agriculture and parallelization", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  parallelize = TRUE
  xa <- aspwbInput(0.75, defaultControl(), examplesoil)
  W_ini <- rlang::duplicate(xa$soil$W)
  for(i in 1:n) x_vec[[i]] = rlang::duplicate(xa)
  # Build runner
  mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
  expect_s4_class(mfmr, "Rcpp_multiple_runner")
  # Run with parallelization
  Sys.setenv(ARMA_NUM_THREADS = "1")
  mfmr$run_day(date, meteo_vec, parallelize)
  expect_s3_class(mfmr$get_output_at(1), "aspwb_day")
  expect_s3_class(mfmr$get_output_at(2), "aspwb_day")
  sa <- mfmr$get_output_at(2)
  # Check that W has changed
  mfmr$update_input_at(2, x_vec[[2]])
  expect_false(all(x_vec[[2]]$soil$W == W_ini))
  # Check that result is equal to a stand-alone run
  xa2 <- rlang::duplicate(xa)
  sa2 <- aspwb_day(xa2, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
  expect_equal(xa2, x_vec[[2]])
  expect_equal(sa2, sa)
})

test_that("multiple_runner initializes and runs correctly for spwb and parallelization", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  parallelize = TRUE
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    ctl <- defaultControl(transpirationMode)
    x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, ctl)
    W_ini <- rlang::duplicate(x1$soil$W)
    for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
    # Build runner
    mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
    expect_s4_class(mfmr, "Rcpp_multiple_runner")
    # Run with parallelization
    Sys.setenv(ARMA_NUM_THREADS = "1")
    mfmr$run_day(date, meteo_vec, parallelize)
    expect_s3_class(mfmr$get_output_at(1), "spwb_day")
    expect_s3_class(mfmr$get_output_at(2), "spwb_day")
    s1 <- mfmr$get_output_at(2)
    # Check that W has changed
    mfmr$update_input_at(2, x_vec[[2]])
    expect_false(all(x_vec[[2]]$soil$W == W_ini))
    # Check that result is equal to a stand-alone run
    x12 <- rlang::duplicate(x1)
    s12 <- spwb_day(x12, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
    expect_equal(x12, x_vec[[2]])
    expect_equal(s12, s1)
  }
})

test_that("multiple_runner initializes and runs correctly for growth and parallelization", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  parallelize = TRUE
  for(transpirationMode in c("Granier", "Sperry", "Sureau")) {
    ctl <- defaultControl(transpirationMode)
    x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, ctl)
    W_ini <- rlang::duplicate(x1$soil$W)
    for(i in 1:n) x_vec[[i]] = rlang::duplicate(x1)
    # Build runner
    mfmr <- new(runners$multiple_runner, x_vec, latitude_vec, elevation_vec, slope_vec, aspect_vec)
    expect_s4_class(mfmr, "Rcpp_multiple_runner")
    # Run with parallelization
    Sys.setenv(ARMA_NUM_THREADS = "1")
    mfmr$run_day(date, meteo_vec, parallelize)
    expect_s3_class(mfmr$get_output_at(1), "growth_day")
    expect_s3_class(mfmr$get_output_at(2), "growth_day")
    s1 <- mfmr$get_output_at(2)
    # Check that W has changed
    mfmr$update_input_at(2, x_vec[[2]])
    expect_false(all(x_vec[[2]]$soil$W == W_ini))
    # Check that result is equal to a stand-alone run
    x12 <- rlang::duplicate(x1)
    s12 <- growth_day(x12, date, meteovec1, 41.82592, 100, 0, 0, modifyInput = TRUE)
    expect_equal(x12, x_vec[[2]])
    expect_equal(s12, s1)

  }
})
