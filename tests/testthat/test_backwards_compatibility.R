library(medfate)

data(exampleforest)
data(SpParamsMED)
data(examplemeteo)
examplemeteo2 <- examplemeteo
row.names(examplemeteo2) <- as.character(examplemeteo2$dates)
examplemeteo2$dates <- NULL
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

oldVer <- "4.8.3"
test_that("spwb can be run with missing or older versions", {
  x_granier <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- oldVer
  s <- spwb(x_granier, 
            examplemeteo[1:10,],
            latitude = 41.82592, elevation = 100)
  expect_s3_class(s, "spwb")
  expect_equal(s$spwbOutput$version, as.character(packageVersion("medfate")))
  x_granier$version <- NULL
  s <- spwb(x_granier, 
            examplemeteo[1:10,],
            latitude = 41.82592, elevation = 100)
  expect_s3_class(s, "spwb")
  expect_equal(s$spwbOutput$version, as.character(packageVersion("medfate")))
})

test_that("spwb_day can be run with missing or older versions", {
  x_granier <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- oldVer
  s <- medfate::spwb_day(x_granier, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_s3_class(s, "spwb_day")
  expect_equal(x_granier$version, as.character(packageVersion("medfate")))
  x_granier <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- NULL
  s <- medfate::spwb_day(x_granier, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_s3_class(s, "spwb_day")
  # If new field is added, then the pointer is new
  # expect_equal(x_granier, as.character(packageVersion("medfate")))
})

test_that("growth can be run with missing or older versions", {
  x_granier <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- oldVer
  s <- growth(x_granier, 
            examplemeteo[1:10,],
            latitude = 41.82592, elevation = 100)
  expect_s3_class(s, "growth")
  expect_equal(s$growthOutput$version, as.character(packageVersion("medfate")))
  x_granier$version <- NULL
  s <- growth(x_granier, 
            examplemeteo[1:10,],
            latitude = 41.82592, elevation = 100)
  expect_s3_class(s, "growth")
  expect_equal(s$growthOutput$version, as.character(packageVersion("medfate")))
})
test_that("growth_day can be run with missing or older versions", {
  x_granier <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- oldVer
  s <- medfate::growth_day(x_granier, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_s3_class(s, "growth_day")
  expect_equal(x_granier$version, as.character(packageVersion("medfate")))
  x_granier <- growthInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  x_granier$version <- NULL
  s <- medfate::growth_day(x_granier, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
  expect_s3_class(s, "growth_day")
  # If new field is added, then the pointer is new
  # expect_equal(x_granier, as.character(packageVersion("medfate")))
})