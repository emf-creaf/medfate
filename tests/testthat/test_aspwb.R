library(medfate)

data(examplemeteo)

control <- defaultControl("Granier")
control$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)
examplemeteo2 <- examplemeteo
row.names(examplemeteo2) <- as.character(examplemeteo2$dates)
examplemeteo2$dates <- NULL
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

test_that("aspwb can be run in example data",{
  x1 <- aspwbInput(0.75, control, examplesoil)
  W_ini <- rlang::duplicate(x1$soil$W)
  expect_s3_class(aspwb_day(x1, 
                            date, meteovec,
                            latitude = 41.82592, elevation = 100, modifyInput = FALSE), "aspwb_day")
  expect_equal(W_ini, x1$soil$W)
  expect_s3_class(aspwb_day(x1, 
                            date, meteovec,
                            latitude = 41.82592, elevation = 100, modifyInput = TRUE), "aspwb_day")
  expect_false(all(W_ini == x1$soil$W))
  expect_s3_class(aspwb(aspwbInput(0.75, control, examplesoil), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "aspwb")
})

