library(medfate)

data(examplemeteo)

control <- defaultControl("Granier")
control$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("aspwb can be run in example data",{
  expect_s3_class(aspwb(aspwbInput(0.75, control, examplesoil), 
                       examplemeteo[1:10,],
                       latitude = 41.82592, elevation = 100), "aspwb")
})

