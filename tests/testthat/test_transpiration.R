library(medfate)

data(exampleforest)

data(SpParamsMED)
data(examplemeteo)

#Define soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE
control_sperry <- defaultControl("Sperry")
control_sperry$verbose <- FALSE
control_sureau <- defaultControl("Sureau")
control_sureau$verbose <- FALSE

#Initialize input
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control_granier)

test_that("transpiration granier can be run ",{
  expect_type(transp_transpirationGranier(x1, examplemeteo, 1, 
                                    latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
                                    modifyInput = FALSE), "list")
})
