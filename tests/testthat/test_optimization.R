library(medfate)

data(exampleforest)
data(SpParamsMED)
data(examplemeteo)
#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("rock optimization can be run in example and empty forests",{
  testthat::expect_vector(utils_rockOptimization(exampleforest, examplesoil, SpParamsMED, defaultControl(), 
                                         meteo = examplemeteo,
                                         latitude = 41.82592, elevation = 100))
  testthat::expect_vector(utils_rockOptimization(emptyforest(), examplesoil, SpParamsMED, defaultControl(), 
                                                 meteo = examplemeteo,
                                                 latitude = 41.82592, elevation = 100))
})
