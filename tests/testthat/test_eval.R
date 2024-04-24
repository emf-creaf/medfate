library(medfate)

data(examplemeteo)
data(exampleforest)
data(SpParamsMED)
examplesoil = soil(defaultSoilParams(4))
control = defaultControl("Granier")
control$verbose = FALSE
data(exampleobs)

test_that("evaluation can be performed",{
  x1 = forest2spwbInput(exampleforest,examplesoil, SpParamsMED, control)
  S1<-spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
  expect_type(evaluation_stats(S1, exampleobs), "double")
  expect_type(evaluation_metric(S1, exampleobs, metric="NSE"), "double")
  expect_s3_class(evaluation_plot(S1, exampleobs), "ggplot")
})
