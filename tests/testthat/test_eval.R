library(medfate)

data(examplemeteo)
data(exampleforest)
data(SpParamsMED)
examplesoil = defaultSoilParams(4)
control = defaultControl("Granier")
control$verbose = FALSE
data(exampleobs)
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)

test_that("evaluation can be performed",{
  expect_type(evaluation_stats(S1, exampleobs), "double")
  expect_type(evaluation_metric(S1, exampleobs, metric="NSE"), "double")
  expect_s3_class(evaluation_plot(S1, exampleobs), "ggplot")
})

test_that("evaluation for first soil layer can be performed",{
  names(exampleobs)[2] <- "SWC.1"
  expect_type(evaluation_stats(S1, exampleobs, type="SWC.1"), "double")
  expect_type(evaluation_metric(S1, exampleobs, type="SWC.1", metric="NSE"), "double")
  expect_s3_class(evaluation_plot(S1, exampleobs, type="SWC.1"), "ggplot")
})