# Default parameters
df_soil <- defaultSoilParams()

# Initializes soil
s <- soil(df_soil)

test_that("infiltration repartition", {
  expect_type(hydrology_infiltrationRepartition(100, s$widths, s$macro), "double")
})