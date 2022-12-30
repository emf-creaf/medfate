library(medfate)

data(exampleforestMED)
data(SpParamsMED)


test_that("Can produce all vertical profiles",{
  expect_s3_class(vprofile_rootDistribution(exampleforestMED, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(exampleforestMED, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforestMED, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_PARExtinction(exampleforestMED, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(exampleforestMED, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(exampleforestMED, SpParamsMED), "ggplot")
})

test_that("Test forest summary",{
  expect_s3_class(summary(exampleforestMED, SpParamsMED), "summary.forest")
})

test_that("Test forest merging",{
  expect_s3_class(forest_mergeTrees(exampleforestMED), "forest")
  expect_s3_class(forest_mergeShrubs(exampleforestMED), "forest")
})

