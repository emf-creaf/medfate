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

test_that("Test stand metrics",{
  expect_type(stand_basalArea(exampleforestMED), "double")
  expect_type(stand_foliarBiomass(exampleforestMED, SpParamsMED), "double")
  expect_type(stand_fuelLoading(exampleforestMED, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(exampleforestMED, SpParamsMED), "double")
  expect_type(stand_LAI(exampleforestMED, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(exampleforestMED), "double")
  expect_type(stand_dominantTreeSpecies(exampleforestMED, SpParamsMED), "character")
  expect_type(stand_dominantTreeHeight(exampleforestMED), "double")
  expect_type(stand_hartBeckingIndex(exampleforestMED), "double")
  expect_type(stand_quadraticMeanTreeDiameter(exampleforestMED), "double")
})
