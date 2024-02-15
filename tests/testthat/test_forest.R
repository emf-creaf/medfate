library(medfate)

data(exampleforest)
data(SpParamsMED)


test_that("Can produce all vertical profiles",{
  expect_s3_class(vprofile_rootDistribution(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_PARExtinction(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(exampleforest, SpParamsMED), "ggplot")
})

test_that("Test forest summary",{
  expect_s3_class(summary(exampleforest, SpParamsMED), "summary.forest")
})

test_that("Test forest merging",{
  expect_s3_class(forest_mergeTrees(exampleforest), "forest")
  expect_s3_class(forest_mergeShrubs(exampleforest), "forest")
})

test_that("Test stand metrics",{
  expect_type(stand_basalArea(exampleforest), "double")
  expect_type(stand_foliarBiomass(exampleforest, SpParamsMED), "double")
  expect_type(stand_fuelLoading(exampleforest, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(exampleforest, SpParamsMED), "double")
  expect_type(stand_LAI(exampleforest, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(exampleforest), "double")
  expect_type(stand_dominantTreeSpecies(exampleforest, SpParamsMED), "character")
  expect_type(stand_dominantTreeHeight(exampleforest), "double")
  expect_type(stand_hartBeckingIndex(exampleforest), "double")
  expect_type(stand_quadraticMeanTreeDiameter(exampleforest), "double")
})
