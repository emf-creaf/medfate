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
  expect_type(stand_shrubVolume(exampleforest, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(exampleforest), "double")
  expect_type(stand_dominantTreeSpecies(exampleforest, SpParamsMED), "character")
  expect_type(stand_treeDensity(exampleforest), "double")
  expect_type(stand_meanTreeHeight(exampleforest), "double")
  expect_type(stand_dominantTreeHeight(exampleforest), "double")
  expect_type(stand_hartBeckingIndex(exampleforest), "double")
  expect_type(stand_quadraticMeanTreeDiameter(exampleforest), "double")
})

test_that("Test stand metrics with missing tree data",{
  f <-exampleforest
  f$treeData<- f$treeData[-c(1:nrow(exampleforest$treeData)), ,drop = FALSE]
  expect_type(stand_basalArea(f), "double")
  expect_type(stand_foliarBiomass(f, SpParamsMED), "double")
  expect_type(stand_fuelLoading(f, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(f, SpParamsMED), "double")
  expect_type(stand_shrubVolume(f, SpParamsMED), "double")
  expect_type(stand_LAI(f, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(f), "double")
  expect_type(stand_dominantTreeSpecies(f, SpParamsMED), "character")
  expect_type(stand_treeDensity(f), "double")
  expect_type(stand_meanTreeHeight(f), "double")
  expect_type(stand_dominantTreeHeight(f), "double")
  expect_type(stand_hartBeckingIndex(f), "double")
  expect_type(stand_quadraticMeanTreeDiameter(f), "double")
})

test_that("Test stand metrics with missing shrub data",{
  f <-exampleforest
  f$shrubData<- f$shrubData[-c(1:nrow(exampleforest$shrubData)), ,drop = FALSE]
  expect_type(stand_basalArea(f), "double")
  expect_type(stand_foliarBiomass(f, SpParamsMED), "double")
  expect_type(stand_fuelLoading(f, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(f, SpParamsMED), "double")
  expect_type(stand_shrubVolume(f, SpParamsMED), "double")
  expect_type(stand_LAI(f, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(f), "double")
  expect_type(stand_dominantTreeSpecies(f, SpParamsMED), "character")
  expect_type(stand_treeDensity(f), "double")
  expect_type(stand_meanTreeHeight(f), "double")
  expect_type(stand_dominantTreeHeight(f), "double")
  expect_type(stand_hartBeckingIndex(f), "double")
  expect_type(stand_quadraticMeanTreeDiameter(f), "double")
})