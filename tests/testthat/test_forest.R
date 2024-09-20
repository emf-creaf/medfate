library(medfate)

data(exampleforest)
data(SpParamsMED)
data(poblet_trees)

test_that("Empty forests can be created",{
  expect_s3_class(emptyforest(), "forest")
})

test_that("Forest object can be build from data frames",{
  f <- emptyforest()
  x <- subset(poblet_trees, Plot.Code=="POBL_CTL")
  sampled_area <- pi*15^2
  mapping_x <- c("Species.name" = "Species", "DBH" = "Diameter.cm")
  species <- c("Erica arborea","Cistus albidus", "Erica arborea", "Cistus albidus", "Cistus albidus")
  H <- c(200,50,100,40,30)
  D1 <- c(140,40,100, 35,30)
  D2 <- D1
  y <- data.frame(species, H, D1, D2)
  mapping_y <- c("Species.name"= "species", "Height" ="H", "D1", "D2")
  
  expect_s3_class(forest_mapWoodyTables(x, y,
                             mapping_x = mapping_x, mapping_y = mapping_y,
                             SpParams = SpParamsMED,
                             plot_size_x = sampled_area, plot_size_y = 4),"forest")
})

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

test_that("Test forest simplification",{
  f <- exampleforest
  f$treeData$Z100 <- c(1500, 1900)
  f$shrubDataZ100 <- 800
  expect_s3_class(forest_reduceToDominant(exampleforest, SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(emptyforest(), SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(f, SpParamsMED), "forest")
  expect_s3_class(forest_mergeTrees(exampleforest), "forest")
  expect_s3_class(forest_mergeTrees(emptyforest()), "forest")
  expect_s3_class(forest_mergeTrees(f), "forest")
  expect_s3_class(forest_mergeShrubs(exampleforest), "forest")
  expect_s3_class(forest_mergeShrubs(emptyforest()), "forest")
  expect_s3_class(forest_mergeShrubs(f), "forest")
  f2 <- exampleforest
  f2$shrubData <- rbind(f2$shrubData, f2$shrubData, f2$shrubData)
  f2$shrubData$ObsID <- c(NA,"1", NA)
  f2$treeData$ObsID <- c("1", NA)
  expect_s3_class(forest_mergeTrees(f2), "forest")
  expect_s3_class(forest_mergeShrubs(f2), "forest")
  expect_s3_class(forest_mergeTrees(f2, keepCohortsWithObsID = TRUE), "forest")
  expect_s3_class(forest_mergeShrubs(f2, keepCohortsWithObsID = TRUE), "forest")
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