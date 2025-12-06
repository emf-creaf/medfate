library(medfate)

data(exampleforest)
data(SpParamsMED)
data(poblet_trees)

exampleforest_minimal <- exampleforest
exampleforest_minimal$shrubData <- exampleforest_minimal$shrubData[numeric(0),,drop = FALSE]
exampleforest_minimal$herbCover <- NULL
exampleforest_minimal$herbHeight <- NULL
exampleforest_minimal$seedBank <- NULL
exampleforest_minimal$seedlingBank <- NULL

emptyforest <- emptyforest()

forest_herbs <- exampleforest
forest_herbs$herbData <- forest_herbs$shrubData

test_that("Empty forests can be created",{
  expect_s3_class(emptyforest(), "forest")
  expect_s3_class(emptyforest(addcolumns = "LAI"), "forest")
  expect_s3_class(emptyforest(addcolumns = c("LAI", "Age")), "forest")
  expect_error(emptyforest(addcolumns = "kk"))
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

test_that("Can produce all vertical profiles with exampleforest",{
  expect_s3_class(vprofile_rootDistribution(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest, SpParamsMED, byCohorts = TRUE), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest, SpParamsMED, byCohorts = TRUE, includeHerbs = TRUE), "ggplot")
  expect_s3_class(vprofile_PARExtinction(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(exampleforest, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(exampleforest, SpParamsMED), "ggplot")
})

test_that("Can produce all vertical profiles with exampleforest2",{
  expect_s3_class(vprofile_rootDistribution(exampleforest2, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(exampleforest2, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest2, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest2, SpParamsMED, byCohorts = TRUE), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest2, SpParamsMED, byCohorts = TRUE, includeHerbs = TRUE), "ggplot")
  expect_s3_class(vprofile_PARExtinction(exampleforest2, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(exampleforest2, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(exampleforest2, SpParamsMED), "ggplot")
})

test_that("Can produce all vertical profiles with minimal forest",{
  expect_s3_class(vprofile_rootDistribution(exampleforest_minimal, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(exampleforest_minimal, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest_minimal, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest_minimal, SpParamsMED, byCohorts = TRUE), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(exampleforest_minimal, SpParamsMED, byCohorts = TRUE, includeHerbs = TRUE), "ggplot")
  expect_s3_class(vprofile_PARExtinction(exampleforest_minimal, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(exampleforest_minimal, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(exampleforest_minimal, SpParamsMED), "ggplot")
})

test_that("Can produce all vertical profiles with herb forest",{
  expect_s3_class(vprofile_rootDistribution(forest_herbs, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_fuelBulkDensity(forest_herbs, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(forest_herbs, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(forest_herbs, SpParamsMED, byCohorts = TRUE), "ggplot")
  expect_s3_class(vprofile_leafAreaDensity(forest_herbs, SpParamsMED, byCohorts = TRUE, includeHerbs = TRUE), "ggplot")
  expect_s3_class(vprofile_PARExtinction(forest_herbs, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_SWRExtinction(forest_herbs, SpParamsMED), "ggplot")
  expect_s3_class(vprofile_windExtinction(forest_herbs, SpParamsMED), "ggplot")
})
test_that("Cannot produce vertical profiles with empty forest",{
  expect_error(vprofile_rootDistribution(emptyforest, SpParamsMED))
  expect_error(vprofile_fuelBulkDensity(emptyforest, SpParamsMED))
  expect_error(vprofile_leafAreaDensity(emptyforest, SpParamsMED))
  expect_error(vprofile_PARExtinction(emptyforest, SpParamsMED))
  expect_error(vprofile_SWRExtinction(emptyforest, SpParamsMED))
  expect_error(vprofile_windExtinction(emptyforest, SpParamsMED))
})

test_that("Test forest summary",{
  expect_s3_class(summary(exampleforest, SpParamsMED), "summary.forest")
  expect_s3_class(summary(exampleforest2, SpParamsMED), "summary.forest")
  expect_s3_class(summary(exampleforest_minimal, SpParamsMED), "summary.forest")
  expect_s3_class(summary(emptyforest, SpParamsMED), "summary.forest")
  expect_s3_class(summary(forest_herbs, SpParamsMED), "summary.forest")
})

test_that("Test forest simplification",{
  f <- exampleforest
  f$treeData$Z100 <- c(1500, 1900)
  f$shrubDataZ100 <- 800
  expect_s3_class(forest_reduceToDominant(exampleforest, SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(exampleforest2, SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(exampleforest_minimal, SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(emptyforest(), SpParamsMED), "forest")
  expect_s3_class(forest_reduceToDominant(f, SpParamsMED), "forest")
  expect_s3_class(forest_mergeTrees(exampleforest), "forest")
  expect_error(forest_mergeTrees(exampleforest2)) # Missing values
  expect_s3_class(forest_mergeTrees(exampleforest_minimal), "forest")
  expect_s3_class(forest_mergeTrees(emptyforest()), "forest")
  expect_s3_class(forest_mergeTrees(emptyforest(addcolumns = "LAI")), "forest")
  expect_s3_class(forest_mergeTrees(f), "forest")
  expect_s3_class(forest_mergeShrubs(exampleforest), "forest")
  expect_error(forest_mergeShrubs(exampleforest2)) # Missing values
  expect_s3_class(forest_mergeShrubs(exampleforest_minimal), "forest")
  expect_s3_class(forest_mergeShrubs(emptyforest()), "forest")
  expect_s3_class(forest_mergeShrubs(emptyforest(addcolumns = "LAI")), "forest")
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

test_that("Test stand metrics with exampleforest",{
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

test_that("Test stand metrics with exampleforest2",{
  expect_type(stand_basalArea(exampleforest2), "double")
  expect_type(stand_foliarBiomass(exampleforest2, SpParamsMED), "double")
  expect_type(stand_fuelLoading(exampleforest2, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(exampleforest2, SpParamsMED), "double")
  expect_type(stand_LAI(exampleforest2, SpParamsMED), "double")
  expect_type(stand_shrubVolume(exampleforest2, SpParamsMED), "double")
  expect_error(stand_dominantTreeDiameter(exampleforest2)) # Missing values
  expect_type(stand_dominantTreeSpecies(exampleforest2, SpParamsMED), "character")
  expect_error(stand_treeDensity(exampleforest2)) # Missing values
  expect_error(stand_meanTreeHeight(exampleforest2)) # Missing values
  expect_error(stand_dominantTreeHeight(exampleforest2)) # Missing values
  expect_error(stand_hartBeckingIndex(exampleforest2)) # Missing values
  expect_error(stand_quadraticMeanTreeDiameter(exampleforest2)) # Missing values
})

test_that("Test stand metrics with empty forest",{
  expect_type(stand_basalArea(emptyforest), "double")
  expect_type(stand_foliarBiomass(emptyforest, SpParamsMED), "double")
  expect_type(stand_fuelLoading(emptyforest, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(emptyforest, SpParamsMED), "double")
  expect_type(stand_LAI(emptyforest, SpParamsMED), "double")
  expect_type(stand_shrubVolume(emptyforest, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(emptyforest), "double")
  expect_type(stand_dominantTreeSpecies(emptyforest, SpParamsMED), "character")
  expect_type(stand_treeDensity(emptyforest), "double")
  expect_type(stand_meanTreeHeight(emptyforest), "double")
  expect_type(stand_dominantTreeHeight(emptyforest), "double")
  expect_type(stand_hartBeckingIndex(emptyforest), "double")
  expect_type(stand_quadraticMeanTreeDiameter(emptyforest), "double")
})

test_that("Test stand metrics with minimal data",{
  expect_type(stand_basalArea(exampleforest_minimal), "double")
  expect_type(stand_foliarBiomass(exampleforest_minimal, SpParamsMED), "double")
  expect_type(stand_fuelLoading(exampleforest_minimal, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(exampleforest_minimal, SpParamsMED), "double")
  expect_type(stand_LAI(exampleforest_minimal, SpParamsMED), "double")
  expect_type(stand_shrubVolume(exampleforest_minimal, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(exampleforest_minimal), "double")
  expect_type(stand_dominantTreeSpecies(exampleforest_minimal, SpParamsMED), "character")
  expect_type(stand_treeDensity(exampleforest_minimal), "double")
  expect_type(stand_meanTreeHeight(exampleforest_minimal), "double")
  expect_type(stand_dominantTreeHeight(exampleforest_minimal), "double")
  expect_type(stand_hartBeckingIndex(exampleforest_minimal), "double")
  expect_type(stand_quadraticMeanTreeDiameter(exampleforest_minimal), "double")
})

test_that("Test stand metrics with herb forest",{
  expect_type(stand_basalArea(forest_herbs), "double")
  expect_type(stand_foliarBiomass(forest_herbs, SpParamsMED), "double")
  expect_type(stand_fuelLoading(forest_herbs, SpParamsMED), "double")
  expect_type(stand_foliarBiomass(forest_herbs, SpParamsMED), "double")
  expect_type(stand_LAI(forest_herbs, SpParamsMED), "double")
  expect_type(stand_shrubVolume(forest_herbs, SpParamsMED), "double")
  expect_type(stand_dominantTreeDiameter(forest_herbs), "double")
  expect_type(stand_dominantTreeSpecies(forest_herbs, SpParamsMED), "character")
  expect_type(stand_treeDensity(forest_herbs), "double")
  expect_type(stand_meanTreeHeight(forest_herbs), "double")
  expect_type(stand_dominantTreeHeight(forest_herbs), "double")
  expect_type(stand_hartBeckingIndex(forest_herbs), "double")
  expect_type(stand_quadraticMeanTreeDiameter(forest_herbs), "double")
})
test_that("Test stand metrics with missing tree data",{
  f <-exampleforest
  f$treeData<- f$treeData[numeric(0), ,drop = FALSE]
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
  f$shrubData<- f$shrubData[numeric(0), ,drop = FALSE]
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


test_that("Test stand metrics with missing tree and shrub data",{
  f <-exampleforest
  f$treeData<- f$treeData[numeric(0), ,drop = FALSE]
  f$shrubData<- f$shrubData[numeric(0), ,drop = FALSE]
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