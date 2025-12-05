library(medfate)

data(exampleforest)
data(SpParamsMED)

exampleforest_minimal <- exampleforest
exampleforest_minimal$shrubData <- exampleforest_minimal$shrubData[numeric(0),,drop = FALSE]
exampleforest_minimal$herbCover <- NULL
exampleforest_minimal$herbHeight <- NULL
exampleforest_minimal$seedBank <- NULL
exampleforest_minimal$seedlingBank <- NULL

emptyforest <- emptyforest()

forest_herbs <- exampleforest
forest_herbs$herbData <- forest_herbs$shrubData

test_that("Plant ID can be retrieved",{
  expect_type(plant_ID(exampleforest, SpParamsMED), "character")
  expect_type(plant_ID(exampleforest_minimal, SpParamsMED), "character")
  expect_type(plant_ID(emptyforest, SpParamsMED), "character")
  expect_type(plant_ID(forest_herbs, SpParamsMED), "character")
})

test_that("Plant species can be retrieved",{
  expect_type(plant_species(exampleforest, SpParamsMED), "integer")
  expect_type(plant_species(exampleforest_minimal, SpParamsMED), "integer")
  expect_type(plant_species(emptyforest, SpParamsMED), "integer")
  expect_type(plant_species(forest_herbs, SpParamsMED), "integer")
})

test_that("Plant species name can be retrieved",{
  expect_type(plant_speciesName(exampleforest, SpParamsMED), "character")
  expect_type(plant_speciesName(exampleforest_minimal, SpParamsMED), "character")
  expect_type(plant_speciesName(emptyforest, SpParamsMED), "character")
  expect_type(plant_speciesName(forest_herbs, SpParamsMED), "character")
})
test_that("Plant height can be retrieved",{
  expect_type(plant_height(exampleforest, SpParamsMED), "double")
  expect_type(plant_height(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_height(emptyforest, SpParamsMED), "double")
  expect_type(plant_height(forest_herbs, SpParamsMED), "double")
})
test_that("Plant density can be retrieved",{
  expect_type(plant_density(exampleforest, SpParamsMED), "double")
  expect_type(plant_density(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_density(emptyforest, SpParamsMED), "double")
  expect_type(plant_density(forest_herbs, SpParamsMED), "double")
})
test_that("Plant cover can be retrieved",{
  expect_type(plant_cover(exampleforest, SpParamsMED), "double")
  expect_type(plant_cover(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_cover(emptyforest, SpParamsMED), "double")
  expect_type(plant_cover(forest_herbs, SpParamsMED), "double")
})