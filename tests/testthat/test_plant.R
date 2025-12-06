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
test_that("Tree basal area can be retrieved",{
  expect_type(plant_basalArea(exampleforest, SpParamsMED), "double")
  expect_type(plant_basalArea(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_basalArea(emptyforest, SpParamsMED), "double")
  expect_type(plant_basalArea(forest_herbs, SpParamsMED), "double")
})

test_that("Larger tree basal area can be retrieved",{
  expect_type(plant_largerTreeBasalArea(exampleforest, SpParamsMED), "double")
  expect_type(plant_largerTreeBasalArea(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_largerTreeBasalArea(emptyforest, SpParamsMED), "double")
  expect_type(plant_largerTreeBasalArea(forest_herbs, SpParamsMED), "double")
})

test_that("Plant individual area can be retrieved",{
  expect_type(plant_individualArea(exampleforest, SpParamsMED), "double")
  expect_type(plant_individualArea(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_individualArea(emptyforest, SpParamsMED), "double")
  expect_type(plant_individualArea(forest_herbs, SpParamsMED), "double")
})

test_that("Plant crown ratio can be retrieved",{
  expect_type(plant_crownRatio(exampleforest, SpParamsMED), "double")
  expect_type(plant_crownRatio(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_crownRatio(emptyforest, SpParamsMED), "double")
  expect_type(plant_crownRatio(forest_herbs, SpParamsMED), "double")
})

test_that("Plant crown base height can be retrieved",{
  expect_type(plant_crownBaseHeight(exampleforest, SpParamsMED), "double")
  expect_type(plant_crownBaseHeight(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_crownBaseHeight(emptyforest, SpParamsMED), "double")
  expect_type(plant_crownBaseHeight(forest_herbs, SpParamsMED), "double")
})

test_that("Plant crown length can be retrieved",{
  expect_type(plant_crownLength(exampleforest, SpParamsMED), "double")
  expect_type(plant_crownLength(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_crownLength(emptyforest, SpParamsMED), "double")
  expect_type(plant_crownLength(forest_herbs, SpParamsMED), "double")
})

test_that("Plant foliar biomass can be retrieved",{
  expect_type(plant_foliarBiomass(exampleforest, SpParamsMED), "double")
  expect_type(plant_foliarBiomass(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_foliarBiomass(emptyforest, SpParamsMED), "double")
  expect_type(plant_foliarBiomass(forest_herbs, SpParamsMED), "double")
})

test_that("Plant fuel loading can be retrieved",{
  expect_type(plant_fuelLoading(exampleforest, SpParamsMED), "double")
  expect_type(plant_fuelLoading(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_fuelLoading(emptyforest, SpParamsMED), "double")
  expect_type(plant_fuelLoading(forest_herbs, SpParamsMED), "double")
})

test_that("Plant equilibrium leaf litter can be retrieved",{
  expect_type(plant_equilibriumLeafLitter(exampleforest, SpParamsMED), "double")
  expect_type(plant_equilibriumLeafLitter(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_equilibriumLeafLitter(emptyforest, SpParamsMED), "double")
  expect_type(plant_equilibriumLeafLitter(forest_herbs, SpParamsMED), "double")
})

test_that("Plant equilibrium small branch litter can be retrieved",{
  expect_type(plant_equilibriumSmallBranchLitter(exampleforest, SpParamsMED), "double")
  expect_type(plant_equilibriumSmallBranchLitter(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_equilibriumSmallBranchLitter(emptyforest, SpParamsMED), "double")
  expect_type(plant_equilibriumSmallBranchLitter(forest_herbs, SpParamsMED), "double")
})


test_that("Plant phytovolume can be retrieved",{
  expect_type(plant_phytovolume(exampleforest, SpParamsMED), "double")
  expect_type(plant_phytovolume(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_phytovolume(emptyforest, SpParamsMED), "double")
  expect_type(plant_phytovolume(forest_herbs, SpParamsMED), "double")
})


test_that("Plant LAI can be retrieved",{
  expect_type(plant_LAI(exampleforest, SpParamsMED), "double")
  expect_type(plant_LAI(exampleforest_minimal, SpParamsMED), "double")
  expect_type(plant_LAI(emptyforest, SpParamsMED), "double")
  expect_type(plant_LAI(forest_herbs, SpParamsMED), "double")
})

test_that("Aboveground data frame can be build",{
  expect_s3_class(forest2aboveground(exampleforest, SpParamsMED), "data.frame")
  expect_s3_class(forest2aboveground(exampleforest_minimal, SpParamsMED), "data.frame")
  expect_s3_class(forest2aboveground(emptyforest, SpParamsMED), "data.frame")
  expect_s3_class(forest2aboveground(forest_herbs, SpParamsMED), "data.frame")
})