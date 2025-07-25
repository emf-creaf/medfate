library(medfate)

data(exampleforest)
data(SpParamsMED)


data(examplemeteo)
#Prepare a two-year meteorological data with half precipitation during 
#the second year
meteo2001 <- examplemeteo
meteo2002 <- examplemeteo
meteo2002$Precipitation <- meteo2002$Precipitation/2
meteo2002$dates <- seq(as.Date("2002-01-01"), 
                       as.Date("2002-12-31"), by="day")
meteo_01_02 <- rbind(meteo2001, meteo2002)
meteo_01_02_B <- meteo_01_02
meteo_01_02_B$dates <- as.Date(meteo_01_02_B$dates)
row.names(meteo_01_02_B) <- NULL

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control <- defaultControl("Granier")
control$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("fordyn can be run in example and empty forests",{
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
  expect_s3_class(fordyn(emptyforest(), examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
})

test_that("fordyn can be run using species codes",{
  f <- exampleforest
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  expect_s3_class(fordyn(f, examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
})
test_that("fordyn can be run using single soilDomains",{
  control_single <- control
  control_single$soilDomains <- "single"
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02, control_single,
                         latitude = 41.82592, elevation = 100), "fordyn")
})

test_that("fordyn can be run using dual soilDomains",{
  control_dual <- control
  control_dual$soilDomains <- "dual"
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02, control_dual,
                         latitude = 41.82592, elevation = 100), "fordyn")
})
test_that("fordyn can be run using partial rhizosphere overlap",{
  control_partial_overlap <- control
  control_partial_overlap$rhizosphereOverlap <- "partial"
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02, control_partial_overlap,
                         latitude = 41.82592, elevation = 100), "fordyn")
})
test_that("fordyn can be run in example and empty forests using management",{
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100,
                         management_function = defaultManagementFunction,
                         management_args = defaultManagementArguments()), "fordyn")
  expect_s3_class(fordyn(emptyforest(), examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100,
                         management_function = defaultManagementFunction,
                         management_args = defaultManagementArguments()), "fordyn")
})

test_that("fordyn can be run using dates as columns",{
  expect_s3_class(fordyn(exampleforest, examplesoil, 
                         SpParamsMED, meteo_01_02_B, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
})
