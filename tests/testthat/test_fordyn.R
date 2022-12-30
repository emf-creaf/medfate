library(medfate)

data(exampleforestMED)
data(SpParamsMED)


data(examplemeteo)
#Prepare a two-year meteorological data with half precipitation during 
#the second year
meteo2001 <- examplemeteo
meteo2002 <- examplemeteo
meteo2002$Precipitation <- meteo2002$Precipitation/2
row.names(meteo2002) <- seq(as.Date("2002-01-01"), 
                            as.Date("2002-12-31"), by="day")
meteo_01_02 <- rbind(meteo2001, meteo2002)

#Load example plot plant data
data(exampleforestMED)

#Default species parameterization
data(SpParamsMED)

#Initialize control parameters
control <- defaultControl("Granier")
control$verbose <- FALSE

#Initialize soil with default soil params (4 layers)
examplesoil <- soil(defaultSoilParams(4))

test_that("fordyn can be run in example and empty forests",{
  expect_s3_class(fordyn(exampleforestMED, examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
  expect_s3_class(fordyn(emptyforest(), examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
})

test_that("fordyn can be run using species codes",{
  f <- exampleforestMED
  f$treeData$Species <- c(148, 168)
  f$shrubData$Species <- 165
  expect_s3_class(fordyn(f, examplesoil, 
                         SpParamsMED, meteo_01_02, control,
                         latitude = 41.82592, elevation = 100), "fordyn")
})

test_that("fordyn can be run in example and empty forests using management",{
  expect_s3_class(fordyn(exampleforestMED, examplesoil, 
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
