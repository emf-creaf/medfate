library(medfate)

data(exampleforest)

forest_herbs <- exampleforest
forest_herbs$herbData <- forest_herbs$shrubData

data(SpParamsMED)
data(examplemeteo)
examplemeteo2 <- examplemeteo
row.names(examplemeteo2) <- as.character(examplemeteo2$dates)
examplemeteo2$dates <- NULL
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

control_granier <- defaultControl("Granier")
control_granier$verbose <- FALSE
control_sperry <- defaultControl("Sperry")
control_sperry$verbose <- FALSE
control_sureau <- defaultControl("Sureau")
control_sureau$verbose <- FALSE


#Initialize soil with default soil params (4 layers)
examplesoil <- defaultSoilParams(4)

test_that("spwb_day can be run with missing meteo data",{
  
  meteovec_mis <- meteovec
  meteovec_mis[4:7] <- NA
  
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  expect_s3_class(suppressWarnings(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau), 
                                    date, meteovec_mis, latitude = 41.82592, elevation = 100, slope=0, aspect=0)), "spwb_day")
  
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["MinTemperature"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                    date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["MaxTemperature"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                 date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
  meteovec_mis2 <- meteovec_mis
  meteovec_mis2["Precipitation"] <- NA
  expect_error(medfate::spwb_day(spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier), 
                                 date, meteovec_mis2, latitude = 41.82592, elevation = 100, slope=0, aspect=0))
})

test_that("spwb_day can be run with truncated root distribution",{
  exampleforest$shrubData$Z50 <- 100
  exampleforest$shrubData$Z95 <- 200
  control_granier$truncateRootDistribution<- TRUE
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_s3_class(spwb_day(x1, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sperry$truncateRootDistribution<- TRUE
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(spwb_day(x2, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sureau$truncateRootDistribution<- TRUE
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_s3_class(spwb_day(x3, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
})

test_that("spwb_day can be run with truncated root distribution and rhizosphere overlap",{
  exampleforest$shrubData$Z50 <- 100
  exampleforest$shrubData$Z95 <- 200
  control_granier$truncateRootDistribution<- TRUE
  control_granier$rhizosphereOverlap <- "partial"
  x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_granier)
  expect_s3_class(spwb_day(x1, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sperry$truncateRootDistribution<- TRUE
  control_sperry$rhizosphereOverlap <- "partial"
  x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sperry)
  expect_s3_class(spwb_day(x2, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
  control_sureau$truncateRootDistribution<- TRUE
  control_sureau$rhizosphereOverlap <- "partial"
  x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control_sureau)
  expect_s3_class(spwb_day(x3, date, meteovec,
                           latitude = 41.82592, elevation = 100, slope=0, aspect=0), "spwb_day")
})



  