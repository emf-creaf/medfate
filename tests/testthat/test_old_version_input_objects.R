library(medfate)

data(exampleforest)
data(SpParamsMED)
data(examplemeteo)
d <- 100
meteovec <- unlist(examplemeteo[d,-1])
date <- as.character(examplemeteo$dates[d])

current_version <- as.character(packageVersion("medfate"))
older_versions <- c("4.8.3","4.8.4")

test_that("spwb(_day) and growth(_day) can be run from stored inputs coming from older versions",{
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  
  # Generate current version objects 
  data(exampleforest)
  data(SpParamsMED)
  examplesoil <- defaultSoilParams(4)
  control <- defaultControl("Granier")
  control$verbose <- FALSE
  x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
  saveRDS(x1, file = paste0("initialized_objects/spwbInput_", current_version,".rds"))
  x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)
  saveRDS(x1, file = paste0("initialized_objects/growthInput_", current_version,".rds"))
  
  # Test older version objects 
  for(ver in older_versions) {
    x1 <- readRDS(paste0("initialized_objects/spwbInput_", ver,".rds"))  
    S1 <- spwb(x1, 
               examplemeteo[1:10,],
               latitude = 41.82592, elevation = 100)
    expect_s3_class(S1, "spwb")
    expect_equal(S1$spwbOutput$version, as.character(packageVersion("medfate")))
    
    x1 <- readRDS(paste0("initialized_objects/growthInput_", ver,".rds")) 
    G1 <- growth(x1, 
                 examplemeteo[1:10,],
                 latitude = 41.82592, elevation = 100)
    expect_s3_class(G1, "growth")
    expect_equal(G1$growthOutput$version, as.character(packageVersion("medfate")))
  }
  
  for(ver in older_versions) {
    x1 <- readRDS(paste0("initialized_objects/spwbInput_", ver,".rds"))  
    s1 <- medfate::spwb_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
    expect_s3_class(s1, "spwb_day")
    
    x1 <- readRDS(paste0("initialized_objects/growthInput_", ver,".rds")) 
    g1 <- medfate::growth_day(x1, date, meteovec, latitude = 41.82592, elevation = 100, slope=0, aspect=0, modifyInput = FALSE)
    expect_s3_class(g1, "growth_day")
  }
})


