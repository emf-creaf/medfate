library(testthat)
library(medfate)

data(exampleforest)
data(SpParamsMED)
data(examplemeteo)

versions <- c("4.8.3","4.8.4")

test_that("spwb and growth can be run from stored inputs coming from older versions",{
  for(ver in versions) {
    x1 <- readRDS(paste0("../initialized_objects/spwbInput_", ver,".rds"))  
    S1 <- spwb(x1, 
               examplemeteo[1:10,],
               latitude = 41.82592, elevation = 100)
    expect_s3_class(S1, "spwb")
    expect_equal(S1$spwbOutput$version, as.character(packageVersion("medfate")))
    
    g1 <- readRDS(paste0("../initialized_objects/growthInput_", ver,".rds")) 
    G1 <- growth(g1, 
                 examplemeteo[1:10,],
                 latitude = 41.82592, elevation = 100)
    expect_s3_class(G1, "growth")
    expect_equal(G1$growthOutput$version, as.character(packageVersion("medfate")))
  }
})
