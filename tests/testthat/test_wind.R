#Default species parameterization
data(SpParamsMED)

#Load example plot plant data
data(exampleforest)

#Canopy height (in m)
h= max(exampleforest$treeData$Height/100) 
d0 = 0.67*h
z0 = 0.08*h

#Height values (cm)
z = seq(50,1000, by=50)
zm = z/100 # (in m)

# Leaf area density
lad = vprofile_leafAreaDensity(exampleforest, SpParamsMED, draw = FALSE,
                               z = c(0,z))

# Effective drag
Cd = 0.2
Cx = Cd*lad

# canopy turbulence model
test_that("wind canopy turbulence can be run ",{
  expect_s3_class(wind_canopyTurbulenceModel(zm, Cx,h,d0,z0), "data.frame")
})
