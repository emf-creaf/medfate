library(medfate)

data(examplemeteo)
data(exampleforest)
data(SpParamsMED)
examplesoil <- defaultSoilParams(4)
control <- defaultControl("Granier")
x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
start_profiler("/tmp/profile_granier.out")
S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
stop_profiler()

control <- defaultControl("Sperry")
x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
start_profiler("/tmp/profile_sperry.out")
S2 <- spwb(x2, examplemeteo, latitude = 41.82592, elevation = 100)
stop_profiler()

control <- defaultControl("Sureau")
x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
start_profiler("/tmp/profile_sureau.out")
S3 <- spwb(x3, examplemeteo, latitude = 41.82592, elevation = 100)
stop_profiler()


# google-pprof --gv /usr/bin/r /tmp/profile_granier.out 
# google-pprof --gv /usr/bin/r /tmp/profile_sperry.out 
# google-pprof --gv /usr/bin/r /tmp/profile_sureau.out 