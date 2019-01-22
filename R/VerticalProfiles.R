vprofile.LeafAreaDensity<-function(x, SpParams, z = NULL, gdd = NA, byCohorts = FALSE,
                                   bySpecies = FALSE, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant.Height(x))/100)*100 , by=10)
  if(!byCohorts) {
    lai = .LAIprofile(z,x, SpParams, gdd)
    if(draw) {
      plot(lai, z[-1], type="l", xlab="Leaf area density (m2/m2)", ylab="Height (cm)")
    }
  } else {
    cohortnames = plant.ID(x, SpParams)
    lai = .LAIdistribution(z,x, SpParams, gdd)
    if(bySpecies) {
      spnames = plant.SpeciesName(x, SpParams)
      lai = t(apply(lai,1, tapply, spnames, sum, na.rm=T))
      cohortnames = colnames(lai)
    } 
    
    matplot(lai, z[-1], 
            type="l", xlab="Leaf area density (m2/m2)", 
            ylab="Height (cm)", lty=1:length(cohortnames), col = 1:length(cohortnames))
    legend("topright", legend = cohortnames, lty=1:length(cohortnames), 
           col = 1:length(cohortnames), bty="n")
  }
  if(draw) invisible(lai)
  else return(lai)
}

vprofile.RootDistribution<-function(x, d = NULL, draw = TRUE) {
  if(is.null(d)){
    zmax = ceiling(max(c(max(x$shrubData$Z), max(x$treeData$Z95)))/100)*100 
    d = rep(10,1+(zmax/10))
    z = seq(0,zmax, by=10)
  } else {
    zmax = sum(d)
    z = numeric(length(d))
    z[1] = 0
    for(i in 2:length(z)) z[i] = z[i-1] + d[i-1]
  }
  cohortnames = plant.ID(x)
  rd = .rootDistribution(d,x)
  matplot(t(rd)*10, z, 
          type="l", xlab="Percentage of fine roots/mm", 
          ylab="Depth (mm)", lty=1:length(cohortnames), col = 1:length(cohortnames),
          ylim=c(zmax,0))
  legend("bottomright", legend = cohortnames, lty=1:length(cohortnames), 
         col = 1:length(cohortnames), bty="n")
  if(draw) invisible(rd)
  else return(rd)
}

vprofile.FuelBulkDensity<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant.Height(x))/100)*100 , by=10)
  wfp = .woodyFuelProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(wfp, z[-1], type="l", xlab="Bulk density (kg/m3)", ylab="Height (cm)")
  }
  if(draw) invisible(wfp)
  else return(wfp)
}
vprofile.PARExtinction<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant.Height(x))/100)*100 , by=10)
  pep = .parExtinctionProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(pep, z, type="l", xlab="Percentage of PAR available", ylab="Height (cm)",
         xlim =c(0,100))
  }
  if(draw) invisible(pep)
  else return(pep)
}
vprofile.SWRExtinction<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant.Height(x))/100)*100 , by=10)
  pep = .swrExtinctionProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(pep, z, type="l", xlab="Percentage of SWR available", ylab="Height (cm)",
         xlim =c(0,100))
  }
  if(draw) invisible(pep)
  else return(pep)
}
vprofile.WindExtinction<-function(x, SpParams, wind20H, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant.Height(x))/100)*100 , by=10)
  fls = fuel.Stratification(x, SpParams, gdd);
  LAIc = fls$canopyLAI;
  canopyHeight = fls$canopyTopHeight;
  wep = .windExtinctionProfile(z, wind20H, LAIc, canopyHeight);
  if(draw) {
    plot(wep, z, type="l", xlab="Wind speed (m/s)", ylab="Height (cm)",
         xlim=c(0,max(wep)))
    abline(v=wind20H, col="gray", lty=2)
    abline(h=canopyHeight, col="gray", lty=2)
  }
  if(draw) invisible(wep)
  else return(wep)
}