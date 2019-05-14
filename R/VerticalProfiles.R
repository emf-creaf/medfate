vprofile_leafAreaDensity<-function(x, SpParams = NULL, z = NULL, gdd = NA, byCohorts = FALSE,
                                   bySpecies = FALSE, draw = TRUE, legend = TRUE, xlim = NULL) {
  if(!(inherits(x,"data.frame") || inherits(x, "forest"))) stop("'x' should be of class 'forest' or 'data.frame'")
  if(inherits(x, "forest")) {
    if(is.null(SpParams)) stop("Please, provide 'SpParams' to calculate leaf area.")
    spnames = plant_speciesName(x, SpParams)
    x = forest2aboveground(x, SpParams, gdd)
  } else {
    if(any(!(c("LAI_expanded", "H", "CR", "SP") %in% names(x)))) {
      stop("Data frame should contain columns 'SP', 'LAI_expanded', 'H' and 'CR'")
    }
    spnames = x$SP
  }

  if(is.null(z)) z = seq(0, ceiling(max(x$H)/100)*100 +10, by=10)
  w = z[2:length(z)]- z[1:(length(z)-1)]
  if(!byCohorts) {
    lai = .LAIprofileVectors(z, x$LAI_expanded, x$H, x$CR)
    lai = 100*lai/w
    if(draw) {
      df = data.frame(lai = c(0,lai), z = z)
      g<-ggplot(df, aes(x=lai, y=z))+
        geom_path()+
        xlab("Leaf Area Density (m2/m3)")+
        ylab("Height (cm)")+
        theme_bw()
      if(!is.null(xlim)) g <- g + xlim(xlim)
    }
  } else {
    cohortnames = row.names(x)
    lai = .LAIdistributionVectors(z, x$LAI_expanded, x$H, x$CR)
    lai = 100*sweep(lai,1,w, "/")
    if(bySpecies) {
      lai = t(apply(lai,1, tapply, spnames, sum, na.rm=T))
      cohortnames = colnames(lai)
    } 
    if(draw) {
      lai = rbind(rep(0, ncol(lai)),lai)
      df = data.frame(lai = as.vector(lai), z = z,
                      Cohort = gl(length(cohortnames), nrow(lai), labels=cohortnames))
      g<-ggplot(df, aes(x=lai, y=z))+
        geom_path(aes(col=Cohort, linetype = Cohort))+
        xlab("Leaf Area Density (m2/m3)")+
        ylab("Height (cm)")+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")+
        theme_bw()
      if(!is.null(xlim)) g <- g + xlim(xlim)
    }
  }
  if(draw) return(g)
  else return(lai)
}

vprofile_rootDistribution<-function(x, SpParams, d = NULL, bySpecies = FALSE, draw = TRUE, xlim = NULL) {
  if(is.null(d)){
    zmax = 0
    if(nrow(x$shrubData)>0) zmax = max(zmax, ceiling(max(x$shrubData$Z95)/100)*100)
    if(nrow(x$treeData)>0) zmax = max(zmax, ceiling(max(x$treeData$Z95)/100)*100)
    d = rep(10,1+(zmax/10))
    z = seq(0,zmax, by=10)
  } else {
    zmax = sum(d)
    z = numeric(length(d))
    z[1] = 0
    for(i in 2:length(z)) z[i] = z[i-1] + d[i-1]
  }
  cohortnames = plant_ID(x)
  rd = .rootDistribution(d,x)
  rownames(rd) = cohortnames
  if(bySpecies) {
    spnames = plant_speciesName(x, SpParams)
    rd = apply(rd,2, tapply, spnames, mean, na.rm=T)
  } 
  if(draw) {
    return(.multiple_x(x=t(rd)*100, y=z/10, xlab="% of fine roots/cm",
                       ylab="Depth (cm)", ylim=c(zmax/10,0), xlim = xlim))
  } else return(rd)
}

vprofile_fuelBulkDensity<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant_height(x))/100)*100 , by=10)
  wfp = .woodyFuelProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(wfp, z[-1], type="l", xlab="Bulk density (kg/m3)", ylab="Height (cm)")
  }
  if(draw) invisible(wfp)
  else return(wfp)
}
vprofile_PARExtinction<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant_height(x))/100)*100 , by=10)
  pep = .parExtinctionProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(pep, z, type="l", xlab="Percentage of PAR available", ylab="Height (cm)",
         xlim =c(0,100))
  }
  if(draw) invisible(pep)
  else return(pep)
}
vprofile_SWRExtinction<-function(x, SpParams, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant_height(x))/100)*100 , by=10)
  pep = .swrExtinctionProfile(z,x, SpParams, gdd)
  if(draw) {
    plot(pep, z, type="l", xlab="Percentage of SWR available", ylab="Height (cm)",
         xlim =c(0,100))
  }
  if(draw) invisible(pep)
  else return(pep)
}
vprofile_windExtinction<-function(x, SpParams, wind20H, z = NULL, gdd = NA, draw = TRUE) {
  if(is.null(z)) z = seq(0, ceiling(max(plant_height(x))/100)*100 , by=10)
  fls = fuel_stratification(x, SpParams, gdd);
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