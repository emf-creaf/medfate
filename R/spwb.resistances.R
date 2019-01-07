spwb.resistances<-function(x, cohort = 1, relative = FALSE, draw = FALSE, cumulative = FALSE, yearAxis = FALSE,  xlab = NULL, ylab=NULL) {
  
  VCroot_kmax = x$spwbInput$below$VCroot_kmax
  VGrhizo_kmax = x$spwbInput$below$VGrhizo_kmax
  VG_nc = x$soilInput$VG_n
  VG_alphac = x$soilInput$VG_alpha
  
  paramsTransp = x$spwbInput$paramsTransp
  VCroot_c = paramsTransp$VCroot_c
  VCroot_d = paramsTransp$VCroot_d
  VCstem_kmax = paramsTransp$VCstem_kmax
  VCstem_c = paramsTransp$VCstem_c
  VCstem_d = paramsTransp$VCstem_d
  VCleaf_kmax = paramsTransp$VCleaf_kmax
  VCleaf_c = paramsTransp$VCleaf_c
  VCleaf_d = paramsTransp$VCleaf_d
  
  psiLeaf = x$LeafPsi
  psiStem = x$StemPsi
  psiRoot = x$RootPsi
  PLCstem = x$PlantStress
  psiRhizo = x$RhizoPsi
  
  nlayers = length(VG_nc)
  psiSoil = x$Soil$psi.1
  if(nlayers>1) psiSoil = cbind(psiSoil, x$Soil$psi.2)
  if(nlayers>2) psiSoil = cbind(psiSoil, x$Soil$psi.3)
  if(nlayers>3) psiSoil = cbind(psiSoil, x$Soil$psi.4)
  if(nlayers>4) psiSoil = cbind(psiSoil, x$Soil$psi.5)
  
  nsteps = nrow(psiSoil)
  resmat = matrix(0, nrow=nsteps, ncol = 4)
  rownames(resmat) = rownames(psiStem)
  colnames(resmat) = c("Rhizosphere", "Root", "Stem", "Leaf")
  for(j in 1:nsteps) {
    rrow  = hydraulics.soilPlantResistances(psiSoil = psiSoil[j,],
                                            psiRhizo = psiRhizo[[cohort]][j,],
                                            psiStem = psiStem[j,cohort],
                                            PLCstem = PLCstem[j,cohort],
                                            psiLeaf = psiLeaf[j,cohort],
                                            VGrhizo_kmax[cohort,],VG_nc,VG_alphac,
                                            VCroot_kmax[cohort,], VCroot_c[cohort],VCroot_d[cohort],
                                            VCstem_kmax[cohort], VCstem_c[cohort],VCstem_d[cohort], 
                                            VCleaf_kmax[cohort], VCleaf_c[cohort],VCleaf_d[cohort])
    if(relative) resmat[j,] = 100*rrow/sum(rrow)
    else resmat[j,] = rrow
  }
  if(draw) {
    cols = c("black", "red", "green", "blue")
    resdraw = resmat
    if(cumulative) {
      resdraw[,2] = resdraw[,1] + resdraw[,2]
      resdraw[,3] = resdraw[,2] + resdraw[,3]
      resdraw[,4] = resdraw[,3] + resdraw[,4]
    }
    dates = as.Date(rownames(x$WaterBalance))
    plotAxes<-function(){
      if(!yearAxis) axis.Date(1, dates)
      else {
        axis(1, at = (0:numYears)*365, labels=FALSE)
        axis(1, at = -182+365*(1:numYears), tick = FALSE, line=FALSE, labels=firstYear:(firstYear+numYears-1))
      }
      axis(2)    
    }
    if(is.null(xlab)) xlab = ifelse(yearAxis,"Year", "Date")  
    if(is.null(ylab)) ylab = ifelse(relative, "Relative resistances (%)", "Resistances")
    ylim = c(0, max(resdraw))
    plot(dates, ylim = ylim, resdraw[,1], type="n", ylab=ylab, 
         xlab=xlab, frame=FALSE, axes=FALSE)
    if(!cumulative) {
      lines(dates, resdraw[,1], lty=1, lwd=1.5, col = cols[1])
      lines(dates, resdraw[,2], lty=2, lwd=1.5, col = cols[2])
      lines(dates, resdraw[,3], lty=3, lwd=1.5, col = cols[3])
      lines(dates, resdraw[,4], lty=4, lwd=1.5, col = cols[4])
      legend("topright", legend = colnames(resdraw), lty=1:4, col = cols, lwd = 1.5, bty="n")
    } else {
      polygon(c(dates[1], dates, dates[length(dates)]),
              c(0, resdraw[,4], 0),
              col = cols[4])
      polygon(c(dates[1], dates, dates[length(dates)]),
              c(0, resdraw[,3], 0),
              col = cols[3])
      polygon(c(dates[1], dates, dates[length(dates)]),
              c(0, resdraw[,2], 0),
              col = cols[2])
      polygon(c(dates[1], dates, dates[length(dates)]),
              c(0, resdraw[,1], 0),
              col = cols[1])
    }
    plotAxes()
    
  }
  if(draw) invisible(resmat)
  else return(resmat)
}