#Draws the supply function (E vs PlantPsi) for the current soil state and plant hydraulic parameters
hydraulics.supplyFunctionPlot<-function(x, soil, type="E") {
  
  TYPES = c("E","dEdP","PsiStem","PsiRoot","PsiRhizo", "Elayers")
  type = match.arg(type,TYPES)  
  
  psic = soil.psi(soil, model="VG")
  VG_nc = soil$VG_n
  VG_alphac = soil$VG_alpha
  cohortnames = row.names(x$cohorts)
  VCroot_kmax = x$below$VCroot_kmax
  VGrhizo_kmax = x$below$VGrhizo_kmax
  pEmb = x$ProportionCavitated
  nlayer = length(psic)
  col = rainbow(nlayer, start = 0.8, end = 0.1)
  
  numericParams = x$control$numericParams
  
  VCroot_c = x$paramsTransp$VCroot_c
  VCroot_d = x$paramsTransp$VCroot_d
  VCstem_kmax = x$paramsTransp$VCstem_kmax
  VCstem_c = x$paramsTransp$VCstem_c
  VCstem_d = x$paramsTransp$VCstem_d
  VCleaf_kmax = x$paramsTransp$VCleaf_kmax
  VCleaf_c = x$paramsTransp$VCleaf_c
  VCleaf_d = x$paramsTransp$VCleaf_d
  ncoh = nrow(x$above)
  l = vector("list", ncoh)
  names(l) = cohortnames
  for(i in 1:ncoh) {
    psiCav = hydraulics.xylemPsi(1.0-pEmb[i], 1.0, VCstem_c[i], VCstem_d[i])
    l[[i]] = hydraulics.supplyFunctionNetwork(psic,
                                          VGrhizo_kmax[i,],VG_nc,VG_alphac,
                                          VCroot_kmax[i,], VCroot_c[i],VCroot_d[i],
                                          VCstem_kmax[i], VCstem_c[i],VCstem_d[i], 
                                          VCleaf_kmax[i], VCleaf_c[i],VCleaf_d[i],
                                          psiCav = psiCav,
                                          minFlow = 0.0, maxNsteps = numericParams$maxNsteps, psiStep = numericParams$psiStep, 
                                          psiMax = numericParams$psiMax, ntrial = numericParams$ntrial,
                                          psiTol = numericParams$psiTol, ETol = numericParams$ETol)
  }
  #Find minimum psi
  minPsi = 0
  for(i in 1:ncoh) {
    minPsi = min(minPsi, min(l[[i]]$PsiLeaf, na.rm = T))
  }
  minPsi = max(minPsi, -40)
  if(type=="E") {
    maxE = 0
    for(i in 1:ncoh) {
      maxE = max(maxE, max(l[[i]]$E, na.rm=T))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        plot(-l[[i]]$PsiLeaf, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
             xlab = "Leaf pressure (-MPa)", 
             ylab = expression(paste("Flow rate    ",(mmolH[2]*O%.%s^{-1}%.%m^{-2}))), 
             col=i)
      } else {
        lines(-l[[i]]$PsiLeaf, l[[i]]$E, lty=i, col=i)
      }
    }
    legend("topleft", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  } 
  else if(type=="dEdP") {
    maxdEdP = 0
    for(i in 1:ncoh) {
      maxdEdP = max(maxdEdP, max(l[[i]]$dEdP, na.rm=T))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        plot(-l[[i]]$PsiLeaf, l[[i]]$dEdP, type="l", ylim=c(0,maxdEdP+0.1), xlim=c(0,-minPsi),
             xlab = "Leaf pressure (-MPa)", 
             ylab = expression(paste("dE/dP  ",(mmol*H[2]*O%.%s^{-1}%.%m^{-2}%.%MPa^{-1}))), 
             col=i)
      } else {
        lines(-l[[i]]$PsiLeaf, l[[i]]$dEdP, lty=i, col=i)
      }
    }
    legend("topright", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  }
  else if(type=="PsiStem") {
    minE = 0
    for(i in 1:ncoh) {
      minE = min(minE, min(l[[i]]$PsiStem, na.rm=T))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        plot(-l[[i]]$PsiLeaf, -l[[i]]$PsiStem, type="l", ylim=c(0,-minE+0.1), xlim=c(0,-minPsi),
             xlab = "Leaf pressure (-MPa)", ylab = "Stem pressure (-MPa)", col=i)
      } else {
        lines(-l[[i]]$PsiLeaf, -l[[i]]$PsiStem, lty=i, col=i)
      }
    }
    abline(h=0, col="gray")
    abline(a=0, b=1, col="gray")
    legend("topleft", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  }
  else if(type=="PsiRoot") {
    minE = 0
    for(i in 1:ncoh) {
      minE = min(minE, min(l[[i]]$PsiRoot, na.rm=T))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        plot(-l[[i]]$PsiLeaf, -l[[i]]$PsiRoot, type="l", ylim=c(0,-minE+0.1), xlim=c(0,-minPsi),
             xlab = "Leaf pressure (-MPa)", ylab = "Root pressure (-MPa)", col=i)
      } else {
        lines(-l[[i]]$PsiLeaf, -l[[i]]$PsiRoot, lty=i, col=i)
      }
    }
    abline(h=0, col="gray")
    abline(a=0, b=1, col="gray")
    legend("topleft", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  }
  else if(type=="Elayers") {
    minE = 0
    maxE = 0
    for(i in 1:ncoh) {
      maxE = max(maxE, max(l[[i]]$Elayers, na.rm=T))
      minE = min(minE, min(l[[i]]$Elayers, na.rm=T))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        matplot(-l[[i]]$PsiLeaf, l[[i]]$Elayers, type="l", lty=i, ylim=c(minE-0.1,maxE+0.1), xlim=c(0,-minPsi),
             xlab = "Leaf pressure (-MPa)", 
             ylab = expression(paste("Flow rate from/to layers   "(mmolH[2]*O%.%s^{-1}%.%m^{-2}))), 
             col = col)
      } else {
        matlines(-l[[i]]$PsiLeaf, l[[i]]$Elayers, lty=i, col = col)
      }
    }
    abline(h=0, col="gray")
    legend("topleft", legend = cohortnames, lty=1:ncoh, bty="n")
    legend("left", legend = paste("Layer", 1:nlayer), lty=1, col=col, bty="n")
  }
  else if(type=="PsiRhizo") {
    minE = 0
    minPsi = 0
    for(i in 1:ncoh) {
      minE = min(minE, min(l[[i]]$PsiRhizo, na.rm=T))
      minPsi = min(minPsi, min(l[[i]]$PsiLeaf))
    }
    for(i in 1:ncoh) {
      if(i==1) {
        matplot(-l[[i]]$PsiLeaf, -l[[i]]$PsiRhizo, type="l", lty=i, ylim=c(0,-minE+0.1), xlim=c(0,-minPsi),
                xlab = "Leaf pressure (-MPa)", ylab = "Rhizosphere pressure (-MPa)", col = col)
      } else {
        matlines(-l[[i]]$PsiLeaf, -l[[i]]$PsiRhizo, lty=i, col = col)
      }
    }
    abline(h=0, col="gray")
    legend("topleft", legend = cohortnames, lty=1:ncoh, bty="n")
    legend("topright", legend = paste("Layer", 1:nlayer), lty=1, col=col, bty="n")
  }
  invisible(l)
}