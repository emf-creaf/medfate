#Draws the supply function (E vs PlantPsi) for the current soil state and plant hydraulic parameters
hydraulics_supplyFunctionPlot<-function(x, soil, draw = TRUE, type="E") {
  
  TYPES = c("E","dEdP","psiStem","psiRoot","psiRhizo", "ERhizo", "resistances")
  type = match.arg(type,TYPES)  
  
  psiSoil = soil_psi(soil, model="VG")
  VG_nc = soil$VG_n
  VG_alphac = soil$VG_alpha
  VCroot_kmax = x$below$VCroot_kmax
  VGrhizo_kmax = x$below$VGrhizo_kmax
  PLCstem = x$PLCstem
  nlayer = length(psiSoil)
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
  cohortnames = row.names(x$cohorts)
  ncoh = length(cohortnames)
  l = vector("list", ncoh)
  names(l) = cohortnames
  for(i in 1:ncoh) {
    VGrhizo_kmaxc = VGrhizo_kmax[i,]
    VCroot_kmaxc = VCroot_kmax[i,]
    psic = psiSoil[VGrhizo_kmaxc>0]
    VGrhizo_kmaxc = VGrhizo_kmaxc[VGrhizo_kmaxc>0]
    VCroot_kmaxc = VCroot_kmaxc[VCroot_kmaxc>0]
    l[[i]] = hydraulics_supplyFunctionNetwork(psic,
                                          VGrhizo_kmaxc,VG_nc,VG_alphac,
                                          VCroot_kmaxc, VCroot_c[i],VCroot_d[i],
                                          VCstem_kmax[i], VCstem_c[i],VCstem_d[i], 
                                          VCleaf_kmax[i], VCleaf_c[i],VCleaf_d[i],
                                          PLCstem = PLCstem,
                                          minFlow = 0.0, maxNsteps = numericParams$maxNsteps, 
                                          ntrial = numericParams$ntrial,
                                          psiTol = numericParams$psiTol, ETol = numericParams$ETol)
  }
  if(draw) {
    minPsi = 0
    psi = numeric(0)
    E = numeric(0)
    dEdP = numeric(0)
    psiStem = numeric(0)
    cohort = character(0)
    xlab = "Leaf pressure (-MPa)"
    for(i in 1:ncoh) {
      minPsi = min(minPsi, min(l[[i]]$psiLeaf, na.rm = T))
      psi = c(psi, -l[[i]]$psiLeaf)
      psiStem = c(psiStem, -l[[i]]$psiStem[,1])
      dEdP = c(dEdP, l[[i]]$dEdP)
      E = c(E, l[[i]]$E)
      cohort = c(cohort, rep(cohortnames[i], length(l[[i]]$E)))
    }
    df = data.frame(psi = psi, E = E, dEdP = dEdP, cohort = cohort)
    if(type=="E") {
      ylab = expression(paste("Flow rate    ",(mmolH[2]*O%.%s^{-1}%.%m^{-2})))
      g<-ggplot(df, aes(x = psi, y=E))+
        geom_path(aes(col=cohort, linetype=cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      return(g)
    } 
    else if(type=="dEdP") {
      ylab = expression(paste("dE/dP  ",(mmol*H[2]*O%.%s^{-1}%.%m^{-2}%.%MPa^{-1})))
      g<-ggplot(df, aes(x = psi, y=dEdP))+
        geom_path(aes(col=cohort, linetype=cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      return(g)
    }
    else if(type=="psiStem") {
      ylab = "Stem pressure (-MPa)"
      g<-ggplot(df, aes(x = psi, y=psiStem))+
        geom_path(aes(col=cohort, linetype=cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      return(g)
    }
    else if(type=="psiRoot") {
      minE = 0
      for(i in 1:ncoh) {
        minE = min(minE, min(l[[i]]$psiRoot, na.rm=T))
      }
      for(i in 1:ncoh) {
        if(i==1) {
          plot(-l[[i]]$psiLeaf, -l[[i]]$psiRoot, type="l", ylim=c(0,-minE+0.1), xlim=c(0,-minPsi),
               xlab = "Leaf pressure (-MPa)", ylab = "Root pressure (-MPa)", col=i)
        } else {
          lines(-l[[i]]$psiLeaf, -l[[i]]$psiRoot, lty=i, col=i)
        }
      }
      abline(h=0, col="gray")
      abline(a=0, b=1, col="gray")
      legend("topleft", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
    }
    else if(type=="ERhizo") {
      minE = 0
      maxE = 0
      for(i in 1:ncoh) {
        maxE = max(maxE, max(l[[i]]$ERhizo, na.rm=T))
        minE = min(minE, min(l[[i]]$ERhizo, na.rm=T))
      }
      for(i in 1:ncoh) {
        if(i==1) {
          matplot(-l[[i]]$psiLeaf, l[[i]]$ERhizo, type="l", lty=i, ylim=c(minE-0.1,maxE+0.1), xlim=c(0,-minPsi),
                  xlab = "Leaf pressure (-MPa)", 
                  ylab = expression(paste("Flow rate from/to layers   "(mmolH[2]*O%.%s^{-1}%.%m^{-2}))), 
                  col = col)
        } else {
          matlines(-l[[i]]$psiLeaf, l[[i]]$ERhizo, lty=i, col = col)
        }
      }
      abline(h=0, col="gray")
      legend("topleft", legend = cohortnames, lty=1:ncoh, bty="n")
      legend("left", legend = paste("Layer", 1:nlayer), lty=1, col=col, bty="n")
    }
    else if(type=="psiRhizo") {
      minE = 0
      minPsi = 0
      for(i in 1:ncoh) {
        minE = min(minE, min(l[[i]]$psiRhizo, na.rm=T))
        minPsi = min(minPsi, min(l[[i]]$psiLeaf))
      }
      for(i in 1:ncoh) {
        if(i==1) {
          matplot(-l[[i]]$psiLeaf, -l[[i]]$psiRhizo, type="l", lty=i, ylim=c(0,-minE+0.1), xlim=c(0,-minPsi),
                  xlab = "Leaf pressure (-MPa)", ylab = "Rhizosphere pressure (-MPa)", col = col)
        } else {
          matlines(-l[[i]]$psiLeaf, -l[[i]]$psiRhizo, lty=i, col = col)
        }
      }
      abline(h=0, col="gray")
      legend("topleft", legend = cohortnames, lty=1:ncoh, bty="n")
      legend("topright", legend = paste("Layer", 1:nlayer), lty=1, col=col, bty="n")
    }
    else if(type=="resistances") {
      for(i in 1:ncoh) {
        nsteps = length(l[[i]]$E)
        resmat = matrix(0, nrow=nsteps, ncol = 4)
        for(j in 1:nsteps) {
          rrow  = hydraulics_soilPlantResistances(psiSoil = psic,
                                                  psiRhizo = l[[i]]$psiRhizo[j,],
                                                  psiStem = l[[i]]$psiStem[j,],
                                                  PLCstem = PLCstem,
                                                  psiLeaf = l[[i]]$psiLeaf[j],
                                                  VGrhizo_kmax[i,],VG_nc,VG_alphac,
                                                  VCroot_kmax[i,], VCroot_c[i],VCroot_d[i],
                                                  VCstem_kmax[i], VCstem_c[i],VCstem_d[i], 
                                                  VCleaf_kmax[i], VCleaf_c[i],VCleaf_d[i])
          resmat[j,] = 100*rrow/sum(rrow)
        }
        if(i==1) {
          plot(-l[[i]]$psiLeaf, resmat[,1], type="l", ylim=c(0,100), xlim=c(0,-minPsi),
               xlab = "Leaf pressure (-MPa)", 
               ylab = expression(paste("Percent resistances")), 
               col=i, lty=1)
          lines(-l[[i]]$psiLeaf, resmat[,2], lty=2, col=i)
          lines(-l[[i]]$psiLeaf, resmat[,3], lty=3, col=i)
          lines(-l[[i]]$psiLeaf, resmat[,4], lty=4, col=i)
        } else {
          lines(-l[[i]]$psiLeaf, resmat[,1], lty=1, col=i)
          lines(-l[[i]]$psiLeaf, resmat[,2], lty=2, col=i)
          lines(-l[[i]]$psiLeaf, resmat[,3], lty=3, col=i)
          lines(-l[[i]]$psiLeaf, resmat[,4], lty=4, col=i)
        }
      }
    }
  }
  if(draw) invisible(l)
  else return(l)
}