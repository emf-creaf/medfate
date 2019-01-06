spwb.resistances<-function(x) {
  
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
  
  nlayers = length(VG_nc)
  psiSoil = x$Soil$psi.1
  if(nlayers>1) psiSoil = cbind(psiSoil, x$Soil$psi.2)
  if(nlayers>2) psiSoil = cbind(psiSoil, x$Soil$psi.3)
  if(nlayers>3) psiSoil = cbind(psiSoil, x$Soil$psi.4)
  if(nlayers>4) psiSoil = cbind(psiSoil, x$Soil$psi.5)
  
  for(i in 1:ncoh) {
    nsteps = nrow(psiSoil)
    resmat = matrix(0, nrow=nsteps, ncol = 4)
    for(j in 1:nsteps) {
      rrow  = hydraulics.soilPlantResistances(psiSoil = psiSoil[j,],
                                              psiRhizo = psiRhizo[j,],
                                              psiStem = psiStem[j,],
                                              PLCstem = PLCstem[j,],
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
  
  print(head(psiSoil))
}