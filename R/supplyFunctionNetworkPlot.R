hydraulics.supplyFunctionNetworkPlot<-function(soil, x, psiCav=0) {
  
  psic = soil$psi
  VG_nc = soil$VG_n
  VG_alphac = soil$VG_alpha
  
  VCroot_kmax = x$below$VCroot_kmax
  VGrhizo_kmax = x$below$VGrhizo_kmax
  
  VCroot_c = x$paramsTransp$VCroot_c
  VCroot_d = x$paramsTransp$VCroot_d
  VCstem_kmax = x$paramsTransp$VCstem_kmax
  VCstem_c = x$paramsTransp$VCstem_c
  VCstem_d = x$paramsTransp$VCstem_d
  
  ncoh = nrow(x$above)
  l = vector("list", ncoh)
  for(i in 1:ncoh) {
    l[[i]] = hydraulics.supplyFunctionNetwork(psic,
                                          VGrhizo_kmax[i,],VG_nc,VG_alphac,
                                          VCroot_kmax[i,], VCroot_c[i],VCroot_d[i],
                                          VCstem_kmax[i], VCstem_c[i],VCstem_d[i], psiCav = psiCav)
  }
  maxE = 0
  minPsi = 0
  for(i in 1:ncoh) {
    maxE = max(maxE, max(l[[i]]$E, na.rm=T))
    minPsi = min(minPsi, min(l[[i]]$PsiPlant))
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylab = "Flow rate")
    } else {
      lines(-l[[i]]$PsiPlant, l[[i]]$E, lty=i)
    }
  }
}
