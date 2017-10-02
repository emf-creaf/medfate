## Plots day transpiration
transp.dayCanopyTranspirationPlot<-function(x, soil, meteo, day, timestep, latitude, elevation, slope=0, aspect=0) {
  
  dctr = transp.dayCanopyTranspiration(x, soil, meteo, day, latitude, elevation, slope, aspect)
  ncoh = length(dctr)
  
  par(mfrow = c(1,2))
  
  l = vector("list", ncoh)
  p = vector("list", ncoh)
  pm = vector("list", ncoh)
  for(i in 1:ncoh) {
    l[[i]] = dctr[[i]]$supply
    p[[i]] = dctr[[i]]$photo[[timestep]]
    pm[[i]] = dctr[[i]]$PM[[timestep]]
  }
  maxE = 0
  maxA = 0
  minPsi = 0
  for(i in 1:ncoh) {
    maxE = max(maxE, max(l[[i]]$E, na.rm=T))
    if(sum(!is.na(p[[i]]$Photosynthesis))>0) maxA = max(maxA, max(p[[i]]$Photosynthesis, na.rm=T))
    minPsi = min(minPsi, min(l[[i]]$PsiPlant))
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylab = "Flow rate")
    } else {
      lines(-l[[i]]$PsiPlant, l[[i]]$E, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  legend("topleft", legend = names(dctr), lty=1:length(dctr), bty = "n")
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, p[[i]]$Photosynthesis, type="l", ylim=c(0,maxA+0.1), xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylab = "Photosynthesis")
    } else {
      lines(-l[[i]]$PsiPlant, p[[i]]$Photosynthesis, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
}