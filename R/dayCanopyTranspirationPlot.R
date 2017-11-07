## Plots day transpiration
transp.dayCanopyTranspirationPlot<-function(x, soil, meteo, day, timestep, latitude, elevation, slope=0, aspect=0) {
  
  dctr = transp.dayCanopyTranspiration(x, soil, meteo, day, latitude, elevation, slope, aspect)
  ncoh = length(dctr)
  
  oldpar = par(mfrow = c(3,2))
  
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
  minTemp = 1000
  maxTemp = -1000
  minVPD = 1000
  maxVPD = -1000
  for(i in 1:ncoh) {
    maxE = max(maxE, max(l[[i]]$E, na.rm=T))
    if(sum(!is.na(p[[i]]$Photosynthesis))>0) maxA = max(maxA, max(p[[i]]$Photosynthesis, na.rm=T))
    minPsi = min(minPsi, min(l[[i]]$PsiPlant))
    minTemp = min(minTemp, min(p[[i]]$LeafTempSL))
    minTemp = min(minTemp, min(p[[i]]$LeafTempSH))
    maxTemp = max(maxTemp, max(p[[i]]$LeafTempSL))
    maxTemp = max(maxTemp, max(p[[i]]$LeafTempSH))
    minVPD = min(minVPD, min(p[[i]]$LeafVPDSL))
    minVPD = min(minVPD, min(p[[i]]$LeafVPDSH))
    maxVPD = max(maxVPD, max(p[[i]]$LeafVPDSL))
    maxVPD = max(maxVPD, max(p[[i]]$LeafVPDSH))
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
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, p[[i]]$LeafTempSL, type="l", xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylim=c(minTemp, maxTemp), ylab = "Leaf temperature (sunlit)")
    } else {
      lines(-l[[i]]$PsiPlant, p[[i]]$LeafTempSL, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, p[[i]]$LeafTempSH, type="l", xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylim=c(minTemp, maxTemp), ylab = "Leaf temperature (shade)")
    } else {
      lines(-l[[i]]$PsiPlant, p[[i]]$LeafTempSH, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, p[[i]]$LeafVPDSL, type="l", xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylim=c(minVPD, maxVPD), ylab = "Leaf VPD (sunlit)")
    } else {
      lines(-l[[i]]$PsiPlant, p[[i]]$LeafVPDSL, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiPlant, p[[i]]$LeafVPDSH, type="l", xlim=c(0,-minPsi),
           xlab = "Plant pressure (-MPa)", ylim=c(minVPD, maxVPD), ylab = "Leaf VPD (shade)")
    } else {
      lines(-l[[i]]$PsiPlant, p[[i]]$LeafVPDSH, lty=i)
    }
    abline(v = -l[[i]]$PsiPlant[pm[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  par(oldpar)
}