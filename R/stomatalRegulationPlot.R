transp.stomatalRegulationPlot<-function(x, soil, meteo, day, timestep, latitude, elevation) {
  
  dctr = transp.stomatalRegulation(x, soil, meteo, day, latitude, elevation)
  ncoh = length(dctr)
  oldpar = par(mar=c(5,5,1,1), mfrow = c(5,2))
  
  l = vector("list", ncoh)
  phsunlit = vector("list", ncoh)
  phshade = vector("list", ncoh)
  pmsunlit = vector("list", ncoh)
  pmshade = vector("list", ncoh)
  
  for(i in 1:ncoh) {
    l[[i]] = dctr[[i]]$supply
    phsunlit[[i]] = dctr[[i]]$photoSunlit[[timestep]]
    phshade[[i]] = dctr[[i]]$photoShade[[timestep]]
    pmsunlit[[i]] = dctr[[i]]$PMSunlit[[timestep]]
    pmshade[[i]] = dctr[[i]]$PMShade[[timestep]]
  }
  maxE = 0
  maxA = 0
  minPsi = 0
  minTemp = 1000
  maxTemp = -1000
  minVPD = 1000
  maxVPD = -1000
  minGw = 1000
  maxGw = -1000
  for(i in 1:ncoh) {
    maxE = max(maxE, max(l[[i]]$E, na.rm=T))
    if(sum(!is.na(phsunlit[[i]]$Photosynthesis))>0) maxA = max(maxA, max(phsunlit[[i]]$Photosynthesis, na.rm=T))
    if(sum(!is.na(phshade[[i]]$Photosynthesis))>0) maxA = max(maxA, max(phshade[[i]]$Photosynthesis, na.rm=T))
    minPsi = min(minPsi, min(l[[i]]$PsiLeaf))
    minTemp = min(minTemp, min(phsunlit[[i]]$LeafTemperature))
    minTemp = min(minTemp, min(phshade[[i]]$LeafTemperature))
    maxTemp = max(maxTemp, max(phsunlit[[i]]$LeafTemperature))
    maxTemp = max(maxTemp, max(phshade[[i]]$LeafTemperature))
    minVPD = min(minVPD, min(phsunlit[[i]]$LeafVPD))
    minVPD = min(minVPD, min(phshade[[i]]$LeafVPD))
    maxVPD = max(maxVPD, max(phsunlit[[i]]$LeafVPD))
    maxVPD = max(maxVPD, max(phshade[[i]]$LeafVPD))
    minGw = min(minGw, min(phsunlit[[i]]$WaterVaporConductance))
    minGw = min(minGw, min(phshade[[i]]$WaterVaporConductance))
    maxGw = max(maxGw, max(phsunlit[[i]]$WaterVaporConductance))
    maxGw = max(maxGw, max(phshade[[i]]$WaterVaporConductance))
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Flow rate sunlit "(mmol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, l[[i]]$E, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  legend("topleft", legend = names(dctr), lty=1:length(dctr), bty = "n")
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab =expression(paste("Flow rate shade "(mmol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, l[[i]]$E, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phsunlit[[i]]$Photosynthesis, type="l", ylim=c(0,maxA+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Photosynthesis sunlit  "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, phsunlit[[i]]$Photosynthesis, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phshade[[i]]$Photosynthesis, type="l", ylim=c(0,maxA+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Photosynthesis shade  "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, phshade[[i]]$Photosynthesis, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phsunlit[[i]]$WaterVaporConductance, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minGw, maxGw), 
           ylab = expression(paste("Leaf sunlit stomatal conductance "(mol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, phsunlit[[i]]$WaterVaporConductance, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phshade[[i]]$WaterVaporConductance, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minGw, maxGw), 
           ylab = expression(paste("Leaf shade stomatal conductance "(mol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$PsiLeaf, phshade[[i]]$WaterVaporConductance, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phsunlit[[i]]$LeafTemperature, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minTemp, maxTemp), 
           ylab = "Leaf sunlit temperature (degrees C)")
    } else {
      lines(-l[[i]]$PsiLeaf, phsunlit[[i]]$LeafTemperature, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phshade[[i]]$LeafTemperature, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minTemp, maxTemp), 
           ylab = "Leaf shade temperature (degrees C)")
    } else {
      lines(-l[[i]]$PsiLeaf, phshade[[i]]$LeafTemperature, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phsunlit[[i]]$LeafVPD, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minVPD, maxVPD), 
           ylab = "Leaf sunlit VPD (kPa)")
    } else {
      lines(-l[[i]]$PsiLeaf, phsunlit[[i]]$LeafVPD, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$PsiLeaf, phshade[[i]]$LeafVPD, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minVPD, maxVPD), 
           ylab = "Leaf shade VPD (kPa)")
    } else {
      lines(-l[[i]]$PsiLeaf, phshade[[i]]$LeafVPD, lty=i)
    }
    abline(v = -l[[i]]$PsiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  par(oldpar)
}