transp_stomatalRegulationPlot<-function(x, soil, meteo, day, timestep, latitude, elevation) {
  
  dctr = transp_transpirationSperry(x, soil, meteo, day, latitude, elevation, stepFunctions = timestep, 
                                    modifyInput = FALSE)
  ncoh = length(dctr$SupplyFunctions)
  oldpar = par(mar=c(5,5,1,1), mfrow = c(5,2))
  
  l = dctr$SupplyFunctions
  phsunlit = dctr$PhotoSunlitFunctions
  phshade = dctr$PhotoShadeFunctions
  pmsunlit = dctr$PMSunlitFunctions
  pmshade = dctr$PMShadeFunctions
  
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
    minPsi = min(minPsi, min(l[[i]]$psiLeaf))
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
      plot(-l[[i]]$psiLeaf, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Flow rate sunlit "(mmol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, l[[i]]$E, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  legend("topleft", legend = names(l), lty=1:ncoh, bty = "n")
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, l[[i]]$E, type="l", ylim=c(0,maxE+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab =expression(paste("Flow rate shade "(mmol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, l[[i]]$E, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phsunlit[[i]]$Photosynthesis, type="l", ylim=c(0,maxA+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Photosynthesis sunlit  "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, phsunlit[[i]]$Photosynthesis, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phshade[[i]]$Photosynthesis, type="l", ylim=c(0,maxA+0.1), xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", 
           ylab = expression(paste("Photosynthesis shade  "(mu*mol*C%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, phshade[[i]]$Photosynthesis, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phsunlit[[i]]$WaterVaporConductance, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minGw, maxGw), 
           ylab = expression(paste("Leaf sunlit stomatal conductance "(mol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, phsunlit[[i]]$WaterVaporConductance, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phshade[[i]]$WaterVaporConductance, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minGw, maxGw), 
           ylab = expression(paste("Leaf shade stomatal conductance "(mol*H[2]*O%.%s^{-1}%.%m^{-2}))))
    } else {
      lines(-l[[i]]$psiLeaf, phshade[[i]]$WaterVaporConductance, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phsunlit[[i]]$LeafTemperature, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minTemp, maxTemp), 
           ylab = "Leaf sunlit temperature (degrees C)")
    } else {
      lines(-l[[i]]$psiLeaf, phsunlit[[i]]$LeafTemperature, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phshade[[i]]$LeafTemperature, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minTemp, maxTemp), 
           ylab = "Leaf shade temperature (degrees C)")
    } else {
      lines(-l[[i]]$psiLeaf, phshade[[i]]$LeafTemperature, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phsunlit[[i]]$LeafVPD, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minVPD, maxVPD), 
           ylab = "Leaf sunlit VPD (kPa)")
    } else {
      lines(-l[[i]]$psiLeaf, phsunlit[[i]]$LeafVPD, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmsunlit[[i]]$iMaxProfit+1], col="orange", lty=i, lwd=2)
  }
  for(i in 1:ncoh) {
    if(i==1) {
      plot(-l[[i]]$psiLeaf, phshade[[i]]$LeafVPD, type="l", xlim=c(0,-minPsi),
           xlab = "Leaf pressure (-MPa)", ylim=c(minVPD, maxVPD), 
           ylab = "Leaf shade VPD (kPa)")
    } else {
      lines(-l[[i]]$psiLeaf, phshade[[i]]$LeafVPD, lty=i)
    }
    abline(v = -l[[i]]$psiLeaf[pmshade[[i]]$iMaxProfit+1], col="gray", lty=i, lwd=2)
  }
  par(oldpar)
}