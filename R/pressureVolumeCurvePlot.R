moisture.pressureVolumeCurvePlot<-function(x, segment="leaf", psiVec =  seq(-0.1, -8.0, by =-0.01)) {
  
  TYPES = c("leaf","stem", "both")
  type = match.arg(segment,TYPES)  
  
  cohortnames = row.names(x$cohorts)
  col = rainbow(2, start = 0.8, end = 0.1)
  
  PVleaf_pi0 = x$paramsWaterStorage$LeafPI0
  PVleaf_eps = x$paramsWaterStorage$LeafEPS
  PVstem_pi0 = x$paramsWaterStorage$StemPI0
  PVstem_eps = x$paramsWaterStorage$StemEPS
  ncoh = nrow(x$above)

  minPsi = min(psiVec)
  if(segment=="leaf") {
    for(i in 1:ncoh) {
      rwcleaf = unlist(lapply(psiVec, moisture.symplasticRWC, PVleaf_pi0[i], PVleaf_eps[i]))*100
      if(i==1) {
        plot(-psiVec, rwcleaf, type="l", xlim=c(0,-minPsi), ylim =c(0, 100),
             xlab = "Leaf pressure (-MPa)", ylab = paste("Leaf relative water content [%]"), col=i)
      } else {
        lines(-psiVec, rwcleaf, lty=i, col=i)
      }
    }
    legend("topright", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  } 
  else if(segment=="stem") {
    for(i in 1:ncoh) {
      rwcstem = unlist(lapply(psiVec, moisture.symplasticRWC, PVstem_pi0[i], PVstem_eps[i]))*100
      if(i==1) {
        plot(-psiVec, rwcstem, type="l", xlim=c(0,-minPsi), ylim =c(0, 100),
             xlab = "Stem pressure (-MPa)", ylab = paste("Stem relative water content [%]"), col=i)
      } else {
        lines(-psiVec, rwcstem, lty=i, col=i)
      }
    }
    legend("topright", legend = cohortnames, lty=1:ncoh, col = 1:ncoh, bty="n")
  } 
  else if(segment=="both") {
    for(i in 1:ncoh) {
      rwcstem = unlist(lapply(psiVec, moisture.symplasticRWC, PVstem_pi0[i], PVstem_eps[i]))*100
      rwcleaf = unlist(lapply(psiVec, moisture.symplasticRWC, PVleaf_pi0[i], PVleaf_eps[i]))*100
      if(i==1) {
        plot(-psiVec, rwcstem, type="l", xlim=c(0,-minPsi), ylim =c(0, 100),
             xlab = "Stem pressure (-MPa)", ylab = paste("Relative water content [%]"), col=1)
        lines(-psiVec, rwcleaf, lty=i, col=2)
      } else {
        lines(-psiVec, rwcstem, lty=i, col=1)
        lines(-psiVec, rwcleaf, lty=i, col=2)
      }
    }
    legend("topright", legend = cohortnames, lty=1:ncoh, col = 1, bty="n")
    legend("bottomleft", legend = c("stem", "leaf"), lty=1, col = 1:2, bty="n")
  } 
}