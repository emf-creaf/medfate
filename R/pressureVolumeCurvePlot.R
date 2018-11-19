moisture.pressureVolumeCurvePlot<-function(x, segment="leaf", 
                                           fraction = "all",
                                           psiVec =  seq(-0.1, -8.0, by =-0.01)) {
  
  TYPES = c("leaf","stem", "both")
  type = match.arg(segment,TYPES)  
  
  cohortnames = row.names(x$cohorts)
  col = rainbow(2, start = 0.8, end = 0.1)
  
  PVleaf_pi0 = x$paramsWaterStorage$LeafPI0
  PVleaf_eps = x$paramsWaterStorage$LeafEPS
  PVleaf_fapo = x$paramsWaterStorage$LeafAF
  PVleaf_c = x$paramsTransp$VCleaf_c
  PVleaf_d = x$paramsTransp$VCleaf_d
  PVstem_pi0 = x$paramsWaterStorage$StemPI0
  PVstem_eps = x$paramsWaterStorage$StemEPS
  PVstem_fapo = x$paramsWaterStorage$StemAF
  PVstem_c = x$paramsTransp$VCstem_c
  PVstem_d = x$paramsTransp$VCstem_d
  ncoh = nrow(x$above)

  minPsi = min(psiVec)
  if(segment=="leaf") {
    for(i in 1:ncoh) {
      rwcsymleaf = unlist(lapply(psiVec, moisture.symplasticRWC, PVleaf_pi0[i], PVleaf_eps[i]))*100
      rwcapoleaf = unlist(lapply(psiVec, moisture.apoplasticRWC, PVleaf_c[i], PVleaf_d[i]))*100
      if(fraction=="all") {
        rwcleaf = rwcapoleaf*PVleaf_fapo[i] + rwcsymleaf*(1-PVleaf_fapo[i])
      } else if(fraction=="symplastic") {
        rwcleaf = rwcsymleaf
      } else if(fraction=="apoplastic") {
        rwcleaf = rwcapoleaf
      }
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
      rwcsymstem = unlist(lapply(psiVec, moisture.symplasticRWC, PVstem_pi0[i], PVstem_eps[i]))*100
      rwcapostem = unlist(lapply(psiVec, moisture.apoplasticRWC, PVstem_c[i], PVstem_d[i]))*100
      if(fraction=="all") {
        rwcstem = rwcapostem*PVstem_fapo[i] + rwcsymstem*(1-PVstem_fapo[i])
      } else if(fraction=="symplastic") {
        rwcstem = rwcsymstem
      } else if(fraction=="apoplastic") {
        rwcstem = rwcapostem
      }
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
      rwcsymleaf = unlist(lapply(psiVec, moisture.symplasticRWC, PVleaf_pi0[i], PVleaf_eps[i]))*100
      rwcapoleaf = unlist(lapply(psiVec, moisture.apoplasticRWC, PVleaf_c[i], PVleaf_d[i]))*100
      rwcsymstem = unlist(lapply(psiVec, moisture.symplasticRWC, PVstem_pi0[i], PVstem_eps[i]))*100
      rwcapostem = unlist(lapply(psiVec, moisture.apoplasticRWC, PVstem_c[i], PVstem_d[i]))*100
      if(fraction=="all") {
        rwcleaf = rwcapoleaf*PVleaf_fapo[i] + rwcsymleaf*(1-PVleaf_fapo[i])
        rwcstem = rwcapostem*PVstem_fapo[i] + rwcsymstem*(1-PVstem_fapo[i])
      } else if(fraction=="symplastic") {
        rwcleaf = rwcsymleaf
        rwcstem = rwcsymstem
      } else if(fraction=="apoplastic") {
        rwcleaf = rwcapoleaf
        rwcstem = rwcapostem
      }
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