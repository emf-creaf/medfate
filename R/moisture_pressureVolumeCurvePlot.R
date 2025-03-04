#' @rdname moisture
#' 
#' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
#' @param segment Segment whose relative water content curve to plot, either \code{"stem"} or \code{"leaf"}
#'                (the latter only available if \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}).
#' @param fraction  Tissue fraction, either \code{"symplastic"}, \code{"apoplastic"} or \code{"all"}.
#' @param psiVec Vector of water potential values to evaluate for the pressure-volume curve.
#' @param speciesNames A flag to indicate the use of species names instead of cohort names in plots.
#' 
#' @keywords internal
moisture_pressureVolumeCurvePlot<-function(x, segment="stem", 
                                           fraction = "all",
                                           psiVec =  seq(-0.1, -8.0, by =-0.01),
                                           speciesNames = FALSE) {
  
  if(!inherits(x, "spwbInput") && !inherits(x, "growthInput")) {
    stop("'x' should be of class 'spwbInput' or 'growthInput'")
  }
  transpirationMode <- x$control$transpirationMode
  if(transpirationMode=="Granier") {
    TYPES = c("stem")
  } else {
    TYPES = c("leaf","stem")
  }
  segment <- match.arg(segment,TYPES)  
  fraction <- match.arg(fraction, c("symplastic", "apoplastic", "all"))  
  
  cohortnames <- row.names(x$cohorts)
  if(speciesNames) cohortnames <- as.character(x$cohorts$Name)
  col <- rainbow(2, start = 0.8, end = 0.1)
  
  ncoh = nrow(x$above)

  minPsi = min(psiVec)
  if(segment=="leaf") {
    PVleaf_pi0 <- x$paramsWaterStorage$LeafPI0
    PVleaf_eps <- x$paramsWaterStorage$LeafEPS
    PVleaf_fapo <- x$paramsWaterStorage$LeafAF
    PVleaf_c <- x$paramsTransp$VCleaf_c
    PVleaf_d <- x$paramsTransp$VCleaf_d
    for(i in 1:ncoh) {
      rwcsymleaf = unlist(lapply(psiVec, moisture_symplasticRWC, PVleaf_pi0[i], PVleaf_eps[i]))*100
      rwcapoleaf = unlist(lapply(psiVec, moisture_apoplasticRWC, PVleaf_c[i], PVleaf_d[i]))*100
      if(fraction=="all") {
        rwcleaf = rwcapoleaf*PVleaf_fapo[i] + rwcsymleaf*(1-PVleaf_fapo[i])
      } else if(fraction=="symplastic") {
        rwcleaf = rwcsymleaf
      } else if(fraction=="apoplastic") {
        rwcleaf = rwcapoleaf
      }
      if(i==1) rwc = matrix(rwcleaf, ncol = 1)
      else rwc = cbind(rwc, rwcleaf)
    }
    colnames(rwc) = cohortnames
    .multiple_y(x = -psiVec, y = rwc, xlab = "Leaf pressure (-MPa)", ylab = "Leaf relative water content [%]", ylim =c(0, 100))
  } 
  else if(segment=="stem") {
    PVstem_pi0 <- x$paramsWaterStorage$StemPI0
    PVstem_eps <- x$paramsWaterStorage$StemEPS
    PVstem_fapo <- x$paramsWaterStorage$StemAF
    PVstem_c <- x$paramsTransp$VCstem_c
    PVstem_d <- x$paramsTransp$VCstem_d
    for(i in 1:ncoh) {
      rwcsymstem = unlist(lapply(psiVec, moisture_symplasticRWC, PVstem_pi0[i], PVstem_eps[i]))*100
      rwcapostem = unlist(lapply(psiVec, moisture_apoplasticRWC, PVstem_c[i], PVstem_d[i]))*100
      if(fraction=="all") {
        rwcstem = rwcapostem*PVstem_fapo[i] + rwcsymstem*(1-PVstem_fapo[i])
      } else if(fraction=="symplastic") {
        rwcstem = rwcsymstem
      } else if(fraction=="apoplastic") {
        rwcstem = rwcapostem
      }
      if(i==1) rwc = matrix(rwcstem, ncol = 1)
      else rwc = cbind(rwc, rwcstem)
    }
    colnames(rwc) = cohortnames
    .multiple_y(x = -psiVec, y = rwc, xlab = "Stem pressure (-MPa)", ylab = "Stem relative water content [%]", ylim =c(0, 100))
  } 
}