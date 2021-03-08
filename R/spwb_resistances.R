spwb_resistances<-function(x, cohort = 1, relative = FALSE, draw = FALSE, 
                           cumulative = FALSE, xlab = NULL, ylab=NULL) {
  
  if(x$spwbInput$control$transpirationMode!="Sperry") {
    stop("Resistances can only be calculated when transpirationMode = 'Sperry'.")
  }
    
  VCroot_kmax = x$spwbInput$belowLayers$VCroot_kmax
  VGrhizo_kmax = x$spwbInput$belowLayers$VGrhizo_kmax
  VG_nc = x$spwbInput$soil$VG_n
  VG_alphac = x$spwbInput$soil$VG_alpha
  
  paramsTranspiration = x$spwbInput$paramsTranspiration
  VCroot_c = paramsTranspiration$VCroot_c
  VCroot_d = paramsTranspiration$VCroot_d
  VCstem_kmax = paramsTranspiration$VCstem_kmax
  VCstem_c = paramsTranspiration$VCstem_c
  VCstem_d = paramsTranspiration$VCstem_d
  VCleaf_kmax = paramsTranspiration$VCleaf_kmax
  VCleaf_c = paramsTranspiration$VCleaf_c
  VCleaf_d = paramsTranspiration$VCleaf_d
  
  LeafPsi = x$Plants$LeafPsiMin
  StemPsi = x$Plants$StemPsi
  RootPsi = x$Plants$RootPsi
  StemPLC = x$Plants$PlantStress
  RhizoPsi = x$Plants$RhizoPsi
  
  nlayers = length(VG_nc)
  psiSoil = x$Soil$psi.1
  if(nlayers>1) psiSoil = cbind(psiSoil, x$Soil$psi.2)
  if(nlayers>2) psiSoil = cbind(psiSoil, x$Soil$psi.3)
  if(nlayers>3) psiSoil = cbind(psiSoil, x$Soil$psi.4)
  if(nlayers>4) psiSoil = cbind(psiSoil, x$Soil$psi.5)
  
  nsteps = nrow(psiSoil)
  resmat = matrix(0, nrow=nsteps, ncol = 4)
  rownames(resmat) = rownames(StemPsi)
  colnames(resmat) = c("Rhizosphere", "Root", "Stem", "Leaf")
  for(j in 1:nsteps) {
    rrow  = hydraulics_soilPlantResistances(psiSoil = psiSoil[j,],
                                            psiRhizo = RhizoPsi[[cohort]][j,],
                                            psiStem = StemPsi[j,cohort],
                                            PLCstem = StemPLC[j,cohort],
                                            psiLeaf = LeafPsi[j,cohort],
                                            VGrhizo_kmax[cohort,],VG_nc,VG_alphac,
                                            VCroot_kmax[cohort,], VCroot_c[cohort],VCroot_d[cohort],
                                            VCstem_kmax[cohort], VCstem_c[cohort],VCstem_d[cohort], 
                                            VCleaf_kmax[cohort], VCleaf_c[cohort],VCleaf_d[cohort])
    if(relative) resmat[j,] = 100*rrow/sum(rrow)
    else resmat[j,] = rrow
  }
  if(draw) {
    if(is.null(ylab)) ylab = ifelse(relative, "Relative resistances (%)", "Resistances")
    if(is.null(xlab)) xlab = ""
    if(!cumulative) {
      g<-.multiple_dynamics(resmat, ylab=ylab, xlab = xlab)
    } else {
      rescum = resmat
      rescum[,2] = rescum[,1] + rescum[,2]
      rescum[,3] = rescum[,2] + rescum[,3]
      rescum[,4] = rescum[,3] + rescum[,4]
      
      df = as.data.frame(rescum)
      df$Date =  as.Date(rownames(x$WaterBalance))
      g<-ggplot(df, aes(x = df$Date))+
        geom_area(aes(y=df$Leaf, fill = "4"))+
        geom_area(aes(y=df$Stem, fill = "3"))+
        geom_area(aes(y=df$Root, fill = "2"))+
        geom_area(aes(y=df$Rhizosphere, fill="1"))+
        scale_fill_manual(name="", values=c("1" = "black", "2"="red",
                                            "3" = "green", "4" = "blue"),
                          labels = c("Rhizosphere", "Root", "Stem", "Leaf"))+
        xlab(xlab)+ylab(ylab)+
        theme_bw()
    }
    return(g)
  } else {
    return(resmat)
  }
}