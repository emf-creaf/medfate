spwb_resistances<-function(x, cohort = 1, relative = FALSE, draw = FALSE, 
                           cumulative = FALSE, xlab = NULL, ylab=NULL) {
  
  if(x$spwbInput$control$transpirationMode!="Sperry") {
    stop("Resistances can only be calculated when transpirationMode = 'Sperry'.")
  }
    
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
  
  psiLeaf = x$LeafPsiMin
  psiStem = x$StemPsi
  psiRoot = x$RootPsi
  PLCstem = x$PlantStress
  psiRhizo = x$RhizoPsi
  
  nlayers = length(VG_nc)
  psiSoil = x$Soil$psi.1
  if(nlayers>1) psiSoil = cbind(psiSoil, x$Soil$psi.2)
  if(nlayers>2) psiSoil = cbind(psiSoil, x$Soil$psi.3)
  if(nlayers>3) psiSoil = cbind(psiSoil, x$Soil$psi.4)
  if(nlayers>4) psiSoil = cbind(psiSoil, x$Soil$psi.5)
  
  nsteps = nrow(psiSoil)
  resmat = matrix(0, nrow=nsteps, ncol = 4)
  rownames(resmat) = rownames(psiStem)
  colnames(resmat) = c("Rhizosphere", "Root", "Stem", "Leaf")
  for(j in 1:nsteps) {
    rrow  = hydraulics_soilPlantResistances(psiSoil = psiSoil[j,],
                                            psiRhizo = psiRhizo[[cohort]][j,],
                                            psiStem = psiStem[j,cohort],
                                            PLCstem = PLCstem[j,cohort],
                                            psiLeaf = psiLeaf[j,cohort],
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