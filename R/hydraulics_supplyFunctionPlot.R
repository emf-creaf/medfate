#Draws the supply function (E vs PlantPsi) for the current soil state and plant hydraulic parameters
#' @rdname hydraulics_supplyfunctions
#' 
#' @param x An object of class \code{\link{spwbInput}}.
#' @param type Plot type for \code{hydraulics_supplyFunctionPlot}, either \code{"E"}, 
#' \code{"ERhizo"}, \code{"StemPsi"}, \code{"RootPsi"} or \code{"dEdP"}).
#' @param draw A flag to indicate whether the supply function should be drawn or just returned.
#' @param speciesNames A flag to indicate the use of species names instead of cohort names in plots.
#' @param ylim Graphical parameter to override function defaults.
#' 
#' @keywords internal
hydraulics_supplyFunctionPlot<-function(x, draw = TRUE, type="E", speciesNames = FALSE, ylim=NULL) {
  
  TYPES = c("E","dEdP","StemPsi","RootPsi","ERhizo")
  type = match.arg(type,TYPES)  
  
  psiSoil = soil_psi(x$soil, model="VG")
  VG_nc = x$soil$VG_n
  VG_alphac = x$soil$VG_alpha
  VCroot_kmax = x$belowLayers$VCroot_kmax
  VGrhizo_kmax = x$belowLayers$VGrhizo_kmax
  StemPLC = x$internalWater$StemPLC
  LeafPLC = x$internalWater$LeafPLC
  nlayer = length(psiSoil)
  col = rainbow(nlayer, start = 0.8, end = 0.1)
  
  numericParams = x$control$numericParams
  
  VCroot_c = x$paramsTransp$VCroot_c
  VCroot_d = x$paramsTransp$VCroot_d
  VCstem_kmax = x$paramsTransp$VCstem_kmax
  VCstem_c = x$paramsTransp$VCstem_c
  VCstem_d = x$paramsTransp$VCstem_d
  VCleaf_kmax = x$paramsTransp$VCleaf_kmax
  VCleaf_c = x$paramsTransp$VCleaf_c
  VCleaf_d = x$paramsTransp$VCleaf_d
  cohortnames = row.names(x$cohorts)
  if(speciesNames) cohortnames = as.character(x$cohorts$Name)
  
  ncoh = length(cohortnames)
  l = vector("list", ncoh)
  names(l) = cohortnames
  for(i in 1:ncoh) {
    VGrhizo_kmaxc = VGrhizo_kmax[i,]
    VCroot_kmaxc = VCroot_kmax[i,]
    psic = psiSoil[VGrhizo_kmaxc>0]
    VGrhizo_kmaxc = VGrhizo_kmaxc[VGrhizo_kmaxc>0]
    VCroot_kmaxc = VCroot_kmaxc[VCroot_kmaxc>0]
    hn = list("numericParams" = x$control$numericParams,
              "stemCavitationEffects" = x$control$stemCavitationEffects,
              "leafCavitationEffects" = x$control$leafCavitationEffects,
              "psisoil" = psiSoil,
               "krhizomax" = VGrhizo_kmaxc, "nsoil" = VG_nc, "alphasoil" = VG_alphac,
               "krootmax" = VCroot_kmaxc, "rootc" = VCroot_c[i], "rootd" = VCroot_d[i],
               "kstemmax" = VCstem_kmax[i], "stemc" = VCstem_c[i], "stemd" = VCstem_d[i], 
               "kleafmax" = VCleaf_kmax[i], "leafc" = VCleaf_c[i], "leafd" = VCleaf_d[i],
               "PLCstem" = StemPLC[i],
               "PLCleaf" = LeafPLC[i])
    l[[i]] = hydraulics_supplyFunctionNetwork(hn)
  }
  if(draw) {
    minPsi = 0
    psi = numeric(0)
    E = numeric(0)
    dEdP = numeric(0)
    psiStem = numeric(0)
    psiRoot = numeric(0)
    cohort = character(0)
    xlab = "Leaf pressure (-MPa)"
    for(i in 1:ncoh) {
      minPsi = min(minPsi, min(l[[i]]$psiLeaf, na.rm = T))
      psi = c(psi, -l[[i]]$psiLeaf)
      psiStem = c(psiStem, -l[[i]]$psiStem[,1])
      psiRoot = c(psiRoot, -l[[i]]$psiRoot)
      dEdP = c(dEdP, l[[i]]$dEdP)
      E = c(E, l[[i]]$E)
      ci = rep(cohortnames[i], length(l[[i]]$E))
      cohort = c(cohort, ci)
      Eri = as.vector(l[[i]]$ERhizo)
      Psiri = as.vector(-l[[i]]$psiRhizo)
      Li = gl(n=ncol(l[[i]]$ERhizo), k=nrow(l[[i]]$ERhizo), labels = paste0("Layer ", 1:ncol(l[[i]]$ERhizo)))
      dfi = data.frame("Psi" = rep(-l[[i]]$psiLeaf, ncol(l[[i]]$ERhizo)), 
                       "PsiRhizo" = Psiri, 
                       "ERhizo" = Eri, 
                       "layer" = Li, "cohort" = ci)         
      if(!exists("dfRhizo")) {
        dfRhizo  = dfi
      } else {
        dfRhizo = rbind(dfRhizo, dfi)
      }
    }
    dfRhizo$cohort = factor(dfRhizo$cohort, levels = cohortnames)
    
    df = data.frame("psi" = psi, 
                    "StemPsi" = psiStem, 
                    "RootPsi" = psiRoot, 
                    "E" = E, "dEdP" = dEdP, 
                    "cohort" = cohort)
    df$cohort = factor(df$cohort, levels = cohortnames)
    if(type=="E") {
      ylab = expression(paste("Flow rate    ",(mmol%.%s^{-1}%.%m^{-2})))
      g<-ggplot(df, aes(x = .data$psi, y=.data$E))+
        geom_path(aes(col=.data$cohort, 
                      linetype=.data$cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    } 
    else if(type=="dEdP") {
      ylab = expression(paste("dE/dP  ",(mmol%.%s^{-1}%.%m^{-2}%.%MPa^{-1})))
      g<-ggplot(df, aes(x = .data$psi, y=.data$dEdP))+
        geom_path(aes(col=.data$cohort, linetype=.data$cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    }
    else if(type=="StemPsi") {
      ylab = "Stem pressure (-MPa)"
      g<-ggplot(df, aes(x = .data$psi, y=.data$StemPsi))+
        geom_path(aes(col=.data$cohort, linetype=.data$cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    }
    else if(type=="RootPsi") {
      ylab = "Root crown pressure (-MPa)"
      g<-ggplot(df, aes(x = .data$psi, y=.data$RootPsi))+
        geom_path(aes(col=.data$cohort, linetype=.data$cohort))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    }
    else if(type=="ERhizo") {
      ylab = expression(paste("Flow rate from/to layers "(mmol%.%s^{-1}%.%m^{-2})))
      g<-ggplot(dfRhizo, aes(x = .data$Psi, y=.data$ERhizo))+
        geom_path(aes(col=.data$cohort, linetype=.data$layer))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")+
        geom_hline(yintercept=0, col="gray")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    }
    else if(type=="PsiRhizo") {
      ylab = "Rhizosphere pressure (-MPa)"
      g<-ggplot(dfRhizo, aes(x = .data$Psi, y=.data$PsiRhizo))+
        geom_path(aes(col=.data$cohort, linetype=.data$layer))+
        scale_color_discrete(name="")+
        scale_linetype_discrete(name="")
      g<-g+xlab(xlab)+ylab(ylab)+theme_bw()
      if(!is.null(ylim)) g<-g+ylim(ylim)
      return(g)
    }
    # else if(type=="resistances") {
    #   for(i in 1:ncoh) {
    #     nsteps = length(l[[i]]$E)
    #     resmat = matrix(0, nrow=nsteps, ncol = 4)
    #     for(j in 1:nsteps) {
    #       rrow  = hydraulics_soilPlantResistances(psiSoil = psic,
    #                                               psiRhizo = l[[i]]$psiRhizo[j,],
    #                                               psiStem = l[[i]]$psiStem[j,],
    #                                               PLCstem = PLCstem,
    #                                               psiLeaf = l[[i]]$psiLeaf[j],
    #                                               VGrhizo_kmax[i,],VG_nc,VG_alphac,
    #                                               VCroot_kmax[i,], VCroot_c[i],VCroot_d[i],
    #                                               VCstem_kmax[i], VCstem_c[i],VCstem_d[i], 
    #                                               VCleaf_kmax[i], VCleaf_c[i],VCleaf_d[i])
    #       resmat[j,] = 100*rrow/sum(rrow)
    #     }
    #     if(i==1) {
    #       plot(-l[[i]]$psiLeaf, resmat[,1], type="l", ylim=c(0,100), xlim=c(0,-minPsi),
    #            xlab = "Leaf pressure (-MPa)", 
    #            ylab = expression(paste("Percent resistances")), 
    #            col=i, lty=1)
    #       lines(-l[[i]]$psiLeaf, resmat[,2], lty=2, col=i)
    #       lines(-l[[i]]$psiLeaf, resmat[,3], lty=3, col=i)
    #       lines(-l[[i]]$psiLeaf, resmat[,4], lty=4, col=i)
    #     } else {
    #       lines(-l[[i]]$psiLeaf, resmat[,1], lty=1, col=i)
    #       lines(-l[[i]]$psiLeaf, resmat[,2], lty=2, col=i)
    #       lines(-l[[i]]$psiLeaf, resmat[,3], lty=3, col=i)
    #       lines(-l[[i]]$psiLeaf, resmat[,4], lty=4, col=i)
    #     }
    #   }
    # }
  }
  return(l)
}
