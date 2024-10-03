.resistances_sim<-function(x, cohort, relative = FALSE) {
  if(inherits(x, c("spwb","pwb"))) {
    input <- x$spwbInput  
  } else {
    input <- x$growthInput
  }
  if(!(input$control$transpirationMode %in% c("Sperry", "Sureau"))) {
    stop("Resistances can only be calculated when transpirationMode = 'Sperry' or 'Sureau'.")
  }
  cn <- row.names(input$cohorts)
  if(!(cohort %in% cn)) stop("'cohort' must be a string identifying a cohort name")
  i_coh <- which(cn==cohort)
  
  VCroot_kmax <- input$belowLayers$VCroot_kmax
  VGrhizo_kmax <- input$belowLayers$VGrhizo_kmax
  VG_nc <- input$soil$VG_n
  VG_alphac <- input$soil$VG_alpha
  
  paramsTranspiration <- input$paramsTranspiration
  VCstem_kmax <- paramsTranspiration$VCstem_kmax
  VCleaf_kmax <- paramsTranspiration$VCleaf_kmax
  if(input$control$transpirationMode=="Sperry") {
    VCroot_c <- paramsTranspiration$VCroot_c
    VCroot_d <- paramsTranspiration$VCroot_d
    VCstem_c <- paramsTranspiration$VCstem_c
    VCstem_d <- paramsTranspiration$VCstem_d
    VCleaf_c <- paramsTranspiration$VCleaf_c
    VCleaf_d <- paramsTranspiration$VCleaf_d
  } else {
    VCroot_P50 <- paramsTranspiration$VCroot_P50
    VCroot_slope <- paramsTranspiration$VCroot_slope
    VCstem_P50 <- paramsTranspiration$VCstem_P50
    VCstem_slope <- paramsTranspiration$VCstem_slope
    VCleaf_P50 <- paramsTranspiration$VCleaf_P50
    VCleaf_slope <- paramsTranspiration$VCleaf_slope
  }
  
  LeafPsi <- x$Plants$LeafPsiMin
  StemPsi <- x$Plants$StemPsi
  RootPsi <- x$Plants$RootPsi
  LeafPLC <- x$Plants$LeafPLC
  StemPLC <- x$Plants$StemPLC
  RhizoPsi <- x$Plants$RhizoPsi
  
  nlayers <- length(VG_nc)
  psiSoil <- x$Soil$Psi[,1:nlayers]

  nsteps <- nrow(psiSoil)
  pathmat <- matrix(0, nrow=nsteps, ncol = 4)
  rownames(pathmat) <- rownames(StemPsi)
  colnames(pathmat) <- c("Rhizosphere", "Root", "Stem", "Leaf")
  rhizomat <- matrix(0, nrow=nsteps, ncol = nlayers)
  rownames(rhizomat) <- rownames(StemPsi)
  colnames(rhizomat) <- 1:nlayers
  rootmat <- rhizomat
  for(j in 1:nsteps) {
    if(input$control$transpirationMode=="Sperry") {
      rrow <- hydraulics_soilPlantResistancesWeibull(psiSoil = psiSoil[j,],
                                                     psiRhizo = RhizoPsi[[i_coh]][j,],
                                                     psiStem = StemPsi[j,i_coh],
                                                     PLCstem = StemPLC[j,i_coh],
                                                     psiLeaf = LeafPsi[j,i_coh],
                                                     PLCleaf = LeafPLC[j,i_coh],
                                                     VGrhizo_kmax[i_coh,],VG_nc,VG_alphac,
                                                     VCroot_kmax[i_coh,], VCroot_c[i_coh],VCroot_d[i_coh],
                                                     VCstem_kmax[i_coh], VCstem_c[i_coh],VCstem_d[i_coh], 
                                                     VCleaf_kmax[i_coh], VCleaf_c[i_coh],VCleaf_d[i_coh])
    } else {
      rrow <- hydraulics_soilPlantResistancesSigmoid(psiSoil = psiSoil[j,],
                                                     psiRhizo = RhizoPsi[[i_coh]][j,],
                                                     psiStem = StemPsi[j,i_coh],
                                                     PLCstem = StemPLC[j,i_coh],
                                                     psiLeaf = LeafPsi[j,i_coh],
                                                     PLCleaf = LeafPLC[j,i_coh],
                                                     VGrhizo_kmax[i_coh,],VG_nc,VG_alphac,
                                                     VCroot_kmax[i_coh,], VCroot_P50[i_coh],VCroot_slope[i_coh],
                                                     VCstem_kmax[i_coh], VCstem_P50[i_coh],VCstem_slope[i_coh], 
                                                     VCleaf_kmax[i_coh], VCleaf_P50[i_coh],VCleaf_slope[i_coh])
    }
    rhizo_row <- rrow[["rhizosphere"]]
    root_row <- rrow[["root"]]
    stem_row <- rrow[["stem"]]
    leaf_row <- rrow[["leaf"]]
    pathmat[j,] <- c(1/sum(1/rhizo_row), 1/sum(1/root_row), stem_row, leaf_row)
    rootmat[j, ] <- root_row
    rhizomat[j, ] <- rhizo_row
    if(relative) {
      pathmat[j,] <- 100*pathmat[j,]/sum(pathmat[j,])
      rootmat[j,] <- 100*rootmat[j,]/sum(rootmat[j,])
      rhizomat[j,] <- 100*rhizomat[j,]/sum(rhizomat[j,])
    }
  }
  return(list("pathway" = pathmat,
              "root" = rootmat,
              "rhizosphere" = rhizomat))
}

#' Soil-plant resistances
#'
#' Calculates and draws rhizosphere, root, stem and leaf resistances for simulation time steps
#' 
#' @param x An object of class \code{\link{spwb}}, \code{\link{pwb}}, \code{\link{growth}} or \code{\link{fordyn}}. 
#' The function only works with the result of simulations with \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}.
#' @param cohort An string indicating the cohort for which resistances are desired.
#' @param relative A boolean flag to indicate that relative percentages are desired as output.
#' @param draw A boolean flag to indicate that a plot should be drawn (only pathway resistances, without discriminating between soil layers).
#' @param cumulative A flag to indicate that drawn series should be cumulative.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' 
#' @details
#' The function makes internal calls to \code{\link{hydraulics_soilPlantResistancesWeibull}} or
#' \code{\link{hydraulics_soilPlantResistancesSigmoid}} depending on the value of \code{transpirationMode} in \code{x}.
#' 
#' @return 
#' If \code{draw = FALSE}, the function returns list with three items:
#' \itemize{
#'  \item{\code{pathway}: A data frame with dates in rows and resistance segments in columns (Rhizosphere, Root, Stem and Leaf). }
#'  \item{\code{root}: A data frame with dates in rows and root resistances for soil layers in columns.}
#'  \item{\code{rhizosphere}: A data frame with dates in rows and rhizosphere resistances for soil layers in columns. }
#'}
#' Values depend on whether \code{relative = TRUE} (percentages) or \code{relative = FALSE} (absolute resistance values). 
#' 
#' If \code{draw = TRUE}, a plot object is returned showing the time series of pathway resistances.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' @author \enc{Léa}{Lea} Veuillen, INRAE-URFM
#' 
#' @seealso \code{\link{waterUseEfficiency}}, \code{\link{droughtStress}}
resistances<-function(x, cohort, relative = FALSE, draw = FALSE, 
                      cumulative = FALSE, xlab = NULL, ylab=NULL) {
  
  if(!inherits(x, c("spwb","pwb","growth", "fordyn"))) {
    stop("'x' should be of class 'spwb', 'pwb', 'growth' or 'fordyn'")
  }
  if(inherits(x, c("spwb","pwb", "growth"))) {
    res <- .resistances_sim(x = x, cohort = cohort, relative = relative)
    pathmat <- res[["pathway"]]
    rootmat <- res[["root"]]
    rhizomat <- res[["rhizosphere"]]
  } else {
    vec_pathway<-vector("list", length(x$GrowthResults))
    vec_root<-vector("list", length(x$GrowthResults))
    vec_rhizo<-vector("list", length(x$GrowthResults))
    for(i in 1:length(x$GrowthResults)) {
      res <- .resistances_sim(x = x$GrowthResults[[i]], cohort = cohort, relative = relative)
      vec_pathway[[i]] <- res[["pathway"]]
      vec_root[[i]] <- res[["root"]]
      vec_rhizo[[i]] <- res[["rhizosphere"]]
    }
    pathmat <- .mergeVectorOfMatrices(vec_pathway)
    rootmat <- .mergeVectorOfMatrices(vec_root)
    rhizomat <- .mergeVectorOfMatrices(vec_rhizo)
  }
  if(draw) {
    if(is.null(ylab)) ylab = ifelse(relative, "Relative resistances (%)", "Resistances")
    if(is.null(xlab)) xlab = ""
    if(!cumulative) {
      g<-.multiple_dynamics(pathmat, ylab=ylab, xlab = xlab)
    } else {
      rescum = pathmat
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
    return(list("pathway" = pathmat, "root" = rootmat, "rhizosphere" = rhizomat))
  }
}
