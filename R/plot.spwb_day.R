#' Plots simulation results for one day
#' 
#' Functions to plot the sub-daily simulation results of \code{\link{spwb_day}}, \code{\link{growth_day}} 
#' or the transpiration calculations of \code{\link{transp_transpirationSperry}} or \code{\link{transp_transpirationSureau}}.
#'
#' @param x An object of class \code{spwb_day}, \code{growth_day} or \code{pwb_day}.
#' @param type The information to be plotted (see details).
#' @param bySpecies Allows aggregating output by species, before drawing plots. Aggregation can involve a sum (as for plant LAI or transpiration) 
#'                  or a LAI-weighted mean (as for plant stress or plant water potential).
#' @param xlim Range of values for x.
#' @param ylim Range of values for y.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param ... Additional parameters for function \code{plot}.
#' 
#' @details The following plots are available for \code{spwb_day} and \code{pwb_day}:
#' \itemize{
#'   \item{\code{"LeafPsi"}: Leaf water potential (for shade and sunlit leaves).}
#'   \item{\code{"LeafPsiAverage"}: Average leaf water potential.}
#'   \item{\code{"RootPsi"}: Root crown water potential.}
#'   \item{\code{"StemPsi"}: Stem water potential.}
#'   \item{\code{"StemPLC"}: (Average) percentage of loss conductance in the stem conduits.}
#'   \item{\code{"StemRWC"}: (Average) relative water content in the stem.}
#'   \item{\code{"LeafRWC"}: Relative water content in the leaf.}
#'   \item{\code{"StemSympRWC"}: (Average) relative water content in the stem symplasm.}
#'   \item{\code{"LeafSympRWC"}: Relative water content in the leaf symplasm.}
#'   \item{\code{"SoilPlantConductance"}: Overall soil plant conductance (calculated as the derivative of the supply function).}
#'   \item{\code{"PlantExtraction"}:  Water extracted from each soil layer.}
#'   \item{\code{"PlantTranspiration"}: Plant cohort transpiration per ground area.}
#'   \item{\code{"TranspirationPerLeaf"}: Plant cohort transpiration per leaf area.}
#'   \item{\code{"PlantGrossPhotosynthesis"}: Plant cohort gross photosynthesis per ground area.}
#'   \item{\code{"GrossPhotosynthesisPerLeaf"}: Plant cohort gross photosynthesis per leaf area.}
#'   \item{\code{"PlantNetPhotosynthesis"}: Plant cohort net photosynthesis per ground area.}
#'   \item{\code{"NetPhotosynthesisPerLeaf"}: Plant cohort net photosynthesis per leaf area.}
#'   \item{\code{"LeafTranspiration"}: Instantaneous transpiration per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafGrossPhotosynthesis"}: Instantaneous gross photosynthesis per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafNetPhotosynthesis"}: Instantaneous net photosynthesis per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafAbsorbedSWR"}: Absorbed short wave radiation per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafAbsorbedPAR"}: Absorbed photosynthetically-active radiation per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafNetLWR"}: Net long wave radiation per leaf area (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafCi"}: Leaf intercellular CO2 concentration (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafIntrinsicWUE"}: Leaf intrinsic water use efficiency, i.e. the ratio between instantaneous photosynthesis and stomatal conductance (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafVPD"}: Leaf vapour pressure deficit (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafStomatalConductance"}: Leaf stomatal conductance to water vapour (differentiates sunlit and shade leaves).}
#'   \item{\code{"LeafTemperature"}: Leaf temperature (differentiates sunlit and shade leaves).}
#'   \item{\code{"Temperature"}: Above-canopy, inside-canopy and soil temperature.}
#'   \item{\code{"CanopyEnergyBalance"}: Canopy energy balance components.}
#'   \item{\code{"SoilEnergyBalance"}: Soil energy balance components.}
#'   \item{\code{"PlantWaterBalance"}: Difference between water extraction from the soil and transpired water per ground area.}
#'   \item{\code{"WaterBalancePerLeaf"}: Difference between water extraction from the soil and transpired water per leaf area.}
#' }
#' And the following plots are additionally available for \code{growth_day}:
#'   \itemize{
#'     \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis rate per dry weight.}
#'     \item{\code{"MaintenanceRespiration"}: Maintenance respiration cost per dry weight.}
#'     \item{\code{"RootExudation"}: Root exudation rate per dry weight.}
#'     \item{\code{"LabileCarbonBalance"}: Labile carbon balance per dry weight.}
#'     \item{\code{"SugarLeaf"}: Sugar concentration in leaves.}
#'     \item{\code{"StarchLeaf"}: Starch concentration in leaves.}
#'     \item{\code{"SugarSapwood"}: Sugar concentration in sapwood.}
#'     \item{\code{"StarchSapwood"}: Starch concentration in sapwood.}
#'     \item{\code{"SugarTransport"}: Phloem sugar transport rate.}
#'   }
#'   
#' @note Only for soil plant water balance simulations using \code{transpirationMode = "Sperry"} or \code{transpirationMode = "Sureau"}. This function can be used to display subdaily dynamics of corresponding to single days on \code{\link{spwb}} runs, if control option \code{subdailyResults} is set to \code{TRUE}. See also option \code{subdaily} in \code{\link{plot.spwb}}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @seealso  \code{\link{spwb_day}}, \code{\link{growth_day}}, \code{\link{plot.spwb}}
#' @return An ggplot object
#' 
#' @examples 
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' #Load example plot plant data
#' data(exampleforest)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Define soil with default soil params (2 layers)
#' examplesoil <- defaultSoilParams(4)
#' 
#' #Switch to 'Sperry' transpiration mode
#' control <- defaultControl("Sperry")
#' 
#' #Simulate one day only
#' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
#' d <- 100
#' date <- examplemeteo$dates[d]
#' meteovec <- unlist(examplemeteo[d,])
#' sd2 <- spwb_day(x2, date, meteovec,
#'               latitude = 41.82592, elevation = 100, slope= 0, aspect = 0)
#' 
#' #Display transpiration for subdaily steps
#' plot(sd2, "PlantTranspiration")
#' 
#' @name plot.spwb_day
plot.spwb_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                        xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  plot.pwb_day(x, type, bySpecies, xlim, ylim, xlab, ylab,...)
}

#' @rdname plot.spwb_day
plot.growth_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                          xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  GROWTH_TYPES = .getSubdailyGROWTHPlotTypes()
  PWB_TYPES = .getSubdailySPWBPlotTypes()
  type = match.arg(type,GROWTH_TYPES)
  if(type %in% PWB_TYPES) {
    return(plot.pwb_day(x, type, bySpecies, xlim, ylim, xlab, ylab,...))
  } else {
    LCBInst = x$LabileCarbonBalanceInst
    OM = LCBInst[[type]]
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
}

#' @rdname plot.spwb_day
plot.pwb_day<-function(x, type="PlantTranspiration", bySpecies = FALSE, 
                       xlim = NULL, ylim=NULL, xlab = NULL, ylab = NULL, ...) {
  if(!("EnergyBalance" %in% names(x))) stop("Plotting function available for transpirationMode = 'Sperry' only.")
  EB = x$EnergyBalance
  Plants = x$Plants
  PlantsInst = x$PlantsInst
  SunlitLeavesInst = x$SunlitLeavesInst
  ShadeLeavesInst = x$ShadeLeavesInst
  type = match.arg(type,.getSubdailySPWBPlotTypes())  
  cohortnames = row.names(x$cohorts)
  timesteps = as.numeric(colnames(x$PlantsInst$LeafPsi))
  if(is.null(xlab)) xlab = "Time step"
  if(type=="LeafPsiAverage") {
    OM = PlantsInst$LeafPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemPsi") {
    OM = PlantsInst$StemPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafSympPsi") {
    OM = PlantsInst$LeafSympPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemSympPsi") {
    OM = PlantsInst$StemSympPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="RootPsi") {
    OM = PlantsInst$RootPsi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafPLC") {
    OM = PlantsInst$LeafPLC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemPLC") {
    OM = PlantsInst$StemPLC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemRWC") {
    OM = PlantsInst$StemRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="StemSympRWC") {
    OM = PlantsInst$StemSympRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafRWC") {
    OM = PlantsInst$LeafRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafSympRWC") {
    OM = PlantsInst$LeafSympRWC*100
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="SoilPlantConductance") {
    OM = PlantsInst$dEdP
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantWaterBalance") {
    OM = PlantsInst$PWB
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab =  .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="WaterBalancePerLeaf") {
    OM = PlantsInst$PWB
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantExtraction") {
    OM = x$ExtractionInst
    nlayers = nrow(OM)
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantTranspiration") {
    OM = PlantsInst$E
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="TranspirationPerLeaf") {
    OM = PlantsInst$E
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantGrossPhotosynthesis") {
    OM = PlantsInst$Ag
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="PlantNetPhotosynthesis") {
    OM = PlantsInst$An
    if(bySpecies) {
      OM = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
    } 
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="GrossPhotosynthesisPerLeaf") {
    OM = PlantsInst$Ag
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="NetPhotosynthesisPerLeaf") {
    OM = PlantsInst$An
    if(bySpecies) {
      m1 = apply(OM,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OM = sweep(m1,1,lai1,"/")
    } else {
      OM = sweep(OM,1,Plants$LAI,"/")
    }
    if(is.null(ylab)) ylab = .getYLab(type)
    return(.multiple_subday_dynamics(t(OM), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafLAI") {
    OM_SL = SunlitLeavesInst$LAI
    OM_SH = ShadeLeavesInst$LAI
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafTranspiration") {
    OM_SL = SunlitLeavesInst$E
    OM_SH = ShadeLeavesInst$E
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafNetPhotosynthesis") {
    OM_SL = SunlitLeavesInst$An
    OM_SH = ShadeLeavesInst$An
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafGrossPhotosynthesis") {
    OM_SL = SunlitLeavesInst$Ag
    OM_SH = ShadeLeavesInst$Ag
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafAbsorbedSWR") {
    OM_SL <- SunlitLeavesInst$Abs_SWR
    OM_SH <- ShadeLeavesInst$Abs_SWR
    if(bySpecies) {
      m1 <- apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$SunlitLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL <- m1/lai1
      m1 <- apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$ShadeLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH <- m1/lai1
    } else {
      OM_SL <- OM_SL/x$SunlitLeavesInst$LAI
      OM_SH <- OM_SH/x$ShadeLeavesInst$LAI
    }
    OM_SL[is.na(OM_SL)] <- 0.0
    OM_SH[is.na(OM_SH)] <- 0.0
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafAbsorbedPAR") {
    OM_SL <- SunlitLeavesInst$Abs_PAR
    OM_SH <- ShadeLeavesInst$Abs_PAR
    if(bySpecies) {
      m1 <- apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$ShadeLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL <- m1/lai1
      m1 <- apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$ShadeLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH <- m1/lai1
    } else {
      OM_SL <- OM_SL/x$ShadeLeavesInst$LAI
      OM_SH <- OM_SH/x$ShadeLeavesInst$LAI
    }
    OM_SL[is.na(OM_SL)] <- 0.0
    OM_SH[is.na(OM_SH)] <- 0.0
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafNetLWR") {
    OM_SL <- SunlitLeavesInst$Net_LWR
    OM_SH <- ShadeLeavesInst$Net_LWR
    if(bySpecies) {
      m1 <- apply(OM_SL,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$ShadeLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL <- m1/lai1
      m1 <- apply(OM_SH,2, tapply, x$cohorts$Name, sum, na.rm=T)
      lai1 <- apply(x$ShadeLeavesInst$LAI,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH <- m1/lai1
    } else {
      OM_SL <- OM_SL/x$ShadeLeavesInst$LAI
      OM_SH <- OM_SH/x$ShadeLeavesInst$LAI
    }
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafPsi") {
    OM_SL = SunlitLeavesInst$Psi
    OM_SH = ShadeLeavesInst$Psi
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafTemperature") {
    OM_SL = SunlitLeavesInst$Temp
    OM_SH = ShadeLeavesInst$Temp
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafVPD") {
    OM_SL = SunlitLeavesInst$VPD
    OM_SH = ShadeLeavesInst$VPD
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafStomatalConductance") {
    OM_SL = SunlitLeavesInst$Gsw
    OM_SH = ShadeLeavesInst$Gsw
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafIntrinsicWUE") {
    OM_SL = SunlitLeavesInst$An/SunlitLeavesInst$Gsw
    OM_SH = ShadeLeavesInst$An/ShadeLeavesInst$Gsw
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="LeafCi") {
    OM_SL = SunlitLeavesInst$Ci
    OM_SH = ShadeLeavesInst$Ci
    if(bySpecies) {
      lai1 = tapply(Plants$LAI, x$cohorts$Name, sum, na.rm=T)
      OMlai = sweep(OM_SL, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SL = sweep(m1,1,lai1,"/")
      OMlai = sweep(OM_SH, 1, Plants$LAI, "*")
      m1 = apply(OMlai,2, tapply, x$cohorts$Name, sum, na.rm=T)
      OM_SH = sweep(m1,1,lai1,"/")
    } 
    if(is.null(ylab)) ylab=.getYLab(type)
    return(.multiple_subday_dynamics_sunlit_shade(t(OM_SL), t(OM_SH), ylab = ylab, ylim = ylim))
  }
  else if(type=="Temperature") {
    
    if(is.null(ylab)) ylab = "Temperature (degrees C)"
    df = data.frame(row.names=timesteps)
    df[["Above-canopy"]] = EB$Temperature$Tatm
    df[["Inside-canopy"]] = EB$Temperature$Tcan
    df[["Soil"]] = EB$Temperature$Tsoil.1
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    return(g)
  } 
  else if(type=="CanopyEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})
    df = data.frame(row.names=timesteps)
    df[["Balance"]] = EB$CanopyEnergyBalance$Ebalcan
    df[["SWR abs."]] = EB$CanopyEnergyBalance$SWRcan 
    df[["LWR net"]] = EB$CanopyEnergyBalance$LWRcan
    df[["Latent heat vaporisation"]] = -EB$CanopyEnergyBalance$LEVcan
    df[["Latent heat snow fusion"]] = -EB$CanopyEnergyBalance$LEFsnow
    df[["Convection can./atm."]] = -EB$CanopyEnergyBalance$Hcan
    df[["Convection soil/can."]] = -EB$SoilEnergyBalance$Hcansoil
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    # cols = c("red","brown","orange", "blue","green", "gray", "dark gray", "black")
    # names(cols) = c("SWR abs. from atm.","LWR abs. from atm.","LWR abs. from soil","LWR emmited", "Latent heat (L)",
    #                 "Convection can./atm.","Convection soil/can.", "Balance")
    return(g)
  }
  else if(type=="SoilEnergyBalance") {
    if(is.null(ylab)) ylab = expression(W%.%m^{-2})    
    df = data.frame(row.names=timesteps)
    df[["Balance"]] = EB$SoilEnergyBalance$Ebalsoil
    df[["SWR abs."]] = EB$SoilEnergyBalance$SWRsoil
    df[["LWR net"]] = EB$SoilEnergyBalance$LWRsoil
    df[["Convection soil/can."]] = EB$SoilEnergyBalance$Hcansoil
    df[["Latent heat vaporisation"]] = -EB$SoilEnergyBalance$LEVsoil
    g<-.multiple_subday_dynamics(as.matrix(df), ylab=ylab, ylim = ylim)
    
    # legend("topright", bty="n", col=c("red","brown","orange", "blue", "gray", "black"), lty=1,
    #        legend=c("SWR abs. from atm.","LWR abs. from atm.", "LWR abs. from canopy","LWR emmited", 
    #                   "Convection soil/can.", "Balance"),...)        
    return(g)
  }
}
