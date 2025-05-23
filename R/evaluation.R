#' Evaluation of simulations results
#' 
#' Functions to compare model predictions against observed values.
#' 
#' @param out An object of class \code{\link{spwb}}, \code{\link{growth}} or \code{\link{pwb}}.
#' @param measuredData A data frame with observed/measured values. Dates should be in row names, whereas columns should be named according to the type of output to be evaluated (see details).
#' @param type A string with the kind of model output to be evaluated. Accepted values are:
#'     \itemize{
#'       \item{\code{"SWC"}: Soil water content (percent volume). See details for specific soil layers.}
#'       \item{\code{"RWC"}: Relative water content (relative to field capacity). See details for specific soil layers.}
#'       \item{\code{"REW"}: Relative extractable water. See details for specific soil layers.}
#'       \item{\code{"ETR"}: Total evapotranspiration.}
#'       \item{\code{"SE+TR"}: Modelled soil evaporation + plant transpiration against observed total evapotranspiration}
#'       \item{\code{"E"}: Transpiration per leaf area}
#'       \item{\code{"LE"}: Latent heat (vaporisation) turbulent flux}
#'       \item{\code{"H"}: Canopy sensible heat turbulent flux}
#'       \item{\code{"GPP"}: Stand-level gross primary productivity}
#'       \item{\code{"LFMC"}: Live fuel moisture content}
#'       \item{\code{"WP"}: Plant water potentials}
#'       \item{\code{"BAI"}: Basal area increment}
#'       \item{\code{"DI"}: Diameter increment}
#'       \item{\code{"DBH"}: Diameter at breast height}
#'       \item{\code{"Height"}: Plant height}
#'     }
#' @param cohort A string of the cohort to be compared (e.g. "T1_68"). If \code{NULL} results for the first cohort will be evaluated.
#' @param temporalResolution A string to indicate the temporal resolution of the model evaluation, which can be "day", "week", "month" or "year". Observed and modelled values are aggregated temporally (using either means or sums) before comparison.
#' @param plotType Plot type to draw, either \code{"dynamics"}, \code{"pointdynamics"}  or \code{"scatter"}.
#' @param metric An evaluation metric:
#'     \itemize{
#'       \item{\code{"MAE"}: Mean absolute error.}
#'       \item{\code{"MAE.rel"}: Mean absolute error in relative terms.}
#'       \item{\code{"r"}: Pearson's linear correlation coefficient.}
#'       \item{\code{"NSE"}: Nash-Sutcliffe model efficiency coefficient.}
#'       \item{\code{"NSE.abs"}: Modified Nash-Sutcliffe model efficiency coefficient (L1 norm) (Legates & McCabe 1999).}
#'       \item{\code{"loglikelihood"}: Logarithm of the likelihood of observing the data given the model predictions, assuming independent Gaussian errors.}
#'     }
#'
#' @details Users should provide the appropriate columns in \code{measuredData}, depending on the type of output to be evaluated:
#' \itemize{
#'   \item{\code{"SWC", "RWC", "REW"}: A column with the same name should be present. By default, the first soil layer is compared. Evaluation can be done for specific soil layers, for example using "RWC.2" for the relative water content of the second layer.}
#'   \item{\code{"ETR"} or \code{"SE+TR"}: A column named \code{"ETR"} should be present, containing stand's evapotranspiration in mm/day (or mm/week, mm/month, etc, depending on the temporal resolution). If \code{type="ETR"} observed values will be compared against modelled evapotranspiration (i.e. sum of transpiration, soil evaporation and interception loss), whereas if \code{type= "SE+TR"} observed values will be compared against the sum of transpiration and soil evaporation only.}
#'   \item{\code{"LE"}: A column named \code{"LE"} should be present containing daily latent heat turbulent flux in MJ/m2.}
#'   \item{\code{"H"}: A column named \code{"H"} should be present containing daily sensible heat turbulent flux in MJ/m2.}
#'   \item{\code{"E"}: For each plant cohort whose transpiration is to be evaluated, a column starting with \code{"E_"} and continuing with a cohort name (e.g. \code{"E_T1_68"}) with transpiration in L/m2/day on a leaf area basis (or L/m2/week, L/m2/month, etc, depending on the temporal resolution).}
#'   \item{\code{"GPP"}: A column named \code{"GPP"} should be present containing daily gross primary productivity in gC/m2.}
#'   \item{\code{"LFMC"}: For each plant cohort whose transpiration is to be evaluated, a column starting with \code{"LFMC_"} and continuing with a cohort name (e.g. \code{"LFMC_T1_68"}) with fuel moisture content as percent of dry weight.}
#'   \item{\code{"WP"}: For each plant cohort whose transpiration is to be evaluated, two columns, one starting with \code{"PD_"} (for pre-dawn) and the other with \code{"MD_"} (for midday), and continuing with a cohort name (e.g. \code{"PD_T1_68"}). They should contain leaf water potential values in MPa. These are compared against sunlit water potentials.}
#'   \item{\code{"BAI"}: For each plant cohort whose growth is to be evaluated, a column starting with \code{"BAI_"} and continuing with a cohort name (e.g. \code{"BAI_T1_68"}) with basal area increment in cm2/day, cm2/week, cm2/month or cm2/year, depending on the temporal resolution.}
#'   \item{\code{"DI"}: For each plant cohort whose growth is to be evaluated, a column starting with \code{"DI_"} and continuing with a cohort name (e.g. \code{"DI_T1_68"}) with basal area increment in cm/day, cm/week, cm/month or cm/year, depending on the temporal resolution.}
#'   \item{\code{"DBH"}: For each plant cohort whose growth is to be evaluated, a column starting with \code{"DBH_"} and continuing with a cohort name (e.g. \code{"DBH_T1_68"}) with DBH values in cm.}
#'   \item{\code{"Height"}: For each plant cohort whose growth is to be evaluated, a column starting with \code{"Height_"} and continuing with a cohort name (e.g. \code{"Height_T1_68"}) with Height values in cm.}
#' }
#' Additional columns may exist with the standard error of measured quantities. These should be named as the referred quantity, followed by \code{"_err"} (e.g. \code{"PD_T1_68_err"}), and are used to draw confidence intervals around observations.
#'  
#' Row names in \code{measuredData} indicate the date of measurement (in the case of days). Alternatively, a column called \code{"dates"} can contain the measurement dates.
#' If measurements refer to months or years, row names should also be in a "year-month-day" format, although with "01" for days and/or months (e.g. "2001-02-01" for february 2001, or "2001-01-01" for year 2001).
#'
#' @return 
#' \itemize{
#'   \item{Function \code{evaluation_table} returns a data frame with dates, observed and predicted values.}
#'   \item{Function \code{evaluation_stats} returns evaluation statistics (a vector or a data frame depending on \code{type}):
#'     \itemize{
#'       \item{\code{Bias}: Mean deviation (positive values correspond to model overestimations).}
#'       \item{\code{Bias.rel}: Bias in relative terms (%).}
#'       \item{\code{MAE}: Mean absolute error.}
#'       \item{\code{MAE.rel}: Mean absolute error in relative terms (%).}
#'       \item{\code{r}: Pearson's linear correlation coefficient.}
#'       \item{\code{NSE}: Nash-Sutcliffe model efficiency coefficient.}
#'       \item{\code{NSE.abs}: Modified Nash-Sutcliffe model efficiency coefficient (L1 norm) (Legates & McCabe 1999).}
#'     }
#'   }
#'   \item{Function \code{evaluation_plot} returns a ggplot object.}
#'   \item{Function \code{evaluation_metric} returns a scalar with the desired metric.}
#' }
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
#' 
#' @references Legates, D.R., McCabe, G.J., 1999. Evaluating the use of “goodness-of-fit” measures in hydrologic and hydroclimatic model validation. Water Resour. Res. 35, 233–241.  
#' 
#' @seealso \code{\link{spwb}}, \code{\link{growth}}, \code{\link{optimization}}, \code{\link{exampleobs}}
#' 
#' @examples 
#' \donttest{
#' #Load example daily meteorological data
#' data(examplemeteo)
#' 
#' #Load example plot plant data
#' data(exampleforest)
#' 
#' #Default species parameterization
#' data(SpParamsMED)
#' 
#' #Define soil with default soil params (4 layers)
#' examplesoil <- defaultSoilParams(4)
#' 
#' #Initialize control parameters
#' control <- defaultControl("Granier")
#' 
#' #Initialize input
#' x1 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
#' 
#' #Call simulation function
#' S1 <- spwb(x1, examplemeteo, latitude = 41.82592, elevation = 100)
#' 
#' #Load observed data (in this case the same simulation results with some added error)  
#' data(exampleobs)
#' 
#' #Evaluation statistics for soil water content
#' evaluation_stats(S1, exampleobs)
#' 
#' #NSE only
#' evaluation_metric(S1, exampleobs, metric="NSE")
#' 
#' #Comparison of temporal dynamics
#' evaluation_plot(S1, exampleobs)
#' 
#' #Loglikelihood value
#' evaluation_metric(S1, exampleobs)
#' }
#' 
#' @name evaluation
evaluation_table<-function(out, measuredData, type = "SWC", cohort = NULL, 
                           temporalResolution = "day") {
  
  # allow dates in column 'dates'
  if("dates" %in% names(measuredData)) {
    row.names(measuredData) <- as.character(measuredData$dates)
  }
  # Check arguments
  temporalResolution = match.arg(temporalResolution, c("day", "week", "month", "year"))
  soil_layers_opts <- c("", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", ".10")
  soil_opts <- c(paste0("SWC", soil_layers_opts),
                 paste0("RWC", soil_layers_opts),
                 paste0("REW", soil_layers_opts))
  if("spwbInput" %in% names(out)) {
    modelInput<-out[["spwbInput"]]
    type = match.arg(type, c(soil_opts,"E", "ETR", "SE+TR", "LE", "H", "WP", "LFMC", "GPP"))
  } else {
    modelInput<- out[["growthInput"]]
    type = match.arg(type, c(soil_opts,"E", "ETR", "SE+TR", "LE", "H", "WP", "LFMC", "GPP", 
                             "BAI", "DI","DBH", "Height"))
  }
  if(substr(type,1,3) %in% c("REW", "SWC", "RWC")) {
    sm = as.data.frame(out$Soil[[substr(type,1,3)]])
    d = rownames(sm)
    if(type=="SWC")  var_mod <- "1"
    else if(type=="REW")  var_mod <- "1"
    else if(type=="RWC")  var_mod <- "1"
    else {
      var_mod  <- substr(type,5, nchar(type))
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = sm[[var_mod]])
    
    if(!(type %in% names(measuredData))) stop(paste0("Column '", type, "' not found in measured data frame."))
    seld = rownames(measuredData) %in% d
    df$Observed[d %in% rownames(measuredData)] = measuredData[seld, type]
    errcolumn = paste0(type, "_err")
    if(errcolumn %in% names(measuredData))  {
      df$obs_lower[d %in% rownames(measuredData)] = df$Observed[d %in% rownames(measuredData)] - 1.96*measuredData[[errcolumn]][seld]
      df$obs_upper[d %in% rownames(measuredData)] = df$Observed[d %in% rownames(measuredData)] + 1.96*measuredData[[errcolumn]][seld]
    }
  } 
  else if(type=="E") {
    pt = out$Plants$Transpiration
    d = rownames(pt)
    LAI = modelInput$above$LAI_live
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = pt[,icoh]/LAI[icoh])
    ## Fill observed values
    obscolumn = paste0("E_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }
  else if(type=="ETR") {
    ET2 = out$WaterBalance$Evapotranspiration
    d = rownames(out$WaterBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = ET2)
    
    if(!("ETR" %in% names(measuredData))) stop(paste0("Column 'ETR' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$ETR[rownames(measuredData) %in% d]
  }
  else if(type=="LE") {
    LEmod = out$EnergyBalance$LEVcan + out$EnergyBalance$LEVsoil
    d = rownames(out$EnergyBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = LEmod)
    if(!("LE" %in% names(measuredData))) stop(paste0("Column 'LE' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$LE[rownames(measuredData) %in% d]
  }
  else if(type=="H") {
    Hmod = out$EnergyBalance$Hcan
    d = rownames(out$EnergyBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = Hmod)
    if(!("H" %in% names(measuredData))) stop(paste0("Column 'H' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$H[rownames(measuredData) %in% d]
  }
  else if(type=="SE+TR") {
    ET1 = out$WaterBalance$SoilEvaporation+out$WaterBalance$Transpiration
    d = rownames(out$WaterBalance)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = ET1)
    
    if(!("ETR" %in% names(measuredData))) stop(paste0("Column 'ETR' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$ETR[rownames(measuredData) %in% d]
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    pdcolumn = paste0("PD_", cohort)
    mdcolumn = paste0("MD_", cohort)
    pderrcolumn = paste0("PD_", cohort, "_err")
    mderrcolumn = paste0("MD_", cohort, "_err")
    icoh = which(allcohnames==cohort)
    
    df = data.frame(Dates = as.Date(d), 
                    PD_obs = NA, MD_obs = NA, PD_obs_lower = NA, PD_obs_upper = NA, MD_obs_lower = NA, MD_obs_upper = NA,
                    PD_mod = wtPD[,icoh], MD_mod = wtMD[,icoh])
    
    seld = rownames(measuredData) %in% d
    if(pdcolumn %in% names(measuredData))  df$PD_obs[d %in% rownames(measuredData)] = measuredData[[pdcolumn]][seld]
    if(mdcolumn %in% names(measuredData))  df$MD_obs[d %in% rownames(measuredData)] = measuredData[[mdcolumn]][seld]
    if(pderrcolumn %in% names(measuredData))  {
      df$PD_obs_lower[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)] - 1.96*measuredData[[pderrcolumn]][seld]
      df$PD_obs_upper[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)] + 1.96*measuredData[[pderrcolumn]][seld]
    } else {
      df$PD_obs_lower[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)]
      df$PD_obs_upper[d %in% rownames(measuredData)] = df$PD_obs[d %in% rownames(measuredData)]
    }
    if(mderrcolumn %in% names(measuredData))  {
      df$MD_obs_lower[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)] - 1.96*measuredData[[mderrcolumn]][seld]
      df$MD_obs_upper[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)] + 1.96*measuredData[[mderrcolumn]][seld]
    } else {
      df$MD_obs_lower[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)]
      df$MD_obs_upper[d %in% rownames(measuredData)] = df$MD_obs[d %in% rownames(measuredData)]
    }
  }
  else if(type=="LFMC") {
    fmc = out$Plants$LFMC
    d = rownames(fmc)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = fmc[,icoh])
    ## Fill observed values
    obscolumn = paste0("LFMC_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }
  else if(type=="GPP") {
    GPP = rowSums(out$Plants$GrossPhotosynthesis)
    d = rownames(out$Plants$GrossPhotosynthesis)
    df = data.frame(Dates = as.Date(d), Observed = NA, Modelled = GPP)
    
    if(!("GPP" %in% names(measuredData))) stop(paste0("Column 'GPP' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData$GPP[rownames(measuredData) %in% d]
  }
  else if(type=="BAI") {
    SAg = out$GrowthMortality$SAgrowth
    d = rownames(SAg)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = SAg[,icoh])
    ## Fill observed values
    obscolumn = paste0("BAI_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }  
  else if(type=="DI") {
    DBH = out$PlantStructure$DBH
    DI = DBH - rbind(modelInput$above$DBH, DBH[-nrow(DBH),])
    d = rownames(DI)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = DI[,icoh])
    ## Fill observed values
    obscolumn = paste0("DI_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }  
  else if(type=="DBH") {
    DBH = out$PlantStructure$DBH
    d = rownames(DBH)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = DBH[,icoh])
    ## Fill observed values
    obscolumn = paste0("DBH_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }  
  else if(type=="Height") {
    Height = out$PlantStructure$Height
    d = rownames(Height)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    df <- data.frame(Dates = as.Date(d), Observed = NA, Modelled = Height[,icoh])
    ## Fill observed values
    obscolumn = paste0("Height_", cohort)
    if(!(obscolumn %in% names(measuredData))) stop(paste0("Column '", obscolumn, "' not found in measured data frame."))
    df$Observed[d %in% rownames(measuredData)] = measuredData[[obscolumn]][rownames(measuredData) %in% d] 
  }  
  if(temporalResolution != "day") {
    d.cut = cut(as.Date(d), breaks=temporalResolution)
    if(type %in% c("SWC", "REW", "LFMC", "DBH", "Height")) {
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      Observed = tapply(df$Observed, d.cut, FUN = mean, na.rm = TRUE),
                      Modelled = tapply(df$Modelled, d.cut, FUN = mean, na.rm = TRUE))
    } else if(type == "WP") {
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      PD_obs = tapply(df$PD_obs, d.cut, FUN = mean, na.rm = TRUE),
                      PD_obs_lower = tapply(df$PD_obs_lower, d.cut, FUN = mean, na.rm = TRUE),
                      PD_obs_upper = tapply(df$PD_obs_upper, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs = tapply(df$MD_obs, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs_lower = tapply(df$MD_obs_lower, d.cut, FUN = mean, na.rm = TRUE),
                      MD_obs_upper = tapply(df$MD_obs_upper, d.cut, FUN = mean, na.rm = TRUE),
                      PD_mod = tapply(df$PD_mod, d.cut, FUN = mean, na.rm = TRUE),
                      MD_mod = tapply(df$MD_mod, d.cut, FUN = mean, na.rm = TRUE))
    } else { # E, ETR, BAI
      df = data.frame(Dates = as.Date(levels(d.cut)),
                      Observed = tapply(df$Observed, d.cut, FUN = sum, na.rm = TRUE),
                      Modelled = tapply(df$Modelled, d.cut, FUN = sum, na.rm = TRUE))
    }
  }
  
  return(df)
}

#' @rdname evaluation
evaluation_stats<-function(out, measuredData, type="SWC", cohort = NULL, 
                           temporalResolution = "day") {
  evalstats<-function(obs, pred) {
    sel_complete = !(is.na(obs) | is.na(pred))
    obs = obs[sel_complete]
    pred = pred[sel_complete]
    E <- pred-obs
    Bias <- mean(E)
    Bias.rel <- 100*Bias/abs(mean(obs))
    MAE <- mean(abs(E))
    MAE.rel <- 100*MAE/abs(mean(obs))
    r<- cor(obs, pred)
    NSE <- 1 - (sum((obs-pred)^2)/sum((obs-mean(obs))^2))
    NSE.abs <- 1 - (sum(abs(obs-pred))/sum(abs(obs-mean(obs))))
    return(c(n = sum(sel_complete), Bias= Bias, Bias.rel= Bias.rel,MAE = MAE, MAE.rel = MAE.rel, r = r, NSE = NSE, NSE.abs = NSE.abs))
  }
  
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution)
  if(type=="WP") {
    eval_res = as.data.frame(rbind(evalstats(df$PD_obs, df$PD_mod),
                                   evalstats(df$MD_obs, df$MD_mod)))
    row.names(eval_res)<-c("Predawn potentials", "Midday potentials")
  } else {
    eval_res = evalstats(df$Observed, df$Modelled)
  }
  return(eval_res)
}

#' @rdname evaluation
evaluation_plot<-function(out, measuredData, type="SWC", cohort = NULL, 
                          temporalResolution = "day",
                          plotType = "dynamics") {
  scatterplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE) {
    g<-ggplot(df, aes(x=.data$Modelled))
    if(err) {
      data_err <- df[!is.na(df$obs_lower) & !is.na(df$obs_upper), , drop = FALSE]
      if(nrow(data_err)>0) {
        g<-g+
          geom_pointrange(aes(y = .data$Observed, ymin = .data$obs_lower, ymax = .data$obs_upper),cex=0.5,
                          data = data_err)
      }
    }
    g<-g + 
      geom_point(aes(y = .data$Observed), cex=0.5)+
      geom_abline(intercept=0, slope=1, col="black")+
      geom_smooth(aes(y = .data$Observed), method="lm", se = FALSE, col="gray", linetype="dashed")+
      xlab(xlab)+
      ylab(ylab)+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }
  dynamicsplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE,
                         str_obs = "Observed", str_mod = "Modelled") {
    g<-ggplot(df, aes(x=.data$Dates))
    if(err) {
      g <- g +          
        geom_ribbon(aes(ymin=.data$obs_lower, ymax=.data$obs_upper), 
                    col="gray", alpha= 0.5)
    }
    g<-g+       
      geom_path(aes(y=.data$Observed, col="Observed"))+
      geom_path(aes(y=.data$Modelled, col="Modelled"))+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", 
                         values=c("Observed"="black", "Modelled"= "red"),
                         labels =c("Observed"=str_obs, "Modelled"=str_mod))+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }
  pointdynamicsplot<-function(df, xlab="", ylab="", title=NULL, err = FALSE,
                              str_obs = "Observed", str_mod = "Modelled") {
    g<-ggplot(df, aes(x=.data$Dates))
    g<- g + geom_path(aes(y=.data$Modelled, col="Modelled"))
    if(err) {
      g <- g + geom_pointrange(aes(x = .data$Dates, y = .data$Observed, ymin = .data$obs_lower, ymax = .data$obs_upper))
    } else {
      g<- g + geom_point(aes(y=.data$Observed, col="Observed"))      
    }
    g<-g+       
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", 
                         values=c("Observed"="black", "Modelled"= "red"),
                         labels =c("Observed"=str_obs, "Modelled"=str_mod))+
      theme_bw()
    if(!is.null(title)) g<-g+labs(title=title)
    return(g)
  }
  # Check arguments
  plotType = match.arg(plotType, c("pointdynamics", "dynamics", "scatter"))
  if("spwbInput" %in% names(out)) {
    modelInput<-out[["spwbInput"]]
  } else {
    modelInput<- out[["growthInput"]]
  }
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution)
  
  if(type %in% paste0("SWC", c("", paste0(".", 1:10)))) {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = expression(paste("Soil moisture ",(m^{3}%.%m^{-3}))),
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = expression(paste("Soil moisture ",(m^{3}%.%m^{-3}))),
                           err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = expression(paste("Measured soil moisture ",(m^{3}%.%m^{-3}))),
                     ylab = expression(paste("Measured soil moisture ",(m^{3}%.%m^{-3}))), 
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  } 
  else if(type %in% paste0("RWC", c("", paste0(".", 1:10)))) {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Relative water content (RWC)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Relative water content (RWC)",
                           err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = "Modelled relative water content (RWC)",
                     ylab = "Measured relative water content (RWC)", 
                     err = (paste0(type,"_err") %in% names(measuredData)))
      
    }
  }
  else if(type %in% paste0("REW", c("", paste0(".", 1:10)))) {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Relative extractable soil water (REW)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Relative extractable soil water (REW)",
                           err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = "Modelled relative extractable soil water (REW)",
                     ylab = "Measured relative extractable soil water (REW)", 
                     err = (paste0(type,"_err") %in% names(measuredData)))
      
    }
  }
  else if(type=="E") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Transpiration per leaf area (l/m2/day)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Transpiration per leaf area (l/m2/day)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = "Modelled transpiration per leaf area (l/m2/day)",
                     ylab = "Measured transpiration per leaf area (l/m2/day)", 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="LFMC") {
    allcohnames <- row.names(modelInput$cohorts)
    spnames <- modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh <- 1
      cohort <- allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh <- which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Fuel moisture content (% of dry weight)", 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Fuel moisture content (% of dry weight)", 
                           title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = "Modelled fuel moisture content (% of dry weight)",
                     ylab = "Measured fuel moisture content (% of dry weight)", 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="BAI") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = paste0("Basal area increment (cm2/", temporalResolution,")"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = paste0("Basal area increment (cm2/", temporalResolution,")"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = paste0("Modelled basal area increment (cm2/", temporalResolution,")"),
                     ylab = paste0("Measured basal area increment (cm2/", temporalResolution,")"), 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="DI") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = paste0("Diameter increment (cm/", temporalResolution,")"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = paste0("Diameter increment (cm/", temporalResolution,")"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = paste0("Modelled diameter increment (cm/", temporalResolution,")"),
                     ylab = paste0("Measured diameter increment (cm/", temporalResolution,")"), 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="DBH") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = paste0("DBH (cm)"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = paste0("DBH (cm)"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = paste0("Modelled DBH (cm)"),
                     ylab = paste0("Measured DBH (cm)"), 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="Height") {
    allcohnames = row.names(modelInput$cohorts)
    spnames = modelInput$cohorts$Name
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = paste0("Height (cm)"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = paste0("Height (cm)"), 
                      title=paste0(cohort , " (",spnames[icoh],")"))
    } else {
      g<-scatterplot(df, 
                     xlab = paste0("Modelled height (cm)"),
                     ylab = paste0("Measured height (cm)"), 
                     title=paste0(cohort , " (",spnames[icoh],")"))
    }
  }
  else if(type=="ETR") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "ETR (mm)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "ETR (mm)",
                           err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = "Modelled ETR (mm)",
                     ylab ="Measured ETR (mm)",
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  }
  else if(type=="LE") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Latent heat (MJ/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Latent heat (MJ/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = "Modelled latent heat (MJ/m2)",
                     ylab ="Measured latent heat (MJ/m2)",
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  }
  else if(type=="H") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "Canpy sensible heat (MJ/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "Canpy sensible heat (MJ/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else{
      g<-scatterplot(df, xlab  = "Modelled sensible heat (MJ/m2)",
                     ylab ="Measured sensible heat (MJ/m2)",
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  }
  else if(type=="SE+TR") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "ETR or SE+TR (mm)", str_obs = "Observed (ETR)", str_mod = "Modelled (SE+TR)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "ETR or SE+TR (mm)", str_obs = "Observed (ETR)", str_mod = "Modelled (SE+TR)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else{
      g<-scatterplot(df, xlab  = "Modelled SE+TR (mm)",
                     ylab ="Measured ETR (mm)",
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  }
  else if(type=="GPP") {
    if(plotType=="dynamics") {
      g<-dynamicsplot(df, ylab = "GPP (gC/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else if(plotType=="pointdynamics") {
      g<-pointdynamicsplot(df, ylab = "GPP (gC/m2)",
                      err = (paste0(type,"_err") %in% names(measuredData)))
    } else {
      g<-scatterplot(df, xlab  = "Modelled GPP (gC/m2)",
                     ylab ="Measured GPP (gC/m2)",
                     err = (paste0(type,"_err") %in% names(measuredData)))
    }
  }
  else if(type=="WP"){
    wtMD = out$Plants$LeafPsiMin
    wtPD = out$Plants$LeafPsiMax
    d = rownames(wtMD)
    spnames = modelInput$cohorts$Name
    allcohnames = row.names(modelInput$cohorts)
    
    if(is.null(cohort)) {
      icoh = 1
      cohort = allcohnames[1] 
      message("Choosing first cohort")
    } else {
      icoh = which(allcohnames==cohort)
    }
    
    if(plotType=="dynamics"  || plotType == "pointdynamics"){
      g<-ggplot(df)+
        geom_path(aes(x=.data$Dates, y=.data$PD_mod, col="Predawn", linetype="Predawn"))+
        geom_path(aes(x=.data$Dates, y=.data$MD_mod, col="Midday", linetype="Midday"))+
        geom_pointrange(aes(x = .data$Dates, y = .data$PD_obs, ymin = .data$PD_obs_lower, ymax = .data$PD_obs_upper, col="Predawn", linetype="Predawn"))+
        geom_pointrange(aes(x = .data$Dates, y = .data$MD_obs, ymin = .data$MD_obs_lower, ymax = .data$MD_obs_upper, col="Midday",linetype="Midday"))+
        scale_color_manual(name="", values=c("Predawn"="blue", "Midday"= "red"))+
        scale_linetype_manual(name="", values=c("Predawn"="dashed", "Midday"= "solid"))+
        labs(title=paste0(cohort , " (",spnames[icoh],")"))+
        xlab("")+
        ylab("Leaf water potential (MPa)")+
        theme_bw()
    } else {
      g<-ggplot(df)+
        geom_abline(intercept=0, slope=1, col="black")+
        geom_pointrange(aes(x = .data$PD_mod, y = .data$PD_obs, ymin = .data$PD_obs_lower, ymax = .data$PD_obs_upper, col="Predawn"))+
        geom_pointrange(aes(x = .data$MD_mod, y = .data$MD_obs, ymin = .data$MD_obs_lower, ymax = .data$MD_obs_upper,col="Midday"))+
        geom_smooth(aes(x = .data$PD_mod, y = .data$PD_obs, col="Predawn"), method="lm", se = FALSE, linetype="dashed")+
        geom_smooth(aes(x = .data$MD_mod, y = .data$MD_obs, col="Midday"), method="lm", se = FALSE, linetype="dashed")+
        scale_color_manual(name="", values=c("Predawn"="blue", "Midday"= "red"))+
        labs(title=paste0(cohort , " (",spnames[icoh],")"))+
        xlab("Modelled leaf water potential (MPa)")+
        ylab("Measured leaf water potential (MPa)")+
        theme_bw()
    }
  }
  return(g)
}

#' @rdname evaluation
evaluation_metric<-function(out, measuredData, type="SWC", cohort=NULL, 
                            temporalResolution = "day",
                            metric = "loglikelihood") {
  df = evaluation_table(out = out, measuredData = measuredData, 
                        type = type, cohort = cohort, 
                        temporalResolution = temporalResolution)
  obs = df$Observed
  pred = df$Modelled
  sel_complete = !(is.na(obs) | is.na(pred))
  obs = obs[sel_complete]
  pred = pred[sel_complete]
  sd <- sd(obs, na.rm=TRUE)
  metric<-match.arg(metric, c("loglikelihood", "NSE", "NSE.abs", "MAE", "MAE.rel", "r"))
  m <- switch(metric,
         "loglikelihood" = sum(dnorm(obs, pred, sd, log=TRUE), na.rm=TRUE),
         "NSE" = 1 - (sum((obs-pred)^2, na.rm=TRUE)/sum((obs-mean(obs, na.rm=TRUE))^2, na.rm=TRUE)),
         "MAE" = mean(abs(pred-obs), na.rm=TRUE),
         "MAE.rel" = 100*mean(abs(pred-obs), na.rm=TRUE)/abs(mean(obs, na.rm=TRUE)),
         "r" = cor(obs, pred),
         "NSE.abs" = 1 - (sum(abs(obs-pred))/sum(abs(obs-mean(obs))))
  )
  return(m)
}
