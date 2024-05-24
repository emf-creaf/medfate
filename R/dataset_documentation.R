#' Example daily meteorology data
#' 
#' Example data set of meteorological input.
#' 
#' @name examplemeteo
#' @docType data
#' 
#' @format
#' A data frame containing daily meteorology of a location in Catalonia (Spain) for year 2001:
#' \describe{
#'   \item{\code{dates}}{Vector of \code{\link{Date}} objects.}
#'    \item{\code{MinTemperature}}{Minimum daily temperature (in degrees Celsius).}
#'    \item{\code{MaxTemperature}}{Maximum daily temperature (in degrees Celsius).}
#'    \item{\code{Precipitation}}{Daily precipitation (in mm of water).}
#'    \item{\code{MinRelativeHumidity}}{Minimum daily relative humidity (in percent).}
#'    \item{\code{MaxRelativeHumidity}}{Maximum daily relative humidity (in percent).}
#'    \item{\code{Radiation}}{Incoming radiation (in MJ/m2).}
#'    \item{\code{WindSpeed}}{Wind speed (in m/s).}
#' }
#' 
#' @source Interpolated from weather station data (Spanish and Catalan meteorology agencies) using package 'meteoland'.
#' 
#' @seealso \code{\link{spwb}}
#' @examples
#'  data(examplemeteo)
#' @keywords data
NULL

#' Data tables with species parameter definition and values for different countries
#' 
#' A data sets of species parameter definition and values, the latter resulting from existing databases, fit to empirical data or expert-based guesses.
#' 
#' @name SpParams
#' @aliases SpParamsDefinition SpParamsMED
#' 
#' @docType data 
#' 
#' @format
#' \itemize{
#'   \item{Data frame \code{SpParamsDefinition} has parameters in rows and columns 'ParameterName', 'ParameterGroup', 'Definition', 'Type' and 'Units'.}
#'   \item{Data frames \code{SpParamsMED} has species or genus as rows and column names equal to parameter names in \code{SpParamsDefinition}.}
#' }
#' @details
#' \code{SpParamsMED} was the official species parameter for package versions up to v.4.0.0, but will not be maintained in the future. Additional 
#' species parameter tables for different countries are distributed via package \code{medfatetraits}, available from GitHub (https://emf-creaf.github.io/medfatetraits/).
#' 
#' @examples
#' data(SpParamsDefinition)
#' data(SpParamsMED)
#' @keywords data
NULL

#' Example observed data
#'
#' Example (fake) data set of variables measured in a plot. 
#'
#' @name exampleobs
#' @docType data
#' 
#' @format
#' A data frame containing daily 'observed' values for year 2001:
#' \describe{
#'   \item{\code{dates}}{Measurement dates.}
#'   \item{\code{SWC}}{Soil moisture content (in m3/m3).}
#'   \item{\code{ETR}}{Total evapotranspiration (mm).}
#'   \item{\code{E_T1_148}}{Transpiration of Pinus halepensis cohort 'T1_148' (L/m2 of leaf area).}
#'   \item{\code{E_T2_168}}{Transpiration of Quercus ilex cohort 'T2_168' (L/m2 of leaf area).}
#'   \item{\code{FMC_T1_148}}{Fuel moisture content of Pinus halepensis cohort 'T1_148'  (in percent).}
#'   \item{\code{FMC_T2_168}}{Fuel moisture content of Quercus ilex cohort 'T2_168' (in percent).}
#'   \item{\code{BAI_T1_148}}{Basal area increment for Pinus halepensis cohort 'T1_148'  (in cm2).}
#'   \item{\code{BAI_T2_168}}{Basal area increment for Quercus ilex cohort 'T2_168' (in cm2).}
#'   \item{\code{DI_T1_148}}{Diameter increment for Pinus halepensis cohort 'T1_148'  (in cm).}
#'   \item{\code{DI_T2_168}}{Diameter increment for Quercus ilex cohort 'T2_168' (in cm).}
#'  }
#'  
#' @source
#' This data set was actually created by running a simulation and adding some gaussian error to the outputs.
#' 
#' @seealso \code{\link{evaluation}}
#' @examples
#' data(exampleobs)
#' @keywords datasets
NULL

#' Standard fuel models (Albini 1976, Scott & Burgan 2005)
#' 
#' Standard fuel models converted to metric system. Copied from package 'Rothermel' (Giorgio Vacchiano, Davide Ascoli).
#' 
#' @name SFM_metric
#' @docType data
#' @format
#'   A data frame including standard fuel models as in Albini (1976) and Scott and Burgan (2005), to be used as input of \code{\link{fire_Rothermel}} function. All values converted to metric format.
#'  \describe{
#'     \item{\code{Fuel_Model_Type}}{A factor with levels \code{D} (for dynamic) or \code{S} (for static).}
#'     \item{\code{Load_1h}}{Loading of 1h fuel class [t/ha].}
#'     \item{\code{Load_10h}}{Loading of 10h fuel class [t/ha].}
#'     \item{\code{Load_100h}}{Loading of 100h fuel class [t/ha]}
#'     \item{\code{Load_Live_Herb}}{Loading of herbaceous fuels [t/ha]}
#'     \item{\code{Load_Live_Woody}}{Loading of woody fuels [t/ha]}
#'     \item{\samp{SA/V_1h}}{Surface area to volume ratio of 1h fuel class [m2/m3]}
#'     \item{\samp{SA/V_10h}}{Surface area to volume ratio of 10h fuel class [m2/m3]}
#'     \item{\samp{SA/V_100h}}{Surface area to volume ratio of 100h fuel class [m2/m3]}
#'     \item{\samp{SA/V_Live_Herb}}{Surface area to volume ratio of herbaceous fuels [m2/m3]}
#'     \item{\samp{SA/V_Live_Woody}}{Surface area to volume ratio of woody fuels [m2/m3]}
#'     \item{\code{Fuel_Bed_Depth}}{Fuel bed depth [cm]}
#'     \item{\code{Mx_dead}}{Dead fuel moisture of extinction [percent]}
#'     \item{\code{Heat_1h}}{Heat content of 1h fuel class [kJ/kg]}
#'     \item{\code{Heat_10h}}{Heat content of 10h fuel class [kJ/kg]}
#'     \item{\code{Heat_100h}}{Heat content of 100h fuel class [kJ/kg]}
#'     \item{\code{Heat_Live_Herb}}{Heat content of herbaceous fuels [kJ/kg]}
#'     \item{\code{Heat_Live_Woody}}{Heat content of woody fuels [kJ/kg]}
#'  }
#' @source
#' Albini, F. A. (1976). Computer-based models of wildland fire behavior: A users' manual. Ogden, UT: US Department of Agriculture, Forest Service, Intermountain Forest and Range Experiment Station.
#'
#' Scott, J., and Burgan, R. E. (2005). A new set of standard fire behavior fuel models for use with Rothermel's surface fire spread model. Gen. Tech. Rep. RMRSGTR-153. Fort Collins, CO: US Department of Agriculture, Forest Service, Rocky Mountain Research Station.
#' 
#' @seealso \code{\link{fire_Rothermel}}
#' @examples
#' data(SFM_metric)
#' @keywords data
NULL


#' Example forest inventory data
#'
#' Example data to illustrate the creation of forest objects from inventory data,
#' coming from a forest inventory survey, used to illustrate the general function \code{\link{forest_mapTreeTable}}:
#' \itemize{
#'  \item \code{poblet_trees} - Data frame with example tree plot data from Poblet, Catalonia (717 observations and 4 variables).
#'    \itemize{
#'      \item Plot.Code - Plot ID (character)
#'      \item Indv.Ref - Tree individual (integer)
#'      \item Species - Species name (character)
#'      \item Diameter.cm - Tree diameter at breast height (cm)
#'    }
#' }
#' @name poblet_trees
#' @docType data
#' @source
#' \itemize{
#'   \item{Data table \code{poblet_trees} corresponds to field data sampled by the Catalan Forest Ownership Center (Centre de la Propietat Forestal; CPF).}
#'  }
#' @seealso \code{\link{forest_mapTreeTable}}
#' @keywords data
NULL