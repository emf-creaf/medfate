#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"
#include "lightextinction_basic.h"
#include "windextinction_c.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "hydraulics_c.h"
#include "forestutils.h"
#include "tissuemoisture_c.h"
#include "carbon.h"
#include "carbon_c.h"
#include "photosynthesis.h"
#include "root.h"
#include "soil.h"
#include "soil_c.h"
#include "transpiration_basic_c.h"
#include "spwb.h"
#include "meteoland/utils_c.hpp"
#include "meteoland/radiation_c.hpp"
#include "meteoland/pet_c.hpp"
using namespace Rcpp;

//' Transpiration modes
//' 
//' High-level sub-models representing transpiration, plant hydraulics, photosynthesis and water relations 
//' within plants. 
//' 
//' Three sub-models are available: 
//' \itemize{
//'   \item{Sub-model in function \code{transp_transpirationGranier} was described in De \enc{Cáceres}{Caceres} et al. (2015), 
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{Sub-model in function \code{transp_transpirationSperry} was described in De \enc{Cáceres}{Caceres} et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).} 
//'   \item{Sub-model in function \code{transp_transpirationSureau} was described for SurEau-Ecos v2.0 model in Ruffault et al. (2022).} 
//' }
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}, built using the 'Granier', 'Sperry' or 'Sureau' transpiration modes.
//' @param meteo A data frame with daily meteorological data series (see \code{\link{spwb}}).
//' @param day An integer to identify a day (row) within the \code{meteo} data frame.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North).
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @return
//' A list with the following elements:
//' \itemize{
//'   \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
//'   \item{\code{"Stand"}: A vector of stand-level variables.}
//'   \item{\code{"Plants"}: A data frame of results for each plant cohort. When using \code{transp_transpirationGranier}, element \code{"Plants"} includes:
//'     \itemize{
//'       \item{\code{"LAI"}: Leaf area index of the plant cohort.}
//'       \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
//'       \item{\code{"AbsorbedSWRFraction"}: Fraction of SWR absorbed by each cohort.}
//'       \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
//'       \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"psi"}: Water potential (in MPa) of the plant cohort (average over soil layers).}
//'       \item{\code{"DDS"}: Daily drought stress \[0-1\] (relative whole-plant conductance).}
//'     }
//'   When using \code{transp_transpirationSperry} or \code{transp_transpirationSureau}, element \code{"Plants"} includes:
//'     \itemize{
//'       \item{\code{"LAI"}: Leaf area index of the plant cohort.}
//'       \item{\code{"LAIlive"}: Leaf area index of the plant cohort, assuming all leaves are unfolded.}
//'       \item{\code{"Extraction"}: Water extracted from the soil (in mm) for each cohort.}
//'       \item{\code{"Transpiration"}: Transpirated water (in mm) corresponding to each cohort.}
//'       \item{\code{"GrossPhotosynthesis"}: Gross photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"NetPhotosynthesis"}: Net photosynthesis (in gC/m2) corresponding to each cohort.}
//'       \item{\code{"RootPsi"}: Minimum water potential (in MPa) at the root collar.}
//'       \item{\code{"StemPsi"}: Minimum water potential (in MPa) at the stem.}
//'       \item{\code{"StemPLC"}: Proportion of conductance loss in stem.}
//'       \item{\code{"LeafPsiMin"}: Minimum (predawn) water potential (in MPa) at the leaf (representing an average leaf).}
//'       \item{\code{"LeafPsiMax"}: Maximum (midday) water potential (in MPa) at the leaf (representing an average leaf).}
//'       \item{\code{"LeafPsiMin_SL"}: Minimum (predawn) water potential (in MPa) at sunlit leaves.}
//'       \item{\code{"LeafPsiMax_SL"}: Maximum (midday) water potential (in MPa) at sunlit leaves.}
//'       \item{\code{"LeafPsiMin_SH"}: Minimum (predawn) water potential (in MPa) at shade leaves.}
//'       \item{\code{"LeafPsiMax_SH"}: Maximum (midday) water potential (in MPa) at shade leaves.}
//'       \item{\code{"dEdP"}: Overall soil-plant conductance (derivative of the supply function).}
//'       \item{\code{"DDS"}: Daily drought stress \[0-1\] (relative whole-plant conductance).}
//'       \item{\code{"StemRWC"}: Relative water content of stem tissue (including symplasm and apoplasm).}
//'       \item{\code{"LeafRWC"}: Relative water content of leaf tissue (including symplasm and apoplasm).}
//'       \item{\code{"LFMC"}: Live fuel moisture content (in percent of dry weight).}
//'       \item{\code{"WaterBalance"}: Plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"Extraction"}: A data frame with mm of water extracted from each soil layer (in columns) by each cohort (in rows). The sum of a given row is equal to the total extraction of the corresponding plant cohort.}
//'   \item{\code{"ExtractionPools"}: A named list with as many elements as plant cohorts, where each element is a matrix data with mm of water extracted from each layer (in columns) of the water pool of each cohort (in rows). The sum of a given matrix is equal to the total extraction of the corresponding plant cohort.}
//' 
//'   The remaining items are only given by \code{transp_transpirationSperry} or \code{transp_transpirationSureau}:
//'   \item{\code{"EnergyBalance"}: A list with the following elements:
//'     \itemize{
//'       \item{\code{"Temperature"}: A data frame with the temperature of the atmosphere ('Tatm'), canopy ('Tcan') and soil ('Tsoil.1', 'Tsoil.2', ...) for each time step.}
//'       \item{\code{"CanopyEnergyBalance"}: A data frame with the components of the canopy energy balance (in W/m2) for each time step.}
//'       \item{\code{"SoilEnergyBalance"}: A data frame with the components of the soil energy balance (in W/m2) for each time step.}
//'     }  
//'   }
//'   \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
//'   \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Data frames for sunlit leaves and shade leaves and the following columns per cohort:
//'     \itemize{
//'       \item{\code{"LAI"}: Cumulative leaf area index of sunlit/shade leaves.}
//'       \item{\code{"Vmax298"}: Average maximum carboxilation rate for sunlit/shade leaves.}
//'       \item{\code{"Jmax298"}: Average maximum electron transport rate for sunlit/shade leaves.}
//'     }  
//'   }
//'   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
//'   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort:
//'     \itemize{
//'       \item{\code{"E"}: A data frame with the cumulative transpiration (mm) for each plant cohort during each time step. }
//'       \item{\code{"Ag"}: A data frame with the cumulative gross photosynthesis (gC/m2) for each plant cohort during each time step. }
//'       \item{\code{"An"}: A data frame with the cumulative net photosynthesis (gC/m2) for each plant cohort during each time step. }
//'       \item{\code{"Sunlitleaves"} and \code{"ShadeLeaves"}: Lists with instantaneous (for each time step) results for sunlit leaves and shade leaves and the following items:
//'         \itemize{
//'           \item{\code{"Abs_SWR"}: A data frame with instantaneous absorbed short-wave radiation (SWR).} 
//'           \item{\code{"Net_LWR"}: A data frame with instantaneous net long-wave radiation (LWR).} 
//'           \item{\code{"An"}: A data frame with instantaneous net photosynthesis (in micromol/m2/s). }
//'           \item{\code{"Ci"}: A data frame with instantaneous intercellular CO2 concentration (in ppm). }
//'           \item{\code{"GW"}: A data frame with instantaneous stomatal conductance (in mol/m2/s). }
//'           \item{\code{"VPD"}: A data frame with instantaneous vapour pressure deficit (in kPa). }
//'           \item{\code{"Temp"}: A data frame with leaf temperature (in degrees Celsius). }
//'           \item{\code{"Psi"}: A data frame with leaf water potential (in MPa). }
//'         }
//'       }
//'       \item{\code{"dEdP"}: A data frame with the slope of the plant supply function (an estimation of whole-plant conductance).}
//'       \item{\code{"RootPsi"}: A data frame with root crown water potential (in MPa) for each plant cohort during each time step.}
//'       \item{\code{"StemPsi"}: A data frame with stem water potential (in MPa) for each plant cohort during each time step.}
//'       \item{\code{"LeafPsi"}: A data frame with leaf (average) water potential (in MPa) for each plant cohort during each time step. }
//'       \item{\code{"StemPLC"}: A data frame with the proportion loss of conductance \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"StemRWC"}: A data frame with the (average) relative water content of stem tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"LeafRWC"}: A data frame with the relative water content of leaf tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"StemSympRWC"}: A data frame with the (average) relative water content of symplastic stem tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"LeafSympRWC"}: A data frame with the relative water content of symplastic leaf tissue \[0-1\] for each plant cohort during each time step. }
//'       \item{\code{"PWB"}: A data frame with plant water balance (extraction - transpiration).}
//'     }
//'   }
//'   \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
//'   \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
//'   \item{\code{"SupplyFunctions"}: If \code{stepFunctions} is not missing, a list of supply functions, photosynthesis functions and profit maximization functions.}
//' }
//' 
//' @references
//' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
//' 
//' De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).
//' 
//' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
//' 
//' Ruffault J, Pimont F, Cochard H, Dupuy JL, Martin-StPaul N (2022) 
//' SurEau-Ecos v2.0: a trait-based plant hydraulics model for simulations of plant water status and drought-induced mortality at the ecosystem level.
//' Geoscientific Model Development 15, 5593-5626 (doi:10.5194/gmd-15-5593-2022).
//' 
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author 
//' \itemize{
//'   \item{Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF}
//'   \item{Nicolas Martin-StPaul, URFM-INRAE}
//' }
//' 
//' @seealso \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Define soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' 
//' #Initialize input
//' x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Granier's model, plant water potential 
//' # and plant stress for a given day
//' t1 <- transp_transpirationGranier(x1, examplemeteo, 1, 
//'                                  latitude = 41.82592, elevation = 100, slope = 0, aspect = 0, 
//'                                  modifyInput = FALSE)
//' 
//' #Switch to 'Sperry' transpiration mode
//' control <- defaultControl("Sperry")
//' 
//' #Initialize input
//' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Sperry's model
//' t2 <- transp_transpirationSperry(x2, examplemeteo, 1, 
//'                                 latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
//'                                 modifyInput = FALSE)
//'                                 
//' #Switch to 'Sureau' transpiration mode
//' control <- defaultControl("Sureau")
//' 
//' #Initialize input
//' x3 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Transpiration according to Sureau model
//' t3 <- transp_transpirationSureau(x3, examplemeteo, 1, 
//'                                   latitude = 41.82592, elevation = 100, slope = 0, aspect = 0,
//'                                   modifyInput = FALSE)
//'                                 
//' @name transp_modes
//' @keywords internal
// [[Rcpp::export("transp_transpirationGranier")]]
List transpirationGranier(List x, DataFrame meteo, int day,
                          double latitude, double elevation, double slope, double aspect, 
                          bool modifyInput = true) {

  // Rcpp::Rcout << "In transpiration Granier\n";
  List control = x["control"];
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  NumericVector Radiation = meteo["Radiation"];
  NumericVector WindSpeed(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector CO2(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) CO2 = meteo["CO2"];
  NumericVector Patm(MinTemperature.length(), NA_REAL);
  if(meteo.containsElementNamed("Patm")) Patm = meteo["Patm"];
  
  if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
  double latrad = latitude * (M_PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  CharacterVector dateStrings = getWeatherDates(meteo);
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = julianDay_c(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));

  double tmin = MinTemperature[day-1];
  double tmax = MaxTemperature[day-1];
  double tday = averageDaylightTemperature_c(tmin, tmax);
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double rad = Radiation[day-1];
  double wind = WindSpeed[day-1];
  double Catm = CO2[day-1];

  double pet = PenmanPET_c(latrad, elevation, slorad, asprad, J, 
                             tmin, tmax, rhmin, rhmax, rad, wind);
  
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];
  
  WeatherInputVector meteovec;
  meteovec.tmax = tmax,
  meteovec.tmin = tmin,
  meteovec.rhmin = rhmin, 
  meteovec.rhmax = rhmax,
  meteovec.tday = tday, 
  meteovec.pet = pet,
  meteovec.Catm = Catm,
  meteovec.Patm = Patm[day-1];
  
  // Rcpp::Rcout << "Building input\n";
  ModelInput x_c = ModelInput(x);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = x_c.soil.getNlayers();
  int numCohorts = x_c.cohorts.SpeciesIndex.size();
  int ncanlayers = x_c.canopy.zlow.size();
  
  //Create comunication structures
  BasicTranspiration_RESULT BTres = BasicTranspiration_RESULT(numCohorts, nlayers);
  
  BasicTranspiration_COMM BTcomm = BasicTranspiration_COMM(numCohorts, ncanlayers, nlayers);
  
  //Perform simulation
  // Rcpp::Rcout << "Transpiration\n";
  transpirationBasic_c(BTres, BTcomm, x_c, 
                       meteovec, elevation);

  //Copy output to Rcpp structures
  // Rcpp::Rcout << "Copying results\n";
  List transpBasic = copyBasicTranspirationResult_c(BTres, x_c);
    
  if(modifyInput) {
    // Modify all state variables of input object from structure
    x_c.copyStateToList(x);
  }
  return(transpBasic);
} 


