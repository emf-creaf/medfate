// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"
#include "communication_structures.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "windextinction_c.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "forestutils.h"
#include "modelInput.h"
#include "photosynthesis.h"
#include "fuelstructure.h"
#include "firebehaviour.h"
#include "tissuemoisture.h"
#include "spwb_day_c.h"
#include "soil.h"
#include "meteoland/utils_c.hpp"
#include "meteoland/radiation_c.hpp"
#include "meteoland/pet_c.hpp"
using namespace Rcpp;


//' Single-day soil-plant water balance
//'
//' Function \code{spwb_day} performs water balance for a single day.
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' @param date Date as string "yyyy-mm-dd".
//' @param meteovec A named numerical vector with weather data. See variable names in parameter \code{meteo} of \code{\link{spwb}}.
//' @param latitude Latitude (in degrees).
//' @param elevation,slope,aspect Elevation above sea level (in m), slope (in degrees) and aspect (in degrees from North). 
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object is allowed to be modified during the simulation.
//' 
//' @details
//' Detailed model description is available in the medfate book. 
//' 
//' Soil-plant water balance simulations allow using different sub-models for bulk soil water flows and different sub-models of transpiration and photosynthesis: 
//' 
//' (1) Sub-models of transpiration and photosynthesis (control parameter \code{transpirationMode}):
//' \itemize{
//'   \item{The sub-model corresponding to 'Granier' transpiration mode is illustrated by internal function \code{\link{transp_transpirationGranier}} and was described in De Caceres et al. (2015),
//'   and implements an approach originally described in Granier et al. (1999).} 
//'   \item{The sub-model corresponding to 'Sperry' transpiration mode is illustrated by internal function \code{\link{transp_transpirationSperry}} and was described in De Caceres et al. (2021), and
//'   implements a modelling approach originally described in Sperry et al. (2017).}  
//'   \item{The sub-model corresponding to 'Sureau' transpiration mode is illustrated by internal function \code{\link{transp_transpirationSureau}} and was described for model SurEau-Ecos v2.0 in Ruffault et al. (2022).} 
//' }
//' Simulations using the 'Sperry' or 'Sureau' transpiration mode are computationally much more expensive than 'Granier' because they
//' include explicit plant hydraulics and canopy/soil energy balance at subdaily time steps.
//' 
//' (2) Sub-models of bulk soil water flows (control parameter \code{soilDomains}; see internal function \code{\link{hydrology_soilWaterBalance}}):
//'  \itemize{
//'    \item{The submodel corresponding to 'buckets' is a multi-bucket model adds/substracts water to each layer and if content is above field capacity the excess percolates to the layer below.}
//'    \item{The submodel corresponding to 'single' is a single-domain solving 1D Richards equation (Bonan et al. 2019).}
//'    \item{The submodel corresponding to 'dual' is a dual-permeability model from MACRO 5.0 (Jarvis et al. 1991; Larsbo et al. 2005)}
//'  }
//'  
//' @return
//' Function \code{spwb_day()} returns a list of class \code{spwb_day} with the 
//' following elements:
//' \itemize{
//'   \item{\code{"cohorts"}: A data frame with cohort information, copied from \code{\link{spwbInput}}.}
//'   \item{\code{"topography"}: Vector with elevation, slope and aspect given as input.} 
//'   \item{\code{"weather"}: A vector with the input weather.}
//'   \item{\code{"WaterBalance"}: A vector of water balance components (rain, snow, net rain, infiltration, ...) for the simulated day, equivalent to one row of 'WaterBalance' object given in \code{\link{spwb}}.}
//'   \item{\code{"EnergyBalance"}: Energy balance of the stand (only returned when \code{transpirationMode = "Sperry"} or  \code{transpirationMode = "Sureau"}; see \code{\link{transp_transpirationSperry}}).}
//'   \item{\code{"Soil"}: A data frame with results for each soil layer:
//'     \itemize{
//'       \item{\code{"Psi"}: Soil water potential (in MPa) at the end of the day.}
//'       \item{\code{"HerbTranspiration"}: Water extracted by herbaceous plants from each soil layer (in mm).}
//'       \item{\code{"HydraulicInput"}: Water entering each soil layer from other layers, transported via plant roots (in mm).}
//'       \item{\code{"HydraulicOutput"}: Water leaving each soil layer (going to other layers or the transpiration stream) (in mm).}
//'       \item{\code{"PlantExtraction"}: Water extracted by woody plants from each soil layer (in mm).}
//'     }
//'   }
//'   \item{\code{"Stand"}: A named vector with with stand values for the simulated day, equivalent to one row of 'Stand' object returned by \code{\link{spwb}}.}
//'   \item{\code{"Plants"}: A data frame of results for each plant cohort (see \code{\link{transp_transpirationGranier}} or \code{\link{transp_transpirationSperry}}).}
//' }
//' 
//' The following additional items are only returned when \code{transpirationMode = "Sperry"} or  \code{transpirationMode = "Sureau"}:
//' \itemize{
//'   \item{\code{"RhizoPsi"}: Minimum water potential (in MPa) inside roots, after crossing rhizosphere, per cohort and soil layer.}
//'   \item{\code{"SunlitLeaves"} and \code{"ShadeLeaves"}: For each leaf type, a data frame with values of LAI, Vmax298 and Jmax298 for leaves of this type in each plant cohort.}
//'   \item{\code{"ExtractionInst"}: Water extracted by each plant cohort during each time step.}
//'   \item{\code{"PlantsInst"}: A list with instantaneous (per time step) results for each plant cohort (see \code{\link{transp_transpirationSperry}}).}
//'   \item{\code{"LightExtinction"}: A list of information regarding radiation balance through the canopy, as returned by function \code{\link{light_instantaneousLightExtinctionAbsortion}}.}
//'   \item{\code{"CanopyTurbulence"}: Canopy turbulence (see \code{\link{wind_canopyTurbulence}}).}
//' }
//'   
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology 213: 77-90 (doi:10.1016/j.agrformet.2015.06.012).
//' 
//' De \enc{Cáceres}{Caceres} M, Mencuccini M, Martin-StPaul N, Limousin JM, Coll L, Poyatos R, Cabon A, Granda V, Forner A, Valladares F, \enc{Martínez}{Martinez}-Vilalta J (2021) Unravelling the effect of species mixing on water use and drought stress in holm oak forests: a modelling approach. Agricultural and Forest Meteorology 296 (doi:10.1016/j.agrformet.2020.108233).
//' 
//' Granier A, \enc{Bréda}{Breda} N, Biron P, Villette S (1999) A lumped water balance model to evaluate duration and intensity of drought constraints in forest stands. Ecol Modell 116:269–283. https://doi.org/10.1016/S0304-3800(98)00205-1.
//' 
//' Jarvis, N.J., Jansson, P‐E., Dik, P.E. & Messing, I. (1991). Modelling water and solute transport in macroporous soil. I. Model description and sensitivity analysis. Journal of Soil Science, 42, 59–70.
//' 
//' Larsbo, M., Roulier, S., Stenemo, F., Kasteel, R. & Jarvis, N. (2005). An Improved Dual‐Permeability Model of Water Flow and Solute Transport in the Vadose Zone. Vadose Zone Journal, 4, 398–406.
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
//' @seealso
//' \code{\link{spwbInput}}, \code{\link{spwb}},  \code{\link{plot.spwb_day}},  
//' \code{\link{growth_day}}  
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
//' #Define soil parameters
//' examplesoil <- defaultSoilParams(4)
//' 
//' # Day to be simulated
//' d <- 100
//' meteovec <- unlist(examplemeteo[d,-1])
//' date <- as.character(examplemeteo$dates[d])
//' 
//' #Simulate water balance one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x1 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd1 <- spwb_day(x1, date, meteovec,  
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0) 
//' 
//' #Simulate water balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x2 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
//' sd2 <-spwb_day(x2, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water balance for one day only (Sureau mode)
//' control <- defaultControl("Sureau")
//' x3 <- spwbInput(exampleforest, examplesoil, SpParamsMED, control)
//' sd3 <-spwb_day(x3, date, meteovec,
//'               latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' 
//' @name spwb_day
// [[Rcpp::export("spwb_day")]]
List spwbDay(List x, CharacterVector date, NumericVector meteovec, 
             double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
             double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
             bool modifyInput = true) {
  
  //Check if input version is lower than current medfate version. If so, try to complete fields
  if(isLowerVersion(x)) spwbInputVersionUpdate(x);

  WeatherInputVector meteovec_c(meteovec);
  ModelInput x_c(x);  

  
  
  // Build communication structures
  int nlayers = x_c.soil.getNlayers();
  int ncanlayers = x_c.canopy.zlow.size();
  int numCohorts = x_c.cohorts.SpeciesIndex.size();
  int ntimesteps = x_c.control.advancedWB.ndailysteps;
  WBCommunicationStructures WBcomm(numCohorts, nlayers, ncanlayers, ntimesteps);

  // Prepare lateral flows
  std::vector<double> lateralFlows_c(nlayers, 0.0);
  NumericVector lateralFlows_mm;
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<lateralFlows_mm.size();l++) {
      lateralFlows_c[l] = lateralFlows_mm[l];
    }
  }
  List l;
  
  if(x_c.control.transpirationMode=="Granier") {
    //Initialises a result
    BasicTranspiration_RESULT BTres(numCohorts, nlayers);
    BasicSPWB_RESULT BSPWBres(BTres, nlayers);
    // Calls simulation
    wb_day_inner_c(BSPWBres, WBcomm, x_c, 
                   as<std::string>(date[0]),
                   meteovec_c, 
                   latitude, elevation, slope, aspect,
                   runon, 
                   lateralFlows_c, waterTableDepth);
    //Copies result
    l = copySPWBResult_c(BSPWBres, x_c);
    //Modifies input list, if required  
    if(modifyInput) x_c.copyStateToList(x);
  } else {
    //Initialises a result
    AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
    AdvancedSPWB_RESULT ASPWBres(ATres, nlayers);
    // Calls simulation
    wb_day_inner_c(ASPWBres, WBcomm, x_c, 
                   as<std::string>(date[0]),
                   meteovec_c, 
                   latitude, elevation, slope, aspect,
                   runon, 
                   lateralFlows_c, waterTableDepth);
    //Copies result
    l = copySPWBResult_c(ASPWBres, x_c);
    //Modifies input list, if required  
    if(modifyInput) {
      x_c.copyStateToList(x); 
    }
  }

  return(l);
}



//' @rdname aspwb
//' @keywords internal
// [[Rcpp::export("aspwb_day")]]
List aspwbDay(List x, CharacterVector date, NumericVector meteovec, 
              double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,
              double runon =  0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
              bool modifyInput = true) {
   
   WeatherInputVector meteovec_c(meteovec);
   AgricultureModelInput x_c(x);  
   int nlayers = x_c.soil.getNlayers();
   
   // Prepare lateral flows
   std::vector<double> lateralFlows_c(nlayers, 0.0);
   NumericVector lateralFlows_mm;
   if(lateralFlows.isNotNull()) {
     lateralFlows_mm = NumericVector(lateralFlows);
     for(int l=0;l<lateralFlows_mm.size();l++) {
       lateralFlows_c[l] = lateralFlows_mm[l];
     }
   }
   
   //Initialises a result
   AgricultureWB_RESULT AgrWBres(nlayers);
   SoilWaterBalance_COMM SWBcomm(nlayers);
   WBCommunicationStructures WBcomm(0, nlayers, 0, 0);
   
   // Calls simulation
   wb_day_inner_c(AgrWBres, WBcomm, x_c, 
                  as<std::string>(date[0]),
                  meteovec_c, 
                  latitude, elevation, slope, aspect,
                  runon, 
                  lateralFlows_c, waterTableDepth);
   //Copies result
   List l = copyAgricultureWBResult_c(AgrWBres, x_c);
   
   if(modifyInput) x_c.copyStateToList(x);
   return(l);
 }
