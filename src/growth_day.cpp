// [[Rcpp::interfaces(r,cpp)]]
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"
#include "carbon.h"
#include "carbon_c.h"
#include "forestutils.h"
#include "forestutils_c.h"
#include "fireseverity_c.h"
#include "firebehaviour.h"
#include "firebehaviour_c.h"
#include "growth_day_c.h"
#include "hydraulics.h"
#include "hydraulics_c.h"
#include "hydrology.h"
#include "lightextinction_basic.h"
#include "lightextinction_advanced.h"
#include "modelInput.h"
#include "root.h"
#include "root_c.h"
#include "soil.h"
#include "spwb.h"
#include "tissuemoisture.h"
#include "tissuemoisture_c.h"
#include "woodformation_c.h"
#include "meteoland/utils_c.hpp"
#include "meteoland/radiation_c.hpp"
#include "meteoland/pet_c.hpp"
using namespace Rcpp;

//' Single-day forest growth
//'
//' Function \code{growth_day} performs water and carbon balance for a single day.
//' 
//' @param x An object of class \code{\link{growthInput}}.
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
//' Forest growth simulations allow using different sub-models for bulk soil water flows and different sub-models of transpiration and photosynthesis (see details in \code{\link{spwb_day}}). 
//' 
//' @return
//' Function \code{growth_day()} returns a list of class \code{growth_day} with the 
//' same elements as \code{\link{spwb_day}} and the following:
//' \itemize{
//'   \item{\code{"CarbonBalance"}: A vector of different stand-level carbon balance components (gross primary production, maintenance respiration, synthesis respiration, net primary production, heterotrophic respiration and net ecosystem exchange), all in g C · m-2.}
//'   \item{\code{"LabileCarbonBalance"}: A data frame with labile carbon balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"GrossPhotosynthesis"}: Daily gross photosynthesis per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"MaintentanceRespiration"}: Daily maintenance respiration per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"GrowthCosts"}: Daily growth costs per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"RootExudation"}: Root exudation per dry weight of living biomass (g gluc · g dry-1).}    
//'     \item{\code{"LabileCarbonBalance"}: Daily labile carbon balance (photosynthesis - maintenance respiration - growth costs - root exudation) per dry weight of living biomass (g gluc · g dry-1).}
//'     \item{\code{"SugarLeaf"}: Sugar concentration (mol·l-1) in leaves.}
//'     \item{\code{"StarchLeaf"}: Starch concentration (mol·l-1) in leaves.}
//'     \item{\code{"SugarSapwood"}: Sugar concentration (mol·l-1) in sapwood.}
//'     \item{\code{"StarchSapwood"}: Starch concentration (mol·l-1) in sapwood.}
//'     \item{\code{"SugarTransport"}:  Average instantaneous rate of carbon transferred between leaves and stem compartments via floem (mol gluc·s-1).}
//'   }
//'   \item{\code{"PlantBiomassBalance"}: A data frame with plant biomass balance results for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"StructuralBiomassBalance"}: Daily structural biomass balance (g dry · m-2).}
//'     \item{\code{"LabileBiomassBalance"}: Daily labile biomass balance (g dry · m-2).}
//'     \item{\code{"PlantBiomassBalance"}: Daily plant biomass balance, i.e. labile change + structural change (g dry · m-2).}
//'     \item{\code{"MortalityBiomassLoss"}: Biomass loss due to mortality (g dry · m-2).}    
//'     \item{\code{"CohortBiomassBalance"}: Daily cohort biomass balance (including mortality) (g dry · m-2).}
//'   }
//'   \item{\code{"PlantStructure"}: A data frame with area and biomass values for compartments of plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LeafBiomass"}: Leaf structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodBiomass"}: Sapwood structural biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootBiomass"}: Fine root biomass (in g dry) for an average individual of each plant cohort.}
//'     \item{\code{"LeafArea"}: Leaf area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"SapwoodArea"}: Sapwood area (in cm2) for an average individual of each plant cohort.}
//'     \item{\code{"FineRootArea"}: Fine root area (in m2) for an average individual of each plant cohort.}
//'     \item{\code{"HuberValue"}: Sapwood area to (target) leaf area (in cm2/m2).}
//'     \item{\code{"RootAreaLeafArea"}: The ratio of fine root area to (target) leaf area (in m2/m2).}
//'     \item{\code{"DBH"}: Diameter at breast height (in cm) for an average individual of each plant cohort.}
//'     \item{\code{"Height"}: Height (in cm) for an average individual of each plant cohort.}
//'   }
//'   \item{\code{"GrowthMortality"}: A data frame with growth and mortality rates for plant cohorts, with elements:}
//'   \itemize{
//'     \item{\code{"LAgrowth"}: Leaf area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"SAgrowth"}: Sapwood area growth rate (in cm2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"FRAgrowth"}: Fine root area growth (in m2·day-1) for an average individual of each plant cohort.}
//'     \item{\code{"StarvationRate"}: Mortality rate from starvation (ind/d-1).}
//'     \item{\code{"DessicationRate"}: Mortality rate from dessication (ind/d-1).}
//'     \item{\code{"MortalityRate"}: Mortality rate (any cause) (ind/d-1).}
//'   }
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
//' @seealso
//' \code{\link{spwb_day}}, \code{\link{growthInput}}, \code{\link{growth}},
//' \code{\link{plot.growth_day}}  
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
//' #Simulate water and carbon balance for one day only (Granier mode)
//' control <- defaultControl("Granier")
//' x4  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd4 <- growth_day(x4, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sperry mode)
//' control <- defaultControl("Sperry")
//' x5  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd5 <- growth_day(x5, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' #Simulate water and carbon balance for one day only (Sureau mode)
//' control <- defaultControl("Sureau")
//' x6  <- growthInput(exampleforest,examplesoil, SpParamsMED, control)
//' sd6 <- growth_day(x6, date, meteovec,
//'                 latitude = 41.82592, elevation = 100, slope=0, aspect=0)
//' 
//' @name growth_day
// [[Rcpp::export("growth_day")]]
List growthDay(List x, CharacterVector date, NumericVector meteovec, 
                 double latitude, double elevation, double slope = NA_REAL, double aspect = NA_REAL,  
                 double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                 bool modifyInput = true) {
  
  //Check if input version is lower than current medfate version. If so, try to complete fields
  if(isLowerVersion(x)) growthInputVersionUpdate(x);
  
  WeatherInputVector meteovec_c(meteovec);
  ModelInput x_c(x);  
  
  // Build communication structures
  int nlayers = x_c.soil.getNlayers();
  int ncanlayers = x_c.canopy.zlow.size();
  int numCohorts = x_c.cohorts.SpeciesIndex.size();
  int ntimesteps = x_c.control.advancedWB.ndailysteps;
  
  // Rcpp::Rcout << "Defining communication structures\n";
  GROWTHCommunicationStructures GROWTHcomm(numCohorts, nlayers, ncanlayers, ntimesteps);
  
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
    BasicGROWTH_RESULT GROWTHres(BSPWBres, numCohorts);
    
    // Calls simulation
    // Rcpp::Rcout << "about to enter growthDay_inner_c\n";
    growthDay_inner_c(GROWTHres, GROWTHcomm, x_c,
                      as<std::string>(date[0]),
                      meteovec_c,
                      latitude, elevation, slope, aspect,
                      runon,
                      lateralFlows_c, waterTableDepth);
    //Copies result
    // Rcpp::Rcout << "about to copy results\n";
    l = copyGROWTHResult_c(GROWTHres, x_c);
  } else if(x_c.control.transpirationMode=="Sperry" || x_c.control.transpirationMode=="Sureau"){
    //Initialises a result
    AdvancedTranspiration_RESULT ATres(numCohorts, nlayers, ncanlayers, ntimesteps);
    AdvancedSPWB_RESULT ASPWBres(ATres, nlayers);
    AdvancedGROWTH_RESULT GROWTHres(ASPWBres, numCohorts, ntimesteps);
    // Calls simulation
    growthDay_inner_c(GROWTHres, GROWTHcomm, x_c,
                    as<std::string>(date[0]),
                    meteovec_c,
                    latitude, elevation, slope, aspect,
                    runon,
                    lateralFlows_c, waterTableDepth);
    //Copies result
    l = copyGROWTHResult_c(GROWTHres, x_c);
  } else {
    throw medfate::MedfateInternalError("Wrong transpiration mode");
  }
  
  //Modifies input list, if required  
  if(modifyInput) {
    // Rcpp::Rcout << "about to update input\n";
    x_c.copyStateToList(x);
  }
  
  return(l);
}
