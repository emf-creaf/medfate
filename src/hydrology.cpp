// [[Rcpp::interfaces(r,cpp)]]

#include <RcppArmadillo.h>
#include <cmath>
#include "soil.h"
#include "soil_c.h"
#include "root.h"
#include "biophysicsutils.h"
#include "biophysicsutils_c.h"
#include "hydraulics.h"
#include "hydraulics_c.h"
#include "hydrology_c.h"
#include "numerical_solving.h"
#include <meteoland.h>
using namespace Rcpp;


//' @param month Month of the year (from 1 to 12).
//' @param prec Precipitation for a given day (mm).
//' @param rainfallIntensityPerMonth A vector with twelve positions with average intensity of rainfall (in mm/h) for each month.
//' 
//' @rdname hydrology_interception
//' @keywords internal
// [[Rcpp::export("hydrology_rainfallIntensity")]]
double rainfallIntensity(int month, double prec, NumericVector rainfallIntensityPerMonth){
   // This makes a copy of the vector
   std::vector<double> Ri_month = as<std::vector<double> >(rainfallIntensityPerMonth);
   return(rainfallIntensity_c(month, prec, Ri_month));
}



//' Bare soil evaporation and herbaceous transpiration
//'
//' Functions:
//' \itemize{
//'   \item{Function \code{hydrology_soilEvaporationAmount} calculates the amount of evaporation from bare soil, following Ritchie (1972).}
//'   \item{Function \code{hydrology_soilEvaporation} calculates the amount of evaporation from bare soil and distributes it among soil layers.}
//'   \item{Function \code{hydrology_herbaceousTranspiration} calculates the amount of transpiration due to herbaceous plants.}
//' }
//' 
//' @param soil An object of class \code{\link{soil}}.
//' @param snowpack The amount of snow (in water equivalents, mm) in the snow pack.
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param pet Potential evapotranspiration for a given day (mm)
//' @param LgroundSWR Percentage of short-wave radiation (SWR) reaching the ground.
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
//' 
//' @return 
//' Function \code{hydrology_soilEvaporationAmount} returns the amount of water evaporated from the soil. 
//' 
//' Function \code{hydrology_soilEvaporation} returns a vector of water evaporated from each soil layer.
//' 
//' @references 
//' Ritchie (1972). Model for predicting evaporation from a row crop with incomplete cover. - Water resources research.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_waterInputs}}, \code{\link{hydrology_infiltration}}
//' 
//' 
//' @name hydrology_soilEvaporation
//' @keywords internal
// [[Rcpp::export("hydrology_soilEvaporation")]]
double soilEvaporation(DataFrame soil, double snowpack, 
                       String soilFunctions, double pet, double LgroundSWR,
                       bool modifySoil = true) {
  Soil soil_c = Soil(soil, soilFunctions);

  double Esoil = soilEvaporation_c(soil_c, 
                                   snowpack, pet, LgroundSWR, modifySoil);
  if(modifySoil) {
    NumericVector W = soil["W"];
    for(int l=0;l<W.size();l++) {
      W[l] = soil_c.getW(l);
    }
  }
                                                              
  return(Esoil);
}

//' @rdname hydrology_soilEvaporation
//' @param LherbSWR Percentage of short-wave radiation (SWR) reaching the herbaceous layer.
//' @param herbLAI Leaf area index of the herbaceous layer.
//' @keywords internal
// [[Rcpp::export("hydrology_herbaceousTranspiration")]]
NumericVector herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                                      DataFrame soil, String soilFunctions, bool modifySoil = true){
  NumericVector widths = soil["widths"];
  NumericVector V = ldrRS_one(50, 500, NA_REAL, widths);
  int nlayers = widths.size();
  NumericVector EherbVec(nlayers,0.0);
  Soil soil_c = Soil(soil, soilFunctions);
  std::vector<double> EherbVec_c = as<std::vector<double> >(EherbVec);
  herbaceousTranspiration_c(EherbVec_c, 
                            soil_c, 
                            pet, LherbSWR, herbLAI,
                            as<std::vector<double> >(V),
                            modifySoil);
  
  // Copy values
  for(int l=0;l<nlayers;l++) EherbVec[l] = EherbVec_c[l];
  if(modifySoil) {
    NumericVector W = soil["W"];
    for(int l=0;l<nlayers;l++) {
      W[l] = soil_c.getW(l);
    }
  }
  return(EherbVec);
}










//' @rdname hydrology_infiltration
//' 
//' @param I Soil infiltration (in mm of water).
//' @param widths Width of soil layers (in mm).
//' @param macro Macroporosity of soil layers (in %).
//' @param a,b Parameters of the extinction function used for water infiltration.
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationRepartition")]]
NumericVector infiltrationRepartition(double I, NumericVector widths, NumericVector macro, 
                                      double a = -0.005, double b = 3.0) {

  int nlayers = widths.size();
  NumericVector Ivec(nlayers, 0.0);
  std::vector<double> Ivec_c = as<std::vector<double> >(Ivec);
  std::vector<double> widths_c = as<std::vector<double> >(widths);
  std::vector<double> macro_c = as<std::vector<double> >(macro);
  infiltrationRepartition_c(I, Ivec_c, widths_c, macro_c, a, b);
  for(int i=0;i<nlayers;i++) {
    Ivec[i] = Ivec_c[i];
  }
  return(Ivec);
}


//' Soil infiltration
//'
//' Soil infiltration functions:
//' \itemize{
//'   \item{Function \code{hydrology_infiltrationBoughton} calculates the amount of water that infiltrates into the topsoil, according to the USDA SCS curve number method (Boughton 1989).}
//'   \item{Function \code{hydrology_infiltrationGreenAmpt} calculates the amount of water that infiltrates into the topsoil, according to the model by Green & Ampt (1911).}
//'   \item{Function \code{hydrology_infiltrationAmount} uses either Green & Ampt (1911) or Boughton (1989) to estimate infiltration.}
//'   \item{Function \code{hydrology_infiltrationRepartition} distributes infiltration among soil layers depending on macroporosity.}
//' }
//' 
//' @param rainfallInput Water from the rainfall event reaching the soil surface (mm)
//' @param soil A list containing the description of the soil (see \code{\link{soil}}).
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param rainfallIntensity rainfall intensity rate (mm/h)
//' @param model Infiltration model, either "GreenAmpt1911" or "Boughton1989"
//' @param K_correction Correction for saturated conductivity, to account for increased infiltration due to macropore presence
//' 
//' 
//' @return 
//' Functions \code{hydrology_infiltrationBoughton}, \code{hydrology_infiltrationGreenAmpt} and \code{hydrology_infiltrationAmount} 
//' return the daily amount of water that infiltrates into the soil (in mm of water). 
//' 
//' Function \code{hydrology_infiltrationRepartition} returns the amount of infiltrated water that reaches each soil layer. 
//' 
//' @references 
//' Boughton (1989). A review of the USDA SCS curve number method. - Australian Journal of Soil Research 27: 511-523.
//' 
//' Green, W.H. and Ampt, G.A. (1911) Studies on Soil Physics, 1: The Flow of Air and Water through Soils. The Journal of Agricultural Science, 4, 1-24. 
//' 
//' @details
//'  When using function \code{hydrology_infiltrationGreenAmpt}, the units of \code{Ksat}, \code{t} and \code{psi_wat} have to be in the same system (e.g. cm/h, h and cm). 
//' 
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_waterInputs}}
//' 
//' @name hydrology_infiltration
//' @keywords internal
//' 
//' @rdname hydrology_infiltration
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationAmount")]]
double infiltrationAmount(double rainfallInput, double rainfallIntensity, DataFrame soil, 
                          String soilFunctions, String model = "GreenAmpt1911", double K_correction = 1.0) {
  if(model!="GreenAmpt1911" && model!="Boughton1989") {
    stop("Wrong infiltration model!");
  }
  Soil soil_c = Soil(soil, soilFunctions);
  std::string model_string = model.get_cstring();
  return(infiltrationAmount_c(rainfallInput, rainfallIntensity, 
                              soil_c, 
                              model_string, K_correction));
}



//' Water vertical inputs
//' 
//' High-level functions to define water inputs into the soil of a stand:
//'  
//' \itemize{
//'   \item{Function \code{hydrology_waterInputs} performs canopy water interception and snow accumulation/melt.}
//'   \item{Function \code{hydrology_snowMelt} estimates snow melt using a simple energy balance, according to Kergoat (1998).}
//' }
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' @param prec Precipitation for the given day (mm)
//' @param pet Potential evapotranspiration for the given day (mm)
//' @param rainfallIntensity Rainfall intensity rate (mm/h).
//' @param tday Average day temperature (ºC).
//' @param rad Solar radiation (in MJ/m2/day).
//' @param elevation Altitude above sea level (m).
//' @param Cm Canopy water storage capacity.
//' @param LgroundPAR Percentage of photosynthetically-active radiation (PAR) reaching the ground.
//' @param LgroundSWR Percentage of short-wave radiation (SWR) reaching the ground.
//' @param modifyInput Boolean flag to indicate that the input \code{x} object should be modified during the simulation.
//' 
//' @details 
//' The function simulates different vertical hydrological processes, which are described separately in other functions. 
//' If \code{modifyInput = TRUE} the function will modify the \code{x} object (including both soil moisture and 
//' the snowpack on its surface) as a result of simulating hydrological processes.
//' 
//' @return 
//' Function \code{hydrology_waterInputs} returns a named vector with the following elements, all in mm:
//' \item{Rain}{Precipitation as rainfall.}
//' \item{Snow}{Precipitation as snow.}
//' \item{Interception}{Rainfall water intercepted by the canopy and evaporated.}
//' \item{Snowmelt}{Snow melted during the day, and added to the water infiltrated.}
//' \item{NetRain}{Rainfall reaching the ground.}
//' 
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//'  Kergoat L. (1998). A model for hydrological equilibrium of leaf area index on a global scale. Journal of Hydrology 212–213: 268–286.
//' 
//' @seealso \code{\link{spwb_day}}, \code{\link{hydrology_rainInterception}}, \code{\link{hydrology_soilEvaporation}}
//' 
//' 
//' @name hydrology_verticalInputs
//' @keywords internal
// [[Rcpp::export("hydrology_waterInputs")]]
NumericVector waterInputs(List x,
                          double prec, double rainfallIntensity,
                          double pet, double tday, double rad, double elevation,
                          double Cm, double LgroundPAR, double LgroundSWR, 
                          bool modifyInput = true) {

  WaterInputs_COMM waterInputs;
  ModelInput x_c = ModelInput(x);
  waterInputs_c(waterInputs, x_c,
                 prec, rainfallIntensity,
                 pet, tday, rad, elevation,
                 Cm, LgroundPAR, LgroundSWR, 
                 modifyInput);
  x["snowpack"] = x_c.snowpack;
  NumericVector WI = NumericVector::create(_["Rain"] = waterInputs.rain, _["Snow"] = waterInputs.snow,
                                           _["Interception"] = waterInputs.interception,
                                           _["NetRain"] = waterInputs.netrain, 
                                           _["Snowmelt"] = waterInputs.melt);
  return(WI);
}


//' @name hydrology_verticalInputs
//' @keywords internal
// [[Rcpp::export("hydrology_agricultureWaterInputs")]]
NumericVector agricultureWaterInputs(List x, 
                                     double prec, double tday, double rad, double elevation,
                                     double LgroundSWR, 
                                     bool modifyInput = true) {
  
  WaterInputs_COMM waterInputs;
  AgricultureModelInput x_c = AgricultureModelInput(x);
  agricultureWaterInputs_c(waterInputs, x_c,
                           prec, tday, rad, elevation,
                           LgroundSWR, 
                           modifyInput);
  x["snowpack"] = x_c.snowpack;
  
  NumericVector WI = NumericVector::create(_["Rain"] = waterInputs.rain, _["Snow"] = waterInputs.snow,
                                           _["Interception"] = waterInputs.interception,
                                           _["NetRain"] = waterInputs.netrain, 
                                           _["Snowmelt"] = waterInputs.melt);
  return(WI);
}



//' Soil water balance
//' 
//' Function \code{hydrology_soilWaterBalance} estimates water balance of soil layers given water inputs/outputs, including the simulation of water movement within the soil.
//' 
//' @param soil Object of class \code{\link{soil}}.
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param rainfallInput Amount of water from rainfall event (after excluding interception), in mm.
//' @param rainfallIntensity Rainfall intensity, in mm/h.
//' @param snowmelt Amount of water originated from snow melt, in mm.
//' @param sourceSink Local source/sink term for each soil layer (from soil evaporation or plant transpiration/redistribution)
//'        as mm/day.
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param lateralFlows Lateral source/sink terms for each soil layer (interflow/to from adjacent locations) as mm/day.
//' @param waterTableDepth Water table depth (in mm). When not missing, capillarity rise will be allowed if lower than total soil depth.
//' @param infiltrationMode Infiltration model, either "GreenAmpt1911" or "Boughton1989"
//' @param infiltrationCorrection Correction for saturated conductivity, to account for increased infiltration due to macropore presence
//' @param soilDomains Either "buckets" (multi-bucket domain), "single" (for single-domain Richards) or "dual" (for dual-permeability model).
//' @param nsteps Number of time steps per day
//' @param max_nsubsteps Maximum number of substeps per time step
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_waterInputs}}, \code{\link{hydrology_infiltration}}
//' 
//' @details
//' The multi-bucket model adds/substracts water to each layer and if content is above field capacity the excess percolates to the layer below.
//' If there is still an excess for the bottom layer, the model will progressively fill upper layers (generating saturation excess if the first layer 
//' becomes saturated). Every day the some layers are over field capacity, the model simulates deep drainage.
//' 
//' The single-domain model simulates water flows by solving Richards's equation using the predictor-corrector method, as described in 
//' Bonan et al. (2019).
//' 
//' The dual-permeability model is an implementation of the model MACRO 5.0 (Jarvis et al. 1991; Larsbo et al. 2005).
//' 
//' Both the multi-bucket and the single-domain model can apply a correction to the infiltration rate to account for macroporosity in infiltration. In the
//' dual-permeability model extra infiltration through macropores is simulated explicitly.
//' 
//' @author 
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' \enc{María González Sanchís}{Maria Gonzalez Sanchis}, UPV-CTFC
//' 
//' @return
//'   Returns a named vector with different elements, depending on \code{soilDomains}. 
//'   
//'   If \code{soilDomains == "buckets"}:
//'   \itemize{
//'     \item{\code{Snowmelt}: Snowmelt input (mm).}
//'     \item{\code{Source/sinks}: Sum of source/sink input across layers (mm).}
//'     \item{\code{Infiltration}: Water infiltrated into the soil (mm).}
//'     \item{\code{InfiltrationExcess}: Excess infiltration in the topmost layer (mm) leading to an increase in runoff.}
//'     \item{\code{SaturationExcess}: Excess saturation in the topmost layer (mm) leading to an increase in runoff.}
//'     \item{\code{Runoff}: Surface runoff generated by saturation excess or infiltration excess (mm).}
//'     \item{\code{DeepDrainage}: Water draining from the bottom layer (mm). This quantity is corrected to close the water balance.}
//'     \item{\code{CapillarityRise}: Water entering the soil via capillarity rise (mm) from the water table, if \code{waterTableDepth} is supplied.}
//'   }
//'   If \code{soilDomains == "single"} the named vector contains the following additional elements:
//'   \itemize{
//'     \item{\code{Correction}: Amount of water (mm) added to deep drainage to correct the water balance.}
//'     \item{\code{VolumeChange}: Change in soil water volume (mm).}
//'     \item{\code{Substep}: Time step of the moisture solving (seconds).}
//'   }
//'  If \code{soilDomains == "dual"} the named vector contains the following additional elements:
//'   \itemize{
//'     \item{\code{Lateral flows}: Sum of water circulating between micropores and macropores, positive when filling micropores (mm).}
//'     \item{\code{InfiltrationMatrix}: Water infiltrated into the soil matrix (mm).}
//'     \item{\code{InfiltrationMacropores}: Water infiltrated into the soil macropore domain (mm).}
//'     \item{\code{InfiltrationExcessMatrix/InfiltrationExcessMacropores}: Excess infiltration in the topmost layer (mm) leading to an increase in runoff.}
//'     \item{\code{SaturationExcessMatrix/SaturationExcessMacropores}: Excess saturation in the topmost layer (mm) leading to an increase in runoff.}
//'     \item{\code{DrainageMatrix}: Water draining from the bottom layer of the matrix domain (mm). This quantity is corrected to close water balance in the micropore domain.}
//'     \item{\code{DrainageMacropores}: Water draining from the bottom layer of the macropore domain (mm). This quantity is corrected to close the water balance in the macropore domain.}
//'     \item{\code{CorrectionMatrix}: Amount of water (mm) added to deep drainage of soil matrix to correct the water balance.}
//'     \item{\code{CorrectionMacropores}: Amount of water (mm) added to deep drainage of macropores to correct the water balance.}
//'     \item{\code{MatrixVolumeChange}: Change in soil water volume in the soil matrix domain (mm).}
//'     \item{\code{MacroporeVolumeChange}: Change in soil water volume in the macropore domain (mm).}
//'   }
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK. 
//' 
//' Jarvis, N.J., Jansson, P‐E., Dik, P.E. & Messing, I. (1991). Modelling water and solute transport in macroporous soil. I. Model description and sensitivity analysis. Journal of Soil Science, 42, 59–70. 
//' 
//' Larsbo, M., Roulier, S., Stenemo, F., Kasteel, R. & Jarvis, N. (2005). An Improved Dual‐Permeability Model of Water Flow and Solute Transport in the Vadose Zone. Vadose Zone Journal, 4, 398–406. 
//' 
//' @examples
//' # Define soil parameters
//' spar <- defaultSoilParams(4)
//' 
//' # Initializes soil hydraulic parameters
//' examplesoil <- soil(spar)
//' 
//' # Water balance in a multi-bucket model
//' hydrology_soilWaterBalance(examplesoil, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
//'                            soilDomains = "buckets", modifySoil = FALSE)
//'                            
//' # Water balance in a single-domain model (Richards equation)
//' hydrology_soilWaterBalance(examplesoil, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
//'                            soilDomains = "single", modifySoil = FALSE)
//'                     
//' # Water balance in a dual-permeability model (MACRO)
//' hydrology_soilWaterBalance(examplesoil, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
//'                            soilDomains = "dual", modifySoil = FALSE)
//'   
//' @name hydrology_soilWaterBalance
//' @keywords internal
// [[Rcpp::export("hydrology_soilWaterBalance")]]
NumericVector soilWaterBalance(DataFrame soil, String soilFunctions, 
                               double rainfallInput, double rainfallIntensity, double snowmelt, NumericVector sourceSink, 
                               double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                               String infiltrationMode = "GreenAmpt1911", double infiltrationCorrection = 5.0, String soilDomains = "buckets", 
                               int nsteps = 24, int max_nsubsteps = 3600, bool modifySoil = true) {
  std::string infiltrationMode_str = infiltrationMode.get_cstring();
  std::string soilDomains_str = soilDomains.get_cstring();
  SoilWaterBalance_COMM SWBcomm = SoilWaterBalance_COMM(soil.nrow());
  SoilWaterBalance_RESULT SWBres;
  Soil soil_c = Soil(soil, soilFunctions);
  std::vector<double> lateralFlows_c(soil.nrow(), 0.0);
  NumericVector lateralFlows_mm;
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<lateralFlows_mm.size();l++) {
      lateralFlows_c[l] = lateralFlows_mm[l];
    }
  }
  std::vector<double> sourceSink_c= Rcpp::as< std::vector<double> >(sourceSink);
  soilWaterBalance_inner_c(SWBres, SWBcomm, soil_c,
                           rainfallInput, rainfallIntensity, snowmelt, sourceSink_c,
                           runon, lateralFlows_c, waterTableDepth,
                           infiltrationMode_str, infiltrationCorrection,
                           soilDomains_str,
                           nsteps, max_nsubsteps);
  NumericVector res;
  if(soilDomains== "buckets") {
    res = NumericVector::create(_["Local source/sinks"] = SWBres.localSourceSinks_mm,
                                _["Lateral source/sinks"] = SWBres.lateralSourceSinks_mm,
                                _["Infiltration"] = SWBres.infiltration_mm,
                                _["InfiltrationExcess"] = SWBres.infiltrationExcess_mm,
                                _["SaturationExcess"] = SWBres.saturationExcess_mm,
                                _["Runoff"] = SWBres.runoff_mm,
                                _["DeepDrainage"] = SWBres.deepDrainage_mm);
  } else if(soilDomains == "single") {
    res = NumericVector::create(_["Local source/sinks"] = SWBres.localSourceSinks_mm,
                                _["Lateral source/sinks"] = SWBres.lateralSourceSinks_mm,
                                _["Infiltration"] = SWBres.infiltration_mm,
                                _["InfiltrationExcess"] = SWBres.infiltrationExcess_mm,
                                _["SaturationExcess"] = SWBres.saturationExcess_mm,
                                _["Runoff"] = SWBres.runoff_mm,
                                _["DeepDrainage"] = SWBres.deepDrainage_mm,
                                _["CapillarityRise"] = SWBres.capillarityRise_mm,
                                _["Correction"] = SWBres.correction_mm,
                                _["VolumeChange"] = SWBres.volumeChange_mm,
                                _["Substeps"] = SWBres.substeps);
  } else if(soilDomains == "dual") {
    res = NumericVector::create(_["Local source/sinks"] = SWBres.localSourceSinks_mm,
                                _["Lateral source/sinks"] = SWBres.lateralSourceSinks_mm,
                                _["Matrix-macropore flow"] = SWBres.matrixMacroporeFlow_mm,
                                _["InfiltrationMatrix"] = SWBres.infiltrationMatrix_mm,
                                _["InfiltrationMacropores"] = SWBres.infiltrationMacropores_mm,
                                _["InfiltrationExcessMatrix"] = SWBres.infiltrationExcessMatrix_mm,
                                _["InfiltrationExcessMacropores"] = SWBres.infiltrationExcessMacropores_mm,
                                _["SaturationExcessMatrix"] = SWBres.saturationExcessMatrix_mm,
                                _["SaturationExcessMacropores"] = SWBres.saturationExcessMacropores_mm,
                                _["DrainageMatrix"] = SWBres.drainageMatrix_mm,
                                _["DrainageMacropores"] = SWBres.drainageMacropores_mm,
                                _["CapillarityMatrix"] = SWBres.capillarityMatrix_mm,
                                _["CapillarityMacropores"] = SWBres.capillarityMacropores_mm,
                                _["CorrectionMatrix"] = SWBres.correctionMatrix_mm,
                                _["CorrectionMacropores"] = SWBres.correctionMacropores_mm,
                                _["MatrixVolumeChange"] = SWBres.matrixVolumeChange_mm,
                                _["MacroporeVolumeChange"] = SWBres.macroporeVolumeChange_mm,
                                _["Infiltration"] = SWBres.infiltration_mm,
                                _["InfiltrationExcess"] = SWBres.infiltrationExcess_mm);
    res.push_back(SWBres.saturationExcess_mm, "SaturationExcess");
    res.push_back(SWBres.runoff_mm, "Runoff");
    res.push_back(SWBres.deepDrainage_mm, "DeepDrainage");
    res.push_back(SWBres.capillarityRise_mm, "CapillarityRise");
    res.push_back(SWBres.correction_mm, "Correction");
    res.push_back(SWBres.volumeChange_mm, "VolumeChange");
    res.push_back((double) SWBres.substeps, "Substeps");
  }
  if(modifySoil) {
    NumericVector W = soil["W"];
    for(int l=0;l<W.size();l++) {
      W[l] = soil_c.getW(l);
    }
  }
  return(res);
}
