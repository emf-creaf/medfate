// [[Rcpp::interfaces(r,cpp)]]

#include <numeric>
#include <math.h>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "photosynthesis.h"
#include "phenology.h"
#include "transpiration.h"
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;


// Soil water balance with simple hydraulic model
List spwbDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
             double rad = NA_REAL, double elevation = NA_REAL, bool verbose = false) {

  //Control parameters
  List control = x["control"];
  bool snowpack = control["snowpack"];
  bool drainage = control["drainage"];
  String soilFunctions = control["soilFunctions"];

  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  int nlayers = W.size();
  

  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIphe.size();


  //Parameters  
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsBase["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsBase["g"]);
  

  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector Phe(numCohorts,0.0);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0, Cm = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(LAIlive[c]>0) Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    else Phe[c]=0.0;
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  double LgroundPAR = exp((-1.0)*s);
  double LgroundSWR = exp((-1.0)*s/1.35);
  
  //Snow pack dynamics and hydrology input
  NumericVector hydroInputs = verticalInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                             Cm, LgroundPAR, LgroundSWR, 
                                             runon,
                                             snowpack, drainage, true);

  //Evaporation from bare soil if there is no snow
  NumericVector EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);

  //Canopy transpiration  
  List transp = transpirationGranier(x, soil, tday, pet, true);
  NumericVector EplantVec(nlayers, 0.0);
  NumericMatrix EplantCoh = Rcpp::as<Rcpp::NumericMatrix>(transp["Extraction"]);
  for(int l=0;l<nlayers;l++) EplantVec[l] = sum(EplantCoh(_,l));
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);
  
  NumericVector psiVec = psi(soil, soilFunctions); //Calculate current soil water potential for output
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, _["Rain"] = hydroInputs["Rain"], _["Snow"] = hydroInputs["Snow"], 
                                           _["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], 
                                           _["Infiltration"] = hydroInputs["Infiltration"], _["Runoff"] = hydroInputs["Runoff"], _["DeepDrainage"] = hydroInputs["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(EplantVec),
                                           _["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, 
                                           _["Cm"] = Cm, _["Lground"] = LgroundPAR);
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = DB, 
                        _["Soil"] = SB,
                        _["Plants"] = Plants);
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}



// Soil water balance with Sperry hydraulic and stomatal conductance models
List spwbDay2(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double slope, double aspect,
             double solarConstant, double delta, 
             double prec, double pet, double er, double runon=0.0, bool verbose = false) {
  
  //Control parameters
  List control = x["control"];
  bool drainage = control["drainage"];
  bool snowpack = control["snowpack"];
  String soilFunctions = control["soilFunctions"];
  int ntimesteps = control["ndailysteps"];

  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  int nlayers = W.size();

  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIlive.size();

  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsBase["Sgdd"]);
  NumericVector albedo = Rcpp::as<Rcpp::NumericVector>(paramsBase["albedo"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsBase["g"]);

  //1. Leaf Phenology: Adjusted leaf area index
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  NumericVector Phe(numCohorts);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0, Cm = 0.0, LAIcellmax = 0.0;
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcellmax += LAIlive[c];
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
  }
  double LgroundPAR = exp((-1.0)*s);
  double LgroundSWR = exp((-1.0)*s/1.35);
  
  //Snow pack dynamics and hydrology input
  NumericVector hydroInputs = verticalInputs(soil, soilFunctions, prec, er, tday, rad, elevation,
                                             Cm, LgroundPAR, LgroundSWR, 
                                             runon,
                                             snowpack, drainage);
  
  
  
  //B.1 - Evaporation from bare soil if there is no snow
  NumericVector EsoilVec = soilEvaporation(soil, soilFunctions, pet, LgroundSWR, true);
  
  //B.2 - Canopy transpiration  
  List transp = transpirationSperry(x, soil,tmin, tmax, rhmin, rhmax, rad, wind, 
                                    latitude, elevation, slope, aspect, 
                                    solarConstant, delta, prec, 
                                    hydroInputs["Interception"], hydroInputs["Snowmelt"], sum(EsoilVec),
                                    verbose, NA_INTEGER, true);

  
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(transp["ExtractionInst"]);
  NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(transp["RhizoPsi"]);
  DataFrame Plants = Rcpp::as<Rcpp::DataFrame>(transp["Plants"]);
  NumericVector Eplant = Plants["Transpiration"];
  
  List PlantsInst = Rcpp::as<Rcpp::List>(transp["PlantsInst"]);
  List EnergyBalance = Rcpp::as<Rcpp::List>(transp["EnergyBalance"]);
  
  //B.3 - Substract evaporated and extracted water from soil moisture 
  NumericVector EplantVec(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  for(int l=0;l<nlayers;l++) {
    for(int n=0;n<ntimesteps;n++) {
      soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
      soilHydraulicOutput[l] += std::max(soilLayerExtractInst(l,n),0.0);
    }
    EplantVec[l] = sum(soilLayerExtractInst(l,_));
  }
  NumericVector psiVec = psi(soil, soilFunctions); //Calculate current soil water potential for output
  
  NumericVector DB = NumericVector::create(_["PET"] = pet,_["Rain"] = hydroInputs["Rain"],_["Snow"] = hydroInputs["Snow"],_["NetRain"] = hydroInputs["NetRain"], _["Snowmelt"] = hydroInputs["Snowmelt"],
                                           _["Runon"] = hydroInputs["Runon"], _["Infiltration"] = hydroInputs["Infiltration"], _["Runoff"] = hydroInputs["Runoff"], _["DeepDrainage"] = hydroInputs["DeepDrainage"],
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(Eplant),
                                           _["HydraulicRedistribution"] = sum(soilHydraulicInput),
                                           _["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, _["Cm"] = Cm, _["Lground"] = LgroundPAR);
  
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = DB, 
                        _["EnergyBalance"] = EnergyBalance,
                        _["Soil"] = SB, 
                        _["RhizoPsi"] = RhizoPsi,
                        _["Plants"] = Plants,
                        _["SunlitLeaves"] = transp["SunlitLeaves"],
                        _["ShadeLeaves"] = transp["ShadeLeaves"],
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["PlantsInst"] = PlantsInst,
                        _["LightExtinction"] = transp["LightExtinction"],
                        _["WindExtinction"] = transp["WindExtinction"]);
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

// [[Rcpp::export("spwb_day")]]
List spwbDay(List x, List soil, CharacterVector date, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
            double latitude, double elevation, double slope, double aspect,  
            double prec, double runon=0.0) {
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
  bool leafPhenology = control["leafPhenology"];
  String transpirationMode = control["transpirationMode"];
  std::string c = as<std::string>(date[0]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  double tday = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  double latrad = latitude * (PI/180.0);
  double asprad = aspect * (PI/180.0);
  double slorad = slope * (PI/180.0);
  double pet = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);

  //Derive doy from date  
  int J0101 = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),1,1);
  int doy = J - J0101+1;
  if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
  if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
  
  //Update phenology
  if(leafPhenology) updateLeaves(x, doy, tday, wind);
  
  double er = erFactor(doy, pet, prec);
  List s;
  if(transpirationMode=="Granier") {
    s = spwbDay1(x,soil, tday, pet, prec, er, runon, rad, elevation, verbose);
  } else {
    s = spwbDay2(x,soil, tmin, tmax, rhmin, rhmax, rad, wind, 
                 latitude, elevation, slope, aspect,
                 solarConstant, delta, prec, pet, er, runon, verbose);
  }
  // Rcout<<"hola4\n";
  return(s);
}

  

IntegerVector order_vector(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}



void checkspwbInput(List x, List soil, String transpirationMode, String soilFunctions) {
  if(!x.containsElementNamed("above")) stop("above missing in spwbInput");
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  if(!above.containsElementNamed("LAI_live")) stop("LAI_live missing in spwbInput$above");
  if(!above.containsElementNamed("CR")) stop("CR missing in spwbInput$above");
  if(!above.containsElementNamed("H")) stop("H missing in spwbInput$above");
  
  if(!x.containsElementNamed("below")) stop("below missing in spwbInput");
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  if(!below.containsElementNamed("V")) stop("V missing in spwbInput$below");
  if(transpirationMode=="Sperry"){
    if(!below.containsElementNamed("VGrhizo_kmax")) stop("VGrhizo_kmax missing in spwbInput$below");
    if(!below.containsElementNamed("VCroot_kmax")) stop("VCroot_kmax missing in spwbInput$below");
  }  
  
  if(!x.containsElementNamed("paramsBase")) stop("paramsBase missing in spwbInput");
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  if(!paramsBase.containsElementNamed("Sgdd")) stop("Sgdd missing in spwbInput$paramsBase");
  if(!paramsBase.containsElementNamed("k")) stop("k missing in spwbInput$paramsBase");
  if(!paramsBase.containsElementNamed("g")) stop("g missing in spwbInput$paramsBase");
  
  if(!x.containsElementNamed("paramsTransp")) stop("paramsTransp missing in spwbInput");
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  if(!paramsTransp.containsElementNamed("pRootDisc")) stop("pRootDisc missing in spwbInput$paramsTransp");
  if(transpirationMode=="Granier") {
    if(!paramsTransp.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("WUE")) stop("WUE missing in spwbInput$paramsTransp");
  } else if(transpirationMode=="Sperry") {
    if(!paramsTransp.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCstem_c")) stop("VCstem_c missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCstem_d")) stop("VCstem_d missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_c")) stop("VCroot_c missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_d")) stop("VCroot_d missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Gwmax")) stop("Gwmax missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Vmax298")) stop("Vmax298 missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Jmax298")) stop("Jmax298 missing in spwbInput$paramsTransp");
  }
  if(transpirationMode=="Sperry") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
  }
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(soilFunctions=="SX") {
    if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
    if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
  }
  if(soilFunctions=="VG") {
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
    if(!soil.containsElementNamed("VG_theta_res")) stop("VG_theta_res missing in soil");
    if(!soil.containsElementNamed("VG_theta_sat")) stop("VG_theta_sat missing in soil");
  }
}

// [[Rcpp::export("spwb_resetInputs")]]
void resetInputs(List x, List soil, List from = R_NilValue, int day = NA_INTEGER) {
  List can = x["canopy"];
  NumericVector W = soil["W"];
  int nlayers = W.size();
  NumericVector Temp = soil["Temp"];
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  
  if(Rf_isNull(from) || from.size()==0) {
    can["gdd"] = 0.0;
    can["Temp"] = NA_REAL;
    for(int i=0;i<nlayers;i++) {
      W[i] = 1.0; //Defaults to soil at field capacity
      Temp[i] = NA_REAL;
    }
    if(transpirationMode=="Sperry") {
      NumericVector psiRoot = Rcpp::as<Rcpp::NumericVector>(x["psiRoot"]);
      NumericMatrix psiStem = Rcpp::as<Rcpp::NumericMatrix>(x["psiStem"]);
      NumericMatrix PLCstem = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
      NumericMatrix RWCsympstem = Rcpp::as<Rcpp::NumericMatrix>(x["RWCsympstem"]);
      NumericVector RWCsympleaf = Rcpp::as<Rcpp::NumericVector>(x["RWCsympleaf"]);
      NumericVector psiLeaf = Rcpp::as<Rcpp::NumericVector>(x["psiLeaf"]);
      NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(x["Einst"]);
      NumericVector Transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
      NumericVector Photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
      for(int i=0;i<PLCstem.nrow();i++) {
        Einst[i] = 0.0;
        psiLeaf[i] = 0.0;
        psiRoot[i] = 0.0;
        RWCsympleaf[i] = 0.0;
        Transpiration[i] = 0.0;
        Photosynthesis[i] = 0.0;
        for(int j=0;j<PLCstem.ncol();j++) {
          psiStem(i,j) = 0.0;
          PLCstem(i,j) = 0.0; 
          RWCsympstem(i,j) = 1.0; 
        }
      }
    } else {
      NumericVector Transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
      NumericVector Photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
      NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(x["PLC"]);
      for(int i=0;i<Transpiration.length();i++) {
        Transpiration[i] = 0.0;
        Photosynthesis[i] = 0.0;
        PLC[i] = 0.0;
      }
    }

  } else {
    if(IntegerVector::is_na(day)) day = 0;
    else day = day-1; //Input will be 1 for first day
    DataFrame DWB = Rcpp::as<Rcpp::DataFrame>(from["WaterBalance"]);
    DataFrame SWB = Rcpp::as<Rcpp::DataFrame>(from["Soil"]);
    NumericVector GDD = DWB["GDD"];
    can["gdd"] = GDD[day];
    can["Temp"] = NA_REAL;
    for(int i=0;i<nlayers;i++) {
      W[i] = Rcpp::as<Rcpp::NumericVector>(SWB[i])[day];
      //TO DO: STORE/RECOVER SOIL LAYER TEMPERATURE?
      Temp[i] = NA_REAL;
    }
    NumericMatrix fromPLC = Rcpp::as<Rcpp::NumericMatrix>(from["PlantStress"]);
    NumericMatrix fromRootPsi = Rcpp::as<Rcpp::NumericMatrix>(from["RootPsi"]);
    NumericMatrix fromLeafPsiMin = Rcpp::as<Rcpp::NumericMatrix>(from["LeafPsiMin"]);
    NumericMatrix fromStemPsi = Rcpp::as<Rcpp::NumericMatrix>(from["StemPsi"]);
    NumericMatrix fromRWCstem = Rcpp::as<Rcpp::NumericMatrix>(from["StemRWC"]);
    NumericMatrix fromRWCleaf = Rcpp::as<Rcpp::NumericMatrix>(from["LeafRWC"]);
    
    NumericVector psiRoot = Rcpp::as<Rcpp::NumericVector>(x["psiRoot"]);
    NumericMatrix psiStem = Rcpp::as<Rcpp::NumericMatrix>(x["psiStem"]);
    NumericVector psiLeaf = Rcpp::as<Rcpp::NumericVector>(x["psiLeaf"]);
    NumericMatrix PLCstem = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
    NumericMatrix RWCsympstem = Rcpp::as<Rcpp::NumericMatrix>(x["RWCsympstem"]);
    NumericVector RWCsympleaf = Rcpp::as<Rcpp::NumericVector>(x["RWCsympleaf"]);
    NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(x["Einst"]);
    NumericVector Transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
    NumericVector Photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
    for(int i=0;i<PLCstem.nrow();i++) {
      Einst[i] = 0.0;
      Transpiration[i] = 0.0;
      Photosynthesis[i] = 0.0;
      psiRoot[i] = fromRootPsi(day,i);
      psiLeaf[i] = fromLeafPsiMin(day,i);
      RWCsympleaf[i] = fromRWCleaf(day,i);
      for(int j=0;j<PLCstem.ncol();j++) {
        psiStem(i,j) = fromStemPsi(day,i);
        PLCstem(i,j) = fromPLC(day,i); 
        RWCsympstem(i,j) = fromRWCstem(day,i); 
      }
    }
  }
  soil["W"] = W;
  soil["Temp"] =Temp;
}

// [[Rcpp::export("spwb")]]
List spwb(List x, List soil, DataFrame meteo, double latitude = NA_REAL, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  checkspwbInput(x, soil, transpirationMode, soilFunctions);
  
  //Store input
  List spwbInput = clone(x);
  List soilInput = clone(soil);
    
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector PET = NumericVector(numDays,0.0);
  if(transpirationMode=="Granier") {
    if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
    PET = meteo["PET"];
    if(control["snowpack"]) {
      if(!meteo.containsElementNamed("Radiation")) stop("If 'snowpack = TRUE', variable 'Radiation' must be provided.");
      else Radiation = meteo["Radiation"];
    }
  } else if(transpirationMode=="Sperry") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
    MinTemperature = meteo["MinTemperature"];
    if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
    MaxTemperature = meteo["MaxTemperature"];
    if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
    MinRelativeHumidity = meteo["MinRelativeHumidity"];
    if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
    MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
    Radiation = meteo["Radiation"];
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  

  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  

  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Detailed subday results
  List subdailyRes(numDays);
  

  //Water balance output variables
  NumericVector GDD(numDays); 
  NumericVector LAIcell(numDays),LAIcelldead(numDays);
  NumericVector Cm(numDays);
  NumericVector Lground(numDays);
  NumericVector Runoff(numDays);
  NumericVector Rain(numDays);
  NumericVector Snow(numDays);
  NumericVector Snowmelt(numDays);
  NumericVector NetRain(numDays);
  NumericVector Interception(numDays);
  NumericVector Infiltration(numDays);
  NumericVector DeepDrainage(numDays);
  NumericVector SoilEvaporation(numDays);
  NumericVector Transpiration(numDays);
  NumericVector PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  
  NumericMatrix Eplantdays(numDays, nlayers);
  NumericMatrix HydrIndays(numDays, nlayers);
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix MLdays(numDays, nlayers);
  NumericVector WaterTable(numDays, NA_REAL);
  NumericVector MLTot(numDays, 0.0);
  NumericVector SWE(numDays, 0.0);
  
  //EnergyBalance output variables
  NumericVector Tatm_mean(numDays, NA_REAL);
  NumericVector Tatm_min(numDays, NA_REAL);
  NumericVector Tatm_max(numDays, NA_REAL);
  NumericVector Tcan_mean(numDays, NA_REAL);
  NumericVector Tcan_min(numDays, NA_REAL);
  NumericVector Tcan_max(numDays, NA_REAL);
  NumericVector Tsoil_mean(numDays, NA_REAL);
  NumericVector Tsoil_min(numDays, NA_REAL);
  NumericVector Tsoil_max(numDays, NA_REAL);
  NumericVector SWRcanin(numDays, NA_REAL);
  NumericVector LWRcanin(numDays, NA_REAL);
  NumericVector LWRcanout(numDays, NA_REAL);
  NumericVector LEcan_heat(numDays, NA_REAL);
  NumericVector LEsoil_heat(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilout(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector LWRsoilcan(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);
  NumericVector RAcan(numDays, NA_REAL);
  NumericVector RAsoil(numDays, NA_REAL);
  
  //Plant output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix dEdP(numDays, numCohorts);
  NumericMatrix LeafPsiMin(numDays, numCohorts);
  NumericMatrix LeafPsiMax(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SH(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SH(numDays, numCohorts);
  NumericMatrix StemPsi(numDays, numCohorts);
  NumericMatrix RootPsi(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  List RhizoPsi(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(meteo.attr("row.names"), seq(1,nlayers)) ;
    RhizoPsi[c] = nm;
  }
  RhizoPsi.attr("names") = above.attr("row.names");
  
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts);
  NumericMatrix LeafRWC(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  NumericMatrix PlantAbsSWR(numDays, numCohorts);
  NumericMatrix PlantAbsLWR(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  
  
  NumericVector Wini = soil["W"];
  Wdays(0,_) = Wini;
  NumericVector initialContent = water(soil, soilFunctions);
  if(verbose) {
    Rcout<<"Initial soil water content (mm): "<< sum(initialContent)<<"\n";
  }

  if(verbose) Rcout << "Performing daily simulations ";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  for(int i=0;i<numDays;i++) {
      if(verbose & (i%10 == 0)) Rcout<<".";//<<i;
      
      double wind = WindSpeed[i];
      if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
      if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
      
      //If DOY == 1 reset PLC (Growth assumed)
      if(DOY[i]==1) {
        if(transpirationMode=="Granier") {
          NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(x["PLC"]);
          for(int j=0;j<PLC.length();j++) PLC[j] = 0.0;
        } else {
          NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
          for(int j=0;j<StemPLC.nrow();j++) for(int k=0;k<StemPLC.ncol();k++)  StemPLC(j,k) = 0.0;
        }
      }

      
      //1. Phenology and leaf fall
      if(leafPhenology) updateLeaves(x, DOY[i], MeanTemperature[i], wind);
      
      //Store GDD
      GDD[i] = canopyParams["gdd"];
      

      
      //2. Water balance and photosynthesis
      if(transpirationMode=="Granier") {
        double er = erFactor(DOY[i], PET[i], Precipitation[i]);
        s = spwbDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], er, 0.0, 
                     Radiation[i], elevation, verbose); //No Runon in simulations for a single cell
      } else if(transpirationMode=="Sperry") {
        int ntimesteps = control["ndailysteps"];
        double tstep = 86400.0/((double) ntimesteps);
        std::string c = as<std::string>(dateStrings[i]);
        int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
        double delta = meteoland::radiation_solarDeclination(J);
        double solarConstant = meteoland::radiation_solarConstant(J);
        double latrad = latitude * (PI/180.0);
        if(NumericVector::is_na(aspect)) aspect = 0.0;
        if(NumericVector::is_na(slope)) slope = 0.0;
        double asprad = aspect * (PI/180.0);
        double slorad = slope * (PI/180.0);
        double tmin = MinTemperature[i];
        double tmax = MaxTemperature[i];
        double rhmin = MinRelativeHumidity[i];
        double rhmax = MaxRelativeHumidity[i];
        double rad = Radiation[i];
        PET[i] = meteoland::penman(latrad, elevation, slorad, asprad, J, tmin, tmax, rhmin, rhmax, rad, wind);
        double er = erFactor(DOY[i], PET[i], Precipitation[i]);
        s = spwbDay2(x, soil, tmin, tmax, 
                     rhmin, rhmax, rad, wind, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, Precipitation[i], PET[i], 
                     er, 0.0, verbose);
        List Plants = Rcpp::as<Rcpp::List>(s["Plants"]);
        List PlantsInst = Rcpp::as<Rcpp::List>(s["PlantsInst"]);
        List ShadeLeaves = Rcpp::as<Rcpp::List>(PlantsInst["ShadeLeaves"]);
        List SunlitLeaves = Rcpp::as<Rcpp::List>(PlantsInst["SunlitLeaves"]);
        
        NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_SWR"]);
        NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_SWR"]);
        NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_LWR"]);
        NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_LWR"]);
        for(int j=0;j<numCohorts;j++) {
          for(int n=0;n<ntimesteps;n++){
            PlantAbsSWR(i,j) += 0.000001*(SWR_SL(j,n)+SWR_SH(j,n))*tstep;
            PlantAbsLWR(i,j) += 0.000001*(LWR_SL(j,n)+LWR_SH(j,n))*tstep;
          }
        }
        List EB = Rcpp::as<Rcpp::List>(s["EnergyBalance"]);
        DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
        DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
        DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
        NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
        NumericVector Tsoil = Rcpp::as<Rcpp::NumericVector>(Tinst["Tsoil.1"]);

        SWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcanin"]))*tstep;
        LWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanin"]))*tstep;
        LWRcanout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanout"]))*tstep;
        LEcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEcan"]))*tstep;
        Hcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
        Ebalcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
        LWRsoilcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRsoilcan"]))*tstep;
        RAcan[i] = sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["RAcan"]))/((double) ntimesteps);
        SWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoilin"]))*tstep;
        LWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilin"]))*tstep;
        LWRsoilout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilout"]))*tstep;
        LEsoil_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEsoil"]))*tstep;
        Hcansoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
        Ebalsoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;
        RAsoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["RAsoil"]))/((double) ntimesteps);
        
        Tatm_min[i] = MinTemperature[i];
        Tatm_max[i] = MaxTemperature[i];
        Tatm_mean[i] = MeanTemperature[i];
        Tcan_min[i] = min(Tcan);
        Tcan_max[i] = max(Tcan);
        Tcan_mean[i] = sum(Tcan)/((double) ntimesteps);
        Tsoil_min[i] = min(Tsoil);
        Tsoil_max[i] = max(Tsoil);
        Tsoil_mean[i] = sum(Tsoil)/((double) ntimesteps);
      }
      List db = s["WaterBalance"];
      Lground[i] = db["Lground"];
      LAIcell[i] = db["LAIcell"];
      LAIcelldead[i] = db["LAIcelldead"];
      Cm[i] = db["Cm"];
      DeepDrainage[i] = db["DeepDrainage"];
      Infiltration[i] = db["Infiltration"];
      Runoff[i] = db["Runoff"];
      Rain[i] = db["Rain"];
      Snow[i] = db["Snow"];
      Snowmelt[i] = db["Snowmelt"];
      NetRain[i] = db["NetRain"];
      PlantExtraction[i] = db["PlantExtraction"];
      if(transpirationMode=="Sperry")  {
        HydraulicRedistribution[i] = db["HydraulicRedistribution"];
      }
      Transpiration[i] = db["Transpiration"];
      SoilEvaporation[i] = db["SoilEvaporation"];
      Interception[i] = Rain[i]-NetRain[i];
      
      List sb = s["Soil"];
      NumericVector psi = sb["psi"];
      psidays(i,_) = psi;
      NumericVector EplantVec = sb["PlantExtraction"];
      SWE[i] = soil["SWE"];
        
      List Plants = s["Plants"];
      PlantLAI(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
      NumericVector EplantCoh = Plants["Transpiration"];
      Eplantdays(i,_) = EplantVec;
      PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["Photosynthesis"]);
      PlantTranspiration(i,_) = EplantCoh;
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
      if(transpirationMode=="Sperry") {
        NumericVector HydrInVec = sb["HydraulicInput"];
        HydrIndays(i,_) = HydrInVec;
        LeafPsiMin(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
        LeafPsiMax(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
        LeafPsiMin_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SL"]);
        LeafPsiMax_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SL"]);
        LeafPsiMin_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SH"]);
        LeafPsiMax_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SH"]);
        RootPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["RootPsi"]); 
        StemPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPsi"]); 
        StemPLC(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
        dEdP(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["dEdP"]); 
        NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(s["RhizoPsi"]);
        for(int c=0;c<numCohorts;c++) {
          NumericMatrix nm = Rcpp::as<Rcpp::NumericMatrix>(RhizoPsi[c]);
          nm(i,_) =  RhizoPsiStep(c,_);
        }
      } else {
        PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);
      }
      EplantCohTot = EplantCohTot + EplantCoh;
      Eplanttot[i] = sum(EplantCoh);
      if(transpirationMode=="Sperry"){
        StemRWC(i,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
        LeafRWC(i,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
      }
      
      if(subdailyResults) {
        subdailyRes[i] = clone(s);
      }
      if(i<(numDays-1)) Wdays(i+1,_) = as<Rcpp::NumericVector>(soil["W"]);
      WaterTable[i] = waterTableDepth(soil, soilFunctions);
  }
  if(verbose) Rcout << "done.\n";
  
  for(int l=0;l<nlayers;l++) {
    MLdays(_,l) = Wdays(_,l)*Water_FC[l]; 
    MLTot = MLTot + MLdays(_,l);
  }
  NumericVector Evapotranspiration = Transpiration+SoilEvaporation + Interception;
  
  if(verbose) {
    NumericVector finalContent = water(soil, soilFunctions);
    Rcout<<"Final soil water content (mm): "<< sum(finalContent)<<"\n";
    Rcout<<"Change in soil water content (mm): "<< sum(finalContent) - sum(initialContent)<<"\n";

    double Precipitationsum = sum(Precipitation);
    double NetRainsum = sum(NetRain);
    double Interceptionsum = sum(Interception);
    double SoilEvaporationsum = sum(SoilEvaporation);
    double Runoffsum  = sum(Runoff);
    double Infiltrationsum  = sum(Infiltration);
    double DeepDrainagesum = sum(DeepDrainage);
    double Transpirationsum = sum(Transpiration);
    
    double wb = Precipitationsum - Interceptionsum - Runoffsum - DeepDrainagesum - SoilEvaporationsum - sum(PlantExtraction);
    Rcout<<"Water balance result (mm): "<< wb<<"\n";
    Rcout<<"Water balance components:\n";
    Rcout<<"  Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    Rcout<<"  Rain (mm) "  <<round(sum(Rain)) <<" Snow (mm) "  <<round(sum(Snow)) <<"\n";
    Rcout<<"  Interception (mm) " << round(Interceptionsum)  <<" Net rainfall (mm) " << round(NetRainsum) <<"\n";
    Rcout<<"  Infiltration (mm) " << round(Infiltrationsum)  <<
      " Runoff (mm) " << round(Runoffsum) <<
        " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
    Rcout<<"  Soil evaporation (mm) " << round(SoilEvaporationsum);
    Rcout<<" Transpiration (mm) "  <<round(Transpirationsum) <<"\n";
    if(transpirationMode =="Sperry") {
      Rcout<<"  Plant extraction from soil (mm) " << round(sum(PlantExtraction));
      Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    }
  }

   DataFrame SWB;
   if(transpirationMode=="Granier") {
     SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                             _["WTD"] = WaterTable,
                             _["SWE"] = SWE, 
                             _["PlantExt"]=Eplantdays,
                             _["psi"]=psidays); 
   } else {
     SWB = DataFrame::create(_["W"]=Wdays, _["ML"]=MLdays,_["MLTot"]=MLTot,
                             _["WTD"] = WaterTable,
                             _["SWE"] = SWE, 
                             _["PlantExt"]=Eplantdays,
                             _["HydraulicInput"] = HydrIndays,
                             _["psi"]=psidays); 
   }
   SWB.attr("row.names") = meteo.attr("row.names") ;
   DataFrame DWB;
   if(transpirationMode=="Granier") {
     DWB = DataFrame::create(_["GDD"] = GDD,
                             _["LAIcell"]=LAIcell, _["LAIcelldead"] = LAIcelldead,  _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                             _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration);
   } else {
     DWB = DataFrame::create(_["GDD"] = GDD,
                             _["LAIcell"]=LAIcell, _["LAIcelldead"] = LAIcelldead,  _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow, 
                             _["NetRain"]=NetRain, _["Snowmelt"] = Snowmelt, _["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                             _["HydraulicRedistribution"] = HydraulicRedistribution);
     
   }
   DWB.attr("row.names") = meteo.attr("row.names") ;

  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPLC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  dEdP.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  RootPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsSWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsLWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantLAI.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  
  DataFrame DEB = DataFrame::create(_["SWRcanin"] = SWRcanin, _["LWRcanin"] = LWRcanin, _["LWRcanout"] = LWRcanout,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebalcan, 
                                    _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = SWRsoilin, _["LWRsoilin"] = LWRsoilin, _["LWRsoilout"] = LWRsoilout,
                                    _["LEsoil"] = LEsoil_heat, _["Ebalsoil"] = Ebalsoil, _["RAcan"] = RAcan, _["RAsoil"] = RAsoil);  
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                   _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DEB.attr("row.names") = meteo.attr("row.names") ;
  DT.attr("row.names") = meteo.attr("row.names") ;
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");
  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("spwbInput") = spwbInput,
                     Named("soilInput") = soilInput,
                     Named("WaterBalance")=DWB, 
                     Named("Soil")=SWB,
                     Named("PlantLAI") = PlantLAI,
                     Named("PlantTranspiration") = PlantTranspiration,
                     Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                     Named("PlantPsi") = PlantPsi, 
                     Named("PlantStress") = PlantStress,
                     Named("subdaily") =  subdailyRes);
  } else {
    CharacterVector ln = CharacterVector(28);
    l = List(28);
    l[0] = latitude;
    ln[0] = "latitude";
    l[1] = topo;
    ln[1] = "topography";
    l[2] = spwbInput;
    ln[2] = "spwbInput";
    l[3] = soilInput;
    ln[3] = "soilInput";
    l[4] = DWB;
    ln[4] = "WaterBalance";
    l[5] = SWB;
    ln[5] = "Soil";
    l[6] = DEB;
    ln[6] = "EnergyBalance";
    l[7] = DT;
    ln[7] = "Temperature";
    l[8] = PlantLAI;
    ln[8] = "PlantLAI";
    l[9] = PlantAbsSWR;
    ln[9] = "PlantAbsorbedSWR";
    l[10] = PlantAbsLWR;
    ln[10] = "PlantAbsorbedLWR";
    l[11] = PlantTranspiration;
    ln[11] = "PlantTranspiration";
    l[12] = PlantPhotosynthesis;
    ln[12] = "PlantPhotosynthesis";
    l[13] = dEdP;
    ln[13] = "dEdP";
    l[14] = LeafPsiMin;
    ln[14] = "LeafPsiMin";
    l[15] = LeafPsiMax;
    ln[15] = "LeafPsiMax";
    l[16] = LeafPsiMin_SL;
    ln[16] = "LeafPsiMin_SL";
    l[17] = LeafPsiMax_SL;
    ln[17] = "LeafPsiMax_SL";
    l[18] = LeafPsiMin_SH;
    ln[18] = "LeafPsiMin_SH";
    l[19] = LeafPsiMax_SH;
    ln[19] = "LeafPsiMax_SH";
    l[20] = LeafRWC;
    ln[20] = "LeafRWC";
    l[21] = StemPsi;
    ln[21] = "StemPsi";
    l[22] = StemPLC;
    ln[22] = "StemPLC";
    l[23] = StemRWC;
    ln[23] = "StemRWC";
    l[24] = RootPsi;
    ln[24] = "RootPsi";
    l[25] = RhizoPsi;
    ln[25] = "RhizoPsi";
    l[26] = PlantStress;
    ln[26] = "PlantStress";
    l[27] = subdailyRes;
    ln[27] = "subdaily";
    l.attr("names") = ln;
  }
  l.attr("class") = CharacterVector::create("spwb","list");
  return(l);
}


// [[Rcpp::export("pwb")]]
List pwb(List x, List soil, DataFrame meteo, NumericMatrix W,
            double latitude = NA_REAL, double elevation = NA_REAL, double slope = NA_REAL, double aspect = NA_REAL, 
            NumericVector canopyEvaporation = NumericVector(0), 
            NumericVector snowMelt = NumericVector(0), 
            NumericVector soilEvaporation = NumericVector(0)) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  String soilFunctions = control["soilFunctions"];
  bool verbose = control["verbose"];
  bool subdailyResults = control["subdailyResults"];
  bool leafPhenology = control["leafPhenology"];
  
  
  
  //Store input
  List spwbInput = clone(x);
  List soilInput = clone(soil);
  
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector Radiation;
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  int numDays = Precipitation.size();
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector WindSpeed(numDays, NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector PET = NumericVector(numDays,0.0);
  if(transpirationMode=="Granier") {
    if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
    PET = meteo["PET"];
  } else if(transpirationMode=="Sperry") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
    MinTemperature = meteo["MinTemperature"];
    if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
    MaxTemperature = meteo["MaxTemperature"];
    if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
    MinRelativeHumidity = meteo["MinRelativeHumidity"];
    if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
    MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
    Radiation = meteo["Radiation"];
  }
  
  if(canopyEvaporation.length()==0) {
    canopyEvaporation = NumericVector(numDays,0.0);
  }
  if(snowMelt.length()==0) {
    snowMelt = NumericVector(numDays,0.0);
  }
  if(soilEvaporation.length()==0) {
    soilEvaporation = NumericVector(numDays,0.0);
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = above.nrow();
  

  
  //Soil input
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = Water_FC.size();
  
  //Detailed subday results
  List subdailyRes(numDays);
  
  //Transpiration output variables
  NumericVector GDD(numDays); 
  NumericVector Transpiration(numDays);
  NumericVector PlantExtraction(numDays);
  NumericVector HydraulicRedistribution(numDays, 0.0);
  NumericMatrix HydrIndays(numDays, nlayers);
  
  //EnergyBalance output variables
  NumericVector Tatm_mean(numDays, NA_REAL);
  NumericVector Tatm_min(numDays, NA_REAL);
  NumericVector Tatm_max(numDays, NA_REAL);
  NumericVector Tcan_mean(numDays, NA_REAL);
  NumericVector Tcan_min(numDays, NA_REAL);
  NumericVector Tcan_max(numDays, NA_REAL);
  NumericVector Tsoil_mean(numDays, NA_REAL);
  NumericVector Tsoil_min(numDays, NA_REAL);
  NumericVector Tsoil_max(numDays, NA_REAL);
  NumericVector SWRcanin(numDays, NA_REAL);
  NumericVector LWRcanin(numDays, NA_REAL);
  NumericVector LWRcanout(numDays, NA_REAL);
  NumericVector LEcan_heat(numDays, NA_REAL);
  NumericVector LEsoil_heat(numDays, NA_REAL);
  NumericVector Hcan_heat(numDays, NA_REAL);
  NumericVector Ebalcan(numDays, NA_REAL);
  NumericVector SWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilin(numDays, NA_REAL);
  NumericVector LWRsoilout(numDays, NA_REAL);
  NumericVector Ebalsoil(numDays, NA_REAL);
  NumericVector LWRsoilcan(numDays, NA_REAL);
  NumericVector Hcansoil(numDays, NA_REAL);
  NumericVector RAcan(numDays, NA_REAL);
  NumericVector RAsoil(numDays, NA_REAL);
  

  //Soil output variables
  NumericMatrix Wdays(numDays, nlayers); //Soil moisture content in relation to field capacity
  NumericMatrix psidays(numDays, nlayers);
  NumericMatrix Eplantdays(numDays, nlayers);
  
  //Plant output variables
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix dEdP(numDays, numCohorts);
  NumericMatrix LeafPsiMin(numDays, numCohorts);
  NumericMatrix LeafPsiMax(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SL(numDays, numCohorts);
  NumericMatrix LeafPsiMin_SH(numDays, numCohorts);
  NumericMatrix LeafPsiMax_SH(numDays, numCohorts);
  NumericMatrix StemPsi(numDays, numCohorts);
  NumericMatrix RootPsi(numDays, numCohorts);
  NumericMatrix StemPLC(numDays, numCohorts);
  List RhizoPsi(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    NumericMatrix nm = NumericMatrix(numDays, nlayers);
    nm.attr("dimnames") = List::create(meteo.attr("row.names"), seq(1,nlayers)) ;
    RhizoPsi[c] = nm;
  }
  RhizoPsi.attr("names") = above.attr("row.names");
  
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix StemRWC(numDays, numCohorts);
  NumericMatrix LeafRWC(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericMatrix PlantPhotosynthesis(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  NumericMatrix PlantAbsSWR(numDays, numCohorts);
  NumericMatrix PlantAbsLWR(numDays, numCohorts);
  NumericMatrix PlantLAI(numDays, numCohorts);
  
  

  if(verbose) Rcout << "Performing daily simulations ";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  for(int i=0;i<numDays;i++) {
    if(verbose & (i%10 == 0)) Rcout<<".";//<<i;
    double wind = WindSpeed[i];
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
    if(wind<0.1) wind = 0.1; //Minimum windspeed abovecanopy
    
    //0. Soil moisture
    soil["W"] = W(i,_);
    Wdays(i,_) = W(i,_);
    psidays(i,_) = psi(soil, soilFunctions); //Get soil water potential
      
      //If DOY == 1 reset PLC (Growth assumed)
      if(DOY[i]==1) {
        if(transpirationMode=="Granier") {
          NumericVector PLC = Rcpp::as<Rcpp::NumericVector>(x["PLC"]);
          for(int j=0;j<PLC.length();j++) PLC[j] = 0.0;
        } else {
          NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
          for(int j=0;j<StemPLC.nrow();j++) for(int k=0;k<StemPLC.ncol();k++)  StemPLC(j,k) = 0.0;
        }
      }
      
      
      //1. Phenology and leaf fall
      if(leafPhenology) updateLeaves(x, DOY[i], MeanTemperature[i], wind);
      
      //Store GDD
      GDD[i] = canopyParams["gdd"];
      
      
    
    int ntimesteps = control["ndailysteps"];
    double tstep = 86400.0/((double) ntimesteps);
    
    //2. transpiration and photosynthesis
    if(transpirationMode=="Granier") {
      s = transpirationGranier(x, soil, MeanTemperature[i], PET[i], true);
    } else if(transpirationMode=="Sperry") {
      std::string c = as<std::string>(dateStrings[i]);
      int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
      double delta = meteoland::radiation_solarDeclination(J);
      double solarConstant = meteoland::radiation_solarConstant(J);
      double tmin = MinTemperature[i];
      double tmax = MaxTemperature[i];
      double rhmin = MinRelativeHumidity[i];
      double rhmax = MaxRelativeHumidity[i];
      double rad = Radiation[i];
      double prec = Precipitation[i];
      
      s = transpirationSperry(x, soil, tmin, tmax, rhmin, rhmax, rad, wind, 
                              latitude, elevation, slope, aspect,
                              solarConstant, delta, prec,
                              canopyEvaporation[i], snowMelt[i], soilEvaporation[i],
                              verbose, NA_INTEGER, true);
      
    }
    List Plants = s["Plants"];
    PlantLAI(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LAI"]);
    NumericVector EplantCoh = Plants["Transpiration"];
    NumericMatrix SoilWaterExtract = s["Extraction"];
    for(int l=0;l<nlayers;l++) {
      Eplantdays(i,l) = sum(SoilWaterExtract(_,l));
    }

    PlantExtraction[i] = sum(SoilWaterExtract);
    Transpiration[i] = sum(EplantCoh);
    NumericVector HydrInVec(nlayers, 0.0);
    PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["Photosynthesis"]);
    PlantTranspiration(i,_) = EplantCoh;
    PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
    
    if(transpirationMode=="Granier") {
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["psi"]);
    }
    else if(transpirationMode=="Sperry")  {
      NumericMatrix soilLayerExtractInst = s["ExtractionInst"];
      for(int l=0;l<nlayers;l++) {
        for(int n=0;n<ntimesteps;n++) {
          HydrInVec[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
        }
      }
      HydraulicRedistribution[i] = sum(HydrInVec);
      
      List PlantsInst = Rcpp::as<Rcpp::List>(s["PlantsInst"]);
      List ShadeLeaves = Rcpp::as<Rcpp::List>(PlantsInst["ShadeLeaves"]);
      List SunlitLeaves = Rcpp::as<Rcpp::List>(PlantsInst["SunlitLeaves"]);
      
      NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_SWR"]);
      NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_SWR"]);
      NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitLeaves["Abs_LWR"]);
      NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeLeaves["Abs_LWR"]);
      for(int j=0;j<numCohorts;j++) {
        for(int n=0;n<ntimesteps;n++){
          PlantAbsSWR(i,j) += 0.000001*(SWR_SL(j,n)+SWR_SH(j,n))*tstep;
          PlantAbsLWR(i,j) += 0.000001*(LWR_SL(j,n)+LWR_SH(j,n))*tstep;
        }
      }
      List EB = Rcpp::as<Rcpp::List>(s["EnergyBalance"]);
      DataFrame Tinst = Rcpp::as<Rcpp::DataFrame>(EB["Temperature"]); 
      DataFrame CEBinst = Rcpp::as<Rcpp::DataFrame>(EB["CanopyEnergyBalance"]); 
      DataFrame SEBinst = Rcpp::as<Rcpp::DataFrame>(EB["SoilEnergyBalance"]); 
      NumericVector Tcan = Rcpp::as<Rcpp::NumericVector>(Tinst["Tcan"]);
      NumericVector Tsoil = Rcpp::as<Rcpp::NumericVector>(Tinst["Tsoil.1"]);
      
      SWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["SWRcanin"]))*tstep;
      LWRcanin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanin"]))*tstep;
      LWRcanout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRcanout"]))*tstep;
      LEcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LEcan"]))*tstep;
      Hcan_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Hcan"]))*tstep;
      Ebalcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["Ebalcan"]))*tstep;
      LWRsoilcan[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["LWRsoilcan"]))*tstep;
      RAcan[i] = sum(Rcpp::as<Rcpp::NumericVector>(CEBinst["RAcan"]))/((double) ntimesteps);
      SWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["SWRsoilin"]))*tstep;
      LWRsoilin[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilin"]))*tstep;
      LWRsoilout[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LWRsoilout"]))*tstep;
      LEsoil_heat[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["LEsoil"]))*tstep;
      Hcansoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Hcansoil"]))*tstep;
      Ebalsoil[i] = 0.000001*sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["Ebalsoil"]))*tstep;
      RAsoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(SEBinst["RAsoil"]))/((double) ntimesteps);
      
      Tatm_min[i] = MinTemperature[i];
      Tatm_max[i] = MaxTemperature[i];
      Tatm_mean[i] = MeanTemperature[i];
      Tcan_min[i] = min(Tcan);
      Tcan_max[i] = max(Tcan);
      Tcan_mean[i] = sum(Tcan)/((double) ntimesteps);
      Tsoil_min[i] = min(Tsoil);
      Tsoil_max[i] = max(Tsoil);
      Tsoil_mean[i] = sum(Tsoil)/((double) ntimesteps);
      HydrIndays(i,_) = HydrInVec;
      LeafPsiMin(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
      LeafPsiMax(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
      LeafPsiMin_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SL"]);
      LeafPsiMax_SL(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SL"]);
      LeafPsiMin_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin_SH"]);
      LeafPsiMax_SH(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax_SH"]);
      RootPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["RootPsi"]); 
      StemPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPsi"]); 
      StemPLC(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["StemPLC"]); 
      dEdP(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["dEdP"]); 
      NumericMatrix RhizoPsiStep = Rcpp::as<Rcpp::NumericMatrix>(s["RhizoPsi"]);
      for(int c=0;c<numCohorts;c++) {
        NumericMatrix nm = Rcpp::as<Rcpp::NumericMatrix>(RhizoPsi[c]);
        nm(i,_) =  RhizoPsiStep(c,_);
      }
      StemRWC(i,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
      LeafRWC(i,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
    } 
    EplantCohTot = EplantCohTot + EplantCoh;
    Eplanttot[i] = sum(EplantCoh);
    
    if(subdailyResults) {
      subdailyRes[i] = clone(s);
    }
  }
  if(verbose) Rcout << "done\n";
  
  if(verbose) {
    double Transpirationsum = sum(Transpiration);
    
    Rcout<<"Transpiration (mm) "  <<round(Transpirationsum);
    if(transpirationMode =="Sperry") {
      Rcout<<" Plant extraction from soil (mm) " << round(sum(PlantExtraction));
      Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    } else {
      Rcout <<"\n";
    }
  }

  
  DataFrame SWB;
  if(transpirationMode=="Granier") {
    SWB = DataFrame::create(_["W"]=Wdays, _["PlantExt"]=Eplantdays, _["psi"]=psidays); 
  } else {
    SWB = DataFrame::create(_["W"]=Wdays, _["PlantExt"]=Eplantdays,
                            _["HydraulicInput"] = HydrIndays,
                            _["psi"]=psidays); 
  }
  SWB.attr("row.names") = meteo.attr("row.names") ;
  
  DataFrame DWB;
  if(transpirationMode=="Granier") {
    DWB = DataFrame::create(_["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration);
  } else {
    DWB = DataFrame::create(_["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                            _["HydraulicRedistribution"] = HydraulicRedistribution);
  }
  DWB.attr("row.names") = meteo.attr("row.names") ;
  
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPLC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  dEdP.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SL.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax_SH.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  RootPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsSWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsLWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantLAI.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  
  DataFrame DEB = DataFrame::create(_["SWRcanin"] = SWRcanin, _["LWRcanin"] = LWRcanin, _["LWRcanout"] = LWRcanout,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebalcan, 
                                      _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = SWRsoilin, _["LWRsoilin"] = LWRsoilin, _["LWRsoilout"] = LWRsoilout,
                                        _["LEsoil"] = LEsoil_heat, _["Ebalsoil"] = Ebalsoil, _["RAcan"] = RAcan, _["RAsoil"] = RAsoil);  
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                     _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DEB.attr("row.names") = meteo.attr("row.names") ;
  DT.attr("row.names") = meteo.attr("row.names") ;
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  
  NumericVector topo = NumericVector::create(elevation, slope, aspect);
  topo.attr("names") = CharacterVector::create("elevation", "slope", "aspect");

  List l;
  if(transpirationMode=="Granier") {
    l = List::create(Named("latitude") = latitude,
                     Named("topography") = topo,
                     Named("spwbInput") = spwbInput,
                     Named("soilInput") = soilInput,
                     Named("WaterBalance")=DWB, 
                     Named("Soil") = SWB,
                     Named("PlantLAI") = PlantLAI,
                     Named("PlantTranspiration") = PlantTranspiration,
                     Named("PlantPhotosynthesis") = PlantPhotosynthesis,
                     Named("PlantPsi") = PlantPsi, 
                     Named("PlantStress") = PlantStress,
                     Named("subdaily") =  subdailyRes);
  } else {
    CharacterVector ln = CharacterVector(28);
    l = List(28);
    l[0] = latitude;
    ln[0] = "latitude";
    l[1] = topo;
    ln[1] = "topography";
    l[2] = spwbInput;
    ln[2] = "spwbInput";
    l[3] = soilInput;
    ln[3] = "soilInput";
    l[4] = DWB;
    ln[4] = "WaterBalance";
    l[5] = SWB;
    ln[5] = "Soil";
    l[6] = DEB;
    ln[6] = "EnergyBalance";
    l[7] = DT;
    ln[7] = "Temperature";
    l[8] = PlantLAI;
    ln[8] = "PlantLAI";
    l[9] = PlantAbsSWR;
    ln[9] = "PlantAbsorbedSWR";
    l[10] = PlantAbsLWR;
    ln[10] = "PlantAbsorbedLWR";
    l[11] = PlantTranspiration;
    ln[11] = "PlantTranspiration";
    l[12] = PlantPhotosynthesis;
    ln[12] = "PlantPhotosynthesis";
    l[13] = dEdP;
    ln[13] = "dEdP";
    l[14] = LeafPsiMin;
    ln[14] = "LeafPsiMin";
    l[15] = LeafPsiMax;
    ln[15] = "LeafPsiMax";
    l[16] = LeafPsiMin_SL;
    ln[16] = "LeafPsiMin_SL";
    l[17] = LeafPsiMax_SL;
    ln[17] = "LeafPsiMax_SL";
    l[18] = LeafPsiMin_SH;
    ln[18] = "LeafPsiMin_SH";
    l[19] = LeafPsiMax_SH;
    ln[19] = "LeafPsiMax_SH";
    l[20] = LeafRWC;
    ln[20] = "LeafRWC";
    l[21] = StemPsi;
    ln[21] = "StemPsi";
    l[22] = StemPLC;
    ln[22] = "StemPLC";
    l[23] = StemRWC;
    ln[23] = "StemRWC";
    l[24] = RootPsi;
    ln[24] = "RootPsi";
    l[25] = RhizoPsi;
    ln[25] = "RhizoPsi";
    l[26] = PlantStress;
    ln[26] = "PlantStress";
    l[27] = subdailyRes;
    ln[27] = "subdaily";
    l.attr("names") = ln;
  }
  l.attr("class") = CharacterVector::create("pwb","list");
  return(l);                    
}