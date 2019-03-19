#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "hydrology.h"
#include "biophysicsutils.h"
#include "forestutils.h"
#include "photosynthesis.h"
#include "phenology.h"
#include "transpiration.h"
#include "tissuemoisture.h"
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*pow(10,-8.0);
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1



// Soil water balance with simple hydraulic model
// [[Rcpp::export(".spwbDay1")]]
List spwbDay1(List x, List soil, double tday, double pet, double prec, double er, double runon=0.0, 
             double rad = NA_REAL, double elevation = NA_REAL, bool verbose = false) {

  //Control parameters
  List control = x["control"];
  bool snowpack = control["snowpack"];
  bool drainage = control["drainage"];
  bool cavitationRefill = control["cavitationRefill"];
  String soilFunctions = control["soilFunctions"];

  
  
  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector psiVec = psi(soil, soilFunctions);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Water_SAT = waterSAT(soil, soilFunctions);
  double swe = soil["SWE"]; //snow pack
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

  //Root distribution input
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(below["V"]);

  //Parameters  
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsBase["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsBase["g"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE"]);
  NumericVector pRootDisc = Rcpp::as<Rcpp::NumericVector>(paramsTransp["pRootDisc"]);
  
  //Communication vectors
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericVector pEmb = Rcpp::as<Rcpp::NumericVector>(x["PLC"]);
  


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
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  LAIdead, H, CR, kPAR);
  NumericVector CohPAR = parcohortC(H, LAIphe, LAIdead, kPAR, CR)/100.0;
  double LgroundPAR = exp((-1.0)*s);
  double LgroundSWR = 1.0 - std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  
  //Snow pack dynamics
  double snow = 0.0, rain=0.0;
  double melt = 0.0;
  if(snowpack) {
    //Turn rain into snow and add it into the snow pack
    if(tday < 0.0) { 
      snow = prec; 
      swe = swe + snow;
    } else {
      rain = prec;
    }
    //Apply snow melting
    if(swe > 0.0) {
      if(NumericVector::is_na(rad)) stop("Missing radiation data for snow melt!");
      if(NumericVector::is_na(elevation)) stop("Missing elevation data for snow melt!");
      double rho = meteoland::utils_airDensity(tday, meteoland::utils_atmosphericPressure(elevation));
      double ten = (86400*tday*rho*1013.86*pow(10,-6.0)/100.0); //ten can be negative if temperature is below zero
      double ren = (rad*LgroundSWR)*(0.1); //90% albedo of snow
      melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
      // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
      swe = std::max(0.0, swe-melt);
    }
    soil["SWE"] = swe;
  } else {
    rain = prec;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(rain>0.0) NetRain = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
  if((NetRain+runon+melt)>0.0) {
    //Interception
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetRain+runon+melt, Water_FC[0]);
    Runoff = (NetRain+runon+melt) - Infiltration;
    //Decide infiltration repartition among layers
    NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
    //Input of the first soil layer is infiltration
    double excess = 0.0;
    double Wn;
    //Update topsoil layer
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0.0) & (Ivec[l]>0.0)) {
        //PROBLEM: THE effect of MACROPOROSITY SHOULD not be affected by layer subdivision
        Wn = W[l]*Water_FC[l] + Ivec[l]; //Update water volume
        if(l<(nlayers-1)) {
          Ivec[l+1] = Ivec[l+1] + std::max(Wn - Water_FC[l],0.0); //update Ivec adding the excess to the infiltrating water (saturated flow)
        } else {
          excess = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer
        }
        W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta (this modifies 'soil')
      } 
    }
    if(drainage) {//Set deep drainage
      DeepDrainage = excess; 
    } else { //Fill to saturation and upwards if needed
      for(int l=(nlayers-1);l>=0;l--) {
        if((dVec[l]>0.0) & (excess>0.0)) {
          Wn = W[l]*Water_FC[l] + excess; //Update water volume
          excess = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
          W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
        }
      }
      if(excess>0.0) { //If soil is completely saturated increase Runoff
        Runoff = Runoff + excess;
      }
    }
  }
  psiVec = psi(soil,soilFunctions); //Update soil water potential


  //Proportion of transpiration that absorbed by each plant cohort (old version)
  // NumericVector PP = CohLight*LAIphe;
  // NumericVector f = PP/std::accumulate(PP.begin(),PP.end(),0.0); 
  // if(LAIcell==0.0) f = rep(0.0,numCohorts); //Avoids NaN values

  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  double Tmax = pet*(-0.006*pow(LAIcell,2.0)+0.134*LAIcell); //From Granier (1999)
  double PETsoil = pet*LgroundSWR;

  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts,0.0);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);

  //Actual plant transpiration
  NumericMatrix EplantCoh(numCohorts, nlayers);
  NumericMatrix PsiRoot(numCohorts, nlayers);
  NumericVector PlantPsi(numCohorts, NA_REAL);
  NumericVector Eplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0);
  NumericVector EplantVec(nlayers, 0.0);
  NumericVector Kl, epc, Vl;
  double WeibullShape=3.0;
  for(int l=0;l<nlayers;l++) {
    Kl = Psi2K(psiVec[l], Psi_Extract, WeibullShape);
 
    //Limit Kl due to previous cavitation
    if(!cavitationRefill) for(int c=0;c<numCohorts;c++) Kl[c] = std::min(Kl[c], 1.0-pEmb[c]);
    //Limit Kl to minimum value for root disconnection
    Vl = V(_,l);
    epc = pmax(TmaxCoh*Kl*Vl,0.0);
    for(int c=0;c<numCohorts;c++) {
      PsiRoot(c,l) = psiVec[l]; //Set initial guess of root potential to soil values
      //If relative conductance is smaller than the value for root disconnection
      //Set root potential to minimum value before disconnection and transpiration from that layer to zero
      if(Kl[c]<pRootDisc[c]) { 
        PsiRoot(c,l) = K2Psi(pRootDisc[c],Psi_Extract[c],WeibullShape);
        Kl[c] = pRootDisc[c]; //So that layer stress does not go below pRootDisc
        epc[c] = 0.0; //Set transpiration from layer to zero
      }
    }
    
    EplantCoh(_,l) = epc;
    Eplant = Eplant + epc;
    EplantVec[l] = std::accumulate(epc.begin(),epc.end(),0.0);
    DDS = DDS + Phe*(Vl*(1.0 - Kl)); //Add stress from the current layer
  }
  for(int c=0;c<numCohorts;c++) {
    PlantPsi[c] = averagePsi(PsiRoot(c,_), V(c,_), WeibullShape, Psi_Extract[c]);
    if(!cavitationRefill) {
      pEmb[c] = std::max(DDS[c], pEmb[c]); //Track current embolism if no refill
      DDS[c] = pEmb[c];
    }
  }

  //Evaporation from bare soil
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  double Esoil = soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
  NumericVector EsoilVec(nlayers,0.0);//Exponential decay to divide bare soil evaporation among layers
  for(int l=0;l<nlayers;l++) {
    double cumAnt = 0.0;
    double cumPost = 0.0;
    for(int l2=0;l2<l;l2++) cumAnt +=dVec[l2];
    cumPost = cumAnt+dVec[l];
    if(l<(nlayers-1)) EsoilVec[l] = Esoil*(exp(-Ksoil*cumAnt)-exp(-Ksoil*cumPost));
    else EsoilVec[l] = Esoil*exp(-Ksoil*cumAnt);
    
  }
  // Rcout<<Esoil<<" "<<  std::accumulate(EsoilVec.begin(),EsoilVec.end(),0.0)<<"\n";

  //Apply decrease in soil layers (psi will be updated next call)
  for(int l=0;l<nlayers;l++) W[l] = W[l] - ((EplantVec[l]+EsoilVec[l])/Water_FC[l]);

  double alpha = std::max(std::min(tday/20.0,1.0),0.0);
  //For comunication with growth and landscape water balance
  for(int c=0;c<numCohorts;c++) {
    transpiration[c] = Eplant[c];
    photosynthesis[c] = alpha*WUE[c]*transpiration[c];
  }

  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  
  NumericVector DB = NumericVector::create(_["PET"] = pet, _["Rain"] = rain, _["Snow"] = snow, _["NetRain"] = NetRain, _["Runon"] = runon, 
                                           _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(EplantVec), _["Transpiration"] = sum(EplantVec),
                                           _["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, 
                                           _["Cm"] = Cm, _["Lground"] = LgroundPAR);
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["PlantExtraction"] = EplantVec, 
                                   _["psi"] = psiVec);
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                             _["Transpiration"] = Eplant, _["psi"] = PlantPsi, _["DDS"] = DDS);
  Plants.attr("row.names") = above.attr("row.names");
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = DB, 
                        _["Soil"] = SB,
                        _["Plants"] = Plants);
  l.attr("class") = CharacterVector::create("spwb.day","list");
  return(l);
}


// Soil water balance with Sperry hydraulic and stomatal conductance models
// [[Rcpp::export(".spwbDay2")]]
List spwbDay2(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation, double solarConstant, double delta, 
             double prec, double er, double runon=0.0, bool verbose = false) {
  
  //Control parameters
  List control = x["control"];
  bool drainage = control["drainage"];
  String soilFunctions = control["soilFunctions"];
  List numericParams = control["numericParams"];
  int ntrial = numericParams["ntrial"];
  int maxNsteps  = numericParams["maxNsteps"];
  double psiTol = numericParams["psiTol"];
  double ETol = numericParams["ETol"];

  bool capacitance = control["capacitance"];
  bool cavitationRefill = control["cavitationRefill"];
  double klat = control["klat"];
  int ntimesteps = control["ndailysteps"];
  int hydraulicCostFunction = control["hydraulicCostFunction"];
  int nStemSegments = control["nStemSegments"];
  double verticalLayerSize = control["verticalLayerSize"];
  double thermalCapacityLAI = control["thermalCapacityLAI"];
  double defaultWindSpeed = control["defaultWindSpeed"];
  
    
  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector psiVec = psi(soil, soilFunctions);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector Water_SAT = waterSAT(soil, soilFunctions);
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector Tsoil = soil["Temp"]; 
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

  //Canopy params
  List canopyParams = Rcpp::as<Rcpp::List>(x["canopy"]);
  
  //Root distribution input
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(below["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VGrhizo_kmax"]);
  
  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(paramsBase["Sgdd"]);
  NumericVector albedo = Rcpp::as<Rcpp::NumericVector>(paramsBase["albedo"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(paramsBase["g"]);

  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector Gwmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gwmin"]);
  NumericVector Gwmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gwmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector Vmax298 = paramsTransp["Vmax298"];
  NumericVector Jmax298 = paramsTransp["Jmax298"];
  NumericVector pRootDisc = Rcpp::as<Rcpp::NumericVector>(paramsTransp["pRootDisc"]);

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]);
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]);

  //Comunication with outside
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericMatrix psiStemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["psiStem"]);
  NumericMatrix PLCstemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
  NumericMatrix RWCsstemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["RWCsympstem"]);
  NumericVector RWCsleafVEC = Rcpp::as<Rcpp::NumericVector>(x["RWCsympleaf"]);
  NumericVector psiLeafVEC = Rcpp::as<Rcpp::NumericVector>(x["psiLeaf"]);
  NumericVector psiRootVEC = Rcpp::as<Rcpp::NumericVector>(x["psiRoot"]);
  NumericMatrix psiRhizoMAT = Rcpp::as<Rcpp::NumericMatrix>(x["psiRhizo"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(x["Einst"]);
  
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  
  double latrad = latitude * (PI/180.0);
  // if(NumericVector::is_na(aspect)) aspect = 0.0;
  // double asprad = aspect * (PI/180.0);
  // if(NumericVector::is_na(slope)) slope = 0.0;
  // double slorad = slope * (PI/180.0);

  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);
  
  //Atmospheric pressure, CO2 concentration
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  double Catm = control["Catm"];
  

  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  
  //Daily cloud cover
  double cloudcover = 0.0;
  if(prec >0.0) cloudcover = 1.0;
  bool clearday = (prec==0);
  
  //1. Leaf Phenology: Adjusted leaf area index
  NumericVector Phe(numCohorts);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0, Cm = 0.0, canopyHeight = 0.0, LAIcellmax = 0.0;
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcellmax += LAIlive[c];
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    Cm += (LAIphe[c]+LAIdead[c])*gRainIntercept[c]; //LAI dead also counts on interception
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  
  NumericVector z(nz+1,0.0);
  NumericVector zmid(nz);
  for(int i=1;i<=nz;i++) {
    z[i] = z[i-1] + verticalLayerSize;
    zmid[i-1] = (verticalLayerSize/2.0) + verticalLayerSize*((double) (i-1));
  }
  NumericMatrix LAIme = LAIdistribution(z, LAIphe, H, CR); //Expanded leaves
  NumericMatrix LAImd = LAIdistribution(z, LAIdead, H, CR); //Dead (standing) leaves
  NumericMatrix LAImx = LAIdistribution(z, LAIlive, H, CR); //Maximum leaf expansion
  
  
  //2. Hydrologic input
  double rain = prec;
  double NetRain = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  double LgroundPAR = exp((-1.0)*s);
  // double propCoverMax = 1.0-exp(-1.0*LAIcellmax);
  if(rain>0.0) {
    //Interception
    NetRain = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetRain+runon, Water_FC[0]);
    Runoff = (NetRain+runon) - Infiltration;
    
    //Decide infiltration repartition among layers
    NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
    //Input of the first soil layer is infiltration
    double excess = 0.0;
    double Wn;
    //Update topsoil layer
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0.0) & (Ivec[l]>0.0)) {
        //PROBLEM: THE effect of MACROPOROSITY SHOULD not be affected by layer subdivision
        Wn = W[l]*Water_FC[l] + Ivec[l]; //Update water volume
        if(l<(nlayers-1)) {
          Ivec[l+1] = Ivec[l+1] + std::max(Wn - Water_FC[l],0.0); //update Ivec adding the excess to the infiltrating water (saturated flow)
        } else {
          excess = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer
        }
        W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta (this modifies 'soil')
      } 
    }
    if(drainage) {//Set deep drainage
      DeepDrainage = excess; 
    } else { //Fill to saturation and upwards if needed
      for(int l=(nlayers-1);l>=0;l--) {
        if((dVec[l]>0.0) & (excess>0.0)) {
          Wn = W[l]*Water_FC[l] + excess; //Update water volume
          excess = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
          W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
        }
      }
      if(excess>0.0) { //If soil is completely saturated increase Runoff
        Runoff = Runoff + excess;
      }
    }
  }
  
  psiVec = psi(soil, soilFunctions); //Update soil water potential
  
  
  //3. Wind extinction profile
  if(NumericVector::is_na(wind)) wind = defaultWindSpeed; //set to default if missing
  NumericVector zWind;
  zWind = windExtinctionCohort(H,CR, wind,LAIcell, canopyHeight);
  double RAcan = aerodynamicResistance(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
  double wind2m = windSpeedMassmanExtinction(200.0, wind, LAIcell, canopyHeight);
  double RAsoil = aerodynamicResistance(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
  
  //4a. Instantaneous direct and diffuse shorwave radiation
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, delta,
                                                        rad, clearday, ntimesteps);
  NumericVector solarElevation = ddd["SolarElevation"]; //in radians
  NumericVector solarHour = ddd["SolarHour"]; //in radians
  NumericVector SWR_direct = ddd["SWR_direct"]; //in kW·m-2
  NumericVector SWR_diffuse = ddd["SWR_diffuse"]; //in kW·m-2
  NumericVector PAR_direct = ddd["PAR_direct"]; //in kW·m-2
  NumericVector PAR_diffuse = ddd["PAR_diffuse"]; //in kW·m-2
  
  //4b. Instantaneous air temperature (above canopy) and longwave radiation
  NumericVector Tatm(ntimesteps), lwdr(ntimesteps), Tcan(ntimesteps, NA_REAL), Tsunrise(ntimesteps);
  NumericVector LEcan_heat(ntimesteps), Hcan_heat(ntimesteps), LWRsoilcan(ntimesteps), LWRcanout(ntimesteps), Ebal(ntimesteps);
  NumericVector LWRsoilout(ntimesteps), Ebalsoil(ntimesteps), Hcansoil(ntimesteps);
  NumericMatrix Tsoil_mat(ntimesteps, nlayers);
  //Daylength in seconds (assuming flat area because we want to model air temperature variation)
  double tauday = meteoland::radiation_daylengthseconds(latrad,0.0,0.0, delta); 
  for(int n=0;n<ntimesteps;n++) {
    //From solar hour (radians) to seconds from sunrise
    Tsunrise[n] = (solarHour[n]*43200.0/PI)+ (tauday/2.0) +(tstep/2.0); 
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern(Tsunrise[n], tmin, tmax, tauday);
    //Longwave sky diffuse radiation (W/m2)
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpatm, cloudcover);
  }
  Tcan[0] = canopyParams["Temp"]; //Take canopy temperature from previous day
  Tsoil_mat(0,_) = Tsoil;

  

  //4c. Light extinction and absortion by time steps
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        kPAR, albedo,
                                                                        ddd,  lwdr,
                                                                        ntimesteps,  "sunshade", 0.1);
  
  List abs_PAR_SL_list = lightExtinctionAbsortion["PAR_SL"];
  List abs_SWR_SL_list = lightExtinctionAbsortion["SWR_SL"];
  List abs_PAR_SH_list = lightExtinctionAbsortion["PAR_SH"];
  List abs_SWR_SH_list = lightExtinctionAbsortion["SWR_SH"];
  List abs_LWR_SL_list = lightExtinctionAbsortion["LWR_SL"];
  List abs_LWR_SH_list = lightExtinctionAbsortion["LWR_SH"];
  NumericVector fsunlit = lightExtinctionAbsortion["fsunlit"];
  NumericVector abs_SWR_can = lightExtinctionAbsortion["SWR_can"];
  NumericVector abs_SWR_soil = lightExtinctionAbsortion["SWR_soil"];
  NumericVector abs_LWR_can = lightExtinctionAbsortion["LWR_can"];
  NumericVector abs_LWR_soil = lightExtinctionAbsortion["LWR_soil"];
  NumericVector emm_LWR_soil(ntimesteps,0.0);
  // double kb = lightExtinctionAbsortion["kb"];  //Proportion of sunlit extinction coefficient
  // double gbf = lightExtinctionAbsortion["gbf"]; //Ground fractions
  // double gdf = lightExtinctionAbsortion["gdf"];


 

  //Hydraulics: determine layers where the plant is connected
  IntegerVector nlayerscon(numCohorts,0);
  NumericMatrix SoilWaterExtract(numCohorts, nlayers);
  NumericMatrix soilLayerExtractInst(nlayers, ntimesteps);
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);
  
  LogicalMatrix layerConnected(numCohorts, nlayers);
  for(int c=0;c<numCohorts;c++) {
    nlayerscon[c] = 0;
    for(int l=0;l<nlayers;l++) {
      SoilWaterExtract(c,l) = 0.0;
      if(V(c,l)>0.0) {
        double pRoot = xylemConductance(psiVec[l], 1.0, VCroot_c[c], VCroot_d[c]); //Relative conductance in the root
        layerConnected(c,l)= (pRoot>=pRootDisc[c]);
        if(layerConnected(c,l)) nlayerscon[c]=nlayerscon[c]+1;
      } else {
        layerConnected(c,l) = false;
      }
    }
    if((nlayerscon[c]==0) & verbose) Rcout<<"D";
  }

  //Hydraulics: build supply functions
  List supply(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    // Copy values from connected layers
    NumericVector Vc = NumericVector(nlayerscon[c]);
    NumericVector VCroot_kmaxc = NumericVector(nlayerscon[c]);
    NumericVector VGrhizo_kmaxc = NumericVector(nlayerscon[c]);
    NumericVector psic = NumericVector(nlayerscon[c]);
    NumericVector VG_nc = NumericVector(nlayerscon[c]);
    NumericVector VG_alphac= NumericVector(nlayerscon[c]);
    int cnt=0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected(c,l)) {
        Vc[cnt] = V(c,l);
        VCroot_kmaxc[cnt] = VCroot_kmax(c,l);
        VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l);
        psic[cnt] = psiVec[l];
        VG_nc[cnt] = VG_n[l];
        VG_alphac[cnt] = VG_alpha[l];
        cnt++;
      }
    }
    // double minFlow = std::max(0.0,1000.0*(Gwmin[c]*(tmin+tmax)/2.0)/Patm);
    // Rcout<<minFlow<<"\n";
    if(nlayerscon[c]>0) {
      NumericVector PLCStemPrev = PLCstemMAT(c,_); //Get row
      // if(!capacitance) {
        supply[c] = supplyFunctionNetwork(psic,
                                          VGrhizo_kmaxc,VG_nc,VG_alphac,
                                          VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                          VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                          VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                          PLCStemPrev,
                                          0.0, maxNsteps, 
                                          ntrial, psiTol, ETol,
                                          0.001); 
      // } else { //Calculate supply function for root system only
      //   supply[c] = supplyFunctionBelowground(psic,
      //                                         VGrhizo_kmaxc,VG_nc,VG_alphac,
      //                                         VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
      //                                         0.0, maxNsteps, psiStep, psiMax, 
      //                                         ntrial, psiTol, ETol); 
      //   
      // }
    } else {
      Rcout<<"D";
    }
  }

  //Transpiration and photosynthesis
  NumericVector psiBk(nlayers);
  for(int l=0;l<nlayers;l++) psiBk[l] = psiVec[l]; //Store initial soil water potential
  NumericMatrix K(numCohorts, nlayers);
  NumericVector Eplant(numCohorts);
  NumericMatrix Rninst(numCohorts,ntimesteps);
  NumericMatrix dEdPinst(numCohorts, ntimesteps);
  NumericMatrix Qinst(numCohorts,ntimesteps);
  NumericMatrix Einst(numCohorts, ntimesteps);
  NumericMatrix Aninst(numCohorts, ntimesteps);
  NumericMatrix PsiLeafinst(numCohorts, ntimesteps);
  NumericMatrix PsiSteminst(numCohorts, ntimesteps);
  NumericMatrix RWCleafinst(numCohorts, ntimesteps);
  NumericMatrix RWCsteminst(numCohorts, ntimesteps);
  NumericMatrix PsiRootinst(numCohorts, ntimesteps);
  NumericMatrix PWBinst(numCohorts, ntimesteps);
  NumericMatrix SWR_SL(numCohorts, ntimesteps);
  NumericMatrix SWR_SH(numCohorts, ntimesteps);
  NumericMatrix LWR_SL(numCohorts, ntimesteps);
  NumericMatrix LWR_SH(numCohorts, ntimesteps);
  NumericMatrix GW_SH(numCohorts, ntimesteps);
  NumericMatrix GW_SL(numCohorts, ntimesteps);
  NumericMatrix VPD_SH(numCohorts, ntimesteps);
  NumericMatrix VPD_SL(numCohorts, ntimesteps);
  NumericMatrix Temp_SH(numCohorts, ntimesteps);
  NumericMatrix Temp_SL(numCohorts, ntimesteps);
  NumericMatrix LAI_SH(numCohorts, ntimesteps);
  NumericMatrix LAI_SL(numCohorts, ntimesteps);
  NumericVector minPsiLeaf(numCohorts,0.0), maxPsiLeaf(numCohorts,-99999.0), minPsiStem(numCohorts, 0.0), minPsiRoot(numCohorts,0.0); //Minimum potentials experienced
  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), 0.0);
  NumericMatrix PLC(numCohorts, ntimesteps);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts);
  NumericVector dEdPm(numCohorts);
  
  
  for(int c=0;c<numCohorts;c++) {
    photosynthesis[c] = 0.0;
    transpiration[c] = 0.0;
  }

  for(int n=0;n<ntimesteps;n++) { //Time loop
    //Long-wave radiation due to canopy temperature
    if(NumericVector::is_na(Tcan[n])) Tcan[n] = Tatm[n]; //If missing take above-canopy air temperature
    if(NumericVector::is_na(Tsoil[0])) {//Initialize Soil temperature (to minimum air temperature) if missing
      for(int l=0;l<nlayers; l++) {
        Tsoil[l] = Tatm[n];
      }
      Tsoil_mat(n,_) = Tsoil; 
    }
    //LWR emmited by the canopy, per ground area
    double LWR_emmcan = 0.95*SIGMA_Wm2*pow(Tcan[n]+273.16,4.0);
    //Soil longwave emmission
    emm_LWR_soil[n] =  0.95*SIGMA_Wm2*pow(Tsoil[0]+273.16,4.0);
    
    //Retrieve radiation absorbed
    NumericVector absPAR_SL = abs_PAR_SL_list[n];
    NumericVector absPAR_SH = abs_PAR_SH_list[n];
    NumericVector absSWR_SL = abs_SWR_SL_list[n];
    NumericVector absSWR_SH = abs_SWR_SH_list[n];
    NumericVector absLWR_SL = abs_LWR_SL_list[n];
    NumericVector absLWR_SH = abs_LWR_SH_list[n];
    
    for(int c=0;c<numCohorts;c++) { //Plant cohort loop

      if(LAIphe[c]>0.0) { //Process transpiration and photosynthesis only if there are some leaves
        SWR_SL(c,n) = absSWR_SL[c];
        SWR_SH(c,n) = absSWR_SH[c];
        LWR_SL(c,n) = absLWR_SL[c];
        LWR_SH(c,n) = absLWR_SH[c];
        
        
        //Constant properties through time steps
        NumericVector Vmax298layer(nz), Jmax298layer(nz);
        NumericVector SLarealayer(nz), SHarealayer(nz);
        NumericVector QSH(nz), QSL(nz), absRadSL(nz), absRadSH(nz);
        double sn =0.0;
        for(int i=(nz-1);i>=0.0;i--) {
          //Effect of nitrogen concentration decay through the canopy
          double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sum(LAIme(_,c)));
          sn+=LAIme(i,c);
          SLarealayer[i] = LAIme(i,c)*fsunlit[i];
          SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
          Vmax298layer[i] = Vmax298[c]*fn;
          Jmax298layer[i] = Jmax298[c]*fn;
        }
        double Vmax298SL= 0.0,Vmax298SH= 0.0,Jmax298SL= 0.0,Jmax298SH= 0.0;
        LAI_SH(c,n) = 0.0; 
        LAI_SL(c,n) = 0.0;
        for(int i=0;i<nz;i++) {
          LAI_SL(c,n) +=SLarealayer[i];
          LAI_SH(c,n) +=SHarealayer[i];
          Vmax298SL +=Vmax298layer[i]*LAIme(i,c)*fsunlit[i];
          Jmax298SL +=Jmax298layer[i]*LAIme(i,c)*fsunlit[i];
          Vmax298SH +=Vmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
          Jmax298SH +=Jmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
        }
        
        NumericVector PLCStemPrev = PLCstemMAT(c,_);
        NumericVector RWCStemPrev = RWCsstemMAT(c,_);
        NumericVector psiStemPrev = psiStemMAT(c,_);
        double psiLeafPrev = psiLeafVEC[c];
        double rwcsleafPrev = RWCsleafVEC[c];
        
        if(nlayerscon[c]>0) {//If the plant is connected to at least one layer build 
          
          List sFunction;
          // Rcout<<c<<" E "<<EinstPrev<<" PR "<< psiRootPrev<<" PL "<<psiLeafPrev<< " PS "<<psiStemPrev[0]<< " "<<rwcsleafPrev<< " "<<RWCStemPrev[0]<<"\n";
          NumericVector Erootcrown;
          NumericVector psiRoot;
          NumericMatrix psiRhizo;
          NumericMatrix ElayersMat;
          NumericVector fittedE, dEdP;
          NumericVector psiLeaf;
          NumericMatrix newPsiStem;
          
          double Gwminc = Gwmin[c];
          
          if(!capacitance) {
            sFunction = supply[c];
            Erootcrown = sFunction["E"];
            psiRoot = sFunction["psiRoot"];
            psiRhizo = Rcpp::as<Rcpp::NumericMatrix>(sFunction["psiRhizo"]);
            psiLeaf = sFunction["psiLeaf"];
            newPsiStem = Rcpp::as<Rcpp::NumericMatrix>(sFunction["psiStem"]);
            ElayersMat = Rcpp::as<Rcpp::NumericMatrix>(sFunction["ERhizo"]);
            fittedE = sFunction["E"];
            dEdP = sFunction["dEdP"];
            //Set minimum conductance to zero to avoid large decreases in water potential to achieve a minimum flow 
            Gwminc = 0.0; 
          } else {
            List RSFunction = supply[c];
            Erootcrown = RSFunction["E"];
            psiRoot = RSFunction["psiRoot"];
            psiRhizo = Rcpp::as<Rcpp::NumericMatrix>(RSFunction["psiRhizo"]);
            ElayersMat = Rcpp::as<Rcpp::NumericMatrix>(RSFunction["ERhizo"]);
            sFunction = supplyFunctionAbovegroundCapacitance(Erootcrown, psiRoot,
                                                             psiStemPrev, PLCStemPrev, 
                                                             psiLeafPrev, 
                                                             VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                                             VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                                             Vsapwood[c], StemAF[c], StemPI0[c], StemEPS[c],
                                                             Vleaf[c], LeafAF[c], LeafPI0[c], LeafEPS[c],
                                                             tstep);
            fittedE = sFunction["E"];
            dEdP = sFunction["dEdP"];
            psiLeaf = sFunction["psiLeaf"];
            newPsiStem = Rcpp::as<Rcpp::NumericMatrix>(sFunction["psiStem"]);
          }
          
          
          if(fittedE.size()>0) {
            //Photosynthesis function for sunlit and shade leaves
            DataFrame photoSunlit = leafPhotosynthesisFunction(fittedE, Catm, Patm,Tcan[n], vpatm, 
                                                               zWind[c], 
                                                               absSWR_SL[c] + LWR_emmcan*LAI_SL(c,n), 
                                                               irradianceToPhotonFlux(absPAR_SL[c]), 
                                                               Vmax298SL, 
                                                               Jmax298SL, 
                                                               Gwminc, Gwmax[c], leafWidth[c], LAI_SL(c,n));
            DataFrame photoShade = leafPhotosynthesisFunction(fittedE, Catm, Patm,Tcan[n], vpatm, 
                                                              zWind[c], 
                                                              absSWR_SH[c] + LWR_emmcan*LAI_SH(c,n), 
                                                              irradianceToPhotonFlux(absPAR_SH[c]),
                                                              Vmax298SH, 
                                                              Jmax298SH, 
                                                              Gwminc, Gwmax[c], leafWidth[c], LAI_SH(c,n));
            
            NumericVector AnSunlit = photoSunlit["NetPhotosynthesis"];
            NumericVector AnShade = photoShade["NetPhotosynthesis"];
            NumericVector GwSunlit = photoSunlit["WaterVaporConductance"];
            NumericVector GwShade = photoShade["WaterVaporConductance"];
            NumericVector VPDSunlit = photoSunlit["LeafVPD"];
            NumericVector VPDShade = photoShade["LeafVPD"];
            NumericVector TempSunlit = photoSunlit["LeafTemperature"];
            NumericVector TempShade = photoShade["LeafTemperature"];
            
            //Profit maximization
            List PMSunlit = profitMaximization(sFunction, photoSunlit,  hydraulicCostFunction, Gwminc, Gwmax[c], VCstem_kmax[c]);
            List PMShade = profitMaximization(sFunction, photoShade,  hydraulicCostFunction, Gwminc,Gwmax[c], VCstem_kmax[c]);
            int iPMSunlit = PMSunlit["iMaxProfit"];
            int iPMShade = PMShade["iMaxProfit"];
            
            // Rcout<<iPMSunlit<<" "<<iPMShade<<"\n";
            //Get leaf status
            GW_SH(c,n)= GwShade[iPMShade];
            GW_SL(c,n)= GwSunlit[iPMSunlit];
            VPD_SH(c,n)= VPDShade[iPMShade];
            VPD_SL(c,n)= VPDSunlit[iPMSunlit];
            Temp_SH(c,n)= TempShade[iPMShade];
            Temp_SL(c,n)= TempSunlit[iPMSunlit];
            
            //Scale photosynthesis
            Aninst(c,n) = AnSunlit[iPMSunlit]*LAI_SL(c,n) + AnShade[iPMShade]*LAI_SH(c,n);
            photosynthesis[c] +=pow(10,-6)*12.01017*Aninst(c,n)*tstep; 
            
            //Average flow from sunlit and shade leaves
            double Eaverage = (fittedE[iPMSunlit]*LAI_SL(c,n) + fittedE[iPMShade]*LAI_SH(c,n))/(LAI_SL(c,n) + LAI_SH(c,n));
            
            
            //Find iPM for  flow corresponding to the  average flow
            double absDiff = 99999999.9;
            int iPM = -1;
            for(int k=0;k<fittedE.size();k++){ //Only check up to the size of fittedE
              double adk = std::abs(fittedE[k]-Eaverage);
              if(adk<absDiff) {
                absDiff = adk;
                iPM = k;
              }
            }
            
            //Calculate transpiration with capacitance effects
            if(iPM==-1) {
              Rcout<<"\n iPM -1! Eaverage="<< Eaverage << " fittedE.size= "<< fittedE.size()<<" iPMSunlit="<< iPMSunlit<< " fittedE[iPMSunlit]="<<fittedE[iPMSunlit]<<" iPMShade="<<iPMShade<<" fittedE[iPMShade]="<<fittedE[iPMShade]<<"\n";
              stop("");
            }
            
            //Scale water extracted from soil to cohort level
            NumericVector Esoilcn(nlayerscon[c],0.0);
            for(int lc=0;lc<nlayerscon[c];lc++) {
              Esoilcn[lc] = ElayersMat(iPM,lc)*0.001*0.01802*LAIphe[c]*tstep; //Scale from flow to water volume in the time step
            }
            
            //Scale from instantaneous flow to water volume in the time step
            Einst(c,n) = fittedE[iPM]*0.001*0.01802*LAIphe[c]*tstep; 
            
            //Store instantaneous total conductance
            dEdPinst(c,n) = dEdP[iPM];
            
            //Balance between extraction and 
            PWBinst(c,n) = sum(Esoilcn) - Einst(c,n);
            
            //Add to daily plant cohort transpiration
            Eplant[c] +=Einst(c,n);
            
            
            //Update symplastic storage and PLC
            psiRootVEC[c] = psiRoot[iPM]; 
            EinstVEC[c] = fittedE[iPM];
            PLC(c,n) = NA_REAL;
            RWCsteminst(c,n) = NA_REAL;
            psiLeafVEC[c] = psiLeaf[iPM];
            RWCsleafVEC[c] = symplasticRelativeWaterContent(psiLeafVEC[c], LeafPI0[c], LeafEPS[c]);
            psiStemMAT(c,_) = newPsiStem(iPM,_);
            for(int i=0;i<nStemSegments;i++) {
              if(!cavitationRefill) PLCstemMAT(c,i) = std::max(PLCstemMAT(c,i), 1.0 - apoplasticRelativeWaterContent(psiStemMAT(c,i), VCstem_c[c], VCstem_d[c]));
              else PLCstemMAT(c,i) = 1.0 - apoplasticRelativeWaterContent(psiStemMAT(c,i), VCstem_c[c], VCstem_d[c]);
              RWCsstemMAT(c,i) = symplasticRelativeWaterContent(psiStemMAT(c,i), StemPI0[c], StemEPS[c]);
            }
            
            
            //Copy transpiration and from connected layers to transpiration from soil layers
            //Copy psiRhizo and from connected layers to psiRhizo from soil layers
            int cl = 0;
            for(int l=0;l<nlayers;l++) {
              if(layerConnected(c,l)) {
                psiRhizoMAT(c,l) = psiRhizo(iPM,cl);
                SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
                soilLayerExtractInst(l,n) += Esoilcn[cl];
                cl++;
              } 
            }
          } else {
            if(verbose) Rcout<<"NS!";
          }
        } else { // If not connected to any soil layer
          if(capacitance) {
            //Sunlit photosynthesis
            double absRadSL = absSWR_SL[c] + LWR_emmcan*LAI_SL(c,n);
            double QSL = irradianceToPhotonFlux(absPAR_SL[c]);
            //Flow (cuticular conductance): iterative procedure to find flow, VPD and leaf temperature given Gwmin
            double ESLbk = std::max(0.0,1000.0*(Gwmin[c]*(meteoland::utils_saturationVP(std::max(0.0,Tcan[n]))-vpatm))/Patm);
            double leafTempSL = leafTemperature(absRadSL/ LAI_SL(c,n), Tcan[n], zWind[c], ESLbk, leafWidth[c]);
            double leafVPDSL = (meteoland::utils_saturationVP(std::max(0.0,leafTempSL))-vpatm);
            double ESL = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSL)/Patm);
            while(std::abs(ESL-ESLbk)> 0.000001) {
              leafTempSL = leafTemperature(absRadSL/ LAI_SL(c,n), Tcan[n], zWind[c], ESL, leafWidth[c]);
              leafVPDSL = meteoland::utils_saturationVP(std::max(0.0,leafTempSL))-vpatm;
              ESLbk = ESL;
              ESL = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSL)/Patm);
            }
            double AgSL = leafphotosynthesis(QSL, Catm, Gwmin[c]/1.6, std::max(0.0,leafTempSL), Vmax298SL, Jmax298SL);
            double AnSL = AgSL - 0.015*VmaxTemp(Vmax298SL, leafTempSL);
            
            //Shade photosynthesis
            double absRadSH = absSWR_SH[c] + LWR_emmcan*LAI_SH(c,n);
            double QSH = irradianceToPhotonFlux(absPAR_SH[c]);
            //Flow (cuticular conductance): iterative procedure to find flow, VPD and leaf temperature given Gwmin
            double ESHbk = std::max(0.0,1000.0*(Gwmin[c]*(meteoland::utils_saturationVP(std::max(0.0,Tcan[n]))-vpatm))/Patm);
            double leafTempSH = leafTemperature(absRadSH/ LAI_SH(c,n), Tcan[n], zWind[c], ESHbk, leafWidth[c]);
            double leafVPDSH = (meteoland::utils_saturationVP(std::max(0.0,leafTempSH))-vpatm);
            double ESH = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSH)/Patm);
            while(std::abs(ESH-ESHbk)> 0.000001) {
              leafTempSH = leafTemperature(absRadSH/ LAI_SH(c,n), Tcan[n], zWind[c], ESH, leafWidth[c]);
              leafVPDSH = meteoland::utils_saturationVP(std::max(0.0,leafTempSH))-vpatm;
              ESHbk = ESH;
              ESH = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSH)/Patm);
            }
            double AgSH = leafphotosynthesis(QSH, Catm, Gwmin[c]/1.6, std::max(0.0,leafTempSH), Vmax298SH, Jmax298SH);
            double AnSH = AgSH - 0.015*VmaxTemp(Vmax298SH, leafTempSH);
            
            //Average flow
            double Eaverage = (ESL*LAI_SL(c,n) + ESH*LAI_SH(c,n))/(LAI_SL(c,n) + LAI_SH(c,n));
            
            //Get leaf status
            GW_SH(c,n)= Gwmin[c];
            GW_SL(c,n)= Gwmin[c];
            VPD_SH(c,n)= leafVPDSH;
            VPD_SL(c,n)= leafVPDSL;
            Temp_SH(c,n)= leafTempSH;
            Temp_SL(c,n)= leafVPDSL;
            
            //Scale photosynthesis
            Aninst(c,n) = AnSL*LAI_SL(c,n) + AnSH*LAI_SH(c,n);
            photosynthesis[c] +=pow(10,-6)*12.01017*Aninst(c,n)*tstep; 
            
            
            Einst(c,n) = Eaverage*0.001*0.01802*LAIphe[c]*tstep; //Scale from instantaneous flow to water volume in the time step
            
            
            //Add to daily plant cohort transpiration
            Eplant[c] +=Einst(c,n);
            
            List sAb = E2psiAbovegroundCapacitanceDisconnected(Eaverage, 
                                                               psiStemPrev, PLCStemPrev, RWCStemPrev, 
                                                               psiLeafPrev, rwcsleafPrev,
                                                               VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                                                                                    VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                                                                                                                         Vsapwood[c], StemAF[c], StemPI0[c], StemEPS[c],
                                                                                                                                                                                    Vleaf[c], LeafAF[c], LeafPI0[c], LeafEPS[c],
                                                                                                                                                                                                                            klat, 
                                                                                                                                                                                                                            tstep);
            
            NumericVector newPsiStem = sAb["psiStem"];
            NumericVector newRWCsympstem = sAb["RWCsympstem"];
            
            
            
            //As it is disconnected, total conductance is equal to leaf hydraulic conductance
            dEdPinst(c,n) = sAb["kleaf"];
            
            //Update symplastic storage and PLC
            psiRootVEC[c] = newPsiStem[0];//Estimate of psiRoot = first stem segment
            psiLeafVEC[c] = sAb["psiLeaf"];
            RWCsleafVEC[c] = sAb["RWCsympleaf"];
            psiStemMAT(c,_) = newPsiStem;
            for(int i=0;i<nStemSegments;i++) {
              if(!cavitationRefill) PLCstemMAT(c,i) = std::max(PLCstemMAT(c,i), 1.0 - apoplasticRelativeWaterContent(psiStemMAT(c,i), VCstem_c[c], VCstem_d[c]));
              else PLCstemMAT(c,i) = 1.0 - apoplasticRelativeWaterContent(psiStemMAT(c,i), VCstem_c[c], VCstem_d[c]);
            }
            // Rcout<< "PLC "<< PLCstemMAT(c,0)<<"\n";
            RWCsstemMAT(c,_) = newRWCsympstem;
          } else {
            GW_SH(c,n)= NA_REAL;
            GW_SL(c,n)= NA_REAL;
            VPD_SH(c,n)= NA_REAL;
            VPD_SL(c,n)= NA_REAL;
            Temp_SH(c,n)= NA_REAL;
            Temp_SL(c,n)= NA_REAL;
            dEdPinst(c,n) = 0.0;
            Einst(c,n) = 0.0;
            Aninst(c,n) = 0.0;
          }
        }
      }


      //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
      PLC(c,n) = PLCstemMAT(c,nStemSegments-1);
      RWCsteminst(c,n) = RWCsstemMAT(c,nStemSegments-1)*(1.0 - StemAF[c]) + (1.0- PLCstemMAT(c,nStemSegments-1))*StemAF[c];
      RWCleafinst(c,n) = RWCsleafVEC[c]*(1.0 - LeafAF[c]) + apoplasticRelativeWaterContent(psiLeafVEC[c], VCleaf_c[c], VCleaf_d[c])*LeafAF[c];
      PsiSteminst(c,n) = psiStemMAT(c, nStemSegments-1); 
      
      PsiLeafinst(c,n) = psiLeafVEC[c]; //Store instantaneous leaf potential
      PsiRootinst(c,n) = psiRootVEC[c]; //Store instantaneous root potential

      //Store the minimum water potential of the day (i.e. mid-day)
      minPsiLeaf[c] = std::min(minPsiLeaf[c],PsiLeafinst(c,n));
      maxPsiLeaf[c] = std::max(maxPsiLeaf[c],PsiLeafinst(c,n));
      minPsiStem[c] = std::min(minPsiStem[c],PsiSteminst(c,n));
      minPsiRoot[c] = std::min(minPsiRoot[c],PsiRootinst(c,n));
      for(int l=0;l<nlayers;l++) {
        minPsiRhizo(c,l) = std::min(minPsiRhizo(c,l),psiRhizoMAT(c,l));
      }
    } //End of cohort loop
    

    
    //CANOPY AND SOIL ENERGY BALANCE
    //Proportion of the canopy exchanging LWR radiation  as the fraction of incoming LWR
    double canLWRexchprop = abs_LWR_can[n]/lwdr[n];
    // Rcout<<canLWRexchprop<<"\n";
    //Latent heat (PROBLEM: does not include interception)
    LEcan_heat[n] = pow(10.0,6.0)*meteoland::utils_latentHeatVaporisation(Tcan[n])*sum(Einst(_,n))/tstep; 
    //Canopy longwave emmission
    LWRcanout[n] = LWR_emmcan*canLWRexchprop;
    //Canopy convective heat exchange
    Hcan_heat[n] = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;
    //Soil-canopy turbulent heat exchange
    Hcansoil[n] = (meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG*(Tcan[n]-Tsoil[0]))/RAsoil;
    //Soil LWR emmission
    LWRsoilout[n] = emm_LWR_soil[n];
    //Soil conductivity
    // NumericVector lambda = layerthermalconductivity(sand, clay, W, Theta_FC);
    // double Ccansoil = lambda[0]*(Tcan[n]-Tsoil[0]);
    //Soil-canopy heat exchange
    LWRsoilcan[n] =  LWRsoilout[n]*canLWRexchprop;
    double G = LWRcanout[n] - LWRsoilcan[n] + Hcansoil[n]; //Only include a fraction equal to absorption
    //Canopy temperature changes
    Ebal[n] = abs_SWR_can[n]+abs_LWR_can[n] - LWRcanout[n] - LEcan_heat[n] - Hcan_heat[n] - G;
    double canopyThermalCapacity = 0.5*(0.8*LAIcellmax + 1.2*LAIcell)*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
    double Tcannext = Tcan[n]+ std::max(-3.0, std::min(3.0, tstep*Ebal[n]/canopyThermalCapacity)); //Avoids changes in temperature that are too fast
    if(n<(ntimesteps-1)) Tcan[n+1] = Tcannext;

    //Soil energy balance
    Ebalsoil[n] = abs_SWR_soil[n] + abs_LWR_soil[n] + LWRcanout[n] + Hcansoil[n] - LWRsoilout[n]; //Here we use all energy escaping to atmosphere
    //Soil temperature changes
    NumericVector soilTchange = soilTemperatureChange(dVec, Tsoil, sand, clay, W, Theta_FC, Ebalsoil[n]);
    for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tstep);
    if(n<(ntimesteps-1)) Tsoil_mat(n+1,_)= Tsoil;

    // Rcout<<n<<", Tatm: "<< Tatm[n]<< " Tcan: "<<Tcan[n]<<" soilT1 "<< Tsoil[0]<<"\n";
    

    //save canopy temperature
    canopyParams["Temp"] = Tcannext;
    
  } //End of timestep loop
  
  
  //Substract extracted water from soil moisture 
  NumericVector soilExtractionBalance(nlayers, 0.0);
  NumericVector soilHydraulicInput(nlayers, 0.0); //Water that entered into the layer across all time steps
  NumericVector soilHydraulicOutput(nlayers, 0.0);  //Water that left the layer across all time steps
  for(int l=0;l<nlayers;l++) {
    for(int n=0;n<ntimesteps;n++) {
      soilHydraulicInput[l] += (-1.0)*std::min(soilLayerExtractInst(l,n),0.0);
      soilHydraulicOutput[l] += std::max(soilLayerExtractInst(l,n),0.0);
    }
    soilExtractionBalance[l] = sum(SoilWaterExtract(_,l));
    W[l] = std::max(W[l]-(soilExtractionBalance[l]/Water_FC[l]),0.0);
  }
  
  //4z. Plant daily drought stress (from root collar mid-day water potential)
  NumericVector SoilExtractCoh(numCohorts,0.0);
  NumericVector DDS(numCohorts, 0.0);
  for(int c=0;c<numCohorts;c++) {
    SoilExtractCoh[c] =  sum(SoilWaterExtract(c,_));
    transpiration[c] = Eplant[c]; 
    PLCm[c] = sum(PLC(c,_))/((double)PLC.ncol());
    RWCsm[c] = sum(RWCsteminst(c,_))/((double)RWCsteminst.ncol());
    RWClm[c] = sum(RWCleafinst(c,_))/((double)RWCleafinst.ncol());
    dEdPm[c] = sum(dEdPinst(c,_))/((double)dEdPinst.ncol());  
    double maxConductance = maximumSoilPlantConductance(VGrhizo_kmax(c,_), VCroot_kmax(c,_), VCstem_kmax[c], VCleaf_kmax[c]);
    DDS[c] = Phe[c]*(1.0 - (dEdPm[c]/maxConductance));
  }
  
  //5. Soil evaporation
  double Rnground = 0.0;
  for(int n=0;n<ntimesteps;n++) {
    Rnground += (abs_SWR_soil[n] + abs_LWR_soil[n] - emm_LWR_soil[n])*tstep; //kJ·m-2 accumulate net soil radiation balance
  }
  Rnground = std::max(0.0,Rnground)/pow(10.0,3.0); //from kJ·m-2 to MJ·m-2
  double PETsoil = meteoland::penmanmonteith(200.0, elevation, tmin, tmax, rhmin, rhmax, Rnground, wind);
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  double Esoil = std::max(0.0,soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil));
  NumericVector EsoilVec(nlayers,0.0);
  //Exponential decay to divide bare soil evaporation among layers
  for(int l=0;l<nlayers;l++) {
    double cumAnt = 0.0;
    double cumPost = 0.0;
    for(int l2=0;l2<l;l2++) cumAnt +=dVec[l2];
    cumPost = cumAnt+dVec[l];
    if(l<(nlayers-1)) EsoilVec[l] = Esoil*(exp(-Ksoil*cumAnt)-exp(-Ksoil*cumPost));
    else EsoilVec[l] = Esoil*exp(-Ksoil*cumAnt);
    W[l] = std::max(W[l]-(EsoilVec[l]/Water_FC[l]),0.0);
  }
  psiVec = psi(soil, soilFunctions); //Update soil water potential
  
  //6. Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  LAIcohort.attr("names") = above.attr("row.names");
    
  DataFrame Tinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                      _["Tatm"] = Tatm, _["Tcan"] = Tcan, _["Tsoil"] = Tsoil_mat);
  DataFrame CEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["SWRcanin"] = abs_SWR_can, _["LWRcanin"] = abs_LWR_can,_["LWRcanout"] = LWRcanout, _["RAcan"] = RAcan,
                                        _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebal);
  DataFrame SEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = abs_SWR_soil, _["LWRsoilin"] = abs_LWR_soil,  _["LWRsoilout"] = LWRsoilout,
                                        _["Ebalsoil"] = Ebalsoil, _["RAsoil"] = RAsoil);
  List EB = List::create(_["Temperature"]=Tinst, _["CanopyEnergyBalance"] = CEBinst, _["SoilEnergyBalance"] = SEBinst);
  SWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GW_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GW_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LAI_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LAI_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  List AbsRadinst = List::create(_["SWR_SH"] = SWR_SH, _["SWR_SL"]=SWR_SL,
                                 _["LWR_SH"] = LWR_SH, _["LWR_SL"] = LWR_SL);
  NumericVector DB = NumericVector::create(_["Rain"] = rain,_["Snow"] = 0.0,_["NetRain"] = NetRain, _["Runon"] = runon, _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                                           _["SoilEvaporation"] = sum(EsoilVec), _["PlantExtraction"] = sum(soilExtractionBalance), _["Transpiration"] = sum(Eplant),
                                           _["HydraulicRedistribution"] = sum(soilHydraulicInput),
                                           _["LAIcell"] = LAIcell, _["LAIcelldead"] = LAIcelldead, _["Cm"] = Cm, _["Lground"] = LgroundPAR);
  
  DataFrame SB = DataFrame::create(_["SoilEvaporation"] = EsoilVec, 
                                   _["HydraulicInput"] = soilHydraulicInput, 
                                   _["HydraulicOutput"] = soilHydraulicOutput, 
                                   _["PlantExtraction"] = soilExtractionBalance, 
                                   _["psi"] = psiVec);
  
  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  dEdPinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PsiLeafinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PsiSteminst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PsiRootinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  RWCleafinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  RWCsteminst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PWBinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  List SoilInst = List::create(_["Extraction"] = soilLayerExtractInst);
  List PlantsInst = List::create(
                                 _["LAIsunlit"] = LAI_SL, _["LAIshade"] = LAI_SH, 
                                 _["AbsRad"] = AbsRadinst, _["E"]=Einst, _["An"]=Aninst,
                                 _["dEdPinst"] = dEdPinst,
                                 _["GWsunlit"] = GW_SL, _["GWshade"] = GW_SH,
                                 _["VPDsunlit"] = VPD_SL, _["VPDshade"] = VPD_SH,
                                 _["Tempsunlit"] = Temp_SL, _["Tempshade"] = Temp_SH,
                                 _["PsiRoot"] = PsiRootinst, 
                                 _["PsiStem"] = PsiSteminst, 
                                 _["PsiLeaf"] = PsiLeafinst, 
                                 _["PLCstem"] = PLC, 
                                 _["RWCstem"] = RWCsteminst,
                                 _["RWCleaf"] = RWCleafinst,
                                 _["PWB"] = PWBinst);
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                             _["Extraction"] = SoilExtractCoh,
                             _["Transpiration"] = Eplant, 
                             _["RootPsi"] = minPsiRoot, 
                             _["StemPsi"] = minPsiStem, 
                             _["StemPLC"] = PLCm, //Average daily stem PLC
                             _["LeafPsiMin"] = minPsiLeaf, 
                             _["LeafPsiMax"] = maxPsiLeaf, 
                             _["dEdP"] = dEdPm,//Average daily soilplant conductance
                             _["DDS"] = DDS, //Daily drought stress is the ratio of average soil plant conductance over its maximum value
                             _["StemRWC"] = RWCsm,
                             _["LeafRWC"] = RWClm);
  Plants.attr("row.names") = above.attr("row.names");
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["WaterBalance"] = DB, 
                        _["EnergyBalance"] = EB,
                        _["Soil"] = SB, 
                        _["SoilInst"] = SoilInst,
                        _["RhizoPsi"] = minPsiRhizo,
                        _["Plants"] = Plants,
                        _["PlantsInst"] = PlantsInst);
  l.attr("class") = CharacterVector::create("spwb.day","list");
  return(l);
}

// [[Rcpp::export("spwb.day")]]
List spwbDay(List x, List soil, CharacterVector date, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
            double latitude, double elevation, double slope, double aspect,  
            double prec, double runon=0.0) {
  //Control parameters
  List control = x["control"];
  bool verbose = control["verbose"];
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

  NumericVector ER = er(IntegerVector::create(doy));
  double er = ER[0];
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  int numCohorts = SP.size();
  
  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];
  double tmean = meteoland::utils_averageDaylightTemperature(tmin, tmax);
  List canopyParams = x["canopy"];
  double gddday = canopyParams["gdd"];
  if((tmean-5.0 < 0.0) & (doy>180)) {
    gddday = 0.0;
  } else if (doy<180){ //Only increase in the first part of the year
    if(tmean-5.0>0.0) gddday = gddday + (tmean-5.0);
  }
  canopyParams["gdd"] = gddday;
  //Update phenological status
  NumericVector phe = leafDevelopmentStatus(Sgdd, gddday);
  for(int j=0;j<numCohorts;j++) {
    LAI_dead[j] *= exp(-1.0*(wind/10.0)); //Decrease dead leaf area according to wind speed
    double LAI_exp_prev= LAI_expanded[j]; //Store previous value
    LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
    LAI_dead[j] += std::max(0.0, LAI_exp_prev-LAI_expanded[j]);//Check increase dead leaf area if expanded leaf area has decreased
  }

  List s;
  if(transpirationMode=="Simple") {
    s = spwbDay1(x,soil, tday, pet, prec, er, runon, rad, elevation, verbose);
  } else {
    if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; 
    if(wind<0.5) wind = 0.5; //Minimum windspeed abovecanopy
    s = spwbDay2(x,soil, tmin, tmax, rhmin, rhmax, rad, wind, latitude, elevation,
                solarConstant, delta, prec, er, runon, verbose);
  }
  // Rcout<<"hola4\n";
  return(s);
}

  
NumericVector getTrackSpeciesTranspiration( NumericVector trackSpecies, NumericVector Eplant, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector Eplantsp(nTrackSpecies, 0.0);
  NumericVector SP = x["SP"];
  int nCoh = SP.size();
  int ts;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        Eplantsp[its] += Eplant[i];
      }
    }
  }
  return(Eplantsp);
}

NumericVector getTrackSpeciesDDS(NumericVector trackSpecies, NumericVector DDS, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector DDSsp(nTrackSpecies, 0.0);
  NumericVector LAI = x["LAI"];
  NumericVector SP = x["SP"];
  int nCoh = LAI.size();
  int ts;
  double laiSum;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    laiSum = 0.0;
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        DDSsp[its] += DDS[i]*LAI[i];
        laiSum +=LAI[i];
      }
    }
    DDSsp = DDSsp/laiSum;
  }
  return(DDSsp);
}

IntegerVector order_vector(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}


// [[Rcpp::export(".spwbgridDay")]]
List spwbgridDay(CharacterVector lct, List xList, List soilList, 
                IntegerVector waterO, List queenNeigh, List waterQ,
                NumericVector tdayVec, NumericVector petVec, NumericVector rainVec, 
                NumericVector erVec, NumericVector radVec, NumericVector elevation,
                NumericVector trackSpecies, double patchsize) {
  int nX = xList.size();
  int nTrackSpecies = trackSpecies.size();
  NumericVector Rain(nX, NA_REAL), Snow(nX, NA_REAL), NetRain(nX,NA_REAL), Runon(nX,0.0), Infiltration(nX,NA_REAL);
  NumericVector Runoff(nX,NA_REAL), DeepDrainage(nX,NA_REAL);
  NumericVector Esoil(nX,NA_REAL), Eplant(nX,NA_REAL);
  NumericMatrix Transpiration(nX, nTrackSpecies), DDS(nX, nTrackSpecies);
  double runoffExport = 0.0;
  
  //A. Subsurface fluxes
  double cellArea = patchsize; //cell size in m2
  double cellWidth = sqrt(patchsize); //cell width in m
  double n = 3.0;
  double K = 7.2; //7.2 m/day
  //1. Calculate water table depth
  NumericVector WTD(nX,NA_REAL); //Water table depth
  NumericVector WaterTableElevation(nX,NA_REAL); //water table elevation (including cell elevation) in meters
  for(int i=0;i<nX;i++){
    if((lct[i]=="wildland") || (lct[i]=="agriculture") ) {
      List x = Rcpp::as<Rcpp::List>(xList[i]);
      List soil = Rcpp::as<Rcpp::List>(soilList[i]);
      List control = x["control"];
      WTD[i] = waterTableDepth(soil, control["soilFunctions"]);
      WaterTableElevation[i] = elevation[i]-(WTD[i]/1000.0);
    }
  }
  //2. Calculate inflow/outflow for each cell (in m3/day)
  NumericVector inflow(nX, 0.0); 
  NumericVector outflow(nX, 0.0); 
  for(int i=0;i<nX;i++){
    if((lct[i]=="wildland") || (lct[i]=="agriculture") ) {
      List soil = Rcpp::as<Rcpp::List>(soilList[i]);
      double D = soil["SoilDepth"]; //Soil depth in mm
      if(WTD[i]<D) {
        double T = ((K*D*0.001)/n)*pow(1.0-(WTD[i]/D),n); //Transmissivity in m2
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[i]);
        //water table slope between target and neighbours
        for(int j=0;j<ni.size();j++) {
          double tanBeta = (WaterTableElevation[i]-WaterTableElevation[ni[j]-1])/cellWidth;
          if(tanBeta>0.0) {
            if((lct[ni[j]-1]=="wildland") || (lct[ni[j]-1]=="agriculture")) { //Only flows to other wildland or agriculture cells
              double qn = tanBeta*T*cellWidth; //flow in m3
              inflow[ni[j]-1] += qn;
              outflow[i] += qn;
            }
          }
        }
      } 
    }
  }
  //3. Apply changes in soil moisture to each cell 
  for(int i=0;i<nX;i++){
    if((lct[i]=="wildland") || (lct[i]=="agriculture") ) {
      double deltaS = 1000.0*((inflow[i]-outflow[i])/cellArea); //change in moisture in mm (L/m2)
      if(deltaS != 0.0) {
        // Rcout<<deltaS<<"_";
        List x = Rcpp::as<Rcpp::List>(xList[i]);
        List soil = Rcpp::as<Rcpp::List>(soilList[i]);
        NumericVector W = soil["W"]; //Access to soil state variable
        NumericVector dVec = soil["dVec"];
        NumericVector macro = soil["macro"];
        NumericVector rfc = soil["rfc"];
        List control = x["control"];
        String soilFunctions = control["soilFunctions"];
        NumericVector Water_FC = waterFC(soil, soilFunctions);
        NumericVector Water_SAT = waterSAT(soil, soilFunctions);
        int nlayers = dVec.length();
        for(int l=(nlayers-1);l>=0;l--) {
          if(dVec[l]>0) {
            double Wn = W[l]*Water_FC[l] + deltaS; //Update water volume
            deltaS = std::max(Wn - Water_SAT[l],0.0); //Update deltaS, using the excess of water over saturation
            W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
          }
        }
        if(deltaS>0) { //If soil is completely saturated increase Runon (return flow) to be processed with vertical flows
          Runon[i] += deltaS;
        }
        // Rcout<<WTD[i]<<"/"<<waterTableDepth(soil,soilFunctions)<<"\n";
      }
    }
  }
  Rcout<<"\n";
  
  //B. Vertical and surface fluxes
  for(int i=0;i<nX;i++) {
    //get next cell in order
    int iCell = waterO[i]-1; //Decrease index!!!!
    if((lct[iCell]=="wildland") || (lct[iCell]=="agriculture") ) {
      List x = Rcpp::as<Rcpp::List>(xList[iCell]);
      List soil = Rcpp::as<Rcpp::List>(soilList[iCell]);
      //Run daily soil water balance for the current cell
      List res = spwbDay1(x, soil, tdayVec[iCell], petVec[iCell], rainVec[iCell], erVec[iCell], 
                          Runon[iCell], radVec[iCell], elevation[iCell]);
      List DB = res["WaterBalance"];
      List SB = res["Soil"];
      List PL = res["Plants"];
      Snow[iCell] = DB["Snow"];
      Rain[iCell] = DB["Rain"];
      NetRain[iCell] = DB["NetRain"];
      Runon[iCell] = DB["Runon"];
      Infiltration[iCell] = DB["Infiltration"];
      Runoff[iCell] = DB["Runoff"];
      DeepDrainage[iCell] = DB["DeepDrainage"];
      Esoil[iCell] = sum(Rcpp::as<Rcpp::NumericVector>(SB["SoilEvaporation"]));
      NumericVector EplantCoh = Rcpp::as<Rcpp::NumericVector>(PL["Transpiration"]);
      NumericVector DDScell = PL["DDS"];
      Eplant[iCell] = sum(EplantCoh);
      if(nTrackSpecies>0) {
        Transpiration(iCell,_) = getTrackSpeciesTranspiration(trackSpecies, EplantCoh, x);
        DDS(iCell,_) = getTrackSpeciesDDS(trackSpecies, DDScell, x);
      }

      //Assign runoff to runon of neighbours
      double ri =  Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri; //decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="rock") {//all Precipitation becomes surface runoff if cell is rock outcrop
      Runoff[iCell] =  Runon[iCell]+rainVec[iCell];
      double ri = Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri;//decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="static") {
      // static cells receive water from other cells or Precipitation
      // but do not export to the atmosphere contribute nor to other cells.
      // Hence, water balance over the landscape is achieved by
      // adding this water to the landscape export via landscape runoff.
      runoffExport += Runon[iCell] + rainVec[iCell];
    }
  }

    
  DataFrame waterBalance = DataFrame::create(_["Rain"] = Rain, _["Snow"] = Snow, _["NetRain"] = NetRain, _["Runon"] = Runon, _["Infiltration"] = Infiltration,
                                   _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                                   _["Esoil"] = Esoil, _["Eplant"] = Eplant);
  return(List::create(_["WaterBalance"] = waterBalance,
                      _["RunoffExport"] = runoffExport,
                      _["Transpiration"] = Transpiration,
                      _["DDS"] = DDS));
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
  if(transpirationMode=="Complex"){
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
  if(transpirationMode=="Simple") {
    if(!paramsTransp.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("WUE")) stop("WUE missing in spwbInput$paramsTransp");
  } else if(transpirationMode=="Complex") {
    if(!paramsTransp.containsElementNamed("VCstem_kmax")) stop("VCstem_kmax missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCstem_c")) stop("VCstem_c missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCstem_d")) stop("VCstem_d missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_c")) stop("VCroot_c missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("VCroot_d")) stop("VCroot_d missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Gwmax")) stop("Gwmax missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Vmax298")) stop("Vmax298 missing in spwbInput$paramsTransp");
    if(!paramsTransp.containsElementNamed("Jmax298")) stop("Jmax298 missing in spwbInput$paramsTransp");
  }
  if(transpirationMode=="Complex") {
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

//
// [[Rcpp::export("spwb.resetInputs")]]
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
    if(transpirationMode=="Complex") {
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
  checkspwbInput(x, soil, transpirationMode, soilFunctions);
  
  //Store input
  List spwbInput = clone(x);
  List soilInput = clone(soil);
    
  //Meteorological input    
  NumericVector MinTemperature, MaxTemperature;
  NumericVector MinRelativeHumidity, MaxRelativeHumidity;
  NumericVector PET;
  NumericVector Radiation, WindSpeed;
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  int numDays = Precipitation.size();
  if(transpirationMode=="Simple") {
    PET = meteo["PET"];
    if(control["snowpack"]) Radiation = meteo["Radiation"];
  } else if(transpirationMode=="Complex") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    MinTemperature = meteo["MinTemperature"];
    MaxTemperature = meteo["MaxTemperature"];
    MinRelativeHumidity = meteo["MinRelativeHumidity"];
    MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    Radiation = meteo["Radiation"];
    WindSpeed = meteo["WindSpeed"];
    PET = NumericVector(numDays);
  }
  CharacterVector dateStrings = meteo.attr("row.names");
  
  IntegerVector DOY = date2doy(dateStrings);
  
  //Canpopy parameters
  List canopyParams = x["canopy"];
  
  //Calculate GDD (taking into account initial sum)
  NumericVector GDD = gdd(DOY, MeanTemperature, 5.0, canopyParams["gdd"]);
  
  NumericVector ER = er(DOY);
  
  
  //Plant input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  NumericVector SP = cohorts["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  int numCohorts = SP.size();
  

  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector Sgdd = paramsBase["Sgdd"];


  //Soil input
  NumericVector W = soil["W"];
  int nlayers = W.size();
  NumericVector Water_FC = waterFC(soil, soilFunctions);

  //Detailed subday results
  List subdailyRes(numDays);
    
  //Water balance output variables
  NumericVector LAIcell(numDays),LAIcelldead(numDays);
  NumericVector Cm(numDays);
  NumericVector Lground(numDays);
  NumericVector Runoff(numDays);
  NumericVector Rain(numDays);
  NumericVector Snow(numDays);
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
  
  
  Wdays(0,_) = W;

  if(verbose) {
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"i:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
  }

  if(verbose) Rcout << "Daily balance:";
  NumericVector Eplanttot(numDays,0.0);
  List s;
  for(int i=0;i<numDays;i++) {
      if(verbose) Rcout<<".";//<<i;
      canopyParams["gdd"] = GDD[i];
      double wind = WindSpeed[i];
      if(NumericVector::is_na(wind)) wind = control["defaultWindSpeed"]; //Default 1 m/s -> 10% of fall every day
      if(wind<0.5) wind = 0.5; //Minimum windspeed abovecanopy
      
      
      //1. Phenology and leaf fall
      NumericVector phe = leafDevelopmentStatus(Sgdd, GDD[i]);
      for(int j=0;j<numCohorts;j++) {
        LAI_dead[j] *= exp(-1.0*(wind/10.0)); //Decrease dead leaf area according to wind speed
        double LAI_exp_prev=0.0;
        if(i>0) LAI_exp_prev= LAI_expanded[j]; //Store previous value
        LAI_expanded[j] = LAI_live[j]*phe[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
        LAI_dead[j] += std::max(0.0, LAI_exp_prev-LAI_expanded[j]);//Check increase dead leaf area if expanded leaf area has decreased
        PlantLAI(i,j) = LAI_expanded[j];
      }
      
      
      //2. Water balance and photosynthesis
      if(transpirationMode=="Simple") {
        s = spwbDay1(x, soil, MeanTemperature[i], PET[i], Precipitation[i], ER[i], 0.0, 
                     Radiation[i], elevation, verbose); //No Runon in simulations for a single cell
      } else if(transpirationMode=="Complex") {
        int ntimesteps = control["ndailysteps"];
        double tstep = 86400.0/((double) ntimesteps);
        std::string c = as<std::string>(dateStrings[i]);
        int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
        double delta = meteoland::radiation_solarDeclination(J);
        double solarConstant = meteoland::radiation_solarConstant(J);
        s = spwbDay2(x, soil, MinTemperature[i], MaxTemperature[i], 
                         MinRelativeHumidity[i], MaxRelativeHumidity[i], Radiation[i], wind, 
                         latitude, elevation, solarConstant, delta, Precipitation[i], ER[i], 0.0, verbose);
        List Plants = Rcpp::as<Rcpp::List>(s["Plants"]);
        List PlantsInst = Rcpp::as<Rcpp::List>(s["PlantsInst"]);
        List AbsRadinst = Rcpp::as<Rcpp::List>(PlantsInst["AbsRad"]);
        
        NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(AbsRadinst["SWR_SL"]);
        NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(AbsRadinst["SWR_SH"]);
        NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(AbsRadinst["LWR_SL"]);
        NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(AbsRadinst["LWR_SH"]);
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
      NetRain[i] = db["NetRain"];
      PlantExtraction[i] = db["PlantExtraction"];
      if(transpirationMode=="Complex")  {
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
      NumericVector EplantCoh = Plants["Transpiration"];
      Eplantdays(i,_) = EplantVec;
      PlantPhotosynthesis(i,_) = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
      PlantTranspiration(i,_) = EplantCoh;
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["DDS"]);
      if(transpirationMode=="Complex") {
        NumericVector HydrInVec = sb["HydraulicInput"];
        HydrIndays(i,_) = HydrInVec;
        LeafPsiMin(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMin"]);
        LeafPsiMax(i,_) = Rcpp::as<Rcpp::NumericVector>(Plants["LeafPsiMax"]);
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
      if(transpirationMode=="Complex"){
        StemRWC(i,_) = as<Rcpp::NumericVector>(Plants["StemRWC"]);
        LeafRWC(i,_) = as<Rcpp::NumericVector>(Plants["LeafRWC"]); 
      }
      
      if(subdailyResults) {
        subdailyRes[i] = clone(s);
      }
      if(i<(numDays-1)) Wdays(i+1,_) = W;
      WaterTable[i] = waterTableDepth(soil, soilFunctions);
  }
  if(verbose) Rcout << "done\n";
  
  for(int l=0;l<nlayers;l++) {
    MLdays(_,l) = Wdays(_,l)*Water_FC[l]; 
    MLTot = MLTot + MLdays(_,l);
  }
  NumericVector Evapotranspiration = Transpiration+SoilEvaporation + Interception;
  
  if(verbose) {
    double Precipitationsum = sum(Precipitation);
    double NetRainsum = sum(NetRain);
    double Interceptionsum = sum(Interception);
    double SoilEvaporationsum = sum(SoilEvaporation);
    double Runoffsum  = sum(Runoff);
    double Infiltrationsum  = sum(Infiltration);
    double DeepDrainagesum = sum(DeepDrainage);
    double Transpirationsum = sum(Transpiration);
    
    Rcout<<"Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    Rcout<<"Rain (mm) "  <<round(sum(Rain)) <<" Snow (mm) "  <<round(sum(Snow)) <<"\n";
    Rcout<<"Interception (mm) " << round(Interceptionsum)  <<" Net rainfall (mm) " << round(NetRainsum) <<"\n";
    Rcout<<"Infiltration (mm) " << round(Infiltrationsum)  <<
      " Runoff (mm) " << round(Runoffsum) <<
        " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
    Rcout<<"Soil evaporation (mm) " << round(SoilEvaporationsum);
    Rcout<<" Transpiration (mm) "  <<round(Transpirationsum) <<"\n";
    if(transpirationMode =="Complex") {
      Rcout<<"Plant extraction from soil (mm) " << round(sum(PlantExtraction));
      Rcout<<" Hydraulic redistribution (mm) " << round(sum(HydraulicRedistribution)) <<"\n";
    }
    for(int l=0;l<nlayers;l++) Rcout << "W"<<(l+1)<<"f:"<< round(100*W[l])/100<<" ";
    Rcout<<"\n";
    Rcout<<"Final soil water content (mm): "<< round(MLTot[numDays-1])<<"\n\n";
    
  }
  if(verbose) Rcout<<"Building SPWB output ...";
  
   DataFrame SWB;
   if(transpirationMode=="Simple") {
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
   DataFrame DWB;
   if(transpirationMode=="Simple") {
     DWB = DataFrame::create(_["GDD"] = GDD,
                             _["LAIcell"]=LAIcell, _["LAIcelldead"] = LAIcelldead,  _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow,
                             _["NetRain"]=NetRain,_["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration);
   } else {
     DWB = DataFrame::create(_["GDD"] = GDD,
                             _["LAIcell"]=LAIcell, _["LAIcelldead"] = LAIcelldead,  _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, 
                             _["Precipitation"] = Precipitation, _["Rain"] = Rain, _["Snow"] = Snow,
                             _["NetRain"]=NetRain,_["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                             _["Evapotranspiration"]=Evapotranspiration,_["Interception"] = Interception, _["SoilEvaporation"]=SoilEvaporation,
                             _["PlantExtraction"] = PlantExtraction, _["Transpiration"]=Transpiration, 
                             _["HydraulicRedistribution"] = HydraulicRedistribution);
     
   }
  
   SWB.attr("row.names") = meteo.attr("row.names") ;
   DWB.attr("row.names") = meteo.attr("row.names") ;

  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafRWC.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  dEdP.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMin.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  LeafPsiMax.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  StemPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  RootPsi.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantPhotosynthesis.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsSWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantAbsLWR.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  PlantLAI.attr("dimnames") = List::create(meteo.attr("row.names"), above.attr("row.names")) ;
  
  DataFrame DEB = DataFrame::create(_["SWRcanin"] = SWRcanin, _["LWRcanin"] = LWRcanin, _["LWRcanout"] = LWRcanout,
                                    _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["LWRsoilcan"] = LWRsoilcan, _["Ebalcan"] = Ebalcan, 
                                    _["Hcansoil"] = Hcansoil, _["SWRsoilin"] = SWRsoilin, _["LWRsoilin"] = LWRsoilin, _["LWRsoilout"] = LWRsoilout,
                                    _["Ebalsoil"] = Ebalsoil, _["RAcan"] = RAcan, _["RAsoil"] = RAsoil);  
  DataFrame DT = DataFrame::create(_["Tatm_mean"] = Tatm_mean, _["Tatm_min"] = Tatm_min, _["Tatm_max"] = Tatm_max,
                                   _["Tcan_mean"] = Tcan_mean, _["Tcan_min"] = Tcan_min, _["Tcan_max"] = Tcan_max,
                                     _["Tsoil_mean"] = Tsoil_mean, _["Tsoil_min"] = Tsoil_min, _["Tsoil_max"] = Tsoil_max);
  DEB.attr("row.names") = meteo.attr("row.names") ;
  DT.attr("row.names") = meteo.attr("row.names") ;
  subdailyRes.attr("names") = meteo.attr("row.names") ;
  List l;
  if(transpirationMode=="Simple") {
    l = List::create(Named("spwbInput") = spwbInput,
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
    CharacterVector ln = CharacterVector(22);
    l = List(22);
    l[0] = spwbInput;
    ln[0] = "spwbInput";
    l[1] = soilInput;
    ln[1] = "soilInput";
    l[2] = DWB;
    ln[2] = "WaterBalance";
    l[3] = SWB;
    ln[3] = "Soil";
    l[4] = DEB;
    ln[4] = "EnergyBalance";
    l[5] = DT;
    ln[5] = "Temperature";
    l[6] = PlantLAI;
    ln[6] = "PlantLAI";
    l[7] = PlantAbsSWR;
    ln[7] = "PlantAbsorbedSWR";
    l[8] = PlantAbsLWR;
    ln[8] = "PlantAbsorbedLWR";
    l[9] = PlantTranspiration;
    ln[9] = "PlantTranspiration";
    l[10] = PlantPhotosynthesis;
    ln[10] = "PlantPhotosynthesis";
    l[11] = dEdP;
    ln[11] = "dEdP";
    l[12] = LeafPsiMin;
    ln[12] = "LeafPsiMin";
    l[13] = LeafPsiMax;
    ln[13] = "LeafPsiMax";
    l[14] = LeafRWC;
    ln[14] = "LeafRWC";
    l[15] = StemPsi;
    ln[15] = "StemPsi";
    l[16] = StemPLC;
    ln[16] = "StemPLC";
    l[17] = StemRWC;
    ln[17] = "StemRWC";
    l[18] = RootPsi;
    ln[18] = "RootPsi";
    l[19] = RhizoPsi;
    ln[19] = "RhizoPsi";
    l[20] = PlantStress;
    ln[20] = "PlantStress";
    l[21] = subdailyRes;
    ln[21] = "subdaily";
    l.attr("names") = ln;
  }
  l.attr("class") = CharacterVector::create("spwb","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}
