#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "carbon.h"
#include "photosynthesis.h"
#include "root.h"
#include "soil.h"
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1
const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double eps_xylem = 1e3; // xylem elastic modulus (1 GPa = 1000 MPa)

//Returns the average soil moisture within the rhizosphere of each cohort
NumericMatrix cohortRhizosphereMoisture(NumericMatrix W, List RHOP) {
  int numCohorts = W.nrow();
  int nlayers = W.ncol();
  NumericMatrix CRM(numCohorts, nlayers);
  CRM.attr("dimnames") = W.attr("dimnames");
  
  for(int coh=0;coh<numCohorts;coh++) {
    NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[coh]);
    for(int l=0;l<nlayers;l++) {
      CRM(coh,l) = 0.0;
      for(int c=0;c<numCohorts;c++) {
        CRM(coh,l) += RHOPcoh(c,l)*W(c,l);
      }
    }
  }
  return(CRM);
}

void rhizosphereMoistureExtraction(NumericMatrix cohExtract, 
                                   NumericVector WaterFC,
                                   NumericMatrix Wpool, 
                                   List RHOP,
                                   NumericVector poolProportions) {
  int numCohorts = Wpool.nrow();
  int nlayers = Wpool.ncol();
  for(int coh=0;coh<numCohorts;coh++) {
    NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[coh]);
    for(int c=0;c<numCohorts;c++) {
      for(int l=0;l<nlayers;l++) {
        Wpool(c,l) -= (RHOPcoh(c,l)*cohExtract(coh,l)/(WaterFC[l]*poolProportions[c]));
      }
    }
  }
}


// [[Rcpp::export("transp_profitMaximization")]]
List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax, 
                        double gainModifier = 1.0, double costModifier = 1.0, String costWater = "dEdP") {
  NumericVector supplyE = supplyFunction["E"];
  NumericVector supplydEdp = supplyFunction["dEdP"];
  NumericVector Ag = photosynthesisFunction["GrossPhotosynthesis"];
  NumericVector leafTemp = photosynthesisFunction["LeafTemperature"];
  NumericVector leafVPD = photosynthesisFunction["LeafVPD"];
  NumericVector Gsw = photosynthesisFunction["Gsw"];
  NumericVector supplyKterm = supplyFunction["kterm"];
  int nsteps = supplydEdp.size();
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  double maxKterm = 0.0, minKterm = 99999999.0;
  double Agmax = 0.0;
  //Find valid limits according to stomatal conductance
  int ini = 0, fin = nsteps-1;
 
  // mindEdp = 0.0; 
  // int imaxdEdp = 0;
  // maxdEdp = supplydEdp[0];
  for(int i=ini;i<fin;i++) {
    // if(supplydEdp[i] > maxdEdp) {
    //   maxdEdp = supplydEdp[i];
    //   imaxdEdp = i;
    // }
    if(costWater=="dEdP") {
      mindEdp = std::min(mindEdp, supplydEdp[i]);
      maxdEdp = std::max(maxdEdp, supplydEdp[i]);
    } else {
      minKterm = std::min(minKterm, supplyKterm[i]);
      maxKterm = std::max(maxKterm, supplyKterm[i]);
    }
    Agmax = std::max(Agmax, Ag[i]);
  }
  
  //Evaluate profit for valid steps
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  for(int i=ini;i<fin;i++) {
    gain[i] = pow(Ag[i]/Agmax, gainModifier);
    if(costWater=="dEdP") {
      cost[i] = pow((maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp), costModifier); 
    }  else {
      cost[i] = pow((maxKterm-supplyKterm[i])/(maxKterm-minKterm), costModifier);
    }
    profit[i] = gain[i]-cost[i];
  }
  
  while((Gsw[ini]<=Gswmin) && (ini<fin)) ini++;
  while((Gsw[fin]>=Gswmax) && (fin>ini)) fin--; 
  
  //Ensure that ini <=fin
  ini = std::min(ini, fin);
  fin = std::max(ini,fin);
  
  int imaxprofit=ini;
  double maxprofit=profit[ini];
  if(fin>ini) {
    for(int i=ini+1;i<=fin;i++){
      if((profit[i]>maxprofit)) {
        maxprofit = profit[i];
        imaxprofit = i;
      }
    }
  }
  // Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
  if((Gsw[imaxprofit] > Gswmax) && (imaxprofit>ini)) {
    Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
    for(int i=0;i<Gsw.size();i++) {
      Rcout<< i << " Gsw "<< Gsw[i] << " supplyE "<< supplyE[i] << " leafT "<< leafTemp[i]<< " leafVPD "<< leafVPD[i]  << "\n";
    }
    stop("Gsw > Gswmax");
  }
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}


List transpirationSperry(List x, double tmin, double tmax, 
                         double tminPrev, double tmaxPrev, double tminNext, 
                         double rhmin, double rhmax, double rad, double wind, 
                  double latitude, double elevation, double slope, double aspect, 
                  double solarConstant, double delta, double prec,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  List numericParams = control["numericParams"];
  int ntrial = numericParams["ntrial"];
  int maxNsteps  = numericParams["maxNsteps"];
  double psiTol = numericParams["psiTol"];
  double ETol = numericParams["ETol"];
  bool capacitance = control["capacitance"];
  bool cochard = control["cochard"];
  String cavitationRefill = control["cavitationRefill"];
  double refillMaximumRate = control["refillMaximumRate"];
  double klatleaf = control["klatleaf"];
  double klatstem = control["klatstem"];
  int ntimesteps = control["ndailysteps"];
  int nsubsteps = control["nsubsteps"];
  String costWater = control["costWater"];
  double costModifier = control["costModifier"];
  double gainModifier = control["gainModifier"];
  bool plantWaterPools = control["plantWaterPools"];
  double verticalLayerSize = control["verticalLayerSize"];
  double windMeasurementHeight  = control["windMeasurementHeight"];
  double thermalCapacityLAI = control["thermalCapacityLAI"];
  bool multiLayerBalance = control["multiLayerBalance"];
  double defaultWindSpeed = control["defaultWindSpeed"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  NumericVector N = Rcpp::as<Rcpp::NumericVector>(above["N"]);
  
  int numCohorts = LAIlive.size();
  
  //Soil input
  List soil = x["soil"];
  NumericVector dVec = soil["dVec"];
  int nlayers = dVec.length();
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector Ws = soil["W"]; //Access to soil state variable
  double SWE = soil["SWE"];
  NumericVector psiSoil = psi(soil, soilFunctions); //Get soil water potential
  
  //Canopy params
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  NumericVector zlow = canopyParams["zlow"];
  NumericVector zmid = canopyParams["zmid"];
  NumericVector zup = canopyParams["zup"];
  NumericVector Tair = canopyParams["Tair"];
  NumericVector VPair = canopyParams["VPair"];
  NumericVector Cair = canopyParams["Cair"];
  int ncanlayers = Tair.size(); //Number of canopy layers
  
  //Root distribution input
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  
  //Water pools
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  NumericMatrix Wrhizo;
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  
  //Base parameters
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector alphaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["alphaSWR"]);
  NumericVector gammaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["gammaSWR"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmax"]);
  NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Plant_kmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  NumericVector VCroot_kmax_sum = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCroot_kmax"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector Vmax298 = paramsTransp["Vmax298"];
  NumericVector Jmax298 = paramsTransp["Jmax298"];

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  

  //Comunication with outside
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector Stem1PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
  NumericVector Stem2PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem2Psi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector NSPLVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["NSPL"]);
  
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double latrad = latitude * (M_PI/180.0);
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);

  //Atmospheric pressure
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  
  double Catm = control["Catm"];
  
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  //If canopy VP is missing or not multilayer initiate it to vpatm
  if(NumericVector::is_na(VPair[0]) || (!multiLayerBalance)){
    for(int i=0;i<ncanlayers;i++) VPair[i] = vpatm;
  }
  //Daily cloud cover
  double cloudcover = 0.0;
  if(prec >0.0) cloudcover = 1.0;
  bool clearday = (prec==0);
  
  
  //1. Leaf Phenology: Adjusted leaf area index
  NumericVector Phe(numCohorts);
  double LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0, canopyHeight = 0.0;
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    if(LAIlive[c]==0.0) Phe[c] = 0.0;
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAIphe[c];
    if((canopyHeight<H[c]) & ((LAIphe[c]+LAIdead[c])>0.0)) canopyHeight = H[c];
  }
  //Create z vector with all layer height limits
  NumericVector z(ncanlayers+1,0.0);
  for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
  //LAI distribution per layer and cohort
  NumericMatrix LAIme = LAIdistributionVectors(z, LAIphe, H, CR); //Expanded leaves
  NumericMatrix LAImd = LAIdistributionVectors(z, LAIdead, H, CR); //Dead (standing) leaves
  NumericMatrix LAImx = LAIdistributionVectors(z, LAIlive, H, CR); //Maximum leaf expansion
  //LAI profile per layer
  NumericVector LAIpx = LAIprofileVectors(z, LAIlive, H, CR);
  NumericVector LAIpe = LAIprofileVectors(z, LAIphe, H, CR);
  NumericVector LAIpd = LAIprofileVectors(z, LAIdead, H, CR);
  NumericVector lad = 100.0*(LAIpe + LAIpd)/verticalLayerSize;
  //3. Wind extinction profile
  if(NumericVector::is_na(wind)) wind = defaultWindSpeed; //set to default if missing
  DataFrame canopyTurbulence = NA_REAL;
  NumericVector zWind(ncanlayers,wind), dU(ncanlayers, 0.0), uw(ncanlayers, 0.0);
  if(canopyHeight>0.0) {
    canopyTurbulence = windCanopyTurbulence(zmid, lad,  canopyHeight, 
                                                      wind, windMeasurementHeight);
    zWind = canopyTurbulence["u"]; 
    dU = Rcpp::as<Rcpp::NumericVector>(canopyTurbulence["du"]);
    uw = canopyTurbulence["uw"];
  } 
  //4a. Instantaneous direct and diffuse shorwave radiation
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, slorad, asprad, delta,
                                                        rad, clearday, ntimesteps);
  NumericVector solarHour = ddd["SolarHour"]; //in radians
  
  //4b. Instantaneous air temperature (above canopy) and longwave radiation
  NumericVector Tatm(ntimesteps), lwdr(ntimesteps), Tcan(ntimesteps, NA_REAL), Tsunrise(ntimesteps);
  NumericVector net_LWR_can(ntimesteps),LEcan_heat(ntimesteps), Hcan_heat(ntimesteps), Ebal(ntimesteps);
  NumericVector net_LWR_soil(ntimesteps), Ebalsoil(ntimesteps), Hcansoil(ntimesteps), LEsoil_heat(ntimesteps);
  NumericMatrix Tsoil_mat(ntimesteps, nlayers);
  NumericMatrix Tcan_mat(ntimesteps, ncanlayers);
  NumericMatrix VPcan_mat(ntimesteps, ncanlayers);
  //Daylength in seconds (assuming flat area because we want to model air temperature variation)
  double tauday = meteoland::radiation_daylengthseconds(latrad,0.0,0.0, delta); 
  for(int n=0;n<ntimesteps;n++) {
    //From solar hour (radians) to seconds from sunrise
    Tsunrise[n] = (solarHour[n]*43200.0/M_PI)+ (tauday/2.0) +(tstep/2.0); 
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern(Tsunrise[n], tmin, tmax, tminPrev, tmaxPrev, tminNext, tauday);
    //Longwave sky diffuse radiation (W/m2)
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpatm, cloudcover);
  }
  if(NumericVector::is_na(Tair[0])) {//If missing initialize canopy profile with atmospheric air temperature 
    for(int i=0;i<ncanlayers;i++) Tair[i] = Tatm[0];
  }
  if(NumericVector::is_na(Tsoil[0])) {//If missing initialize soil temperature with atmospheric air temperature 
    for(int l=0;l<nlayers; l++) Tsoil[l] = Tatm[0];
  }
  //Take initial canopy air temperature from previous day
  Tcan[0] = sum(Tair*LAIpx)/sum(LAIpx);
  for(int j=0;j<ncanlayers; j++) {
    Tcan_mat(0,j) = Tair[j];
    VPcan_mat(0,j) = VPair[j];
  }
  //Take temperature soil vector 
  Tsoil_mat(0,_) = Tsoil; 
  
  
  
  //4c. Light extinction and absortion by time steps
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        kPAR, alphaSWR, gammaSWR,
                                                                        ddd, 
                                                                        ntimesteps, 0.1);
  List sunshade = lightExtinctionAbsortion["sunshade"];
  List abs_PAR_SL_COH_list = sunshade["PAR_SL"];
  List abs_PAR_SH_COH_list = sunshade["PAR_SH"];
  List abs_SWR_SL_COH_list = sunshade["SWR_SL"];
  List abs_SWR_SH_COH_list = sunshade["SWR_SH"];
  List multilayer = lightExtinctionAbsortion["multilayer"];
  List abs_SWR_SL_ML_list = multilayer["SWR_SL"];
  List abs_SWR_SH_ML_list = multilayer["SWR_SH"];
  NumericVector fsunlit = lightExtinctionAbsortion["fsunlit"];
  NumericVector abs_SWR_can = lightExtinctionAbsortion["SWR_can"];
  NumericVector abs_SWR_soil = lightExtinctionAbsortion["SWR_soil"];

  NumericVector LAI_SL(numCohorts,0.0);
  NumericVector LAI_SH(numCohorts,0.0);
  NumericVector Vmax298SL(numCohorts,0.0);
  NumericVector Vmax298SH(numCohorts,0.0);
  NumericVector Jmax298SL(numCohorts,0.0);
  NumericVector Jmax298SH(numCohorts,0.0);

  for(int c=0;c<numCohorts;c++) {
    // Rcout<<"cohort "<<c<<":\n";
    //Constant properties through time steps
    NumericVector Vmax298layer(ncanlayers), Jmax298layer(ncanlayers);
    NumericVector SLarealayer(ncanlayers), SHarealayer(ncanlayers);
    double sn =0.0;
    for(int i=(ncanlayers-1);i>=0.0;i--) {
      //Effect of nitrogen concentration decay through the canopy
      double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sum(LAIme(_,c)));
      // Rcout<<" l"<<i<<" fsunlit: "<< fsunlit[i]<<" lai: "<< LAIme(i,c)<<" fn: "<< fn <<"\n";
      sn+=LAIme(i,c);
      SLarealayer[i] = LAIme(i,c)*fsunlit[i];
      SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
      Vmax298layer[i] = Vmax298[c]*fn;
      Jmax298layer[i] = Jmax298[c]*fn;
    }
    for(int i=0;i<ncanlayers;i++) {
      LAI_SL[c] +=SLarealayer[i];
      LAI_SH[c] +=SHarealayer[i];
      Vmax298SL[c] +=Vmax298layer[i]*LAIme(i,c)*fsunlit[i];
      Jmax298SL[c] +=Jmax298layer[i]*LAIme(i,c)*fsunlit[i];
      Vmax298SH[c] +=Vmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
      Jmax298SH[c] +=Jmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
    }
    
  }
  // double kb = lightExtinctionAbsortion["kb"];  //Proportion of sunlit extinction coefficient
  // double gbf = lightExtinctionAbsortion["gbf"]; //Ground fractions
  // double gdf = lightExtinctionAbsortion["gdf"];
  
  //Determine canopy vertical layer corresponding to cohort canopy, sunlit and shade leaves for each cohort
  IntegerVector iLayerCohort(numCohorts), iLayerSunlit(numCohorts), iLayerShade(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    double num = 0.0, den = 0.0, numsl=0.0, densl =0.0, numsh = 0.0, densh=0.0;
    for(int i=0;i<ncanlayers;i++) {
      num += LAIme(i,c)*zmid[i];
      den += LAIme(i,c);
      numsl += LAIme(i,c)*zmid[i]*fsunlit[i];
      densl += LAIme(i,c)*fsunlit[i];
      numsh += LAIme(i,c)*zmid[i]*(1.0 - fsunlit[i]);
      densh += LAIme(i,c)*(1.0-fsunlit[i]);
    }
    double hc_sl = numsl/densl;
    double hc_sh = numsh/densh;
    double hc  = num/den;
    for(int i=0;i<ncanlayers;i++) {
      if((hc > zlow[i]) & (hc <=zup[i])) iLayerCohort[c] = i;
      if((hc_sl > zlow[i]) & (hc_sl <=zup[i])) iLayerSunlit[c] = i;
      if((hc_sh > zlow[i]) & (hc_sh <=zup[i])) iLayerShade[c] = i;
    }
    // Rcout << c << " "<< hc_sl<<" "<< iLayerSunlit[c]<< " "<< hc_sh<<" "<< iLayerShade[c]<<"\n";
  }
  
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
        layerConnected(c,l)= true;
        if(layerConnected(c,l)) nlayerscon[c]=nlayerscon[c]+1;
      } else {
        layerConnected(c,l) = false;
      }
    }
    if((nlayerscon[c]==0) & verbose) Rcout<<"D";
  }
  
  //Average sap fluidity
  double sapFluidityDay = 1.0/waterDynamicViscosity((tmin+tmax)/2.0);
  
  //Hydraulics: supply functions
  List soil_c;
  if(plantWaterPools) {
    soil_c= clone(soil); //Clone soil
    //Calculate average rhizosphere moisture, including rhizosphere overlaps
    Wrhizo = cohortRhizosphereMoisture(Wpool, RHOP);
  }
  List supply(numCohorts);
  List supplyAboveground(numCohorts);
  supply.attr("names") = above.attr("row.names");
  for(int c=0;c<numCohorts;c++) {
    if(plantWaterPools) { 
      //Copy rhizosphere moisture to soil moisture
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wrhizo(c,l);
      //Update soil water potential from pool moisture
      psiSoil = psi(soil_c,soilFunctions); 
    }
    
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
        VCroot_kmaxc[cnt] = sapFluidityDay*VCroot_kmax(c,l);
        VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l);
        psic[cnt] = psiSoil[l];
        VG_nc[cnt] = VG_n[l];
        VG_alphac[cnt] = VG_alpha[l];
        cnt++;
      }
    }
    // double minFlow = std::max(0.0,1000.0*(Gwmin[c]*(tmin+tmax)/2.0)/Patm);
    // Rcout<<minFlow<<"\n";
    if(nlayerscon[c]>0) {
      //Build supply function networks 
      if(!capacitance) {
        supply[c] = supplyFunctionNetwork(psic,
                                          VGrhizo_kmaxc,VG_nc,VG_alphac,
                                          VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                                          sapFluidityDay*VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                          sapFluidityDay*VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                          NumericVector::create(StemPLCVEC[c],StemPLCVEC[c]), 
                                          0.0, maxNsteps, 
                                          ntrial, psiTol, ETol, 0.001); 
      } else {
        supply[c] = supplyFunctionNetworkStem1(psic,
                                               VGrhizo_kmaxc,VG_nc,VG_alphac,
                                               sapFluidityDay*VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                               sapFluidityDay*VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                               0.0, //StemPLCVEC[c],
                                               0.0, maxNsteps, 
                                               ntrial, psiTol, ETol, 0.001); 
        
      }
    } else {
      stop("Plant cohort not connected to any soil layer!");
    }
  }
  //Sugar conc in sapwood and leaf of each cohort
  NumericVector sugarLeaf(numCohorts, 0.0);
  NumericVector sugarSapwood(numCohorts, 0.0);
  for(int c=0;c<numCohorts;c++) {
    sugarLeaf[c] = sugarConcentration(LeafPI0[c],20.0, nonSugarConcentration);
    sugarSapwood[c] = sugarConcentration(StemPI0[c],20.0, nonSugarConcentration);
  }
  
  //Transpiration and photosynthesis
  NumericVector psiBk(nlayers);
  for(int l=0;l<nlayers;l++) psiBk[l] = psiSoil[l]; //Store initial soil water potential
  NumericMatrix K(numCohorts, nlayers);
  NumericVector Eplant(numCohorts, 0.0), Anplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericMatrix Rninst(numCohorts,ntimesteps);
  NumericMatrix dEdPInst(numCohorts, ntimesteps);
  NumericMatrix Qinst(numCohorts,ntimesteps);
  NumericMatrix Einst(numCohorts, ntimesteps);
  NumericMatrix Aninst(numCohorts, ntimesteps), Aginst(numCohorts, ntimesteps);
  NumericMatrix LeafPsiInst(numCohorts, ntimesteps), StemPsiInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympPsiInst(numCohorts, ntimesteps), StemSympPsiInst(numCohorts, ntimesteps);
  NumericMatrix LeafRWCInst(numCohorts, ntimesteps), StemRWCInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympRWCInst(numCohorts, ntimesteps), StemSympRWCInst(numCohorts, ntimesteps);
  NumericMatrix RootPsiInst(numCohorts, ntimesteps);
  NumericMatrix PWBinst(numCohorts, ntimesteps);
  NumericMatrix E_SL(numCohorts, ntimesteps);
  NumericMatrix E_SH(numCohorts, ntimesteps);
  NumericMatrix An_SL(numCohorts, ntimesteps), Ag_SL(numCohorts, ntimesteps);
  NumericMatrix An_SH(numCohorts, ntimesteps), Ag_SH(numCohorts, ntimesteps);
  NumericMatrix Psi_SL(numCohorts, ntimesteps);
  NumericMatrix Psi_SH(numCohorts, ntimesteps);
  NumericMatrix Ci_SL(numCohorts, ntimesteps);
  NumericMatrix Ci_SH(numCohorts, ntimesteps);
  NumericMatrix SWR_SL(numCohorts, ntimesteps);
  NumericMatrix SWR_SH(numCohorts, ntimesteps);
  NumericMatrix PAR_SL(numCohorts, ntimesteps);
  NumericMatrix PAR_SH(numCohorts, ntimesteps);
  NumericMatrix LWR_SL(numCohorts, ntimesteps);
  NumericMatrix LWR_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SL(numCohorts, ntimesteps);
  NumericMatrix VPD_SH(numCohorts, ntimesteps);
  NumericMatrix VPD_SL(numCohorts, ntimesteps);
  NumericMatrix Temp_SH(numCohorts, ntimesteps);
  NumericMatrix Temp_SL(numCohorts, ntimesteps);
  NumericVector minLeafPsi(numCohorts,0.0), maxLeafPsi(numCohorts,-99999.0); 
  NumericVector maxGSW_SL(numCohorts,-99999.0), maxGSW_SH(numCohorts,-99999.0); 
  NumericVector minGSW_SL(numCohorts,99999.0), minGSW_SH(numCohorts,99999.0); 
  NumericVector maxTemp_SL(numCohorts,-99999.0), maxTemp_SH(numCohorts,-99999.0); 
  NumericVector minTemp_SL(numCohorts,99999.0), minTemp_SH(numCohorts,99999.0); 
  NumericVector minLeafPsi_SL(numCohorts,0.0), maxLeafPsi_SL(numCohorts,-99999.0); 
  NumericVector minLeafPsi_SH(numCohorts,0.0), maxLeafPsi_SH(numCohorts,-99999.0);
  NumericVector minStemPsi(numCohorts, 0.0), minRootPsi(numCohorts,0.0); //Minimum potentials experienced
  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  if(numCohorts>0) std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), 0.0);
  NumericMatrix PLC(numCohorts, ntimesteps);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  NumericVector dEdPm(numCohorts);
  NumericVector PWB(numCohorts,0.0);
  
  List outPhotoSunlit(numCohorts);
  List outPhotoShade(numCohorts);
  List outPMSunlit(numCohorts);
  List outPMShade(numCohorts);
  outPhotoSunlit.attr("names") = above.attr("row.names");
  outPhotoShade.attr("names") = above.attr("row.names");
  outPMSunlit.attr("names") = above.attr("row.names");
  outPMShade.attr("names") = above.attr("row.names");
  
  List lwrExtinctionList(ntimesteps);
  
  for(int n=0;n<ntimesteps;n++) { //Time loop
    //Longwave radiation
    List lwrExtinction = longwaveRadiationSHAW(LAIme, LAImd, LAImx, 
                                               lwdr[n], Tsoil[0], Tair);
    lwrExtinctionList[n] = lwrExtinction;
    net_LWR_soil[n] = lwrExtinction["Lnet_ground"];
    net_LWR_can[n]= lwrExtinction["Lnet_canopy"];
    NumericMatrix Lnet_cohort_layer = lwrExtinction["Lnet_cohort_layer"];
    
    //Retrieve radiation absorbed for the current time step
    NumericVector absPAR_SL_COH = abs_PAR_SL_COH_list[n];
    NumericVector absPAR_SH_COH = abs_PAR_SH_COH_list[n];
    NumericVector absSWR_SL_COH = abs_SWR_SL_COH_list[n];
    NumericVector absSWR_SH_COH = abs_SWR_SH_COH_list[n];
    NumericMatrix absSWR_SL_ML = abs_SWR_SL_ML_list[n];
    NumericMatrix absSWR_SH_ML = abs_SWR_SH_ML_list[n];

    for(int c=0;c<numCohorts;c++) { //Plant cohort loop
      //Current osmotic potentials
      double leafpi0 = osmoticWaterPotential(sugarLeaf[c], Tair[iLayerCohort[c]], nonSugarConcentration);
      double stempi0 = osmoticWaterPotential(sugarSapwood[c], Tair[iLayerCohort[c]], nonSugarConcentration);
      
      
      //default values
      dEdPInst(c,n) = 0.0;
      Einst(c,n) = 0.0;
      Aginst(c,n) = 0.0;
      Aninst(c,n) = 0.0;
      
      if(LAIphe[c]>0.0) { //Process transpiration and photosynthesis only if there are some leaves
        PAR_SL(c,n) = absPAR_SL_COH[c];
        PAR_SH(c,n) = absPAR_SH_COH[c];
        SWR_SL(c,n) = absSWR_SL_COH[c];
        SWR_SH(c,n) = absSWR_SH_COH[c];
        // for(int j=0;j<ncanlayers;j++) Rcout<< n << " "<< c<< " "<<j<<" " << Lnet_cohort_layer(j,c)<<"\n";
        // Rcout<< n << " "<< c<< " LAIsl: " << LAI_SL[c]<< " LAIsh: " << LAI_SH[c]<< " LWRnet: "<< sum(Lnet_cohort_layer(_,c))<<" "<< sum(Lnet_cohort_layer(_,c)*fsunlit)<< " "<<sum(Lnet_cohort_layer(_,c)*(1.0 - fsunlit))<<"\n";
        LWR_SL(c,n) = sum(Lnet_cohort_layer(_,c)*fsunlit);
        LWR_SH(c,n) = sum(Lnet_cohort_layer(_,c)*(1.0 - fsunlit));
        
        //NumericVector PLCStemPrev = NumericVector::create(StemPLCVEC[c],StemPLCVEC[c]);
        //NumericVector psiStemPrev = NumericVector::create(Stem1PsiVEC[c],Stem2PsiVEC[c]);
        //double LeafPsiPrev = LeafPsiVEC[c];

        // Rcout<<c<<" E "<<EinstPrev<<" PR "<< RootPsiPrev<<" PL "<<LeafPsiPrev<< " PS "<<psiStemPrev[0]<< " "<<rwcsleafPrev<< " "<<RWCStemPrev[0]<<"\n";
        NumericVector fittedE, dEdP;
        NumericVector LeafPsi, psiRootCrown;
        
        //Retrieve supply functions
        List sFunctionBelow, sFunctionAbove;
        if(!capacitance) {
          sFunctionAbove = supply[c];
          sFunctionBelow = supply[c];
        } else {
          double psiPLCStem = apoplasticWaterPotential(1.0-StemPLCVEC[c], VCstem_c[c], VCstem_d[c]);
          double psiRootCrownFake = std::min(0.0,E2psiXylemUp(EinstVEC[c], Stem1PsiVEC[c],VCstem_kmax[c]*2.0, VCstem_c[c], VCstem_d[c], psiPLCStem));
          if(NumericVector::is_na(psiRootCrownFake)) psiRootCrownFake = 0.0;
          double psiFineRootFake= std::min(0.0,E2psiXylemUp(EinstVEC[c], psiRootCrownFake,VCroot_kmax_sum[c], VCroot_c[c], VCroot_d[c]));
          if(NumericVector::is_na(psiFineRootFake)) psiFineRootFake = 0.0;
          // Rcout<< c << " EinstVEC[c] "<< EinstVEC[c] << " Stem1PsiVEC[c] "<< Stem1PsiVEC[c]<<" psiFineRootFake "<< psiFineRootFake << " psiRootCrownFake "<< psiRootCrownFake<<"\n";
          // sFunctionAbove = supplyAboveground[c];
          double sapFluidityBelow = 1.0/waterDynamicViscosity(Tsoil[0]);
          double sapFluidityAbove = 1.0/waterDynamicViscosity(Tair[iLayerCohort[c]]);
          sFunctionAbove = supplyFunctionFineRootLeaf(psiFineRootFake,
                                                      sapFluidityBelow*VCroot_kmax_sum[c], VCroot_c[c], VCroot_d[c],
                                                      sapFluidityAbove*VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                                      sapFluidityAbove*VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                                      StemPLCVEC[c], 
                                                      0.0, maxNsteps, 
                                                      ETol, 0.001);
          
          sFunctionBelow = supply[c];
        }

        //Determine turgor loss point (as proxy of stomatal closure)
        double psiTlp = turgorLossPoint(leafpi0, LeafEPS[c]);
        
        //Retrieve transpiration, LeafPsi and dEdP vectors
        fittedE = sFunctionAbove["E"];
        dEdP = sFunctionAbove["dEdP"];
        LeafPsi = sFunctionAbove["psiLeaf"];
        
        //Get info from sFunctionAbove
        psiRootCrown = sFunctionAbove["psiRootCrown"];
        
        if(fittedE.size()>0) {
          //Photosynthesis function for sunlit and shade leaves
          DataFrame photoSunlit = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerSunlit[c]], Patm,
                                                             Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], 
                                                             zWind[iLayerSunlit[c]], 
                                                             SWR_SL(c,n), LWR_SL(c,n), 
                                                             irradianceToPhotonFlux(PAR_SL(c,n)), 
                                                             NSPLVEC[c]*Vmax298SL[c], 
                                                             NSPLVEC[c]*Jmax298SL[c], 
                                                             leafWidth[c], LAI_SL[c]);
          DataFrame photoShade = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerShade[c]], Patm,
                                                            Tair[iLayerShade[c]], VPair[iLayerShade[c]], 
                                                            zWind[iLayerShade[c]], 
                                                            SWR_SH(c,n), LWR_SH(c,n), 
                                                            irradianceToPhotonFlux(PAR_SH(c,n)),
                                                            NSPLVEC[c]*Vmax298SH[c], 
                                                            NSPLVEC[c]*Jmax298SH[c], 
                                                            leafWidth[c], LAI_SH[c]);
          
          NumericVector AgSunlit = photoSunlit["GrossPhotosynthesis"];
          NumericVector AgShade = photoShade["GrossPhotosynthesis"];
          NumericVector AnSunlit = photoSunlit["NetPhotosynthesis"];
          NumericVector AnShade = photoShade["NetPhotosynthesis"];
          NumericVector GswSunlit = photoSunlit["Gsw"];
          NumericVector GswShade = photoShade["Gsw"];
          NumericVector CiSunlit = photoSunlit["Ci"];
          NumericVector CiShade = photoShade["Ci"];
          NumericVector VPDSunlit = photoSunlit["LeafVPD"];
          NumericVector VPDShade = photoShade["LeafVPD"];
          NumericVector TempSunlit = photoSunlit["LeafTemperature"];
          NumericVector TempShade = photoShade["LeafTemperature"];
          
          
          //Profit maximization
          List PMSunlit, PMShade;
          int iPMSunlit = 0, iPMShade = 0;
          
          if(!cochard) { //Pure Sperry model
            PMSunlit = profitMaximization(sFunctionAbove, photoSunlit,  Gswmin[c], Gswmax[c], gainModifier, costModifier, costWater);
            PMShade = profitMaximization(sFunctionAbove, photoShade,  Gswmin[c],Gswmax[c], gainModifier, costModifier, costWater);
            iPMSunlit = PMSunlit["iMaxProfit"];
            iPMShade = PMShade["iMaxProfit"];
          } else {
            if(LeafPsi[c] < psiTlp) {  //Is leaf turgor zero
              iPMSunlit = 0;
              iPMShade  = 0;
              for(int j=0;j<(GswSunlit.size()-1);j++) if(GswSunlit[j]<Gswmin[c]) iPMSunlit++;
              for(int j=0;j<(GswShade.size()-1);j++) if(GswShade[j]<Gswmin[c]) iPMShade++;
            } else {
              PMSunlit = profitMaximization(sFunctionAbove, photoSunlit,  Gswmin[c], Gswmax[c], gainModifier, costModifier, costWater);
              PMShade = profitMaximization(sFunctionAbove, photoShade,  Gswmin[c],Gswmax[c], gainModifier, costModifier, costWater);
              iPMSunlit = PMSunlit["iMaxProfit"];
              iPMShade = PMShade["iMaxProfit"];
            }
          }
          
          //Store?
          if(!IntegerVector::is_na(stepFunctions)) {
            if(n==stepFunctions) {
              outPhotoSunlit[c] = photoSunlit;
              outPhotoShade[c] = photoShade;
              outPMSunlit[c] = PMSunlit;
              outPMShade[c] = PMShade;
            }
          }
          // Rcout<<iPMSunlit<<" "<<iPMShade <<" "<<GwSunlit[iPMSunlit]<<" "<<GwShade[iPMShade]<<" "<<fittedE[iPMSunlit]<<" "<<fittedE[iPMShade]<<"\n";
          //Get leaf status
          E_SH(c,n) = fittedE[iPMShade];
          E_SL(c,n) = fittedE[iPMSunlit];
          Psi_SH(c,n) = LeafPsi[iPMShade];
          Psi_SL(c,n) = LeafPsi[iPMSunlit];
          An_SH(c,n) = AnShade[iPMShade];
          An_SL(c,n) = AnSunlit[iPMSunlit];
          Ag_SH(c,n) = AgShade[iPMShade];
          Ag_SL(c,n) = AgSunlit[iPMSunlit];
          Ci_SH(c,n) = CiShade[iPMShade];
          Ci_SL(c,n) = CiSunlit[iPMSunlit];
          GSW_SH(c,n)= GswShade[iPMShade];
          GSW_SL(c,n)= GswSunlit[iPMSunlit];
          VPD_SH(c,n)= VPDShade[iPMShade];
          VPD_SL(c,n)= VPDSunlit[iPMSunlit];
          Temp_SH(c,n)= TempShade[iPMShade];
          Temp_SL(c,n)= TempSunlit[iPMSunlit];
          
          //Scale photosynthesis
          double Agsum = AgSunlit[iPMSunlit]*LAI_SL[c] + AgShade[iPMShade]*LAI_SH[c];
          double Ansum = AnSunlit[iPMSunlit]*LAI_SL[c] + AnShade[iPMShade]*LAI_SH[c];
          Aginst(c,n) = (1e-6)*12.01017*Agsum*tstep;
          Aninst(c,n) = (1e-6)*12.01017*Ansum*tstep;
          
          //Average flow from sunlit and shade leaves
          double Eaverage = (fittedE[iPMSunlit]*LAI_SL[c] + fittedE[iPMShade]*LAI_SH[c])/(LAI_SL[c] + LAI_SH[c]);
          
          
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
          if(iPM==-1) {
            Rcout<<"\n iPM -1! Eaverage="<< Eaverage << " fittedE.size= "<< fittedE.size()<<" iPMSunlit="<< iPMSunlit<< " fittedE[iPMSunlit]="<<fittedE[iPMSunlit]<<" iPMShade="<<iPMShade<<" fittedE[iPMShade]="<<fittedE[iPMShade]<<"\n";
            stop("");
          }

          //Store instantaneous total conductance
          dEdPInst(c,n) = dEdP[iPM];
          
          //Store instantaneous flow and leaf water potential
          EinstVEC[c] = fittedE[iPM];
          LeafPsiVEC[c] = LeafPsi[iPM];
          RootCrownPsiVEC[c] = psiRootCrown[iPM]; 
          
          //Scale from instantaneous flow to water volume in the time step
          Einst(c,n) = fittedE[iPM]*0.001*0.01802*LAIphe[c]*tstep; 
          
          NumericVector Esoilcn(nlayerscon[c],0.0);
          NumericVector ElayersVEC(nlayerscon[c],0.0);
          
          
          //Get info from sFunctionBelow (this will be different depending on wether capacitance is considered)
          NumericMatrix ERhizo = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["ERhizo"]);
          NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["psiRhizo"]);
          
          if(!capacitance) {
            //Store steady state stem and rootcrown and root surface water potential values
            NumericMatrix newStemPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionAbove["psiStem"]);
            Stem1PsiVEC[c] = newStemPsi(iPM,0); 
            Stem2PsiVEC[c] = newStemPsi(iPM,1);
            for(int lc=0;lc<nlayerscon[c];lc++) {
              ElayersVEC[lc] = ERhizo(iPM,lc)*tstep; //Scale according to the time step
            }
            //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
            int cl = 0;
            for(int l=0;l<nlayers;l++) {
              if(layerConnected(c,l)) {
                RhizoPsiMAT(c,l) = RhizoPsi(iPM,cl);
                cl++;
              } 
            }
            StemSympPsiVEC[c] = Stem1PsiVEC[c]; //Stem symplastic compartment coupled with apoplastic compartment
            LeafSympPsiVEC[c] = LeafPsiVEC[c]; //Leaf symplastic compartment coupled with apoplastic compartment
            
            // Store the PLC corresponding to stem1 water potential
            if(cavitationRefill!="total") {
              StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
            } else { //Immediate refilling
              StemPLCVEC[c] = 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
            }
            
          } else {
            //Store steady state stem2 water potential
            NumericVector newStemPsi2 = Rcpp::as<Rcpp::NumericVector>(sFunctionAbove["psiStem2"]);
            Stem2PsiVEC[c] = newStemPsi2[iPM];
            
            NumericVector newStemPsi1 = Rcpp::as<Rcpp::NumericVector>(sFunctionBelow["psiStem1"]);

            int iPMB = -1;
            
            //TO DO: stem segment water balance
            //Estimate current apoplastic and symplastic volumes
            // NOTE: Vsapwood and Vleaf are in l·m-2
            double VLeafSymp_mmolmax = 1000.0*((Vleaf[c]*(1.0-LeafAF[c]))/0.018); //mmol·m-2
            double VStemSymp_mmolmax = 1000.0*((Vsapwood[c]*(1.0-StemAF[c]))/0.018); //mmol·m-2
            //Substract from maximum apoplastic compartment embolized conduits
            double VStemApo_mmolmax = 1000.0*((Vsapwood[c]*StemAF[c])/0.018); //mmol·m-2
            double RWCLeafSymp = symplasticRelativeWaterContent(LeafSympPsiVEC[c], leafpi0, LeafEPS[c]); //mmol·m-2
            double RWCStemSymp = symplasticRelativeWaterContent(StemSympPsiVEC[c], stempi0, StemEPS[c]); //mmol·m-2
            double VLeafSymp_mmol = VLeafSymp_mmolmax * RWCLeafSymp;
            double VStemSymp_mmol = VStemSymp_mmolmax * RWCStemSymp;
            double Vcav = 0.0;
            //Perform water balance
            // Rcout<<"\n"<<c<<" Before - iPM " << iPM<< " EinstVEC[c]: "<< EinstVEC[c]<<" Vol: "<<VStemApo_mmol<<" RWC:"<< RWCStemApo <<" Psi: "<< Stem1PsiVEC[c]<< " LeafPsiVEC[c]: "<<LeafPsiVEC[c]<<"\n";
            for(double scnt=0.0; scnt<tstep;scnt += 1.0) {
              //Find flow corresponding to Stem1PsiVEC[c]
              //Find iPM for water potential corresponding to the current water potential
              double absDiff = 99999999.9;
              iPMB = -1;
              for(int k=0;k<newStemPsi1.size();k++){ //Only check up to the size of fittedE
                double adk = std::abs(newStemPsi1[k]-Stem1PsiVEC[c]);
                if(adk<absDiff) {
                  absDiff = adk;
                  iPMB = k;
                }
              }
              if(iPMB==-1) {
                Rcout<<"\n Stem1PsiVEC[c]="<< Stem1PsiVEC[c] << " newStemPsi1.size= "<< newStemPsi1.size()<<"\n";
                stop("iPMB = -1");
              }
              // Stem1PsiVEC[c] = newStemPsi1[iPMB];
                
              //Add flow from soil to ElayersVEC
              for(int lc=0;lc<nlayerscon[c];lc++) ElayersVEC[lc] += ERhizo(iPMB,lc); 
              
              //Calculate stem and leaf lateral flows
              double Flatstem = (StemSympPsiVEC[c] - Stem1PsiVEC[c])*klatstem;
              double Flatleaf = (LeafSympPsiVEC[c] - LeafPsiVEC[c])*klatleaf;


              //Leaf symplastic water balance
              VLeafSymp_mmol += (-Flatleaf);
              RWCLeafSymp = VLeafSymp_mmol/VLeafSymp_mmolmax;
              LeafSympPsiVEC[c] = symplasticWaterPotential(std::min(1.0,RWCLeafSymp), leafpi0, LeafEPS[c]);
              if(NumericVector::is_na(LeafSympPsiVEC[c]))  LeafSympPsiVEC[c] = -40.0;
              
              //Stem symplastic water balance
              VStemSymp_mmol += (-Flatstem);
              RWCStemSymp = VStemSymp_mmol/VStemSymp_mmolmax;
              StemSympPsiVEC[c] = symplasticWaterPotential(std::min(1.0,RWCStemSymp), stempi0, StemEPS[c]);
              if(NumericVector::is_na(StemSympPsiVEC[c]))  StemSympPsiVEC[c] = -40.0;
              
              //Stem apoplastic water balance
              double Vchange = (Flatstem + sum(ERhizo(iPMB,_)) - (EinstVEC[c] - Flatleaf)) + Vcav;
              
              Stem1PsiVEC[c] = Stem1PsiVEC[c] + eps_xylem*(Vchange/VStemApo_mmolmax);
              
              // VStemApo_mmol += (Flatstem + sum(ERhizo(iPMB,_)) - (EinstVEC[c] - Flatleaf));
              // RWCStemApo = VStemApo_mmol/VStemApo_mmolmax;
              // Stem1PsiVEC[c] = apoplasticWaterPotential(std::min(1.0,RWCStemApo), VCstem_c[c], VCstem_d[c]);
              // if(NumericVector::is_na(Stem1PsiVEC[c]))  Stem1PsiVEC[c] = -40.0;
              

              //Recalculate PLC and calculate volume corresponding to new cavitation
              double plc_old = StemPLCVEC[c];
              if(cavitationRefill!="total") {
                StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
                Vcav = VStemApo_mmolmax*(StemPLCVEC[c]-plc_old);
              } else { //Immediate refilling
                StemPLCVEC[c] = 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
                Vcav = 0.0;
              }

              // if((c==2) & (n==1)) {
              //   Rcout<< iPMB<<"  sum(ERhizo(iPMB,_) "<<  sum(ERhizo(iPMB,_)) << " Flatstem: "<<Flatstem<<" Flatleaf: "<<Flatleaf<<" VStemApo_mmol: "<<VStemApo_mmol<<" Psi: "<< Stem1PsiVEC[c]<<" StemSympPsiVEC: "<< StemSympPsiVEC[c]<<" LeafSympPsiVEC: "<< LeafSympPsiVEC[c]<<"\n";
              //   if(scnt>10.0) stop("");
              // }
            }
            
            // Rcout<<c<<" after - EinstVEC: "<<EinstVEC[c] << " RWCStemApo: " << RWCStemApo << "  Stem1PsiVEC:"<< Stem1PsiVEC[c]<<" StemSympPsiVEC: "<< StemSympPsiVEC[c]<<" LeafSympPsiVEC: "<< LeafSympPsiVEC[c] <<"\n";

            //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
            int cl = 0;
            for(int l=0;l<nlayers;l++) {
              if(layerConnected(c,l)) {
                RhizoPsiMAT(c,l) = RhizoPsi(iPMB,cl);
                cl++;
              } 
            }
          }
          
          //Scale soil water extracted from leaf to cohort level
          for(int lc=0;lc<nlayerscon[c];lc++) {
            Esoilcn[lc] = ElayersVEC[lc]*0.001*0.01802*LAIphe[c]; //Scale from flow to water volume in the time step
          }
          
          //Balance between extraction and transpiration
          PWBinst(c,n) = sum(Esoilcn) - Einst(c,n);
          
          //Add step transpiration to daily plant cohort transpiration
          Eplant[c] += Einst(c,n);
          Anplant[c] += Aninst(c,n);
          Agplant[c] += Aginst(c,n);
          //Add PWB
          PWB[c] += PWBinst(c,n); 
          

          
          //Copy transpiration and from connected layers to transpiration from soil layers
          int cl = 0;
          for(int l=0;l<nlayers;l++) {
            if(layerConnected(c,l)) {
              SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
              soilLayerExtractInst(l,n) += Esoilcn[cl];
              cl++;
            } 
          }
          
        } else {
          if(verbose) Rcout<<"NS!";
          Psi_SH(c,n) = NA_REAL;
          Psi_SL(c,n) = NA_REAL;
          GSW_SH(c,n)= NA_REAL;
          GSW_SL(c,n)= NA_REAL;
          VPD_SH(c,n)= NA_REAL;
          VPD_SL(c,n)= NA_REAL;
          Temp_SH(c,n)= NA_REAL;
          Temp_SL(c,n)= NA_REAL;
        }        
      } else if(N[c]>0.0) { //Cohorts with living individuals but no LAI should be in equilibrium with soil (i.e. no transpiration)
        List sFunctionBelow = supply[c];
        NumericVector  psiRootCrown = sFunctionBelow["psiRootCrown"];
        RootCrownPsiVEC[c] = psiRootCrown[0];
        if(!capacitance) {
          NumericVector  psiStem1 = sFunctionBelow["psiStem"];
          Stem1PsiVEC[c] = psiStem1[0];
          StemSympPsiVEC[c] = psiStem1[0];
          NumericVector LeafPsi = sFunctionBelow["psiLeaf"];
          LeafPsiVEC[c] = LeafPsi[0];
          LeafSympPsiVEC[c] = LeafPsi[0];
        } else {
          NumericVector  psiStem1 = sFunctionBelow["psiStem1"];
          Stem1PsiVEC[c] = psiStem1[0];
          StemSympPsiVEC[c] = psiStem1[0];
          LeafPsiVEC[c] = Stem1PsiVEC[c];
          LeafSympPsiVEC[c] = StemSympPsiVEC[c];
        }
      }
      
      if(N[c]>0.0) {
        //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
        PLC(c,n) = StemPLCVEC[c];
        StemSympRWCInst(c,n) = symplasticRelativeWaterContent(StemSympPsiVEC[c], stempi0, StemEPS[c]);
        LeafSympRWCInst(c,n) = symplasticRelativeWaterContent(LeafSympPsiVEC[c], leafpi0, LeafEPS[c]);
        StemRWCInst(c,n) = StemSympRWCInst(c,n)*(1.0 - StemAF[c]) + apoplasticRelativeWaterContent(Stem1PsiVEC[c], VCstem_c[c], VCstem_d[c])*StemAF[c];
        LeafRWCInst(c,n) = LeafSympRWCInst(c,n)*(1.0 - LeafAF[c]) + apoplasticRelativeWaterContent(LeafPsiVEC[c], VCleaf_c[c], VCleaf_d[c])*LeafAF[c];
        StemPsiInst(c,n) = Stem1PsiVEC[c]; 
        LeafPsiInst(c,n) = LeafPsiVEC[c]; //Store instantaneous (average) leaf potential
        RootPsiInst(c,n) = RootCrownPsiVEC[c]; //Store instantaneous root crown potential
        LeafSympPsiInst(c,n) = LeafSympPsiVEC[c];
        StemSympPsiInst(c,n) = StemSympPsiVEC[c];
        
        //Store the minimum water potential of the day (i.e. mid-day)
        minGSW_SL[c] = std::min(minGSW_SL[c], GSW_SL(c,n));
        minGSW_SH[c] = std::min(minGSW_SH[c], GSW_SH(c,n));
        maxGSW_SL[c] = std::max(maxGSW_SL[c], GSW_SL(c,n));
        maxGSW_SH[c] = std::max(maxGSW_SH[c], GSW_SH(c,n));
        minTemp_SL[c] = std::min(minTemp_SL[c], Temp_SL(c,n));
        minTemp_SH[c] = std::min(minTemp_SH[c], Temp_SH(c,n));
        maxTemp_SL[c] = std::max(maxTemp_SL[c], Temp_SL(c,n));
        maxTemp_SH[c] = std::max(maxTemp_SH[c], Temp_SH(c,n));
        minLeafPsi_SL[c] = std::min(minLeafPsi_SL[c],Psi_SL(c,n));
        minLeafPsi_SH[c] = std::min(minLeafPsi_SH[c],Psi_SH(c,n));
        maxLeafPsi_SL[c] = std::max(maxLeafPsi_SL[c],Psi_SL(c,n));
        maxLeafPsi_SH[c] = std::max(maxLeafPsi_SH[c],Psi_SH(c,n));
        minLeafPsi[c] = std::min(minLeafPsi[c],LeafPsiInst(c,n));
        maxLeafPsi[c] = std::max(maxLeafPsi[c],LeafPsiInst(c,n));
        minStemPsi[c] = std::min(minStemPsi[c],StemPsiInst(c,n));
        minRootPsi[c] = std::min(minRootPsi[c],RootPsiInst(c,n));
        for(int l=0;l<nlayers;l++) {
          minPsiRhizo(c,l) = std::min(minPsiRhizo(c,l),RhizoPsiMAT(c,l));
        }
      }
    } //End of cohort loop
    
    //CANOPY AND SOIL ENERGY BALANCE
    
    //Soil latent heat (soil evaporation)
    //Latent heat (snow fusion) as J/m2/s
    double soilEvapStep = abs_SWR_soil[n]*(soilEvaporation/sum(abs_SWR_soil));
    double snowMeltStep = abs_SWR_soil[n]*(snowMelt/sum(abs_SWR_soil));
    if(sum(abs_SWR_soil)==0.0) { // avoid zero sums
      soilEvapStep = 0.0; 
      snowMeltStep = 0.0;
    }
    if(SWE>0.0) {
      abs_SWR_soil[n] = 0.0; //Set SWR absorbed by soil to zero (for energy balance) if snow pack is present 
    }
    double LEsoilevap = (1e6)*meteoland::utils_latentHeatVaporisation(Tsoil[0])*soilEvapStep/tstep;
    double LEsnow = (1e6)*(snowMeltStep*0.33355)/tstep; // 0.33355 = latent heat of fusion
    LEsoil_heat[n] = LEsoilevap + LEsnow;
    
    
    //Canopy evaporation (mm) in the current step
    double canEvapStep = canopyEvaporation*(abs_SWR_can[n]/sum(abs_SWR_can));

    //Canopy convective heat exchange
    double RAcan = aerodynamicResistance(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
    Hcan_heat[n] = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;
    
    if(!multiLayerBalance) {//Canopy balance assuming a single layer
      //Soil-canopy turbulent heat exchange
      double wind2m = windSpeedMassmanExtinction(200.0, wind, LAIcell, canopyHeight);
      double RAsoil = aerodynamicResistance(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
      Hcansoil[n] = (meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG*(Tcan[n]-Tsoil[0]))/RAsoil;
      //Latent heat (evaporation + transpiration)
      double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tcan[n])*(sum(Einst(_,n)) + canEvapStep)/tstep;
      LEcan_heat[n] = LEwat; 
      //Canopy temperature changes
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEcan_heat[n] - Hcan_heat[n] - Hcansoil[n];
      double canopyAirThermalCapacity = meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG;
      double canopyThermalCapacity =  canopyAirThermalCapacity + (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
      double Tcannext = Tcan[n]+ std::max(-3.0, std::min(3.0, tstep*Ebal[n]/canopyThermalCapacity)); //Avoids changes in temperature that are too fast
      if(n<(ntimesteps-1)) Tcan[n+1] = Tcannext;
      for(int i=0;i<ncanlayers;i++) Tair[i] = Tcannext;
      
      
      //Soil energy balance including exchange with canopy
      Ebalsoil[n] = abs_SWR_soil[n] + net_LWR_soil[n] + Hcansoil[n] - LEsoil_heat[n]; //Here we use all energy escaping to atmosphere
      
      
      //Soil temperature changes
      NumericVector soilTchange = temperatureChange(dVec, Tsoil, sand, clay, Ws, Theta_FC, Ebalsoil[n]);
      for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tstep);
      if(n<(ntimesteps-1)) Tsoil_mat(n+1,_)= Tsoil;
      
    } else { //Multilayer canopy balance
      double moistureAtm = 0.622*(vpatm/Patm)*meteoland::utils_airDensity(Tatm[n],Patm);
      double CO2Atm = 0.409*Catm*44.01; //mg/m3
        
      double tsubstep = tstep/((double) nsubsteps); 
      double maxTchange = 3.0/((double) nsubsteps);
      double maxMoistureChange = 0.001/((double)nsubsteps); //=0.16 kPa per step
      double maxCO2Change = 180.0/((double)nsubsteps); //= 10 ppm per step
      double deltaZ = (verticalLayerSize/100.0); //Vertical layer size in m
      DataFrame LWR_layer = Rcpp::as<Rcpp::DataFrame>(lwrExtinction["LWR_layer"]);
      NumericVector LWRnet_layer = LWR_layer["Lnet"];
      Ebal[n] = 0.0;
      LEcan_heat[n] = 0.0;
      NumericVector Tairnext(ncanlayers), LElayer(ncanlayers), absSWRlayer(ncanlayers), Rnlayer(ncanlayers), Hleaflayer(ncanlayers);
      NumericVector layerThermalCapacity(ncanlayers);
      NumericVector moistureET(ncanlayers), rho(ncanlayers), moistureLayer(ncanlayers), moistureLayernext(ncanlayers);
      NumericVector CO2An(ncanlayers), CO2Layer(ncanlayers), CO2Layernext(ncanlayers);
      for(int i=0;i<ncanlayers;i++) {
        rho[i] = meteoland::utils_airDensity(Tair[i],Patm);
        absSWRlayer[i] = sum(absSWR_SL_ML(i,_)) + sum(absSWR_SH_ML(i,_));
        //Radiation balance
        Rnlayer[i] = absSWRlayer[i] + LWRnet_layer[i];
        NumericVector pLayer = LAIme(i,_)/LAIphe; //Proportion of each cohort LAI in layer i
        //Instantaneous layer transpiration
        //from mmolH2O/m2/s to kgH2O/m2/s
        double ElayerInst = 0.001*0.01802*sum(LAIme(i,_)*(E_SL(_,n)*fsunlit[i] + E_SH(_,n)*(1.0-fsunlit[i])));
        //Assumes Layers contribute to evaporation proportionally to their LAI fraction
        double layerEvapInst = (canEvapStep/tstep)*(LAIpe[i]/LAIcellexpanded);
        //Estimate instantaneous mgCO2/m2 absorption for the layer, taking into account the proportion of sunlit and shade leaves of each cohort
        //from micro.molCO2/m2/s to mgCO2/m2/s
        double Anlayer =(1e-3)*44.01*sum(LAIme(i,_)*(An_SL(_,n)*fsunlit[i] + An_SH(_,n)*(1.0-fsunlit[i])));
        // 1000.0*(44.01/12.0)*sum(Aninst(_,n)*pLayer); 
        double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tair[i])*(ElayerInst + layerEvapInst);
        LElayer[i] = LEwat; //Energy spent in vaporisation
        LEcan_heat[n] = LElayer[i];
        // layerThermalCapacity[i] = (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI/((double) ncanlayers);
        layerThermalCapacity[i] =  (0.5*(0.8*LAIpx[i] + 1.2*LAIpe[i]) + LAIpd[i])*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
        
        moistureLayer[i] = 0.622*(VPair[i]/Patm)*rho[i]; //kg water vapour/m3
        moistureET[i] = (ElayerInst + layerEvapInst)/(deltaZ); //kg water vapour /m3/s
        
        CO2Layer[i] = 0.409*Cair[i]*44.01; //mg/m3
        CO2An[i] = -1.0*Anlayer/(deltaZ); //mg/m3/s
        // Rcout<<n<< " "<< i<< " - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<< " Tini: "<< Tair[i]<<"\n";
        // Rcout<<n<< " "<< i<< " - moistureET: "<<moistureET[i]<<" moistureLayer: "<<moistureLayer[i]<<" CO2An: "<<CO2An[i]<< " CO2Layer: "<< CO2Layer[i]<<"\n";
      }
      //Add soil moisture evaporation
      moistureET[0] += soilEvapStep/(deltaZ*tstep); //kg/m3/s
      
      
      for(int s=0;s<nsubsteps;s++) {
        double RAsoil = aerodynamicResistance(200.0, std::max(zWind[0],1.0)); //Aerodynamic resistance to convective heat transfer from soil
        double Hcansoils = Cp_JKG*meteoland::utils_airDensity(Tair[0],Patm)*(Tair[0]-Tsoil[0])/RAsoil;
        // double Hcan_heats = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tair[ncanlayers-1]-Tatm[n]))/RAcan;
        for(int i=0;i<ncanlayers;i++) {
          double deltaH = 0.0;
          double deltaMoisture = 0.0;
          double deltaCO2 = 0.0;
          // double Hlayers = 0.0;
          //Add turbulent heat flow (positive gradient when temperature is larger above)
          if(i==0) { //Lower layer
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i])*uw[i])/(deltaZ*dU[i]);
            deltaH -= Hcansoils;
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i])*uw[i])/(deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i])*uw[i])/(deltaZ*dU[i]);
          } else if((i > 0) & (i<(ncanlayers-1))) { //Intermediate layers
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          } else if(i==(ncanlayers-1)){ //Upper layer
            deltaH -= (Cp_JKG*rho[i]*(Tatm[n] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureAtm - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Atm - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          }
          Hleaflayer[i] = 0.0;
          for(int c=0;c<numCohorts;c++) {
            double gHa = 0.189*pow(std::max(zWind[i],0.1)/(leafWidth[c]*0.0072), 0.5);
            double Hsunlit = 2.0*Cp_Jmol*rho[i]*(Temp_SL(c, n)-Tair[i])*gHa;
            double Hshade = 2.0*Cp_Jmol*rho[i]*(Temp_SH(c, n)-Tair[i])*gHa;
            // Rcout<<c<<" " << Hsunlit<< " "<<Hshade<<" \n";
            Hleaflayer[i] +=(Hsunlit*fsunlit[i] + Hshade*(1.0-fsunlit[i]))*LAIme(i,c);
          }
          double EbalLayer = Rnlayer[i] - LElayer[i] + Hleaflayer[i] + deltaH; 
          //Instantaneous changes in temperature due to internal energy balance
          double deltaT = EbalLayer/(rho[i]*Cp_JKG + layerThermalCapacity[i]); 
          
          // if(s==0) Rcout<<n<< " "<< i<< " "<< s <<" - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<<" H: "<<Hlayers<< " Ebal: "<<EbalLayer<< " LTC: " << rholayer*Cp_JKG + layerThermalCapacity[i]<< " Tini: "<< Tair[i]<< " deltaT: "<<deltaT<<"\n";
          Tairnext[i] = Tair[i] +  std::max(-1.0*maxTchange, std::min(maxTchange, tsubstep*deltaT)); //Avoids changes in temperature that are too fast
          //Changes in water vapour
          moistureLayernext[i] = moistureLayer[i] + std::max(-1.0*maxMoistureChange, std::min(maxMoistureChange, tsubstep*(moistureET[i]+ deltaMoisture)));
          // if(i==0) Rcout<<n<< " "<< i<< " "<< s <<" - moisture: "<<moistureLayer[i]<<" delta: "<<deltaMoisture<<" next: "<< moistureLayernext[i]<<"\n";
          //Changes in CO2
          CO2Layernext[i] = CO2Layer[i] + std::max(-1.0*maxCO2Change, std::min(maxCO2Change, tsubstep*(CO2An[i] + deltaCO2)));
        }
        for(int i=0;i<ncanlayers;i++) {
          Tair[i] = Tairnext[i]; 
          moistureLayer[i] = moistureLayernext[i];
          CO2Layer[i] = CO2Layernext[i];
        }
        
        //Soil energy balance including exchange with canopy
        double Ebalsoils =  abs_SWR_soil[n] + net_LWR_soil[n]  - LEsoil_heat[n] + Hcansoils; //Here we use all energy escaping to atmosphere
        Ebalsoil[n] +=Ebalsoils;
        Hcansoil[n] +=Hcansoils;
        //Soil temperature changes
        NumericVector soilTchange = temperatureChange(dVec, Tsoil, sand, clay, Ws, Theta_FC, Ebalsoils);
        for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tsubstep);
      }
      Hcansoil[n] = Hcansoil[n]/((double) nsubsteps);
      Ebalsoil[n] = Ebalsoil[n]/((double) nsubsteps);
      for(int i=0;i<ncanlayers;i++) {
        VPair[i] = ((moistureLayer[i]/rho[i])*Patm)/0.622;
        Cair[i] = CO2Layer[i]/(0.409*44.01);
        // Rcout<< n << " "<<i << " - " << moistureLayer[i]<< " "<< VPair[i]<<"\n";
      }
      // stop("kk");
      // Rcout<< n << " "<<abs_SWR_can[n]<< " = "<< sum(absSWRlayer)<<" = "<< (sum(absSWR_SL_ML) + sum(absSWR_SH_ML))<<" "<<net_LWR_can[n]<< " = "<< sum(LWRnet_layer)<<"\n";
      //Canopy energy balance
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEcan_heat[n] - Hcan_heat[n] - Hcansoil[n];
      if(n<(ntimesteps-1)) {
        Tcan[n+1] = sum(Tair*LAIpx)/sum(LAIpx); 
        Tsoil_mat(n+1,_)= Tsoil;
      }
    }
    if(n<(ntimesteps-1)) for(int i=0;i<ncanlayers;i++) {
      Tcan_mat(n+1,i) = Tair[i];
      VPcan_mat(n+1,i) = VPair[i];
    }
  } //End of timestep loop

  //4z. Plant daily drought stress (from root collar mid-day water potential)
  NumericVector SoilExtractCoh(numCohorts,0.0);
  NumericVector DDS(numCohorts, 0.0);
  for(int c=0;c<numCohorts;c++) {
    SoilExtractCoh[c] =  sum(SoilWaterExtract(c,_));
    PLCm[c] = sum(PLC(c,_))/((double)PLC.ncol());
    RWCsm[c] = sum(StemRWCInst(c,_))/((double)StemRWCInst.ncol());
    RWClm[c] = sum(LeafRWCInst(c,_))/((double)LeafRWCInst.ncol());
    RWCssm[c] = sum(StemSympRWCInst(c,_))/((double)StemSympRWCInst.ncol());
    RWClsm[c] = sum(LeafSympRWCInst(c,_))/((double)LeafSympRWCInst.ncol());
    dEdPm[c] = sum(dEdPInst(c,_))/((double)dEdPInst.ncol());  
    double maxConductance = maximumSoilPlantConductance(VGrhizo_kmax(c,_), VCroot_kmax(c,_), VCstem_kmax[c], VCleaf_kmax[c]);
    DDS[c] = Phe[c]*(1.0 - (dEdPm[c]/(sapFluidityDay*maxConductance)));
    
    if(cavitationRefill=="rate") {
      double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
      double r = refillMaximumRate*std::max(0.0, (StemSympPsiVEC[c] + 1.5)/1.5);
      StemPLCVEC[c] = std::max(0.0, StemPLCVEC[c] - (r/SAmax));
    }
  }
  
  
  //B.3 - Substract  extracted water from soil moisture 
  if(modifyInput){
    for(int l=0;l<nlayers;l++) {
      Ws[l] = std::max(Ws[l] - (sum(soilLayerExtractInst(l,_))/Water_FC[l]),0.0);
    } 
    if(plantWaterPools) {
      rhizosphereMoistureExtraction(SoilWaterExtract, Water_FC,
                                    Wpool, RHOP,
                                    poolProportions);
      // for(int l=0;l<nlayers;l++) {
      //   double Ws2 = 0.0;
      //   for(int c=0;c<numCohorts;c++) Ws2 +=Wpool(c,l)*poolProportions[c];
      // 
      //   Rcout<<l<<": "<< Ws[l]<< " = " << Ws2<<"\n";
      // }
    } else { //copy soil to the pools of all cohorts
      for(int c=0;c<numCohorts;c++) {
        for(int l=0;l<nlayers;l++) {
          Wpool(c,l) = Ws[l];
        }
      }
    }
  }
  //6. Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  LAIcohort.attr("names") = above.attr("row.names");
  
  DataFrame Tinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                      _["Tatm"] = Tatm, _["Tcan"] = Tcan, _["Tsoil"] = Tsoil_mat);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  DataFrame CEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["SWRcan"] = abs_SWR_can, _["LWRcan"] = net_LWR_can,
                                        _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["Ebalcan"] = Ebal);
  DataFrame SEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["Hcansoil"] = Hcansoil, _["LEsoil"] = LEsoil_heat, _["SWRsoil"] = abs_SWR_soil, _["LWRsoil"] = net_LWR_soil,
                                        _["Ebalsoil"] = Ebalsoil);
  List EB = List::create(_["Temperature"]=Tinst, _["CanopyEnergyBalance"] = CEBinst, _["SoilEnergyBalance"] = SEBinst,
                         _["TemperatureLayers"] = NA_REAL, _["VaporPressureLayers"] = NA_REAL);
  if(multiLayerBalance) {
    EB["TemperatureLayers"] = Tcan_mat;
    EB["VaporPressureLayers"] = VPcan_mat;
  }
  E_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  E_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));

  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  dEdPInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  RootPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aginst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PWBinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  for(int c=0;c<numCohorts;c++) {
    Vmax298SL[c] = Vmax298SL[c]/LAI_SL[c];
    Jmax298SL[c] = Jmax298SL[c]/LAI_SL[c];
    Vmax298SH[c] = Vmax298SH[c]/LAI_SH[c];
    Jmax298SH[c] = Jmax298SH[c]/LAI_SH[c];
  }
  DataFrame Sunlit = DataFrame::create(
    _["LAI"] = LAI_SL, 
    _["Vmax298"] = Vmax298SL,
    _["Jmax298"] = Jmax298SL,
    _["LeafPsiMin"] = minLeafPsi_SL, 
    _["LeafPsiMax"] = maxLeafPsi_SL, 
    _["GSWMin"] = minGSW_SL,
    _["GSWMax"] = maxGSW_SL,
    _["TempMin"] = minTemp_SL,
    _["TempMax"] = maxTemp_SL  
  );
  DataFrame Shade = DataFrame::create(
    _["LAI"] = LAI_SH, 
    _["Vmax298"] = Vmax298SH,
    _["Jmax298"] = Jmax298SH,
    _["LeafPsiMin"] = minLeafPsi_SH, 
    _["LeafPsiMax"] = maxLeafPsi_SH, 
    _["GSWMin"] = minGSW_SH,
    _["GSWMax"] = maxGSW_SH,
    _["TempMin"] = minTemp_SH,
    _["TempMax"] = maxTemp_SH  
  );
  Sunlit.attr("row.names") = above.attr("row.names");
  Shade.attr("row.names") = above.attr("row.names");
  
  List ShadeInst = List::create(
    _["Abs_SWR"] = SWR_SH,
    _["Abs_PAR"]=PAR_SH,
    _["Net_LWR"] = LWR_SH,
    _["Ag"] = Ag_SH,
    _["An"] = An_SH,
    _["Ci"] = Ci_SH,
    _["E"] = E_SH,
    _["Gsw"] = GSW_SH,
    _["VPD"] = VPD_SH,
    _["Temp"] = Temp_SH,
    _["Psi"] = Psi_SH);
  List SunlitInst = List::create(
    _["Abs_SWR"]=SWR_SL,
    _["Abs_PAR"]=PAR_SL,
    _["Net_LWR"] = LWR_SL,
    _["Ag"] = Ag_SL,
    _["An"] = An_SL,
    _["Ci"] = Ci_SL,
    _["E"] = E_SL,
    _["Gsw"] = GSW_SL,
    _["VPD"] = VPD_SL,
    _["Temp"] = Temp_SL,
    _["Psi"] = Psi_SL);
  
  List PlantsInst = List::create(
    _["E"]=Einst, _["Ag"]=Aginst, _["An"]=Aninst,
    _["dEdP"] = dEdPInst,
    _["RootPsi"] = RootPsiInst, 
    _["StemPsi"] = StemPsiInst,
    _["LeafPsi"] = LeafPsiInst,
    _["StemSympPsi"] = StemSympPsiInst,
    _["LeafSympPsi"] = LeafSympPsiInst,
    _["StemPLC"] = PLC, 
    _["StemRWC"] = StemRWCInst,
    _["LeafRWC"] = LeafRWCInst,
    _["StemSympRWC"] = StemSympRWCInst,
    _["LeafSympRWC"] = LeafSympRWCInst,
    _["PWB"] = PWBinst);
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["Extraction"] = SoilExtractCoh,
                                       _["Transpiration"] = Eplant,
                                       _["GrossPhotosynthesis"] = Agplant,
                                       _["NetPhotosynthesis"] = Anplant,
                                       _["RootPsi"] = minRootPsi, 
                                       _["StemPsi"] = minStemPsi, 
                                       _["StemPLC"] = PLCm, //Average daily stem PLC
                                       _["LeafPsiMin"] = minLeafPsi, 
                                       _["LeafPsiMax"] = maxLeafPsi, 
                                       _["dEdP"] = dEdPm,//Average daily soilplant conductance
                                       _["DDS"] = DDS, //Daily drought stress is the ratio of average soil plant conductance over its maximum value
                                       _["StemRWC"] = RWCsm,
                                       _["LeafRWC"] = RWClm,
                                       _["StemSympRWC"] = RWCssm,
                                       _["LeafSympRWC"] = RWClsm,
                                       _["WaterBalance"] = PWB);
  Plants.attr("row.names") = above.attr("row.names");
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell, 
                                              _["LAIlive"] = LAIcelllive, 
                                              _["LAIexpanded"] = LAIcellexpanded, 
                                              _["LAIdead"] = LAIcelldead);
  
  List l;
  if(!IntegerVector::is_na(stepFunctions)){
    l = List::create(_["cohorts"] = clone(cohorts),
                     _["EnergyBalance"] = EB,
                     _["Stand"] = Stand,
                     _["Extraction"] = SoilWaterExtract,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["RhizoPsi"] = minPsiRhizo,
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["PlantsInst"] = PlantsInst,
                     _["SunlitLeavesInst"] = SunlitInst,
                     _["ShadeLeavesInst"] = ShadeInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["LWRExtinction"] = lwrExtinctionList,
                     _["CanopyTurbulence"] = canopyTurbulence,
                     _["SupplyFunctions"] = supply,
                     _["PhotoSunlitFunctions"] = outPhotoSunlit,
                     _["PhotoShadeFunctions"] = outPhotoShade,
                     _["PMSunlitFunctions"] = outPMSunlit,
                     _["PMShadeFunctions"] = outPMShade);
  } else {
    l = List::create(_["cohorts"] = clone(cohorts),
                     _["EnergyBalance"] = EB,
                     _["Extraction"] = SoilWaterExtract,
                     _["RhizoPsi"] = minPsiRhizo,
                     _["Stand"] = Stand,
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["PlantsInst"] = PlantsInst,
                     _["SunlitLeavesInst"] = SunlitInst,
                     _["ShadeLeavesInst"] = ShadeInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["LWRExtinction"] = lwrExtinctionList,
                     _["CanopyTurbulence"] = canopyTurbulence,
                     _["SupplyFunctions"] = supply);
    
  } 
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}

// [[Rcpp::export("transp_transpirationSperry")]]
List transpirationSperry(List x, DataFrame meteo, int day,
                        double latitude, double elevation, double slope, double aspect,
                        double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                        int stepFunctions = NA_INTEGER, 
                        bool modifyInput = true) {
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
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  CharacterVector dateStrings = meteo.attr("row.names");
  NumericVector WindSpeed(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  std::string c = as<std::string>(dateStrings[day-1]);
  double prec = Precipitation[day-1];
  double rad = Radiation[day-1];
  double tmax = MaxTemperature[day-1];
  double tmin = MinTemperature[day-1];
  double tmaxPrev = tmax;
  double tminPrev = tmin;
  double tminNext = tmin;
  if(day>1) {
    tmaxPrev = MaxTemperature[day-2];
    tminPrev = MinTemperature[day-2];
  }
  if(day<(MaxTemperature.length()-1)) tminNext = MinTemperature[day];
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);

  return(transpirationSperry(x, tmin, tmax, tminPrev, tmaxPrev, tminNext, rhmin, rhmax, rad, wind, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, prec,
                     canopyEvaporation, snowMelt, soilEvaporation,
                     false, stepFunctions, 
                     modifyInput));
} 


List transpirationGranier(List x, double tday, double pet, 
                          bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
  String cavitationRefill = control["cavitationRefill"];
  String soilFunctions = control["soilFunctions"];
  double verticalLayerSize = control["verticalLayerSize"];
  bool plantWaterPools = control["plantWaterPools"];

  //Soil water at field capacity
  List soil = x["soil"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  
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
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);

  //Water pools
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  NumericMatrix Wrhizo;
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  //Parameters  
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  NumericVector Psi_Critic = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Critic"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE"]);
  NumericVector pRootDisc = Rcpp::as<Rcpp::NumericVector>(paramsTransp["pRootDisc"]);
  NumericVector Tmax_LAI(numCohorts, 0.134);
  NumericVector Tmax_LAIsq(numCohorts, -0.006);
  if(paramsTransp.containsElementNamed("Tmax_LAI")) {
    Tmax_LAI = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAI"]);
    Tmax_LAIsq = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Tmax_LAIsq"]);
  }
  //Communication vectors
  //Comunication with outside
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector PlantPsi = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]));
  NumericVector StemPLC = clone(Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]));
  
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector Phe(numCohorts,0.0);
  double s = 0.0, LAIcell = 0.0, canopyHeight = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0,LAIcelldead = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(LAIlive[c]>0) Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    else Phe[c]=0.0;
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
    LAIcellexpanded +=LAIphe[c];
    LAIcelllive += LAIlive[c];
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  NumericVector z(nz+1,0.0);
  NumericVector zmid(nz);
  for(int i=1;i<=nz;i++) {
    z[i] = z[i-1] + verticalLayerSize;
    zmid[i-1] = (verticalLayerSize/2.0) + verticalLayerSize*((double) (i-1));
  }
  
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(z, LAIphe,  LAIdead, H, CR, kPAR);

  
  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  NumericVector Tmax = pet*(Tmax_LAIsq*pow(LAIcell,2.0)+ Tmax_LAI*LAIcell); //From Granier (1999)
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts,0.0);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  //Soil input
  NumericVector psiSoil = psi(soil,soilFunctions); //Update soil water potential
  List soil_c;
  if(plantWaterPools) {
    soil_c= clone(soil); //Clone soil
    //Calculate average rhizosphere moisture, including rhizosphere overlaps
    Wrhizo = cohortRhizosphereMoisture(Wpool, RHOP);
  }
  int nlayers = Wpool.ncol();
  NumericMatrix EplantCoh(numCohorts, nlayers);
  NumericMatrix RootPsi(numCohorts, nlayers);
  NumericVector Eplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0);
  NumericVector Kl, epc, Vl;
  double WeibullShape=3.0;
  for(int c=0;c<numCohorts;c++) {
    if(plantWaterPools) { 
      //Copy rhizosphere moisture to soil moisture
      NumericVector W_c = soil_c["W"];
      for(int l=0;l<nlayers;l++) W_c[l] = Wrhizo(c,l);
      //Update soil water potential from pool moisture
      psiSoil = psi(soil_c,soilFunctions); 
    }
    for(int l=0;l<nlayers;l++) {
      double Klc = Psi2K(psiSoil[l], Psi_Extract[c], WeibullShape);
      //Limit Kl due to previous cavitation
      if(cavitationRefill!="total") {
        Klc = std::min(Klc, 1.0-StemPLC[c]); 
      }
      double epc = std::max(TmaxCoh[c]*Klc*V(c,l),0.0);
      RootPsi(c,l) = psiSoil[l]; //Set initial guess of root potential to soil values
      //If relative conductance is smaller than the value for root disconnection
      //Set root potential to minimum value before disconnection and transpiration from that layer to zero
      if(Klc<pRootDisc[c]) { 
        RootPsi(c,l) = K2Psi(pRootDisc[c],Psi_Extract[c],WeibullShape);
        Klc = pRootDisc[c]; //So that layer stress does not go below pRootDisc
        epc = 0.0; //Set transpiration from layer to zero
      }
      EplantCoh(c,l) = epc;
      Eplant[c] = Eplant[c] + epc;
      DDS[c] = DDS[c] + Phe[c]*(V(c,l)*(1.0 - Klc)); //Add stress from the current layer
    }
  }

  for(int c=0;c<numCohorts;c++) {
    PlantPsi[c] = averagePsi(RootPsi(c,_), V(c,_), WeibullShape, Psi_Extract[c]);
    if(cavitationRefill!="total") {
      StemPLC[c] = std::max(1.0 - Psi2K(PlantPsi[c],Psi_Critic[c],WeibullShape), StemPLC[c]); //Track current embolism if no refill
    } else {
      StemPLC[c] = 1.0 - Psi2K(PlantPsi[c],Psi_Critic[c],WeibullShape);
    }
    Agplant[c] = WUE[c]*Eplant[c];
  }
  
  
  if(modifyInput) {
    internalWater["StemPLC"] = StemPLC;
    internalWater["PlantPsi"] = PlantPsi;
  }
  //Modifies input soil
  if(modifyInput) {
    NumericVector Ws = soil["W"];
    for(int l=0;l<nlayers;l++) Ws[l] = Ws[l] - (sum(EplantCoh(_,l))/Water_FC[l]); 
    if(plantWaterPools) {
      rhizosphereMoistureExtraction(EplantCoh, Water_FC,
                                    Wpool, RHOP,
                                    poolProportions);
    } else { //copy soil to the pools of all cohorts
      for(int c=0;c<numCohorts;c++) {
        for(int l=0;l<nlayers;l++) {
          Wpool(c,l) = Ws[l];
        }
      }
    }
  }
  
  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell,
                                              _["LAIlive"] = LAIcelllive, 
                                              _["LAIexpanded"] = LAIcellexpanded, 
                                              _["LAIdead"] = LAIcelldead);
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["AbsorbedSWRFraction"] = CohASWRF, 
                                       _["Transpiration"] = Eplant, 
                                       _["GrossPhotosynthesis"] = Agplant,
                                       _["PlantPsi"] = PlantPsi, 
                                       _["DDS"] = DDS,
                                       _["StemPLC"] = StemPLC);
  Plants.attr("row.names") = above.attr("row.names");
  EplantCoh.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = EplantCoh);
  return(l);
}


// [[Rcpp::export("transp_transpirationGranier")]]
List transpirationGranier(List x, DataFrame meteo, int day,
                          bool modifyInput = true) {
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
  NumericVector PET = meteo["PET"];
  double pet = PET[day-1];
  double tday = MeanTemperature[day-1];
  return(transpirationGranier(x, tday, pet, modifyInput));
} 

