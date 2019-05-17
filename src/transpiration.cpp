#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "photosynthesis.h"
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1

// [[Rcpp::export("transp_profitMaximization")]]
List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, int type, double Gwmin, double Gwmax, double kleafmax = NA_REAL) {
  NumericVector supplyKterm = supplyFunction["kterm"];
  NumericVector supplyE = supplyFunction["E"];
  NumericVector supplydEdp = supplyFunction["dEdP"];
  NumericVector Ag = photosynthesisFunction["Photosynthesis"];
  NumericVector Gw = photosynthesisFunction["WaterVaporConductance"];
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
    mindEdp = std::min(mindEdp, supplydEdp[i]);
    maxdEdp = std::max(maxdEdp, supplydEdp[i]);
    minKterm = std::min(minKterm, supplyKterm[i]);
    maxKterm = std::max(maxKterm, supplyKterm[i]);
    Agmax = std::max(Agmax, Ag[i]);
  }
  //Evaluate profit for valid steps
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  for(int i=ini;i<fin;i++) {
    if(type==1) {
      // if(i<imaxdEdp) {
      //   cost[i]= 0.0; 
      // } else {
      //   cost[i] = std::max(0.0,(maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp));  
      // }
      cost[i] = (maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp);  
    }
    else {
      cost[i] = (maxKterm-supplyKterm[i])/(maxKterm - minKterm);
    }
    gain[i] = Ag[i]/Agmax;
    profit[i] = gain[i]-cost[i];
  }
  
  while((Gw[ini]<=Gwmin) & (ini<fin)) ini +=1; 
  while((Gw[fin]>=Gwmax) & (fin>ini)) fin -=1; 
  
  int imaxprofit=ini;
  double maxprofit=profit[ini];
  for(int i=ini+1;i<fin;i++){
    if((profit[i]>maxprofit)) {
      maxprofit = profit[i];
      imaxprofit = i;
    }
  }
  // Rcout<<ini<< " "<< fin<<" Gwmx= "<<Gwmax<<" Gwmin "<<Gwmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}


List transpirationSperry(List x, List soil, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
                  double latitude, double elevation, double slope, double aspect, 
                  double solarConstant, double delta, double prec,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
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
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIlive.size();
  
  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector psiVec = psi(soil, soilFunctions); //Get soil water potential
  int nlayers = psiVec.length();
  
  //Canopy params
  List canopyParams = Rcpp::as<Rcpp::List>(x["canopy"]);
  
  //Root distribution input
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(below["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VGrhizo_kmax"]);
  
  
  //Base parameters
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector albedo = Rcpp::as<Rcpp::NumericVector>(paramsBase["albedo"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  
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
  for(int c=0;c<numCohorts;c++) { //Reset photosynthesis and transpiration
    photosynthesis[c] = 0.0;
    transpiration[c] = 0.0;
  }
  NumericMatrix psiStemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["psiStem"]);
  NumericMatrix PLCstemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
  NumericMatrix RWCsstemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["RWCsympstem"]);
  NumericVector RWCsleafVEC = Rcpp::as<Rcpp::NumericVector>(x["RWCsympleaf"]);
  NumericVector psiLeafVEC = Rcpp::as<Rcpp::NumericVector>(x["psiLeaf"]);
  NumericVector psiRootVEC = Rcpp::as<Rcpp::NumericVector>(x["psiRoot"]);
  NumericMatrix psiRhizoMAT = Rcpp::as<Rcpp::NumericMatrix>(x["psiRhizo"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(x["Einst"]);
  
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double latrad = latitude * (PI/180.0);
  double asprad = aspect * (PI/180.0);
  double slorad = slope * (PI/180.0);
  
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
  double LAIcell = 0.0, LAIcelldead = 0.0, canopyHeight = 0.0, LAIcellmax = 0.0;
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcellmax += LAIlive[c];
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  NumericVector z(nz+1,0.0);
  NumericVector zmid(nz);
  for(int i=1;i<=nz;i++) {
    z[i] = z[i-1] + verticalLayerSize;
    zmid[i-1] = (verticalLayerSize/2.0) + verticalLayerSize*((double) (i-1));
  }
  NumericMatrix LAIme = LAIdistributionVectors(z, LAIphe, H, CR); //Expanded leaves
  NumericMatrix LAImd = LAIdistributionVectors(z, LAIdead, H, CR); //Dead (standing) leaves
  NumericMatrix LAImx = LAIdistributionVectors(z, LAIlive, H, CR); //Maximum leaf expansion
  
  //3. Wind extinction profile
  if(NumericVector::is_na(wind)) wind = defaultWindSpeed; //set to default if missing
  NumericVector zWind;
  zWind = windExtinctionCohort(H,CR, wind,LAIcell, canopyHeight);
  zWind.attr("names") = above.attr("row.names");
  double RAcan = aerodynamicResistance(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
  double wind2m = windSpeedMassmanExtinction(200.0, wind, LAIcell, canopyHeight);
  double RAsoil = aerodynamicResistance(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
  
  //4a. Instantaneous direct and diffuse shorwave radiation
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, slorad, asprad, delta,
                                                        rad, clearday, ntimesteps);
  NumericVector solarHour = ddd["SolarHour"]; //in radians
  
  //4b. Instantaneous air temperature (above canopy) and longwave radiation
  NumericVector Tatm(ntimesteps), lwdr(ntimesteps), Tcan(ntimesteps, NA_REAL), Tsunrise(ntimesteps);
  NumericVector LEcan_heat(ntimesteps), Hcan_heat(ntimesteps), LWRsoilcan(ntimesteps), LWRcanout(ntimesteps), Ebal(ntimesteps);
  NumericVector LWRsoilout(ntimesteps), Ebalsoil(ntimesteps), Hcansoil(ntimesteps), LEsoil_heat(ntimesteps);
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
  
  NumericVector LAI_SL(numCohorts,0.0);
  NumericVector LAI_SH(numCohorts,0.0);
  NumericVector Vmax298SL(numCohorts,0.0);
  NumericVector Vmax298SH(numCohorts,0.0);
  NumericVector Jmax298SL(numCohorts,0.0);
  NumericVector Jmax298SH(numCohorts,0.0);

  for(int c=0;c<numCohorts;c++) {
    // Rcout<<"cohort "<<c<<":\n";
    //Constant properties through time steps
    NumericVector Vmax298layer(nz), Jmax298layer(nz);
    NumericVector SLarealayer(nz), SHarealayer(nz);
    double sn =0.0;
    for(int i=(nz-1);i>=0.0;i--) {
      //Effect of nitrogen concentration decay through the canopy
      double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sum(LAIme(_,c)));
      // Rcout<<" l"<<i<<" fsunlit: "<< fsunlit[i]<<" lai: "<< LAIme(i,c)<<" fn: "<< fn <<"\n";
      sn+=LAIme(i,c);
      SLarealayer[i] = LAIme(i,c)*fsunlit[i];
      SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
      Vmax298layer[i] = Vmax298[c]*fn;
      Jmax298layer[i] = Jmax298[c]*fn;
    }
    for(int i=0;i<nz;i++) {
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
  supply.attr("names") = above.attr("row.names");
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
      supply[c] = supplyFunctionNetwork(psic,
                                        VGrhizo_kmaxc,VG_nc,VG_alphac,
                                        VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                                                          VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                                                                                                               VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c],
                                                                                                                                                    PLCStemPrev, 0.0, maxNsteps, 
                                                                                                                                                    ntrial, psiTol, ETol, 0.001); 
    } else {
      Rcout<<"D";
    }
  }
  
  
  //Transpiration and photosynthesis
  NumericVector psiBk(nlayers);
  for(int l=0;l<nlayers;l++) psiBk[l] = psiVec[l]; //Store initial soil water potential
  NumericMatrix K(numCohorts, nlayers);
  NumericVector Eplant(numCohorts, 0.0), Anplant(numCohorts, 0.0);
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
  NumericMatrix An_SL(numCohorts, ntimesteps);
  NumericMatrix An_SH(numCohorts, ntimesteps);
  NumericMatrix Psi_SL(numCohorts, ntimesteps);
  NumericMatrix Psi_SH(numCohorts, ntimesteps);
  NumericMatrix Ci_SL(numCohorts, ntimesteps);
  NumericMatrix Ci_SH(numCohorts, ntimesteps);
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
  NumericVector minPsiLeaf(numCohorts,0.0), maxPsiLeaf(numCohorts,-99999.0); 
  NumericVector minPsiLeaf_SL(numCohorts,0.0), maxPsiLeaf_SL(numCohorts,-99999.0); 
  NumericVector minPsiLeaf_SH(numCohorts,0.0), maxPsiLeaf_SH(numCohorts,-99999.0);
  NumericVector minPsiStem(numCohorts, 0.0), minPsiRoot(numCohorts,0.0); //Minimum potentials experienced
  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), 0.0);
  NumericMatrix PLC(numCohorts, ntimesteps);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts);
  NumericVector dEdPm(numCohorts);
  
  List outPhotoSunlit(numCohorts);
  List outPhotoShade(numCohorts);
  List outPMSunlit(numCohorts);
  List outPMShade(numCohorts);
  outPhotoSunlit.attr("names") = above.attr("row.names");
  outPhotoShade.attr("names") = above.attr("row.names");
  outPMSunlit.attr("names") = above.attr("row.names");
  outPMShade.attr("names") = above.attr("row.names");
  
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
                                                                    absSWR_SL[c] + LWR_emmcan*LAI_SL[c], 
                                                                    irradianceToPhotonFlux(absPAR_SL[c]), 
                                                                    Vmax298SL[c], 
                                                                    Jmax298SL[c], 
                                                                    Gwminc, Gwmax[c], leafWidth[c], LAI_SL[c]);
            DataFrame photoShade = leafPhotosynthesisFunction(fittedE, Catm, Patm,Tcan[n], vpatm, 
                                                              zWind[c], 
                                                                   absSWR_SH[c] + LWR_emmcan*LAI_SH[c], 
                                                                   irradianceToPhotonFlux(absPAR_SH[c]),
                                                                   Vmax298SH[c], 
                                                                   Jmax298SH[c], 
                                                                   Gwminc, Gwmax[c], leafWidth[c], LAI_SH[c]);
            
            NumericVector AnSunlit = photoSunlit["NetPhotosynthesis"];
            NumericVector AnShade = photoShade["NetPhotosynthesis"];
            NumericVector GwSunlit = photoSunlit["WaterVaporConductance"];
            NumericVector GwShade = photoShade["WaterVaporConductance"];
            NumericVector CiSunlit = photoSunlit["Ci"];
            NumericVector CiShade = photoShade["Ci"];
            NumericVector VPDSunlit = photoSunlit["LeafVPD"];
            NumericVector VPDShade = photoShade["LeafVPD"];
            NumericVector TempSunlit = photoSunlit["LeafTemperature"];
            NumericVector TempShade = photoShade["LeafTemperature"];
            

            //Profit maximization
            List PMSunlit = profitMaximization(sFunction, photoSunlit,  hydraulicCostFunction, Gwminc, Gwmax[c], VCstem_kmax[c]);
            List PMShade = profitMaximization(sFunction, photoShade,  hydraulicCostFunction, Gwminc,Gwmax[c], VCstem_kmax[c]);
            int iPMSunlit = PMSunlit["iMaxProfit"];
            int iPMShade = PMShade["iMaxProfit"];
            
            //Store?
            if(!IntegerVector::is_na(stepFunctions)) {
              if(n==stepFunctions) {
                outPhotoSunlit[c] = photoSunlit;
                outPhotoShade[c] = photoShade;
                outPMSunlit[c] = PMSunlit;
                outPMShade[c] = PMShade;
              }
            }
            // Rcout<<iPMSunlit<<" "<<iPMShade<<"\n";
            //Get leaf status
            Psi_SH(c,n) = psiLeaf[iPMShade];
            Psi_SL(c,n) = psiLeaf[iPMSunlit];
            An_SH(c,n) = AnShade[iPMShade];
            An_SL(c,n) = AnSunlit[iPMSunlit];
            Ci_SH(c,n) = CiShade[iPMShade];
            Ci_SL(c,n) = CiSunlit[iPMSunlit];
            GW_SH(c,n)= GwShade[iPMShade];
            GW_SL(c,n)= GwSunlit[iPMSunlit];
            VPD_SH(c,n)= VPDShade[iPMShade];
            VPD_SL(c,n)= VPDSunlit[iPMSunlit];
            Temp_SH(c,n)= TempShade[iPMShade];
            Temp_SL(c,n)= TempSunlit[iPMSunlit];
            
            //Scale photosynthesis
            double Ansum = AnSunlit[iPMSunlit]*LAI_SL[c] + AnShade[iPMShade]*LAI_SH[c];
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
            Anplant[c] +=Aninst(c,n);
            
            
            //Update symplastic storage and PLC
            psiRootVEC[c] = psiRoot[iPM]; 
            EinstVEC[c] = fittedE[iPM];
            PLC(c,n) = NA_REAL;
            RWCsteminst(c,n) = NA_REAL;
            psiLeafVEC[c] = psiLeaf[iPM];
            RWCsleafVEC[c] = symplasticRelativeWaterContent(psiLeafVEC[c], LeafPI0[c], LeafEPS[c]);
            psiStemMAT(c,_) = newPsiStem(iPM,_);
            for(int i=0;i<nStemSegments;i++) {
              double psiPLC; // Store the PLC corresponding to the upstream water potential
              if(i==0) psiPLC = psiRootVEC[c];
              else psiPLC = psiStemMAT(c, i-1);
              if(!cavitationRefill) PLCstemMAT(c,i) = std::max(PLCstemMAT(c,i), 1.0 - xylemConductance(psiPLC, 1.0, VCstem_c[c], VCstem_d[c]));
              else PLCstemMAT(c,i) = 1.0 - xylemConductance(psiPLC, 1.0, VCstem_c[c], VCstem_d[c]);
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
            double leafTempSL = leafTemperature(absRadSL/ LAI_SL[c], Tcan[n], zWind[c], ESLbk, leafWidth[c]);
            double leafVPDSL = (meteoland::utils_saturationVP(std::max(0.0,leafTempSL))-vpatm);
            double ESL = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSL)/Patm);
            while(std::abs(ESL-ESLbk)> 0.000001) {
              leafTempSL = leafTemperature(absRadSL/ LAI_SL[c], Tcan[n], zWind[c], ESL, leafWidth[c]);
              leafVPDSL = meteoland::utils_saturationVP(std::max(0.0,leafTempSL))-vpatm;
              ESLbk = ESL;
              ESL = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSL)/Patm);
            }
            NumericVector LP = leafphotosynthesis(QSL, Catm, Gwmin[c]/1.6, std::max(0.0,leafTempSL), Vmax298SL[c], Jmax298SL[c]);
            double CiSL = LP[0];
            double AgSL = LP[1];
            double AnSL = AgSL - 0.015*VmaxTemp(Vmax298SL[c]/LAI_SL[c], leafTempSL);
            
            //Shade photosynthesis
            double absRadSH = absSWR_SH[c] + LWR_emmcan*LAI_SH(c,n);
            double QSH = irradianceToPhotonFlux(absPAR_SH[c]);
            //Flow (cuticular conductance): iterative procedure to find flow, VPD and leaf temperature given Gwmin
            double ESHbk = std::max(0.0,1000.0*(Gwmin[c]*(meteoland::utils_saturationVP(std::max(0.0,Tcan[n]))-vpatm))/Patm);
            double leafTempSH = leafTemperature(absRadSH/ LAI_SH[c], Tcan[n], zWind[c], ESHbk, leafWidth[c]);
            double leafVPDSH = (meteoland::utils_saturationVP(std::max(0.0,leafTempSH))-vpatm);
            double ESH = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSH)/Patm);
            while(std::abs(ESH-ESHbk)> 0.000001) {
              leafTempSH = leafTemperature(absRadSH/ LAI_SH[c], Tcan[n], zWind[c], ESH, leafWidth[c]);
              leafVPDSH = meteoland::utils_saturationVP(std::max(0.0,leafTempSH))-vpatm;
              ESHbk = ESH;
              ESH = std::max(0.0,1000.0*(Gwmin[c]*leafVPDSH)/Patm);
            }
            LP = leafphotosynthesis(QSH, Catm, Gwmin[c]/1.6, std::max(0.0,leafTempSH), Vmax298SH[c], Jmax298SH[c]);
            double CiSH = LP[0];
            double AgSH = LP[1];
            double AnSH = AgSH - 0.015*VmaxTemp(Vmax298SH[c]/LAI_SH[c], leafTempSH);
            
            //Average flow
            double Eaverage = (ESL*LAI_SL[c] + ESH*LAI_SH[c])/(LAI_SL[c] + LAI_SH[c]);
            
            //Get leaf status
            An_SH(c,n) = AgSH;
            An_SL(c,n) = AnSL;
            Ci_SH(c,n) = CiSH;
            Ci_SL(c,n) = CiSL;
            GW_SH(c,n)= Gwmin[c];
            GW_SL(c,n)= Gwmin[c];
            VPD_SH(c,n)= leafVPDSH;
            VPD_SL(c,n)= leafVPDSL;
            Temp_SH(c,n)= leafTempSH;
            Temp_SL(c,n)= leafVPDSL;
            
            //Scale photosynthesis
            double Ansum = AnSL*LAI_SL[c] + AnSH*LAI_SH[c];
            Aninst(c,n) =(1e-6)*12.01017*Ansum*tstep; 

            Einst(c,n) = Eaverage*0.001*0.01802*LAIphe[c]*tstep; //Scale from instantaneous flow to water volume in the time step
            
            
            //Add to daily plant cohort transpiration
            Eplant[c] +=Einst(c,n);
            Anplant[c] +=Aninst(c,n);
            
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
            
            
            
            //As it is disconnected, leaf and shade leaves have the same wp
            Psi_SH(c,n) =  sAb["psiLeaf"];
            Psi_SL(c,n) =  sAb["psiLeaf"];
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
            Psi_SH(c,n) = NA_REAL;
            Psi_SL(c,n) = NA_REAL;
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
      
      PsiLeafinst(c,n) = psiLeafVEC[c]; //Store instantaneous (average) leaf potential
      PsiRootinst(c,n) = psiRootVEC[c]; //Store instantaneous root crown potential
      
      //Store the minimum water potential of the day (i.e. mid-day)
      minPsiLeaf_SL[c] = std::min(minPsiLeaf_SL[c],Psi_SL(c,n));
      minPsiLeaf_SH[c] = std::min(minPsiLeaf_SH[c],Psi_SH(c,n));
      maxPsiLeaf_SL[c] = std::max(maxPsiLeaf_SL[c],Psi_SL(c,n));
      maxPsiLeaf_SH[c] = std::max(maxPsiLeaf_SH[c],Psi_SH(c,n));
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
    //Latent heat (evaporation of intercepted + transpiration)
    double canEvapStep = abs_SWR_can[n]*(canopyEvaporation/sum(abs_SWR_can));
    double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tcan[n])*(sum(Einst(_,n)) + canEvapStep)/tstep;
    LEcan_heat[n] = LEwat; 
    //Canopy longwave emmission
    LWRcanout[n] = LWR_emmcan*canLWRexchprop;
    //Canopy convective heat exchange
    Hcan_heat[n] = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;
    
    //Soil latent heat (soil evaporation)
    //Latent heat (snow fusion) as J/m2/s
    double soilEvapStep = abs_SWR_soil[n]*(soilEvaporation/sum(abs_SWR_soil));
    double LEsoilevap = (1e6)*meteoland::utils_latentHeatVaporisation(Tsoil[0])*soilEvapStep/tstep;
    double snowMeltStep = abs_SWR_soil[n]*(snowMelt/sum(abs_SWR_soil));
    double LEsnow = (1e6)*(snowMeltStep*0.33355)/tstep; // 0.33355 = latent heat of fusion
    LEsoil_heat[n] = LEsoilevap + LEsnow;
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
    Ebalsoil[n] = abs_SWR_soil[n] + abs_LWR_soil[n] + LWRcanout[n] + Hcansoil[n] - LEsoil_heat[n] - LWRsoilout[n]; //Here we use all energy escaping to atmosphere
    //Soil temperature changes
    NumericVector soilTchange = temperatureChange(dVec, Tsoil, sand, clay, W, Theta_FC, Ebalsoil[n]);
    for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tstep);
    if(n<(ntimesteps-1)) Tsoil_mat(n+1,_)= Tsoil;
    
    // Rcout<<n<<", Tatm: "<< Tatm[n]<< " Tcan: "<<Tcan[n]<<" soilT1 "<< Tsoil[0]<<"\n";
    
    
    //save canopy temperature
    canopyParams["Temp"] = Tcannext;
    
  } //End of timestep loop

  //4z. Plant daily drought stress (from root collar mid-day water potential)
  NumericVector SoilExtractCoh(numCohorts,0.0);
  NumericVector DDS(numCohorts, 0.0);
  for(int c=0;c<numCohorts;c++) {
    SoilExtractCoh[c] =  sum(SoilWaterExtract(c,_));
    transpiration[c] = Eplant[c]; 
    photosynthesis[c] = Anplant[c];
    PLCm[c] = sum(PLC(c,_))/((double)PLC.ncol());
    RWCsm[c] = sum(RWCsteminst(c,_))/((double)RWCsteminst.ncol());
    RWClm[c] = sum(RWCleafinst(c,_))/((double)RWCleafinst.ncol());
    dEdPm[c] = sum(dEdPinst(c,_))/((double)dEdPinst.ncol());  
    double maxConductance = maximumSoilPlantConductance(VGrhizo_kmax(c,_), VCroot_kmax(c,_), VCstem_kmax[c], VCleaf_kmax[c]);
    DDS[c] = Phe[c]*(1.0 - (dEdPm[c]/maxConductance));
  }
  
  
  //B.3 - Substract  extracted water from soil moisture 
  if(modifyInput){
    for(int l=0;l<nlayers;l++) {
      W[l] = std::max(W[l] - (sum(soilLayerExtractInst(l,_))/Water_FC[l]),0.0);
    } 
  }
  
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
                                        _["Hcansoil"] = Hcansoil, _["LEsoil"] = LEsoil_heat, _["SWRsoilin"] = abs_SWR_soil, _["LWRsoilin"] = abs_LWR_soil,  _["LWRsoilout"] = LWRsoilout,
                                        _["Ebalsoil"] = Ebalsoil, _["RAsoil"] = RAsoil);
  List EB = List::create(_["Temperature"]=Tinst, _["CanopyEnergyBalance"] = CEBinst, _["SoilEnergyBalance"] = SEBinst);
  Psi_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GW_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GW_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));

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
  for(int c=0;c<numCohorts;c++) {
    Vmax298SL[c] = Vmax298SL[c]/LAI_SL[c];
    Jmax298SL[c] = Jmax298SL[c]/LAI_SL[c];
    Vmax298SH[c] = Vmax298SH[c]/LAI_SH[c];
    Jmax298SH[c] = Jmax298SH[c]/LAI_SH[c];
  }
  DataFrame Sunlit = DataFrame::create(
    _["LAI"] = LAI_SL, 
    _["Vmax298"] = Vmax298SL,
    _["Jmax298"] = Jmax298SL
  );
  DataFrame Shade = DataFrame::create(
    _["LAI"] = LAI_SH, 
    _["Vmax298"] = Vmax298SH,
    _["Jmax298"] = Jmax298SH
  );
  Sunlit.attr("row.names") = above.attr("row.names");
  Shade.attr("row.names") = above.attr("row.names");
  
  List ShadeInst = List::create(
    _["Abs_SWR"] = SWR_SH,
    _["Abs_LWR"] = LWR_SH,
    _["An"] = An_SH,
    _["Ci"] = Ci_SH,
    _["GW"] = GW_SH,
    _["VPD"] = VPD_SH,
    _["Temp"] = Temp_SH,
    _["Psi"] = Psi_SH);
  List SunlitInst = List::create(
    _["Abs_SWR"]=SWR_SL,
    _["Abs_LWR"] = LWR_SL,
    _["An"] = An_SL,
    _["Ci"] = Ci_SL,
    _["GW"] = GW_SL,
    _["VPD"] = VPD_SL,
    _["Temp"] = Temp_SL,
    _["Psi"] = Psi_SL);
  
  List PlantsInst = List::create(
    _["E"]=Einst, _["An"]=Aninst,
    _["SunlitLeaves"] = SunlitInst,
    _["ShadeLeaves"] = ShadeInst,
    _["dEdPinst"] = dEdPinst,
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
                                       _["Photosynthesis"] = Anplant,
                                       _["RootPsi"] = minPsiRoot, 
                                       _["StemPsi"] = minPsiStem, 
                                       _["StemPLC"] = PLCm, //Average daily stem PLC
                                       _["LeafPsiMin"] = minPsiLeaf, 
                                       _["LeafPsiMax"] = maxPsiLeaf, 
                                       _["LeafPsiMin_SL"] = minPsiLeaf_SL, 
                                       _["LeafPsiMax_SL"] = maxPsiLeaf_SL, 
                                       _["LeafPsiMin_SH"] = minPsiLeaf_SH, 
                                       _["LeafPsiMax_SH"] = maxPsiLeaf_SH, 
                                       _["dEdP"] = dEdPm,//Average daily soilplant conductance
                                       _["DDS"] = DDS, //Daily drought stress is the ratio of average soil plant conductance over its maximum value
                                       _["StemRWC"] = RWCsm,
                                       _["LeafRWC"] = RWClm);
  Plants.attr("row.names") = above.attr("row.names");
  
  List l;
  if(!IntegerVector::is_na(stepFunctions)){
    l = List::create(_["cohorts"] = clone(cohorts),
                     _["EnergyBalance"] = EB,
                     _["Extraction"] = SoilWaterExtract,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["RhizoPsi"] = minPsiRhizo,
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["PlantsInst"] = PlantsInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["WindExtinction"] = zWind,
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
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["PlantsInst"] = PlantsInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["WindExtinction"] = zWind,
                     _["SupplyFunctions"] = supply);
    
  } 
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}

// [[Rcpp::export("transp_transpirationSperry")]]
List transpirationSperry(List x, List soil, DataFrame meteo, int day,
                        double latitude, double elevation, double slope, double aspect,
                        double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                        int stepFunctions = NA_INTEGER, bool modifyInput = true) {
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
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);

  return(transpirationSperry(x,soil, tmin, tmax, rhmin, rhmax, rad, wind, 
                     latitude, elevation, slope, aspect,
                     solarConstant, delta, prec,
                     canopyEvaporation, snowMelt, soilEvaporation,
                     false, stepFunctions, modifyInput));
} 


List transpirationGranier(List x, List soil, double tday, double pet, bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
  bool cavitationRefill = control["cavitationRefill"];
  String soilFunctions = control["soilFunctions"];
  
  //Soil input
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector psiSoil = psi(soil,soilFunctions); //Update soil water potential
  NumericVector Water_FC = waterFC(soil, soilFunctions);
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
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Psi_Extract"]);
  NumericVector WUE = Rcpp::as<Rcpp::NumericVector>(paramsTransp["WUE"]);
  NumericVector pRootDisc = Rcpp::as<Rcpp::NumericVector>(paramsTransp["pRootDisc"]);
  
  //Communication vectors
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericVector pEmb = clone(Rcpp::as<Rcpp::NumericVector>(x["PLC"]));
  
  
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector Phe(numCohorts,0.0);
  double s = 0.0, LAIcell = 0.0, LAIcelldead = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(LAIlive[c]>0) Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    else Phe[c]=0.0;
    s += (kPAR[c]*(LAIphe[c]+LAIdead[c]));
    LAIcell += LAIphe[c]+LAIdead[c];
    LAIcelldead += LAIdead[c];
  }
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  LAIdead, H, CR, kPAR);

  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  double Tmax = pet*(-0.006*pow(LAIcell,2.0)+0.134*LAIcell); //From Granier (1999)
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts,0.0);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  NumericMatrix EplantCoh(numCohorts, nlayers);
  NumericMatrix PsiRoot(numCohorts, nlayers);
  NumericVector PlantPsi(numCohorts, NA_REAL);
  NumericVector Eplant(numCohorts, 0.0), Anplant(numCohorts, 0.0);
  NumericVector DDS(numCohorts, 0.0);
  NumericVector Kl, epc, Vl;
  double WeibullShape=3.0;
  for(int l=0;l<nlayers;l++) {
    Kl = Psi2K(psiSoil[l], Psi_Extract, WeibullShape);
    
    //Limit Kl due to previous cavitation
    if(!cavitationRefill) for(int c=0;c<numCohorts;c++) Kl[c] = std::min(Kl[c], 1.0-pEmb[c]);
    //Limit Kl to minimum value for root disconnection
    Vl = V(_,l);
    epc = pmax(TmaxCoh*Kl*Vl,0.0);
    for(int c=0;c<numCohorts;c++) {
      PsiRoot(c,l) = psiSoil[l]; //Set initial guess of root potential to soil values
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
    DDS = DDS + Phe*(Vl*(1.0 - Kl)); //Add stress from the current layer
  }
  
  double alpha = std::max(std::min(tday/20.0,1.0),0.0);
  for(int c=0;c<numCohorts;c++) {
    PlantPsi[c] = averagePsi(PsiRoot(c,_), V(c,_), WeibullShape, Psi_Extract[c]);
    if(!cavitationRefill) {
      pEmb[c] = std::max(DDS[c], pEmb[c]); //Track current embolism if no refill
      DDS[c] = pEmb[c];
    }
    Anplant[c] = alpha*WUE[c]*Eplant[c];
  }
  
  
  if(modifyInput) {
    for(int c=0;c<numCohorts;c++) {
      transpiration[c] = Eplant[c];
      photosynthesis[c] = Anplant[c];
    }
    x["PLC"] = pEmb;
  }
  
  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  
  if(modifyInput) {
    for(int l=0;l<nlayers;l++) {
      W[l] = W[l] - (sum(EplantCoh(_,l))/Water_FC[l]); 
    }
  }
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["Transpiration"] = Eplant, 
                                       _["Photosynthesis"] = Anplant,
                                       _["psi"] = PlantPsi, _["DDS"] = DDS);
  Plants.attr("row.names") = above.attr("row.names");
  EplantCoh.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Plants"] = Plants,
                        _["Extraction"] = EplantCoh);
  return(l);
}


// [[Rcpp::export("transp_transpirationGranier")]]
List transpirationGranier(List x, List soil, DataFrame meteo, int day,
                          bool modifyInput = true) {
  if(!meteo.containsElementNamed("MeanTemperature")) stop("Please include variable 'MeanTemperature' in weather input.");
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  if(!meteo.containsElementNamed("PET")) stop("Please include variable 'PET' in weather input.");
  NumericVector PET = meteo["PET"];
  double pet = PET[day-1];
  double tday = MeanTemperature[day-1];
  return(transpirationGranier(x,soil, tday, pet, modifyInput));
} 

