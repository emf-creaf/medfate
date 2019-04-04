#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "photosynthesis.h"
#include "soil.h"
#include <Rcpp.h>
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;

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


// [[Rcpp::export("transp_stomatalRegulation")]]
List stomatalRegulation(List x, List soil, DataFrame meteo, int day,
                        double latitude, double elevation) {
  
  //Extract meteo
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector MinTemperature = meteo["MinTemperature"];
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  NumericVector Radiation = meteo["Radiation"];
  NumericVector WindSpeed = meteo["WindSpeed"];
  
  CharacterVector dateStrings = meteo.attr("row.names");
  IntegerVector DOY = date2doy(dateStrings);
  
  NumericVector GDD = gdd(DOY, MeanTemperature, 5.0);
  
  
  double rad = Radiation[day-1];
  double rain = Precipitation[day-1];
  double tmin = MinTemperature[day-1];
  double tmax = MaxTemperature[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double rhmax = MaxRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
    
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  
  //Control parameters
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  List numericParams = control["numericParams"];
  int ntrial = numericParams["ntrial"];
  int maxNsteps  = numericParams["maxNsteps"];
  double psiTol = numericParams["psiTol"];
  double ETol = numericParams["ETol"];
  
  int ntimesteps = control["ndailysteps"];
  int hydraulicCostFunction = control["hydraulicCostFunction"];
  double verticalLayerSize = control["verticalLayerSize"];
  
  
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector albedo = Rcpp::as<Rcpp::NumericVector>(paramsBase["albedo"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  
  //Soil input
  NumericVector W = soil["W"];
  NumericVector dVec = soil["dVec"];
  NumericVector psiVec = psi(soil, soilFunctions);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  NumericVector om = soil["om"];
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
  
  
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  NumericVector photosynthesis = Rcpp::as<Rcpp::NumericVector>(x["Photosynthesis"]);
  NumericMatrix PLCstemMAT = Rcpp::as<Rcpp::NumericMatrix>(x["PLCstem"]);
  
  //Canopy params
  List canopyParams = Rcpp::as<Rcpp::List>(x["canopy"]);
  
  //Root distribution input
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(below["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VGrhizo_kmax"]);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  
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
  
  
  //Leaf phenology and the adjusted leaf area index
  double canopyHeight = 0.0, LAIcell= 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(canopyHeight<H[c]) canopyHeight = H[c];
    LAIcell += (LAIphe[c]+LAIdead[c]);
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  
  double latrad = latitude * (PI/180.0);

  NumericVector z(nz+1,0.0);
  NumericVector zmid(nz);
  for(int i=1;i<=nz;i++) {
    z[i] = z[i-1] + verticalLayerSize;
    zmid[i-1] = (verticalLayerSize/2.0) + verticalLayerSize*((double) (i-1));
  }
  NumericMatrix LAIme = LAIdistribution(z, LAIphe, H, CR); //Expanded leaves
  NumericMatrix LAImd = LAIdistribution(z, LAIdead, H, CR); //Dead (standing) leaves
  NumericMatrix LAImx = LAIdistribution(z, LAIlive, H, CR); //Maximum leaf expansion
  

  //Day length (latitude in radians), atmospheric pressure, CO2 concentration
  double tauday = meteoland::radiation_daylengthseconds(latrad,  0.0,0.0, delta); 
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  double Catm = 386.0;
  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);
  
  //Daily average water vapor pressure
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  
  //Daily cloud cover
  double cloudcover = 0.0;
  if(rain >0.0) cloudcover = 1.0;
  bool clearday = (rain==0);
  
  //Instantaneous direct and diffuse shorwave radiation
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, delta,
                                                        rad, clearday, ntimesteps);
  NumericVector solarElevation = ddd["SolarElevation"]; //in radians
  NumericVector solarHour = ddd["SolarHour"]; //in radians
  NumericVector SWR_direct = ddd["SWR_direct"]; //in kW·m-2
  NumericVector SWR_diffuse = ddd["SWR_diffuse"]; //in kW·m-2
  NumericVector PAR_direct = ddd["PAR_direct"]; //in kW·m-2
  NumericVector PAR_diffuse = ddd["PAR_diffuse"]; //in kW·m-2
  
  //Instantaneous air temperature and longwave radiation
  NumericVector Tatm(ntimesteps), lwdr(ntimesteps), Tcan(ntimesteps, NA_REAL);
  for(int n=0;n<ntimesteps;n++) {
    //From solar hour (radians) to seconds from sunrise
    double t_sunrise = (solarHour[n]*43200.0/PI)+ (tauday/2.0) +(tstep/2.0); 
    
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern(t_sunrise, tmin, tmax, tauday);
    //Longwave sky diffuse radiation
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpatm, cloudcover);
  }
  Tcan[0] = canopyParams["Temp"]; //Take canopy temperature from previous day
  

  //Light extinction and absortion by time steps
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        kPAR, albedo,
                                                                        ddd, lwdr,
                                                                        ntimesteps,  "sunshade", 0.1);
  List abs_PAR_SL_list = lightExtinctionAbsortion["PAR_SL"];
  List abs_PAR_SH_list = lightExtinctionAbsortion["PAR_SH"];
  List abs_SWR_SL_list = lightExtinctionAbsortion["SWR_SL"];
  List abs_SWR_SH_list = lightExtinctionAbsortion["SWR_SH"];
  List abs_LWR_SL_list = lightExtinctionAbsortion["LWR_SL"];
  List abs_LWR_SH_list = lightExtinctionAbsortion["LWR_SH"];
  NumericVector fsunlit = lightExtinctionAbsortion["fsunlit"];

  NumericMatrix Einst(numCohorts, ntimesteps);
  NumericMatrix Aninst(numCohorts, ntimesteps);
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
  

  //Wind extinction profile
  NumericVector zWind;
  if(NumericVector::is_na(wind)) wind = 2.0; //Default wind speed when missing
  zWind = windExtinctionCohort(H,CR, wind,LAIcell, canopyHeight);
  
  
  NumericVector Vc, VCroot_kmaxc, VGrhizo_kmaxc, psic, VG_nc,VG_alphac;
  
  List cohort_list(numCohorts);
  cohort_list.attr("names") = above.attr("row.names");
  for(int c=0;c<numCohorts;c++) {
    //Determine to which layers is plant connected
    LogicalVector layerConnected(nlayers, false);
    int nlayersc = 0;
    for(int l=0;l<nlayers;l++) {
      if(V(c,l)>0.0) {
        double pRoot = xylemConductance(psiVec[l], 1.0, VCroot_c[c], VCroot_d[c]); //Relative conductance in the root
        layerConnected[l]= (pRoot>=pRootDisc[c]);
        if(layerConnected[l]) nlayersc++;
      }
    }
    // Rcout<<nlayersc;
    Vc = NumericVector(nlayersc);
    VCroot_kmaxc = NumericVector(nlayersc);
    VGrhizo_kmaxc = NumericVector(nlayersc);
    psic = NumericVector(nlayersc);
    VG_nc = NumericVector(nlayersc);
    VG_alphac= NumericVector(nlayersc);
    int cnt=0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected[l]) {
        Vc[cnt] = V(c,l);
        VCroot_kmaxc[cnt] = VCroot_kmax(c,l);
        VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l);
        psic[cnt] = psiVec[l];
        VG_nc[cnt] = VG_n[l];
        VG_alphac[cnt] = VG_alpha[l];
        cnt++;
      }
    }
    
    double minFlow = 1000.0*(Gwmin[c]*(tmin+tmax)/2.0)/Patm;
    NumericVector PLCStemPrev = PLCstemMAT(c,_); //Get row
    List supply = supplyFunctionNetwork(psic,
                                        VGrhizo_kmaxc,VG_nc,VG_alphac,
                                        VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                        VCstem_kmax[c], VCstem_c[c],VCstem_d[c], 
                                        VCleaf_kmax[c], VCleaf_c[c],VCleaf_d[c], 
                                        PLCStemPrev,                                    
                                        minFlow, maxNsteps,
                                        ntrial, psiTol, ETol);
    NumericVector fittedE = supply["E"];
    
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

    List photoSunlit(ntimesteps);
    List photoShade(ntimesteps);
    List PMSunlit(ntimesteps);
    List PMShade(ntimesteps);
    
    photosynthesis[c] = 0.0;
    transpiration[c] = 0.0;
    
      
    for(int n=0;n<ntimesteps;n++) {

      //Long-wave radiation due to canopy temperature
      if(NumericVector::is_na(Tcan[n])) Tcan[n] = Tatm[n]; //If missing take above-canopy air temperature

      //LWR emmited by the canopy, per ground area
      double LWR_emmcan = 0.95*SIGMA_Wm2*pow(Tcan[n]+273.16,4.0);
      
      // Rcout<< n<< " sh "<<solarHour[n]<< " t_sr " << t_sunrise<< " Tatm "<<Tatm<<"\n";
      NumericVector absPAR_SL = abs_PAR_SL_list[n];
      NumericVector absPAR_SH = abs_PAR_SH_list[n];
      NumericVector absSWR_SL = abs_SWR_SL_list[n];
      NumericVector absSWR_SH = abs_SWR_SH_list[n];
      NumericVector absLWR_SL = abs_LWR_SL_list[n];
      NumericVector absLWR_SH = abs_LWR_SH_list[n];
      
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
      SWR_SL(c,n) = absSWR_SL[c];
      SWR_SH(c,n) = absSWR_SH[c];
      LWR_SL(c,n) = absLWR_SL[c];
      LWR_SH(c,n) = absLWR_SH[c];
      

      //Photosynthesis function for sunlit and shade leaves
      DataFrame psl = leafPhotosynthesisFunction(fittedE, Catm, Patm,Tcan[n], vpatm, 
                                                    zWind[c], 
                                                         absSWR_SL[c] + LWR_emmcan*LAI_SL(c,n), 
                                                         irradianceToPhotonFlux(absPAR_SL[c]), 
                                                         Vmax298SL, 
                                                         Jmax298SL, 
                                                         Gwmin[c], Gwmax[c], leafWidth[c], LAI_SL(c,n));
      DataFrame psh = leafPhotosynthesisFunction(fittedE, Catm, Patm,Tcan[n], vpatm, 
                                                   zWind[c], 
                                                        absSWR_SH[c] + LWR_emmcan*LAI_SH(c,n), 
                                                        irradianceToPhotonFlux(absPAR_SH[c]),
                                                        Vmax298SH, 
                                                        Jmax298SH, 
                                                        Gwmin[c], Gwmax[c], leafWidth[c], LAI_SH(c,n));

      NumericVector AnSunlit = psl["NetPhotosynthesis"];
      NumericVector AnShade = psh["NetPhotosynthesis"];
      NumericVector GwSunlit = psl["WaterVaporConductance"];
      NumericVector GwShade = psh["WaterVaporConductance"];
      NumericVector VPDSunlit = psl["LeafVPD"];
      NumericVector VPDShade = psh["LeafVPD"];
      NumericVector TempSunlit = psl["LeafTemperature"];
      NumericVector TempShade = psh["LeafTemperature"];

      //Profit maximization
      List pmsl = profitMaximization(supply, psl,  hydraulicCostFunction, Gwmin[c],Gwmax[c], VCstem_kmax[c]);
      List pmsh = profitMaximization(supply, psh,  hydraulicCostFunction, Gwmin[c],Gwmax[c], VCstem_kmax[c]);

      PMSunlit[n] =pmsl;
      PMShade[n] = pmsh;
      photoSunlit[n] = psl;
      photoShade[n] = psh;
      
      int iPMSunlit = pmsl["iMaxProfit"];
      int iPMShade = pmsh["iMaxProfit"];
      
      GW_SH(c,n)= GwShade[iPMShade];
      GW_SL(c,n)= GwSunlit[iPMSunlit];
      VPD_SH(c,n)= VPDShade[iPMShade];
      VPD_SL(c,n)= VPDSunlit[iPMSunlit];
      Temp_SH(c,n)= TempShade[iPMShade];
      Temp_SL(c,n)= TempSunlit[iPMSunlit];
      
      //Scale photosynthesis
      Aninst(c,n) = AnSunlit[iPMSunlit]*LAI_SL(c,n) + AnShade[iPMShade]*LAI_SH(c,n);
      photosynthesis[c] +=(1e-6)*12.01017*Aninst(c,n)*tstep; 
      
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
      
      //Scale from instantaneous flow to water volume in the time step
      Einst(c,n) = fittedE[iPM]*0.001*0.01802*LAIphe[c]*tstep; 
      transpiration[c] +=Einst(c,n);
    }

    //for each cohort, table with time steps in rows and PM results in columns for sunlit and shade leaves
    //table with cohorts in rows and daily PM results in columns, with units scaled by leaf area
    cohort_list[c] = List::create(_["supply"]=supply,
                                  _["photoSunlit"]=photoSunlit,
                                  _["photoShade"]=photoShade,
                                  _["PMSunlit"] = PMSunlit,
                                  _["PMShade"] = PMShade);
  }
  Rcout<<"done\n";
  
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
  
  List PlantsInst = List::create(
    _["LAIsunlit"] = LAI_SL, _["LAIshade"] = LAI_SH, 
    _["AbsRad"] = AbsRadinst, _["E"]=Einst, _["An"]=Aninst,
    _["GWsunlit"] = GW_SL, _["GWshade"] = GW_SH,
    _["VPDsunlit"] = VPD_SL, _["VPDshade"] = VPD_SH,
    _["Tempsunlit"] = Temp_SL, _["Tempshade"] = Temp_SH);
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["PlantsInst"] = PlantsInst,
                        _["CohortDetails"] = cohort_list);
  return(l);
} 

// [[Rcpp::export("transp_Granier")]]
List transpGranier(List x, List soil, double tday, double pet) {
  //Control parameters
  List control = x["control"];
  bool cavitationRefill = control["cavitationRefill"];
  String soilFunctions = control["soilFunctions"];
  
  //Soil input
  NumericVector psiSoil = psi(soil,soilFunctions); //Update soil water potential
  
  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIphe.size();
  int nlayers = psiSoil.length();
  
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
  NumericVector Eplant(numCohorts, 0.0);
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
  for(int c=0;c<numCohorts;c++) {
    PlantPsi[c] = averagePsi(PsiRoot(c,_), V(c,_), WeibullShape, Psi_Extract[c]);
    if(!cavitationRefill) {
      pEmb[c] = std::max(DDS[c], pEmb[c]); //Track current embolism if no refill
      DDS[c] = pEmb[c];
    }
  }
  
  double alpha = std::max(std::min(tday/20.0,1.0),0.0);
  //For comunication with growth and landscape water balance
  for(int c=0;c<numCohorts;c++) {
    transpiration[c] = Eplant[c];
    photosynthesis[c] = alpha*WUE[c]*transpiration[c];
  }
  
  //Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) LAIcohort[c]= LAIphe[c];
  
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["Transpiration"] = Eplant, _["psi"] = PlantPsi, _["DDS"] = DDS);
  Plants.attr("row.names") = above.attr("row.names");
  EplantCoh.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  List l = List::create(_["Plants"] = Plants,
                        _["Uptake"] = EplantCoh);
  return(l);
}
