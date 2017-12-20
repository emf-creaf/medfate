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

const double SIGMA_W = 5.67*pow(10,-8.0);

List profitMaximization1(List supplyFunction, List photosynthesisFunction, double Gwmax) {
  NumericVector supplydEdp = supplyFunction["dEdP"];
  NumericVector Ag = photosynthesisFunction["Photosynthesis"];
  NumericVector Gw = photosynthesisFunction["WaterVaporConductance"];
  int nsteps = supplydEdp.size();
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  double Agmax = 0.0;
  int nvalidsteps = 0;
  while((Gw[nvalidsteps]< Gwmax) & (nvalidsteps<nsteps)) {
    nvalidsteps++;
  }
  for(int i=0;i<nsteps;i++) {
    maxdEdp = std::max(maxdEdp, supplydEdp[i]);
    mindEdp =  std::min(mindEdp, supplydEdp[i]);
    Agmax = std::max(Agmax, Ag[i]);
  }
  //Evaluate profit for valid steps
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  for(int i=0;i<nvalidsteps;i++) {
    cost[i] = (maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp);
    gain[i] = Ag[i]/Agmax;
    profit[i] = gain[i]-cost[i];
  }
  
  int imaxprofit=0;
  double maxprofit=profit[0];
  for(int i=(imaxprofit+1);i<nvalidsteps;i++){
    if((profit[i]>maxprofit)) {
      maxprofit = profit[i];
      imaxprofit = i;
    }
  }
  // if(Gw[imaxprofit] == Gwmax) {
  //   Rcout<<"GX";
  // }
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}

List profitMaximization2(List supplyFunction, List photosynthesisFunction, double Gwmax, double kleafmax) {
  NumericVector supplyKterm = supplyFunction["kterm"];
  NumericVector Ag = photosynthesisFunction["Photosynthesis"];
  NumericVector Gw = photosynthesisFunction["WaterVaporConductance"];
  int nsteps = supplyKterm.size();
  double Agmax = 0.0;
  int nvalidsteps = 0;
  while((Gw[nvalidsteps]< Gwmax) & (nvalidsteps<nsteps)) {
    nvalidsteps++;
  }
    
  for(int i=0;i<nsteps;i++) {
    Agmax = std::max(Agmax, Ag[i]);
  }
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  //Evaluate profit for valid steps
  for(int i=0;i<nvalidsteps;i++) {
      cost[i] = (kleafmax-supplyKterm[i])/(kleafmax-0.0);
      gain[i] = Ag[i]/Agmax;
      profit[i] = gain[i]-cost[i];
  }
  
  int imaxprofit=0;
  double maxprofit=profit[0];
  for(int i=0;i<nvalidsteps;i++){
    if((profit[i]>maxprofit)) {
      maxprofit = profit[i];
      imaxprofit = i;
    }
  }
  // if(Gw[imaxprofit] == Gwmax) {
  //   Rcout<<"GX";
  // }
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}

// [[Rcpp::export("transp.profitMaximization")]]
List profitMaximization(List supplyFunction, List photosynthesisFunction, int type, double Gwmax, double kleafmax = NA_REAL) {
  if(type==1) return(profitMaximization1(supplyFunction, photosynthesisFunction, Gwmax));
  return(profitMaximization2(supplyFunction, photosynthesisFunction, Gwmax, kleafmax));
}


// [[Rcpp::export("transp.dayCanopyTranspiration")]]
List dayCanopyTranspiration(List x, List soil, DataFrame meteo, int day,
                            double latitude, double elevation, double slope = 0.0, double aspect = 0.0) {
  
  //Extract meteo
  IntegerVector DOY = meteo["DOY"];
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector MeanTemperature = meteo["MeanTemperature"];
  NumericVector MinTemperature = meteo["MinTemperature"];
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  NumericVector Radiation = meteo["Radiation"];
  NumericVector WindSpeed = meteo["WindSpeed"];
  CharacterVector dateStrings = meteo.attr("row.names");
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
  List numericParams = control["numericParams"];
  double psiStep = numericParams["psiStep"];
  double psiMax = numericParams["psiMax"];
  int ntrial = numericParams["ntrial"];
  int maxNsteps  = numericParams["maxNsteps"];
  double psiTol = numericParams["psiTol"];
  double ETol = numericParams["ETol"];
  
  String canopyMode = Rcpp::as<Rcpp::String>(control["canopyMode"]);
  bool cavitationRefill = control["cavitationRefill"];
  int ntimesteps = control["ndailysteps"];
  int hydraulicCostFunction = control["hydraulicCostFunction"];
  double verticalLayerSize = control["verticalLayerSize"];
  
  
  DataFrame paramsBase = Rcpp::as<Rcpp::DataFrame>(x["paramsBase"]);
  NumericVector albedo = Rcpp::as<Rcpp::NumericVector>(paramsBase["albedo"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsBase["k"]);
  
  //Soil input
  NumericVector W = soil["W"];
  NumericVector psi = soil["psi"];
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  NumericVector om = soil["om"];
  NumericVector Water_FC = waterFC(soil);
  int nlayers = W.size();
  
  //Vegetation input
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  int numCohorts = LAIlive.size();
  NumericVector pEmb = Rcpp::as<Rcpp::NumericVector>(x["ProportionCavitated"]);
  
  //Canopy params
  List canopyParams = Rcpp::as<Rcpp::List>(x["canopy"]);
  
  //Root distribution input
  List below = Rcpp::as<Rcpp::List>(x["below"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(below["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(below["VGrhizo_kmax"]);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  
  //Transpiration params
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTransp"]);
  NumericVector leafWidth = paramsTransp["LeafWidth"];
  NumericVector Vmax298 = paramsTransp["Vmax298"];
  NumericVector Jmax298 = paramsTransp["Jmax298"];
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector pRootDisc = Rcpp::as<Rcpp::NumericVector>(paramsTransp["pRootDisc"]);
  NumericVector Gwmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gwmin"]);
  NumericVector Gwmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gwmax"]);
  
  
  //Leaf phenology and the adjusted leaf area index
  double canopyHeight = 0.0, LAIcell= 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(canopyHeight<H[c]) canopyHeight = H[c];
    LAIcell += (LAIphe[c]+LAIdead[c]);
  }
  int nz = ceil(canopyHeight/verticalLayerSize); //Number of vertical layers
  
  double latrad = latitude * (PI/180.0);
  if(NumericVector::is_na(aspect)) aspect = 0.0;
  double asprad = aspect * (PI/180.0);
  if(NumericVector::is_na(slope)) slope = 0.0;
  double slorad = slope * (PI/180.0);
  
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
  double tauday = meteoland::radiation_daylengthseconds(latrad,  slorad,asprad, delta); 
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  double Catm = 386.0;
  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);
  
  //Daily average water vapor pressure
  double vpa = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  
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
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpa, cloudcover);
  }
  Tcan[0] = canopyParams["Temp"]; //Take canopy temperature from previous day
  
  //Light extinction and absortion by time steps
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        kPAR, albedo,
                                                                        ddd, lwdr,
                                                                        ntimesteps,  canopyMode);
  List abs_PAR_SL_list = lightExtinctionAbsortion["PAR_SL"];
  List abs_PAR_SH_list = lightExtinctionAbsortion["PAR_SH"];
  List abs_SWR_SL_list = lightExtinctionAbsortion["SWR_SL"];
  List abs_SWR_SH_list = lightExtinctionAbsortion["SWR_SH"];
  List abs_LWR_SL_list = lightExtinctionAbsortion["LWR_SL"];
  List abs_LWR_SH_list = lightExtinctionAbsortion["LWR_SH"];
  NumericVector fsunlit = lightExtinctionAbsortion["fsunlit"];

  

  //Wind extinction profile
  NumericVector zWind;
  if(NumericVector::is_na(wind)) wind = 2.0; //Default wind speed when missing
  if(canopyMode=="multilayer") {
    zWind = windExtinctionProfile(zmid, wind, LAIcell, canopyHeight);
  } else if(canopyMode=="sunshade"){
    zWind = windExtinctionCohort(H,CR, wind,LAIcell, canopyHeight);
  }

  
  NumericVector Vc, VCroot_kmaxc, VGrhizo_kmaxc, psic, VG_nc,VG_alphac;
  
  List cohort_list(numCohorts);
  cohort_list.attr("names") = above.attr("row.names");
  
  for(int c=0;c<numCohorts;c++) {
    //Determine to which layers is plant connected
    LogicalVector layerConnected(nlayers, false);
    int nlayersc = 0;
    for(int l=0;l<nlayers;l++) {
      if(V(c,l)>0.0) {
        double pRoot = xylemConductance(psi[l], 1.0, VCroot_c[c], VCroot_d[c]); //Relative conductance in the root
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
        psic[cnt] = psi[l];
        VG_nc[cnt] = VG_n[l];
        VG_alphac[cnt] = VG_alpha[l];
        cnt++;
      }
    }
    
    double psiCav = 0.0;
    if(!cavitationRefill) {
      psiCav = xylemPsi(1.0-pEmb[c], 1.0, VCstem_c[c], VCstem_d[c]);//find water potential corresponding to this percentage of conductance loss
      // Rcout<< c <<" "<<psiCav<<"\n";
    }
    double minFlow = 1000.0*(Gwmin[c]*(tmin+tmax)/2.0)/Patm;
    List supplyNetwork = supplyFunctionNetwork(psic,
                                          VGrhizo_kmaxc,VG_nc,VG_alphac,
                                          VCroot_kmaxc, VCroot_c[c],VCroot_d[c],
                                          VCstem_kmax[c], VCstem_c[c],VCstem_d[c], 
                                          minFlow, psiCav, maxNsteps, psiStep, psiMax , ntrial, psiTol, ETol);
    
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
    double SLarea = 0.0, SHarea = 0.0;
    if(canopyMode=="sunshade"){
      for(int i=0;i<nz;i++) {
        SLarea +=SLarealayer[i];
        SHarea +=SHarealayer[i];
        Vmax298SL +=Vmax298layer[i]*LAIme(i,c)*fsunlit[i];
        Jmax298SL +=Jmax298layer[i]*LAIme(i,c)*fsunlit[i];
        Vmax298SH +=Vmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
        Jmax298SH +=Jmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
      }
    }
    List photo(ntimesteps);
    List PM(ntimesteps);
    for(int n=0;n<ntimesteps;n++) {
      
      //Long-wave radiation due to canopy temperature
      if(NumericVector::is_na(Tcan[n])) Tcan[n] = Tatm[n]; //If missing take above-canopy air temperature
      double lwcan = SIGMA_W*pow(Tcan[n]+273.16,4.0);
      
      // Rcout<< n<< " sh "<<solarHour[n]<< " t_sr " << t_sunrise<< " Tatm "<<Tatm<<"\n";
      if(canopyMode=="multilayer"){
        //Retrieve Light extinction
        NumericMatrix absPAR_SL = abs_PAR_SL_list[n];
        NumericMatrix absPAR_SH = abs_PAR_SH_list[n];
        NumericMatrix absSWR_SL = abs_SWR_SL_list[n];
        NumericMatrix absSWR_SH = abs_SWR_SH_list[n];
        NumericMatrix absLWR_SL = abs_LWR_SL_list[n];
        NumericMatrix absLWR_SH = abs_LWR_SH_list[n];
        for(int i=0;i<nz;i++) {
          QSL[i] = irradianceToPhotonFlux(absPAR_SL(i,c));
          QSH[i] = irradianceToPhotonFlux(absPAR_SH(i,c));
          //Add short wave and long wave radiation
          absRadSL[i] = absSWR_SL(i,c) + absLWR_SL(i,c) + SLarealayer[i]*0.97*lwcan;
          absRadSH[i] = absSWR_SH(i,c)+ absLWR_SH(i,c) + SHarealayer[i]*0.97*lwcan;
        }
        photo[n] = multilayerPhotosynthesisFunction(supplyNetwork, Catm, Patm,Tcan[n], vpa, 
                                                 SLarealayer, SHarealayer, 
                                                 zWind, absRadSL, absRadSH, 
                                                 QSL, QSH,
                                                 Vmax298layer,Jmax298layer,
                                                 Gwmin[c], Gwmax[c], leafWidth[c]);
      } else if(canopyMode=="sunshade"){
        //Retrieve Light extinction
        NumericVector absPAR_SL = abs_PAR_SL_list[n];
        NumericVector absPAR_SH = abs_PAR_SH_list[n];
        NumericVector absSWR_SL = abs_SWR_SL_list[n];
        NumericVector absSWR_SH = abs_SWR_SH_list[n];
        NumericVector absLWR_SL = abs_LWR_SL_list[n];
        NumericVector absLWR_SH = abs_LWR_SH_list[n];
        // Rcout<< "cohort " << c << " step "<< n << " "<<absSWR_SL[c]<<" " << absLWR_SL[c] << " "<< absSWR_SH[c] << " "<< absLWR_SH[c]<<"\n";
        photo[n] = sunshadePhotosynthesisFunction(supplyNetwork, Catm, Patm,Tcan[n], vpa, 
                                                  SLarea, SHarea, 
                                                  zWind[c], 
                                                  absSWR_SL[c] + absLWR_SL[c] + SLarea*0.97*lwcan, 
                                                  absSWR_SH[c] + absLWR_SH[c] + SHarea*0.97*lwcan, 
                                               irradianceToPhotonFlux(absPAR_SL[c]), irradianceToPhotonFlux(absPAR_SH[c]),
                                               Vmax298SL, Vmax298SH,
                                               Jmax298SL, Jmax298SH,
                                               Gwmin[c], Gwmax[c], leafWidth[c]);
        // Rcout<<n<<", "<< c <<": "<<lwcan<< " "<<absSWR_SL[c]<<" "<< absLWR_SL[c] << " "<<SLarea*0.97*lwcan<< " "<<absSWR_SH[c]<<" "<< absLWR_SH[c] << " "<<SHarea*0.97*lwcan<<"\n";
      }

      //Profit maximization
      PM[n] = profitMaximization(supplyNetwork, photo[n],  hydraulicCostFunction, VCstem_kmax[c]);
      

    }
    cohort_list[c] = List::create(_["supply"]=supplyNetwork,
                                  _["photo"]=photo,
                                  _["PM"] = PM);
  }
  return(cohort_list);
} 