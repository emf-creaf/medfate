// [[Rcpp::interfaces(r,cpp)]]
#include <cmath>
#include <vector>
#include <RcppArmadillo.h>
#include <meteoland.h>
#include "modelInput_c.h"
#include "communication_structures_c.h"
#include "soil_c.h"
#include "hydraulics_c.h"
#include "hydrology_c.h"
#include "biophysicsutils_c.h"
#include "numerical_solving_c.h"


/*=============================================================================
 * Implementation of hydrology routines using C++ code
 *=============================================================================*/


//' @rdname hydrology_soilEvaporation
//' 
//' @param DEF Water deficit in the (topsoil) layer.
//' @param PETs Potential evapotranspiration at the soil surface.
//' @param Gsoil Gamma parameter (maximum daily evaporation).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_soilEvaporationAmount")]]
double soilEvaporationAmount_c(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t + 1.0)-sqrt(t)), PETs);
  return(Esoil);
}

double soilEvaporation_c(Soil& soil,  
                         double snowpack, 
                         double pet, double LgroundSWR,
                         bool modifySoil = true) {
  
  
  double Esoil = 0.0;
  if(snowpack == 0.0) {
    double PETsoil = pet*(LgroundSWR/100.0);
    double Gsoil = 0.5; //TO DO, implement pedotransfer functions for Gsoil
    double W0 = soil.getW(0);
    double water_FC0 = soil.getWaterFC(0);
    double psi0 = soil.getPsi(0);
    // Allow evaporation only if water potential is higher than -2 MPa
    if(psi0 > -2.0) Esoil = soilEvaporationAmount_c((water_FC0*(1.0 - W0)), PETsoil, Gsoil);
    if(modifySoil){
      soil.setW(0, W0 - (Esoil/water_FC0));
    }
  }
  return(Esoil);
}


void herbaceousTranspiration_c(std::vector<double>& EherbVec, 
                               Soil& soil, 
                               double pet, double LherbSWR, 
                               double herbLAI,
                               const std::vector<double> V,
                               bool modifySoil = true){
  int nlayers = soil.getNlayers();
  for(int i=0;i<nlayers;i++) EherbVec[i] = 0.0;
  if(!std::isnan(herbLAI)) {
    double Tmax_herb = pet*(LherbSWR/100.0)*(0.134*herbLAI - 0.006*pow(herbLAI, 2.0));
    double psi0 = soil.getPsi(0);
    for(int l=0;l<nlayers;l++) {
      EherbVec[l] = V[l]*Tmax_herb*Psi2K_c(psi0, -1.5, 2.0); 
      if(modifySoil) {
        soil.setW(l, soil.getW(l) - (EherbVec[l]/soil.getWaterFC(l)));
      }
    }
  }
}

//Old defaults
//ERconv=0.05, ERsyn = 0.2
//New defaults
//Rconv = 5.6, Rsyn = 1.5

// Rainfall intensity calculation
double rainfallIntensity_c(int month, double prec, const std::vector<double>& rainfallIntensityPerMonth){
  double Ri_month = rainfallIntensityPerMonth[month - 1];
  double Ri = std::max(prec/24.0,Ri_month);
  return(Ri); 
}


// [[Rcpp::export(".hydrology_interceptionGashDay")]]
double interceptionGashDay_c(double Rainfall, double Cm, double p, double ER=0.05) {
  double I = 0.0;
  double PG = (-Cm/(ER*(1.0-p)))*log(1.0-ER); //Rainfall need to saturate the canopy
  if(Cm==0.0 || p==1.0) PG = 0.0; //Avoid NAs
  if(Rainfall>PG) {
    I = (1-p)*PG + (1-p)*ER*(Rainfall-PG);
  } else {
    I = (1-p)*Rainfall;
  }
  return(I);
}

// [[Rcpp::export(".hydrology_interceptionLiuDay")]]
double interceptionLiuDay_c(double Rainfall, double Cm, double p, double ER=0.05){
  double I = Cm*(1.0 - exp(-1.0*(Rainfall)*((1.0 - p)/Cm)))*(1.0 - (ER/(1.0 - p))) + (ER*Rainfall);
  return(I);
}

void infiltrationRepartition_c(double I, 
                               std::vector<double> &Ivec, 
                               const std::vector<double> &widths, 
                               const std::vector<double> &macro, 
                               double a = -0.005, double b = 3.0) {
  int nlayers = widths.size();
  double pi = 0.0;
  double z1 = 0.0;
  double p1 = 1.0;
  for(int i=0;i<nlayers;i++) {
    double ai = a*pow(1.0-macro[i],b);
    if(i<(nlayers-1)) {
      pi = p1*(1.0-exp(ai*widths[i]));
    } else {
      pi = p1;
    }
    p1 = p1*exp(ai*widths[i]);
    z1 = z1 + widths[i];
    Ivec[i] = I*pi;
  }
}


//' @rdname hydrology_infiltration
//' 
//' @param input A numeric vector of (daily) water input (in mm of water).
//' @param Ssoil Soil water storage capacity (can be referred to topsoil) (in mm of water).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationBoughton")]]
double infiltrationBoughton_c(double input, double Ssoil) {
  double I = 0;
  if(input>0.2*Ssoil) {
    I = input-(pow(input-0.2*Ssoil,2.0)/(input+0.8*Ssoil));
  } else {
    I = input;
  }
  return(I);
}

double fGreenAmpt_c(double x, double t, double psi_w, double Ksat, double delta_theta) {
  double f = Ksat*t + std::abs(psi_w)*delta_theta*log(1.0 + (x/(std::abs(psi_w)*delta_theta))) - x;
  return(f);
}
double fGreenAmptDer_c(double x, double t, double psi_w, double Ksat, double delta_theta) {
  double fder = (log(1.0 + (x/(std::abs(psi_w)*delta_theta)))/(1.0+ (x/(std::abs(psi_w)*delta_theta))))-1.0;
  return(fder);
}

//' @rdname hydrology_infiltration
//' 
//' @param t Time of the infiltration event
//' @param psi_w Matric potential at the wetting front
//' @param Ksat hydraulic conductivity at saturation
//' @param theta_sat volumetric content at saturation
//' @param theta_dry volumetric content at the dry side of the wetting front
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationGreenAmpt")]]
double infitrationGreenAmpt_c(double t, double psi_w, double Ksat, double theta_sat, double theta_dry) {
   double delta_theta = theta_sat - theta_dry;
   double x,x1,e,fx,fx1;
   x1 = 0.0;//initial guess
   e = 0.001; // accuracy in mm
   int cnt = 0;
   int mxiter = 100;
   do {
     x=x1; /*make x equal to the last calculated value of  x1*/
     fx=fGreenAmpt_c(x, t, psi_w, Ksat, delta_theta);            //simplifying f(x)to fx
     fx1=fGreenAmptDer_c(x, t, psi_w, Ksat, delta_theta);            //simplifying fprime(x) to fx1
     x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
     cnt++;
   } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
   return(x);
 }


double infiltrationAmount_c(double rainfallInput, double rainfallIntensity, Soil& soil, 
                            std::string model = "GreenAmpt1911", double K_correction = 1.0) {
  double infiltration = 0.0;
  if(model=="GreenAmpt1911") {
    double theta_dry0 = soil.getTheta(0);
    ClappHornberger cp = soil.getClappHornberger();
    double t = std::min(24.0, rainfallInput/rainfallIntensity); // time in hours
    double psi_w = cp.psi_sat_cm*((2.0*cp.b + 3.0)/(2*cp.b + 6.0));
    double K_sat_0 = K_correction*soil.getKsat(0)/(24.0*cmdTOmmolm2sMPa); // from mmolH20*m-2*MPa-1*s-1 to cm_h
    infiltration = infitrationGreenAmpt_c(t, psi_w, K_sat_0, cp.theta_sat, theta_dry0);
  } else if(model=="Boughton1989") {
    infiltration = infiltrationBoughton_c(rainfallInput, soil.getWaterFC(0));
  }
  infiltration = std::min(infiltration, rainfallInput);
  return(infiltration);
}

//' @rdname hydrology_verticalInputs
//' 
//' @param tday Average day temperature (ºC).
//' @param rad Solar radiation (in MJ/m2/day).
//' @param elevation Altitude above sea level (m).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_snowMelt")]]
double snowMelt_c(double tday, double rad, double LgroundSWR, double elevation) {
  //missing data checks
  if(std::isnan(rad)) Rcpp::stop("Missing radiation data for snow melt!\n");
  if(std::isnan(elevation)) Rcpp::stop("Missing elevation data for snow melt!\n");
  double rho = meteoland::utils_airDensity(tday, meteoland::utils_atmosphericPressure(elevation));
  double ten = (86400.0*tday*rho*1013.86*1e-6/100.0); //ten can be negative if temperature is below zero
  double ren = (rad*(LgroundSWR/100.0))*(0.1); //90% albedo of snow
  double melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
  return(melt);
}


void waterInputs_c(std::vector<double>& waterInputs,
                   ModelInput& x,
                   double prec, double rainfallIntensity,
                   double pet, double tday, double rad, double elevation,
                   double Cm, double LgroundPAR, double LgroundSWR, 
                   bool modifyInput = true) {
  
  //Soil input
  std::string soilFunctions = x.control.commonWB.soilFunctions; 
  std::string interceptionMode = x.control.commonWB.interceptionMode;
  
  double swe = x.snowpack; //snow pack
  double er = pet/(24.0*rainfallIntensity);
  
  //Snow pack dynamics
  double snow = 0.0, rain=0.0;
  double melt = 0.0;
  //Turn rain into snow and add it into the snow pack
  if(tday < 0.0) { 
    snow = prec; 
    swe = swe + snow;
  } else {
    rain = prec;
  }
  //Apply snow melting
  if(swe > 0.0) {
    melt = std::min(swe, snowMelt_c(tday, rad, LgroundSWR, elevation));
    // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
    swe = swe-melt;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Interception = 0.0;
  if(rain>0.0)  {
    if(interceptionMode=="Gash1995") {
      Interception = interceptionGashDay_c(rain,Cm,LgroundPAR/100.0,er);
    } else if(interceptionMode =="Liu2001") {
      Interception = interceptionLiuDay_c(rain,Cm,LgroundPAR/100.0,er);
    } else {
      throw medfate::MedfateInternalError("Wrong interception model!");
    }
    NetRain = rain - Interception; 
  }
  waterInputs[0] = rain;
  waterInputs[1] = snow;
  waterInputs[2] = Interception;
  waterInputs[3] = NetRain;
  waterInputs[4] = melt;
  if(modifyInput) {
    x.snowpack = swe;
  }
}

//Imbitition from macropores to micropores following Larsbo et al. 2005, eq. 6-7
//Diffusion equation with gradients in water content as driving force
//Returns m3/m3/s
double microporeImbibitionRate_c(double theta_b, double theta_micro, 
                                 double D_theta_b, double D_theta_micro,
                                 double S_macro) {
  double G_f = 3.0;//geometry factor
  double gamma_w = 0.1; //Scaling factor
  double D_w = ((D_theta_b + D_theta_micro)/2.0)*S_macro;//Effective water diffusivity
  double d = 9.35*1e-3; //Effective diffusion pathlength in m
  double S_w = std::max(0.0, ((G_f*D_w*gamma_w)/(d*d))*(theta_b - theta_micro));
  // Rcout<< "D_w " << D_w << " S_w "<< S_w<<"\n";
  return(S_w);
}


double rootFindingMacropores_c(double S_t, double K_up, double Ksat_ms, double Ksat_b_ms, double kin_exp,
                               double e_macro, double lambda, double dZ_m, double sourceSink_macro_m3s, double tstep, 
                               int Nmax = 100) {
  double a = 0.0;
  //function at 'a'
  double Ka = (Ksat_ms - Ksat_b_ms)*pow(a, kin_exp);
  double f_a = S_t - a + (tstep/(e_macro*lambda*dZ_m))*((K_up - Ka) + sourceSink_macro_m3s);
  double b = a + 1.0;
  //function at 'b'
  double Kb = (Ksat_ms - Ksat_b_ms)*pow(b, kin_exp);
  double f_b = S_t - b + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kb) + sourceSink_macro_m3s);
  while(f_b > 0.0) {
    b = b + 1.0;
    Kb = (Ksat_ms - Ksat_b_ms)*pow(b, kin_exp);
    f_b = S_t - b + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kb) + sourceSink_macro_m3s);
    if(b>10.0) Rcpp::stop("Could not find appropriate bounds for macropore circulation");
  }
  // Rcout<< "Ka "<< Ka <<" f_a "<< f_a<<" Kb "<< Kb <<" f_b "<< f_b<<"\n";
  
  int N = 1;
  double c = NA_REAL;
  double Kc = NA_REAL;
  double f_c = NA_REAL;
  double tol = 1e-7;
  bool has_to_stop = false;
  bool found = false;
  while(!has_to_stop) {
    c = (a + b)/2.0; // new midpoint
    //function at 'c'
    Kc = (Ksat_ms - Ksat_b_ms)*pow(c, kin_exp);
    f_c = S_t - c + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kc) + sourceSink_macro_m3s);
    // Rcout<<N << " std::abs((b - a)/2.0 "<< std::abs((b - a)/2.0) << " c "<< c <<" f_c "<< f_c<<"\n";
    if((f_c == 0) || (std::abs((b - a)/2.0) < tol)) { // solution found
      found = true;
      has_to_stop = true;
    }
    if(((f_c > 0) && (f_a > 0)) || ((f_c < 0) && (f_a < 0)) ) { //new interval
      a = c;
      f_a = f_c;
    } else {
      b = c;
      f_b = f_c;
    }
    N = N + 1; //Increment step counter
    if(N == Nmax) {
      has_to_stop = true;
    }
  }
  if(!found) {
    Rcpp::stop("Not found");
  }
  return(c);
}


void soilWaterBalance_inner_c(SoilWaterBalance_RESULT &SWBres, SoilWaterBalance_COMM &SWBcomm, Soil &soil, 
                              double rainfallInput, double rainfallIntensity, double snowmelt, const std::vector<double> &sourceSink, 
                              double runon, const std::vector<double> &lateralFlows, double waterTableDepth,
                              std::string infiltrationMode = "GreenAmpt1911", double infiltrationCorrection = 5.0, 
                              std::string soilDomains = "buckets", 
                              int nsteps = 24, int max_nsubsteps = 3600) {
  
  if((soilDomains!="single") && (soilDomains!="dual")  && (soilDomains!="buckets") ) throw medfate::MedfateInternalError("Unrecognized soilDomain value");

  const std::vector<double> widths = soil.getWidths();
  const std::vector<double> macro = soil.getMacro();
  const std::vector<double> rfc = soil.getRFC();
  
  int nlayers = soil.getNlayers();
  double soilDepth = std::accumulate(widths.begin(), widths.end(), 0.0);
  
  bool micropore_imbibition = true;

  double mm_2_m3 = 0.001;
  double m3_2_mm = 1.0/mm_2_m3;
  double mm_day_2_m3_s = mm_2_m3*(1.0/86400.0);//From mm/day = l/m2/day = dm3/day to m3/m2/s

  //Infiltration
  double K_correction = 1.0;
  if(soilDomains!="dual") K_correction = infiltrationCorrection;
  double infiltration_matrix_mm = infiltrationAmount_c(rainfallInput, rainfallIntensity, soil,
                                                       infiltrationMode, K_correction);
  double infiltration_macropores_mm = 0.0;
  double infiltration_excess_matrix_mm = rainfallInput - infiltration_matrix_mm;
  double infiltration_target_macropores_mm = infiltration_excess_matrix_mm;
  double infiltration_excess_macropores_mm = 0.0;
  double matrix_correction_mm = 0.0;
  double macropore_correction_mm = 0.0;

  
  //Add snow-melt and runon to infiltration_matrix_mm
  double snowmelt_mm = 0.0;
  if(!std::isnan(snowmelt)) {
    snowmelt_mm = snowmelt;
    infiltration_matrix_mm += snowmelt_mm;
  }
  double runon_mm = 0.0;
  if(!std::isnan(runon)) {
    runon_mm = runon;
    infiltration_matrix_mm += runon_mm;
  }

  //Copy sinks
  std::vector<double> source_sink_def_mm(nlayers, 0.0);
  for(int l=0;l<nlayers;l++) {
    source_sink_def_mm[l] = sourceSink[l];
  }

  //Add infiltration to matrix def source/sinks
  if(soilDomains!="dual") {
    std::vector<double> IVec(nlayers, 0.0); 
    infiltrationRepartition_c(infiltration_matrix_mm,
                              IVec,
                              widths, macro);
    for(int l=0;l<nlayers;l++) {
      source_sink_def_mm[l] += IVec[l];
    }
  } else {
    source_sink_def_mm[0] += infiltration_matrix_mm;
  }
  //Add lateral flows to matrix def source/sinks
  for(int l=0;l<nlayers;l++) {
    if(!std::isnan(lateralFlows[l])) source_sink_def_mm[l] += lateralFlows[l];
  }
  
  //Main soil water balance calculation
  if(soilDomains == "buckets") {
    const std::vector<double> Ksat = soil.getKsat();
    
    std::vector<double> Water_FC(nlayers);
    std::vector<double> Water_SAT(nlayers);
    std::vector<double> W(nlayers);
    for(int l=0;l<nlayers;l++) {
      Water_FC[l] = soil.getWaterFC(l); //Field capacity in mm
      Water_SAT[l] = soil.getWaterSAT(l); //Saturation in mm
      W[l] = soil.getW(l); //Current water content (relative to field capacity)
    }
    double drainage_matrix_mm = 0.0;
    double saturation_excess_matrix_mm = 0.0;
    double Wn;
    for(int l=0;l<nlayers;l++) {
      if((widths[l]>0.0)) {
        Wn = W[l]*Water_FC[l] + source_sink_def_mm[l]; //Update water volume
        W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta
        // Rcout<< source_sink_def_mm[l]<< " " << W[l] <<"\n";
        if(l<(nlayers-1)) {
          //update source_sink_def_mm adding the excess to the infiltrating water (saturated flow)
          source_sink_def_mm[l+1] = source_sink_def_mm[l+1] + std::max(Wn - Water_FC[l],0.0);
        } else {
          saturation_excess_matrix_mm = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer using field capacity
          // Rcout << drainage_matrix_mm << "\n";
        }
      }
    }
    //If there still excess fill layers over field capacity
    if((saturation_excess_matrix_mm>0.0)) {
      for(int l=(nlayers-1);l>=0;l--) {
        if((widths[l]>0.0) && (saturation_excess_matrix_mm>0.0)) {
          Wn = W[l]*Water_FC[l] + saturation_excess_matrix_mm; //Update water volume
          saturation_excess_matrix_mm = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
          W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta here no upper
        }
      }
    }
    //If there is still room for additional drainage (water in macropores accumulated from previous days)
    double head = 0.0;
    for(int l=0;l<nlayers;l++) { //Add mm over field capacity
      head += Water_FC[l]*std::max(W[l] - 1.0, 0.0);
    }
    if(head>0.0) {
      // Saturated vertical hydraulic conductivity (mm/day)
      double cmdTOmmolm2sMPa = 655.2934;
      double Kdrain = 10.0*(Ksat[nlayers-1]*(1.0 - (rfc[nlayers-1]/100.0))/cmdTOmmolm2sMPa);
      double maxDrainage = Kdrain;
      // Rcout<<head<< " "<< maxDrainage <<"\n";
      for(int l=0;l<nlayers;l++) {
        if(maxDrainage>0.0) {
          double Wn = W[l]*Water_FC[l];
          double toDrain = std::min(std::max(Wn - Water_FC[l], 0.0), maxDrainage);
          if(toDrain > 0.0) {
            drainage_matrix_mm +=toDrain;
            maxDrainage -=toDrain;
            Wn -= toDrain;
            W[l] = std::max(0.0, Wn/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
          }
        }
      }
    }
    for(int l=0;l<nlayers;l++) {
      soil.setW(l, W[l]);
    }
    SWBres.localSourceSinks_mm = std::accumulate(sourceSink.begin(), sourceSink.end(), 0.0);
    SWBres.lateralSourceSinks_mm = std::accumulate(lateralFlows.begin(), lateralFlows.end(), 0.0);
    SWBres.infiltration_mm = infiltration_matrix_mm;
    SWBres.infiltrationExcess_mm = infiltration_excess_matrix_mm;
    SWBres.saturationExcess_mm = saturation_excess_matrix_mm;
    SWBres.runoff_mm = infiltration_excess_matrix_mm + saturation_excess_matrix_mm;
    SWBres.deepDrainage_mm = drainage_matrix_mm;
    SWBres.capillarityRise_mm = 0.0;
  } else if(soilDomains == "dual" || soilDomains == "single") {
    //Initialize matrix-macrore flows (positive in the direction of matrix)
    std::vector<double> matrix_macropore_flows_mm(nlayers, 0.0);

    //Set time steps
    double tstep = 86400.0/((double) nsteps);
    double tsubstep = tstep;
    double halftsubstep = tsubstep/2.0;
    double rainfallIntensity_step = rainfallIntensity*24.0/((double) nsteps); //mm/step

    std::vector<double> source_sink_def_m3s(nlayers, 0.0);
    for(int l=0;l<nlayers;l++) source_sink_def_m3s[l] =  source_sink_def_mm[l]*mm_day_2_m3_s;
    std::vector<double> matrixImbibition_m3s(nlayers, 0.0);
    std::vector<double> matrixExcess_m3s(nlayers, 0.0);
    std::vector<double> saturated_matrix_correction_m3s(nlayers, 0.0);
    std::vector<double> saturated_macropore_correction_m3s(nlayers, 0.0);

    std::vector<double> dZ_m = SWBcomm.dZ_m;
    for(int l=0;l<nlayers;l++) dZ_m[l] = widths[l]*0.001; //mm to m
    std::vector<double>  dZUp = SWBcomm.dZUp;
    std::vector<double> dZDown = SWBcomm.dZDown;
    std::vector<double> lambda = SWBcomm.lambda;

    //Estimate layer interfaces
    for(int l=0;l<nlayers;l++) {
      lambda[l] = 1.0 - (rfc[l]/100.0);
      if(l==0) { //first layer
        dZUp[l] = dZ_m[0]/2.0;
      } else {
        dZUp[l] = (dZ_m[l - 1]/2.0) + (dZ_m[l]/2.0);
      }
      if(l<(nlayers - 1)) {
        dZDown[l] = (dZ_m[l]/2.0) + (dZ_m[l + 1]/2.0);
      } else { //last layer
        dZDown[l] = dZ_m[l]/2.0;
      }
    }

    std::vector<double> prop_saturated(nlayers, 0.0);
    int num_saturated = 0;
    double freeDrainage = true;
    if(!std::isnan(waterTableDepth)) {
      freeDrainage = (waterTableDepth > soilDepth);
      double sZ = 0.0;
      for(int i=0;i<nlayers;i++){
        prop_saturated[i] = std::min(1.0, std::max(0.0,((sZ + widths[i]) - waterTableDepth)/widths[i]));
        if(prop_saturated[i]==1.0) num_saturated++;
        sZ += widths[i];
        // Rcout << i <<" "<< prop_saturated[i] << "\n";
      }
    }
    // Rcout << " num_saturated " << num_saturated << "\n";

    //boundary condition of water table (if freeDrainage = FALSE)
    double Psi_bc = 0.0;
    double Psi_quasi_sat = -0.0000001;
    //Boundary water potential for dual porosity
    double Psi_b = -0.1*mTOMPa; // 10 cm = 0.1 m
    double kin_exp = 2.23; //Kinematic exponent reflecting macropore size distribution and tortuosity

    //Retrieve VG parameters
    const std::vector<double>  Ksat_ori = soil.getKsat();
    const std::vector<double>  n = soil.getVG_n();
    const std::vector<double>  alpha = soil.getVG_alpha();
    const std::vector<double>  theta_res = soil.getVG_theta_res();
    const std::vector<double>  theta_sat = soil.getVG_theta_sat();
    
    //Get initial theta and psi state
    std::vector<double>  Theta(nlayers, 0.0);
    for(int l=0;l<nlayers;l++) {
      Theta[l] = soil.getTheta(l);
    }
    
    //Estimate Psi, C, K

    //Microporosity or single domain
    std::vector<double> theta_micro = SWBcomm.theta_micro;
    std::vector<double> theta_b = SWBcomm.theta_b;
    std::vector<double> theta_macro = SWBcomm.theta_macro;
    std::vector<double> theta_sat_fict = SWBcomm.theta_sat_fict;
    std::vector<double> Ksat_b = SWBcomm.Ksat_b;
    std::vector<double> Ksat_b_ms = SWBcomm.Ksat_b_ms;
    std::vector<double> Ksat = SWBcomm.Ksat;
    std::vector<double> Ksat_ms = SWBcomm.Ksat_ms;
    std::vector<double> Psi = SWBcomm.Psi;
    std::vector<double> K = SWBcomm.K;
    std::vector<double> C = SWBcomm.C;
    std::vector<double> Psi_m = SWBcomm.Psi_m;
    std::vector<double> K_ms = SWBcomm.K_ms;
    std::vector<double> Kbc = SWBcomm.Kbc;
    std::vector<double> Kbc_ms = SWBcomm.Kbc_ms;
    std::vector<double> C_m = SWBcomm.C_m;
    for(int l=0;l<nlayers;l++) {
      theta_micro[l] = 0.0;
      theta_b[l] = 0.0;
      theta_macro[l] = 0.0;
      theta_sat_fict[l] = 0.0;
      Ksat_b[l] = 0.0;
      Ksat_b_ms[l] = 0.0;
      Ksat[l] = 0.0;
      Psi[l] = 0.0;
      K[l] = 0.0;
      C[l] = 0.0;
      Psi_m[l] = 0.0;
      K_ms[l] = 0.0;
      Kbc[l] = 0.0;
      Kbc_ms[l] = 0.0;
      C_m[l] = 0.0;
    }

    //Macroporosity domain
    std::vector<double> S_macro, e_macro, Kmacro_ms;
    if(soilDomains=="dual") {
      S_macro = SWBcomm.S_macro;
      e_macro = SWBcomm.e_macro;
      Kmacro_ms = SWBcomm.Kmacro_ms;
      for(int l=0;l<nlayers;l++) {
        S_macro[l] = 0.0;
        e_macro[l] = 0.0;
        Kmacro_ms[l] = 0.0;
      }
    }

    std::vector<double> waterFluidity = SWBcomm.waterFluidity;
    for(int l=0;l<nlayers;l++) {
      waterFluidity[l] = 1.0;
      double temp = soil.getTemp(l);
      if(!std::isnan(temp)) {
        if(temp>0) {
          waterFluidity[l] = 1.0/waterDynamicViscosity_c(temp);
        } else {
          waterFluidity[l] = 0.0;
        }
      }
      Ksat[l] = Ksat_ori[l]*lambda[l];//Multiply K for the space available for water movement
      Ksat_ms[l] = 0.01*waterFluidity[l]*Ksat[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      if(soilDomains=="single") {
        Psi[l] = theta2psiVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]);
        C[l] = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
        K[l] = waterFluidity[l]*psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
        Kbc[l] = waterFluidity[l]*psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_bc);
      } else {
        theta_sat_fict[l] = theta_sat[l] - macro[l];
        //Matching theta point and theta partitioning
        //This ensures theta_b < theta_sat - macro & e_macro > macro
        theta_b[l] = psi2thetaVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat_fict[l], Psi_b);
        e_macro[l] = theta_sat[l] - theta_b[l];
        // Rcout<<e_macro[l]<< " "<< macro[l]<<"\n";
        Ksat_b[l] = psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_b);
        theta_micro[l] = std::min(Theta[l], theta_b[l]);
        theta_macro[l] = std::max(Theta[l] - theta_b[l], 0.0);
        //Water potential, conductivity and capacitance in the micropore domain according to the modified retention
        Psi[l] = theta2psiVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat_fict[l], theta_micro[l]);
        C[l] = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi[l]);
        K[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores_c(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                              Psi[l], Psi_b);
        Kbc[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores_c(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                                Psi_bc, Psi_b);
        //Effective saturation and conductivity in the macropore domain
        S_macro[l] = theta_macro[l]/e_macro[l];
        Ksat_b_ms[l] = 0.01*waterFluidity[l]*Ksat_b[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
        Kmacro_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_macro[l], kin_exp);
      }
      Psi_m[l]= Psi[l]/mTOMPa; // MPa to m
      K_ms[l] = 0.01*K[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      Kbc_ms[l] = 0.01*Kbc[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      C_m[l] = C[l]*mTOMPa; //From MPa-1 to m-1
    }
    //Initial soil volume
    double Vini_mm = 0.0;
    double Vini_micro_mm = 0.0;
    double Vini_macro_mm = 0.0;
    double Vini_step_mm = 0.0;
    double Vini_step_micro_mm = 0.0;
    double Vini_step_macro_mm = 0.0;
    double Vfin_mm = 0.0;
    double Vfin_micro_mm = 0.0;
    double Vfin_macro_mm = 0.0;

    if(soilDomains=="dual"){
      for(int l=0;l<nlayers;l++) {
        Vini_micro_mm += theta_micro[l]*widths[l]*lambda[l];
        Vini_macro_mm += theta_macro[l]*widths[l]*lambda[l];
      }
      Vini_step_micro_mm = Vini_micro_mm;
      Vini_step_macro_mm = Vini_macro_mm;
    } else {
      for(int l=0;l<nlayers;l++) {
        Vini_mm += Theta[l]*widths[l]*lambda[l];
      }
      Vini_step_mm = Vini_mm;
    }
    double Vini0_mm = Vini_mm;
    double Vini0_macro_mm = Vini_macro_mm;
    double Vini0_micro_mm = Vini_micro_mm;

    //Temporary step variables
    double drainage_matrix_m3 = 0.0, drainage_macropores_m3 = 0.0;
    double capillarity_matrix_m3 = 0.0;
    double capillarity_macropores_m3 = 0.0;
    double saturation_excess_matrix_mm = 0.0;
    double saturation_excess_macropores_mm = 0.0;
    double K_up= 0.0, K_down= 0.0;

    double C_step_05, K_step_05;

    std::vector<double> a = SWBcomm.a;
    std::vector<double> b = SWBcomm.b;
    std::vector<double> c = SWBcomm.c;
    std::vector<double> d = SWBcomm.d;
    std::vector<double> e = SWBcomm.e;
    std::vector<double> f = SWBcomm.f;
    std::vector<double> Psi_step_t1(nlayers), Psi_step_t05(nlayers);

    std::vector<double> K_step_ms05 = SWBcomm.K_step_ms05;
    std::vector<double> C_step_m05 = SWBcomm.C_step_m05;
    std::vector<double> C_step = SWBcomm.C_step;
    std::vector<double> C_step_m = SWBcomm.C_step_m;
    std::vector<double> K_step_ms = SWBcomm.K_step_ms;
    std::vector<double> K_step = SWBcomm.K_step;
    std::vector<double> Psi_step = SWBcomm.Psi_step;
    std::vector<double> Psi_step_m = SWBcomm.Psi_step_m;
    std::vector<double> S_macro_step = SWBcomm.S_macro_step;
    std::vector<double> Kmacro_step_ms = SWBcomm.Kmacro_step_ms;
    std::vector<double> theta_macro_step = SWBcomm.theta_macro_step;
    std::vector<double> theta_micro_step = SWBcomm.theta_micro_step;
    for(int l=0;l<nlayers;l++) {
      K_step_ms05[l] = 0.0;
      C_step_m05[l] = 0.0;
      C_step[l] = 0.0;
      C_step_m[l] = 0.0;
      K_step_ms[l] = 0.0;
      K_step[l] = 0.0;
      Psi_step[l] = 0.0;
      Psi_step_m[l] = 0.0;
      S_macro_step[l] = 0.0;
      Kmacro_step_ms[l] = 0.0;
      theta_macro_step[l] = 0.0;
      theta_micro_step[l] = 0.0;
    }
    std::vector<double> finalSourceSinks_m3s = SWBcomm.finalSourceSinks_m3s;
    std::vector<double> capill_below = SWBcomm.capill_below;
    std::vector<double> drain_above = SWBcomm.drain_above;
    std::vector<double> drain_below = SWBcomm.drain_below;
    std::vector<double> lateral_flows_step_mm = SWBcomm.lateral_flows_step_mm;
    for(int l=0;l<nlayers;l++) {
      finalSourceSinks_m3s[l] = 0.0;
      capill_below[l] = 0.0;
      drain_above[l] = 0.0;
      drain_below[l] = 0.0;
      lateral_flows_step_mm[l] = 0.0;
    }

    double drainage_matrix_step_m3 = 0.0;
    double drainage_macropores_step_m3 = 0.0;
    double capillarity_matrix_step_m3 = 0.0;
    double capillarity_macropores_step_m3 = 0.0;
    double saturation_excess_matrix_step_mm = 0.0;
    double saturation_excess_macropores_step_mm = 0.0;
    double infiltration_excess_macropores_step_mm = 0.0;
    double infiltration_macropores_step_mm = 0.0;
    double matrix_correction_step_mm = 0.0;
    double macropore_correction_step_mm = 0.0;
    double infiltration_remaining_macropores_step_mm = 0.0;
    double infiltration_target_macropores_step_mm = 0.0;

    int total_nsubsteps = 0;
    int nsubsteps = 1;
    bool cont = true;

    for(int s =0;s<nsteps;s++) {

      //Target infiltration for macropores
      infiltration_target_macropores_step_mm = std::min(infiltration_target_macropores_mm, rainfallIntensity_step);
      infiltration_target_macropores_mm -= infiltration_target_macropores_step_mm;

      nsubsteps = 1;
      cont = true;

      while(cont) {
        saturation_excess_macropores_step_mm = 0.0;
        infiltration_excess_macropores_step_mm = 0.0;
        saturation_excess_matrix_step_mm = 0.0;
        drainage_matrix_step_m3 = 0.0;
        capillarity_matrix_step_m3 = 0.0;
        capillarity_macropores_step_m3 = 0.0;
        drainage_macropores_step_m3 = 0.0;
        infiltration_macropores_step_mm = 0.0;
        for(int l=0;l<nlayers;l++) {
          //Reset comunication variables
          lateral_flows_step_mm[l] = 0.0;
          matrixImbibition_m3s[l] = 0.0;
          matrixExcess_m3s[l] = 0.0;
          //Copy to step variables
          C_step[l] = C[l];
          C_step_m[l] = C_m[l];
          K_step[l] = K[l];
          K_step_ms[l] = K_ms[l];
          Psi_step[l] = Psi[l];
          Psi_step_m[l] = Psi_m[l];
          if(soilDomains=="dual") {
            S_macro_step[l] = S_macro[l];
            Kmacro_step_ms[l] = Kmacro_ms[l];
            theta_macro_step[l] = theta_macro[l];
          }
          theta_micro_step[l] = theta_micro[l];
        }

        tsubstep = tstep/((double) nsubsteps);
        halftsubstep = tsubstep/2.0;

        infiltration_remaining_macropores_step_mm = infiltration_target_macropores_step_mm;


        for(int ss=0;ss<nsubsteps;ss++) {
          total_nsubsteps++;
          //Correction for saturation
          for(int l=0;l<nlayers;l++) {
            saturated_matrix_correction_m3s[l] = 0.0;
            saturated_macropore_correction_m3s[l] = 0.0;
            if(prop_saturated[l]>0.0) {
              double theta_l= 0.0, quasi_sat_theta_l = 0.0;
              if(soilDomains =="single") {
                theta_l = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
                quasi_sat_theta_l = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_quasi_sat);
                saturated_matrix_correction_m3s[l] = std::max(0.0, ((prop_saturated[l]*quasi_sat_theta_l - theta_l)*lambda[l]*widths[l]*0.001)/tsubstep);
              } else {
                theta_l = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
                quasi_sat_theta_l = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_quasi_sat);
                saturated_matrix_correction_m3s[l] = std::max(0.0, ((prop_saturated[l]*theta_b[l] - theta_l)*lambda[l]*widths[l]*0.001)/tsubstep);
                saturated_macropore_correction_m3s[l] = std::max(0.0, ((1.0 - S_macro_step[l])*(prop_saturated[l]*quasi_sat_theta_l - theta_b[l])*lambda[l]*widths[l]*0.001)/tsubstep);
              }
              saturated_matrix_correction_m3s[l] = std::min(Ksat_ms[l], saturated_matrix_correction_m3s[l]);
              saturated_macropore_correction_m3s[l] = std::min(Ksat_ms[l], saturated_macropore_correction_m3s[l]);
              // Rcout<< s << " "<< nsubsteps<< " " << l << " theta " << theta << " quasi_sat_theta " << quasi_sat_theta <<  " Ksat_ms "<< Ksat_ms[l]<< " sat micro corr: " << saturated_matrix_correction_m3s[l]<< " sat macro corr: " << saturated_macropore_correction_m3s[l] <<"\n";
              saturated_matrix_correction_m3s[l] = std::max(0.0, saturated_matrix_correction_m3s[l] - source_sink_def_m3s[l] - matrixImbibition_m3s[l]);
              capillarity_matrix_step_m3 += tsubstep*saturated_matrix_correction_m3s[l];
            }
          }
          if(soilDomains=="dual"){
            if(micropore_imbibition) {
              for(int l=0;l<nlayers;l++) {
                //Update imbibition rate (m3s)
                double Ksat_fict_ms = psi2kVanGenuchtenMicropores_c(Ksat_b_ms[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                                    0.0, Psi_b);
                double D_theta_b_m2s = psi2DVanGenuchten_c(Ksat_fict_ms, n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                           Psi_b);
                double D_theta_micro_m2s = psi2DVanGenuchten_c(Ksat_fict_ms, n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                               Psi_step[l]);
                //Calculate imbibition rate in m3/m3/s = s-1
                double imbibitionRate = microporeImbibitionRate_c(theta_b[l], theta_micro_step[l],
                                                                  D_theta_b_m2s, D_theta_micro_m2s,
                                                                  S_macro_step[l]);
                matrixImbibition_m3s[l] = dZ_m[l]*lambda[l]*imbibitionRate;
                lateral_flows_step_mm[l] += imbibitionRate*widths[l]*lambda[l]*tsubstep; //From m3/m3/s to mm/step
              }
            }
          }

          //Psi-based solution of the Richards equation using implicit solution for psi
          //but with explicit linearization for K and C (pp. 126, Bonan)

          //0. Sum source/sinks
          for(int l=0;l<nlayers;l++) finalSourceSinks_m3s[l] = source_sink_def_m3s[l] + matrixImbibition_m3s[l] + saturated_matrix_correction_m3s[l];

          //A. Predictor sub-step
          for(int l=0;l<nlayers;l++) {

            if(l==0) { //first layer
              K_up = 0.0;
              // K_down = 0.5*(K_step_ms[l] + K_step_ms[l+1]);
              K_down = K_step_ms[l];
              drain_below[l] = K_down;
              if(prop_saturated[l]==1.0) {
                drain_below[l] = 0.0;
              }
              capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_step_m[l+1]);
              a[l] = 0.0;
              c[l] = -1.0*K_down/dZDown[l];
              b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - c[l];
              d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l] - drain_below[l] + finalSourceSinks_m3s[l];
            } else if(l<(nlayers - 1)) {
              // K_up = 0.5*(K_step_ms[l-1] + K_step_ms[l]);
              // K_down = 0.5*(K_step_ms[l] + K_step_ms[l+1]);
              K_up = K_step_ms[l-1];
              K_down = K_step_ms[l];
              drain_above[l] = K_up;
              drain_below[l] = K_down;
              if(prop_saturated[l]==1.0) {
                drain_above[l] = 0.0;
                drain_below[l] = 0.0;
              }
              capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_step_m[l+1]);
              a[l] = -1.0*K_up/dZUp[l];
              c[l] = -1.0*K_down/dZDown[l];
              b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - a[l] - c[l];
              d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l] +  drain_above[l] - drain_below[l] + finalSourceSinks_m3s[l];
            } else { // last layer
              // K_up = 0.5*(K_step_ms[l-1] + K_step_ms[l]);
              K_up = K_step_ms[l-1];
              K_down = K_step_ms[l];
              drain_below[l] = K_down;
              drain_above[l] = K_up;
              if(prop_saturated[l]==1.0) {
                drain_above[l] = 0.0;
                drain_below[l] = 0.0;
              }
              capill_below[l] = 0.0;
              if(!freeDrainage) {
                capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_bc);
              }
              a[l] = -1.0*K_up/dZUp[l];
              c[l] = 0.0;
              b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - a[l];
              d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l] + drain_above[l] - drain_below[l] + capill_below[l] + finalSourceSinks_m3s[l];
            }
          }
          //TRIDIAGONAL SOLVING
          tridiagonalSolving_c(a,b,c,d,e,f, Psi_step_t05);
          //MODIFY UNITS OF OUTPUT PSI
          for(int l=0;l<nlayers;l++) Psi_step_t05[l] = Psi_step_t05[l]*mTOMPa; // m to MPa

          //Calculate K and C at t05
          for(int l=0;l<nlayers;l++) {
            if(soilDomains=="single") {
              C_step_05 = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step_t05[l]);
              K_step_05 = waterFluidity[l]*psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step_t05[l]);
            } else {
              C_step_05 = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step_t05[l]);
              K_step_05 = waterFluidity[l]*psi2kVanGenuchtenMicropores_c(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                                         Psi_step_t05[l], Psi_b);
            }
            K_step_ms05[l] = 0.01*K_step_05/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
            C_step_m05[l] = C_step_05*mTOMPa; //From MPa-1 to m-1
          }
          //B. Corrector sub-step also Crank-Nicolson
          for(int l=0;l<nlayers;l++) {
            if(l==0) { //first layer
              K_up = 0.0;
              // K_down = 0.5*(K_step_ms05[l] + K_step_ms05[l+1]);
              K_down = K_step_ms05[l];
              drain_below[l] = K_down;
              if(prop_saturated[l]==1.0) {
                drain_below[l] = 0.0;
              }
              capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_step_m[l+1]);
              a[l] = 0.0;
              c[l] = -1.0*K_down/(2.0*dZDown[l]);
              b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - c[l];
              d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] - c[l]*(Psi_step_m[l] - Psi_step_m[l+1])  - drain_below[l] + finalSourceSinks_m3s[l];
            } else if(l<(nlayers - 1)) {
              // K_up = 0.5*(K_step_ms05[l-1] + K_step_ms05[l]);
              // K_down = 0.5*(K_step_ms05[l] + K_step_ms05[l+1]);
              K_up = K_step_ms05[l-1];
              K_down = K_step_ms05[l];
              drain_above[l] = K_up;
              drain_below[l] = K_down;
              if(prop_saturated[l] == 1.0) {
                drain_above[l] = 0.0;
                drain_below[l] = 0.0;
              }
              capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_step_m[l+1]);
              a[l] = -1.0*K_up/(2.0*dZUp[l]);
              c[l] = -1.0*K_down/(2.0*dZDown[l]);
              b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - a[l] - c[l];
              d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] + a[l]*(Psi_step_m[l - 1] - Psi_step_m[l]) - c[l]*(Psi_step_m[l] - Psi_step_m[l+1]) + drain_above[l] - drain_below[l] + finalSourceSinks_m3s[l];
            } else { // last layer
              // K_up = 0.5*(K_step_ms05[l-1] + K_step_ms05[l]);
              K_up = K_step_ms05[l-1];
              K_down = K_step_ms05[l];
              drain_below[l] = K_down;
              drain_above[l] = K_up;
              if(prop_saturated[l]==1.0) {
                drain_above[l] = 0.0;
                drain_below[l] = 0.0;
              }
              capill_below[l] = 0.0;
              //Changes the boundary conditions allowing capillarity from layer below
              if(!freeDrainage) {
                capill_below[l] = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_bc);
              }
              a[l] = -1.0*K_up/(2.0*dZUp[l]);
              c[l] = 0.0;
              b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - a[l];
              d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] + a[l]*(Psi_step_m[l - 1] - Psi_step_m[l]) + capill_below[l] + drain_above[l] - drain_below[l] + finalSourceSinks_m3s[l];
            }
          }
          //TRIDIAGONAL SOLVING
          tridiagonalSolving_c(a,b,c,d, e, f, Psi_step_t1);
          //MODIFY UNITS OF OUTPUT PSI
          for(int l=0;l<nlayers;l++) Psi_step_t1[l] = Psi_step_t1[l]*mTOMPa; // m to MPa

          //calculate drainage (m3)
          if(freeDrainage) {
            drainage_matrix_step_m3 += drain_below[nlayers -1]*tsubstep;
          } else {
            if(num_saturated < nlayers) {
              drainage_matrix_step_m3 += drain_below[nlayers - num_saturated -1]*tsubstep;
              capillarity_matrix_step_m3 += capill_below[nlayers - num_saturated - 1]*tsubstep;
            }
          }

          //Update Psi and theta
          double res_mm = 0.0;
          for(int l=(nlayers-1);l>=0;l--) {
            // Rcout<<" step "<<s<<" layer " <<l<< " "<< Psi_step[l]<< " to " << Psi_step_t1[l]<<"\n";
            Psi_step[l] = Psi_step_t1[l];
            //If there is residue from below, update Psi
            double new_theta, quasi_sat_theta;
            //Calculate new and sat theta depending on soilDomains
            if(soilDomains=="single") {
              new_theta = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
              quasi_sat_theta = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_quasi_sat);
            } else {
              new_theta = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              quasi_sat_theta = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_quasi_sat);
            }
            //Correct for positive psi
            if(Psi_step[l] > 0.0) {
              if(soilDomains=="single") {
                new_theta = 2.0*theta_sat[l] - new_theta;
              } else {
                new_theta = 2.0*theta_sat_fict[l] - new_theta;
              }
            }
            // Add previous residue
            if(res_mm > 0.0) {
              // Rcout << l << " residue " << res_mm << "\n";
              new_theta += res_mm/(widths[l]*lambda[l]);
              res_mm = 0.0;
            }
            //Manage oversaturation (set res_mm to layer above)
            if(new_theta > quasi_sat_theta) {
              res_mm = std::abs(new_theta - quasi_sat_theta)*widths[l]*lambda[l];
              Psi_step[l] = Psi_quasi_sat;
            } else { // Update psi
              if(soilDomains=="single") {
                Psi_step[l] = theta2psiVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], new_theta);
              } else {
                Psi_step[l] = theta2psiVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l],new_theta);
              }
            }
            // Rcout<<" step "<<s<<" layer " <<l<< " final "<< Psi_step[l]<<"\n";
            //If dual model, update theta and manage excess to macropores
            if((soilDomains=="dual")) {
              //Update theta_micro for the next substep
              theta_micro_step[l] = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              //If needed, add exceeding moisture to the macropores and correct Psi
              double excess_theta_step = 0.0, excess_to_macro= 0.0;
              if(theta_micro_step[l] > theta_b[l]) {
                //Maximum macropore capacity
                double C_macro_step = (1.0 - S_macro_step[l])*(theta_sat[l] - theta_b[l]);
                excess_theta_step = theta_micro_step[l] - theta_b[l];
                double excess_to_macro = std::min(excess_theta_step, C_macro_step);
                // Rcout<< s<<" "<< l << " theta_micro "<< theta_micro_step[l] <<"  "<<theta_b[l]<< " "<< excess_theta_step << "\n";
                res_mm += (excess_theta_step - excess_to_macro)*widths[l]*lambda[l]; //Set remaining to upper layer
                theta_micro_step[l] = theta_b[l];
                Psi_step[l] = Psi_b;
              }
              if(excess_to_macro>0.0){
                lateral_flows_step_mm[l] -= excess_to_macro*widths[l]*lambda[l]; //negative flow (in mm/step)
                matrixExcess_m3s[l] = excess_to_macro*widths[l]*lambda[l]*0.001/tsubstep; //Source of water flowing into macropores (m3/s)
                // Rcout<< excess_theta_step <<" to macropores in "<<l<<"\n";
              } else {
                matrixExcess_m3s[l] = 0.0;
              }
            }
          }
          //Generate saturation excess if there was some residue in the top layer
          saturation_excess_matrix_step_mm += res_mm;


          //Update (micropore) capacitances and conductances for next substep
          for(int l=0;l<nlayers;l++) {
            if(soilDomains=="single") {
              C_step[l] = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
              K_step[l] = waterFluidity[l]*psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
            } else {
              C_step[l] = psi2cVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              K_step[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores_c(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                                         Psi_step[l], Psi_b);
            }
            Psi_step_m[l]= Psi_step[l]/mTOMPa; // MPa to m
            K_step_ms[l] = 0.01*K_step[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
            C_step_m[l] = C_step[l]*mTOMPa; //From MPa-1 to m-1
          }

          if(soilDomains=="dual") {
            //Solve macropore domain by the bisection method
            double infiltration_macropores_substep_mm = 0.0;
            for(int l=0;l<nlayers;l++) {
              //Set imbibition as sink for macropore
              double sourceSink_macro_m3s = matrixExcess_m3s[l] - matrixImbibition_m3s[l];

              //Infiltration into macropores
              if(l==0) {
                if(infiltration_remaining_macropores_step_mm > 0.0) {
                  //Maximum infiltration in this time step
                  double infiltration_macropores_substep_m3s = infiltration_remaining_macropores_step_mm*((double) nsteps)*mm_day_2_m3_s; //From mm/step to m3s
                  sourceSink_macro_m3s += infiltration_macropores_substep_m3s;
                  infiltration_macropores_substep_mm = infiltration_remaining_macropores_step_mm/((double) nsubsteps);
                  // Rcout << "infiltration step " << infiltration_macropores_substep_mm<<"\n";
                  infiltration_remaining_macropores_step_mm -= infiltration_macropores_substep_mm;
                }
              }
              //Correct saturated correction with source/sinks
              // saturated_macropore_correction_m3s[l] = std::max(0.0, saturated_macropore_correction_m3s[l] - sourceSink_macro_m3s);
              capillarity_macropores_step_m3 += tsubstep*saturated_macropore_correction_m3s[l];
              //Add to source/sinks
              sourceSink_macro_m3s += saturated_macropore_correction_m3s[l];
              if(l==0) { //first layer
                K_up = 0.0;
              } else {
                K_up = Kmacro_step_ms[l-1]; //Get last updated value
              }
              double ksat_i = Ksat_ms[l];
              double ksat_b_i = Ksat_b_ms[l];
              //If layer is below the water table, gravitational fluxes are set to zero
              if(prop_saturated[l]==1.0) {
                K_up = 0.0;
                ksat_i = 0.0;
                ksat_b_i = 0.0;
              }
              //Find solution for half substep
              double S_t1 = rootFindingMacropores_c(S_macro_step[l], K_up, ksat_i, ksat_b_i, kin_exp,
                                                    e_macro[l], lambda[l], dZ_m[l], sourceSink_macro_m3s, tsubstep);
              // Rcout << "S " << S_macro_step[l] << " source/sink step " << sourceSink_macro_m3s<< " S1 "<< S_t1<<"\n";
              //
              //Update macropore conductances for next step (sets K_up for next layer)
              Kmacro_step_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_t1, kin_exp);
              //Update S_macro_step using full substep
              S_macro_step[l] = S_t1;
              // S_macro_step[l] = S_macro_step[l] + (tsubstep/(e_macro[l]*lambda[l]*dZ_m[l]))*((K_up - Kmacro_step_ms[l]) + sourceSink_macro_m3s);
            }
            //Drainage
            if(freeDrainage) {
              drainage_macropores_step_m3 += Kmacro_step_ms[nlayers-1]*tsubstep;
            } else {
              if(num_saturated<nlayers) {
                double flow = Kmacro_step_ms[nlayers - num_saturated -1];
                // Rcout<< " ss " << ss << " drainage flow " << flow << "\n";
                drainage_macropores_step_m3 += flow*tsubstep;
              }

            }
            //Update theta_macro for the next step
            double res_mm = 0.0;
            for(int l=(nlayers-1);l>=0;l--) {
              theta_macro_step[l] = S_macro_step[l]*e_macro[l] + res_mm/(widths[l]*lambda[l]);
              // Rcout << " layer "<< l <<" res_mm " << res_mm << " theta_macro "<< theta_macro_step[l] <<"\n";
              // If oversaturation of macroporosity occurs
              if(theta_macro_step[l] > e_macro[l]) {
                res_mm = widths[l]*(theta_macro_step[l] - e_macro[l])*lambda[l]; //residue in mm for layers above
                // Rcout << " layer "<< l <<" residue "<< res_mm<<"\n";
                theta_macro_step[l] = e_macro[l];
                S_macro_step[l] = 1.0; //Correct S_macro to saturation
                Kmacro_step_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_macro_step[l], kin_exp);
                if(l==0) { // If there has been infiltration, add it back to infiltration target
                  double min_diff = std::min(res_mm, infiltration_macropores_substep_mm);
                  infiltration_macropores_substep_mm -= min_diff;
                  infiltration_remaining_macropores_step_mm += min_diff;
                  res_mm -= min_diff;
                  // Rcout << " layer "<< l <<" res_mm " << res_mm << " min_diff "<< min_diff << " infiltration "<< infiltration_macropores_substep_mm<<"\n";
                }
              } else {
                res_mm = 0.0;
              }
            }
            //Generate saturation excess if there was some residue in the top layer
            saturation_excess_macropores_step_mm += res_mm;
            // Store infiltration
            infiltration_macropores_step_mm += infiltration_macropores_substep_mm;
            // Rcout<< "infiltration step " << infiltration_macropores_step_mm<<"\n";
          }

        }

        //Update theta, W and final volumes
        Vfin_mm = 0.0;
        Vfin_micro_mm = 0.0;
        Vfin_macro_mm = 0.0;
        for(int l=0;l<nlayers;l++) {
          if(soilDomains=="single") {
            Theta[l] = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
          } else {
            Theta[l] = theta_micro_step[l] + theta_macro_step[l];
            Vfin_micro_mm += theta_micro_step[l]*widths[l]*lambda[l];
            Vfin_macro_mm += theta_macro_step[l]*widths[l]*lambda[l];
          }
          Vfin_mm += Theta[l]*widths[l]*lambda[l];
        }

        double max_abs_correction_step_mm;
        double sum_source_sink_def_m3s = std::accumulate(source_sink_def_m3s.begin(), source_sink_def_m3s.end(), 0.0);
        if(soilDomains=="single") {
          // Rcout << s << " "<< nsubsteps<< " sat_ex " << saturation_excess_matrix_mm << " dr " << drainage_matrix_step_m3*m3_2_mm << " cap " << capillarity_matrix_step_m3*m3_2_mm<<"\n";
          double balance_step_mm = m3_2_mm*(tstep*sum_source_sink_def_m3s  + capillarity_matrix_step_m3 - drainage_matrix_step_m3) - saturation_excess_matrix_step_mm;
          matrix_correction_step_mm = balance_step_mm + Vini_step_mm - Vfin_mm;
          max_abs_correction_step_mm = std::abs(matrix_correction_step_mm);
          // Rcout << s << " "<< nsubsteps<<" ini "<< Vini_step_mm <<" fin " <<  Vfin_mm << " dif "<< Vfin_mm - Vini_step_mm
          //       << " dra " << m3_2_mm*drainage_matrix_step_m3 <<" cap " << m3_2_mm*capillarity_matrix_step_m3 <<"  bal "<< balance_step_mm<< " corr "<< matrix_correction_step_mm<<"\n";
        } else {
          //Add remaining target to infiltration excess
          infiltration_excess_macropores_step_mm = infiltration_remaining_macropores_step_mm;
          double sum_lateral_flows_step_mm = 0.0;
          for(int l=0;l<nlayers;l++) sum_lateral_flows_step_mm += lateral_flows_step_mm[l];
          double balance_micro_step_mm = m3_2_mm*(tstep*sum_source_sink_def_m3s + capillarity_matrix_step_m3 - drainage_matrix_step_m3) + sum_lateral_flows_step_mm - saturation_excess_matrix_step_mm;
          double balance_macro_step_mm = infiltration_macropores_step_mm + m3_2_mm*(capillarity_macropores_step_m3 - drainage_macropores_step_m3) - sum_lateral_flows_step_mm - saturation_excess_macropores_step_mm;
          // Correct possible mismatch between balance and volume change
          matrix_correction_step_mm = balance_micro_step_mm + Vini_step_micro_mm - Vfin_micro_mm;
          macropore_correction_step_mm = balance_macro_step_mm + Vini_step_macro_mm - Vfin_macro_mm;
          // Rcout << s << " "<< nsubsteps<<" micro ini "<< Vini_step_micro_mm <<" fin " <<  Vfin_micro_mm << " dif "<< Vfin_micro_mm - Vini_step_micro_mm << "  bal "<< balance_micro_step_mm<< " corr "<< matrix_correction_step_mm<< "\n";
          // Rcout << s << " "<< nsubsteps<<" macro ini "<< Vini_step_macro_mm <<" fin " <<  Vfin_macro_mm << " dif "<< Vfin_macro_mm - Vini_step_macro_mm <<
          //   " inf " << infiltration_macropores_step_mm<< " sat exc " << saturation_excess_macropores_step_mm << " dra " << m3_2_mm*drainage_macropores_step_m3 <<" cap " << m3_2_mm*capillarity_macropores_step_m3 <<
          //   " lat " << sum(lateral_flows_step_mm) << " bal "<< balance_macro_step_mm<<  " corr "<< macropore_correction_step_mm<< "\n";
          max_abs_correction_step_mm = std::max(std::abs(matrix_correction_step_mm), std::abs(macropore_correction_step_mm));
        }
        if((max_abs_correction_step_mm > 0.001)) {
          nsubsteps = nsubsteps*2;
          if(nsubsteps >= max_nsubsteps) {
            nsubsteps = max_nsubsteps;
            cont = false;
          }
        } else {
          cont = false;
        }
      }


      //Add drainage correction
      drainage_matrix_step_m3 += std::max(0.0, matrix_correction_step_mm*mm_2_m3);
      capillarity_matrix_step_m3 += std::max(0.0, -1.0*matrix_correction_step_mm*mm_2_m3);
      if(soilDomains=="dual") {
        drainage_macropores_step_m3 += std::max(0.0, macropore_correction_step_mm*mm_2_m3);
        capillarity_macropores_step_m3 += std::max(0.0, -1.0*macropore_correction_step_mm*mm_2_m3);
      }

      //Add to totals
      drainage_matrix_m3 += drainage_matrix_step_m3;
      matrix_correction_mm += matrix_correction_step_mm;
      saturation_excess_matrix_mm += saturation_excess_matrix_step_mm;
      capillarity_matrix_m3 += capillarity_matrix_step_m3;
      if(soilDomains=="dual") {
        drainage_macropores_m3 += drainage_macropores_step_m3;
        capillarity_macropores_m3 += capillarity_macropores_step_m3;
        infiltration_macropores_mm +=infiltration_macropores_step_mm;
        macropore_correction_mm += macropore_correction_step_mm;
        saturation_excess_macropores_mm += saturation_excess_macropores_step_mm;
        infiltration_excess_macropores_mm += infiltration_excess_macropores_step_mm;
      }

      for(int l=0;l<nlayers;l++) {
        matrix_macropore_flows_mm[l] += lateral_flows_step_mm[l];
        //Copy back step variables
        C[l] = C_step[l];
        C_m[l] = C_step_m[l];
        K[l] = K_step[l];
        K_ms[l] = K_step_ms[l];
        Psi[l] = Psi_step[l];
        Psi_m[l] = Psi_step_m[l];
        // Rcout << s << " "<< l << " psi " << Psi[l] << " K " << K[l]  << " theta " << Theta[l] <<"\n";
        if(soilDomains=="dual") {
          S_macro[l] = S_macro_step[l];
          Kmacro_ms[l] = Kmacro_step_ms[l];
          theta_macro[l] = theta_macro_step[l];
        }
        theta_micro[l] = theta_micro_step[l];
        // Rcout << s << " "<< l << " th_mi " << theta_micro[l] << " K_mi " << K_ms[l] << " th_ma " << theta_macro[l] << " K_ma " << Kmacro_ms[l]<< " mat_im " << matrixImbibition_m3s[l] << " mat_ex " << matrixExcess_m3s[l] <<"\n";
      }
      //Prepare for next step
      Vini_step_macro_mm = Vfin_macro_mm;
      Vini_step_micro_mm = Vfin_micro_mm;
      Vini_step_mm = Vfin_mm;
    }

    double drainage_matrix_mm = drainage_matrix_m3*1000.0; //m3/m2 to mm/m2
    double capillarity_matrix_mm = capillarity_matrix_m3*1000.0; //m3/m2 to mm/m2
    double capillarity_macropores_mm = capillarity_macropores_m3*1000.0; //m3/m2 to mm/m2


    for(int l=0; l< nlayers;l++) {
      soil.setW(l,Theta[l]/soil.getThetaFC(l));
      //   Rcout << "Final "<<l<< " W " << W[l]<<" Theta " << Theta[l]<< " theta_sat "<< theta_sat[l] <<" theta_micro "<< theta_micro[l] << " theta_b "<< theta_b[l] << " theta_macro "<< theta_macro[l] << " S " << S_macro[l]<<"\n";
    }

    //Correct overestimation of capillarity by using deep drainage
    double dif_mm = std::min(capillarity_matrix_mm,  drainage_matrix_mm);
    capillarity_matrix_mm -= dif_mm;
    drainage_matrix_mm  -= dif_mm;

    //Output
    if(soilDomains=="dual") {
      double drainage_macropores_mm = drainage_macropores_m3*1000.0; //m3/m2 to mm/m2
      double dif_mm = std::min(capillarity_macropores_mm,  drainage_macropores_mm);
      capillarity_macropores_mm -= dif_mm;
      drainage_macropores_mm  -= dif_mm;
      double runoff_mm = infiltration_target_macropores_mm + infiltration_excess_macropores_mm + saturation_excess_matrix_mm + saturation_excess_macropores_mm;
      SWBres.localSourceSinks_mm = std::accumulate(sourceSink.begin(), sourceSink.end(), 0.0);
      SWBres.lateralSourceSinks_mm = std::accumulate(lateralFlows.begin(), lateralFlows.end(), 0.0);
      SWBres.runoff_mm = runoff_mm;
      SWBres.matrixMacroporeFlow_mm = std::accumulate(matrix_macropore_flows_mm.begin(), matrix_macropore_flows_mm.end(), 0.0);;
      SWBres.infiltrationMatrix_mm = infiltration_matrix_mm;
      SWBres.infiltrationMacropores_mm = infiltration_macropores_mm;
      SWBres.infiltration_mm = infiltration_matrix_mm + infiltration_macropores_mm;
      SWBres.infiltrationExcessMatrix_mm = infiltration_excess_matrix_mm;
      SWBres.infiltrationExcessMacropores_mm = infiltration_excess_macropores_mm;
      SWBres.infiltrationExcess_mm = infiltration_target_macropores_mm + infiltration_excess_macropores_mm;
      SWBres.saturationExcessMatrix_mm = saturation_excess_matrix_mm;
      SWBres.saturationExcessMacropores_mm = saturation_excess_macropores_mm;
      SWBres.saturationExcess_mm = saturation_excess_matrix_mm + saturation_excess_macropores_mm;
      SWBres.drainageMatrix_mm = drainage_matrix_mm;
      SWBres.drainageMacropores_mm = drainage_macropores_mm;
      SWBres.deepDrainage_mm = drainage_matrix_mm + drainage_macropores_mm;
      SWBres.runoff_mm = runoff_mm;
      SWBres.capillarityMatrix_mm = capillarity_matrix_mm;
      SWBres.capillarityMacropores_mm = capillarity_macropores_mm;
      SWBres.capillarityRise_mm = capillarity_matrix_mm + capillarity_macropores_mm;
      SWBres.correctionMatrix_mm = matrix_correction_mm;
      SWBres.correctionMacropores_mm = macropore_correction_mm;
      SWBres.correction_mm = matrix_correction_mm + macropore_correction_mm;
      SWBres.matrixVolumeChange_mm = Vfin_micro_mm - Vini0_micro_mm;
      SWBres.macroporeVolumeChange_mm = Vfin_macro_mm - Vini0_macro_mm;
      SWBres.volumeChange_mm = Vfin_micro_mm - Vini0_micro_mm  + Vfin_macro_mm - Vini0_macro_mm;
      SWBres.substeps = total_nsubsteps;
    } else if(soilDomains=="single") {
      double runoff_mm = infiltration_excess_matrix_mm  + saturation_excess_matrix_mm;
      SWBres.localSourceSinks_mm = std::accumulate(sourceSink.begin(), sourceSink.end(), 0.0);
      SWBres.lateralSourceSinks_mm = std::accumulate(lateralFlows.begin(), lateralFlows.end(), 0.0);
      SWBres.infiltration_mm = infiltration_matrix_mm;
      SWBres.infiltrationExcess_mm = infiltration_excess_matrix_mm;
      SWBres.saturationExcess_mm = saturation_excess_matrix_mm;
      SWBres.runoff_mm = runoff_mm;
      SWBres.deepDrainage_mm = drainage_matrix_mm;
      SWBres.capillarityRise_mm = capillarity_matrix_mm;
      SWBres.correction_mm = matrix_correction_mm;
      SWBres.volumeChange_mm = Vfin_mm - Vini0_mm;
      SWBres.substeps = total_nsubsteps;
    }
  }
}

 