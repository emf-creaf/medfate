// [[Rcpp::interfaces(r,cpp)]]
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <RcppArmadillo.h>
#include "modelInput_c.h"
#include "lowlevel_structures_c.h"
#include "soil_c.h"
#include "hydraulics_c.h"
#include "hydrology_c.h"
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
#include "numerical_solving_c.h"



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
  // Fix B5: Exact derivative of Ksat*t + A*log(1 + x/A) - x with respect to x is 1/(1 + x/A) - 1
  double A = std::abs(psi_w)*delta_theta;
  double fder = (1.0/(1.0 + (x/A))) - 1.0;
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
  if(Ksat <= 0.0 || t <= 0.0) return(0.0);
  double delta_theta = theta_sat - theta_dry;
  if(delta_theta <= 1e-12 || std::abs(psi_w) <= 1e-12) return(Ksat*t);
  double x,x1,e,fx,fx1;
  x1 = std::max(Ksat*t, 1e-6); // non-zero initial guess (since f'(0) == 0)
  e = 0.001; // accuracy in mm
  int cnt = 0;
  int mxiter = 100;
  do {
    x=x1; /*make x equal to the last calculated value of  x1*/
    fx=fGreenAmpt_c(x, t, psi_w, Ksat, delta_theta);            //simplifying f(x)to fx
    fx1=fGreenAmptDer_c(x, t, psi_w, Ksat, delta_theta);        //simplifying fprime(x) to fx1
    x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
    if(x1 <= 0.0) x1 = 0.5*x;
    cnt++;
  } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
  return(x1);
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
  double rho = airDensity_c(tday, atmosphericPressure_c(elevation));
  double ten = (86400.0*tday*rho*1013.86*1e-6/100.0); //ten can be negative if temperature is below zero
  double ren = (rad*(LgroundSWR/100.0))*(0.1); //90% albedo of snow
  double melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
  return(melt);
}


void waterInputs_c(WaterInputs_COMM& waterInputs,
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
  waterInputs.rain = rain;
  waterInputs.snow = snow;
  waterInputs.interception = Interception;
  waterInputs.netrain = NetRain;
  waterInputs.melt = melt;
  if(modifyInput) {
    x.snowpack = swe;
  }
}

void agricultureWaterInputs_c(WaterInputs_COMM& waterInputs,
                              AgricultureModelInput& x,
                              double prec, double tday, double rad, double elevation,
                              double LgroundSWR, 
                              bool modifyInput = true) {
  
  double swe = x.snowpack; //snow pack
  
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
    NetRain = rain - Interception; 
  }
  waterInputs.rain = rain;
  waterInputs.snow = snow;
  waterInputs.interception = Interception;
  waterInputs.netrain = NetRain;
  waterInputs.melt = melt;
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
  // Fix B4 & P3: Exact 0 short-circuit when dry/draining + Newton-Raphson solver
  double eff_cap = e_macro*lambda*dZ_m;
  if(eff_cap <= 0.0) return(0.0);
  double alpha_dt = tstep/eff_cap;
  double R = S_t + alpha_dt*(K_up + sourceSink_macro_m3s);
  if(R <= 0.0) return(0.0);
  double dK = Ksat_ms - Ksat_b_ms;
  if(dK <= 0.0) return(R);
  double beta = alpha_dt*dK;
  double c = (S_t > 1e-6) ? S_t : std::min(R, 1.0);
  double tol = 1e-10;
  for(int N = 0; N < Nmax; N++) {
    if(c < 0.0) c = 0.0;
    double c_km1 = (c > 0.0) ? pow(c, kin_exp - 1.0) : 0.0;
    double f_c = c + beta*(c*c_km1) - R;
    if(std::abs(f_c) < tol) return(c);
    double df_c = 1.0 + beta*kin_exp*c_km1;
    double step = f_c/df_c;
    c -= step;
    if(std::abs(step) < tol) return(c > 0.0 ? c : 0.0);
  }
  return(c > 0.0 ? c : 0.0);
}

// Fused Van Genuchten evaluation (computes theta, C, and K in 2 pow() calls instead of 11)
inline void evalVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat,
                               double Ksat, double waterFluidity, double psi,
                               bool is_dual_micropore, double Ksat_b, double Psi_b,
                               double& theta_out, double& C_out, double& K_out) {
  if(psi >= 0.0) {
    theta_out = theta_sat;
    C_out = 0.0;
    K_out = waterFluidity*(is_dual_micropore ? Ksat_b : Ksat);
    return;
  }
  double m = 1.0 - (1.0/n);
  double abs_psi = std::abs(psi);
  double u = alpha*abs_psi;
  double u_n1 = pow(u, n - 1.0);
  double x = u*u_n1;
  double Se = pow(1.0 + x, -m);
  double w = u_n1*Se; // exact algebraic identity: (1 - Se^(1/m))^m == u^(n-1)*Se
  theta_out = theta_res + (theta_sat - theta_res)*Se;
  C_out = (theta_sat - theta_res)*n*m*alpha*w/(1.0 + x);
  if(is_dual_micropore && (psi >= Psi_b)) {
    K_out = waterFluidity*Ksat_b;
  } else {
    double one_minus_w = 1.0 - w;
    K_out = waterFluidity*Ksat*sqrt(Se)*one_minus_w*one_minus_w;
  }
}


/*=============================================================================
 * Implementation of hydrology routines using C++ code
 * (Modified with Celia et al. 1990 Mixed theta-psi Modified Picard Solver)
 *=============================================================================*/

void soilWaterBalance_inner_c(SoilWaterBalance_RESULT &SWBres, SoilWaterBalance_COMM &SWBcomm, Soil &soil, 
                              double rainfallInput, double rainfallIntensity, double snowmelt, const std::vector<double> &sourceSink, 
                              double runon, const std::vector<double> &lateralFlows, double waterTableDepth,
                              std::string infiltrationMode = "GreenAmpt1911", double infiltrationCorrection = 5.0, 
                              std::string soilDomains = "buckets", 
                              int nsteps = 24, int max_nsubsteps = 3600, double richards_tolerance = 0.001) {
  
  if((soilDomains!="single") && (soilDomains!="dual")  && (soilDomains!="buckets") ) throw medfate::MedfateInternalError("Unrecognized soilDomain value");
  
  const std::vector<double>& widths = soil.getWidths();
  const std::vector<double>& macro = soil.getMacro();
  const std::vector<double>& rfc = soil.getRFC();
  const std::vector<double>& Ksat_ori = soil.getKsat();
  
  
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
  std::vector<double>& source_sink_def_mm = SWBcomm.source_sink_def_mm;
  for(int l=0;l<nlayers;l++) {
    source_sink_def_mm[l] = sourceSink[l];
  }
  
  //Add infiltration to matrix def source/sinks
  std::vector<double>& IVec_mm = SWBcomm.IVec_mm;
  if(soilDomains!="dual") {
    infiltrationRepartition_c(infiltration_matrix_mm,
                              IVec_mm,
                              widths, macro);
    for(int l=0;l<nlayers;l++) {
      source_sink_def_mm[l] += IVec_mm[l];
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
        if(l<(nlayers-1)) {
          //update source_sink_def_mm adding the excess to the infiltrating water (saturated flow)
          source_sink_def_mm[l+1] = source_sink_def_mm[l+1] + std::max(Wn - Water_FC[l],0.0);
        } else {
          saturation_excess_matrix_mm = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer using field capacity
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
      double Kdrain = 10.0*(Ksat_ori[nlayers-1]*(1.0 - (rfc[nlayers-1]/100.0))/cmdTOmmolm2sMPa);
      double maxDrainage = Kdrain;
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
    
    //Set time steps
    double tstep = 86400.0/((double) nsteps);
    double tsubstep = tstep;
    double rainfallIntensity_step = rainfallIntensity*24.0/((double) nsteps); //mm/step
    
    std::vector<double>& matrix_macropore_flows_mm = SWBcomm.matrix_macropore_flows_mm;
    std::vector<double>& source_sink_def_m3s = SWBcomm.source_sink_def_m3s;
    std::vector<double>& matrixImbibition_m3s = SWBcomm.matrixImbibition_m3s;
    std::vector<double>& matrixExcess_m3s = SWBcomm.matrixExcess_m3s;
    std::vector<double>& saturated_matrix_correction_m3s = SWBcomm.saturated_matrix_correction_m3s;
    std::vector<double>& saturated_macropore_correction_m3s = SWBcomm.saturated_macropore_correction_m3s;
    for(int l=0;l<nlayers;l++) {
      matrix_macropore_flows_mm[l] = 0.0;
      source_sink_def_m3s[l] =  source_sink_def_mm[l]*mm_day_2_m3_s; 
      matrixImbibition_m3s[l] = 0.0;
      matrixExcess_m3s[l] = 0.0;
      saturated_matrix_correction_m3s[l] = 0.0;
      saturated_macropore_correction_m3s[l] = 0.0;
    }
    
    //Estimate layer interfaces
    std::vector<double>& dZ_m = SWBcomm.dZ_m;
    std::vector<double>& dZUp = SWBcomm.dZUp;
    std::vector<double>& dZDown = SWBcomm.dZDown;
    std::vector<double>& lambda = SWBcomm.lambda;
    for(int l=0;l<nlayers;l++) {
      dZ_m[l] = widths[l]*0.001; //mm to m
      lambda[l] = 1.0 - (rfc[l]/100.0);
    }
    for(int l=0;l<nlayers;l++) {
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
    
    std::vector<double>& prop_saturated = SWBcomm.prop_saturated;
    int num_saturated = 0;
    bool freeDrainage = true;
    for(int i=0;i<nlayers;i++) prop_saturated[i] = 0.0;
    if(!std::isnan(waterTableDepth)) {
      freeDrainage = (waterTableDepth > soilDepth);
      double sZ = 0.0;
      for(int i=0;i<nlayers;i++){
        prop_saturated[i] = std::min(1.0, std::max(0.0,((sZ + widths[i]) - waterTableDepth)/widths[i]));
        if(prop_saturated[i]==1.0) num_saturated++;
        sZ += widths[i];
      }
    }
    
    //boundary condition of water table (if freeDrainage = FALSE)
    double Psi_bc = 0.0;
    double Psi_quasi_sat = -0.0000001;
    //Boundary water potential for dual porosity
    double Psi_b = -0.1*mTOMPa; // 10 cm = 0.1 m
    double kin_exp = 2.23; //Kinematic exponent reflecting macropore size distribution and tortuosity
    
    //Retrieve VG parameters
    const std::vector<double>&  n = soil.getVG_n();
    const std::vector<double>&  alpha = soil.getVG_alpha();
    const std::vector<double>&  theta_res = soil.getVG_theta_res();
    const std::vector<double>&  theta_sat = soil.getVG_theta_sat();
    
    //Get initial theta and psi state
    std::vector<double>& Theta = SWBcomm.Theta;
    std::vector<double>& Psi = SWBcomm.Psi;
    for(int l=0;l<nlayers;l++) {
      Theta[l] = soil.getTheta(l);
      Psi[l] = soil.getPsi(l);
    }
    
    //Microporosity or single domain
    std::vector<double>& theta_micro = SWBcomm.theta_micro;
    std::vector<double>& theta_b = SWBcomm.theta_b;
    std::vector<double>& theta_macro = SWBcomm.theta_macro;
    std::vector<double>& theta_sat_fict = SWBcomm.theta_sat_fict;
    std::vector<double>& Ksat_b = SWBcomm.Ksat_b;
    std::vector<double>& Ksat_b_ms = SWBcomm.Ksat_b_ms;
    std::vector<double>& Ksat = SWBcomm.Ksat;
    std::vector<double>& Ksat_ms = SWBcomm.Ksat_ms;
    std::vector<double>& K = SWBcomm.K;
    std::vector<double>& C = SWBcomm.C;
    std::vector<double>& Psi_m = SWBcomm.Psi_m;
    std::vector<double>& K_ms = SWBcomm.K_ms;
    std::vector<double>& Kbc = SWBcomm.Kbc;
    std::vector<double>& Kbc_ms = SWBcomm.Kbc_ms;
    std::vector<double>& C_m = SWBcomm.C_m;
    std::vector<double>& S_macro = SWBcomm.S_macro;
    std::vector<double>& e_macro = SWBcomm.e_macro;
    std::vector<double>& Kmacro_ms = SWBcomm.Kmacro_ms;
    
    for(int l=0;l<nlayers;l++) {
      theta_micro[l] = 0.0;
      theta_b[l] = 0.0;
      theta_macro[l] = 0.0;
      theta_sat_fict[l] = 0.0;
      Ksat_b[l] = 0.0;
      Ksat_b_ms[l] = 0.0;
      Ksat[l] = 0.0;
      K[l] = 0.0;
      C[l] = 0.0;
      Psi_m[l] = 0.0;
      K_ms[l] = 0.0;
      Kbc[l] = 0.0;
      Kbc_ms[l] = 0.0;
      C_m[l] = 0.0;
      if(soilDomains=="dual") {
        S_macro[l] = 0.0;
        e_macro[l] = 0.0;
        Kmacro_ms[l] = 0.0;
      }
    }
    
    // Hoisted layer-invariant quantities (Fix P2)
    std::vector<double> quasi_sat_theta_vec(nlayers, 0.0);
    std::vector<double> Ksat_fict_ms_vec(nlayers, 0.0);
    std::vector<double> D_theta_b_m2s_vec(nlayers, 0.0);
    
    std::vector<double>& waterFluidity = SWBcomm.waterFluidity;
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
        theta_sat_fict[l] = theta_sat[l];
        theta_micro[l] = Theta[l];
        Psi[l] = theta2psiVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]);
        double th_dummy;
        evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Ksat[l], waterFluidity[l],
                           Psi[l], false, 0.0, Psi_b, th_dummy, C[l], K[l]);
        evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Ksat[l], waterFluidity[l],
                           Psi_bc, false, 0.0, Psi_b, th_dummy, th_dummy, Kbc[l]);
        quasi_sat_theta_vec[l] = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_quasi_sat);
      } else {
        theta_sat_fict[l] = theta_sat[l] - macro[l];
        //Matching theta point and theta partitioning
        //This ensures theta_b < theta_sat - macro & e_macro > macro
        theta_b[l] = psi2thetaVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat_fict[l], Psi_b);
        e_macro[l] = theta_sat[l] - theta_b[l];
        Ksat_b[l] = psi2kVanGenuchten_c(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_b);
        theta_micro[l] = std::min(Theta[l], theta_b[l]);
        theta_macro[l] = std::max(Theta[l] - theta_b[l], 0.0);
        //Water potential, conductivity and capacitance in the micropore domain according to the modified retention
        Psi[l] = theta2psiVanGenuchten_c(n[l],alpha[l],theta_res[l], theta_sat_fict[l], theta_micro[l]);
        double th_dummy;
        evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Ksat[l], waterFluidity[l],
                           Psi[l], true, Ksat_b[l], Psi_b, th_dummy, C[l], K[l]);
        evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Ksat[l], waterFluidity[l],
                           Psi_bc, true, Ksat_b[l], Psi_b, th_dummy, th_dummy, Kbc[l]);
        //Effective saturation and conductivity in the macropore domain
        S_macro[l] = theta_macro[l]/e_macro[l];
        Ksat_b_ms[l] = 0.01*waterFluidity[l]*Ksat_b[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
        Kmacro_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_macro[l], kin_exp);
        quasi_sat_theta_vec[l] = psi2thetaVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_quasi_sat);
        Ksat_fict_ms_vec[l] = psi2kVanGenuchtenMicropores_c(Ksat_b_ms[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 0.0, Psi_b);
        D_theta_b_m2s_vec[l] = psi2DVanGenuchten_c(Ksat_fict_ms_vec[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_b);
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
    
    std::vector<double>& a = SWBcomm.a;
    std::vector<double>& b = SWBcomm.b;
    std::vector<double>& c = SWBcomm.c;
    std::vector<double>& d = SWBcomm.d;
    std::vector<double>& e = SWBcomm.e;
    std::vector<double>& f = SWBcomm.f;
    std::vector<double>& C_step = SWBcomm.C_step;
    std::vector<double>& C_step_m = SWBcomm.C_step_m;
    std::vector<double>& K_step_ms = SWBcomm.K_step_ms;
    std::vector<double>& K_step = SWBcomm.K_step;
    std::vector<double>& Psi_step = SWBcomm.Psi_step;
    std::vector<double>& Psi_step_m = SWBcomm.Psi_step_m;
    std::vector<double>& Psi_step_t1 = SWBcomm.Psi_step_t1;
    
    std::vector<double>& S_macro_step = SWBcomm.S_macro_step;
    std::vector<double>& Kmacro_step_ms = SWBcomm.Kmacro_step_ms;
    std::vector<double>& theta_macro_step = SWBcomm.theta_macro_step;
    std::vector<double>& theta_micro_step = SWBcomm.theta_micro_step;
    std::vector<double>& finalSourceSinks_m3s = SWBcomm.finalSourceSinks_m3s;
    std::vector<double>& capill_below = SWBcomm.capill_below;
    std::vector<double>& drain_above = SWBcomm.drain_above;
    std::vector<double>& drain_below = SWBcomm.drain_below;
    std::vector<double>& lateral_flows_step_mm = SWBcomm.lateral_flows_step_mm;
    for(int l=0;l<nlayers;l++) {
      C_step[l] = 0.0;
      C_step_m[l] = 0.0;
      K_step_ms[l] = 0.0;
      K_step[l] = 0.0;
      Psi_step[l] = 0.0;
      Psi_step_m[l] = 0.0;
      Psi_step_t1[l] = 0.0;
      S_macro_step[l] = 0.0;
      Kmacro_step_ms[l] = 0.0;
      theta_macro_step[l] = 0.0;
      theta_micro_step[l] = 0.0;
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
    
    double sum_source_sink_def_m3s = std::accumulate(source_sink_def_m3s.begin(), source_sink_def_m3s.end(), 0.0);
    int total_nsubsteps = 0;
    const int max_picard = 15;
    const double picard_theta_tol = richards_tolerance*1.0e-3;
    
    for(int s =0;s<nsteps;s++) {
      
      //Target infiltration for macropores
      infiltration_target_macropores_step_mm = std::min(infiltration_target_macropores_mm, rainfallIntensity_step);
      infiltration_target_macropores_mm -= infiltration_target_macropores_step_mm;
      
      saturation_excess_macropores_step_mm = 0.0;
      infiltration_excess_macropores_step_mm = 0.0;
      saturation_excess_matrix_step_mm = 0.0;
      drainage_matrix_step_m3 = 0.0;
      capillarity_matrix_step_m3 = 0.0;
      capillarity_macropores_step_m3 = 0.0;
      drainage_macropores_step_m3 = 0.0;
      infiltration_macropores_step_mm = 0.0;
      for(int l=0;l<nlayers;l++) {
        lateral_flows_step_mm[l] = 0.0;
        matrixImbibition_m3s[l] = 0.0;
        matrixExcess_m3s[l] = 0.0;
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
      
      tsubstep = tstep;
      infiltration_remaining_macropores_step_mm = infiltration_target_macropores_step_mm;
      
      // 1. Update micropore imbibition BEFORE water-table saturation correction (Fix B6)
      if(soilDomains=="dual") {
        if(micropore_imbibition) {
          for(int l=0;l<nlayers;l++) {
            if(S_macro_step[l] > 1e-12 && theta_micro_step[l] < theta_b[l]) {
              double D_theta_micro_m2s = psi2DVanGenuchten_c(Ksat_fict_ms_vec[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l],
                                                             Psi_step[l]);
              double imbibitionRate = microporeImbibitionRate_c(theta_b[l], theta_micro_step[l],
                                                                D_theta_b_m2s_vec[l], D_theta_micro_m2s,
                                                                S_macro_step[l]);
              double max_from_macro = (S_macro_step[l]*e_macro[l])/tsubstep;
              double max_to_micro = (theta_b[l] - theta_micro_step[l])/tsubstep;
              imbibitionRate = std::min(imbibitionRate, std::min(max_from_macro, max_to_micro));
              matrixImbibition_m3s[l] = dZ_m[l]*lambda[l]*imbibitionRate;
              lateral_flows_step_mm[l] += imbibitionRate*widths[l]*lambda[l]*tsubstep;
            } else {
              matrixImbibition_m3s[l] = 0.0;
            }
          }
        }
      }
      
      // 2. Correction for saturation
      for(int l=0;l<nlayers;l++) {
        saturated_matrix_correction_m3s[l] = 0.0;
        saturated_macropore_correction_m3s[l] = 0.0;
        if(prop_saturated[l]>0.0) {
          double theta_l = theta_micro_step[l];
          double quasi_sat_theta_l = quasi_sat_theta_vec[l];
          if(soilDomains =="single") {
            saturated_matrix_correction_m3s[l] = std::max(0.0, ((prop_saturated[l]*quasi_sat_theta_l - theta_l)*lambda[l]*widths[l]*0.001)/tsubstep);
          } else {
            saturated_matrix_correction_m3s[l] = std::max(0.0, ((prop_saturated[l]*theta_b[l] - theta_l)*lambda[l]*widths[l]*0.001)/tsubstep);
            saturated_macropore_correction_m3s[l] = std::max(0.0, (std::max(0.0, 1.0 - S_macro_step[l])*(prop_saturated[l]*quasi_sat_theta_l - theta_b[l])*lambda[l]*widths[l]*0.001)/tsubstep);
          }
          saturated_matrix_correction_m3s[l] = std::min(Ksat_ms[l], saturated_matrix_correction_m3s[l]);
          saturated_macropore_correction_m3s[l] = std::min(Ksat_ms[l], saturated_macropore_correction_m3s[l]);
          saturated_matrix_correction_m3s[l] = std::max(0.0, saturated_matrix_correction_m3s[l] - source_sink_def_m3s[l] - matrixImbibition_m3s[l]);
          capillarity_matrix_step_m3 += tsubstep*saturated_matrix_correction_m3s[l];
        }
      }
      
      // 3. Sum source/sinks
      for(int l=0;l<nlayers;l++) finalSourceSinks_m3s[l] = source_sink_def_m3s[l] + matrixImbibition_m3s[l] + saturated_matrix_correction_m3s[l];
      
      // 4. Celia, Bouloutas & Zarba (1990) Mixed (theta-psi) Modified Picard Iteration
      for(int pic = 0; pic < max_picard; pic++) {
        total_nsubsteps++;
        for(int l=0;l<nlayers;l++) {
          if(soilDomains=="single") {
            evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Ksat[l], waterFluidity[l],
                               Psi_step[l], false, 0.0, Psi_b, theta_micro_step[l], C_step[l], K_step[l]);
          } else {
            evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Ksat[l], waterFluidity[l],
                               Psi_step[l], true, Ksat_b[l], Psi_b, theta_micro_step[l], C_step[l], K_step[l]);
          }
          K_step_ms[l] = 0.01*K_step[l]/(86400.0*cmdTOmmolm2sMPa);
          C_step_m[l] = std::max(C_step[l]*mTOMPa, 1.0e-7);
        }
        
        for(int l=0;l<nlayers;l++) {
          double storage_coef = (lambda[l]*dZ_m[l])/tsubstep;
          double c_term = storage_coef*C_step_m[l];
          double mass_defect = storage_coef*(theta_micro_step[l] - theta_micro[l]);
          
          if(l==0) { //first layer
            K_up = 0.0;
            K_down = K_step_ms[l];
            drain_below[l] = K_down;
            if(prop_saturated[l]==1.0) {
              drain_below[l] = 0.0;
            }
            a[l] = 0.0;
            c[l] = (nlayers > 1) ? (-1.0*K_down/dZDown[l]) : 0.0;
            double bc_cond = (!freeDrainage && num_saturated < nlayers && nlayers == 1) ? (K_down/dZDown[l]) : 0.0;
            b[l] = c_term - c[l] + bc_cond;
            d[l] = c_term*Psi_step_m[l] - mass_defect - drain_below[l] + bc_cond*(Psi_bc/mTOMPa) + finalSourceSinks_m3s[l];
          } else if(l<(nlayers - 1)) {
            K_up = K_step_ms[l-1];
            K_down = K_step_ms[l];
            drain_above[l] = K_up;
            drain_below[l] = K_down;
            if(prop_saturated[l]==1.0) {
              drain_above[l] = 0.0;
              drain_below[l] = 0.0;
            }
            a[l] = -1.0*K_up/dZUp[l];
            c[l] = -1.0*K_down/dZDown[l];
            b[l] = c_term - a[l] - c[l];
            d[l] = c_term*Psi_step_m[l] - mass_defect + drain_above[l] - drain_below[l] + finalSourceSinks_m3s[l];
          } else { // last layer
            K_up = K_step_ms[l-1];
            K_down = K_step_ms[l];
            drain_below[l] = K_down;
            drain_above[l] = K_up;
            if(prop_saturated[l]==1.0) {
              drain_above[l] = 0.0;
              drain_below[l] = 0.0;
            }
            a[l] = -1.0*K_up/dZUp[l];
            c[l] = 0.0;
            double bc_cond = (!freeDrainage && num_saturated < nlayers) ? (K_down/dZDown[l]) : 0.0;
            b[l] = c_term - a[l] + bc_cond;
            d[l] = c_term*Psi_step_m[l] - mass_defect + drain_above[l] - drain_below[l] + bc_cond*(Psi_bc/mTOMPa) + finalSourceSinks_m3s[l];
          }
        }
        //TRIDIAGONAL SOLVING
        tridiagonalSolving_c(a,b,c,d,e,f, Psi_step_t1);
        
        // Conservative Celia moisture update & convergence check
        double max_dtheta = 0.0;
        for(int l=0;l<nlayers;l++) {
          double psi_new_m = Psi_step_t1[l];
          double dpsi_m = psi_new_m - Psi_step_m[l];
          double dtheta = C_step_m[l]*dpsi_m;
          if(std::abs(dtheta) > max_dtheta) max_dtheta = std::abs(dtheta);
          theta_micro_step[l] += dtheta;
          Psi_step_m[l] = psi_new_m;
          Psi_step[l] = psi_new_m*mTOMPa;
        }
        if(max_dtheta < picard_theta_tol) {
          break;
        }
      }
      
      //calculate drainage and bottom boundary capillary flux (m3) (Fix B7)
      if(freeDrainage) {
        drainage_matrix_step_m3 += drain_below[nlayers -1]*tsubstep;
      } else {
        if(num_saturated < nlayers) {
          drainage_matrix_step_m3 += drain_below[nlayers - num_saturated -1]*tsubstep;
          double cap_bot_m3s = -1.0*(K_step_ms[nlayers - 1]/dZDown[nlayers - 1])*(Psi_step_m[nlayers - 1] - (Psi_bc/mTOMPa));
          capillarity_matrix_step_m3 += cap_bot_m3s*tsubstep;
        }
      }
      
      //Update Psi and theta (upward oversaturation cascade)
      double res_mm = 0.0;
      for(int l=(nlayers-1);l>=0;l--) {
        double new_theta = theta_micro_step[l];
        double quasi_sat_theta = quasi_sat_theta_vec[l];
        // Add previous residue
        if(res_mm > 0.0) {
          new_theta += res_mm/(widths[l]*lambda[l]);
          res_mm = 0.0;
        }
        if(soilDomains=="single") {
          if(new_theta > quasi_sat_theta) {
            res_mm = std::abs(new_theta - quasi_sat_theta)*widths[l]*lambda[l];
            new_theta = quasi_sat_theta;
            Psi_step[l] = Psi_quasi_sat;
          } else {
            Psi_step[l] = theta2psiVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], new_theta);
          }
          theta_micro_step[l] = new_theta;
          Theta[l] = new_theta;
        } else {
          theta_micro_step[l] = new_theta;
          //If needed, add exceeding moisture to the macropores and correct Psi (Fix B1 & B9)
          double excess_theta_step = 0.0, excess_to_macro = 0.0;
          if(theta_micro_step[l] > theta_b[l]) {
            double C_macro_step = std::max(0.0, 1.0 - S_macro_step[l])*e_macro[l];
            excess_theta_step = theta_micro_step[l] - theta_b[l];
            excess_to_macro = std::min(excess_theta_step, C_macro_step);
            res_mm += (excess_theta_step - excess_to_macro)*widths[l]*lambda[l];
            theta_micro_step[l] = theta_b[l];
            Psi_step[l] = Psi_b;
          } else {
            Psi_step[l] = theta2psiVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], theta_micro_step[l]);
          }
          if(excess_to_macro>0.0){
            lateral_flows_step_mm[l] -= excess_to_macro*widths[l]*lambda[l]; //negative flow (in mm/step)
            matrixExcess_m3s[l] = excess_to_macro*widths[l]*lambda[l]*0.001/tsubstep; //Source of water flowing into macropores (m3/s)
          } else {
            matrixExcess_m3s[l] = 0.0;
          }
        }
        double th_dummy;
        if(soilDomains=="single") {
          evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat[l], Ksat[l], waterFluidity[l],
                             Psi_step[l], false, 0.0, Psi_b, th_dummy, C_step[l], K_step[l]);
        } else {
          evalVanGenuchten_c(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Ksat[l], waterFluidity[l],
                             Psi_step[l], true, Ksat_b[l], Psi_b, th_dummy, C_step[l], K_step[l]);
        }
        Psi_step_m[l]= Psi_step[l]/mTOMPa; // MPa to m
        K_step_ms[l] = 0.01*K_step[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
        C_step_m[l] = C_step[l]*mTOMPa; //From MPa-1 to m-1
      }
      //Generate saturation excess if there was some residue in the top layer
      saturation_excess_matrix_step_mm += res_mm;
      
      if(soilDomains=="dual") {
        //Solve macropore domain (Fix B3, B4, B10)
        double infiltration_macropores_substep_mm = 0.0;
        for(int l=0;l<nlayers;l++) {
          //Set imbibition as sink for macropore
          double sourceSink_macro_m3s = matrixExcess_m3s[l] - matrixImbibition_m3s[l];
          
          //Infiltration into macropores (Fix B3: full step target attempted)
          if(l==0) {
            if(infiltration_remaining_macropores_step_mm > 0.0) {
              infiltration_macropores_substep_mm = infiltration_remaining_macropores_step_mm;
              double infiltration_macropores_substep_m3s = (infiltration_macropores_substep_mm*mm_2_m3)/tsubstep;
              sourceSink_macro_m3s += infiltration_macropores_substep_m3s;
              infiltration_remaining_macropores_step_mm = 0.0;
            }
          }
          capillarity_macropores_step_m3 += tsubstep*saturated_macropore_correction_m3s[l];
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
          double S_t1 = rootFindingMacropores_c(S_macro_step[l], K_up, ksat_i, ksat_b_i, kin_exp,
                                                e_macro[l], lambda[l], dZ_m[l], sourceSink_macro_m3s, tsubstep);
          Kmacro_step_ms[l] = (S_t1 > 0.0) ? ((ksat_i - ksat_b_i)*pow(S_t1, kin_exp)) : 0.0;
          S_macro_step[l] = S_t1;
        }
        //Drainage
        if(freeDrainage) {
          drainage_macropores_step_m3 += Kmacro_step_ms[nlayers-1]*tsubstep;
        } else {
          if(num_saturated<nlayers) {
            double flow = Kmacro_step_ms[nlayers - num_saturated -1];
            drainage_macropores_step_m3 += flow*tsubstep;
          }
        }
        //Update theta_macro for the next step
        double res_mm = 0.0;
        for(int l=(nlayers-1);l>=0;l--) {
          theta_macro_step[l] = S_macro_step[l]*e_macro[l] + res_mm/(widths[l]*lambda[l]);
          // If oversaturation of macroporosity occurs
          if(theta_macro_step[l] > e_macro[l]) {
            res_mm = widths[l]*(theta_macro_step[l] - e_macro[l])*lambda[l]; //residue in mm for layers above
            theta_macro_step[l] = e_macro[l];
            S_macro_step[l] = 1.0; //Correct S_macro to saturation
            Kmacro_step_ms[l] = (prop_saturated[l] == 1.0) ? 0.0 : (Ksat_ms[l] - Ksat_b_ms[l]);
            if(l==0) { // If there has been infiltration, add it back to infiltration target
              double min_diff = std::min(res_mm, infiltration_macropores_substep_mm);
              infiltration_macropores_substep_mm -= min_diff;
              infiltration_remaining_macropores_step_mm += min_diff;
              res_mm -= min_diff;
            }
          } else {
            // Fix B10: Keep S_macro_step[l] and Kmacro_step_ms[l] synchronized when residue is absorbed
            if(res_mm > 0.0) {
              S_macro_step[l] = theta_macro_step[l]/e_macro[l];
              Kmacro_step_ms[l] = (prop_saturated[l] == 1.0 || S_macro_step[l] <= 0.0)
                ? 0.0
              : ((Ksat_ms[l] - Ksat_b_ms[l])*pow(S_macro_step[l], kin_exp));
            }
            res_mm = 0.0;
          }
        }
        //Generate saturation excess if there was some residue in the top layer
        saturation_excess_macropores_step_mm += res_mm;
        // Store infiltration
        infiltration_macropores_step_mm += infiltration_macropores_substep_mm;
      }
      
      //Update theta, W and final volumes
      Vfin_mm = 0.0;
      Vfin_micro_mm = 0.0;
      Vfin_macro_mm = 0.0;
      for(int l=0;l<nlayers;l++) {
        if(soilDomains=="single") {
          Theta[l] = theta_micro_step[l];
        } else {
          Theta[l] = theta_micro_step[l] + theta_macro_step[l];
          Vfin_micro_mm += theta_micro_step[l]*widths[l]*lambda[l];
          Vfin_macro_mm += theta_macro_step[l]*widths[l]*lambda[l];
        }
        Vfin_mm += Theta[l]*widths[l]*lambda[l];
      }
      
      if(soilDomains=="single") {
        double balance_step_mm = m3_2_mm*(tstep*sum_source_sink_def_m3s + capillarity_matrix_step_m3 - drainage_matrix_step_m3) - saturation_excess_matrix_step_mm;
        matrix_correction_step_mm = balance_step_mm + Vini_step_mm - Vfin_mm;
      } else {
        infiltration_excess_macropores_step_mm = infiltration_remaining_macropores_step_mm;
        double sum_lateral_flows_step_mm = 0.0;
        for(int l=0;l<nlayers;l++) sum_lateral_flows_step_mm += lateral_flows_step_mm[l];
        double balance_micro_step_mm = m3_2_mm*(tstep*sum_source_sink_def_m3s + capillarity_matrix_step_m3 - drainage_matrix_step_m3) + sum_lateral_flows_step_mm - saturation_excess_matrix_step_mm;
        double balance_macro_step_mm = infiltration_macropores_step_mm + m3_2_mm*(capillarity_macropores_step_m3 - drainage_macropores_step_m3) - sum_lateral_flows_step_mm - saturation_excess_macropores_step_mm;
        matrix_correction_step_mm = balance_micro_step_mm + Vini_step_micro_mm - Vfin_micro_mm;
        macropore_correction_step_mm = balance_macro_step_mm + Vini_step_macro_mm - Vfin_macro_mm;
      }
      
      // Fix B11: Do NOT dump matrix_correction_step_mm or macropore_correction_step_mm
      // into drainage or capillarity (Celia's mixed form is intrinsically mass-conservative).
      
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
        if(soilDomains=="dual") {
          S_macro[l] = S_macro_step[l];
          Kmacro_ms[l] = Kmacro_step_ms[l];
          theta_macro[l] = theta_macro_step[l];
        }
        theta_micro[l] = theta_micro_step[l];
      }
      //Prepare for next step
      Vini_step_macro_mm = Vfin_macro_mm;
      Vini_step_micro_mm = Vfin_micro_mm;
      Vini_step_mm = Vfin_mm;
    }
    
    double drainage_matrix_mm = drainage_matrix_m3*1000.0; //m3/m2 to mm/m2
    double capillarity_matrix_mm = capillarity_matrix_m3*1000.0; //m3/m2 to mm/m2
    double capillarity_macropores_mm = capillarity_macropores_m3*1000.0; //m3/m2 to mm/m2
    
    
    // This effectively changes soil moisture states in the soil object
    for(int l=0; l< nlayers;l++) {
      soil.setW(l,Theta[l]/soil.getThetaFC(l));
    }
    
    //Correct overestimation of capillarity by using deep drainage
    double dif_mm = std::min(capillarity_matrix_mm,  drainage_matrix_mm);
    capillarity_matrix_mm -= dif_mm;
    drainage_matrix_mm  -= dif_mm;
    //Correct overestimation of capillarity by using saturation excess
    dif_mm = std::min(capillarity_matrix_mm,  saturation_excess_matrix_mm);
    capillarity_matrix_mm -= dif_mm;
    saturation_excess_matrix_mm  -= dif_mm;
    
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
      SWBres.matrixMacroporeFlow_mm = std::accumulate(matrix_macropore_flows_mm.begin(), matrix_macropore_flows_mm.end(), 0.0);
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
      SWBres.capillarityMatrix_mm = capillarity_matrix_mm;
      SWBres.capillarityMacropores_mm = capillarity_macropores_mm;
      SWBres.capillarityRise_mm = capillarity_matrix_mm + capillarity_macropores_mm;
      SWBres.correctionMatrix_mm = matrix_correction_mm;
      SWBres.correctionMacropores_mm = macropore_correction_mm;
      SWBres.correction_mm = matrix_correction_mm + macropore_correction_mm;
      SWBres.matrixVolumeChange_mm = Vfin_micro_mm - Vini0_micro_mm;
      SWBres.macroporeVolumeChange_mm = Vfin_macro_mm - Vini0_macro_mm;
      SWBres.volumeChange_mm = Vfin_micro_mm - Vini0_micro_mm  + Vfin_macro_mm - Vini0_macro_mm;
      SWBres.flowBalance_mm = (SWBres.infiltration_mm + SWBres.capillarityRise_mm + SWBres.localSourceSinks_mm + SWBres.lateralSourceSinks_mm) - (SWBres.saturationExcess_mm + SWBres.deepDrainage_mm);
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
      SWBres.flowBalance_mm = (SWBres.infiltration_mm + SWBres.capillarityRise_mm + SWBres.localSourceSinks_mm + SWBres.lateralSourceSinks_mm) - (SWBres.saturationExcess_mm + SWBres.deepDrainage_mm);
      SWBres.substeps = total_nsubsteps;
    }
  }
}