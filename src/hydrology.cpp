// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include "root.h"
#include "hydraulics.h"
#include "communication_structures.h"
#include "biophysicsutils.h"
#include "numerical_solving.h"
#include <meteoland.h>
using namespace Rcpp;

//Old defaults
//ERconv=0.05, ERsyn = 0.2
//New defaults
//Rconv = 5.6, Rsyn = 1.5

//' @param month Month of the year (from 1 to 12).
//' @param prec Precipitation for a given day (mm).
//' @param rainfallIntensityPerMonth A vector with twelve positions with average intensity of rainfall (in mm/h) for each month.
//' 
//' @rdname hydrology_interception
//' @keywords internal
// [[Rcpp::export("hydrology_rainfallIntensity")]]
double rainfallIntensity(int month, double prec, NumericVector rainfallIntensityPerMonth){
  double Ri_month = rainfallIntensityPerMonth[month - 1];
  double Ri = std::max(prec/24.0,Ri_month);
  return(Ri);
}

// [[Rcpp::export(".hydrology_interceptionGashDay")]]
double interceptionGashDay(double Rainfall, double Cm, double p, double ER=0.05) {
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
double interceptionLiuDay(double Rainfall, double Cm, double p, double ER=0.05){
  double I = Cm*(1.0 - exp(-1.0*(Rainfall)*((1.0 - p)/Cm)))*(1.0 - (ER/(1.0 - p))) + (ER*Rainfall);
  return(I);
}

//' @rdname hydrology_soilEvaporation
//' 
//' @param DEF Water deficit in the (topsoil) layer.
//' @param PETs Potential evapotranspiration at the soil surface.
//' @param Gsoil Gamma parameter (maximum daily evaporation).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_soilEvaporationAmount")]]
double soilEvaporationAmount(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t + 1.0)-sqrt(t)), PETs);
  return(Esoil);
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
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector widths = soil["widths"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector psiSoil = psi(soil, soilFunctions);
  double Esoil = 0.0;
  if(snowpack == 0.0) {
    double PETsoil = pet*(LgroundSWR/100.0);
    double Gsoil = 0.5; //TO DO, implement pedotransfer functions for Gsoil
    // Allow evaporation only if water potential is higher than -2 MPa
    if(psiSoil[0] > -2.0) Esoil = soilEvaporationAmount((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
    if(modifySoil){
      W[0] = W[0] - (Esoil/Water_FC[0]);
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
  if(NumericVector::is_na(herbLAI)) return(0.0);
  double Tmax_herb = pet*(LherbSWR/100.0)*(0.134*herbLAI - 0.006*pow(herbLAI, 2.0));
  NumericVector widths = soil["widths"];
  NumericVector W = soil["W"];
  int nlayers = widths.size();
  NumericVector psiSoil = psi(soil, soilFunctions);
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector EherbVec(nlayers,0.0);
  NumericVector V = ldrRS_one(50, 500, NA_REAL, widths);
  for(int l=0;l<nlayers;l++) {
    EherbVec[l] = V[l]*Tmax_herb*Psi2K(psiSoil[0], -1.5, 2.0); 
    if(modifySoil) {
      W[l] = W[l] - (EherbVec[l]/Water_FC[l]);
    }
  }
  return(EherbVec);
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
//' @param input A numeric vector of (daily) water input (in mm of water).
//' @param Ssoil Soil water storage capacity (can be referred to topsoil) (in mm of water).
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
// [[Rcpp::export("hydrology_infiltrationBoughton")]]
double infiltrationBoughton(double input, double Ssoil) {
  double I = 0;
  if(input>0.2*Ssoil) {
    I = input-(pow(input-0.2*Ssoil,2.0)/(input+0.8*Ssoil));
  } else {
    I = input;
  }
  return(I);
}

double fGreenAmpt(double x, double t, double psi_w, double Ksat, double delta_theta) {
  double f = Ksat*t + std::abs(psi_w)*delta_theta*log(1.0 + (x/(std::abs(psi_w)*delta_theta))) - x;
  return(f);
}
double fGreenAmptDer(double x, double t, double psi_w, double Ksat, double delta_theta) {
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
double infitrationGreenAmpt(double t, double psi_w, double Ksat, double theta_sat, double theta_dry) {
  double delta_theta = theta_sat - theta_dry;
  double x,x1,e,fx,fx1;
  x1 = 0.0;//initial guess
  e = 0.001; // accuracy in mm
  int cnt = 0;
  int mxiter = 100;
  do {
    x=x1; /*make x equal to the last calculated value of  x1*/
    fx=fGreenAmpt(x, t, psi_w, Ksat, delta_theta);            //simplifying f(x)to fx
    fx1=fGreenAmptDer(x, t, psi_w, Ksat, delta_theta);            //simplifying fprime(x) to fx1
    x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
    cnt++;
  } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
  return(x);
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
  int nlayers = widths.length();
  NumericVector Pvec = NumericVector(nlayers,0.0);
  NumericVector Ivec = NumericVector(nlayers,0.0);
  double z1 = 0.0;
  double p1 = 1.0;
  for(int i=0;i<nlayers;i++) {
    double ai = a*pow(1.0-macro[i],b);
    if(i<(nlayers-1)) {
      Pvec[i] = p1*(1.0-exp(ai*widths[i]));
    } else {
      Pvec[i] = p1;
    }
    p1 = p1*exp(ai*widths[i]);
    z1 = z1 + widths[i];
    Ivec[i] = I*Pvec[i];
  }
  return(Ivec);
}


//' @rdname hydrology_infiltration
//' 
//' @param rainfallInput Water from the rainfall event reaching the soil surface (mm)
//' @param soil A list containing the description of the soil (see \code{\link{soil}}).
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param rainfallIntensity rainfall intensity rate (mm/h)
//' @param model Infiltration model, either "GreenAmpt1911" or "Boughton1989"
//' @param K_correction Correction for saturated conductivity, to account for increased infiltration due to macropore presence
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationAmount")]]
double infiltrationAmount(double rainfallInput, double rainfallIntensity, DataFrame soil, 
                          String soilFunctions, String model = "GreenAmpt1911", double K_correction = 1.0) {
  double infiltration = 0.0;
  if(model=="GreenAmpt1911") {
    NumericVector clay = soil["clay"];
    NumericVector sand = soil["sand"];
    NumericVector bd = soil["bd"];
    NumericVector Ksat = soil["Ksat"];
    String usda = USDAType(clay[0], sand[0]);
    NumericVector cp = campbellParamsClappHornberger(usda);
    NumericVector theta_dry = theta(soil, soilFunctions);
    double t = std::min(24.0, rainfallInput/rainfallIntensity); // time in hours
    double b = cp["b"];
    double psi_w = cp["psi_sat_cm"]*((2.0*b + 3.0)/(2*b + 6.0));
    double theta_sat = cp["theta_sat"];
    double K_sat_0 = K_correction*Ksat[0]/(24.0*cmdTOmmolm2sMPa); // from mmolH20*m-2*MPa-1*s-1 to cm_h
    infiltration = infitrationGreenAmpt(t, psi_w, K_sat_0, theta_sat, theta_dry[0]);
  } else if(model=="Boughton1989") {
    NumericVector Water_FC = waterFC(soil, soilFunctions);
    infiltration = infiltrationBoughton(rainfallInput, Water_FC[0]);
  } else {
    stop("Wrong infiltration model!");
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
double snowMelt(double tday, double rad, double LgroundSWR, double elevation) {
  if(NumericVector::is_na(rad)) stop("Missing radiation data for snow melt!");
  if(NumericVector::is_na(elevation)) stop("Missing elevation data for snow melt!");
  double rho = meteoland::utils_airDensity(tday, meteoland::utils_atmosphericPressure(elevation));
  double ten = (86400.0*tday*rho*1013.86*1e-6/100.0); //ten can be negative if temperature is below zero
  double ren = (rad*(LgroundSWR/100.0))*(0.1); //90% albedo of snow
  double melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
  return(melt);
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
  //Soil input
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  String interceptionMode = control["interceptionMode"];
  double swe = x["snowpack"]; //snow pack
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
    melt = std::min(swe, snowMelt(tday, rad, LgroundSWR, elevation));
    // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
    swe = swe-melt;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Interception = 0.0;
  if(rain>0.0)  {
    if(interceptionMode=="Gash1995") {
      Interception = interceptionGashDay(rain,Cm,LgroundPAR/100.0,er);
    } else if(interceptionMode =="Liu2001") {
      Interception = interceptionLiuDay(rain,Cm,LgroundPAR/100.0,er);
    } else {
      stop("Wrong interception model!");
    }
    NetRain = rain - Interception; 
  }
  if(modifyInput) {
    x["snowpack"] = swe;
  }
  NumericVector WI = NumericVector::create(_["Rain"] = rain, _["Snow"] = snow,
                                           _["Interception"] = Interception,
                                           _["NetRain"] = NetRain, 
                                           _["Snowmelt"] = melt);
  return(WI);
}


//Imbitition from macropores to micropores following Larsbo et al. 2005, eq. 6-7
//Diffusion equation with gradients in water content as driving force
//Returns m3/m3/s
double microporeImbibitionRate(double theta_b, double theta_micro, 
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

double rootFindingMacropores(double S_t, double K_up, double Ksat_ms, double Ksat_b_ms, double kin_exp,
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
    if(b>10.0) stop("Could not find appropriate bounds for macropore circulation");
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
    stop("Not found");
  }
  return(c);
}


NumericVector soilWaterBalance_inner(List SWBcommunication, DataFrame soil, String soilFunctions, 
                                     double rainfallInput, double rainfallIntensity, double snowmelt, NumericVector sourceSink, 
                                     double runon = 0.0, Nullable<NumericVector> lateralFlows = R_NilValue, double waterTableDepth = NA_REAL,
                                     String infiltrationMode = "GreenAmpt1911", double infiltrationCorrection = 5.0, String soilDomains = "buckets", 
                                     int nsteps = 24, int max_nsubsteps = 3600, bool modifySoil = true) {
  
  if((soilDomains!="single") && (soilDomains!="dual")  && (soilDomains!="buckets") ) stop("Unrecognized soilDomain value");
  if(!modifySoil) soil = clone(soil);
  
  NumericVector widths = soil["widths"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector W = soil["W"];
  int nlayers = soil.nrow();
  
  bool micropore_imbibition = true;
  
  double mm_2_m3 = 0.001;
  double m3_2_mm = 1.0/mm_2_m3;
  double mm_day_2_m3_s = mm_2_m3*(1.0/86400.0);//From mm/day = l/m2/day = dm3/day to m3/m2/s
  
  //Infiltration
  double K_correction = 1.0;
  if(soilDomains!="dual") K_correction = infiltrationCorrection;
  double infiltration_matrix_mm = infiltrationAmount(rainfallInput, rainfallIntensity, soil, 
                                                     soilFunctions, infiltrationMode, K_correction);
  double infiltration_macropores_mm = 0.0;
  double infiltration_excess_matrix_mm = rainfallInput - infiltration_matrix_mm;
  double infiltration_target_macropores_mm = infiltration_excess_matrix_mm;
  double infiltration_excess_macropores_mm = 0.0;
  double matrix_correction_mm = 0.0;
  double macropore_correction_mm = 0.0;
  
  
  
  //Add snow-melt and runon to infiltration_matrix_mm
  double snowmelt_mm = 0.0;
  if(!NumericVector::is_na(snowmelt)) {
    snowmelt_mm = snowmelt;
    infiltration_matrix_mm += snowmelt_mm; 
  }
  double runon_mm = 0.0;
  if(!NumericVector::is_na(runon)) {
    runon_mm = runon;
    infiltration_matrix_mm += runon_mm; 
  }
  
  //Copy sinks
  NumericVector source_sink_def_mm(nlayers);
  for(int l=0;l<nlayers;l++) {
    source_sink_def_mm[l] = sourceSink[l];
  }
  
  //Add infiltration to matrix def source/sinks
  if(soilDomains!="dual") {
    NumericVector IVec = infiltrationRepartition(infiltration_matrix_mm, widths, macro);
    for(int l=0;l<nlayers;l++) {
      source_sink_def_mm[l] += IVec[l];
    }
  } else {
    source_sink_def_mm[0] += infiltration_matrix_mm;
  }
  //Add lateral flows to matrix def source/sinks
  NumericVector lateralFlows_mm(0.0, nlayers);
  if(lateralFlows.isNotNull()) {
    lateralFlows_mm = NumericVector(lateralFlows);
    for(int l=0;l<nlayers;l++) {
      source_sink_def_mm[l] += lateralFlows_mm[l];
    }
  }
  
  NumericVector res;
  if(soilDomains == "buckets") {
    NumericVector Ksat = soil["Ksat"];
    NumericVector W = clone(Rcpp::as<Rcpp::NumericVector>(soil["W"])); //Access to soil state variable
    NumericVector Water_FC = waterFC(soil, soilFunctions);
    NumericVector Water_SAT = waterSAT(soil, soilFunctions);
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
    if(modifySoil) {
      NumericVector Ws = soil["W"];
      for(int l=0;l<nlayers;l++) Ws[l] = W[l];
    }
    res = NumericVector::create(_["Local source/sinks"] = sum(sourceSink),
                                _["Lateral source/sinks"] = sum(lateralFlows_mm),
                                _["Infiltration"] = infiltration_matrix_mm,
                                _["InfiltrationExcess"] = infiltration_excess_matrix_mm,
                                _["SaturationExcess"] = saturation_excess_matrix_mm,
                                _["Runoff"] = saturation_excess_matrix_mm + infiltration_excess_matrix_mm,
                                _["DeepDrainage"] = drainage_matrix_mm,
                                _["CapillarityRise"] = 0.0);
  } else if(soilDomains == "dual" || soilDomains == "single") {
    //Initialize matrix-macrore flows (positive in the direction of matrix)
    NumericVector matrix_macropore_flows_mm(nlayers, 0.0);
    
    //Set time steps
    double tstep = 86400.0/((double) nsteps);
    double tsubstep = tstep;
    double halftsubstep = tsubstep/2.0;
    double rainfallIntensity_step = rainfallIntensity*24.0/((double) nsteps); //mm/step
    
    NumericVector source_sink_def_m3s =  source_sink_def_mm*mm_day_2_m3_s;
    NumericVector matrixImbibition_m3s(nlayers, 0.0);
    NumericVector matrixExcess_m3s(nlayers, 0.0);
    NumericVector saturated_matrix_correction_m3s(nlayers, 0.0);
    NumericVector saturated_macropore_correction_m3s(nlayers, 0.0);
    
    
    NumericVector dZ_m = SWBcommunication[SOILWBCOM_dZ_m];
    for(int l=0;l<nlayers;l++) dZ_m[l] = widths[l]*0.001; //mm to m
    
    NumericVector dZUp = SWBcommunication[SOILWBCOM_dZUp];
    NumericVector dZDown = SWBcommunication[SOILWBCOM_dZDown];
    NumericVector lambda = SWBcommunication[SOILWBCOM_lambda];
    
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
    
    NumericVector prop_saturated(nlayers, 0.0);
    int num_saturated = 0;
    double freeDrainage = true;
    if(!NumericVector::is_na(waterTableDepth)) {
      freeDrainage = (waterTableDepth > sum(widths));
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
    NumericVector Ksat_ori = soil["Ksat"];
    NumericVector n =soil["VG_n"];
    NumericVector alpha = soil["VG_alpha"];
    NumericVector theta_res = soil["VG_theta_res"];
    NumericVector theta_sat = soil["VG_theta_sat"];
    NumericVector Tsoil = soil["Temp"];
    
    //Estimate Theta, Psi, C, K
    NumericVector Theta = theta(soil, "VG");
    NumericVector Theta_FC = thetaFC(soil, "VG");
    
    //Microporosity or single domain
    NumericVector theta_micro = SWBcommunication[SOILWBCOM_theta_micro];
    NumericVector theta_b = SWBcommunication[SOILWBCOM_theta_b];
    NumericVector theta_macro = SWBcommunication[SOILWBCOM_theta_macro];
    NumericVector theta_sat_fict = SWBcommunication[SOILWBCOM_theta_sat_fict];
    NumericVector Ksat_b = SWBcommunication[SOILWBCOM_Ksat_b];
    NumericVector Ksat_b_ms = SWBcommunication[SOILWBCOM_Ksat_b_ms];
    NumericVector Ksat = SWBcommunication[SOILWBCOM_Ksat];
    NumericVector Ksat_ms = SWBcommunication[SOILWBCOM_Ksat_ms];
    NumericVector Psi = SWBcommunication[SOILWBCOM_Psi];
    NumericVector K = SWBcommunication[SOILWBCOM_K];
    NumericVector C = SWBcommunication[SOILWBCOM_C];
    NumericVector Psi_m = SWBcommunication[SOILWBCOM_Psi_m];
    NumericVector K_ms = SWBcommunication[SOILWBCOM_K_ms];
    NumericVector Kbc = SWBcommunication[SOILWBCOM_Kbc];
    NumericVector Kbc_ms = SWBcommunication[SOILWBCOM_Kbc_ms];
    NumericVector C_m = SWBcommunication[SOILWBCOM_C_m];
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
    NumericVector S_macro, e_macro, Kmacro_ms;
    if(soilDomains=="dual") {
      S_macro = SWBcommunication[SOILWBCOM_S_macro];
      e_macro = SWBcommunication[SOILWBCOM_e_macro];
      Kmacro_ms = SWBcommunication[SOILWBCOM_Kmacro_ms];
      for(int l=0;l<nlayers;l++) {
        S_macro[l] = 0.0;
        e_macro[l] = 0.0;
        Kmacro_ms[l] = 0.0;
      }
    }
    
    NumericVector waterFluidity = SWBcommunication[SOILWBCOM_waterFluidity];
    for(int l=0;l<nlayers;l++) {
      waterFluidity[l] = 1.0;
      if(!NumericVector::is_na(Tsoil[l])) {
        if(Tsoil[l]>0) {
          waterFluidity[l] = 1.0/waterDynamicViscosity(Tsoil[l]); 
        } else {
          waterFluidity[l] = 0.0;
        }
      }
      Ksat[l] = Ksat_ori[l]*lambda[l];//Multiply K for the space available for water movement
      Ksat_ms[l] = 0.01*waterFluidity[l]*Ksat[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      if(soilDomains=="single") {
        Psi[l] = theta2psiVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]);
        C[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
        K[l] = waterFluidity[l]*psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
        Kbc[l] = waterFluidity[l]*psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_bc);
      } else {
        theta_sat_fict[l] = theta_sat[l] - macro[l];
        //Matching theta point and theta partitioning
        //This ensures theta_b < theta_sat - macro & e_macro > macro
        theta_b[l] = psi2thetaVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat_fict[l], Psi_b);
        e_macro[l] = theta_sat[l] - theta_b[l];
        // Rcout<<e_macro[l]<< " "<< macro[l]<<"\n";
        Ksat_b[l] = psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_b);
        theta_micro[l] = std::min(Theta[l], theta_b[l]);
        theta_macro[l] = std::max(Theta[l] - theta_b[l], 0.0);
        //Water potential, conductivity and capacitance in the micropore domain according to the modified retention
        Psi[l] = theta2psiVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat_fict[l], theta_micro[l]);
        C[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi[l]);
        K[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
                                                            Psi[l], Psi_b);
        Kbc[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
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
    
    NumericVector a = SWBcommunication[SOILWBCOM_a];
    NumericVector b = SWBcommunication[SOILWBCOM_b];
    NumericVector c = SWBcommunication[SOILWBCOM_c];
    NumericVector d = SWBcommunication[SOILWBCOM_d];
    NumericVector e = SWBcommunication[SOILWBCOM_e];
    NumericVector f = SWBcommunication[SOILWBCOM_f];
    NumericVector Psi_step_t1(nlayers), Psi_step_t05(nlayers);
    
    NumericVector K_step_ms05 = SWBcommunication[SOILWBCOM_K_step_ms05];
    NumericVector C_step_m05 = SWBcommunication[SOILWBCOM_C_step_m05];
    NumericVector C_step = SWBcommunication[SOILWBCOM_C_step];
    NumericVector C_step_m = SWBcommunication[SOILWBCOM_C_step_m];
    NumericVector K_step_ms = SWBcommunication[SOILWBCOM_K_step_ms];
    NumericVector K_step = SWBcommunication[SOILWBCOM_K_step];
    NumericVector Psi_step = SWBcommunication[SOILWBCOM_Psi_step];
    NumericVector Psi_step_m = SWBcommunication[SOILWBCOM_Psi_step_m];
    NumericVector S_macro_step = SWBcommunication[SOILWBCOM_S_macro_step];
    NumericVector Kmacro_step_ms = SWBcommunication[SOILWBCOM_Kmacro_step_ms];
    NumericVector theta_macro_step = SWBcommunication[SOILWBCOM_theta_macro_step];
    NumericVector theta_micro_step = SWBcommunication[SOILWBCOM_theta_micro_step];
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
    NumericVector finalSourceSinks_m3s = SWBcommunication[SOILWBCOM_finalSourceSinks_m3s];
    NumericVector capill_below = SWBcommunication[SOILWBCOM_capill_below];
    NumericVector drain_above = SWBcommunication[SOILWBCOM_drain_above];
    NumericVector drain_below = SWBcommunication[SOILWBCOM_drain_below];
    NumericVector lateral_flows_step_mm = SWBcommunication[SOILWBCOM_lateral_flows_step_mm];
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
                theta_l = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
                quasi_sat_theta_l = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_quasi_sat);
                saturated_matrix_correction_m3s[l] = std::max(0.0, ((prop_saturated[l]*quasi_sat_theta_l - theta_l)*lambda[l]*widths[l]*0.001)/tsubstep);
              } else {
                theta_l = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
                quasi_sat_theta_l = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_quasi_sat);
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
                double Ksat_fict_ms = psi2kVanGenuchtenMicropores(Ksat_b_ms[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
                                                                  0.0, Psi_b);
                double D_theta_b_m2s = psi2DVanGenuchten(Ksat_fict_ms, n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
                                                         Psi_b);
                double D_theta_micro_m2s = psi2DVanGenuchten(Ksat_fict_ms, n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
                                                             Psi_step[l]);
                //Calculate imbibition rate in m3/m3/s = s-1
                double imbibitionRate = microporeImbibitionRate(theta_b[l], theta_micro_step[l], 
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
          tridiagonalSolving(a,b,c,d,e,f, Psi_step_t05);
          //MODIFY UNITS OF OUTPUT PSI
          for(int l=0;l<nlayers;l++) Psi_step_t05[l] = Psi_step_t05[l]*mTOMPa; // m to MPa
          
          //Calculate K and C at t05
          for(int l=0;l<nlayers;l++) {
            if(soilDomains=="single") {
              C_step_05 = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step_t05[l]);
              K_step_05 = waterFluidity[l]*psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step_t05[l]);
            } else {
              C_step_05 = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step_t05[l]);
              K_step_05 = waterFluidity[l]*psi2kVanGenuchtenMicropores(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
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
          tridiagonalSolving(a,b,c,d, e, f, Psi_step_t1);
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
              new_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
              quasi_sat_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_quasi_sat);
            } else {
              new_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              quasi_sat_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_quasi_sat);
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
                Psi_step[l] = theta2psiVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], new_theta);
              } else {
                Psi_step[l] = theta2psiVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l],new_theta);
              }
            }
            // Rcout<<" step "<<s<<" layer " <<l<< " final "<< Psi_step[l]<<"\n";
            //If dual model, update theta and manage excess to macropores
            if((soilDomains=="dual")) {
              //Update theta_micro for the next substep
              theta_micro_step[l] = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
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
              C_step[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
              K_step[l] = waterFluidity[l]*psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
            } else {
              C_step[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              K_step[l] = waterFluidity[l]*psi2kVanGenuchtenMicropores(Ksat_b[l], n[l], alpha[l], theta_res[l], theta_sat_fict[l], 
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
              double S_t1 = rootFindingMacropores(S_macro_step[l], K_up, ksat_i, ksat_b_i, kin_exp,
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
            Theta[l] = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
          } else {
            Theta[l] = theta_micro_step[l] + theta_macro_step[l];
            Vfin_micro_mm += theta_micro_step[l]*widths[l]*lambda[l];
            Vfin_macro_mm += theta_macro_step[l]*widths[l]*lambda[l];
          }
          Vfin_mm += Theta[l]*widths[l]*lambda[l];
        }
        
        double max_abs_correction_step_mm; 
        if(soilDomains=="single") {
          // Rcout << s << " "<< nsubsteps<< " sat_ex " << saturation_excess_matrix_mm << " dr " << drainage_matrix_step_m3*m3_2_mm << " cap " << capillarity_matrix_step_m3*m3_2_mm<<"\n";
          double balance_step_mm = m3_2_mm*(tstep*sum(source_sink_def_m3s)  + capillarity_matrix_step_m3 - drainage_matrix_step_m3) - saturation_excess_matrix_step_mm;
          matrix_correction_step_mm = balance_step_mm + Vini_step_mm - Vfin_mm;
          max_abs_correction_step_mm = std::abs(matrix_correction_step_mm);
          // Rcout << s << " "<< nsubsteps<<" ini "<< Vini_step_mm <<" fin " <<  Vfin_mm << " dif "<< Vfin_mm - Vini_step_mm
          //       << " dra " << m3_2_mm*drainage_matrix_step_m3 <<" cap " << m3_2_mm*capillarity_matrix_step_m3 <<"  bal "<< balance_step_mm<< " corr "<< matrix_correction_step_mm<<"\n";
        } else {
          //Add remaining target to infiltration excess
          infiltration_excess_macropores_step_mm = infiltration_remaining_macropores_step_mm;
          double sum_lateral_flows_step_mm = 0.0;
          for(int l=0;l<nlayers;l++) sum_lateral_flows_step_mm += lateral_flows_step_mm[l];
          double balance_micro_step_mm = m3_2_mm*(tstep*sum(source_sink_def_m3s) + capillarity_matrix_step_m3 - drainage_matrix_step_m3) + sum_lateral_flows_step_mm - saturation_excess_matrix_step_mm;
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
      W[l] = Theta[l]/Theta_FC[l];
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
      res = NumericVector::create(_["Local source/sinks"] = sum(sourceSink),
                                  _["Lateral source/sinks"] = sum(lateralFlows_mm),
                                  _["Matrix-macropore flow"] = sum(matrix_macropore_flows_mm),
                                  _["InfiltrationMatrix"] = infiltration_matrix_mm,
                                  _["InfiltrationMacropores"] = infiltration_macropores_mm,
                                  _["InfiltrationExcessMatrix"] = infiltration_target_macropores_mm,
                                  _["InfiltrationExcessMacropores"] = infiltration_excess_macropores_mm,
                                  _["SaturationExcessMatrix"] = saturation_excess_matrix_mm,
                                  _["SaturationExcessMacropores"] = saturation_excess_macropores_mm,
                                  _["DrainageMatrix"] = drainage_matrix_mm,
                                  _["DrainageMacropores"] = drainage_macropores_mm,
                                  _["CapillarityMatrix"] = capillarity_matrix_mm,
                                  _["CapillarityMacropores"] = capillarity_macropores_mm,
                                  _["CorrectionMatrix"] = matrix_correction_mm,
                                  _["CorrectionMacropores"] = macropore_correction_mm,
                                  _["MatrixVolumeChange"] = Vfin_micro_mm - Vini0_micro_mm,
                                  _["MacroporeVolumeChange"] = Vfin_macro_mm - Vini0_macro_mm,
                                  _["Infiltration"] = infiltration_matrix_mm + infiltration_macropores_mm,
                                  _["InfiltrationExcess"] = infiltration_target_macropores_mm + infiltration_excess_macropores_mm);
      res.push_back(saturation_excess_macropores_mm + saturation_excess_matrix_mm, "SaturationExcess");
      res.push_back(runoff_mm, "Runoff");
      res.push_back(drainage_macropores_mm + drainage_matrix_mm, "DeepDrainage");
      res.push_back(capillarity_matrix_mm + capillarity_macropores_mm, "CapillarityRise");
      res.push_back(matrix_correction_mm + macropore_correction_mm, "Correction");
      res.push_back(Vfin_micro_mm - Vini0_micro_mm  + Vfin_macro_mm - Vini0_macro_mm, "VolumeChange");
      res.push_back((double) total_nsubsteps, "Substeps");
    } else if(soilDomains=="single") {
      double runoff_mm = infiltration_excess_matrix_mm  + saturation_excess_matrix_mm;
      res = NumericVector::create(_["Local source/sinks"] = sum(sourceSink),
                                  _["Lateral source/sinks"] = sum(lateralFlows_mm),
                                  _["Infiltration"] = infiltration_matrix_mm,
                                  _["InfiltrationExcess"] = infiltration_excess_matrix_mm,
                                  _["SaturationExcess"] = saturation_excess_matrix_mm,
                                  _["Runoff"] = runoff_mm,
                                  _["DeepDrainage"] = drainage_matrix_mm,
                                  _["CapillarityRise"] = capillarity_matrix_mm,
                                  _["Correction"] = matrix_correction_mm,
                                  _["VolumeChange"] = Vfin_mm - Vini0_mm,
                                  _["Substeps"] = (double) total_nsubsteps);
    } 
  }
  
  return(res);
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
  List SWBcommunication = communicationSoilWaterBalance(soil.nrow());
  NumericVector res = soilWaterBalance_inner(SWBcommunication, soil, soilFunctions, 
                                             rainfallInput, rainfallIntensity, snowmelt, sourceSink, 
                                             runon, lateralFlows, waterTableDepth,
                                             infiltrationMode, infiltrationCorrection, soilDomains, 
                                             nsteps, max_nsubsteps, modifySoil);
  return(res);
}
