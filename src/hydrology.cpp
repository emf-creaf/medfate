// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include "root.h"
#include "hydraulics.h"
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
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}, \code{\link{hydrology_infiltration}}
//' 
//' 
//' @name hydrology_soilEvaporation
// [[Rcpp::export("hydrology_soilEvaporation")]]
double soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                       bool modifySoil = true) {
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector psiSoil = psi(soil, soilFunctions);
  double Esoil = 0.0;
  double swe = soil["SWE"]; //snow pack
  if(swe == 0.0) {
    double PETsoil = pet*(LgroundSWR/100.0);
    double Gsoil = soil["Gsoil"];
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
// [[Rcpp::export("hydrology_herbaceousTranspiration")]]
NumericVector herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                               List soil, String soilFunctions, bool modifySoil = true){
  if(NumericVector::is_na(herbLAI)) return(0.0);
  double Tmax_herb = pet*(LherbSWR/100.0)*(0.134*herbLAI - 0.006*pow(herbLAI, 2.0));
  NumericVector dVec = soil["W"];
  int nlayers = dVec.size();
  NumericVector psiSoil = psi(soil, soilFunctions);
  NumericVector W = soil["W"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector EherbVec(nlayers,0.0);
  NumericVector V = ldrRS_one(50, 500, dVec);
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
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}
//' 
//' @name hydrology_infiltration
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
//' @param dVec Width of soil layers (in mm).
//' @param macro Macroporosity of soil layers (in \%).
//' @param a,b Parameters of the extinction function used for water infiltration.
//' 
// [[Rcpp::export("hydrology_infiltrationRepartition")]]
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro, 
                                      double a = -0.005, double b = 3.0) {
  int nlayers = dVec.length();
  NumericVector Pvec = NumericVector(nlayers,0.0);
  NumericVector Ivec = NumericVector(nlayers,0.0);
  double z1 = 0.0;
  double p1 = 1.0;
  for(int i=0;i<nlayers;i++) {
    double ai = a*pow(1.0-macro[i],b);
    if(i<(nlayers-1)) {
      Pvec[i] = p1*(1.0-exp(ai*dVec[i]));
    } else {
      Pvec[i] = p1;
    }
    p1 = p1*exp(ai*dVec[i]);
    z1 = z1 + dVec[i];
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
//' 
// [[Rcpp::export("hydrology_infiltrationAmount")]]
double infiltrationAmount(double rainfallInput, double rainfallIntensity, List soil, 
                          String soilFunctions, String model = "GreenAmpt1911") {
  double infiltration = 0.0;
  if(model=="GreenAmpt1911") {
    NumericVector clay = soil["clay"];
    NumericVector sand = soil["sand"];
    String usda = USDAType(clay[0], sand[0]);
    NumericVector cp = campbellParamsClappHornberger(usda);
    NumericVector theta_dry = theta(soil, soilFunctions);
    double t = rainfallInput/rainfallIntensity; // time in hours
    double b = cp["b"];
    double psi_w = cp["psi_sat_cm"]*((2.0*b + 3.0)/(2*b + 6.0));
    double theta_sat = cp["theta_sat"];
    double K_sat = cp["K_sat_cm_h"];
    infiltration = infitrationGreenAmpt(t, psi_w, K_sat, theta_sat, theta_dry[0]);
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



//' Soil vertical inputs
//' 
//' High-level functions to define water inputs into the soil:
//'  
//' \itemize{
//'   \item{Function \code{hydrology_soilWaterInputs} performs canopy water interception and snow accumulation/melt.}
//'   \item{Function \code{hydrology_snowMelt} estimates snow melt using a simple energy balance, according to Kergoat (1998).}
//' }
//' 
//' @param soil A list containing the description of the soil (see \code{\link{soil}}).
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param prec Precipitation for the given day (mm)
//' @param pet Potential evapotranspiration for the given day (mm)
//' @param rainfallIntensity Rainfall intensity rate (mm/h).
//' @param tday Average day temperature (ºC).
//' @param rad Solar radiation (in MJ/m2/day).
//' @param elevation Altitude above sea level (m).
//' @param Cm Canopy water storage capacity.
//' @param LgroundPAR Percentage of photosynthetically-active radiation (PAR) reaching the ground.
//' @param LgroundSWR Percentage of short-wave radiation (SWR) reaching the ground.
//' @param runon Surface water amount running on the target area from upslope (in mm).
//' @param snowpack Boolean flag to indicate the simulation of snow accumulation and melting.
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
//' @details 
//' The function simulates different vertical hydrological processes, which are described separately in other functions. 
//' If \code{modifySoil = TRUE} the function will modify the \code{soil} object (including both soil moisture and 
//' the snowpack on its surface) as a result of simulating hydrological processes.
//' 
//' @return 
//' Function \code{hydrology_soilWaterInputs} returns a named vector with the following elements, all in mm:
//' \item{Rain}{Precipitation as rainfall.}
//' \item{Snow}{Precipitation as snow.}
//' \item{Interception}{Rainfall water intercepted by the canopy and evaporated.}
//' \item{NetRain}{Rainfall reaching the ground (throughfall).}
//' \item{Snowmelt}{Snow melted during the day, and added to the water infiltrated.}
//' \item{Runon}{Surface water amount running on the target area from upslope.}
//' \item{RainfallInput}{Rainfall input, including runon and net rain.}
//' \item{TotalInput}{Total soil input, including runon, snowmelt and net rain.}
//' 
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//'  Kergoat L. (1998). A model for hydrological equilibrium of leaf area index on a global scale. Journal of Hydrology 212–213: 268–286.
//' 
//' @seealso \code{\link{spwb_day}}, \code{\link{hydrology_rainInterception}}, \code{\link{hydrology_soilEvaporation}}
//' 
//' @param interceptionMode Infiltration model, either "Gash1995" or "Liu2001".
//' 
//' @name hydrology_verticalInputs
// [[Rcpp::export("hydrology_soilWaterInputs")]]
NumericVector soilWaterInputs(List soil, String soilFunctions, String interceptionMode,
                              double prec, double rainfallIntensity,
                              double pet, double tday, double rad, double elevation,
                              double Cm, double LgroundPAR, double LgroundSWR, 
                              double runon = 0.0,
                              bool snowpack = true, bool modifySoil = true) {
  //Soil input
  double swe = soil["SWE"]; //snow pack
  double er = pet/(24.0*rainfallIntensity);
  
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
      melt = std::min(swe, snowMelt(tday, rad, LgroundSWR, elevation));
      // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
      swe = swe-melt;
    }
  } else {
    rain = prec;
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
  if(modifySoil) {
    soil["SWE"] = swe;
  }
  NumericVector WI = NumericVector::create(_["Rain"] = rain, _["Snow"] = snow,
                                           _["Interception"] = Interception,
                                           _["NetRain"] = NetRain, 
                                           _["Snowmelt"] = melt,
                                           _["Runon"] = runon,
                                           _["RainfallInput"] = runon+NetRain,
                                           _["TotalInput"] = runon+melt+NetRain);
  return(WI);
}

// NumericVector soilInfiltrationPercolation(List soil, String soilFunctions, 
//                                           double waterInput,
//                                           bool modifySoil = true) {
//   //Soil input
//   NumericVector W = clone(Rcpp::as<Rcpp::NumericVector>(soil["W"])); //Access to soil state variable
//   NumericVector dVec = soil["dVec"];
//   NumericVector macro = soil["macro"];
//   NumericVector rfc = soil["rfc"];
//   NumericVector Water_FC = waterFC(soil, soilFunctions);
//   NumericVector Water_SAT = waterSAT(soil, soilFunctions);
//   int nlayers = W.size();
//   
// 
//   //Hydrologic input
//   double Infiltration= 0.0, Runoff= 0.0;
//   if(waterInput>0.0) {
//     //Interception
//     //Net Runoff and infiltration
//     Infiltration = infiltrationBoughton(waterInput, Water_FC[0]);
//     Runoff = waterInput - Infiltration;
//     //Decide infiltration repartition among layers
//     NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
//     //Input of the first soil layer is infiltration
//     double infiltrationExcess = 0.0;
//     double Wn;
//     for(int l=0;l<nlayers;l++) {
//       if((dVec[l]>0.0) && (Ivec[l]>0.0)) {
//         Wn = W[l]*Water_FC[l] + Ivec[l]; //Update water volume
//         if(l<(nlayers-1)) {
//           Ivec[l+1] = Ivec[l+1] + std::max(Wn - Water_FC[l],0.0); //update Ivec adding the excess to the infiltrating water (saturated flow)
//           W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta
//         } else {
//           W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta
//           infiltrationExcess = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer using field capacity
//         }
//       } 
//     } 
//     //If there still excess fill layers over field capacity
//     if((infiltrationExcess>0.0)) {
//       for(int l=(nlayers-1);l>=0;l--) {
//         if((dVec[l]>0.0) && (infiltrationExcess>0.0)) {
//           Wn = W[l]*Water_FC[l] + infiltrationExcess; //Update water volume
//           infiltrationExcess = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
//           W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta here no upper
//         }
//       }
//       //If soil is completely saturated increase surface Runoff and decrease Infiltration
//       if(infiltrationExcess>0.0) { 
//         Runoff = Runoff + infiltrationExcess;
//         Infiltration = Infiltration - infiltrationExcess;
//       }
//     } 
//   }
  // //If there is still room for additional drainage (water in macropores accumulated from previous days)
  // double head = 0.0;
  // for(int l=0;l<nlayers;l++) { //Add mm over field capacity
  //   if((l<(nlayers-1)) || rockyLayerDrainage) {
  //     head += Water_FC[l]*std::max(W[l] - 1.0, 0.0);
  //   }
  // }
  // // Rcout<<head<<"\n";
  // if(head>0.0) {
  //   double maxDrainage = head*Kdrain;
  //   for(int l=0;l<nlayers;l++) {
  //     if(maxDrainage>0.0) {
  //       double Wn = W[l]*Water_FC[l];
  //       double toDrain = std::min(std::max(Wn - Water_FC[l], 0.0), maxDrainage);
  //       if((l==(nlayers-1)) && (rfc[l] >= 95.0) && (!rockyLayerDrainage)) { //Prevent drainage for last rocky layer if not allowed
  //         toDrain = 0.0;
  //       }
  //       if(toDrain > 0.0) {
  //         DeepDrainage +=toDrain;
  //         maxDrainage -=toDrain;
  //         Wn -= toDrain;
  //         W[l] = std::max(0.0, Wn/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
  //       }
  //     }
  //   }
  // }
//   if(modifySoil) {
//     NumericVector Ws = soil["W"];
//     for(int l=0;l<nlayers;l++) Ws[l] = W[l];
//   }
//   NumericVector DB = NumericVector::create(_["Infiltration"] = Infiltration, 
//                                            _["Runoff"] = Runoff);
//   return(DB);
// }


NumericVector tridiagonalSolving(NumericVector a, NumericVector b, NumericVector c, NumericVector d) {
  int n = a.size();
  NumericVector e(n), f(n), u(n);
  
  //Forward steps
  double e_prev = 0.0;
  double f_prev = 0.0;
  for(int i=0;i<n;i++) {
    e[i] = c[i]/(b[i] - a[i]*e_prev);
    f[i] = (d[i] - a[i]*f_prev)/(b[i] - a[i]*e_prev);
    // Rcout<<i<< " "<< e[i]<< " "<< f[i]<<"\n";
    e_prev = e[i];
    f_prev = f[i];
  }
  //Backward steps
  u[n-1] = f[n-1];
  for(int i = (n - 2);i>=0;i--) {
    u[i] = f[i] - e[i]*u[i + 1];
  }  
  return(u);
}


//' Soil flows
//' 
//' Function \code{hydrology_soilFlows} estimates water movement within the soil according to Richards equation.
//' 
//' @param soil Object of class \code{\link{soil}}.
//' @param sourceSink Source/sink term for each soil layer (from snowmelt, soil evaporation or plant transpiration/redistribution)
//'        as mm/day.
//' @param nsteps  Number of time steps per day
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}, \code{\link{hydrology_infiltration}}
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @return
//'   Returns the water draining from the bottom layer.
//'   
//' @name hydrology_soilFlows
// [[Rcpp::export("hydrology_soilFlows")]]
double soilFlows(List soil, NumericVector sourceSink, int nsteps = 24,
                 bool modifySoil = true) {
  
  NumericVector dVec = soil["dVec"];
  double mm_day_2_m3_s = 0.001*(1.0/86400.0);//From mm/day = l/day = dm3/day to m3/m2/s
  NumericVector sourceSink_m3s =  sourceSink*mm_day_2_m3_s;
  NumericVector dZ_m = dVec*0.001; //mm to m
  NumericVector rfc = soil["rfc"];

  double maxSource =  max(abs(sourceSink));
  if(maxSource > 10) nsteps = 48;
  if(maxSource > 20) nsteps = 96;
  if(maxSource > 30) nsteps = 144;
  if(maxSource > 40) nsteps = 192;
  
  double tstep = 86400.0/((double) nsteps);
  double halftstep = tstep/2.0;
  
  //Estimate layer interfaces
  int nlayers = dVec.size();
  NumericVector dZUp(nlayers), dZDown(nlayers), lambda(nlayers);
  for(int l=0;l<nlayers;l++) {
    lambda[l] = 1.0 - rfc[l]/100.0;
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
  
  //Retrieve VG parameters
  NumericVector Ksat =soil["Ksat"];
  NumericVector n =soil["VG_n"];
  NumericVector alpha = soil["VG_alpha"];
  NumericVector theta_res = soil["VG_theta_res"];
  NumericVector theta_sat = soil["VG_theta_sat"];
  
  //Estimate Theta, Psi, C, K
  NumericVector Theta = theta(soil, "VG");
  NumericVector Psi(nlayers), K(nlayers), C(nlayers);
  NumericVector Psi_m(nlayers), K_ms(nlayers), C_m(nlayers), K_ms05(nlayers), C_m05(nlayers);
  for(int l=0;l<nlayers;l++) {
    Psi[l] = theta2psiVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], Theta[l]); 
    C[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
    K[l] = psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
    Psi_m[l]= Psi[l]/mTOMPa; // MPa to m
    K_ms[l] = 0.01*K[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
    C_m[l] = C[l]*mTOMPa; //From MPa-1 to m-1
  }
  

  double drainage = 0.0;
  double K_up, K_down;
  NumericVector a(nlayers), b(nlayers), c(nlayers), d(nlayers);
  //Psi-based solution of the Richards equation using implicit solution for psi
  //but with explicit linearization for K and C (pp. 126, Bonan)
  for(int s =0;s<nsteps;s++) {
    //A. Predictor step
    for(int l=0;l<nlayers;l++) {
      if(l==0) { //first layer
        K_up = K_ms[l];
        K_down = 0.5*(K_ms[l] + K_ms[l+1]);
        a[l] = 0.0;
        c[l] = -1.0*K_down/dZDown[l];
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep) - c[l];
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep)*Psi_m[l]  - K_down + sourceSink_m3s[l];
      } else if(l<(nlayers - 1)) {
        K_up = 0.5*(K_ms[l-1] + K_ms[l]);
        K_down = 0.5*(K_ms[l] + K_ms[l+1]);
        a[l] = -1.0*K_up/dZUp[l];
        c[l] = -1.0*K_down/dZDown[l];
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep) - a[l] - c[l];
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep)*Psi_m[l] + K_up - K_down + sourceSink_m3s[l];
      } else { // last layer
        K_up = 0.5*(K_ms[l-1] + K_ms[l]);
        K_down = K_ms[l];
        a[l] = -1.0*K_up/dZUp[l];
        c[l] = 0.0;
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep) - a[l] - c[l];
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/halftstep)*Psi_m[l] + K_up - K_down + sourceSink_m3s[l];
      }
    }
    NumericVector Psi_m_t05 = tridiagonalSolving(a,b,c,d);
    NumericVector Psi_t05 = Psi_m_t05*mTOMPa; // m to MPa
    //Calculate K and C at t05
    for(int l=0;l<nlayers;l++) {
      double C_05 = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_t05[l]);
      double K_05 = psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi_t05[l]);
      K_ms05[l] = 0.01*K_05/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      C_m05[l] = C_05*mTOMPa; //From MPa-1 to m-1
    }
    //B. Corrector step also Crank-Nicolson
    for(int l=0;l<nlayers;l++) {
      if(l==0) { //first layer
        K_up = K_ms05[l];
        K_down = 0.5*(K_ms05[l] + K_ms05[l+1]);
        a[l] = 0.0;
        c[l] = -1.0*K_down/(2.0*dZDown[l]);
        b[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep) - c[l];
        d[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep)*Psi_m[l] - c[l]*(Psi_m[l] - Psi_m[l+1])  - K_down + sourceSink_m3s[l];
      } else if(l<(nlayers - 1)) {
        K_up = 0.5*(K_ms05[l-1] + K_ms05[l]);
        K_down = 0.5*(K_ms05[l] + K_ms05[l+1]);
        a[l] = -1.0*K_up/(2.0*dZUp[l]);
        c[l] = -1.0*K_down/(2.0*dZDown[l]);
        b[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep) - a[l] - c[l];
        d[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep)*Psi_m[l] + a[l]*(Psi_m[l - 1] - Psi_m[l]) - c[l]*(Psi_m[l] - Psi_m[l+1]) + K_up - K_down + sourceSink_m3s[l];
      } else { // last layer
        K_up = 0.5*(K_ms05[l-1] + K_ms05[l]);
        K_down = K_ms05[l];
        a[l] = -1.0*K_up/(2.0*dZUp[l]);
        c[l] = 0.0;
        b[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep) - a[l];
        d[l] = (lambda[l]*C_m05[l]*dZ_m[l]/tstep)*Psi_m[l] + a[l]*(Psi_m[l - 1] - Psi_m[l])  + K_up - K_down + sourceSink_m3s[l];
      }
    }
    NumericVector Psi_m_t1 = tridiagonalSolving(a,b,c,d);
    NumericVector Psi_t1 = Psi_m_t1*mTOMPa; // m to MPa
    //calculate free drainage (m3)
    drainage += std::max(0.0, K_ms05[nlayers -1]*tstep);
    //Update psi, capacitances and conductances for next step
    for(int l=0;l<nlayers;l++) {
      // Rcout<<" step "<<s<<" layer " <<l<< " "<< Psi[l]<< " to " << Psi_t1[l]<<"\n";
      Psi[l] = Psi_t1[l];
      C[l] = psi2cVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
      K[l] = psi2kVanGenuchten(Ksat[l], n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
      Psi_m[l]= Psi[l]/mTOMPa; // MPa to m
      K_ms[l] = 0.01*K[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
      C_m[l] = C[l]*mTOMPa; //From MPa-1 to m-1
    }
  }
  double drainage_mm = drainage*1000.0; //m3/m2 to mm/m2
  if(modifySoil) {
    NumericVector W = soil["W"];
    for(int l=0;l<nlayers;l++) {
      double theta_fc = psi2thetaVanGenuchten(n[l],alpha[l],theta_res[l], theta_sat[l], -0.033);
      Theta[l] = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi[l]);
      W[l] = Theta[l]/theta_fc;
    }
  }
  return(drainage_mm);
}