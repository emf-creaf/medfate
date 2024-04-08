// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include "root.h"
#include "hydraulics.h"
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
    NumericVector bd = soil["bd"];
    double bdsoil = 2.73; //Density of soil particles
    double bdref = 1.2; //Reference bulk density for Ksat
    String usda = USDAType(clay[0], sand[0]);
    NumericVector cp = campbellParamsClappHornberger(usda);
    NumericVector theta_dry = theta(soil, soilFunctions);
    double t = std::min(24.0, rainfallInput/rainfallIntensity); // time in hours
    double b = cp["b"];
    double psi_w = cp["psi_sat_cm"]*((2.0*b + 3.0)/(2*b + 6.0));
    double theta_sat = cp["theta_sat"];
    double K_sat = cp["K_sat_cm_h"];
    K_sat = K_sat*std::pow((bdsoil - bd[0])/(bdsoil - bdref),3.0);
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
  double S_w = std::max(0.0, ((G_f*D_w*gamma_w)/pow(d, 2.0))*(theta_b - theta_micro));
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

//' Soil flows
//' 
//' Function \code{hydrology_soilFlows} estimates water balance of soil layers given water inputs/outputs, including the simulation of water movement within the soil.
//' 
//' @param soil Object of class \code{\link{soil}}.
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param rainfallInput Amount of water from rainfall event (after excluding interception), in mm.
//' @param rainfallIntensity Rainfall intensity, in mm/h.
//' @param snowmelt Amount of water originated from snow melt, in mm.
//' @param sourceSink Source/sink term for each soil layer (from soil evaporation or plant transpiration/redistribution)
//'        as mm/day.
//' @param infiltrationMode Infiltration model, either "GreenAmpt1911" or "Boughton1989"
//' @param soilDomains Either "single" (for single-domain) or "dual" (for dual-permeability).
//' @param freeDrainage Boolean flag to indicate that lower boundary condition is free drainage.
//' @param nsteps  Number of time steps per day
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
//' @seealso  \code{\link{spwb}}, \code{\link{hydrology_soilWaterInputs}}, \code{\link{hydrology_infiltration}}
//' 
//' @author 
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' \enc{María González Sanchís}{Maria Gonzalez Sanchis}, UPV-CTFC
//' 
//' @return
//'   Returns a named vector with different elements, depending on \code{soilDomains}. If
//'   \code{soilDomains == "single"}:
//'   \itemize{
//'     \item{\code{Snowmelt}: Snowmelt input (mm).}
//'     \item{\code{Source/sinks}: Sum of source/sink input across layers (mm).}
//'     \item{\code{Infiltration}: Water infiltrated into the soil (mm).}
//'     \item{\code{Runoff}: Surface runoff generated (mm).}
//'     \item{\code{DeepDrainage}: Water draining from the bottom layer (mm). This quantity is corrected to close the water balance.}
//'     \item{\code{VolumeChange}: Change in soil water volume (mm).}
//'     \item{\code{UncorrectedWaterBalance}: Uncorrected balance of inputs/outputs (mm).}
//'   }
//'  If \code{soilDomains == "dual"} the named vector contains the following additional elements:
//'   \itemize{
//'     \item{\code{Lateral flows}: Sum of water circulating between micropores and macropores, positive when filling micropores (mm).}
//'     \item{\code{InfiltrationMicropores}: Water infiltrated into the soil matrix (mm).}
//'     \item{\code{InfiltrationMacropores}: Water infiltrated into the soil macropore domain (mm).}
//'     \item{\code{DrainageMicropores}: Water draining from the bottom layer of the micropore domain (mm). This quantity is corrected to close water balance in the micropore domain.}
//'     \item{\code{DrainageMacropores}: Water draining from the bottom layer of the macropore domain (mm). This quantity is corrected to close the water balance in the macropore domain.}
//'     \item{\code{MicroporeVolumeChange}: Change in soil water volume in the micropore domain (mm).}
//'     \item{\code{UncorrectedMicroporeBalance}: Uncorrected balance of inputs/outputs in the micropore domain (mm).}
//'     \item{\code{MacroporeVolumeChange}: Change in soil water volume in the macropore domain (mm).}
//'     \item{\code{UncorrectedMacroporeBalance}: Uncorrected balance of inputs/outputs in the macropore domain (mm).}
//'   }
//'   
//' @examples
//' # Initialize soil example
//' examplesoil <- soil(defaultSoilParams(4))
//' 
//' # Water balance in a single-domain simulation (Richards equation)
//' hydrology_soilWaterBalance(examplesoil, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
//'                            soilDomains = "single", modifySoil = FALSE)
//'                     
//' # Water balance in a dual-permeability model (MACRO)
//' hydrology_soilWaterBalance(examplesoil, "VG", 10, 5, 0, c(-1,-1,-1,-1), 
//'                            soilDomains = "dual", modifySoil = FALSE)
//'   
//' @name hydrology_soilWaterBalance
// [[Rcpp::export("hydrology_soilWaterBalance")]]
NumericVector soilWaterBalance(List soil, String soilFunctions, 
                               double rainfallInput, double rainfallIntensity, double snowmelt, NumericVector sourceSink, 
                               String infiltrationMode = "GreenAmpt1911", String soilDomains = "single", bool freeDrainage = true, 
                               int nsteps = 24, bool modifySoil = true) {
  
  if((soilDomains!="single") && (soilDomains!="dual")) stop("Unrecognized soilDomain value");
  if(!modifySoil) soil = clone(soil);
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector W = soil["W"];
  int nlayers = dVec.size();
  
  bool micropore_excess = true;
  bool micropore_imbibition = true;
  
  double mm_2_m3 = 0.001;
  double m3_2_mm = 1.0/mm_2_m3;
  double mm_day_2_m3_s = mm_2_m3*(1.0/86400.0);//From mm/day = l/m2/day = dm3/day to m3/m2/s
  
  //Infiltration
  double infiltration_matrix_mm = infiltrationAmount(rainfallInput, rainfallIntensity, soil, 
                                                     soilFunctions, infiltrationMode);
  double infiltration_macropores_mm = 0.0;
  double infiltration_matrix_excess_mm = rainfallInput - infiltration_matrix_mm;
  double runoff_mm = 0.0;
  double snowmelt_mm = snowmelt;
  double correction_mm = 0.0;
  double micropore_correction_mm = 0.0;
  double macropore_correction_mm = 0.0;
  
  //Initialize lateral flows (positive in the direction of matrix)
  NumericVector lateral_flows_mm(nlayers, 0.0);
  
  //Copy sinks
  NumericVector source_sink_def_mm(nlayers);
  for(int l=0;l<nlayers;l++) {
    source_sink_def_mm[l] = sourceSink[l];
  }
  //Add snow-melt to source/sinks
  source_sink_def_mm[0] += snowmelt_mm;
  if(soilDomains=="single") {
    runoff_mm = infiltration_matrix_excess_mm;
    NumericVector IVec = infiltrationRepartition(infiltration_matrix_mm, dVec, macro);
    //Add infiltration to matrix def source/sinks
    for(int l=0;l<nlayers;l++) {
      source_sink_def_mm[l] += IVec[l];
    }
  } else {
    //Add matrix infiltration to first layer
    source_sink_def_mm[0] += infiltration_matrix_mm;
  }

  double tstep = 86400.0/((double) nsteps);
  int nsubsteps = 1;
  int max_nsubsteps = 3600;
  double tsubstep = tstep;
  
  NumericVector source_sink_def_m3s =  source_sink_def_mm*mm_day_2_m3_s;
  NumericVector matrixImbibition_m3s(nlayers, 0.0);
  NumericVector matrixExcess_m3s(nlayers, 0.0);
  NumericVector dZ_m = dVec*0.001; //mm to m

  //Estimate layer interfaces
  NumericVector dZUp(nlayers), dZDown(nlayers), lambda(nlayers);
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
  
  //boundary condition of water table (if freeDrainage = FALSE)
  double Psi_bc = 0.0;
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
  NumericVector theta_micro(nlayers), theta_b(nlayers), theta_macro(nlayers), theta_sat_fict(nlayers), Ksat_b(nlayers), Ksat_b_ms(nlayers);
  NumericVector Ksat(nlayers), Ksat_ms(nlayers);
  NumericVector Psi(nlayers), K(nlayers), C(nlayers);
  NumericVector Psi_m(nlayers), K_ms(nlayers), Kbc(nlayers), Kbc_ms(nlayers), C_m(nlayers);
  //Macroporosoty domain
  NumericVector S_macro(nlayers), e_macro(nlayers), Kmacro_ms(nlayers);
  NumericVector waterFluidity(nlayers, 1.0);
  for(int l=0;l<nlayers;l++) {
    if(!NumericVector::is_na(Tsoil[l])) {
      if(Tsoil[l]>0) {
        waterFluidity[l] = 1.0/waterDynamicViscosity(Tsoil[l]); 
      } else {
        waterFluidity[l] = 0.0;
      }
    }
    Ksat[l] = Ksat_ori[l]*lambda[l];//Multiply K for the space available for water movement
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
      Ksat_ms[l] = 0.01*waterFluidity[l]*Ksat[l]/(86400.0*cmdTOmmolm2sMPa); //mmolH20*m-2*MPa-1*s-1 to m*s-1
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
      Vini_micro_mm += theta_micro[l]*dVec[l]*lambda[l];
      Vini_macro_mm += theta_macro[l]*dVec[l]*lambda[l];
    }
    Vini_step_micro_mm = Vini_micro_mm;
    Vini_step_macro_mm = Vini_macro_mm;
  } else {
    for(int l=0;l<nlayers;l++) {
      Vini_mm += Theta[l]*dVec[l]*lambda[l];
    }
    Vini_step_mm = Vini_mm;
  }
  double Vini0_mm = Vini_mm;
  double Vini0_macro_mm = Vini_macro_mm;
  double Vini0_micro_mm = Vini_micro_mm;
  
  double drainage_matrix_m3 = 0.0, drainage_macropores_m3 = 0.0;
  double saturation_excess_mm = 0.0;
  double saturation_excess_micropores_mm = 0.0;
  double saturation_excess_macropores_mm = 0.0;
  double K_up, K_down;
  NumericVector a(nlayers), b(nlayers), c(nlayers), d(nlayers);
  //Temporary step variables
  NumericVector C_step(nlayers), C_step_m(nlayers), K_step_ms(nlayers), K_step(nlayers), Psi_step(nlayers), Psi_step_m(nlayers);
  NumericVector S_macro_step(nlayers), Kmacro_step_ms(nlayers), theta_macro_step(nlayers), theta_micro_step(nlayers);
  for(int s =0;s<nsteps;s++) {
    double drainage_matrix_step_m3;
    double drainage_macropores_step_m3;
    double saturation_excess_step_mm;
    double saturation_excess_micropores_step_mm;
    double saturation_excess_macropores_step_mm;
    double infiltration_macropores_step_mm;
    double correction_step_mm;
    double macropore_correction_step_mm;
    double micropore_correction_step_mm;
    double infiltration_matrix_excess_step_mm; 
    NumericVector lateral_flows_step_mm(nlayers);
    
    bool cont = true;
    while(cont) {
      infiltration_matrix_excess_step_mm = infiltration_matrix_excess_mm;
      saturation_excess_macropores_step_mm = 0.0;
      saturation_excess_micropores_step_mm = 0.0;
      saturation_excess_step_mm = 0.0;
      drainage_matrix_step_m3 = 0.0;
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
        S_macro_step[l] = S_macro[l];
        Kmacro_step_ms[l] = Kmacro_ms[l];
        theta_macro_step[l] = theta_macro[l];
        theta_micro_step[l] = theta_micro[l];
      }
      
      int nsubsteps_day = (nsteps*nsubsteps);
      double rainfallIntensity_substep = rainfallIntensity*24.0/((double) nsubsteps_day); //mm/substep
      tsubstep = tstep/((double) nsubsteps);
      double halftsubstep = tsubstep/2.0;
      
      for(int ss=0;ss<nsubsteps;ss++) {
        //Psi-based solution of the Richards equation using implicit solution for psi
        //but with explicit linearization for K and C (pp. 126, Bonan)
        //A. Predictor sub-step
        for(int l=0;l<nlayers;l++) {
          if(l==0) { //first layer
            K_up = K_step_ms[l];
            K_down = 0.5*(K_step_ms[l] + K_step_ms[l+1]);
            a[l] = 0.0;
            c[l] = -1.0*K_down/dZDown[l];
            b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - c[l];
            d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l]  - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l];
          } else if(l<(nlayers - 1)) {
            K_up = 0.5*(K_step_ms[l-1] + K_step_ms[l]);
            K_down = 0.5*(K_step_ms[l] + K_step_ms[l+1]);
            a[l] = -1.0*K_up/dZUp[l];
            c[l] = -1.0*K_down/dZDown[l];
            b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - a[l] - c[l];
            d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l] + K_up - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l];
          } else { // last layer
            K_up = 0.5*(K_step_ms[l-1] + K_step_ms[l]);
            a[l] = -1.0*K_up/dZUp[l];
            double drain_bc = 0.0;
            if(freeDrainage) {
              K_down = K_step_ms[l];
            } else {
              K_down = K_step_ms[l]; // 0.5*(K_step_ms[l] + Kbc_ms[l]);
              drain_bc = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_bc);
            }
            c[l] = 0.0;
            b[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep) - a[l];
            d[l] = (lambda[l]*C_step_m[l]*dZ_m[l]/halftsubstep)*Psi_step_m[l] + K_up - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l] + drain_bc;
          }
        }
        NumericVector Psi_step_m_t05 = tridiagonalSolving(a,b,c,d);
        NumericVector Psi_step_t05 = Psi_step_m_t05*mTOMPa; // m to MPa
        //Calculate K and C at t05
        double C_step_05, K_step_05;
        NumericVector K_step_ms05(nlayers), C_step_m05(nlayers);
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
            K_up = K_step_ms05[l];
            K_down = 0.5*(K_step_ms05[l] + K_step_ms05[l+1]);
            a[l] = 0.0;
            c[l] = -1.0*K_down/(2.0*dZDown[l]);
            b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - c[l];
            d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] - c[l]*(Psi_step_m[l] - Psi_step_m[l+1])  - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l];
          } else if(l<(nlayers - 1)) {
            K_up = 0.5*(K_step_ms05[l-1] + K_step_ms05[l]);
            K_down = 0.5*(K_step_ms05[l] + K_step_ms05[l+1]);
            a[l] = -1.0*K_up/(2.0*dZUp[l]);
            c[l] = -1.0*K_down/(2.0*dZDown[l]);
            b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - a[l] - c[l];
            d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] + a[l]*(Psi_step_m[l - 1] - Psi_step_m[l]) - c[l]*(Psi_step_m[l] - Psi_step_m[l+1]) + K_up - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l];
          } else { // last layer
            K_up = 0.5*(K_step_ms05[l-1] + K_step_ms05[l]);
            a[l] = -1.0*K_up/(2.0*dZUp[l]);
            double drain_bc = 0.0;
            if(freeDrainage) {
              K_down = K_step_ms05[l];
            } else {
              K_down = K_step_ms05[l]; //0.5*(K_step_ms05[l] + Kbc_ms[l]);
              drain_bc = -1.0*K_down/dZDown[l]*(Psi_step_m[l] - Psi_bc);
            }
            c[l] = 0.0;
            b[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep) - a[l] - c[l];
            d[l] = (lambda[l]*C_step_m05[l]*dZ_m[l]/tsubstep)*Psi_step_m[l] + a[l]*(Psi_step_m[l - 1] - Psi_step_m[l])  + K_up - K_down + source_sink_def_m3s[l] + matrixImbibition_m3s[l] + drain_bc;
          }
        }
        NumericVector Psi_step_m_t1 = tridiagonalSolving(a,b,c,d);
        NumericVector Psi_step_t1 = Psi_step_m_t1*mTOMPa; // m to MPa
        //calculate drainage (m3)
        if(freeDrainage) {
          drainage_matrix_step_m3 += K_step_ms05[nlayers -1]*tsubstep;
        } else {
          //Calculate outward/inward flux between last layer and boundary condition
          K_down = K_step_ms05[nlayers-1];
          double flow = K_down/dZDown[nlayers-1]*(Psi_step_m[nlayers-1] - Psi_bc) + K_down;
          // Rcout<< dZDown[nlayers-1]<<" "<< (Psi_step_m[nlayers-1] - Psi_bc)<<" "<<K_down<<" " << flow <<"\n";
          drainage_matrix_step_m3 += flow*tsubstep;
        }
        
        //Update Psi and theta
        for(int l=0;l<nlayers;l++) {
          // Rcout<<" step "<<s<<" layer " <<l<< " "<< Psi_step[l]<< " to " << Psi_step_t1[l]<<"\n";
          Psi_step[l] = Psi_step_t1[l];
          
          if((soilDomains=="single")) {
            if((Psi_step[l] > -0.0001)) { //oversaturation generates saturation excess
              double new_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], Psi_step[l]);
              double sat_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat[l], -0.0001);
              double res_mm = std::abs(sat_theta - new_theta)*dVec[l]*lambda[l];
              saturation_excess_step_mm += res_mm;
              Psi_step[l] = -0.0001;
            }
          } else {
            if((Psi_step[l] > -0.0001)) { //oversaturation generates micropore saturation excess
              double new_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
              double sat_theta = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], -0.0001);
              double res_mm = std::abs(sat_theta - new_theta)*dVec[l]*lambda[l];
              saturation_excess_micropores_step_mm += res_mm;
              Psi_step[l] = -0.0001;
            }
            //Update theta_micro for the next substep
            theta_micro_step[l] = psi2thetaVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], Psi_step[l]);
            //If needed, add exceeding moisture to the macropores and correct Psi
            double excess_theta_step = 0.0;
            if(theta_micro_step[l] > theta_b[l]) {
              //Maximum macropore capacity 
              double C_macro_step = (1.0 - S_macro_step[l])*(theta_sat[l] - theta_b[l]);
              double excess_theta_step = std::min(theta_micro_step[l] - theta_b[l], C_macro_step);
              theta_micro_step[l] -= excess_theta_step;
              Psi_step[l] = theta2psiVanGenuchten(n[l], alpha[l], theta_res[l], theta_sat_fict[l], theta_micro_step[l]);
            } 
            if(excess_theta_step>0.0){
              lateral_flows_step_mm[l] -= excess_theta_step*dVec[l]*lambda[l]; //negative flow (in mm/step)
              matrixExcess_m3s[l] = excess_theta_step*dVec[l]*lambda[l]*0.001/tsubstep; //Source of water flowing into macropores (m3/s)
              // Rcout<< excess_theta <<" to macropores in "<<l<<"\n";
            } else {
              matrixExcess_m3s[l] = 0.0;
            }
          }
        }
        
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
            if(l==0) { //first layer
              K_up = 0.0;
              //Infiltration into macropores
              if(infiltration_matrix_excess_step_mm > 0.0) {
                //Maximum infiltration in this time step
                infiltration_macropores_substep_mm = std::min(infiltration_matrix_excess_step_mm, rainfallIntensity_substep);
                // Rcout << "infiltration step " << infiltration_macropores_substep_mm<<"\n";
                double infiltration_macropores_substep_m3s = infiltration_macropores_substep_mm*((double) nsubsteps_day)*mm_day_2_m3_s; //From mm/substep to m3s
                sourceSink_macro_m3s += infiltration_macropores_substep_m3s;
                infiltration_matrix_excess_step_mm -= infiltration_macropores_substep_mm;
              }
            } else {
              K_up = Kmacro_step_ms[l-1]; //Get last updated value
            }
            //Find solution for half substep
            double S_t1 = rootFindingMacropores(S_macro_step[l], K_up, Ksat_ms[l], Ksat_b_ms[l], kin_exp,
                                                e_macro[l], lambda[l], dZ_m[l], sourceSink_macro_m3s, tsubstep);
            // Rcout << "S " << S_macro_step[l] << " source/sink step " << sourceSink_macro_m3s<< " S1 "<< S_t1<<"\n";
            // 
            //Update macropore conductances for next step (sets K_up for next layer)
            Kmacro_step_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_t1, kin_exp);
            //Update S_macro_step using full substep
            S_macro_step[l] = S_t1;
            // S_macro_step[l] = S_macro_step[l] + (tsubstep/(e_macro[l]*lambda[l]*dZ_m[l]))*((K_up - Kmacro_step_ms[l]) + sourceSink_macro_m3s);
            
            if(micropore_imbibition) {
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
              lateral_flows_step_mm[l] += imbibitionRate*dVec[l]*lambda[l]*tsubstep; //From m3/m3/s to mm/step
            }
          }
          
          //Update theta_macro for the next step
          double res_mm = 0.0;
          for(int l=(nlayers-1);l>=0;l--) {
            theta_macro_step[l] = S_macro_step[l]*e_macro[l] + res_mm/(dVec[l]*lambda[l]);
            // If oversaturation of macroporosity occurs
            if(theta_macro_step[l] > e_macro[l]) {
              res_mm = dVec[l]*(theta_macro_step[l] - e_macro[l])*lambda[l]; //residue in mm for layers above
              // Rcout << " layer "<< l <<" residue "<< res_mm<<"\n";
              theta_macro_step[l] = e_macro[l];
              S_macro_step[l] = 1.0; //Correct S_macro to saturation
              Kmacro_step_ms[l] = (Ksat_ms[l] - Ksat_b_ms[l])*pow(S_macro_step[l], kin_exp);
              if(l==0) {
                infiltration_macropores_substep_mm = infiltration_macropores_substep_mm - res_mm; 
                // Rcout << " layer "<< l <<" infiltration "<< infiltration_macropores_substep_mm<<"\n";
              }
            } else {
              res_mm = 0.0;
            }
          }
          //Generate saturation excess if there was some residue in the top layer
          saturation_excess_macropores_step_mm += res_mm;
          //Drainage
          drainage_macropores_step_m3 += Kmacro_step_ms[nlayers-1]*tsubstep;
          infiltration_macropores_step_mm += infiltration_macropores_substep_mm;
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
          Vfin_micro_mm += theta_micro_step[l]*dVec[l]*lambda[l];
          Vfin_macro_mm += theta_macro_step[l]*dVec[l]*lambda[l];
        }
        Vfin_mm += Theta[l]*dVec[l]*lambda[l];
        W[l] = Theta[l]/Theta_FC[l];
      }
      
      double max_abs_correction_step_mm; 
      if(soilDomains=="single") {
        double balance_step_mm = m3_2_mm*(tstep*sum(source_sink_def_m3s) - drainage_matrix_step_m3) - saturation_excess_step_mm;
        correction_step_mm = balance_step_mm + Vini_step_mm - Vfin_mm;
        max_abs_correction_step_mm = correction_step_mm;
        // Rcout << s << " "<< nsubsteps<<" ini "<< Vini_step_mm <<" fin " <<  Vfin_mm << " dif "<< Vfin_mm - Vini_step_mm << "  bal "<< balance_step_mm<< " corr "<< correction_step_mm<< " sat_excess "<< saturation_excess_step_mm<<"\n";
      } else {
        double balance_micro_step_mm = m3_2_mm*(tstep*sum(source_sink_def_m3s) - drainage_matrix_step_m3) + sum(lateral_flows_step_mm) - saturation_excess_micropores_step_mm;
        double balance_macro_step_mm = infiltration_macropores_step_mm - m3_2_mm*drainage_macropores_step_m3 - sum(lateral_flows_step_mm);
        // Correct possible mismatch between balance and volume change
        micropore_correction_step_mm = balance_micro_step_mm + Vini_step_micro_mm - Vfin_micro_mm;
        macropore_correction_step_mm = balance_macro_step_mm + Vini_step_macro_mm - Vfin_macro_mm;
        // Rcout << s << " "<< nsubsteps<<" micro ini "<< Vini_step_micro_mm <<" fin " <<  Vfin_micro_mm << " dif "<< Vfin_micro_mm - Vini_step_micro_mm << "  bal "<< balance_micro_step_mm<< " corr "<< micropore_correction_step_mm<< "\n";
        // Rcout << s << " "<< nsubsteps<<" macro ini "<< Vini_step_macro_mm <<" fin " <<  Vfin_macro_mm << " dif "<< Vfin_macro_mm - Vini_step_macro_mm <<
        //   " inf " << infiltration_macropores_step_mm<< " sat exc " << saturation_excess_macropores_step_mm << " dra " << m3_2_mm*drainage_macropores_step_m3 <<
        //   " lat " << sum(lateral_flows_step_mm) << " bal "<< balance_macro_step_mm<<  " corr "<< macropore_correction_step_mm<< "\n";
        max_abs_correction_step_mm = std::max(std::abs(micropore_correction_step_mm), std::abs(macropore_correction_step_mm));
      }
      if(max_abs_correction_step_mm > 0.1) {
        nsubsteps = nsubsteps*2;
        if(nsubsteps >= max_nsubsteps) {
          nsubsteps = max_nsubsteps;
          cont = false;
        }
      } else {
        cont = false;
        if(max_abs_correction_step_mm < 0.0001) {
          if(nsubsteps > 1) {
            nsubsteps = nsubsteps/2; 
          }
        }
      }
    }
    //Copy final values
    //Add correction
    if(soilDomains=="single") {
      drainage_matrix_step_m3 += correction_step_mm*mm_2_m3;
    } else {
      drainage_matrix_step_m3 += micropore_correction_step_mm*mm_2_m3;
      drainage_macropores_step_m3 += macropore_correction_step_mm*mm_2_m3;
    }
    
    //Add to totals
    drainage_matrix_m3 += drainage_matrix_step_m3;
    if(soilDomains=="single") {
      correction_mm += correction_step_mm;
      saturation_excess_mm += saturation_excess_step_mm;
    } else {
      //Copy back infiltration excess
      infiltration_matrix_excess_mm = infiltration_matrix_excess_step_mm;
      //add to totals
      drainage_macropores_m3 += drainage_macropores_step_m3;
      infiltration_macropores_mm +=infiltration_macropores_step_mm;
      micropore_correction_mm += micropore_correction_step_mm;
      macropore_correction_mm += macropore_correction_step_mm;
      saturation_excess_micropores_mm += saturation_excess_micropores_step_mm;
      saturation_excess_macropores_mm += saturation_excess_macropores_step_mm;
    }
    
    for(int l=0;l<nlayers;l++) {
      lateral_flows_mm[l] += lateral_flows_step_mm[l];
      //Copy back step variables
      C[l] = C_step[l];
      C_m[l] = C_step_m[l];
      K[l] = K_step[l];
      K_ms[l] = K_step_ms[l];
      Psi[l] = Psi_step[l];
      Psi_m[l] = Psi_step_m[l];
      S_macro[l] = S_macro_step[l];
      Kmacro_ms[l] = Kmacro_step_ms[l];
      theta_macro[l] = theta_macro_step[l];
      theta_micro[l] = theta_micro_step[l];
      // Rcout << s << " "<< l << " th_mi " << theta_micro[l] << " K_mi " << K_ms[l] << " th_ma " << theta_macro[l] << " K_ma " << Kmacro_ms[l]<< " mat_im " << matrixImbibition_m3s[l] << " mat_ex " << matrixExcess_m3s[l] <<"\n";
    }
    //Prepare for next step
    Vini_step_macro_mm = Vfin_macro_mm;
    Vini_step_micro_mm = Vfin_micro_mm;
    Vini_step_mm = Vfin_mm;
  }
  
  double drainage_matrix_mm = drainage_matrix_m3*1000.0; //m3/m2 to mm/m2
  
  //Output
  NumericVector res;
  if(soilDomains=="dual") {
    double drainage_macropores_mm = drainage_macropores_m3*1000.0; //m3/m2 to mm/m2
    //Correct infiltration and runoff according to saturation excess
    infiltration_matrix_mm = infiltration_matrix_mm - saturation_excess_micropores_mm;
    runoff_mm = runoff_mm  + saturation_excess_micropores_mm + saturation_excess_macropores_mm;
    res = NumericVector::create(_["Snowmelt"] = snowmelt_mm,
                                _["Source/sinks"] = sum(sourceSink),
                                _["Lateral flows"] = sum(lateral_flows_mm),
                                _["InfiltrationMicropores"] = infiltration_matrix_mm,
                                _["InfiltrationMacropores"] = infiltration_macropores_mm,
                                _["SaturationExcessMicropores"] = saturation_excess_micropores_mm,
                                _["SaturationExcessMacropores"] = saturation_excess_macropores_mm,
                                _["DrainageMicropores"] = drainage_matrix_mm,
                                _["DrainageMacropores"] = drainage_macropores_mm,
                                _["DrainageMicroporeCorrection"] = micropore_correction_mm,
                                _["DrainageMacroporeCorrection"] = macropore_correction_mm,
                                _["MicroporeVolumeChange"] = Vfin_micro_mm - Vini0_micro_mm,
                                _["MacroporeVolumeChange"] = Vfin_macro_mm - Vini0_macro_mm,
                                _["SaturationExcess"] = saturation_excess_macropores_mm + saturation_excess_micropores_mm,
                                _["Infiltration"] = infiltration_matrix_mm + infiltration_macropores_mm,
                                _["Runoff"] = runoff_mm,
                                _["DeepDrainage"] = drainage_macropores_mm + drainage_matrix_mm,
                                _["VolumeChange"] = Vfin_micro_mm - Vini0_micro_mm  + Vfin_macro_mm - Vini0_macro_mm,
                                _["Substep"] = tsubstep);
  } else {
    //Correct infiltration and runoff according to saturation excess
    infiltration_matrix_mm = infiltration_matrix_mm - saturation_excess_mm;
    runoff_mm = runoff_mm  + saturation_excess_mm;
    res = NumericVector::create(_["Snowmelt"] = snowmelt_mm,
                                _["Source/sinks"] = sum(sourceSink),
                                _["Infiltration"] = infiltration_matrix_mm,
                                _["SaturationExcess"] = saturation_excess_mm,
                                _["Runoff"] = runoff_mm,
                                _["DeepDrainage"] = drainage_matrix_mm,
                                _["DrainageCorrection"] = correction_mm,
                                _["VolumeChange"] = Vfin_mm - Vini0_mm,
                                _["Substep"] = tsubstep);
  } 
  return(res);
}