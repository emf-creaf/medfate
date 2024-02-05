// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include "hydraulics.h"
#include <meteoland.h>
using namespace Rcpp;

//Old defaults
//ERconv=0.05, ERsyn = 0.2
//New defaults
//Rconv = 5.6, Rsyn = 1.5

//' @param doy Day of the year.
//' @param pet Potential evapotranspiration for a given day (mm).
//' @param prec Precipitation for a given day (mm).
//' @param Rconv,Rsyn Rainfall rate for convective storms and synoptic storms, respectively, in mm/h.
//' 
//' @rdname hydrology_interception
// [[Rcpp::export("hydrology_erFactor")]]
double erFactor(int doy, double pet, double prec, double Rconv = 5.6, double Rsyn = 1.5){
  double Ri = 0.0; //mm/h
  if((doy<=120) || (doy>=335)) {
    Ri = std::max(prec/24.0,Rsyn);
  } else {
    Ri = std::max(prec/24.0,Rconv);
  }
  double Ei =pet/24.0;
  return(Ei/Ri);
}

// [[Rcpp::export(".hydrology_interceptionGashDay")]]
double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05) {
  double I = 0.0;
  double PG = (-Cm/(ER*(1.0-p)))*log(1.0-ER); //Precipitation need to saturate the canopy
  if(Cm==0.0 || p==1.0) PG = 0.0; //Avoid NAs
  if(Precipitation>PG) {
    I = (1-p)*PG + (1-p)*ER*(Precipitation-PG);
  } else {
    I = (1-p)*Precipitation;
  }
  return(I);
}


//' @rdname hydrology_soil
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

//' @rdname hydrology_soil
//' 
//' @param soil An object of class \code{\link{soil}}.
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param pet Potential evapotranspiration for a given day (mm)
//' @param LgroundSWR Percentage of short-wave radiation (SWR) reaching the ground.
//' @param modifySoil Boolean flag to indicate that the input \code{soil} object should be modified during the simulation.
//' 
// [[Rcpp::export("hydrology_soilEvaporation")]]
NumericVector soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                              bool modifySoil = true) {
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector psiSoil = psi(soil, soilFunctions);
  int nlayers = W.size();
  NumericVector EsoilVec(nlayers,0.0);
  double swe = soil["SWE"]; //snow pack
  if(swe == 0.0) {
    double PETsoil = pet*(LgroundSWR/100.0);
    double Gsoil = soil["Gsoil"];
    double Ksoil = soil["Ksoil"];
    double Esoil = 0.0;
    // Allow evaporation only if water potential is less than -2 MPa
    if(psiSoil[0] > -2.0) Esoil = soilEvaporationAmount((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
    for(int l=0;l<nlayers;l++) {
      double cumAnt = 0.0;
      double cumPost = 0.0;
      for(int l2=0;l2<l;l2++) cumAnt +=dVec[l2];
      cumPost = cumAnt+dVec[l];
      //Exponential decay to divide bare soil evaporation among layers
      if(l<(nlayers-1)) EsoilVec[l] = Esoil*(exp(-Ksoil*cumAnt)-exp(-Ksoil*cumPost));
      else EsoilVec[l] = Esoil*exp(-Ksoil*cumAnt);
      if(modifySoil) W[l] = W[l] - ((EsoilVec[l])/Water_FC[l]);
    }
  }
  return(EsoilVec);
}

//' @rdname hydrology_soil
//' @param LherbSWR Percentage of short-wave radiation (SWR) reaching the herbaceous layer.
//' @param herbLAI Leaf area index of the herbaceous layer.
// [[Rcpp::export("hydrology_herbaceousTranspiration")]]
double herbaceousTranspiration(double pet, double LherbSWR, double herbLAI, 
                               List soil, String soilFunctions, bool modifySoil = true){
  if(NumericVector::is_na(herbLAI)) return(0.0);
  double Tmax_herb = pet*(LherbSWR/100.0)*(0.134*herbLAI - 0.006*pow(herbLAI, 2.0));
  NumericVector psiSoil = psi(soil, soilFunctions);
  NumericVector W = soil["W"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  double Eherb = Tmax_herb*Psi2K(psiSoil[0], -1.5, 2.0);
  if(modifySoil) {
    W[0] = W[0] - (Eherb/Water_FC[0]);
  }
  return(Eherb);
}

// [[Rcpp::export(".hydrology_infiltrationAmount")]]
double infiltrationAmount(double input, double Ssoil) {
  double I = 0;
  if(input>0.2*Ssoil) {
    I = input-(pow(input-0.2*Ssoil,2.0)/(input+0.8*Ssoil));
  } else {
    I = input;
  }
  return(I);
}

/**
 * Calculates infiltrated water that goes to each layer
 */
//' @rdname hydrology_soil
//' 
//' @param I Soil infiltration (in mm of water).
//' @param dVec Width of soil layers (in mm).
//' @param macro Macroporosity of soil layers (in \%).
//' @param a,b Parameters of the extinction function used for water infitration.
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

//' @rdname hydrology_soil
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
//' High-level functions for hydrological processes. Function \code{hydrology_soilWaterInputs} performs 
//' canopy water interception and snow accumulation/melt. Function \code{hydrology_soilInfiltrationPercolation} 
//' performs soil infiltration and percolation from the input given by the previous function.
//' 
//' @param soil A list containing the description of the soil (see \code{\link{soil}}).
//' @param soilFunctions Soil water retention curve and conductivity functions, either 'SX' (for Saxton) or 'VG' (for Van Genuchten).
//' @param prec Precipitation for a given day (mm)
//' @param er The ratio of evaporation rate to rainfall rate.
//' @param tday Average day temperature (ºC).
//' @param rad Solar radiation (in MJ/m2/day).
//' @param elevation Altitude above sea level (m).
//' @param Cm Canopy water storage capacity.
//' @param LgroundPAR Percentage of photosynthetically-acvive radiation (PAR) reaching the ground.
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
//' \item{NetRain}{Rainfall reaching the ground.}
//' \item{Snowmelt}{Snow melted during the day, and added to the water infiltrated.}
//' \item{Runon}{Surface water amount running on the target area from upslope.}
//' \item{Input}{Total soil input, including runon, snowmelt and net rain.}
//' 
//' Function \code{hydrology_soilInfiltrationPercolation} returns a named vector with the following elements, all in mm:
//' \item{Infiltration}{Water infiltrated into the soil (i.e. throughfall + runon + snowmelt - runoff).}
//' \item{Runoff}{Surface water leaving the target area.}
//' \item{DeepDrainage}{Water leaving the target soil towards the water table.}
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwb_day}}, \code{\link{hydrology_rainInterception}}, \code{\link{hydrology_soilEvaporation}}
//' 
//' @name hydrology_verticalInputs
// [[Rcpp::export("hydrology_soilWaterInputs")]]
NumericVector soilWaterInputs(List soil, String soilFunctions, double prec, double er, double tday, double rad, double elevation,
                             double Cm, double LgroundPAR, double LgroundSWR, 
                             double runon = 0.0,
                             bool snowpack = true, bool modifySoil = true) {
  //Soil input
  double swe = soil["SWE"]; //snow pack

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
    Interception = interceptionGashDay(rain,Cm,LgroundPAR/100.0,er);
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
                                           _["Input"] = runon+melt+NetRain);
  return(WI);
}

//' @rdname hydrology_verticalInputs
//' 
//' @param waterInput Soil water input for a given day (mm).
//' @param rockyLayerDrainage Boolean flag to indicate the simulation of drainage from rocky layers (> 95\% of rocks).
//' 
// [[Rcpp::export("hydrology_soilInfiltrationPercolation")]]
NumericVector soilInfiltrationPercolation(List soil, String soilFunctions, 
                                          double waterInput,
                                          bool rockyLayerDrainage = true, bool modifySoil = true) {
  //Soil input
  NumericVector W = clone(Rcpp::as<Rcpp::NumericVector>(soil["W"])); //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  double Kdrain = soil["Kdrain"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Water_SAT = waterSAT(soil, soilFunctions);
  int nlayers = W.size();
  

  //Hydrologic input
  double Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(waterInput>0.0) {
    //Interception
    //Net Runoff and infiltration
    Infiltration = infiltrationAmount(waterInput, Water_FC[0]);
    Runoff = waterInput - Infiltration;
    //Decide infiltration repartition among layers
    NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
    //Input of the first soil layer is infiltration
    double percolationExcess = 0.0;
    double Wn;
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0.0) && (Ivec[l]>0.0)) {
        Wn = W[l]*Water_FC[l] + Ivec[l]; //Update water volume
        if(l<(nlayers-1)) {
          Ivec[l+1] = Ivec[l+1] + std::max(Wn - Water_FC[l],0.0); //update Ivec adding the excess to the infiltrating water (saturated flow)
          W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta (this modifies 'soil')
        } else {
          if((rfc[l]<95.0) || rockyLayerDrainage) { //If not a rock layer or rocky layer drainage is allowed
            W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta (this modifies 'soil')
            percolationExcess = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer using field capacity
          } else {
            W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil')
            percolationExcess = std::max(Wn - Water_SAT[l],0.0); //Set excess of the bottom layer using saturation
          }
        }
      } 
    } 
    //If there still excess fill layers over field capacity
    if(percolationExcess>0.0) {
      for(int l=(nlayers-1);l>=0;l--) {
        if((dVec[l]>0.0) && (percolationExcess>0.0)) {
          Wn = W[l]*Water_FC[l] + percolationExcess; //Update water volume
          percolationExcess = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
          W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
        }
      }
      //If soil is completely saturated increase surface Runoff
      if(percolationExcess>0.0) { 
        Runoff = Runoff + percolationExcess;
      }
    }
  }
  //If there is still room for additional drainage (water in macropores accumulated from previous days)
  double head = 0.0;
  for(int l=0;l<nlayers;l++) { //Add mm over field capacity
    if((l<(nlayers-1)) || rockyLayerDrainage) {
      head += Water_FC[l]*std::max(W[l] - 1.0, 0.0);
    }
  }
  // Rcout<<head<<"\n";
  if(head>0.0) {
    double maxDrainage = head*Kdrain;
    for(int l=0;l<nlayers;l++) {
      if(maxDrainage>0.0) {
        double Wn = W[l]*Water_FC[l];
        double toDrain = std::min(std::max(Wn - Water_FC[l], 0.0), maxDrainage);
        if((l==(nlayers-1)) && (rfc[l] >= 95.0) && (!rockyLayerDrainage)) { //Prevent drainage for last rocky layer if not allowed
          toDrain = 0.0;
        }
        if(toDrain > 0.0) {
          DeepDrainage +=toDrain;
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
  NumericVector DB = NumericVector::create(_["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage);
  return(DB);
}


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


//' @rdname hydrology_soil
//' 
//' @param sourceSink Source/sink term for each soil layer (from snowmelt, soil evaporation or plant transpiration/redistribution)
//'        as mm/day.
//' @param nsteps  Number of time steps per day
//' 
// [[Rcpp::export("hydrology_soilFlows")]]
double soilFlows(List soil, NumericVector sourceSink, int nsteps = 24, bool modifySoil = true) {
  
  NumericVector dVec = soil["dVec"];
  double mm_day_2_m3_s = 0.001*(1.0/86400.0);//From mm/day = l/day = dm3/day to m3/s
  NumericVector sourceSink_m3s =  sourceSink*mm_day_2_m3_s;
  NumericVector dZ_m = dVec*0.001; //mm to m
  NumericVector rfc = soil["rfc"];

  double tstep = 86400.0/((double) nsteps);
  
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
  NumericVector Psi_m(nlayers), K_ms(nlayers), C_m(nlayers);
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
    //build tridiagonal terms
    for(int l=0;l<nlayers;l++) {
      if(l==0) { //first layer
        K_up = K_ms[l];
        K_down = 0.5*(K_ms[l] + K_ms[l+1]);
        a[l] = 0.0;
        c[l] = -1.0*K_down/dZDown[l];
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep) - c[l];
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep)*Psi_m[l]  - K_down - sourceSink_m3s[l];
      } else if(l<(nlayers - 1)) {
        K_up = 0.5*(K_ms[l-1] + K_ms[l]);
        K_down = 0.5*(K_ms[l] + K_ms[l+1]);
        a[l] = -1.0*K_up/dZUp[l];
        c[l] = -1.0*K_down/dZDown[l];
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep) - a[l] - c[l];
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep)*Psi_m[l] + K_up - K_down - sourceSink_m3s[l];
      } else { // last layer
        K_up = 0.5*(K_ms[l-1] + K_ms[l]);
        K_down = K_ms[l];
        a[l] = -1.0*K_up/dZUp[l];
        b[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep) - a[l];
        c[l] = 0.0;
        d[l] = (lambda[l]*C_m[l]*dZ_m[l]/tstep)*Psi_m[l] + K_up - K_down - sourceSink_m3s[l];
      }
    }
    NumericVector Psi_m_t1 = tridiagonalSolving(a,b,c,d);
    NumericVector Psi_t1 = Psi_m_t1*mTOMPa; // m to MPa
    //calculate free drainage (m3)
    drainage += std::max(0.0, K_ms[nlayers -1]*tstep);
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