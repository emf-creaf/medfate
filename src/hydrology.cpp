// [[Rcpp::interfaces(r,cpp)]]

#include <Rcpp.h>
#include "soil.h"
#include <meteoland.h>
using namespace Rcpp;
//Old defaults
//ERconv=0.05, ERsyn = 0.2
//New defaults
//Rconv = 5.6, Rsyn = 1.5
// [[Rcpp::export("hydrology_erFactor")]]
double erFactor(int doy, double pet, double prec, double Rconv = 5.6, double Rsyn = 1.5){
  double Ri = 0.0; //mm/h
  if((doy<=120)|(doy>=335)) {
    Ri = std::max(prec/24.0,Rsyn);
  } else {
    Ri = std::max(prec/24.0,Rconv);
  }
  double Ei =pet/24.0;
  return(Ei/Ri);
}

// [[Rcpp::export("hydrology_soilEvaporationAmount")]]
double soilEvaporationAmount(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t+1)-sqrt(t)), PETs);
  return(Esoil);
}

// [[Rcpp::export("hydrology_soilEvaporation")]]
NumericVector soilEvaporation(List soil, String soilFunctions, double pet, double LgroundSWR,
                              bool modifySoil = true) {
  NumericVector W = soil["W"]; //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  int nlayers = W.size();
  NumericVector EsoilVec(nlayers,0.0);
  double swe = soil["SWE"]; //snow pack
  if(swe == 0.0) {
    double PETsoil = pet*LgroundSWR;
    double Gsoil = soil["Gsoil"];
    double Ksoil = soil["Ksoil"];
    double Esoil = soilEvaporationAmount((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
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
// [[Rcpp::export("hydrology_infiltrationRepartition")]]
NumericVector infiltrationRepartition(double I, NumericVector dVec, NumericVector macro) {
  int nlayers = dVec.length();
  NumericVector Ivec = NumericVector(nlayers,0.0);
  double macro_mean = sum(macro)/((double) nlayers);
  double a = -0.01*pow(1.0-macro_mean,3.0);
  double z1 = 0.0;
  for(int i=0;i<nlayers;i++) {
    if(i<(nlayers-1)) {
      Ivec[i] = I*(exp(a*z1)-exp(a*(z1 + dVec[i])));
    } else {
      Ivec[i] = I*exp(a*z1);
    }
    z1 = z1 + dVec[i];
  }
  return(Ivec);
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

// [[Rcpp::export("hydrology_verticalInputs")]]
NumericVector verticalInputs(List soil, String soilFunctions, double prec, double er, double tday, double rad, double elevation,
                    double Cm, double LgroundPAR, double LgroundSWR, 
                    double runon = 0.0,
                    bool snowpack = true, bool drainage = true, bool modifySoil = true) {
  //Soil input
  NumericVector W = clone(Rcpp::as<Rcpp::NumericVector>(soil["W"])); //Access to soil state variable
  NumericVector dVec = soil["dVec"];
  NumericVector macro = soil["macro"];
  NumericVector rfc = soil["rfc"];
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Water_SAT = waterSAT(soil, soilFunctions);
  double swe = soil["SWE"]; //snow pack
  int nlayers = W.size();
  
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
      if(NumericVector::is_na(rad)) stop("Missing radiation data for snow melt!");
      if(NumericVector::is_na(elevation)) stop("Missing elevation data for snow melt!");
      double rho = meteoland::utils_airDensity(tday, meteoland::utils_atmosphericPressure(elevation));
      double ten = (86400.0*tday*rho*1013.86*1e-6/100.0); //ten can be negative if temperature is below zero
      double ren = (rad*LgroundSWR)*(0.1); //90% albedo of snow
      melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
      // Rcout<<" swe: "<< swe<<" temp: "<<ten<< " rad: "<< ren << " melt : "<< melt<<"\n";
      swe = std::max(0.0, swe-melt);
    }
  } else {
    rain = prec;
  }
  
  //Hydrologic input
  double NetRain = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  double interception = interceptionGashDay(rain,Cm,LgroundPAR,er);
  if(rain>0.0) NetRain = rain - interception;
  if((NetRain+runon+melt)>0.0) {
    //Interception
    //Net Runoff and infiltration
    Infiltration = infiltrationAmount(NetRain+runon+melt, Water_FC[0]);
    Runoff = (NetRain+runon+melt) - Infiltration;
    //Decide infiltration repartition among layers
    NumericVector Ivec = infiltrationRepartition(Infiltration, dVec, macro);
    //Input of the first soil layer is infiltration
    double excess = 0.0;
    double Wn;
    //Update topsoil layer
    for(int l=0;l<nlayers;l++) {
      if((dVec[l]>0.0) & (Ivec[l]>0.0)) {
        Wn = W[l]*Water_FC[l] + Ivec[l]; //Update water volume
        if(l<(nlayers-1)) {
          Ivec[l+1] = Ivec[l+1] + std::max(Wn - Water_FC[l],0.0); //update Ivec adding the excess to the infiltrating water (saturated flow)
        } else {
          excess = std::max(Wn - Water_FC[l],0.0); //Set excess of the bottom layer
        }
        W[l] = std::max(0.0,std::min(Wn, Water_FC[l])/Water_FC[l]); //Update theta (this modifies 'soil')
      } 
    }
    if(drainage) {//Set deep drainage
      DeepDrainage = excess; 
    } else { //Fill to saturation and upwards if needed
      for(int l=(nlayers-1);l>=0;l--) {
        if((dVec[l]>0.0) & (excess>0.0)) {
          Wn = W[l]*Water_FC[l] + excess; //Update water volume
          excess = std::max(Wn - Water_SAT[l],0.0); //Update excess, using the excess of water over saturation
          W[l] = std::max(0.0,std::min(Wn, Water_SAT[l])/Water_FC[l]); //Update theta (this modifies 'soil') here no upper
        }
      }
      if(excess>0.0) { //If soil is completely saturated increase Runoff
        Runoff = Runoff + excess;
      }
    }
  }
  if(modifySoil) {
    soil["SWE"] = swe;
    soil["W"] = W;
  }
  NumericVector DB = NumericVector::create(_["Rain"] = rain, _["Snow"] = snow,
                                           _["Interception"] = interception,
                                           _["Throughfall"] = NetRain, 
                                           _["Snowmelt"] = melt,
                                           _["Runon"] = runon, 
                                           _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage);
  return(DB);
}