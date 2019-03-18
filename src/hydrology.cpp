#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".er")]]
NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2){
  int nDays = DOY.size();
  NumericVector ER=rep(0.0,nDays);
  for(int i=0;i<nDays;i++){
    if((DOY[i]<=120)|(DOY[i]>=335)) {
      ER[i] = ERsyn;
    } else {
      ER[i] = ERconv;
    }
  }
  return(ER);
  
}


// [[Rcpp::export("hydrology.soilEvaporation")]]
double soilevaporation(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t+1)-sqrt(t)), PETs);
  return(Esoil);
}
// [[Rcpp::export(".infiltrationDay")]]
double infiltrationDay(double input, double Ssoil) {
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
// [[Rcpp::export("hydrology.infiltrationRepartition")]]
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


// [[Rcpp::export(".interceptionGashDay")]]
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
