#include <cmath>
#include "woodformation_c.h"

////// Effect of temperature (on metabolic rate and microtubule stability)
double _microT_c(double Tc, double inflection, double scale=5.0){
  double out = 1.0/(1.0+std::exp((-Tc+inflection)*scale));
  return(out);
}
double _metR_c(double Tc, double DHa, double DSd, double DHd){
  double Tk = Tc-T0;
  double out = Tk*std::exp(-DHa/(Rn*Tk)) / (1.0+std::exp(DSd/Rn*(1.0-(DHd/(DSd*Tk)))));
  return(out);
}

//' @rdname woodformation
//' @keywords internal
// [[Rcpp::export("woodformation_temperatureEffect")]]
double temperature_function_c(double Tc, double Y_T=5.0, double DHa=87.5e3, double DSd=1.09e3, double DHd=333e3){
  double out = _metR_c(Tc, DHa, DSd, DHd);
  out = out/_metR_c(30.0, DHa, DSd, DHd); // the output is equal to 1 at 30 degC
  out = out*_microT_c(Tc, Y_T);
  // out = 1;
  return(out);
}


//// Convert osmotic potential to osmolyte quantity and back
double _pi2n_c(double pi, double V, double Tc){
  double n = -pi*V/(Rn*(Tc-T0));
  return(n);
}

double _n2pi_c(double n, double V, double Tc){
  double pi = -n*Rn*(Tc-T0)/V;
  return(pi);
}

////// Cell expansion model
//' @rdname woodformation
//' @keywords internal
// [[Rcpp::export("woodformation_relativeExpansionRate")]]
double relative_expansion_rate_c(double psi, double Tc, double pi, double phi, double Y_P, double Y_T){
  double out = phi*(psi-pi-Y_P);
  if(out<0.0) out=0.0;
  out = out*temperature_function_c(Tc,Y_T);
  return(out);
}

////// Cell division model
double _divide_c(double psi, double Tc,
               double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
               double Y_P=0.05, double Y_T=5.0){
  // Default parameters from Cabon et al New Phytologist (2020)
  
  double r; //  Cell relative growth rate
  double P; // Cell production rate
  double pi_Tcorr = _n2pi_c(_pi2n_c(pi0,1.0,Tref),1.0,Tc);
  r = relative_expansion_rate_c(psi, Tc, pi_Tcorr, phi0, Y_P, Y_T);
  P = r/std::log(2.0)*Nc;
  
  return(P);
}

