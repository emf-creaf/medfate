#include <Rcpp.h>
using namespace Rcpp;

/// Constants
double Rn = 8.314; // The perfect gas constant
double T0 = -273.15; // Absolute 0 temperature in degC
double Tref = 15; // Reference temperature in degC


List initialize_ring(double phi0=0.13, double pi0=-0.8, double CRD0=8.3){
  
  IntegerVector formation_dog(0);
  int l = formation_dog.size();
  NumericVector phi(l, phi0);
  NumericVector pi(l, pi0);
  NumericVector CRD(l, CRD0);
  
  IntegerVector dog(0);
  NumericVector P;
  
  DataFrame cells = DataFrame::create(_["formation_dog"] = formation_dog,
                                      _["phi"] = phi,
                                      _["pi"] = pi,
                                      _["CRD"] = CRD);
  
  DataFrame divisions = DataFrame::create(_["dog"] = dog,
                                          _["P"] = P);
  
  List ring = List::create(_["divisions"] = divisions,
                           _["cells"] = cells);
  
  
  return ring;
}

////// Effect of temperature (on metabolic rate and microtubule stability)
double _microT(double Tc, double inflection, double scale=5){
  double out = 1/(1+exp((-Tc+inflection)*scale));
  return out;
}
double _metR(double Tc, double DHa, double DSd, double DHd){
  double Tk = Tc-T0;
  double out = Tk*exp(-DHa/(Rn*Tk)) / (1+exp(DSd/Rn*(1-(DHd/(DSd*Tk)))));
  return out;
}
double T_fun(double Tc, double Y_T=8, double DHa=87.5e3, double DSd=1.09e3, double DHd=333e3){
  double out = _metR(Tc, DHa, DSd, DHd);
  out = out/_metR(Tref, DHa, DSd, DHd); // the output is equal to 1 at Tref degC
  out = out*_microT(Tc, Y_T);
  // out = 1;
  return out;
}


//// Convert osmotic potential to osmolyte quantity and back
double _pi2n(double pi, double V, double Tc){
  double n = -pi*V/(Rn*(Tc-T0));
  return n;
}

double _n2pi(double n, double V, double Tc){
  double pi = -n*Rn*(Tc-T0)/V;
  return pi;
}

////// Cell expansion model
double _r(double psi, double Tc, double pi, double phi, double Y_P, double Y_T){
  double out = phi*(psi-pi-Y_P);
  if(out<0) out=0;
  out = out*T_fun(Tc,Y_T);
  return out;
}

////// Cell division model
double _divide(double psi, double Tc,
               double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
               double Y_P=0.05, double Y_T=8){
  // Default parameters from Cabon et al New Phytologist (2020)
  
  double r; //  Cell relative growth rate
  double P; // Cell production rate
  double pi_Tcorr = _n2pi(_pi2n(pi0,1,Tref),1,Tc);
  r = _r(psi, Tc, pi_Tcorr, phi0, Y_P, Y_T);
  P = r/log(2)*Nc;
  
  return(P);
}
DataFrame _expand_cell(double psi, double Tc,
                  double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                  double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){
  // default parameters from Cabon et al. New Phytologist 2020. h is different because of potential error in the ref value.
  
  double n = _pi2n(pi0, CRD0, Tref);
  pi0 = _n2pi(n, CRD0, Tc); // updates the value of pi0 which is given at Tref for the current temperature Tc
  
  // Calculate relative volume expansion rate
  double r = _r(psi, Tc, pi0, phi0, Y_P, Y_T);
  
  // Variable update
  double CRD1 = CRD0*(1+r); // cell diameter (volume) increment
  double pi1 = _n2pi(n, CRD1, Tref); // pi is returned at Tref in order to be consistent with input
  double phi1 = phi0 + phi0*(s*r - h*T_fun(Tc, -999)); //changes in cell wall properties. Hardening (thickening and lignification) is temperature sensitive but not threshold prone because lignification does not need microtubules
  if(phi1<0) {
    phi1=0;
  } else {}
  
  // return outputs
  return(DataFrame::create(_["phi"]=phi1,
                           _["pi"]=pi1,
                           _["CRD"]=CRD1));
}


void _expand_ring(List ring, double psi, double Tc, 
                  double Y_P=0.05, double Y_T=5, double h=0.043*1.8, double s=1.8){
  
  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  IntegerVector formation_dog = cells["formation_dog"];
  int l = cells.nrow();
  
  for(int i=0; i<l; i++){
    DataFrame temp = _expand_cell(psi, Tc,
                             phi[i], pi[i], CRD[i],
                             Y_P, Y_T, h, s);
    
    phi[i] = temp["phi"];
    pi[i] = temp["pi"];
    CRD[i] = temp["CRD"];
  }
  
  cells["phi"] = phi;
  cells["pi"] = pi;
  cells["CRD"] = CRD;
  
  ring["cells"] = cells;
  // return ring;
}

// [[Rcpp::export]]
void grow_ring(List ring, double psi, double Tc,
                double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                double Y_P=0.05, double Y_T=8, double h=0.043*1.8, double s=1.8){
  
  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  DateVector formation_dog = cells["formation_dog"];
  
  DataFrame divisions = as<DataFrame>(ring["divisions"]);
  NumericVector P = divisions["P"];
  IntegerVector P_dog = divisions["dog"];
  
  int dog = P_dog.size()+1;
  
  // Calculate cell production
  double P_i = _divide(psi, Tc, Nc, phi0, pi0, Y_P, Y_T);
  // Update the cell production vector
  P.push_back(P_i);
  P_dog.push_back(dog);
  divisions = DataFrame::create(_["dog"] = P_dog,
                                _["P"] = P);
  
  
  // Calculate the whole number of cells formed at the current and previous timestep
  double Pnew = sum(P); int Pnew_int = floor(Pnew);
  double Pold = Pnew-P_i; int Pold_int = floor(Pold);
  // Add a new value in the cell expansion vectors for each whole cell number increment
  for (int i=0; i < Pnew_int-Pold_int; i++){
    // Update cell expansion vectors
    phi.push_back(phi0);
    pi.push_back(pi0);
    CRD.push_back(CRD0);
    formation_dog.push_back(dog);
  }
  
  // Create new cells data frame
  cells = DataFrame::create(_["formation_dog"] = formation_dog,
                            _["phi"] = phi,
                            _["pi"] = pi,
                            _["CRD"] = CRD);
  // Update ring object
  ring["divisions"] = divisions;
  ring["cells"] = cells;
  // Calculate cell expansion
  _expand_ring(ring, psi, Tc, Y_P, Y_T, h, s);
}
