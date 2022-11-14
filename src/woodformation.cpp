#include <Rcpp.h>
using namespace Rcpp;

/// Constants
double Rn = 8.314; // The perfect gas constant
double T0 = -273.15; // Absolute 0 temperature in degC
double Tref = 15.0; // Reference temperature in degC

//' Wood formation
//' 
//' Functions to initialize and expand a ring of tracheids to simulate secondary growth.
//' 
//' @param ring An object of class \code{\link{ring}} returned by function \code{woodformation_initRing}.
//' @param psi Water potential (in MPa).
//' @param Tc Temperature in Celsius.
//' @param Nc Number of active cells in the cambium.
//' @param phi0 Initial value of cell extensibility (in MPa-1 day-1)
//' @param pi0 Initial value of cell osmotic potential (in MPa)
//' @param CRD0 Initial value of cell radial diameter
//' @param Y_P Turgor pressure yield threshold (in MPa)
//' @param Y_T Temperature yield threshold (in Celsius)
//' @param h Cell wall hardening coefficient (in day-1)
//' @param s Cell wall softening coefficient (unitless)
//' @param pi Osmotic potential (in MPa)
//' @param phi Cell extensibility (in MPa-1 day-1)
//' @param DHa,DSd,DHd Enthalpy of activation, enthalpy difference and entropy difference between the catalytically active and inactive states of the enzymatic system (Parent et al. 2010).
//'
//' @return 
//' Function \code{woodformation_initRing()} returns a list of class 'ring', 
//' that is a list containing a data frame \code{cells} and two vectors: \code{P} and \code{SA}. 
//' Dataframe \code{cells} contains the columns "formation_date", "phi", "pi" and "CRD" and as many rows as dates processed. 
//' Vectors \code{P} and \code{SA} contain, respectively, the number of cells produced and the sapwood area 
//' corresponding to the ring of cells (assuming a tangencial radius of 20 micrometers). 
//' 
//' Function \code{woodformation_growRing()} modifies the input 'ring' object according to the environmental conditions given as input.
//' 
//' Function \code{woodformation_relativeExpansionRate()} returns a numeric scalar with the relative expansion rate. 
//' 
//' Function \code{woodformation_temperatureEffect()} returns a scalar between 0 and 1 reflecting the temperature effect on tissue formation rate. 
//' 
//' Function \code{woodformation_relativeGrowthRate} returns the annual growth rate, relative to cambium perimeter, estimated from initial and final diameter values.
//' 
//' @note Code modified from package xylomod by Antoine Cabon, available at GitHub
//' 
//' @author
//' Antoine Cabon, CTFC
//' 
//' Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Cabon A, \enc{Fernández-de-Uña}{Fernandez-de-Una} L, Gea-Izquierdo G, Meinzer FC, Woodruff DR, \enc{Martínez-Vilalta}{Martinez-Vilalta} J, De \enc{Cáceres}{Caceres} M. 2020a. Water potential control of turgor-driven tracheid enlargement in Scots pine at its xeric distribution edge. New Phytologist 225: 209–221.
//' 
//' Cabon A, Peters RL, Fonti P, \enc{Martínez-Vilalta}{Martinez-Vilalta}  J, De \enc{Cáceres}{Caceres} M. 2020b. Temperature and water potential co-limit stem cambial activity along a steep elevational gradient. New Phytologist: nph.16456.
//' 
//' Parent, B., O. Turc, Y. Gibon, M. Stitt, and F. Tardieu. 2010. Modelling temperature-compensated physiological rates, based on the co-ordination of responses to temperature of developmental processes. Journal of Experimental Botany 61:2057–2069.
//' 
//' @seealso \code{\link{growth}}
//' 
//' @name woodformation
// [[Rcpp::export("woodformation_initRing")]]
List initialize_ring(){
  
  IntegerVector formation(0);
  NumericVector phi(0);
  NumericVector pi(0);
  NumericVector CRD(0);
  
  IntegerVector dog(0);
  NumericVector P(0),SA(0);
  
  DataFrame cells = DataFrame::create(_["formation"] = formation,
                                      _["phi"] = phi,
                                      _["pi"] = pi,
                                      _["CRD"] = CRD);
  
  List ring = List::create(_["P"] = P,
                           _["SA"] = SA,
                           _["cells"] = cells);
  
  
  return(ring);
}

////// Effect of temperature (on metabolic rate and microtubule stability)
double _microT(double Tc, double inflection, double scale=5.0){
  double out = 1.0/(1.0+exp((-Tc+inflection)*scale));
  return(out);
}
double _metR(double Tc, double DHa, double DSd, double DHd){
  double Tk = Tc-T0;
  double out = Tk*exp(-DHa/(Rn*Tk)) / (1.0+exp(DSd/Rn*(1.0-(DHd/(DSd*Tk)))));
  return(out);
}

//' @rdname woodformation
// [[Rcpp::export("woodformation_temperatureEffect")]]
double temperature_function(double Tc, double Y_T=5.0, double DHa=87.5e3, double DSd=1.09e3, double DHd=333e3){
  double out = _metR(Tc, DHa, DSd, DHd);
  out = out/_metR(30.0, DHa, DSd, DHd); // the output is equal to 1 at 30 degC
  out = out*_microT(Tc, Y_T);
  // out = 1;
  return(out);
}


//// Convert osmotic potential to osmolyte quantity and back
double _pi2n(double pi, double V, double Tc){
  double n = -pi*V/(Rn*(Tc-T0));
  return(n);
}

double _n2pi(double n, double V, double Tc){
  double pi = -n*Rn*(Tc-T0)/V;
  return(pi);
}

////// Cell expansion model
//' @rdname woodformation
// [[Rcpp::export("woodformation_relativeExpansionRate")]]
double relative_expansion_rate(double psi, double Tc, double pi, double phi, double Y_P, double Y_T){
  double out = phi*(psi-pi-Y_P);
  if(out<0.0) out=0.0;
  out = out*temperature_function(Tc,Y_T);
  return(out);
}

////// Cell division model
double _divide(double psi, double Tc,
               double Nc = 8.85, double phi0=0.13, double pi0=-0.8,
               double Y_P=0.05, double Y_T=5.0){
  // Default parameters from Cabon et al New Phytologist (2020)
  
  double r; //  Cell relative growth rate
  double P; // Cell production rate
  double pi_Tcorr = _n2pi(_pi2n(pi0,1.0,Tref),1.0,Tc);
  r = relative_expansion_rate(psi, Tc, pi_Tcorr, phi0, Y_P, Y_T);
  P = r/log(2.0)*Nc;
  
  return(P);
}
List _expand_cell(double psi, double Tc,
                  double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                  double Y_P=0.05, double Y_T=5.0, double h=0.043*1.8, double s=1.8){
  // default parameters from Cabon et al. New Phytologist 2020. h is different because of potential error in the ref value.
  
  double n = _pi2n(pi0, CRD0, Tref);
  pi0 = _n2pi(n, CRD0, Tc); // updates the value of pi0 which is given at Tref for the current temperature Tc
  
  // Calculate relative volume expansion rate
  double r = relative_expansion_rate(psi, Tc, pi0, phi0, Y_P, Y_T);
  
  // Variable update
  double CRD1 = CRD0*(1.0+r); // cell diameter (volume) increment
  double pi1 = _n2pi(n, CRD1, Tref); // pi is returned at Tref in order to be consistent with input
  double phi1 = phi0 + phi0*(s*r - h*(_metR(Tc, 87.5e3, 1.09e3, 333e3)/_metR(Tref, 87.5e3, 1.09e3, 333e3))); //changes in cell wall properties. Hardening (thickening and lignification) is temperature sensitive but not threshold prone because lignification does not need microtubules
  if(phi1<0.0) {
    phi1=0.0;
  }
  
  // return outputs
  return(List::create(_["phi"]=phi1,
                      _["pi"]=pi1,
                      _["CRD"]=CRD1));
}


void _expand_ring(List ring, double psi, double Tc, 
                  double Y_P=0.05, double Y_T=5.0, double h=0.043*1.8, double s=1.8){
  
  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  IntegerVector formation = cells["formation"];
  int l = cells.nrow();
  
  for(int i=0; i<l; i++){
    List temp = _expand_cell(psi, Tc,
                             phi[i], pi[i], CRD[i],
                             Y_P, Y_T, h, s);
    
    phi[i] = temp["phi"];
    pi[i] = temp["pi"];
    CRD[i] = temp["CRD"];
  }
}

//' @rdname woodformation
// [[Rcpp::export("woodformation_growRing")]]
void grow_ring(List ring, double psi, double Tc,
                double Nc=8.85, double phi0=0.13, double pi0=-0.8, double CRD0=8.3,
                double Y_P=0.05, double Y_T=5.0, double h=0.043*1.8, double s=1.8){
  
  DataFrame cells = as<DataFrame>(ring["cells"]);
  NumericVector phi = cells["phi"];
  NumericVector pi = cells["pi"];
  NumericVector CRD = cells["CRD"];
  IntegerVector formation = cells["formation"];
  
  NumericVector P = as<NumericVector>(ring["P"]);
  NumericVector SA = as<NumericVector>(ring["SA"]);
  
  int dog = P.size()+1;
  
  // Calculate cell production
  double P_i = _divide(psi, Tc, Nc, phi0, pi0, Y_P, Y_T);
  // Update the cell production vector
  P.push_back(P_i);

  
  // Calculate the whole number of cells formed at the current and previous timestep
  double Pnew = sum(P); int Pnew_int = floor(Pnew);
  double Pold = Pnew-P_i; int Pold_int = floor(Pold);
  // Add a new value in the cell expansion vectors for each whole cell number increment
  for (int i=0; i < Pnew_int-Pold_int; i++){
    // Update cell expansion vectors
    phi.push_back(phi0);
    pi.push_back(pi0);
    CRD.push_back(CRD0);
    formation.push_back(dog);
  }
  // double SAprev = 0.0;
  // if(SA.length()>0) SAprev = SA[SA.length()-1];
  double SAnew = 0.0;
  for(int i=0;i<CRD.length();i++) {
    double w = std::min(1.0, ((double)(dog-formation[i]))/20.0);
    SAnew += w*CRD[i]*20.0; //20 micras diametro tangencial 
  }
  SA.push_back(SAnew);
  
  // Create new cells data frame
  cells = DataFrame::create(_["formation"] = formation,
                            _["phi"] = phi,
                            _["pi"] = pi,
                            _["CRD"] = CRD);
  // Update ring object
  ring["P"] = P;
  ring["SA"] = SA;
  ring["cells"] = cells;
  // Calculate cell expansion
  _expand_ring(ring, psi, Tc, Y_P, Y_T, h, s);
}

