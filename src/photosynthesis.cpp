#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils.h"
#include <meteoland.h>
using namespace Rcpp;

const double R_gas = 8.314; //(J/mol/ºK) Universal gas constant
const double O2_conc = 209.0; //mmol*mol-1 (Collatz et al. 2001)
const double quantumYield = 0.3; //mol photon * mol-1 electron
const double lightResponseCurvature = 0.9;




/**
 * Species-independent photosynthesis terms depenent of leaf temperature
 *  From:
 *  Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, and S. P. Long. 2001. Improved temperature response functions for models of Rubisco-limited photosynthesis. 
 *  Plant, Cell and Environment 24:253–259.
 *  
 *  Tleaf - Leaf temperature (ºC)
 *  Oi - Oxigen concentration (mmol*mol-1)
 */
//Compensation point (micromol * mol-1)
// [[Rcpp::export("photo_GammaTemp")]]
double gammaTemp(double Tleaf) {return(42.75*exp((37830*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));} 
//Michaelis-Menten coefficients of Rubisco for Carbon (micromol * mol-1)
double KcTemp(double Tleaf) {return(404.9*exp((79430*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}
//Michaelis-Menten coefficients of Rubisco for Oxigen (mmol * mol-1)
double KoTemp(double Tleaf) {return(278.4*exp((36380*(Tleaf-25.0))/(298.0*R_gas*(Tleaf+273))));}
// [[Rcpp::export("photo_KmTemp")]]
double KmTemp(double Tleaf, double Oi = 209.0) {
  double Kc = KcTemp(Tleaf);
  double Ko = KoTemp(Tleaf);  
  return(Kc*(1.0+(Oi/Ko)));
}

/**
 * Temperature correction of Vmax298
 * 
 * From:
 *    Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. 
 *    Plant, Cell and Environment 25:1205–1210.
 *    
 * Eq.1 with parameters from Table 2
 * 
 *  Tleaf - Leaf temperature (ºC)
 *  Vmax298 - maximum carboxylation rate at 298ºK (ie. 25 ºC) (micromol*s-1*m-2)
 */
// [[Rcpp::export("photo_VmaxTemp")]]
double VmaxTemp(double Vmax298, double Tleaf) {
  double Ha = 73637.0; //Energy of activation J * mol-1
  double Hd = 149252.0; //Energy of deactivation J * mol-1
  double Sv = 486.0;  //Entropy term J * mol-1 * K-1
  double C = 1.0+exp((Sv*298.2-Hd)/(R_gas*298.2));
  return(Vmax298*(C*exp((Ha/(R_gas*298.2))*(1.0-298.2/(Tleaf+273.2))))/(1.0+exp((Sv*Tleaf-Hd)/(R_gas*(Tleaf+273.2)))));
}

/**
 * Temperature correction of Jmax298
 * 
 * From:
 *    Leuning, R. 2002. Temperature dependence of two parameters in a photosynthesis model. 
 *    Plant, Cell and Environment 25:1205–1210.
 *    
 * Eq.1 with parameters from Table 2
 * 
 *  Tleaf - Leaf temperature (ºC)
 *  Jmax298 - maximum electron transport rate at 298ºK (ie. 25 ºC) (micromol*s-1*m-2)
 */
// [[Rcpp::export("photo_JmaxTemp")]]
double JmaxTemp(double Jmax298, double Tleaf) {
  double Ha = 50300.0; //Energy of activation J * mol-1
  double Hd = 152044.0; //Energy of deactivation J * mol-1
  double Sv = 495.0;  //Entropy term J * mol-1 * K-1
  double C = 1.0+exp((Sv*298.2-Hd)/(R_gas*298.2));
  return(Jmax298*(C*exp((Ha/(R_gas*298.2))*(1.0-298.2/(Tleaf+273.2))))/(1.0+exp((Sv*Tleaf-Hd)/(R_gas*(Tleaf+273.2)))));
}


/**
 * Calculates electron-limited photosynthesis (Farquhar et al. 1980)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Jmax - maximum electron transport rate per leaf area (micromol*s-1*m-2)
 * 
 * return units: micromol*s-1*m-2
 */
// [[Rcpp::export("photo_electronLimitedPhotosynthesis")]]
double electronLimitedPhotosynthesis(double Q, double Ci, double GT, double Jmax) {
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);
  return((J/4.0)*((Ci-GT)/(Ci+2.0*GT)));
}
double electronLimitedPhotosynthesisDerivative(double Q, double Ci, double GT, double Jmax){
  double J = ((quantumYield*Q+Jmax)-sqrt(pow(quantumYield*Q+Jmax, 2.0) - 4.0*lightResponseCurvature*quantumYield*Q*Jmax))/(2.0*lightResponseCurvature);
  return((J/4.0)*((3.0*GT)/pow(Ci+2.0*GT,2.0)));
}

/**
 * Calculates rubisco-limited photosynthesis (Farquhar et al. 1980)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1)
 * Vmax - maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2)
 * 
 * return units: micromol*s-1*m-2
 */
// [[Rcpp::export("photo_rubiscoLimitedPhotosynthesis")]]
double rubiscoLimitedPhotosynthesis(double Ci, double GT, double Km, double Vmax) {
  return(Vmax *(Ci-GT)/(Ci+Km));
}
double rubiscoLimitedPhotosynthesisDerivative(double Ci, double GT, double Km, double Vmax) {
  return(Vmax *(Km+GT)/pow(Ci+Km,2.0));
}
/**
 * Calculates photosynthesis (Farquhar et al. 1980/Collatz et al 1991)
 * 
 * Ci - CO2 internal concentration (micromol * mol-1)
 * GT - CO2 saturation point corrected by temperature (micromol * mol-1)
 * Km = Kc*(1.0+(Oi/Ko)) - Michaelis-Menten term corrected by temperature (in micromol * mol-1)
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Jmax - maximum electron transport rate per leaf area (micromol*s-1*m-2) corrected  by temperature
 * Vmax - maximum Rubisco carboxylation rate per leaf area (micromol*s-1*m-2) corrected  by temperature
 * 
 * return units: micromol*s-1*m-2
 */
double photosynthesis_Ci(double Q, double Ci, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis(Q, Ci, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis(Ci, GT, Km, Vmax);
  return(std::max(0.0,(Je+Jc-sqrt(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc))/(2.0*0.98)));
}

// Auxiliary functions for Newton-Raphson
double f(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  return(photosynthesis_Ci(Q,x, GT, Km, Vmax, Jmax)-(Gc*(Ca-x)));
}
double fder(double x, double Q, double Ca, double Gc, double GT, double Km, double Vmax, double Jmax) {
  double Je = electronLimitedPhotosynthesis(Q, x, GT, Jmax);
  double dJe = electronLimitedPhotosynthesisDerivative(Q, x, GT, Jmax);
  double Jc = rubiscoLimitedPhotosynthesis(x, GT, Km, Vmax);
  double dJc = rubiscoLimitedPhotosynthesisDerivative(x, GT, Km, Vmax);
  double dA1 = (1.0/(2.0*0.98))*(dJe+dJc-(0.5*pow(pow(Je+Jc,2.0)-4.0*0.98*Je*Jc,-0.5)*(2.0*Je*dJe+2.0*Jc*dJc+(2.0-4.0*0.98)*(dJe*Jc + dJc*Je))));
  double dA2 = -Gc;
  return(dA1-dA2);
}

/**
 * Calculates photosynthesis (Farquhar et al. 1980/Collatz et al 1991)
 * 
 * Q - Active photon flux density (micromol * s-1 * m-2)
 * Ca - CO2 air concentration (micromol * mol-1)
 * Gc - CO2 stomatal conductance (mol * s-1 * m-2)
 * Tleaf - Leaf temperature (ºC)
 * Jmax298 - maximum electron transport rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2) 
 * Vmax298 - maximum Rubisco carboxylation rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2) 
 * 
 * return units: micromol*s-1*m-2
 */
// [[Rcpp::export("photo_photosynthesis")]]
NumericVector leafphotosynthesis(double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298, bool verbose=false) {
  //Corrections per leaf temperature
  double GT = gammaTemp(Tleaf);
  double Km = KmTemp(Tleaf, O2_conc);
  double Vmax = VmaxTemp(Vmax298, Tleaf);
  double Jmax = JmaxTemp(Jmax298, Tleaf);
  double x,x1,e,fx,fx1;
  x1 = 0.0;//initial guess
  e = 0.001; // accuracy in micromol * mol-1
  int cnt = 0;
  int mxiter = 100;
  if(verbose) Rcout <<"x{i}"<<"    "<<"x{i+1}"<<"        "<<"|x{i+1}-x{i}|\n";                
  do {
    x=x1; /*make x equal to the last calculated value of                             x1*/
    fx=f(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying f(x)to fx
    fx1=fder(x, Q, Catm, Gc, GT, Km, Vmax, Jmax);            //simplifying fprime(x) to fx1
    x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
    cnt++;
    if(verbose) Rcout<<x<<"     "<<x1<<"           "<<std::abs(x1-x)<<"\n";        
  } while ((std::abs(x1-x)>=e) & (cnt < mxiter));
  double A = photosynthesis_Ci(Q,x1,GT,Km,Vmax,Jmax);
  NumericVector res = NumericVector::create(x1, A);
  res.attr("names") = CharacterVector::create("Ci", "A");
  return(res);
}


// [[Rcpp::export("photo_leafPhotosynthesisFunction")]]
DataFrame leafPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                             double absRad, double Q, double Vmax298, double Jmax298, 
                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  double Gwdiff, Gbw;
  for(int i=0;i<nsteps;i++){
    leafTemp[i] = leafTemperature(absRad/refLeafArea, Tair, u, E[i], leafWidth);
    leafVPD[i] = std::max(0.0,leafVapourPressure(leafTemp[i], psiLeaf[i]) - vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw[i]  =  Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    // NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gw[i]/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    Gbw = 0.397*pow(u/(leafWidth*0.0072), 0.5); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
    Gsw[i]  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    Ci[i] = LP[0];
    Ag[i] = LP[1];
    An[i] = Ag[i] - 0.015*VmaxTemp(Vmax298/refLeafArea, leafTemp[i]);
  }
  return(DataFrame::create(Named("LeafTemperature") = leafTemp,
                      Named("LeafVPD") = leafVPD,
                      Named("Gsw") = Gsw,
                      Named("Ci") = Ci,
                      Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}
// [[Rcpp::export("photo_leafPhotosynthesisFunction2")]]
DataFrame leafPhotosynthesisFunction2(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                     double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                     double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  double Gwdiff, Gbw;
  // Rcout <<refLeafArea<<" "<<SWRabs/refLeafArea<< " "<< LWRnet/refLeafArea << "  "<<Tair<< " "<<u<<"\n";
  for(int i=0;i<nsteps;i++){
    leafTemp[i] = leafTemperature2(SWRabs/refLeafArea, LWRnet/refLeafArea, Tair, u, E[i], leafWidth);
    leafVPD[i] = std::max(0.0,leafVapourPressure(leafTemp[i], psiLeaf[i]) - vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw[i]  =  Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    // NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gw[i]/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    Gbw = 0.397*pow(u/(leafWidth*0.0072), 0.5); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
    Gsw[i]  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    Ci[i] = LP[0];
    Ag[i] = LP[1];
    An[i] = Ag[i] - 0.015*VmaxTemp(Vmax298/refLeafArea, leafTemp[i]);
  }
  return(DataFrame::create(Named("LeafTemperature") = leafTemp,
                           Named("LeafVPD") = leafVPD,
                           Named("Gsw") = Gsw,
                           Named("Ci") = Ci,
                           Named("GrossPhotosynthesis") = Ag,
                           Named("NetPhotosynthesis") = An));
}


/**
 * Calculates gross/net canopy photosynthesis function, considering a multilayer canopy 
 * and sunlit/shade leaves.
 * (Farquhar et al. 1980/Collatz et al 1991/De Pury and Farquhar)
 * 
 * supplyFunction - Hydraulic supply function
 * 
 * Catm - CO2 air concentration (micromol * mol-1)
 * Patm - Air pressure (kPa)
 * Tair - Air temperature (ºC) - changes through the day and from one day to the other
 * vpa - Air actual vapour pressure (kPa)
 * 
 * SLarea, SHarea - leaf area index of sunlit/shade leaves for each canopy layer
 * u - Wind speed (m/s) for each canopy layer
 * absRadSL, absRadSH - instantaneous absorbed radiation (W·m-2) per unit of sunlit/shade leaf area, for each canopy layer
 * QSL, QSH - Active photon flux density (micromol * s-1 * m-2) per unit of sunlit/shade leaf area, for each canopy layer
 * Vmax298 - maximum Rubisco carboxylation rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2), for each canopy layer
 * Jmax298 - maximum electron transport rate per leaf area at 298 ºK (i.e. 25 ºC) (micromol*s-1*m-2), for each canopy layer
 * 
 * return units: micromol*s-1*m-2
 */

// [[Rcpp::export("photo_sunshadePhotosynthesisFunction")]]
DataFrame sunshadePhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, 
                                  double SLarea, double SHarea,
                                  double u, double absRadSL, double absRadSH,
                                  double QSL, double QSH, 
                                  double Vmax298SL, double Vmax298SH, 
                                  double Jmax298SL, double Jmax298SH, 
                                  double leafWidth = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector Ag(nsteps,0.0), An(nsteps,0.0);
  NumericVector leafCiSL(nsteps,0.0), leafCiSH(nsteps,0.0);
  NumericVector leafTSL(nsteps,0.0), leafTSH(nsteps,0.0);
  NumericVector leafVPDSL(nsteps,0.0), leafVPDSH(nsteps,0.0);
  // Rcout<<"ws "<<u<<" tair "<< Tair<< " SLarea "<< SLarea << " SHarea "<< SHarea<< " absRadSL"<< absRadSL<< " absRadSH "<< absRadSH<< " QSL "<<QSL<<" QSH "<<QSH<<"\n";
  double leafT,Gsw, Gwdiff, Gbw, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    //Sunlit leaves
    Ag[i]=0.0;
    An[i]=0.0;
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature(absRadSL/SLarea, Tair, u, E[i], leafWidth);
    leafTSL[i]= leafT;
    leafVPDSL[i] = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSL[i];
    // Gw = Gw*SLarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSL[i]; //Transform flow from mmol to mol
    Gbw = 0.397*pow(u/(leafWidth*0.0072), 0.5); // mol (boundary layer conductance)
    Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
    Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSL>0.0) {
      NumericVector LP = leafphotosynthesis(QSL, Catm, Gwdiff/1.6, leafT, Vmax298SL, Jmax298SL);//Call photosynthesis with aggregated values
      leafCiSL[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp(Vmax298SL, leafT);
      Ag[i]+=Agj;
      An[i]+=Anj;
    }
    //SHADE leaves
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature(absRadSH/SHarea, Tair, u, E[i], leafWidth);
    leafTSH[i]= leafT;
    leafVPDSH[i] = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
    // leafVPDSH[i] = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSH[i];
    // Gw = Gw*SHarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSH[i]; //Transform flow from mmol to mol
    Gbw = 0.397*pow(u/(leafWidth*0.0072), 0.5); // mol
    Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
    Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSH>0.0) {
      NumericVector LP = leafphotosynthesis(QSH, Catm, Gwdiff/1.6, leafT, Vmax298SH, Jmax298SH); //Call photosynthesis with aggregated values
      leafCiSH[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp(Vmax298SH, leafT);
      Ag[i]+=Agj;
      An[i]+=Anj;
    }
    
  }
  return(DataFrame::create(Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An,
                      Named("LeafCiSL") = leafCiSL,
                      Named("LeafCiSH") = leafCiSH,
                      Named("LeafTempSL") = leafTSL,
                      Named("LeafTempSH") = leafTSH,
                      Named("LeafVPDSL") = leafVPDSL,
                      Named("LeafVPDSH") = leafVPDSH));
}

// [[Rcpp::export("photo_multilayerPhotosynthesisFunction")]]
DataFrame multilayerPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, 
                                           double Catm, double Patm, double Tair, double vpa, 
                                  NumericVector SLarea, NumericVector SHarea,
                                  NumericVector u, NumericVector absRadSL, NumericVector absRadSH,
                                  NumericVector QSL, NumericVector QSH, 
                                  NumericVector Vmax298, NumericVector Jmax298, 
                                  double leafWidth = 1.0, bool verbose = false) {
  int nsteps = E.size();
  int nlayers = SLarea.size();
  NumericVector Ag(nsteps,0.0), An(nsteps,0.0);
  double leafT,leafVPD, Gwdiff, Gsw, Gbw, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    Ag[i]=0.0;
    An[i]=0.0;
    for(int j=0;j<nlayers;j++) {
      //Sunlit leaves
      leafT = leafTemperature(absRadSL[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbw = 0.397*pow(u[j]/(leafWidth*0.0072), 0.5); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
      Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
      if(QSL[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSL[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp(Vmax298[j], leafT);
        //From A per leaf area to A per ground area
        Ag[i]+=Agj*SLarea[j];
        An[i]+=Anj*SLarea[j];
      }
      //SHADE leaves
      leafT = leafTemperature(absRadSH[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbw = 0.397*pow(u[j]/(leafWidth*0.0072), 0.5); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbw); //Diffusive conductance cannot be lower than boundary layer conductance
      Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
      if(QSH[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSH[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp(Vmax298[j], leafT);
        Ag[i]+=Agj*SHarea[j];
        An[i]+=Anj*SHarea[j];
      }
    }
  }
  return(DataFrame::create(Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}



 