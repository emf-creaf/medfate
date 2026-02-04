#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "biophysicsutils_c.h"
#include "photosynthesis_c.h"
#include <meteoland.h>
using namespace Rcpp;



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
//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_photosynthesis")]]
NumericVector leafphotosynthesis(double Q, double Catm, double Gc, double Tleaf, double Vmax298, double Jmax298, bool verbose=false) {
  Photo photo;
  leafphotosynthesis_inner_c(photo, Q, Catm, Gc, Tleaf, Vmax298, Jmax298);
  NumericVector res = NumericVector::create(photo.Ci, photo.A);
  res.attr("names") = CharacterVector::create("Ci", "A");
  return(res);
}



// From Baldocchi D (1994). An analytical solution for the coupled leaf photosynthesis and stomatal conductance models. Tree Physiology 14: 1069-1079 
//' @rdname photo
//' @param Gsw_AC_slope Slope of the An/C vs Gsw relationship 
//' @param Gsw_AC_intercept Intercept of the An/C vs Gsw relationship 
//' @keywords internal
// [[Rcpp::export("photo_photosynthesisBaldocchi")]]
NumericVector photosynthesisBaldocchi(double Q, 
                                      double Catm, 
                                      double Tleaf, 
                                      double u,
                                      double Vmax298, 
                                      double Jmax298, 
                                      double leafWidth,
                                      double Gsw_AC_slope,
                                      double Gsw_AC_intercept) {
  BaldocchiPhoto photoOut;
  photosynthesisBaldocchi_inner_c(photoOut, Q, Catm, Tleaf, u,Vmax298,Jmax298, leafWidth,Gsw_AC_slope,Gsw_AC_intercept);
  NumericVector res = {photoOut.Gsw, 
                       photoOut.Cs,
                       photoOut.Ci,
                       photoOut.An,
                       photoOut.Ag};
  res.attr("names") = CharacterVector::create("Gsw", "Cs" ,"Ci", "An", "Ag");
  return(res);
}

//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_leafPhotosynthesisFunction")]]
DataFrame leafPhotosynthesisFunction(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                             double absRad, double Q, double Vmax298, double Jmax298, 
                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  double Gwdiff, Gbound;
  for(int i=0;i<nsteps;i++){
    leafTemp[i] = leafTemperature_c(absRad/refLeafArea, Tair, u, E[i], leafWidth);
    leafVPD[i] = std::max(0.0,leafVapourPressure_c(leafTemp[i], psiLeaf[i]) - vpa);
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPD[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary_c(u, leafWidth); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive resistance cannot be smaller than the boundary layer resistance
    Gsw[i]  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbound))); //Determine stomatal conductance after accounting for leaf boundary conductance
    NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp[i]), Vmax298/refLeafArea, Jmax298/refLeafArea);
    Ci[i] = LP[0];
    Ag[i] = LP[1];
    An[i] = Ag[i] - 0.015*VmaxTemp_c(Vmax298/refLeafArea, leafTemp[i]);
  }
  return(DataFrame::create(Named("LeafTemperature") = leafTemp,
                      Named("LeafVPD") = leafVPD,
                      Named("Gsw") = Gsw,
                      Named("Ci") = Ci,
                      Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}


NumericVector leafPhotosynthesisOneFunction2(double E, double psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                             double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                             double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  double leafTemp, leafVPD, Gsw, Ci;
  double Ag, An;
  double Gwdiff, Gbound;
  leafTemp = leafTemperature2_c(SWRabs/refLeafArea, LWRnet/refLeafArea, Tair, u, E, leafWidth);
  leafVPD = std::max(0.0,leafVapourPressure_c(leafTemp, psiLeaf) - vpa);
  // Separates diffusive conductance into stomatal and boundary layer conductance
  Gwdiff = Patm*(E/1000.0)/leafVPD; //Transform flow from mmol to mol
  Gbound = gLeafBoundary_c(u, leafWidth); // mol boundary layer conductance
  Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be smaller than the boundary layer resistances
  Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbound))); //Determine stomatal conductance after accounting for leaf boundary conductance
  NumericVector LP = leafphotosynthesis(Q/refLeafArea, Catm, Gwdiff/1.6, std::max(0.0,leafTemp), Vmax298/refLeafArea, Jmax298/refLeafArea);
  Ci = LP[0];
  Ag = LP[1];
  An = Ag - 0.015*VmaxTemp_c(Vmax298/refLeafArea, leafTemp);
  return(NumericVector::create(Named("LeafTemperature") = leafTemp,
                               Named("LeafVPD") = leafVPD,
                               Named("Gsw") = Gsw,
                               Named("Ci") = Ci,
                               Named("GrossPhotosynthesis") = Ag,
                               Named("NetPhotosynthesis") = An));
}


//' @rdname photo
//' @keywords internal
// [[Rcpp::export("photo_leafPhotosynthesisFunction2")]]
DataFrame leafPhotosynthesisFunction2(NumericVector E, NumericVector psiLeaf, double Catm, double Patm, double Tair, double vpa, double u, 
                                     double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                                     double leafWidth = 1.0, double refLeafArea = 1.0, bool verbose = false) {
  int nsteps = E.size();
  NumericVector leafTemp(nsteps);
  NumericVector leafVPD(nsteps);
  NumericVector Gsw(nsteps), Ci(nsteps);
  NumericVector Ag(nsteps), An(nsteps);
  for(int i=0;i<nsteps;i++){
    NumericVector lpf = leafPhotosynthesisOneFunction2(E[i], psiLeaf[i], Catm, Patm, Tair, vpa, u,
                                                       SWRabs, LWRnet, Q, Vmax298, Jmax298,
                                                       leafWidth, refLeafArea, verbose);
    leafTemp[i] = lpf["LeafTemperature"];
    leafVPD[i] = lpf["LeafVPD"];
    Gsw[i] = lpf["Gsw"];
    Ci[i] = lpf["Ci"];
    Ag[i] = lpf["GrossPhotosynthesis"];
    An[i] = lpf["NetPhotosynthesis"];
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
//' @rdname photo
//' @keywords internal
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
  double leafT, Gwdiff, Gbound, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    //Sunlit leaves
    Ag[i]=0.0;
    An[i]=0.0;
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature_c(absRadSL/SLarea, Tair, u, E[i], leafWidth);
    leafTSL[i]= leafT;
    leafVPDSL[i] = std::max(0.0,leafVapourPressure_c(leafT, psiLeaf[i]) - vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSL[i];
    // Gw = Gw*SLarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSL[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary_c(u, leafWidth); // mol (boundary layer conductance)
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSL>0.0) {
      NumericVector LP = leafphotosynthesis(QSL, Catm, Gwdiff/1.6, leafT, Vmax298SL, Jmax298SL);//Call photosynthesis with aggregated values
      leafCiSL[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp_c(Vmax298SL, leafT);
      Ag[i]+=Agj;
      An[i]+=Anj;
    }
    //SHADE leaves
    //From rad per ground area to rad per leaf area
    leafT = leafTemperature_c(absRadSH/SHarea, Tair, u, E[i], leafWidth);
    leafTSH[i]= leafT;
    leafVPDSH[i] = std::max(0.0,leafVapourPressure_c(leafT, psiLeaf[i]) - vpa);
    // leafVPDSH[i] = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
    // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
    // Gw = Patm*(E[i]/1000.0)/leafVPDSH[i];
    // Gw = Gw*SHarea; //From Gw per leaf area to Gw per ground area
    // Separates diffusive conductance into stomatal and boundary layer conductance
    Gwdiff = Patm*(E[i]/1000.0)/leafVPDSH[i]; //Transform flow from mmol to mol
    Gbound = gLeafBoundary_c(u, leafWidth); // mol boundary layer conductance
    Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
    // Gsw  = std::abs(1.0/((1.0/Gwdiff) - (1.0/Gbw))); //Determine stomatal conductance after accounting for leaf boundary conductance
    Gwdiff = Gwdiff*SLarea; //From Gwdiff per leaf area to Gwdiff per ground area
    if(QSH>0.0) {
      NumericVector LP = leafphotosynthesis(QSH, Catm, Gwdiff/1.6, leafT, Vmax298SH, Jmax298SH); //Call photosynthesis with aggregated values
      leafCiSH[i] = LP[0];
      Agj = LP[1];
      Anj = Agj - 0.015*VmaxTemp_c(Vmax298SH, leafT);
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

//' @rdname photo
//' @keywords internal
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
  double leafT,leafVPD, Gwdiff, Gbound, Agj, Anj;
  for(int i=0;i<nsteps;i++){
    Ag[i]=0.0;
    An[i]=0.0;
    for(int j=0;j<nlayers;j++) {
      //Sunlit leaves
      leafT = leafTemperature_c(absRadSL[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure_c(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbound = gLeafBoundary_c(u[j], leafWidth); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
      if(QSL[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSL[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp_c(Vmax298[j], leafT);
        //From A per leaf area to A per ground area
        Ag[i]+=Agj*SLarea[j];
        An[i]+=Anj*SLarea[j];
      }
      //SHADE leaves
      leafT = leafTemperature_c(absRadSH[j], Tair, u[j], E[i], leafWidth);
      leafVPD = std::max(0.0,leafVapourPressure_c(leafT, psiLeaf[i]) - vpa);
      // leafVPD = std::max(0.0,meteoland::utils_saturationVP(leafT)-vpa);
      // Assumes infinite boundary layer conductance (well-coupling) so that stomatal conductance equals diffusive conductance
      // Gw = Patm*(E[i]/1000.0)/leafVPD;
      // Separates diffusive conductance into stomatal and boundary layer conductance
      Gwdiff = Patm*(E[i]/1000.0)/leafVPD; //Transform flow from mmol to mol
      Gbound = gLeafBoundary_c(u[j], leafWidth); // mol boundary layer conductance
      Gwdiff = std::min(Gwdiff, Gbound); //Diffusive conductance cannot be lower than boundary layer conductance
      if(QSH[j]>0.0) {
        NumericVector LP = leafphotosynthesis(QSH[j], Catm, Gwdiff/1.6, leafT, Vmax298[j], Jmax298[j]);
        Agj = LP[1];
        Anj = Agj - 0.015*VmaxTemp_c(Vmax298[j], leafT);
        Ag[i]+=Agj*SHarea[j];
        An[i]+=Anj*SHarea[j];
      }
    }
  }
  return(DataFrame::create(Named("GrossPhotosynthesis") = Ag,
                      Named("NetPhotosynthesis") = An));
}



 
