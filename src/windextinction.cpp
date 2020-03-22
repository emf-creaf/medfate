#include <Rcpp.h>
#include <math.h> 
#include "fuelstructure.h"
using namespace Rcpp;


/**
 * Albini & Baughman (1979)
 * wind20H: wind speed (m/s) measured at 20 feet (6 m) over the canopy
 * H: canopy height
 * 
 * 1.181102 (m) = 0.36*3.2808399 (ft/m)
 * 0.4265092 (m) = 0.13*3.2808399 (ft/m)
 */
// [[Rcpp::export(".windSpeedAtCanopyHeight")]]
double windSpeedAtCanopyHeight(double wind20H, double canopyHeight) {
  return(1.01857*wind20H/log((20+1.181102*canopyHeight)/(0.4265092*canopyHeight)));
}


/**
 * Unsheltered midflame wind speed
 * Albini & Baughman (1979), for surface fires (no canopy over the flames)
 * Andrews, P.L., 2012. Modeling wind adjustment factor and midflame wind speed for Rothermel’s surface fire spread model. USDA For. Serv. - Gen. Tech. Rep. RMRS-GTR 1–39.
 * wind20H: wind speed (m/s) measured at 20 feet (6 m) over the canopy
 * fuelBedHeight: Fuel bed height above the ground (m)
 * Assumes flameHeightOverFuelBed (i.e. Height of flame over the fuel bed (m)) is equal to fuelbedheihgt
 */
// [[Rcpp::export(".unshelteredMidflameWindSpeed")]]
double unshelteredMidflameWindSpeed(double wind20H, double fuelBedHeight) {
  double num = 1.83*wind20H;
  double den = log((20.0 + 1.181102*fuelBedHeight)/(0.4265092*fuelBedHeight));
  return(num/den);
}

/**
 * Sheltered midflame windspeed (for surface fire)
 * Andrews, P.L., 2012. Modeling wind adjustment factor and midflame wind speed for Rothermel’s surface fire spread model. USDA For. Serv. - Gen. Tech. Rep. RMRS-GTR 1–39.
 * wind20H: wind speed (m/s) measured at 20 feet (6 m) over the canopy
 * topCanopyHeight: canopy height (m)
 */
// [[Rcpp::export(".shelteredMidflameWindSpeed")]]
double shelteredMidflameWindSpeed(double wind20H, double crownFillProportion, double topCanopyHeight) {
  double num = 0.55*wind20H;
  double den = sqrt(crownFillProportion*topCanopyHeight)*log((20.0 + 1.181102*topCanopyHeight)/(0.4265092*topCanopyHeight));
  return(num/den);
}

/**
 * Midflame windspeed (for surface fire) adjustment factor
 * Andrews, P.L., 2012. Modeling wind adjustment factor and midflame wind speed for Rothermel’s surface fire spread model. USDA For. Serv. - Gen. Tech. Rep. RMRS-GTR 1–39.
 * topShrubHeight: canopy height of shrub stratum (m)
 * topCanopyHeight: top canopy height (m)
 * bottomCanopyHeight: bottom canopy height (m)
 * canopyCover: canopy cover (%)
 * 
 */
// [[Rcpp::export("fuel_windAdjustmentFactor")]]
double windAdjustmentFactor(double topShrubHeight, double bottomCanopyHeight, double topCanopyHeight, double canopyCover){
  double crownFillProportion = ((topCanopyHeight-bottomCanopyHeight)/topCanopyHeight)*(canopyCover/300.0);
  if(NumericVector::is_na(topCanopyHeight)) crownFillProportion=0.0;
  double WAF = NA_REAL;
  //unsheltered vs sheltered (5% of crownFillProportion)
  if(crownFillProportion<0.05) WAF = unshelteredMidflameWindSpeed(10.0, topShrubHeight)/10.0;
  else WAF = shelteredMidflameWindSpeed(10.0, crownFillProportion, topCanopyHeight)/10.0;
  return(WAF);
}

/**
 * Albini & Baughman (1979)
 * z: Height above the ground (m)
 * wind20H: wind speed (m/s) measured at 20 feet (6 m) over the canopy
 * H: canopy height
 */
// [[Rcpp::export(".windSpeedAtHeightOverCanopy")]]
double windSpeedAtHeightOverCanopy(double z, double wind20H, double canopyHeight) {
  if(canopyHeight<=0.0) canopyHeight = 0.00000001; //to avoid NaNs
  double Uh = windSpeedAtCanopyHeight(wind20H, canopyHeight);
//  double z0 = ; // friction length in m
  return(Uh*log((z-0.64*canopyHeight)/(0.13*canopyHeight))/log(0.36/0.13));
}

/**
 * From Massman (1987). 
 * A comparative study of some mathematical models of the mean wind structure and aerodynamic drag of plant canopies. 
 * Boundary-Layer Meteorol. 40, 179–197
 */
double windSpeedMassmanExtinction(double z, double wind20H, double LAIc, double canopyHeight) {
  if(canopyHeight<200.0) canopyHeight = 200.0; //Minimum canopy height of 2 m
  double Uh = windSpeedAtCanopyHeight(wind20H, canopyHeight);
  double beta = 4.0*0.2*LAIc/(0.16*pow(1.5,2.0));
  return(Uh*pow(cosh(beta*(z/canopyHeight))/cosh(beta),0.5));
}

/**
 * Following Lalic et al. (2003), to calculate wind inside the canopy (for surface fires under canopy)
 */
double windSpeedLalicExtinction(double z, double wind20H, double LAIc, double canopyHeight, double crownBaseHeight) {
  double Uh = windSpeedAtCanopyHeight(wind20H, canopyHeight);
  if(z<crownBaseHeight) z = crownBaseHeight;
  double beta = 4.0*0.2*LAIc/(0.16*pow(1.5,2.0));
  return(Uh*pow(cosh(beta*(z-crownBaseHeight)/canopyHeight)/cosh(beta*(1.0-(crownBaseHeight/canopyHeight))),3.5));
}


/**
 * Wind extinction profile (m/s)
 */
// [[Rcpp::export(".windExtinctionProfile")]]
NumericVector windExtinctionProfile(NumericVector z, double wind20H, double LAIc, double canopyHeight) {
  int n = z.size();
  NumericVector wp(n);
  for(int i =0;i<n;i++) {
    if(z[i]>canopyHeight) {
      wp[i] = windSpeedAtHeightOverCanopy(z[i], wind20H, canopyHeight);
    } else {
      wp[i] = windSpeedMassmanExtinction(z[i], wind20H, LAIc, canopyHeight);
    }
  }
  return(wp);
}

NumericVector windExtinctionCohort(NumericVector H, NumericVector CR, double wind20H, double LAIc, double canopyHeight) {
  int n = H.size();
  NumericVector wp(n);
  for(int i =0;i<n;i++) {
    wp[i] = windSpeedMassmanExtinction(H[i]*(1.0-(1.0-CR[i])/2.0), wind20H, LAIc, canopyHeight);
  }
  return(wp);
}


double aerodynamicResistance(double canopyHeight, double wind) {
  if(canopyHeight<200.0) canopyHeight = 200.0; //Minimum canopy height of 2 m
  canopyHeight = canopyHeight/100.0;
  double d= 2/3*canopyHeight;
  double zom = 0.123*canopyHeight;
  double zoh = 0.1*zom;
  double zm = 2.0+canopyHeight;
  double zh = zm;
  return(log((zm-d)/zom)*log((zh-d)/zoh)/(pow(0.41,2.0)*wind));
}

