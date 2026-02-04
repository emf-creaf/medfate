//
// CODE COMES FROM METEOLAND, BUT IS REPEATED HERE TO AVOID RCPP STRUCTURES
//
#include "radiation_c.h"

const double GSC_MJminm2 = 0.0820; //solar constant in MJ/m2/min
const double GSC_kWm2 = 1.361; //in kW/m2
const double SIGMA_MJ_day = 4.903*pow(10,-9.0); //Stefan-Boltzmann constant MJ/K^4/m2/day
const double SIGMA_Wm2 = 5.67*pow(10,-8.0); //Stefan-Boltzmann constant J/s/K^4/m2
/**
 * Computes Julian Day from a given date.
 * Julian Day, or number of days since January 1, 4713 BCE at noon UTC.
 *
 * Code from package 'insol' (J. G. Corripio)
 * Danby, J. M. Eqn. 6.16.4 in Fundamentals of Celestial Mechanics, 2nd ed. Richmond, VA: Willmann-Bell, p. 207, 1988.
 * http://scienceworld.wolfram.com/astronomy/JulianDate.html
 */
int julianDay_c(int year, int month, int day) {
  int jd = 367 * year - (int)((7 * (year + (month + 9)/12))/4) + (int)((275 * month)/9) + day + 1721013.5 + 0.5;
  return(jd);
}
/**
 * Returns the solar declination (in radians)
 * Code from package 'insol' (J. G. Corripio)
 * http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
 * Meeus, J. 1999. Astronomical Algorithms. Willmann-Bell, Richmond, Virginia, USA.
 *
 * J - Julian days (integer), number of days since January 1, 4713 BCE at noon UTC.
 */
double solarDeclination_c(int J) {
  // Original calculation
  //double delta = 0.4093 * sin(2.0 * M_PI/365.0 * ((double)J) - 1.405); // solar declination in radians
  // New calculation (Bourges 1985)
  // double D = (2.0*M_PI/365.25)*(((double) J)-79.346);
  // double delta = 0.3723+23.2567*sin(D)-0.758*cos(D)+0.1149*sin(2.0*D)+0.3656*cos(2.0*D)-0.1712*sin(3.0*D)+0.00201*cos(3.0*D);
  // return(delta*(M_PI/180.0));
  //
  double jdc = (((double)J) - 2451545.0)/36525.0;
  double sec = 21.448 - jdc * (46.815 + jdc * (0.00059 - jdc * (0.001813)));
  double e0 = 23.0 + (26.0 + (sec/60.0))/60.0;
  double oblcorr = e0 + 0.00256 * cos((M_PI/180.0)*(125.04 - 1934.136 * jdc));
  double l0 = 280.46646 + jdc * (36000.76983 + jdc * (0.0003032));
  l0 = (double)((((int)l0) - 360*(((int)l0)/360))%360);
  double gmas = (M_PI/180.0)*(357.52911 + jdc * (35999.05029 - 0.0001537 * jdc));
  double seqcent = sin(gmas) * (1.914602 - jdc * (0.004817 + 1.4e-05 * jdc)) + sin(2 * gmas) * (0.019993 - 0.000101 * jdc) + sin(3 * gmas) * 0.000289;
  double suntl = l0 + seqcent;
  double sal = suntl - 0.00569 - 0.00478 * sin((M_PI/180.0)*(125.04 - 1934.136 * jdc));
  double delta = asin(sin((M_PI/180.0)*(oblcorr)) * sin((M_PI/180.0)*(sal)));
  return(delta);
}
 
/**
 * Calculates the earth radius vector (from package 'insol') to correct solar constant
 * Code from package 'insol' (J. G. Corripio)
 */
double solarConstant_c(int J) {
  double jdc = (((double)J) - 2451545.0)/36525.0; //julian century
  double ecc = 0.016708634 - jdc * (4.2037e-05 + 1.267e-07 * jdc);
  double gmasr = (M_PI/180.0)*(357.52911 + jdc * (35999.05029 - 0.0001537 * jdc));
  double seqc = sin(gmasr) * (1.914602 - jdc * (0.004817 + 1.4e-05 *
                    jdc)) + sin(2.0 * gmasr) * (0.019993 - 0.000101 * jdc) +
                    sin(3.0 * gmasr) * 0.000289;
  double star = gmasr + (M_PI/180.0)*seqc;
  double sunrv = (1.000001018 * (1.0 - pow(ecc,2.0)))/(1.0 + ecc * cos(star));
  return(GSC_kWm2*(1.0/sunrv));
}

/**
 * Returns the sunrise and sunset hours in hour angle (radians)
 *
 * latrad - Latitude of actual slope, in radians
 * slorad - Inclination of slope, in radians above horizontal
 * asprad - Azimuth of slope (aspect), in radians from north
 * delta - Solar declination, in radians
 */
void sunRiseSet_c(SunriseSet& ss, double latrad, double slorad, double asprad, double delta){
  double L1 = asin(cos(slorad)*sin(latrad)+sin(slorad)*cos(latrad)*cos(asprad)); //latitude on equivalent slope
  double den = cos(slorad)*cos(latrad)-sin(slorad)*sin(latrad)*cos(asprad);
  double L2;
  if(den<0) {
    L2 = atan((sin(slorad)*sin(asprad))/den)+M_PI;
  } else {
    L2 = atan((sin(slorad)*sin(asprad))/den);
  }
  double T = acos(std::max(std::min(-tan(L1)*tan(delta),1.0),-1.0));
  double T7 = T-L2; //hour angle of sunset on equivalent slope
  double T6 = -T-L2; //hour angle of sunrise on equivalent slope
  double T1 = acos(std::max(std::min(-tan(latrad)*tan(delta),1.0),-1.0));  //hour angle of sunset on horizontal surface
  double T0 = - T1; //hour angle of sunrise on horizontal surface
  double T3 = std::min(T1,T7); //hour angle of sunset on slope
  double T2 = std::max(T0,T6); //hour angle of sunrise on slope
  
  ss.sunrise = T2;
  ss.set = T3;
}

/**
 * Calculates solar elevation angle (in radians)
 *
 * latrad - Latitude (in radians)
 * delta - Solar declination (in radians)
 * hrad - Solar hour angle (in radians)
 */
double solarElevation_c(double latrad, double delta, double hrad) {
  double sinb = sin(latrad)*sin(delta)+cos(latrad)*cos(delta)*cos(hrad);
  return(asin(sinb));
}
 
/**
 *  Calculates day length (in hours)
 *
 *  latrad - Latitude (in radians)
 *  slorad - Inclination of slope, in radians above horizontal
 *  asprad - Azimuth of slope (aspect), in radians from north
 *  delta - Solar declination, in radians
 */
double daylength_c(double latrad, double slorad, double asprad, double delta) {
  SunriseSet v;
  sunRiseSet_c(v, latrad, slorad, asprad, delta); //solar hour at sunrise and sunset, in radians
  double N = 12.0/M_PI*(v.set-v.sunrise);
  return(std::max(N,0.0));
}
/**
 *  Calculates day length (in seconds)
 *
 *  latrad - Latitude (in radians)
 *  slorad - Inclination of slope, in radians above horizontal
 *  asprad - Azimuth of slope (aspect), in radians from north
 *  delta - Solar declination, in radians
 */
double daylengthseconds_c(double latrad, double slorad, double asprad, double delta) {
  return(daylength_c(latrad,slorad, asprad, delta)*3600.0);
}
 
/**
 * Returns instant radiation (in kW/m2). From Granier & Ohmura (1968)
 * Garnier, B.J., Ohmura, A., 1968. A method of calculating the direct shortwave radiation income of slopes. J. Appl. Meteorol. 7: 796-800
 *
 *  solarConstant - Solar constant (in kW/m2)
 *  latrad - latitude of the slope, in radians
 *  A - Azimuth of slope, in radians from north
 *  Zx - Zenith angle of the vector normal to the slope (= slope?), in radians
 *  delta - Solar declination, in radians
 *  H - hour angle measured from solar noon, in radians
 *
 */
double RpotInstant_c(double solarConstant, double latrad, double slorad, double asprad, double delta, double H) {
  double t1 = (sin(latrad)*cos(H))*(-1.0*cos(asprad)*sin(slorad)) -1.0*sin(H)*(sin(asprad)*sin(slorad)) + (cos(latrad)*cos(H))*cos(slorad);
  double t2 = cos(latrad)*(cos(asprad)*sin(slorad)) + sin(latrad)*cos(slorad);
  return(solarConstant*(t1*cos(delta)+ t2*sin(delta)));
}

 
/**
 * One-day potential solar radiation (in MJ/m2). From Granier & Ohmura (1968)
 * Garnier, B.J., Ohmura, A., 1968. A method of calculating the direct shortwave radiation income of slopes. J. Appl. Meteorol. 7: 796-800
 *
 *  solarConstant - Solar constant (in kW/m2)
 *  latrad - latitude of the slope, in radians
 *  asprad - Azimuth of slope, in radians from north
 *  slorad - Zenith angle of the vector normal to the slope (= slope) in radians
 *  J - Julian day (1 to 365)
 *  step - Number of seconds step (for integration of instant radiation)
 */
double RpotDay_c(double solarConstant, double latrad,  double slorad, double asprad, double delta) {
  double step = 600.0; //10 min = 600 s step
  double hradstep = step*(5.0/1200.0)*(M_PI/180.0);
  double Rpot = 0.0;
  SunriseSet srs;
  sunRiseSet_c(srs, latrad, slorad, asprad, delta);
  for(double hrad=srs.sunrise;hrad<srs.set;hrad+=hradstep) {
    Rpot += step*std::max(0.0, RpotInstant_c(solarConstant, latrad, slorad, asprad, delta, hrad));
  }
  return(Rpot/1000.0);
}

/**
 * One-day solar radiation (in MJ/m2). From Thornton & Running (1999)
 *
 *  solarConstant - Solar constant (in kW/m2)
 *  latrad - latitude of the slope, in radians
 *  elevation - Elevation from sea level (in meters)
 *  slorad - Zenith angle of the vector normal to the slope (= slope?) in radians
 *  asprad - Azimuth of slope, in radians from north
 *  delta - Solar declination (in radians)
 *  diffTemp - Difference between maximum and minimum temperature (in degrees Celsius)
 *  vpa - Vapor pressure (in kPa)
 *  precipitation - Precipitation (in mm of water)
 *
 *  Bibliography:
 *  Thornton, P.E., Running, S.W., 1999. An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity, and precipitation. Agric. For. Meteorol. 93, 211–228. doi:10.1016/S0168-1923(98)00126-9
 *  Garnier, B.J., Ohmura, A., 1968. A method of calculating the direct shortwave radiation income of slopes. J. Appl. Meteorol. 7: 796-800
 *  Pearcy, R.W., Ehleringer, J.R., Mooney, H.A., Rundel, P.W., 1991. Plant Physiological Ecology: Field Methods and Instrumentation. Chapman & Hall, New York 455 pp.
 */
double RDay_c(double solarConstant, double latrad, double elevation, double slorad, double asprad, double delta,
              double diffTemp, double diffTempMonth, double vpa, double precipitation) {
  //Number of seconds step (for integration of instant potential radiation and maximum transmittance)
  double step = 600.0; //10 min = 600 s step
  double Rpot = 0.0;
  double RpotFlat = 0.0; //Rpot for a flat surface (needed for diffuse radiation)
  double B = 0.031+0.201*exp(-0.185*diffTempMonth);
  double Tfmax = 1.0-0.9*exp(-B*pow(diffTemp,1.5));
  if(std::isnan(precipitation)) precipitation = 0.0;
  if(precipitation>0.0) Tfmax *=0.75; //Correction for wet days
  double Ttmax = 0.0, RpotInst = 0.0, RpotInstFlat = 0.0, costheta=0.0;
  double pratio = pow(1.0 -2.2569e-5*elevation,5.2553); // from Pearcy et al. 1991
  double hradstep = step*(5.0/1200.0)*(M_PI/180.0);
  SunriseSet srsFlat;
  sunRiseSet_c(srsFlat, latrad, 0.0, 0.0, delta);
  for(double hrad=srsFlat.sunrise;hrad<srsFlat.set;hrad+=hradstep) {
    RpotInstFlat = std::max(0.0,RpotInstant_c(solarConstant, latrad, 0.0, 0.0, delta, hrad));
    RpotFlat += step*RpotInstFlat;
    //Solar zenith angle = 90º - solar elevation angle
    //cos(solar zenith angle) = sin(solar elevation angle)
    costheta = sin(solarElevation_c(latrad,delta, hrad));
    if(costheta>0.0) {
      Ttmax +=step*RpotInstFlat*pow(0.87,pratio*(1.0/costheta));
    }
  }
  SunriseSet srsSlope;
  sunRiseSet_c(srsSlope, latrad, slorad, asprad, delta);
  for(double hrad=srsSlope.sunrise;hrad<srsSlope.set;hrad+=hradstep) {
    RpotInst = std::max(0.0,RpotInstant_c(solarConstant, latrad, slorad, asprad, delta, hrad));
    Rpot += step*RpotInst;
  }
  Ttmax = (Ttmax/RpotFlat) -6.1e-2*vpa; //vpa in kPa
  // Rcout<<Rpot<<" "<<RpotFlat<<" "<<Ttmax<<" "<<Tfmax<<"\n";
  if(RpotFlat==0.0) return(0.0);
  double dirdifCorrected = (Rpot/1000.0)*Ttmax*Tfmax;
  double difCorrected = 0.3*(RpotFlat/1000.0)*(1.0-Ttmax*Tfmax); //Campbell, G.S., & Norman, J.M. 1998. AN INTRODUCTION TO ENVIRONMENTAL BIOPHYSICS.: 2nd edition.
  return(std::max(dirdifCorrected, difCorrected)); //Radiation in MJ/m2
}

 
/**
 * Calculates instantaneous direct and diffuse radiation (in kW) from daily global radiation
 *
 *  solarConstant - Solar constant (in kW/m2)
 *  latrad - Latitude (in radians)
 *  slorad - Zenith angle of the vector normal to the slope (= slope?) in radians
 *  asprad - Azimuth of slope, in radians from north
 *  delta - Solar declination (in radians)
 *  hrad - Solar hour (in radians)
 *  R_p - Potential daily solar radiation (MJ·m-2)
 *  R_s - Global daily solar radiation (MJ·m-2)
 *  clearday - Boolean to indicate that is a clearsky (TRUE) vs overcast (FALSE)
 *
 * Spitters, C.J.T., Toussaint, H.A.J.M. & Goudriaan, J. (1986). Separating the diffuse and direct components of global radiation and its implications for modeling canopy photosynthesis. I. Components of incoming radiation.
 * Agricultural and Forest Meteoroloogy, 38, 231–242.
 */
void directDiffuseInstant_c(DirectDiffuseInstant_RESULT& res,
                            double solarConstant, double latrad, double slorad, double asprad, double delta,
                            double hrad, double R_s,
                            double R_p_flat, double Rpotinst_flat, double R_p_topo, double Rpotinst_topo,
                            bool clearday) {
  // Rcout<< slorad<<" "<<R_p_topo<<"\n";
  //Solar elevation (for corrections)
  double beta = solarElevation_c(latrad, delta, hrad);
  
  //Estimation of SgSo ratio (transmittance) assuming flat surface.
  double SgSoday = R_s/R_p_flat;
  double SdfSgday = NA_REAL;
  if(SgSoday<0.07){
    SdfSgday = 1.0;
  } else if(SgSoday<0.35) {
    SdfSgday = 1.0 - 2.3*pow(SgSoday-0.07,2.0);
  } else if(SgSoday<0.75) {
    SdfSgday = 1.33 -1.46*SgSoday;
  } else {
    SdfSgday = 0.23;
  }
  double SdfSgday2 = SdfSgday;
  //If clear day (e.g. not rainy) modify for circumsolar part of diffuse radiation
  if(clearday) {
    SdfSgday2 = SdfSgday/(1.0+(1.0-pow(SdfSgday,2.0))*pow(cos(M_PI/4.0-beta),2.0)*pow(cos(beta),3.0));
  }
  
  double PARday = R_s*0.5; //Daily PAR radiation (MJ)
  double SdfSdPAR = (1.0+0.3*(1.0-pow(SdfSgday,2.0)))*SdfSgday2;
  double Sdfday = SdfSgday2*R_s; //MJ Diffuse daily radiation
  double Sdrday = R_s - Sdfday; //MJ Direct daily radiation
  double Sdrinst = (Sdrday*1000.0)*(Rpotinst_topo/(R_p_topo*1000.0));//kW Direct light is affected by topography
  double Sdfinst = (Sdfday*1000.0)*(Rpotinst_flat/(R_p_flat*1000.0));//kW Diffuse light not affected by topography
  if(R_p_topo==0.0) {
    Sdrinst = 0.0;
  }
  if(R_p_flat==0.0) {
    Sdfinst = 0.0;
  }
  double Sginst = Sdfinst + Sdrinst;
  double SdfdayPAR = SdfSdPAR*PARday;
  double SginstPAR = Sginst*0.5;
  double SdfinstPAR = std::min((SdfdayPAR*1000.0)*(Rpotinst_flat/(R_p_flat*1000.0)), SginstPAR);//kW
  if(R_p_flat==0.0) {
    SdfinstPAR = 0.0;
  }
  double SdrinstPAR = SginstPAR-SdfinstPAR;
  
  res.SolarElevation = beta;
  res.Rpot= Rpotinst_topo;
  res.Rpot_flat = Rpotinst_flat;
  res.Rg = Sginst;
  res.SWR_direct = Sdrinst;
  res.SWR_diffuse = Sdfinst;
  res.PAR_direct = SdrinstPAR;
  res.PAR_diffuse = SdfinstPAR;
}

/**
 * Calculates daily variation of direct and diffuse radiation (in kW) from daily global radiation
 *
 *  solarConstant - Solar constant (in kW/m2)
 *  latrad - Latitude (in radians)
 *  slorad - Zenith angle of the vector normal to the slope (= slope?) in radians
 *  asprad - Azimuth of slope, in radians from north
 *  delta - Solar declination (in radians)
 *  R_s - Global daily solar radiation (MJ·m-2)
 *  clearday - Boolean to indicate that is a clearsky (TRUE) vs overcast (FALSE)
 *  nsteps - Number of steps to divide the day
 *
 * Spitters, C.J.T., Toussaint, H.A.J.M. & Goudriaan, J. (1986). Separating the diffuse and direct components of global radiation and its implications for modeling canopy photosynthesis. I. Components of incoming radiation.
 * Agricultural and Forest Meteoroloogy, 38, 231–242.
 */
void directDiffuseDay_c(DirectDiffuseDay_RESULT& res, 
                        double solarConstant, double latrad, double slorad, double asprad, double delta,
                        double R_s, bool clearday) {
  int nsteps = res.PAR_diffuse.size();
    // NumericVector Rpot(nsteps),Rpot_flat(nsteps), Rg(nsteps), SWR_direct(nsteps), SWR_diffuse(nsteps), PAR_direct(nsteps), PAR_diffuse(nsteps);
    // NumericVector Hrad(nsteps), beta(nsteps);
  double R_p_flat = RpotDay_c(solarConstant, latrad, 0.0, 0.0, delta);
  double R_p_topo = R_p_flat;
  if(slorad>0.0) {
    R_p_topo = RpotDay_c(solarConstant, latrad, slorad, asprad, delta);
  }
  SunriseSet srs; 
  sunRiseSet_c(srs, latrad, slorad, asprad, delta);
  DirectDiffuseInstant_RESULT ddi;
  for(int i=0;i<nsteps;i++) {
    res.SolarHour[i] = -M_PI + (((double)i)+0.5)*(2.0*M_PI/((double)nsteps));
    double Rpotinst_flat = std::max(0.0,RpotInstant_c(solarConstant, latrad, 0.0, 0.0, delta, res.SolarHour[i]));//kW
    double Rpotinst_topo = Rpotinst_flat;
    if(slorad>0.0) {
      Rpotinst_topo = 0.0;
      if(res.SolarHour[i] >= srs.sunrise && res.SolarHour[i]< srs.set) Rpotinst_topo = std::max(0.0,RpotInstant_c(solarConstant, latrad, slorad, asprad, delta, res.SolarHour[i]));//kW
    }
    directDiffuseInstant_c(ddi, solarConstant, latrad, slorad, asprad, delta,
                           res.SolarHour[i], R_s, R_p_flat, Rpotinst_flat, R_p_topo, Rpotinst_topo, clearday);
    res.SolarElevation[i] = ddi.SolarElevation;
    res.Rpot[i] = Rpotinst_topo;
    res.Rpot_flat[i] = Rpotinst_flat;
    res.Rg[i] = ddi.Rg;
    res.SWR_direct[i] = ddi.SWR_direct;
    res.SWR_diffuse[i] = ddi.SWR_diffuse;
    res.PAR_direct[i] = ddi.PAR_direct;
    res.PAR_diffuse[i] = ddi.PAR_diffuse;
  }
}
/**
 *  Estimates instantaneous thermal (long-wave) radiation from the sky (i.e. counterradiation)
 *  Units: W/m2
 *
 *  Campbell, G.S., & Norman, J.M. 1998. AN INTRODUCTION TO ENVIRONMENTAL BIOPHYSICS.: 2nd edition.
 *
 *  Tair - Air temperature (ºC)
 *  vpa - water vapour pressure (kPa)
 *  c - fraction of sky covered by clouds
 */
double skyLongwaveRadiation_c(double Tair, double vpa, double c) {
  double Tkelv = Tair+273.0;
  double bb = SIGMA_Wm2*pow(Tkelv,4.0); //eq. 10.7 (Stefan-Boltzmann law for blackbodies)
  double eac = 1.72 * pow(vpa/Tkelv,1.0/7.0); // eq. 10.10
  double ea = (1.0 - 0.84*c)*eac + 0.84*c; //eq.10.12
  return(ea*bb); // eq. 10.9 Energy emitted by gray bodies (no dependence of emissivity on wavelength)
}
Rcpp::DataFrame copyDirectDiffuseDayResult_c(const DirectDiffuseDay_RESULT& ddd) {
  Rcpp::DataFrame res = Rcpp::DataFrame::create(Rcpp::Named("SolarHour") = Rcpp::wrap(ddd.SolarHour),
                                                Rcpp::Named("SolarElevation") = Rcpp::wrap(ddd.SolarElevation),
                                                Rcpp::Named("Rpot") = Rcpp::wrap(ddd.Rpot),
                                                Rcpp::Named("Rpot_flat") = Rcpp::wrap(ddd.Rpot_flat),
                                                Rcpp::Named("Rg") = Rcpp::wrap(ddd.Rg),
                                                Rcpp::Named("SWR_direct") = Rcpp::wrap(ddd.SWR_direct),
                                                Rcpp::Named("SWR_diffuse") = Rcpp::wrap(ddd.SWR_diffuse),
                                                Rcpp::Named("PAR_direct") = Rcpp::wrap(ddd.PAR_direct),
                                                Rcpp::Named("PAR_diffuse") = Rcpp::wrap(ddd.PAR_diffuse));
  return(res);
}