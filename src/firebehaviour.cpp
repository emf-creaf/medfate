// [[Rcpp::depends(meteoland)]]
#include <Rcpp.h>
#include <numeric>
#include <math.h>
#include "spwb.h"
#include "forestutils.h"
#include "fuelstructure.h"
#include "fuelmoisture.h"
#include "windextinction.h"
#include <meteoland.h>
using namespace Rcpp;

const double V = 1116.0; //Latent heat of vaporisation of water in Btu/lb
const double etaK = 0.42; //Mineral damping coefficient
const double Mx_dead_ca = 0.25; //Moisture of extinction for dead canopy in proportion
const double Mx_dead_sh = 0.25; //Moisture of extinction for dead shrub in proportion
const double Mx_dead_he = 0.25; //Moisture of extinction for dead herbs in proportion
const double Mx_wo = 0.25; //Moisture of extinction for woody in proportion
const double Mx_li = 0.25; //Moisture of extinction for litter in proportion
const double Mx_live_ca = 1.80; //Moisture of extinction for live canopy in proportion
const double Mx_live_sh = 1.80; //Moisture of extinction for live shrub in proportion
const double Mx_live_he = 1.20; //Moisture of extinction for live herb in proportion
const double BMU = 352.0; //Benchmark midflame windspeed (ft/min)
const double VS = 900.0; //Vertical stack velocity (ft/min)
const double B = 1.2; //Exponential response of wind coefficient to windspeed (Sandberg et al. 2007)


/**
 * Vector addition from polar coordinates (length, angles in radians)
 * Angles are measured from the y-axis (north)
 */
NumericVector vectorAddition(NumericVector v1, NumericVector v2) {
  //add coordinates
  double x = v1[0]*sin(v1[1])+v2[0]*sin(v2[1]);
  double y = v1[0]*cos(v1[1])+v2[0]*cos(v2[1]);
  //  Rcout << x << " "<< y<<"\n";
  return(NumericVector::create(sqrt(pow(x,2.0)+pow(y,2.0)), atan2(x, y)));
}

/**
 * From van Wagner (1977) model
 * CBH - crown base height (in m)
 * M - Crown foliar moisture (in percent)
 */
// [[Rcpp::export(".criticalFirelineIntensity")]]
double criticalFirelineIntensity(double CBH, double M) {
  return(pow(0.010*(CBH)*(460.0+25.9*M),1.5));
}

/**
 * FCCS
 * 
 *  FCCSpropsSI - Dataframe with fuel properties
 *  MliveSI - Moisture of live materials (in percent of dry weight) for canopy, shrub, and herb strata
 *  MdeadSI - Moisture of dead materials (in percent of dry weight) for canopy, shrub, herb, woody and litter strata
 *  windSpeedSI - Wind speed (m/s) at 20 ft (6 m) over vegetation (default 11 m/s = 40 km/h)
 *  slope - Slope (in degrees)
 *  
 *  Default moisture, slope and windspeed values are benchmark conditions used 
 *  to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire
 *  
 *  Strata indices:  0 - Canopy, 1 - Shrub, 2- Herb, 3 - Woody, 4 - Litter
 */
// [[Rcpp::export("fire_FCCS")]]
List FCCSbehaviour(DataFrame FCCSpropsSI,
          NumericVector MliveSI = NumericVector::create(90, 90, 60), 
          NumericVector MdeadSI = NumericVector::create(6, 6, 6, 6, 6), 
          double slope = 0.0, double windSpeedSI = 11.0) {
  //Extract vectors
  NumericVector wSI = FCCSpropsSI["w"];
  NumericVector cover = FCCSpropsSI["cover"];
  NumericVector hbcSI = FCCSpropsSI["hbc"];
  NumericVector htcSI = FCCSpropsSI["htc"];
  NumericVector deltaSI = FCCSpropsSI["delta"];
  NumericVector rhopSI = FCCSpropsSI["rhop"];
  NumericVector PVSI = FCCSpropsSI["PV"];
  NumericVector beta = FCCSpropsSI["beta"]; //unitless
  NumericVector betarel = FCCSpropsSI["betarel"]; //unitless
  NumericVector etabetarel = FCCSpropsSI["etabetarel"]; //unitless
  NumericVector sigmaSI = FCCSpropsSI["sigma"];
  NumericVector pDead = FCCSpropsSI["pDead"];
  NumericVector FAI = FCCSpropsSI["FAI"]; //unitless
  NumericVector hdefSI = FCCSpropsSI["h"];
  NumericVector etaF = FCCSpropsSI["etaF"]; //unitless
  NumericVector RVSI = FCCSpropsSI["RV"];
  //Rescale variables to British units
  NumericVector Mlive = MliveSI/100.0; //from percent to proportions
  NumericVector Mdead = MdeadSI/100.0; //from percent to proportions
  NumericVector w = wSI * 0.204918; //from kg/m2 to lb/ft2
  NumericVector delta = deltaSI * 3.2808399; //from m to ft
  NumericVector rhop = rhopSI * 0.06242796; //from kg/m3 to lb/ft3
  NumericVector PV = PVSI * 3.2808399; //from m3/m2 to ft3/ft2
  NumericVector sigma = sigmaSI * 0.3048; //from m2/m3 to ft2/ft3
  NumericVector hdef = hdefSI * 0.429922614; //from kJ/kg to Btu/lb
  NumericVector RV = RVSI * 3.2808399; //from m3/m2 to ft3/ft2

  //Heat content (after correcting for moisture in live fuels)
  NumericVector h(5);
  h[0] = hdef[0] - Mlive[0]*V;
  h[1] = hdef[1] - Mlive[1]*V;
  h[2] = hdef[2] - Mlive[2]*V;
  h[3] = hdef[3];
  h[4] = hdef[4];
  // Rcout<<"H: "<<h[0]<<" "<<h[1]<<" "<<h[2]<<" "<<h[3]<< " "<<h[4]<<"\n";

  
  //Rothermel's A
  double A_ca = 133*pow(sigma[0], -0.7913);
  double A_sh = 133*pow(sigma[1], -0.7913);
  double A_he = 133*pow(sigma[2], -0.7913);
  double A_wo = 1.0;
  double A_li = 1.0;
  // Rcout<<"A: "<<A_sh<<" "<<A_he<<" "<<A_wo<<" "<<A_li<<"\n";
  
  //Maximum reaction velocity (in min-1)
  double Gmax_ca = 15.0;
  double Gmax_sh = 9.495 *(sigma[1]/sigma[3]);
  double Gmax_he = 9.495 *(sigma[2]/sigma[3]);
  double Gmax_wo = 9.495;
  double Gmax_li = 15.0;
  // Rcout<<"Gmax: "<<Gmax_sh<<" "<<Gmax_he<<" "<<Gmax_wo<<" "<<Gmax_li<<"\n";
  
  //Moisture damping coefficient
  NumericVector etaMlive(3);
  NumericVector etaMdead(3);
  NumericVector etaM(5);
  etaMlive[0] = (1.0 - 2.59*(Mlive[0]/Mx_live_ca))+(5.11*pow(Mlive[0]/Mx_live_ca,2.0))-(3.52*pow(Mlive[0]/Mx_live_ca,3.0));
  etaMlive[1] = (1.0 - 2.59*(Mlive[1]/Mx_live_sh))+(5.11*pow(Mlive[1]/Mx_live_sh,2.0))-(3.52*pow(Mlive[1]/Mx_live_sh,3.0));
  etaMlive[2] = (1.0 - 2.59*(Mlive[2]/Mx_live_he))+(5.11*pow(Mlive[2]/Mx_live_he,2.0))-(3.52*pow(Mlive[2]/Mx_live_he,3.0));
  etaMdead[0] = (1.0 - 2.59*(Mdead[0]/Mx_dead_ca))+(5.11*pow(Mdead[0]/Mx_dead_ca,2.0))-(3.52*pow(Mdead[0]/Mx_dead_ca,3.0));
  etaMdead[1] = (1.0 - 2.59*(Mdead[1]/Mx_dead_sh))+(5.11*pow(Mdead[1]/Mx_dead_sh,2.0))-(3.52*pow(Mdead[1]/Mx_dead_sh,3.0));
  etaMdead[2] = (1.0 - 2.59*(Mdead[2]/Mx_dead_he))+(5.11*pow(Mdead[2]/Mx_dead_he,2.0))-(3.52*pow(Mdead[2]/Mx_dead_he,3.0));
  etaM[0] = std::max(0.0,etaMlive[0]*(1.0 - pDead[0]) + etaMdead[0]*pDead[0]);
  etaM[1] = std::max(0.0,etaMlive[1]*(1.0 - pDead[1]) + etaMdead[1]*pDead[1]);
  etaM[2] = std::max(0.0,etaMlive[2]*(1.0 - pDead[2]) + etaMdead[2]*pDead[2]);
  etaM[3] = std::max(0.0,(1.0 - 2.59*(Mdead[3]/Mx_wo))+(5.11*pow(Mdead[3]/Mx_wo,2.0))-(3.52*pow(Mdead[3]/Mx_wo,3.0)));
  etaM[4] = std::max(0.0,(1.0 - 2.59*(Mdead[4]/Mx_li))+(5.11*pow(Mdead[4]/Mx_li,2.0))-(3.52*pow(Mdead[4]/Mx_li,3.0)));
  // Rcout<<"etaM: "<<etaM[0]<<" "<<etaM[1]<<" "<<etaM[2]<<" "<<etaM[3]<<" "<<etaM[4]<<"\n";
  
  //Reaction intensity
  double I_r_sh = pow(etabetarel[1], A_sh)*Gmax_sh*w[1]*h[1]*etaM[1]*etaK*etaF[1];
  double I_r_he = pow(etabetarel[2], A_he)*Gmax_he*w[2]*h[2]*etaM[2]*etaK*etaF[2];
  double I_r_wo = pow(etabetarel[3], A_wo)*Gmax_wo*w[3]*h[3]*etaM[3]*etaK*etaF[3];
  double I_r_li = pow(etabetarel[4], A_li)*Gmax_li*w[4]*h[4]*etaM[4]*etaK*etaF[4];
  double I_r_surf = I_r_sh+ I_r_he+I_r_wo+I_r_li;
  double I_r_litter = I_r_li;
  
  //Propagating flux ratio (unitless)
  double delta_surf_heatsink = (RV[1]*delta[1]+RV[2]*delta[2]+RV[3]*delta[3]+RV[4]*delta[4])/(RV[1]+RV[2]+RV[3]+RV[4]);
  double xi_surf = 0.03 + 2.5*std::min(0.06, (RV[1]+RV[2]+RV[3]+RV[4])/delta_surf_heatsink);
  double xi_litter = 0.03 + 2.5*std::min(0.06, RV[4]/delta[4]);
  
  //Heat of preignition (Btu/lb)
  double Qig_live_ca = 250.0 + (V*(Mlive[0]));
  double Qig_live_sh = 250.0 + (V*(Mlive[1]));
  double Qig_live_he = 250.0 + (V*(Mlive[2]));
  double Qig_dead_ca = 250.0;
  double Qig_dead_sh = 250.0;
  double Qig_dead_he = 250.0;
  double Qig_ca = Qig_live_ca*(1.0-pDead[0])+Qig_dead_ca*pDead[0];
  double Qig_sh = Qig_live_sh*(1.0-pDead[1])+Qig_dead_sh*pDead[1];
  double Qig_he = Qig_live_he*(1.0-pDead[2])+Qig_dead_he*pDead[2];
  double Qig_wo = 250.0;
  double Qig_li = 250.0;
  
  //Heat sink (Btu/ft3)
  double q_sh = 0.0, q_he=0.0, q_wo=0.0, q_li=0.0;
  if(delta[1]>0.0) q_sh = etabetarel[1]*RV[1]*rhop[1]*Qig_sh/std::min(1.0,delta[1]);
  if(delta[2]>0.0) q_he = etabetarel[2]*RV[2]*rhop[2]*Qig_he/std::min(1.0,delta[2]);
  if(delta[3]>0.0) q_wo = etabetarel[3]*RV[3]*rhop[3]*Qig_wo/std::min(1.0,delta[3]);
  if(delta[4]>0.0) q_li = etabetarel[4]*RV[4]*rhop[4]*Qig_li/std::min(1.0,delta[4]);
  double q_surf = q_sh + q_he + q_wo + q_li;
  double q_litter = q_li;
  
  //Wind extinction
  //Andrews, P.L., 2012. Modeling wind adjustment factor and midflame wind speed for Rothermel’s surface fire spread model. USDA For. Serv. - Gen. Tech. Rep. RMRS-GTR 1–39.
  double crownFillProportion = ((htcSI[0]-hbcSI[0])/htcSI[0])*(cover[0]/300.0);
  // Rcout<<crownFillProportion<<"\n";
  if(NumericVector::is_na(htcSI[0])) crownFillProportion=0.0;
  double midflameWindSpeedSI = NA_REAL;
  //unsheltered vs sheltered (5% of crownFillProportion)
  //from m/s to mph
  if(crownFillProportion<0.05) midflameWindSpeedSI = unshelteredMidflameWindSpeed(windSpeedSI, delta_surf_heatsink/3.2808399);
  else midflameWindSpeedSI = shelteredMidflameWindSpeed(windSpeedSI, crownFillProportion, htcSI[0]);
  double crownfireWindSpeedSI = windSpeedAtCanopyHeight(windSpeedSI,htcSI[0]);
  double midflameWindSpeed = midflameWindSpeedSI *2.23693629; 
  double crownfireWindSpeed = crownfireWindSpeedSI *2.23693629; 
  
  //Wind coefficient
  double E = 0.55 - 0.2*((FAI[1]+FAI[2])/(FAI[1]+FAI[2]+FAI[3]));
  double phi_wind_surf = 8.8*pow(betarel[1], -1.0*E)*pow(88.0*midflameWindSpeed/BMU,B);
  double phi_wind_litter = 8.8*pow(betarel[4], -1.0*E)*pow(88.0*midflameWindSpeed/BMU,B);
  double phi_wind = ((1.0 - (I_r_litter/I_r_surf))*phi_wind_surf) + ((I_r_litter/I_r_surf)*phi_wind_litter);
  // Rcout<< "Wind: "<< E<<" "<<phi_wind_surf<<" "<<phi_wind_litter<<"\n";
    
  //Slope coefficient
  double phi_slope_surf = 5.275*pow(slope/100.0,2.0)*pow(beta[1]+beta[2]+beta[3],-0.3);
  double phi_slope_litter = 5.275*pow(slope/100.0,2.0)*pow(beta[4],-0.3);
  double phi_slope = ((1.0 - (I_r_litter/I_r_surf))*phi_slope_surf) + ((I_r_litter/I_r_surf)*phi_slope_litter);
  
  //Rate of spread of a surface fire  (ft/min)
  double ros_surf = I_r_surf * xi_surf * (1.0 + phi_wind + phi_slope)/q_surf;
  //Rate of spread of a litter fire  (ft/min)
  double ros_litter = I_r_litter * xi_litter * (1.0 + phi_wind + phi_slope)/q_litter;
  //Maximum rate of spread calculated from wind and slope (ft/min)
  double windslopecap = 88.0*midflameWindSpeed * (1+ phi_slope);
  //Final rate of spread (ft/min)
  double ros = std::min(windslopecap, std::max(ros_surf, ros_litter));
  // Rcout<<"ros: "<<ros<<"\n";
  
  //Reaction thickness (ft)
  double RT_ca = std::min(0.0028, 2.0/sigma[0]);
  double RT_sh = std::min(0.0028, 2.0/sigma[1]);
  double RT_he = std::min(0.0028, 2.0/sigma[2]);
  double RT_wo = std::min(0.0028, 2.0/sigma[3]);
  double RT_li = std::min(0.0028, 2.0/sigma[4]);
  
  //Residence time (min)
  double t_r = 192.0*((I_r_sh*RT_sh) + (I_r_he*RT_he)+(I_r_wo*RT_wo)+(I_r_li*RT_li))/I_r_surf;
  // Rcout<<"T_r: "<<t_r<<"\n";
  
  //Fireline intensity (Btu/ft/min)
  double I_b = I_r_surf*ros*t_r;
  // Rcout<<"I_b: "<<I_b<<"\n";
  
  //Flame length (ft)
  double FL = 0.45*pow(I_b/60.0, 0.46);
  // Rcout<<"FL: "<<FL<<"\n";
  
  //Surface fire potentials
  double RP = 0.0, SP = 0.0, FP = 0.0;
  if(!NumericVector::is_na(I_r_surf)) RP = std::min(9.0,0.08*pow(I_r_surf, 0.5));
  if(!NumericVector::is_na(ros)) SP = std::min(9.0, 2.5*pow(ros, 0.5));
  if(!NumericVector::is_na(FL)) FP = std::min(9.0, 2.5*pow(FL, 0.5));
  double SFP = std::min(9.0, std::max(SP,FP));
  
  //Canopy gap
  double GAP_ca = hbcSI[0]-htcSI[1];
  
  //Van Wagner's critical fireline intensity (in Btu/ft/s)
  double I_c = 0.288894658*criticalFirelineIntensity(GAP_ca, 100.0*Mlive[0]); //express GAP in m and M in percent before calling Van Wagner
  // Rcout<<"GAP (m) "<<GAP_ca<<" I_b: "<< (I_b/60.0)<<" I_c: "<< I_c<<"\n";

  //Canopy windspeed adjustment factor
  double U = 88.0*crownfireWindSpeed; // conversion to ft/min
  double WAF = (U / sqrt(pow(U,2.0) + pow(VS, 2.0)))/(BMU / sqrt(pow(BMU,2.0) + pow(VS, 2.0)));
  // Rcout<<"WAF "<<WAF<<"\n";
  
  //Efficiency of crown-to-crown transfer
  double TCq = pow(std::max(0.0, cover[0]*WAF-40.0),0.3)/pow(100.0*WAF-40.0, 0.3);
  
  //Threshold for FAIc
  double Ac = 2.6296; //if sigma_c <= 2000 ft2/ft3
  if(sigma[0] > 2000) Ac = 3.2868; //if sigma_c > 2000 ft2/ft3
  double TFAIc = Ac*exp(-0.0019*U);
  // Rcout << " FAIc "<< FAI[0] << " TFAI/3pi "<< (TFAIc/(3.0*PI))<<"\n";

  //Crown fire spread rate
  double I_r_ca  = pow(etabetarel[0], A_ca)*Gmax_ca*w[0]*h[0]*etaM[0]*etaK*etaF[0];
  double I_r_crown = I_r_surf+I_r_ca;
  double xi_crown = 1.0 - exp(-1.0*(FAI[0]/(4.0*delta[0])));
  double q_ca = (0.5*FAI[0]*RT_ca*rhop[0]*Qig_ca)/((cover[0]/100.0)*delta[0]);
  double q_crown = q_surf + q_ca;
  double ros_crown = I_r_crown * xi_crown*WAF/q_ca;
  
  //Crown fire residence time (min)
  double t_r_crown = 192.0*RT_ca;
  
  //Crown fire fireline intensity (Btu/ft/min)
  double I_b_crown = I_r_crown*ros_crown*t_r_crown;

  //Crown flame length (ft)
  double FL_crown = 0.45*pow(I_b_crown/60.0, 0.46);
  
  //Crown fire porentials
  double IC = 0.0;
  double RC = 0.0;
  double TC = 0.0;
  double ic_ratio = NA_REAL;
  if(FAI[0]>0.0) {
    if((!NumericVector::is_na(I_b)) & (!NumericVector::is_na(I_c))) {
      ic_ratio = (I_b/60.0)/I_c;
      IC = std::min(9.0, 4.0*pow(ic_ratio,0.2));
    }
  } 
  if(FAI[0]> (TFAIc/(3.0*PI))) {
    if(!NumericVector::is_na(TCq)) TC = std::min(9.0, 10.0*TCq);
  }
  if(!NumericVector::is_na(ros_crown)) RC = std::min(9.0, 1.0*pow(ros_crown, 0.5));

  // double CFP = 0.4286*(IC+(TC/3.0)+RC);
  double CFP = pow(IC*RC, 0.5);
  //   double AC = pow(IC*TC*RC,1.0/3.0);
//   double CFP = std::max(IC, AC);
  
  //Rescale output to metric units
  ros_surf = ros_surf * 0.3048; //ft/min to m/min
  ros_litter = ros_litter * 0.3048; //ft/min to m/min
  windslopecap = windslopecap* 0.3048; //ft/min to m/min
  ros = ros * 0.3048; //ft/min to m/min
  ros_crown = ros_crown* 0.3048; //ft/min to m/min
  FL = FL * 0.3048; //ft to m
  FL_crown = FL_crown * 0.3048; //ft to m
  //From Kennard, Relationship Between Flame Length and Fireline Intensity — Forest Encyclopedia 
  //Add half of mean canopy top height to FL of Byram
  if(!NumericVector::is_na(htcSI[0])) FL_crown = FL_crown + 0.5*htcSI[0]; 
  I_r_surf = I_r_surf*11.3484; //Btu/ft2/min to kJ/m2/min
  I_r_litter = I_r_litter*11.3484; //Btu/ft2/min to kJ/m2/min
  I_r_ca = I_r_ca*11.3484; //Btu/ft2/min to kJ/m2/min
  I_r_crown = I_r_crown*11.3484; //Btu/ft2/min to kJ/m2/min
  q_surf = q_surf* 37.2589458; //Btu/ft3 to kJ/m3
  q_litter = q_litter* 37.2589458;//Btu/ft3 to kJ/m3
  q_ca = q_ca* 37.2589458;//Btu/ft3 to kJ/m3
  q_crown = q_crown* 37.2589458;//Btu/ft3 to kJ/m3
  I_b = I_b * 0.0576911555; //Btu/ft/min to kW/m
  I_b_crown = I_b_crown * 0.0576911555; //Btu/ft/min to kW/m
  
  //Build output list
    
  List firePotentials=List::create(_["RP"] = RP,
                                   _["SP"] = SP,
                                   _["FP"] = FP,
                                   _["SFP"] = SFP,
                                   _["IC"] = IC,
                                   _["TC"] = TC,
                                   _["RC"] = RC,
                                   _["CFP"] = CFP);

  List crownFire=List::create(_["I_R_canopy [kJ/m2/min]"]=I_r_ca,
                              _["I_R_crown [kJ/m2/min]"]=I_r_crown,
                              _["q_canopy [kJ/m2]"]=q_ca,
                              _["q_crown [kJ/m2]"]=q_crown,
                              _["xi_crown"]=xi_crown,
                              _["canopy_WindSpeed [m/s]"] = crownfireWindSpeedSI,
                              _["WAF"]=WAF,
                              _["ROS_crown [m/min]"]=ros_crown,
                              _["I_b_crown [kW/m]"] = I_b_crown,
                              _["Ic_ratio"] = ic_ratio,
                              _["FL_crown [m]"]=FL_crown);
    
  List surfaceFire=List::create(_["midflame_WindSpeed [m/s]"] = midflameWindSpeedSI,
                                _["phi_wind"] = phi_wind,
                           _["phi_slope"] = phi_slope,
                           _["I_R_surf [kJ/m2/min]"]=I_r_surf,
                           _["I_R_litter [kJ/m2/min]"]=I_r_litter,
                           _["q_surf [kJ/m2]"]=q_surf,
                           _["q_litter [kJ/m2]"]=q_litter,
                           _["xi_surf"]=xi_surf,
                           _["xi_litter"]=xi_litter,
                           _["ROS_surf [m/min]"]=ros_surf,
                           _["ROS_litter [m/min]"]=ros_litter,
                           _["ROS_windslopecap [m/min]"]=windslopecap,
                           _["ROS [m/min]"]=ros,
                           _["I_b [kW/m]"] = I_b,
                           _["FL [m]"]=FL);
  return(List::create(_["SurfaceFire"] = surfaceFire, _["CrownFire"] = crownFire, _["FirePotentials"] = firePotentials));
}

/** ROTHERMEL
 *  Recodified from function 'ros' in package 'Rothermel' (Vacchiano & Ascoli)
 *  Fuel classes are: 1-hour, 10-hour, 100-hour, live herbs and live woody
 * 
 *  modeltype: 'S'(tatic) or 'D'(ynamic) 
 *  wSI: vector of fuel load (t/ha) for five fuel classes
 *  sSI: vector of surface-to-volume ratio (m2/m3) for five fuel classes
 *  delta: a value of fuel bed depth (cm)
 *  mx_dead: a value of dead fuel moisture of extinction (percent)
 *  hSI: a vector of heat content (kJ/kg) for five fuel classes
 *  mSI: a vector of percent moisture on a dry weight basis (percent) for five fuel classes
 *  u: a value of windspeed (m/s) at midflame height. 
 *  windDir: wind direction (in degrees from north). North means blowing from north to south
 *  slope: a value of site slope (in degrees). Negative values are possible.
 *  aspect: aspect (in degrees from north)
 * 
 */
// [[Rcpp::export("fire_Rothermel")]]
List rothermel(String modeltype, NumericVector wSI, NumericVector sSI, double delta, double mx_dead,
                  NumericVector hSI, NumericVector mSI, double u, double windDir, double slope, double aspect) {
  //Rescale variables to units of rothermel model
  NumericVector m = mSI/100.0; //from percent to proportions
  NumericVector w = wSI*0.02048161; //from t/ha to lbs/ft2
  NumericVector s = sSI/3.281; //from m-1 to ft-1
  NumericVector h = hSI*0.429922614; //from kJ/kg to ?
  delta = delta*0.0328084; //from cm to feet
  mx_dead = mx_dead/100.0;//from percent to proportions
  u=u*196.8504; //from m/s to ft/min
  slope = tan(slope*180.0/PI)/100.0; //from degrees to proportions

  // transfer partially cured herbaceous fuels to dead
  if (modeltype=="D") {
    double kt = 0.0;
    if((m[3]>=0.3) & (m[3]<1.2)) kt = (1.20-m[3])/0.9;
    // weighting SAV from transferred cured herbaceous
    double f1=w[0]*s[0]/32.0;
    double f4 = w[3]*kt*s[3]/32.0;
    s[0]=(f1*s[0]+f4*s[3])/(f1+f4);        
    if((f1+f4)==0.0) s[0] = 0.0;
    w[0]=w[0]+w[3]*kt;
    w[3]=w[3]-w[3]*kt;    
  }
  
  //Constants 
  double rhop= 32.0; // = 513*0.0624279606 Scott and Burgan (2005)
  double st=0.0555;
  double se=0.01;
  
  //Area fractions
  NumericVector a = s*w;
  a = a/rhop;
  double a_dead=a[0]+a[1]+a[2];
  double a_live=a[3]+a[4];
  double a_tot=(a_dead+a_live);
  NumericVector f = NumericVector::create(a[0]/a_dead,a[1]/a_dead,a[2]/a_dead,a[3]/a_live,a[4]/a_live);
  if(a_live==0.0) {
    f[3] = 0.0;
    f[4] = 0.0;
  }
  if(a_dead==0.0) {
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
  }
  double f_dead=a_dead/a_tot;
  double f_live=a_live/a_tot;
  if(a_tot==0.0) {
    f_dead = 0.0;
    f_live = 0.0;
  }
  
  //net (weighted) fuel loadings
  NumericVector wn=w*(1.0-st); // Albini 1976
  double wn_dead=f[0]*wn[0]+f[1]*wn[1]+f[2]*wn[2];
  double wn_live=wn[3]+wn[4]; // corrects models w/ 2 live fuel classes  (undocumented)  
  //weighted fuel moisture
  double mf_dead=f[0]*m[0]+f[1]*m[1]+f[2]*m[2];
  double mf_live=f[3]*m[3]+f[4]*m[4];  
  //weighted SAV ratio
  double sigma_dead=f[0]*s[0]+f[1]*s[1]+f[2]*s[2];
  double sigma_live=f[3]*s[3]+f[4]*s[4];
  double sigma_tot=(f_dead*sigma_dead+f_live*sigma_live); //characteristic SAV
  
  //weighted heat content
  double h_dead=f[0]*h[0]+f[1]*h[1]+f[2]*h[2];
  double h_live=f[3]*h[3]+f[4]*h[4];
  //mean packing ratio for fuel complex
  double beta=(1.0/delta)*(w[0]/rhop+w[1]/rhop+w[2]/rhop+w[3]/rhop+w[4]/rhop);
  if(delta==0.0) beta = 0.0;
  
  //live fuel moisture of extinction
  double Wden = (w[3]*exp(-500.0/s[3])+w[4]*exp(-500.0/s[4]));
  double W=(w[0]*exp(-138.0/s[0])+w[1]*exp(-138.0/s[1])+w[2]*exp(-138.0/s[2]))/Wden;
  if(Wden==0.0) W = 0.0;
  double den = (w[0]*exp(-138/s[0])+w[1]*exp(-138.0/s[1])+w[2]*exp(-138.0/s[2]));
  double mfpd=(w[0]*m[0]*exp(-138.0/s[0])+w[1]*m[1]*exp(-138.0/s[1])+w[2]*m[2]*exp(-138.0/s[2]))/den;
  if(den==0.0) mfpd = 0.0;
  double mx_live=std::max(2.9*W*(1.0-mfpd/mx_dead)-0.226, mx_dead);
  if(Wden==0) mx_live = mx_dead;
//  Rcout<<mx_dead<<", "<< mfpd<<", "<< mx_live<<"\n";

  //damping coefficients
  double ns=0.174*pow(se,-0.19);
  if(se==0.0) ns = 0.0;
  double nm_dead=1-2.59*(mf_dead/mx_dead)+5.11*pow(mf_dead/mx_dead,2.0)-3.52*pow(mf_dead/mx_dead,3.0);
  double nm_live=1-2.59*(mf_live/mx_live)+5.11*pow(mf_live/mx_live,2.0)-3.52*pow(mf_live/mx_live,3.0);
  if(mx_dead==0.0) nm_dead = 0.0;
  if(mx_live==0.0) nm_live = 0.0;
  if(mf_dead>mx_dead) nm_dead = 0.0;  
  if(mf_live>mx_live) nm_live = 0.0;  
//  Rcout<<nm_dead<<", "<< nm_live<<"\n";
  
  // Andrews 2013 pag.E
  //optimum packing ratio
  double beta_op=3.348*pow(sigma_tot,-0.8189);
  if(sigma_tot==0.0) beta_op = 0.0;
//  Rcout<<beta_op<<"\n";

  //relative packing ratio
  double rpr=beta/beta_op; 
  if(beta_op==0.0) rpr = 0.0;
//  Rcout<<rpr<<"\n";

  //maximum reaction velocity
  double gamma_max=pow(sigma_tot,1.5)/(495.0+0.0594*pow(sigma_tot,1.5));
  if(sigma_tot==0.0) gamma_max = 0.0;
//  Rcout<<gamma_max<<"\n";
  
  //reaction intensity (in kJ/ft2/min)
  double sum_dead=(wn_dead*h_dead*nm_dead*ns);
  double sum_live=(wn_live*h_live*nm_live*ns);
  //alternate formulation from Albini (1976)
  double A=133.0*pow(sigma_tot,-0.7913);
  if(sigma_tot==0.0) A = 0.0;
  double ir_dead=gamma_max*pow(rpr*exp(1.0-rpr), A)*sum_dead; //*f.dead removed by Frandsen 73
  double ir_live=gamma_max*pow(rpr*exp(1.0-rpr), A)*sum_live; //*f.live removed by Frandsen 73
  if(rpr==0.0) {
    ir_dead=0.0;
    ir_live=0.0;
  }
  double ir=ir_dead+ir_live;
//  Rcout<<sum_dead<<", "<<sum_live<<"\n";
//  Rcout<<ir<<"\n";
  
  //propagating flux ratio
  double xi=pow(192+0.2595*sigma_tot,-1.0)*exp((0.792+0.681*sqrt(sigma_tot))*(beta+0.1));  
//  Rcout<<xi<<"\n";

//wind coefficient  
  double C=7.47*exp(-0.133*pow(sigma_tot,0.55));
  if(sigma_tot==0.0) C = 0.0;
  double B=0.02526*pow(sigma_tot,0.54);
  if(sigma_tot==0.0) B = 0.0;
  double E=0.715*exp(-0.000359*sigma_tot);
  double fw=C*pow(u,B)*pow(rpr,-E);
  if(rpr==0.0) fw = 0.0;
  
  //slope coefficient
  double fs=5.275*pow(beta,-0.3)*pow(slope,2.0);
  if(beta==0.0) fs = 0.0;
  
  //calculate 'effective wind' direction and spread
  NumericVector vS = NumericVector::create(fs, (aspect*PI/180.0)+PI); //add 180 degrees for pushing fire in the direction opposite of the aspect
  NumericVector vW = NumericVector::create(fw, (windDir*PI/180.0)+PI); //add 180 degrees for blowing 'in the direction of'
  NumericVector vEff = vectorAddition(vS,vW);
//  Rcout << " fw " << fw << " dir " << (windDir+3.141592)*(180.0/3.141592) << " fs " << fs << " dir " << (slopeDir+3.141592)*(180.0/3.141592) << "res: " << vEff[0] <<","<<vEff[1]*(180.0/3.141592)<<"\n";
  
  //resulting effect and virtual wind speed
  double phi = vEff[0];
  double vws = pow(phi/(C*pow(rpr,-E)), 1.0/B);
  if((B==0.0) |(rpr==0.0)) vws = 0.0;
  vws = vws/196.8504; //from ft/min to m/s
//  Rcout<<vws<<"\n";

  //for heat sink (denominator ROS)
  double rhob=(w[0]+w[1]+w[2]+w[3]+w[4])/delta; 
  if(delta==0.0) rhob = 0.0;

  NumericVector qig=250+1116*m;
  for(int i=0;i<5;i++) if(w[i]==0.0) qig[i]=0.0;
  
  double eps=f_dead*(f[0]*qig[0]*exp(-138.0/s[0])+f[1]*qig[1]*exp(-138.0/s[1])+f[2]*qig[2]*exp(-138.0/s[2]))+f_live*(f[3]*qig[3]*exp(-138.0/s[3])+f[4]*qig[4]*exp(-138.0/s[4]));
//  Rcout<<eps<<"\n";

  //ROS
  double r = (ir*xi*(1+phi))/(rhob*eps);
  if(rhob==0.0) r = 0.0;
  if(eps==0.0) r = 0.0;
//  Rcout<<r<<"\n";  
  r = 0.3048*r; //from feet/min to m/min

  //Limit the rate of spread to be smaller or equal the effective wind speed
  //Andrews, P.L., Cruz, M.G., Rothermel, R.C., 2013. Examination of the wind speed limit function in the Rothermel surface fire spread model. Int. J. Wildl. Fire 22, 959–969. doi:10.1071/WF12122
  r = std::min(r, vws*60.0);
  
  
  double tr = 1259.843/(sigma_tot*3.281); //Residence time (in min) Anderson (1969)
  double Ib = (ir*0.1893)*tr*r;//Fireline intensity (kW/m)
  
  //return values for rothermel function
  List output=List::create(Named("Characteristic dead fuel moisture [%]")=mf_dead*100.0,
                           Named("Characteristic live fuel moisture [%]")=mf_live*100.0,
                           Named("Live fuel moisture of extinction [%]")=mx_live*100.0,
                           Named("Characteristic SAV [m2/m3]")=sigma_tot*3.281,
                           Named("Bulk density [kg/m3]") = rhob*16.0184634,
                           Named("Packing ratio [dimensionless]")=beta,
                           Named("Relative packing ratio [dimensionless]")=rpr,
                           Named("Dead fuel Reaction intensity [kW/m2]") = ir_dead* 0.1893,
                           Named("Live fuel Reaction intensity [kW/m2]") = ir_live* 0.1893,
                           Named("Reaction intensity [kW/m2]") = ir*0.1893,
                           Named("Fireline intensity [kW/m]") = Ib,
                           Named("Wind factor [dimensionless]") = fw,
                           Named("Slope factor [dimensionless]") = fs,
                           Named("Slope-wind vector") = vEff,
                           Named("Virtual wind speed [m/s]") = vws,
                           Named("Heat source [kW/m2]") = (ir*xi*(1+phi))*0.1893,
                           Named("Heat sink [kJ/m3]") = (rhob*eps)*37.25894580781,
                           Named("ROS [m/min]") = r);
  return(output);
}


/** FIRE BEHAVIOUR
 * 
 *  x: An object of class 'forest'
 *  soil: An object of class 'soil'
 *  slope: A value of site slope (in degrees). Negative values are possible.
 *  aspect: aspect (in degrees from north)
 *  meteo: A data frame with meteorological variables for each day.
 *  SpParams: A data frame with species parameters
 *  surfaceToVolumeRatios: A vector of surface-to-volume ratio (m2/m3) for five fuel classes
 *  heatContent: A vector of heat content (kJ/kg) for five fuel classes
 *  deadFuelMoistureExtinction: a value of dead fuel moisture of extinction (percent)
 *  control: A list with default parameter values
 */
// List fb(List x, List soil, double latitude, double slope, double aspect, DataFrame meteo, DataFrame SpParams, 
//         DataFrame FuelModelParams, List control) {
//   bool verbose = control["verbose"];
//   String liveFMCmode = control["liveFMCmode"];
//   bool useModelForLive = control["useModelForLive"]; 
//   if(verbose) Rcout<<"Initializing";
//   
//   //Prepare swb input
//   DataFrame swbInput = forest2swbInput(x, soil, SpParams, control);
//   NumericVector pBole = swbInput["pBole"];
//   NumericVector H = swbInput["H"];
//   NumericVector Sgdd = swbInput["Sgdd"];
//   IntegerVector SP = swbInput["SP"];
//   NumericVector LAI = swbInput["LAI"];
//   
//   NumericVector cohLoading = cohortFuel(x, SpParams);
// 
//   //Initialize other cohort-based variables
//   int numCohorts = LAI.size();
//   NumericVector LAIphe(numCohorts), cohFMC(numCohorts), cohLoadingPhe(numCohorts);
//   
//   //Meteorological variables
//   IntegerVector DOY = meteo["DOY"]; // day of the year
//   NumericVector Precipitation = meteo["Precipitation"]; // in mm of water (L/m2)
//   NumericVector MeanTemperature = meteo["MeanTemperature"]; // in degrees Celsius
//   NumericVector MinTemperature = meteo["MinTemperature"]; // in degrees Celsius
//   NumericVector MaxTemperature = meteo["MaxTemperature"]; // in degrees Celsius
//   NumericVector RHmean = meteo["MeanRelativeHumidity"]; // in percentage
//   NumericVector RHmin = meteo["MinRelativeHumidity"]; // in percentage
//   NumericVector RHmax = meteo["MaxRelativeHumidity"]; // in percentage
//   NumericVector WS = meteo["WindSpeed"]; // in m/s
//   NumericVector WD = meteo["WindDirection"]; // in radians from north
//   NumericVector Rn = meteo["Radiation"];  //in MJ/m2
//   NumericVector PET = meteo["PET"]; //in mm of water
// 
//   NumericVector GDD = gdd(DOY, MeanTemperature, 5.0);
//   NumericVector ER = er(DOY);
// 
//   int numDays = Precipitation.size();
//   
//   //Output FS variables
//   NumericVector fuelbedHeight(numDays), fine1hLoading(numDays), 
//                 coarse10hLoading(numDays), coarse100hLoading(numDays), 
//                 herbaceousLoading(numDays), woodyLoading(numDays),
//                 canopyBaseHeight(numDays), canopyTopHeight(numDays),
//                 canopyLength(numDays), canopyBulkDensity(numDays),
//                 canopyLAI(numDays);
// 
//   //Output FM variables
//   NumericVector fuelWind(numDays), fine1hFMC(numDays), coarse10hFMC(numDays), 
//                 coarse100hFMC(numDays), canopyFMC(numDays), fuelbedWoodyFMC(numDays),
//                 fuelbedHerbaceousFMC(numDays);
// 
//   //Output FB variables
//   IntegerVector burningType(numDays);
//   NumericVector firelineIntensity(numDays), spreadRate(numDays), surfaceSpreadRate(numDays), 
//                 midflameWind(numDays), virtualWindSpeed(numDays), spreadDirection(numDays), I0(numDays);
//   
//   //Initialize fuel moisture
//   NumericVector fMoisture(5);
//   fMoisture[0] = x["dead1hourMoisture"];
//   fMoisture[1] = x["dead10hourMoisture"];
//   fMoisture[2] = x["dead100hourMoisture"];
//   
//   NumericVector fcafternoonprev; //Stores fuel conditions of afternoon of previous day
//   if(verbose) Rcout<<" - Daily loop ";
//   for(int i=0;i<numDays;i++) {
//     //SOIL WATER BALANCE
//     if(verbose) Rcout<<".";
//     List s = swbDay1(swbInput, soil, GDD[i], PET[i], Precipitation[i], ER[i], 0.0); //No Runon in simulations for a single cell    
//     double Lground = s["Lground"]; //Percent light at the ground level
// 
//     //FUEL STRUCTURE
//     NumericVector Phe = leafDevelopmentStatus(Sgdd, GDD[i]);
//     for(int c=0;c<numCohorts;c++) {
//       if(Sgdd[c]==0.0) Phe[c]=1.0;
//       LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology 
//       cohLoadingPhe[c] = cohLoading[c]*Phe[c]; //LAI modified by phenology 
//     }
//     //Loadings and FMC by cohort
//     cohFMC = cohortFuelMoistureContent(s, swbInput, SpParams);
//     
//     List fs = fuelStructure(x, SpParams, FuelModelParams, GDD[i], 30, 5000, 0.01, useModelForLive);
//     fuelbedHeight[i] = fs["fuelbedHeight [cm]"];
//     NumericVector fBiomass = fs["fuelbedLoading [Mg/ha]"];
//     NumericVector surfaceToVolumeRatios = fs["fuelbedSAV [m2/m3]"];
//     fine1hLoading[i] = fBiomass[0];
//     coarse10hLoading[i] = fBiomass[1];
//     coarse100hLoading[i] = fBiomass[2];
//     herbaceousLoading[i] = fBiomass[3];
//     woodyLoading[i] = fBiomass[4];
//     canopyBaseHeight[i] = fs["canopyBaseHeight [cm]"];
//     canopyTopHeight[i] = fs["canopyTopHeight [cm]"];
//     canopyLength[i] = fs["canopyLength [cm]"];
//     canopyBulkDensity[i] = fs["canopyBulkDensity [kg/m3]"];
//     canopyLAI[i] = fs["canopyLAI [dimensionless]"];
// 
//     //WIND EXTINCTION
//     //Check NA
//     if(NumericVector::is_na(WS[i])) WS[i] = 2.0; //Default value
//     if(NumericVector::is_na(WD[i])) WD[i] = 180.0; //Wind from south
//     
//     //Calculates midflameWind
//     if(NumericVector::is_na(canopyBaseHeight[i])) { //NO canopy over fire
//        midflameWind[i] = unshelteredMidflameWindSpeed(WS[i], fuelbedHeight[i]/100.0); 
//        fuelWind[i] = windSpeedAtCanopyHeight(WS[i], fuelbedHeight[i]/100.0);
//     } else { //Wind extinction due to crowns
//        midflameWind[i] = windSpeedMassmanExtinction(canopyBaseHeight[i]/100.0,WS[i], canopyLAI[i], canopyTopHeight[i]/100.0);
//        fuelWind[i] = midflameWind[i];
//     }    
// 
//     //FUEL MOISTURE
//     double netPrec = s["NetPrec"];
//     double rainDuration = std::min(24.0,Precipitation[i]/6.35);
//     
//     //Number of hours of insolation
//     NumericVector srs = meteoland::radiation_sunRiseSet(latitude*(PI/180.0),  slope*(PI/180.0), aspect*(PI/180.0),DOY[i]);
//     double n = std::max(0.0, (srs[1] - srs[0])*12.0/PI); 
//     if(n<=0.0) n = 0.0;
//     
//     //Dead fuel moisture
//     double radSec = Rn[i]*(Lground/100.0)*(1000000.0/(n*3600.0)); //Transform to MJ/m2.day to W/m2 assuming n hours of equal sunlight
//     if(n<=0.0) radSec = 0.0;
// 
//     NumericVector fcmean = fuelConditions(MeanTemperature[i],RHmean[i], radSec, fuelWind[i]);
//     NumericVector fcmorning = fuelConditions(MinTemperature[i],RHmax[i], radSec, fuelWind[i]);
//     NumericVector fcafternoon = fuelConditions(MaxTemperature[i],RHmin[i], radSec, fuelWind[i]);
//     if(i==0) fcafternoonprev = fcafternoon; 
//     fine1hFMC[i] = fine1hday(fMoisture[0], fcmean[0], fcmean[1], fuelWind[i], netPrec);
//     coarse10hFMC[i] = coarse10hday(fMoisture[1], 
//                                    (fcafternoonprev[0] +fcmorning[0])/2.0, (fcafternoonprev[1] + fcmorning[1])/2.0,
//                                    (fcmorning[0]+fcafternoon[0])/2.0, (fcmorning[1]+fcafternoon[1])/2.0, rainDuration);
//     coarse100hFMC[i] = coarse100hday(fMoisture[2], 
//                                      fcmorning[0], fcmorning[1],
//                                      fcafternoon[0], fcafternoon[1], n, rainDuration);
//     //Store for next day (and for fire behaviour)
//     fMoisture[0] = fine1hFMC[i]; 
//     fMoisture[1] =  coarse10hFMC[i]; 
//     fMoisture[2] = coarse100hFMC[i]; 
//     fcafternoonprev = fcafternoon;
//       
//     //Live fuel moisture
//     NumericVector psi = s["psiVec"];
//     if(liveFMCmode == "constant") {
//       fMoisture[3] = 200.0; 
//       fMoisture[4] = 100.0; 
//       canopyFMC[i] = 100.0;
//     } else if(liveFMCmode =="ffmc") {
//       fMoisture[3] = fine1hFMC[i]; 
//       fMoisture[4] = fine1hFMC[i]; 
//       canopyFMC[i] = fine1hFMC[i];
//     } else if(liveFMCmode =="swb") {
//       //Level of physiological activity is defined in terms of growth degree days
//       double act = leafDevelopmentStatus(200.0, GDD[i]);
//       //In winter, herbaceous FMC is modelled as fine (1h) dead fuels
//       //When physiologically active, FMC is modelled as function of water content in the topsoil
//       fMoisture[3] = act*std::max(30.0, 200.0+20*(psi[0]/1000.0))+(1.0-act)*fine1hFMC[i]; 
//       
//       fMoisture[4] = fuelbedLiveFuelMoisture(fuelbedHeight[i], cohFMC, cohLoadingPhe, H, pBole); 
//       canopyFMC[i] = canopyLiveFuelMoisture(canopyBaseHeight[i], canopyTopHeight[i], cohFMC, cohLoadingPhe, H, pBole); 
//     } else {
//       stop("Wrong 'liveFMCmode'");
//     }
//     
//     fuelbedHerbaceousFMC[i] = fMoisture[3];
//     fuelbedWoodyFMC[i] = fMoisture[4];
// 
// 
//    CharacterVector FMcode = x["FuelModelCode"];
//    CharacterVector models = FuelModelParams.attr("row.names");
//    NumericVector heatContent(5);
//    double deadFuelMoistureExtinction = 30.0; //Default value
//    for(int j=0;j<models.length();j++) {
//      if(models[j]==FMcode[0]) {
//        heatContent[0] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_1h"])[j];
//        heatContent[1] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_10h"])[j];
//        heatContent[2] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_100h"])[j];
//        heatContent[3] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Herb"])[j];
//        heatContent[4] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Woody"])[j];
//        deadFuelMoistureExtinction = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Mx_dead"])[j];
//      }
//    }
//    
//     //FIRE BEHAVIOUR
//     //1. Calculate surface fire characteristics
//     //1.2 Call to Rothermel model to obtain surface fire spread rate and intensity
//     List rtm = rothermel("D", fBiomass, surfaceToVolumeRatios, fuelbedHeight[i],
//                         deadFuelMoistureExtinction, heatContent, fMoisture, 
//                         midflameWind[i], WD[i], 
//                         slope, aspect);
//     surfaceSpreadRate[i] = rtm["ROS [m/min]"];
//     spreadRate[i] = surfaceSpreadRate[i];
//     firelineIntensity[i] = rtm["Fireline intensity [kW/m]"];
//     virtualWindSpeed[i] = rtm["Virtual wind speed [m/s]"];
//     NumericVector slopeWindVec = rtm["Slope-wind vector"];
//     
//     //1.3 Determine whether crowning will occur  (Van Wagner 1989)
//     if(!NumericVector::is_na(canopyBaseHeight[i])) {
//       I0[i] = criticalFirelineIntensity(canopyBaseHeight[i]/100.0, canopyFMC[i]); //critical fireline intensity (MW/m)
//       //Is crown fire (active or passive)?
//       if(firelineIntensity[i] >= I0[i]) { 
//         double RAC = 3.0/canopyBulkDensity[i]; // Threshold for active crown spread rate (m/min)
//         double Ei = 1.0;
//         NumericVector fBiomass10 = NumericVector::create(0.38, 0.092, 0.230,0.0, 0.092); //lb/square feet
//         fBiomass10 = fBiomass10*(1.0/(0.3048*0.3048))*0.45359237*(10.0); // transform to tons/hectare
//         List rtm10 = rothermel("D", fBiomass10, surfaceToVolumeRatios, 30.48,
//                         25, heatContent, fMoisture, WS[i]*0.4, WD[i], 
//                         slope, aspect);
//         double R10 = rtm10["ROS [m/min]"]; //Spread rate with fuel model 10 (Rothermel 1991)
//         double R0 = (I0[i]/firelineIntensity[i])*surfaceSpreadRate[i]; //critical surface spread rate (m/s) associated with critical fireline intensity
//         double ac = -log(0.1)/(0.9*(RAC-R0)); //Scaling coefficient
//         double cFractionBurn = 1.0 - exp(-ac*(surfaceSpreadRate[i]-R0)); //Fraction of canopy burn with crowning
//         double Rcmax = (3.34)*R10*Ei; // Maximum crown spread rate (m/s)
//         double Rcactual = surfaceSpreadRate[i] + cFractionBurn*(Rcmax - surfaceSpreadRate[i]); //Actual crown spread rate (m/s)
//         if(Rcactual >= RAC) {
//           burningType[i] = 4;
//           spreadRate[i] = Rcactual;
//         } else {
//           burningType[i] = 3;
//         }
//         //Modify fireline intensity for crown fire
//         firelineIntensity[i] = 300.0 * ((firelineIntensity[i]/(300.0*surfaceSpreadRate[i])) + cFractionBurn*canopyBulkDensity[i]*(canopyTopHeight[i]-canopyBaseHeight[i])*0.01)*spreadRate[i];        
//       } else {
//         burningType[i] = 2; // Surface fire
//       }
//     } else {
//       burningType[i] = 0; // Shrubland (open) fire
//     }
// 
//     spreadDirection[i] = slopeWindVec[1]*(180.0/3.141592);  
//     if(spreadDirection[i]<0) spreadDirection[i] +=360.0; 
//   }
//   
//   if(verbose) Rcout<<"Building output ...";  
//   Rcpp::DataFrame FS = DataFrame::create(_["fuelbedHeight"] = fuelbedHeight,
//                       _["fine1hLoading"] = fine1hLoading,
//                       _["coarse10hLoading"] = coarse10hLoading,
//                       _["coarse100hLoading"] = coarse100hLoading,
//                       _["herbaceousLoading"] = herbaceousLoading,
//                       _["woodyLoading"] = woodyLoading,
//                       _["canopyBaseHeight"] = canopyBaseHeight,
//                       _["canopyTopHeight"] = canopyTopHeight,
//                       _["canopyLength"] = canopyLength,
//                       _["canopyBulkDensity"] = canopyBulkDensity,
//                       _["canopyLAI"] = canopyLAI);
//   Rcpp::DataFrame FM = DataFrame::create(_["fuelWind"]=fuelWind, 
//                                          _["fine1hFMC"]=fine1hFMC, 
//                                          _["coarse10hFMC"]=coarse10hFMC,
//                                          _["coarse100hFMC"]=coarse100hFMC,
//                                          _["liveHerbaceousFMC"] = fuelbedHerbaceousFMC,
//                                          _["liveWoodyFMC"] = fuelbedWoodyFMC,
//                                          _["canopyFMC"] = canopyFMC);
// 
//   Rcpp::DataFrame FB = DataFrame::create(_["burningType"]=burningType,
//                                          _["midflameWind"] = midflameWind,
//                                          _["virtualWindSpeed"]=virtualWindSpeed,
//                                          _["spreadDirection"]=spreadDirection,
//                                          _["surfaceSpreadRate"]=surfaceSpreadRate,
//                                          _["intensityCrowningThreshold"]=I0,
//                                          _["spreadRate"]=spreadRate,
//                                          _["firelineIntensity"]=firelineIntensity);
//   
//   FS.attr("row.names") = meteo.attr("row.names");
//   FM.attr("row.names") = meteo.attr("row.names");
//   FB.attr("row.names") = meteo.attr("row.names");
//   List l = List::create(Named("FuelStructure")=FS, Named("FuelMoisture")=FM, Named("Behaviour")=FB);
//   l.attr("class") = CharacterVector::create("fb","list");
//   if(verbose) Rcout<<"done.\n";
//   return(l);  
// }
