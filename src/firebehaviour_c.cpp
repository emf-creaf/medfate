#include "RcppArmadillo.h"
#include "modelInput_c.h"
#include "lowlevel_structures_c.h"
#include "firebehaviour_c.h"
#include "fuelstructure_c.h"
#include "windextinction_c.h"
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"
using namespace Rcpp;

Rcpp::NumericVector copyFCCSResult_c(const FCCS_RESULT& fccs) {
  Rcpp::NumericVector fireHazard = Rcpp::NumericVector::create(
    Rcpp::Named("Loading_overstory [kg/m2]") = fccs.Loading_overstory,
    Rcpp::Named("Loading_understory [kg/m2]") = fccs.Loading_understory,
    Rcpp::Named("CFMC_overstory [%]") = fccs.CFMC_overstory,
    Rcpp::Named("CFMC_understory [%]") = fccs.CFMC_understory,
    Rcpp::Named("DFMC [%]") = fccs.DFMC,
    Rcpp::Named("ROS_surface [m/min]") = fccs.ROS_surface,
    Rcpp::Named("I_b_surface [kW/m]") = fccs.I_b_surface,
    Rcpp::Named("t_r_surface [s]") = fccs.t_r_surface,
    Rcpp::Named("FL_surface [m]") = fccs.FL_surface,
    Rcpp::Named("Ic_ratio") = fccs.Ic_ratio,
    Rcpp::Named("ROS_crown [m/min]") = fccs.ROS_crown,
    Rcpp::Named("I_b_crown [kW/m]") = fccs.I_b_crown,
    Rcpp::Named("t_r_crown [s]") = fccs.t_r_crown,
    Rcpp::Named("FL_crown [m]") = fccs.FL_crown,
    Rcpp::Named("SFP") = fccs.SFP,
    Rcpp::Named("CFP") = fccs.CFP);
  return(fireHazard);
}

/**
 * From van Wagner (1977) model
 * CBH - crown base height (in m)
 * M - Crown foliar moisture (in percent)
 */
// [[Rcpp::export(".criticalFirelineIntensity")]]
double criticalFirelineIntensity_c(double CBH, double M) {
  return(pow(0.010*(CBH)*(460.0+25.9*M),1.5));
}

void FCCSbehaviour_c(FCCSBehaviour_RESULT& res,
                     const InternalFCCS& FCCSpropsSI,
                     std::vector<double>& MliveSI, 
                     std::vector<double>& MdeadSI, 
                     double slope, double windSpeedSI) {
  //Extract vectors
  const std::vector<double>& wSI = FCCSpropsSI.w;
  const std::vector<double>&  cover = FCCSpropsSI.cover;
  const std::vector<double>&  hbcSI = FCCSpropsSI.hbc;
  const std::vector<double>&  htcSI = FCCSpropsSI.htc;
  const std::vector<double>&  deltaSI = FCCSpropsSI.delta;
  const std::vector<double>&  rhopSI = FCCSpropsSI.rhop;
  const std::vector<double>&  PVSI = FCCSpropsSI.PV;
  const std::vector<double>&  beta = FCCSpropsSI.beta; //unitless
  const std::vector<double>&  betarel = FCCSpropsSI.betarel; //unitless
  const std::vector<double>&  etabetarel = FCCSpropsSI.etabetarel; //unitless
  const std::vector<double>&  sigmaSI = FCCSpropsSI.sigma;
  const std::vector<double>&  pDead = FCCSpropsSI.pDead;
  const std::vector<double>&  FAI = FCCSpropsSI.FAI; //unitless
  const std::vector<double>&  hdefSI = FCCSpropsSI.h;
  const std::vector<double>&  RVSI = FCCSpropsSI.RV;
  const std::vector<double>&  ActFMC = FCCSpropsSI.ActFMC;
  
  //Replace fuel moisture if available
  if(!std::isnan(ActFMC[0])) MliveSI[0] = ActFMC[0];
  if(!std::isnan(ActFMC[1])) MliveSI[1] = ActFMC[1];
  if(!std::isnan(ActFMC[2])) MliveSI[2] = ActFMC[2];
  
  //Rescale variables to British units
  double* Mlive = new double[5];
  double* Mdead = new double[5];
  double* w = new double[5];
  double* delta = new double[5];
  double* rhop = new double[5];
  double* PV = new double[5];
  double* sigma = new double[5];
  double* hdef = new double[5];
  double* RV = new double[5];
  for(int i=0;i<5;i++) {
    Mlive[i] = MliveSI[i]/100.0; //from percent to proportions
    Mdead[i] = MdeadSI[i]/100.0; //from percent to proportions
    w[i] = wSI[i] * 0.204918; //from kg/m2 to lb/ft2
    delta[i] = deltaSI[i] * 3.2808399; //from m to ft
    rhop[i] = rhopSI[i] * 0.06242796; //from kg/m3 to lb/ft3
    PV[i] = PVSI[i] * 3.2808399; //from m3/m2 to ft3/ft2
    sigma[i] = sigmaSI[i] * 0.3048; //from m2/m3 to ft2/ft3
    hdef[i] = hdefSI[i] * 0.429922614; //from kJ/kg to Btu/lb
    RV[i] = RVSI[i] * 3.2808399; //from m3/m2 to ft3/ft2
  }
  
  //Heat content (after correcting for moisture in live fuels)
  double* h = new double[5];
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
  double* etaMlive = new double[3];
  double* etaMdead = new double[3];
  double* etaM = new double[5];
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
  double I_r_sh = pow(etabetarel[1], A_sh)*Gmax_sh*w[1]*h[1]*etaM[1]*etaK;
  double I_r_he = pow(etabetarel[2], A_he)*Gmax_he*w[2]*h[2]*etaM[2]*etaK;
  double I_r_wo = pow(etabetarel[3], A_wo)*Gmax_wo*w[3]*h[3]*etaM[3]*etaK;
  double I_r_li = pow(etabetarel[4], A_li)*Gmax_li*w[4]*h[4]*etaM[4]*etaK;
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
  if(std::isnan(htcSI[0])) crownFillProportion=0.0;
  double midflameWindSpeedSI = NA_REAL;
  //unsheltered vs sheltered (5% of crownFillProportion)
  //from m/s to mph
  if(crownFillProportion<0.05) midflameWindSpeedSI = unshelteredMidflameWindSpeed_c(windSpeedSI, delta_surf_heatsink/3.2808399);
  else midflameWindSpeedSI = shelteredMidflameWindSpeed_c(windSpeedSI, crownFillProportion, htcSI[0]);
  double crownfireWindSpeedSI = windSpeedAtCanopyHeight_c(windSpeedSI,htcSI[0]);
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
  res.potentials.RP = 0.0;
  res.potentials.SP = 0.0;
  res.potentials.FP = 0.0;
  if(!std::isnan(I_r_surf)) res.potentials.RP = std::min(9.0,0.08*pow(I_r_surf, 0.5));
  if(!std::isnan(ros)) res.potentials.SP = std::min(9.0, 2.5*pow(ros, 0.5));
  if(!std::isnan(FL)) res.potentials.FP = std::min(9.0, 2.5*pow(FL, 0.5));
  res.potentials.SFP = std::min(9.0, std::max(res.potentials.SP,res.potentials.FP));
  
  //Canopy gap
  double GAP_ca = hbcSI[0]-htcSI[1];
  
  //Van Wagner's critical fireline intensity (in Btu/ft/s)
  double I_c = 0.288894658*criticalFirelineIntensity_c(GAP_ca, 100.0*Mlive[0]); //express GAP in m and M in percent before calling Van Wagner
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
  // Rcout << " FAIc "<< FAI[0] << " TFAI/3pi "<< (TFAIc/(3.0*M_PI))<<"\n";
  
  //Crown fire spread rate
  double I_r_ca  = pow(etabetarel[0], A_ca)*Gmax_ca*w[0]*h[0]*etaM[0]*etaK;
  double I_r_crown = I_r_surf+I_r_ca;
  double xi_crown = 1.0 - exp(-1.0*(FAI[0]/(4.0*delta[0])));
  double q_ca = (0.5*FAI[0]*RT_ca*rhop[0]*Qig_ca)/((cover[0]/100.0)*delta[0]);
  double q_crown = q_surf + q_ca;
  double ros_crown = I_r_crown * xi_crown*WAF/q_ca;
  
  //Crown fire residence time (min)
  double t_r_crown = 192.0*RT_ca;
  // Rcout<<"T_r_crown: "<<t_r_crown<<"\n";
  
  //Crown fire fireline intensity (Btu/ft/min)
  double I_b_crown = I_r_crown*ros_crown*t_r_crown;
  
  //Crown flame length (ft)
  double FL_crown = 0.45*pow(I_b_crown/60.0, 0.46);
  
  //Crown fire porentials
  res.potentials.IC = 0.0;
  res.potentials.RC = 0.0;
  res.potentials.TC = 0.0;
  double ic_ratio = NA_REAL;
  if(FAI[0]>0.0) {
    if((!std::isnan(I_b)) && (!std::isnan(I_c))) {
      ic_ratio = (I_b/60.0)/I_c;
      res.potentials.IC = std::min(9.0, 4.0*pow(ic_ratio,0.2));
    }
  } 
  if(FAI[0]> (TFAIc/(3.0*M_PI))) {
    if(!std::isnan(TCq)) res.potentials.TC = std::min(9.0, 10.0*TCq);
  }
  if(!std::isnan(ros_crown)) res.potentials.RC = std::min(9.0, 1.0*pow(ros_crown, 0.5));
  
  // double CFP = 0.4286*(IC+(TC/3.0)+RC);
  res.potentials.CFP = pow(res.potentials.IC*res.potentials.RC, 0.5);
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
  if(!std::isnan(htcSI[0])) FL_crown = FL_crown + 0.5*htcSI[0]; 
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
  

  //Copy crown fire values to structure
  res.crown.I_R_canopy = I_r_ca;
  res.crown.I_R_crown = I_r_crown;
  res.crown.q_canopy = q_ca;
  res.crown.q_crown = q_crown;
  res.crown.xi_crown = xi_crown;
  res.crown.canopy_WindSpeed = crownfireWindSpeedSI;
  res.crown.WAF = WAF;
  res.crown.ROS_crown = ros_crown;
  res.crown.I_b_crown = I_b_crown;
  res.crown.t_r_crown = t_r_crown*60.0; //from min to sec
  res.crown.Ic_ratio = ic_ratio;
  res.crown.FL_crown = FL_crown;

  
  //Copy surface fire values to structure
  res.surface.midflame_WindSpeed = midflameWindSpeedSI;
  res.surface.phi_wind = phi_wind;
  res.surface.phi_slope = phi_slope;
  res.surface.I_R_surf =I_r_surf;
  res.surface.I_R_litter = I_r_litter;
  res.surface.q_surf =q_surf;
  res.surface.q_litter=q_litter;
  res.surface.xi_surf=xi_surf;
  res.surface.xi_litter=xi_litter;
  res.surface.ROS_surf=ros_surf;
  res.surface.ROS_litter=ros_litter;
  res.surface.ROS_windslopecap=windslopecap;
  res.surface.ROS=ros;
  res.surface.I_b = I_b;
  res.surface.t_r = t_r*60.0; //From minutes to sec
  res.surface.FL=FL;
    
  //Free memory for dynamic vectors
  delete[] h;
  delete[] etaMlive;
  delete[] etaMdead;
  delete[] etaM;
  delete[] Mlive;
  delete[] Mdead;
  delete[] w;
  delete[] delta;
  delete[] rhop;
  delete[] PV;
  delete[] sigma;
  delete[] hdef;
  delete[] RV;
}


void fccsHazard_c(FCCSBehaviour_RESULT& FCCSbehres, FCCS_RESULT& FCCSres, ModelInput& x, 
                  const WeatherInputVector& meteovec, 
                  const std::vector<double>& LFMC, 
                  const std::vector<double>& PLC, 
                  const double slope) {
  
  double fireHazardStandardWind = x.control.fireHazard.standardWind;
  double fireHazardStandardDFMC = x.control.fireHazard.standardDFMC; 
  
  InternalFCCS& FCCSprops = x.internalFCCS;

  
  double fm_dead = NA_REAL;
  if(!NumericVector::is_na(fireHazardStandardDFMC)) {
    fm_dead = fireHazardStandardDFMC;
  } else {
    // Estimate moisture of dead fine fuels (Resco de Dios et al. 2015)
    double vp = averageDailyVapourPressure_c(meteovec.tmin, meteovec.tmax, meteovec.rhmin, meteovec.rhmax);
    double D = std::max(0.0, saturationVapourPressure_c(meteovec.tmax) - vp);
    fm_dead = 5.43 + 52.91*exp(-0.64*D); 
  }
  double wind = meteovec.wind;
  if(!std::isnan(fireHazardStandardWind)) wind = fireHazardStandardWind;
  
  int numCohorts = x.above.H.size();
  //Calculate cohort canopy moisture to the average of canopy live and dead fuels, considering that a fraction of LAI is dead
  //proportionally to stem PLC (Ruffault et al. 2023)
  //Correct loading for phenology
  std::vector<double> cohLoading(numCohorts);
  std::vector<double> canopyFMC(numCohorts);
  for(int i=0;i<numCohorts;i++){
    canopyFMC[i] = (LFMC[i]*(1.0 - PLC[i]) + fm_dead*PLC[i]);
    if(x.above.LAI_live[i]>0.0) cohLoading[i] = x.above.Loading[i]*(x.above.LAI_expanded[i]/x.above.LAI_live[i]);
    else cohLoading[i] = 0.0;
  }
  
  //Average canopy moisture in the crown and surface layers
  if(FCCSprops.w[0] > 0.0) FCCSprops.ActFMC[0] = layerFuelAverageParameter_c(200.0, 10000.0, canopyFMC, cohLoading, x.above.H, x.above.CR);
  if(FCCSprops.w[1] > 0.0) FCCSprops.ActFMC[1] = layerFuelAverageParameter_c(0.0, 200.0, canopyFMC, cohLoading, x.above.H, x.above.CR);
  
  std::vector<double> MdeadSI = {fm_dead, fm_dead, fm_dead, fm_dead, fm_dead}; 
  std::vector<double> MliveSI = {90.0, 90.0, 60.0}; //Default values (not actually used if ActFMC is non-missing)
  
  FCCSbehaviour_c(FCCSbehres,
                  FCCSprops,
                  MliveSI, 
                  MdeadSI, 
                  slope, wind);
  
  //Copy results to fireBehavior vector
  FCCSres.Loading_overstory = FCCSprops.w[0];
  FCCSres.Loading_understory =  FCCSprops.w[1];
  FCCSres.CFMC_overstory =  FCCSprops.ActFMC[0];
  FCCSres.CFMC_understory = FCCSprops.ActFMC[1];
  FCCSres.DFMC = fm_dead;
  FCCSres.ROS_surface = FCCSbehres.surface.ROS;
  FCCSres.I_b_surface = FCCSbehres.surface.I_b;
  FCCSres.t_r_surface = FCCSbehres.surface.t_r;
  FCCSres.FL_surface = FCCSbehres.surface.FL;
  FCCSres.Ic_ratio = FCCSbehres.crown.Ic_ratio;
  FCCSres.ROS_crown = FCCSbehres.crown.ROS_crown;
  FCCSres.I_b_crown = FCCSbehres.crown.I_b_crown;
  FCCSres.t_r_crown = FCCSbehres.crown.t_r_crown;
  FCCSres.FL_crown = FCCSbehres.crown.FL_crown;
  FCCSres.SFP = FCCSbehres.potentials.SFP;
  FCCSres.CFP = FCCSbehres.potentials.CFP;
}


Rcpp::List copyFCCSBehaviour_Result_c(const FCCSBehaviour_RESULT& fccs) {
  List firePotentials=List::create(_["RP"] = fccs.potentials.RP,
                                   _["SP"] = fccs.potentials.SP,
                                   _["FP"] = fccs.potentials.FP,
                                   _["SFP"] = fccs.potentials.SFP,
                                   _["IC"] = fccs.potentials.IC,
                                   _["TC"] = fccs.potentials.TC,
                                   _["RC"] = fccs.potentials.RC,
                                   _["CFP"] = fccs.potentials.CFP);
  
  List crownFire=List::create(_["I_R_canopy [kJ/m2/min]"]= fccs.crown.I_R_canopy,
                              _["I_R_crown [kJ/m2/min]"]=fccs.crown.I_R_crown,
                              _["q_canopy [kJ/m2]"]=fccs.crown.q_canopy,
                              _["q_crown [kJ/m2]"]=fccs.crown.q_crown,
                              _["xi_crown"]=fccs.crown.xi_crown,
                              _["canopy_WindSpeed [m/s]"] = fccs.crown.canopy_WindSpeed,
                              _["WAF"]=fccs.crown.WAF,
                              _["ROS_crown [m/min]"]=fccs.crown.ROS_crown,
                              _["I_b_crown [kW/m]"] = fccs.crown.I_b_crown,
                              _["t_r_crown [s]"] = fccs.crown.t_r_crown,
                              _["Ic_ratio"] = fccs.crown.Ic_ratio,
                              _["FL_crown [m]"]=fccs.crown.FL_crown);
  
  List surfaceFire=List::create(_["midflame_WindSpeed [m/s]"] = fccs.surface.midflame_WindSpeed,
                                _["phi_wind"] = fccs.surface.phi_wind,
                                _["phi_slope"] = fccs.surface.phi_slope,
                                _["I_R_surf [kJ/m2/min]"]=fccs.surface.I_R_surf,
                                _["I_R_litter [kJ/m2/min]"]=fccs.surface.I_R_litter,
                                _["q_surf [kJ/m2]"]=fccs.surface.q_surf,
                                _["q_litter [kJ/m2]"]=fccs.surface.q_litter,
                                _["xi_surf"]=fccs.surface.xi_surf,
                                _["xi_litter"]=fccs.surface.xi_litter,
                                _["ROS_surf [m/min]"]=fccs.surface.ROS_surf,
                                _["ROS_litter [m/min]"]=fccs.surface.ROS_litter,
                                _["ROS_windslopecap [m/min]"]=fccs.surface.ROS_windslopecap,
                                _["ROS [m/min]"]=fccs.surface.ROS,
                                _["I_b [kW/m]"] = fccs.surface.I_b,
                                _["t_r [s]"] = fccs.surface.t_r, 
                                _["FL [m]"]=fccs.surface.FL);
  return(List::create(_["SurfaceFire"] = surfaceFire, _["CrownFire"] = crownFire, _["FirePotentials"] = firePotentials));
}

