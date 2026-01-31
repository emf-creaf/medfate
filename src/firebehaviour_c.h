#include <RcppArmadillo.h>
#include "modelInput_c.h"
#include "communication_structures_c.h"

#ifndef FIREBEHAVIOUR_C_H
#define FIREBEHAVIOUR_C_H

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

struct FCCS_RESULT {
  double Loading_overstory, Loading_understory, CFMC_overstory, CFMC_understory;
  double DFMC, ROS_surface, I_b_surface, t_r_surface;
  double FL_surface, Ic_ratio, ROS_crown, I_b_crown;
  double t_r_crown, FL_crown, SFP, CFP;
};
Rcpp::NumericVector copyFCCSResult_c(const FCCS_RESULT& fccs);

struct FCCSFirePotentials {
  double RP, SP, FP, SFP, IC, TC, RC, CFP;
};
struct FCCSCrownFire {
  double I_R_canopy, I_R_crown, q_canopy, q_crown, xi_crown;
  double canopy_WindSpeed, WAF, ROS_crown, I_b_crown, t_r_crown, Ic_ratio, FL_crown;
};
struct FCCSSurfaceFire {
  double midflame_WindSpeed, phi_wind, phi_slope;
  double I_R_surf, I_R_litter, q_surf, q_litter;
  double xi_surf, xi_litter;
  double ROS_surf, ROS_litter, ROS_windslopecap, ROS;
  double I_b, t_r, FL;
};
struct FCCSBehaviour_RESULT {
  FCCSSurfaceFire surface;
  FCCSCrownFire crown;
  FCCSFirePotentials potentials;
};
Rcpp::List copyFCCSBehaviour_Result_c(const FCCSBehaviour_RESULT& fccs);

double criticalFirelineIntensity_c(double CBH, double M);

void FCCSbehaviour_c(FCCSBehaviour_RESULT& res,
                     const InternalFCCS& FCCSpropsSI,
                     std::vector<double>& MliveSI, 
                     std::vector<double>& MdeadSI, 
                     double slope, double windSpeedSI);

void fccsHazard_c(FCCSBehaviour_RESULT& FCCSbehres, FCCS_RESULT& FCCSres, ModelInput& x, 
                  const WeatherInputVector& meteovec, 
                  const std::vector<double>& LFMC, 
                  const std::vector<double>& PLC, 
                  const double slope);
#endif
