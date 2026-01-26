#include "medfate.h"
#include "Rcpp.h"

#ifndef SOIL_C_H
#define SOIL_C_H

//*
//* Clapp-Hornberger parameters structure
//* 
struct ClappHornberger{
  double theta_sat;
  double psi_sat_cm;
  double b;
  double K_sat_cm_h;
  ClappHornberger() : theta_sat(0.0), psi_sat_cm(0.0), b(0.0), K_sat_cm_h(0.0) {}
  ClappHornberger(std::string soilType) {
    if(soilType=="Sand") {theta_sat=0.395; psi_sat_cm=-12.1; b = 4.05; K_sat_cm_h=63.36;}
    else if(soilType=="Loamy sand") {theta_sat=0.410; psi_sat_cm=-9.1;b = 4.38; K_sat_cm_h=56.28;}
    else if(soilType=="Sandy loam") {theta_sat=0.435; psi_sat_cm=-21.8; b = 4.90; K_sat_cm_h=12.48;}
    else if(soilType=="Silt loam") {theta_sat=0.485; psi_sat_cm=-78.6; b = 5.30; K_sat_cm_h=2.59;}
    else if(soilType=="Loam") {theta_sat=0.451; psi_sat_cm=-47.8; b = 5.39; K_sat_cm_h=2.50;}
    else if(soilType=="Silt") {theta_sat=0.485; psi_sat_cm=-78.6; b = 5.30; K_sat_cm_h=2.59;} // EQUAL TO SILT LOAM
    else if(soilType=="Sandy clay loam") {theta_sat=0.420; psi_sat_cm=-29.9; b = 7.12; K_sat_cm_h=2.27;}
    else if(soilType=="Silty clay loam") {theta_sat=0.477; psi_sat_cm=-35.6; b = 7.75; K_sat_cm_h=0.61;}
    else if(soilType=="Clay loam") {theta_sat=0.476; psi_sat_cm=-63.0; b = 8.52; K_sat_cm_h=0.88;}
    else if(soilType=="Sandy clay") {theta_sat=0.426; psi_sat_cm=-15.3; b = 10.4; K_sat_cm_h=0.38;}
    else if(soilType=="Silty clay") {theta_sat=0.492; psi_sat_cm=-49.0; b = 10.4; K_sat_cm_h=0.37;}
    else if(soilType=="Clay") {theta_sat=0.482; psi_sat_cm=-40.5; b = 11.4; K_sat_cm_h=0.46;}
  }
};

//*
//* Soil structure to hold soil properties and state variables
//* 
class Soil {
  private:
  int nlayers;
  std::string model;
  ClappHornberger clapp_hornberger;
  std::vector<double> widths;
  std::vector<double> clay;
  std::vector<double> sand;
  std::vector<double> om;
  std::vector<double> nitrogen;
  std::vector<double> ph;
  std::vector<double> bd;
  std::vector<double> rfc;
  std::vector<double> macro;
  std::vector<double> Ksat;
  std::vector<double> VG_alpha;
  std::vector<double> VG_n;
  std::vector<double> VG_theta_res;
  std::vector<double> VG_theta_sat;
  std::vector<std::string> usda_type;
  std::vector<double> theta_SAT;
  std::vector<double> theta_FC;
  std::vector<double> W;
  std::vector<double> psi;
  std::vector<double> theta;
  std::vector<double> Temp;
  public:
    Soil();
    Soil(int nlayersIn, 
             std::string& modelIn,
             std::vector<double>& widthsIn,
             std::vector<double>& clayIn,
             std::vector<double>& sandIn,
             std::vector<double>& omIn,
             std::vector<double>& nitrogenIn,
             std::vector<double>& phIn,
             std::vector<double>& bdIn,
             std::vector<double>& rfcIn,
             std::vector<double>& macroIn,
             std::vector<double>& KsatIn,
             std::vector<double>& VG_alphaIn,
             std::vector<double>& VG_nIn,
             std::vector<double>& VG_theta_resIn,
             std::vector<double>& VG_theta_satIn,
             std::vector<std::string>& usda_typeIn,
             std::vector<double>& theta_SATIn,
             std::vector<double>& theta_FCIn,
             std::vector<double>& WIn,
             std::vector<double>& psiIn,
             std::vector<double>& thetaIn,
             std::vector<double>& TempIn,
             ClappHornberger& clapp_hornbergerIn);
    Soil(Rcpp::DataFrame x, Rcpp::String model);
    
    double getW(int layer);
    int getNlayers();
    std::string getModel();
    ClappHornberger getClappHornberger();
    double getWidth(int layer);  
    double getClay(int layer);
    double getSand(int layer);
    double getOM(int layer);
    double getNitrogen(int layer);
    double getPH(int layer);
    double getBD(int layer);
    double getRFC(int layer);
    double getMacro(int layer);
    double getKsat(int layer);
    double getVG_alpha(int layer);
    double getVG_n(int layer);
    double getVG_theta_res(int layer);
    double getVG_theta_sat(int layer);
    std::string getUSDAType(int layer);
    double getThetaSAT(int layer);
    double getThetaFC(int layer);
    double getPsi(int layer);
    double getTheta(int layer);
    double getTemp(int layer);
    double getWaterSAT(int layer);
    double getWaterFC(int layer);
    void setPsi(int layer, double value);
    void setTheta(int layer, double value);
    void setW(int layer, double value);
    void setTemp(int layer, double value);
    
};

/**
 * Conversion factor from conductivity in cm·day-1 to molH20·m-2·MPa-1·s-1
 *  1 day = 86400 sec
 *  1 mol H20 = 18.01528 g
 */
const double cmdTOmmolm2sMPa = 655.2934;//100.0/(18.01528*86400.0*0.00009804139432); 
/**
 * Conversion factor from cm to MPa
 */
const double cmTOMPa = 0.00009804139432;
/**
 * Conversion factor from m to MPa
 */
const double mTOMPa = 0.009804139432;//1/9.804139*0.000001; 

const double minimumPsi = -40.0;
const double fieldCapacityPsi = -0.033; // -33 kPa

std::string USDAType_c(double clay, double sand);

/**
 * Saxton et al. (1986) pedotransfer functions
 */
double saturatedConductivitySaxton_c(double clay, double sand, double bd, double om, bool mmol);
double unsaturatedConductivitySaxton_c(double theta, double clay, double sand, double bd, double om, bool mmol);
double theta2psiSaxton_c(double clay, double sand, double theta, double om);
double psi2thetaSaxton_c(double clay, double sand, double psi, double om);
double thetaSATSaxton_c(double clay, double sand, double om);

/**
 * Van Genuchten-Mualem pedotransfer functions
 */
double theta2psiVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double theta);
double psi2thetaVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchten_c(double ksat, double n, double alpha, double theta_res, double theta_sat, double psi);
double psi2kVanGenuchtenMicropores_c(double k_b, double n, double alpha, double theta_res, double theta_sat, 
                                   double psi, double psi_b);
double psi2DVanGenuchten_c(double k_sat, double n, double alpha, double theta_res, double theta_sat, 
                         double psi);
double psi2cVanGenuchten_c(double n, double alpha, double theta_res, double theta_sat, double psi);
#endif
