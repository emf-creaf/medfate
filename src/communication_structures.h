#include <Rcpp.h>

#ifndef COMMUNICATION_STRUCTURES_H
#define COMMUNICATION_STRUCTURES_H
#endif
using namespace Rcpp;

const int SOILCOM_dZ_m = 0;
const int SOILCOM_dZUp = 1;
const int SOILCOM_dZDown = 2;
const int SOILCOM_lambda = 3;
const int SOILCOM_theta_micro = 4;
const int SOILCOM_theta_b = 5;
const int SOILCOM_theta_macro = 6;
const int SOILCOM_theta_sat_fict = 7;
const int SOILCOM_Ksat_b = 8;
const int SOILCOM_Ksat_b_ms = 9;
const int SOILCOM_Ksat = 10;
const int SOILCOM_Ksat_ms = 11;
const int SOILCOM_Psi = 12;
const int SOILCOM_K = 13;
const int SOILCOM_C = 14;
const int SOILCOM_Psi_m = 15;
const int SOILCOM_K_ms = 16;
const int SOILCOM_Kbc = 17;
const int SOILCOM_Kbc_ms = 18;
const int SOILCOM_C_m = 19;
const int SOILCOM_S_macro = 20;
const int SOILCOM_e_macro = 21;
const int SOILCOM_Kmacro_ms = 22;
const int SOILCOM_waterFluidity = 23;
const int SOILCOM_a = 24;
const int SOILCOM_b = 25;
const int SOILCOM_c = 26;
const int SOILCOM_d = 27;
const int SOILCOM_e = 28;
const int SOILCOM_f = 29;
const int SOILCOM_K_step_ms05 = 30;
const int SOILCOM_C_step_m05 = 31;
const int SOILCOM_C_step = 32;
const int SOILCOM_C_step_m = 33;
const int SOILCOM_K_step_ms = 34;
const int SOILCOM_K_step = 35;
const int SOILCOM_Psi_step = 36;
const int SOILCOM_Psi_step_m = 37;
const int SOILCOM_S_macro_step = 38;
const int SOILCOM_Kmacro_step_ms = 39;
const int SOILCOM_theta_macro_step = 40;
const int SOILCOM_theta_micro_step = 41;
const int SOILCOM_finalSourceSinks_m3s = 42;
const int SOILCOM_capill_below = 43;
const int SOILCOM_drain_above = 44;
const int SOILCOM_drain_below = 45;
const int SOILCOM_lateral_flows_step_mm = 46;


List instanceCommunicationStructures(List x, String model);
List generalCommunicationStructures(int numCohorts, int nlayers, int ncanlayers, int ntimesteps,
                                    String model);
List communicationSoilWaterBalance(int nlayers);

List basicTranspirationCommunicationOutput(int numCohorts, int nlayers);
List advancedTranspirationCommunicationOutput(int numCohorts, int nlayers, int ncanlayers, int ntimesteps);
List copyBasicTranspirationOutput(List btc, List x);
List copyAdvancedTranspirationOutput(List atc, List x);
List copyBasicSPWBOutput(List boc, List x);
List copyAdvancedSPWBOutput(List aoc, List x);
List copyBasicGROWTHOutput(List boc, List x);
List copyAdvancedGROWTHOutput(List aoc, List x);
List copyModelOutput(List internalCommunication, List x, String model);
