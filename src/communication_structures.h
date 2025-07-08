#include <Rcpp.h>

#ifndef COMMUNICATION_STRUCTURES_H
#define COMMUNICATION_STRUCTURES_H
#endif
using namespace Rcpp;



List instanceCommunicationStructures(List x, String model);
List generalCommunicationStructures(int numCohorts, int nlayers, int ncanlayers, int ntimesteps,
                                    String model);

const int SOILWBCOM_dZ_m = 0;
const int SOILWBCOM_dZUp = 1;
const int SOILWBCOM_dZDown = 2;
const int SOILWBCOM_lambda = 3;
const int SOILWBCOM_theta_micro = 4;
const int SOILWBCOM_theta_b = 5;
const int SOILWBCOM_theta_macro = 6;
const int SOILWBCOM_theta_sat_fict = 7;
const int SOILWBCOM_Ksat_b = 8;
const int SOILWBCOM_Ksat_b_ms = 9;
const int SOILWBCOM_Ksat = 10;
const int SOILWBCOM_Ksat_ms = 11;
const int SOILWBCOM_Psi = 12;
const int SOILWBCOM_K = 13;
const int SOILWBCOM_C = 14;
const int SOILWBCOM_Psi_m = 15;
const int SOILWBCOM_K_ms = 16;
const int SOILWBCOM_Kbc = 17;
const int SOILWBCOM_Kbc_ms = 18;
const int SOILWBCOM_C_m = 19;
const int SOILWBCOM_S_macro = 20;
const int SOILWBCOM_e_macro = 21;
const int SOILWBCOM_Kmacro_ms = 22;
const int SOILWBCOM_waterFluidity = 23;
const int SOILWBCOM_a = 24;
const int SOILWBCOM_b = 25;
const int SOILWBCOM_c = 26;
const int SOILWBCOM_d = 27;
const int SOILWBCOM_e = 28;
const int SOILWBCOM_f = 29;
const int SOILWBCOM_K_step_ms05 = 30;
const int SOILWBCOM_C_step_m05 = 31;
const int SOILWBCOM_C_step = 32;
const int SOILWBCOM_C_step_m = 33;
const int SOILWBCOM_K_step_ms = 34;
const int SOILWBCOM_K_step = 35;
const int SOILWBCOM_Psi_step = 36;
const int SOILWBCOM_Psi_step_m = 37;
const int SOILWBCOM_S_macro_step = 38;
const int SOILWBCOM_Kmacro_step_ms = 39;
const int SOILWBCOM_theta_macro_step = 40;
const int SOILWBCOM_theta_micro_step = 41;
const int SOILWBCOM_finalSourceSinks_m3s = 42;
const int SOILWBCOM_capill_below = 43;
const int SOILWBCOM_drain_above = 44;
const int SOILWBCOM_drain_below = 45;
const int SOILWBCOM_lateral_flows_step_mm = 46;
List communicationSoilWaterBalance(int nlayers);

const int SOILEBCOM_dZ_m = 0;
const int SOILEBCOM_dZUp = 1;
const int SOILEBCOM_dZDown = 2;
const int SOILEBCOM_Zup = 3;
const int SOILEBCOM_Zdown = 4;
const int SOILEBCOM_Zcent = 5;
const int SOILEBCOM_a = 6;
const int SOILEBCOM_b = 7;
const int SOILEBCOM_c = 8;
const int SOILEBCOM_d = 9;
const int SOILEBCOM_e = 10;
const int SOILEBCOM_f = 11;
const int SOILEBCOM_k_up = 12;
const int SOILEBCOM_k_down = 13;
const int SOILEBCOM_tempch = 14;
List communicationSoilEnergyBalance(int nlayers);

List basicTranspirationCommunicationOutput(int numCohorts, int nlayers);
List advancedTranspirationCommunicationOutput(int numCohorts, int nlayers, int ncanlayers, int ntimesteps);
List copyBasicTranspirationOutput(List btc, List x);
List copyAdvancedTranspirationOutput(List atc, List x);
List copyBasicSPWBOutput(List boc, List x);
List copyAdvancedSPWBOutput(List aoc, List x);
List copyBasicGROWTHOutput(List boc, List x);
List copyAdvancedGROWTHOutput(List aoc, List x);
List copyModelOutput(List internalCommunication, List x, String model);
