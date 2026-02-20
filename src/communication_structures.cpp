// [[Rcpp::interfaces(r,cpp)]]
#include <RcppArmadillo.h>
#include "communication_structures.h"
using namespace Rcpp;

List communicationLongWaveRadiation(int ncanlayers) {
  NumericVector Lup(ncanlayers, NA_REAL), Ldown(ncanlayers, NA_REAL), Lnet(ncanlayers, NA_REAL);
  NumericVector tau(ncanlayers, NA_REAL), sumTauComp(ncanlayers, NA_REAL);
  DataFrame LWR_layer = DataFrame::create(_["tau"] = tau,
                                          _["sumTauComp"] = sumTauComp,
                                          _["Ldown"] = Ldown, 
                                          _["Lup"] = Lup,
                                          _["Lnet"] = Lnet);
  List lwr_struct = List::create(_["LWR_layer"] = LWR_layer,
                                 _["Ldown_ground"] = NA_REAL,
                                 _["Lup_ground"] = NA_REAL,
                                 _["Lnet_ground"] = NA_REAL,
                                 _["Ldown_canopy"] = NA_REAL,
                                 _["Lup_canopy"] = NA_REAL,
                                 _["Lnet_canopy"] = NA_REAL,
                                 _["Lnet_cohort_layer"] = NA_REAL);
  return(lwr_struct);
}

DataFrame communicationCanopyTurbulence(int ncanlayers) {
  DataFrame output = DataFrame::create(Named("zmid") = NumericVector(ncanlayers, NA_REAL),
                                       Named("u") = NumericVector(ncanlayers, NA_REAL),
                                       Named("du") = NumericVector(ncanlayers, NA_REAL),
                                       Named("epsilon") = NumericVector(ncanlayers, NA_REAL),
                                       Named("k") = NumericVector(ncanlayers, NA_REAL),
                                       Named("uw") = NumericVector(ncanlayers, NA_REAL));
  return(output);
}


List communicationSoilWaterBalance(int nlayers) {
  int ncol = 50;
  List out(ncol);
  CharacterVector colnames(ncol);
  for(int i = 0; i<ncol; i++) out[i] = NumericVector(nlayers, NA_REAL);
  colnames[SOILWBCOM_dZ_m] = "dZ_m";
  colnames[SOILWBCOM_dZUp] = "dZUp";
  colnames[SOILWBCOM_dZDown] = "dZDown";
  colnames[SOILWBCOM_lambda] = "lambda";
  colnames[SOILWBCOM_theta_micro] = "theta_micro";
  colnames[SOILWBCOM_theta_b] = "theta_b";
  colnames[SOILWBCOM_theta_macro] = "theta_macro";
  colnames[SOILWBCOM_theta_sat_fict] = "theta_sat_fict";
  colnames[SOILWBCOM_Ksat_b] = "Ksat_b";
  colnames[SOILWBCOM_Ksat_b_ms] = "Ksat_b_ms";
  colnames[SOILWBCOM_Ksat] = "Ksat";
  colnames[SOILWBCOM_Ksat_ms] = "Ksat_ms";
  colnames[SOILWBCOM_Psi] = "Psi";
  colnames[SOILWBCOM_K] = "K";
  colnames[SOILWBCOM_C] = "C";
  colnames[SOILWBCOM_Psi_m] = "Psi_m";
  colnames[SOILWBCOM_K_ms] = "K_ms";
  colnames[SOILWBCOM_Kbc] = "Kbc";
  colnames[SOILWBCOM_Kbc_ms] = "Kbc_ms";
  colnames[SOILWBCOM_C_m] = "C_m";
  colnames[SOILWBCOM_S_macro] = "S_macro";
  colnames[SOILWBCOM_e_macro] = "e_macro";
  colnames[SOILWBCOM_Kmacro_ms] = "Kmacro_ms";
  colnames[SOILWBCOM_waterFluidity] = "waterFluidity";
  colnames[SOILWBCOM_a] = "a";
  colnames[SOILWBCOM_b] = "b";
  colnames[SOILWBCOM_c] = "c";
  colnames[SOILWBCOM_d] = "d";
  colnames[SOILWBCOM_e] = "e";
  colnames[SOILWBCOM_f] = "f";
  colnames[SOILWBCOM_K_step_ms05] = "K_step_ms05";
  colnames[SOILWBCOM_C_step_m05] = "C_step_m05";
  colnames[SOILWBCOM_C_step] = "C_step";
  colnames[SOILWBCOM_C_step_m] = "C_step_m";
  colnames[SOILWBCOM_K_step_ms] = "K_step_ms";
  colnames[SOILWBCOM_K_step] = "K_step";
  colnames[SOILWBCOM_Psi_step] = "Psi_step";
  colnames[SOILWBCOM_Psi_step_m] = "Psi_step_m";
  colnames[SOILWBCOM_S_macro_step] = "S_macro_step";
  colnames[SOILWBCOM_Kmacro_step_ms] = "Kmacro_step_ms";
  colnames[SOILWBCOM_theta_macro_step] = "theta_macro_step";
  colnames[SOILWBCOM_theta_micro_step] = "theta_micro_step";
  colnames[SOILWBCOM_finalSourceSinks_m3s] = "finalSourceSinks_m3s";
  colnames[SOILWBCOM_capill_below] = "capill_below";
  colnames[SOILWBCOM_drain_above] = "drain_above";
  colnames[SOILWBCOM_drain_below] = "drain_below";
  colnames[SOILWBCOM_lateral_flows_step_mm] = "lateral_flows_step_mm";
  out.attr("names") = colnames;
  return(out);
}

List communicationSoilEnergyBalance(int nlayers) {
  int ncol = 15;
  List out(ncol);
  CharacterVector colnames(ncol);
  for(int i = 0; i<ncol; i++) out[i] = NumericVector(nlayers, NA_REAL);
  colnames[SOILEBCOM_dZ_m] = "dZ_m";
  colnames[SOILEBCOM_dZUp] = "dZUp";
  colnames[SOILEBCOM_dZDown] = "dZDown";
  colnames[SOILEBCOM_Zup] = "Zup";
  colnames[SOILEBCOM_Zdown] = "Zdown";
  colnames[SOILEBCOM_Zcent] = "Zcent";
  colnames[SOILEBCOM_a] = "a";
  colnames[SOILEBCOM_b] = "b";
  colnames[SOILEBCOM_c] = "c";
  colnames[SOILEBCOM_d] = "d";
  colnames[SOILEBCOM_e] = "e";
  colnames[SOILEBCOM_f] = "f";
  colnames[SOILEBCOM_k_up] = "k_up";
  colnames[SOILEBCOM_k_down] = "k_down";
  colnames[SOILEBCOM_tempch] = "tempch";
  out.attr("names") = colnames;
  return(out);
}
NumericVector communicationFireHazard() {
  NumericVector fireHazard = NumericVector::create(
    _["Loading_overstory [kg/m2]"] = NA_REAL,
    _["Loading_understory [kg/m2]"] = NA_REAL,
    _["CFMC_overstory [%]"] = NA_REAL,
    _["CFMC_understory [%]"] = NA_REAL,
    _["DFMC [%]"] = NA_REAL,
    _["ROS_surface [m/min]"] = NA_REAL,
    _["I_b_surface [kW/m]"] = NA_REAL,
    _["t_r_surface [s]"] = NA_REAL,
    _["FL_surface [m]"] = NA_REAL,
    _["Ic_ratio"] = NA_REAL,
    _["ROS_crown [m/min]"] = NA_REAL,
    _["I_b_crown [kW/m]"] = NA_REAL,
    _["t_r_crown [s]"] = NA_REAL,
    _["FL_crown [m]"] = NA_REAL,
    _["SFP"] = NA_REAL,
    _["CFP"] = NA_REAL);
  return(fireHazard);
}


NumericVector communicationSnagDecomposition() {
  NumericVector output = NumericVector::create(_["transfer_surface_active"] = 0.0,
                                               _["transfer_surface_slow"] = 0.0,
                                               _["flux_respiration"] = 0.0);
  return(output);
}


NumericVector communicationLitterDecomposition() {
  NumericVector output = NumericVector::create(_["transfer_surface_active"] = 0.0,
                                               _["transfer_surface_slow"] = 0.0,
                                               _["transfer_soil_active"] = 0.0,
                                               _["transfer_soil_slow"] = 0.0,
                                               _["flux_respiration"] = 0.0);
  return(output);
}


// Creates list with the following matrices:
// pools:
//   \itemize{
//     \item{\code{xi}: Environmental scalar matrix.}
//     \item{\code{A}: Carbon transfer matrix.} 
//     \item{\code{pathf}: Fractional carbon flow from pool j to pool i.} 
//     \item{\code{respf}: Fractional respiration loss for carbon flow from pool j to pool i.} 
// }
List communicationDecomposition() {
  int npool = 7;
  NumericVector snagDecompositionOutput = communicationSnagDecomposition();
  NumericVector litterDecompositionOutput = communicationLitterDecomposition();
  NumericVector xi(npool), K(npool);
  NumericMatrix A(npool, npool);
  NumericMatrix pathf(npool, npool);
  NumericMatrix respf(npool, npool);
  for(int i = 0; i< npool; i++) {
    xi[i] = 0.0;
    K[i] = 0.0;
    for(int j = 0; j< npool; j++) {
      A(i,j) = 0.0;
      pathf(i,j) = 0.0;
      respf(i,j) = 0.0;
    }
  }
  List l = List::create(_["sdo"] = snagDecompositionOutput,
                        _["ldo"] = litterDecompositionOutput,
                        _["xi"] = xi,
                        _["K"] = K,
                        _["Kmix"] = 0.0,
                        _["K_s21"] = 0.0,
                        _["A"] = A,
                        _["pathf"] = pathf,
                        _["respf"] = respf);
  return(l);
}

DataFrame communicationCarbonCompartments(int numCohorts) {
  DataFrame df = DataFrame::create(
    _["LeafStorageVolume"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStorageVolume"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStarchMaximumConcentration"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStarchMaximumConcentration"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStarchCapacity"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStarchCapacity"] = NumericVector(numCohorts, NA_REAL),
    _["LeafStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TwigStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TwigLivingStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["DeadLeafStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["DeadTwigStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["SapwoodLivingStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["HeartwoodStructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["AbovegroundWoodBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["BelowgroundWoodBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["FineRootBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["LabileBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["StructuralBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["DeadBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TotalLivingBiomass"] = NumericVector(numCohorts, NA_REAL),
    _["TotalBiomass"] = NumericVector(numCohorts, NA_REAL)
  );
  return(df);
}

DataFrame copyDataFrame(DataFrame comm, int numRows) {
  CharacterVector colnames = comm.attr("names");
  int n = colnames.size();
  List out(n); 
  for(int i = 0;i<n;i++) {
    String nameCol = colnames[i];
    NumericVector vecIn = comm[nameCol];
    NumericVector vecOut(numRows);
    for(int c = 0;c<numRows;c++) {
      vecOut[c] = vecIn[c];
    }
    out[i] = vecOut;
  }
  out.attr("names") = clone(colnames);
  DataFrame dfout(out);
  return(dfout);
}
NumericMatrix copyNumericMatrix(NumericMatrix comm, int rows, int cols) {
  NumericMatrix out(rows, cols);
  for(int r=0;r<rows;r++) {
    for(int c=0;c<cols; c++) {
      out(r, c) = comm(r, c);
    }
  }
  return(out);
}
List copyList(List comm, int n) {
  List out(n);
  for(int i=0;i<n;i++) {
    out[i] = clone(as<List>(comm[i]));
  }
  return(out);
}
List copyCommunicationLongWaveRadiation(List clwr, int ncanlayers) {
  List lwr_struct = List::create(_["LWR_layer"] = copyDataFrame(as<DataFrame>(clwr["LWR_layer"]), ncanlayers),
                                 _["Ldown_ground"] = clwr["Ldown_ground"],
                                 _["Lup_ground"] = clwr["Lup_ground"],
                                 _["Lnet_ground"] = clwr["Lnet_ground"],
                                 _["Ldown_canopy"] = clwr["Ldown_canopy"],
                                 _["Lup_canopy"] = clwr["Lup_canopy"],
                                 _["Lnet_canopy"] = clwr["Lnet_canopy"],
                                 _["Lnet_cohort_layer"] = clwr["Lnet_cohort_layer"]);
  return(lwr_struct);
}
List copyBasicTranspirationOutput(List btc, List x) {
  List control = x["control"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = above.nrow();
  
  NumericMatrix ExtractionComm = btc["Extraction"];
  NumericMatrix Extraction = copyNumericMatrix(ExtractionComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
  Extraction.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  List ExtractionPoolsComm = btc["ExtractionPools"];
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCohComm = ExtractionPoolsComm[c];
      NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = above.attr("row.names");
  }
  
  NumericVector StandComm = btc["Stand"];
  NumericVector Stand = clone(StandComm);
  
  DataFrame PlantsComm = Rcpp::as<Rcpp::DataFrame>(btc["Plants"]);
  DataFrame Plants = copyDataFrame(PlantsComm, numCohorts);
  Plants.attr("row.names") = above.attr("row.names");
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["Extraction"] = Extraction,
                        _["ExtractionPools"] = ExtractionPools);
  return(l);
}
List copyPlantsInstOutput(List PlantsInstComm, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  NumericMatrix Einst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["E"]), numCohorts, ntimesteps);
  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Aginst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["Ag"]), numCohorts, ntimesteps);
  Aginst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Aninst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["An"]), numCohorts, ntimesteps);
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix dEdPInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["dEdP"]), numCohorts, ntimesteps);
  dEdPInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafPsi"]), numCohorts, ntimesteps);
  LeafPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemPsi"]), numCohorts, ntimesteps);
  StemPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix RootPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["RootPsi"]), numCohorts, ntimesteps);
  RootPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafSympPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafSympPsi"]), numCohorts, ntimesteps);
  LeafSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemSympPsiInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemSympPsi"]), numCohorts, ntimesteps);
  StemSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafSympRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafSympRWC"]), numCohorts, ntimesteps);
  LeafSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemSympRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemSympRWC"]), numCohorts, ntimesteps);
  StemSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemPLC = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemPLC"]), numCohorts, ntimesteps);
  StemPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafPLC = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafPLC"]), numCohorts, ntimesteps);
  LeafPLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LeafRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["LeafRWC"]), numCohorts, ntimesteps);
  LeafRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix StemRWCInst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["StemRWC"]), numCohorts, ntimesteps);
  StemRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix PWBinst = copyNumericMatrix(as<NumericMatrix>(PlantsInstComm["PWB"]), numCohorts, ntimesteps);
  PWBinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  
  List PlantsInst = List::create(
    _["E"]=Einst, _["Ag"]=Aginst, _["An"]=Aninst,
    _["dEdP"] = dEdPInst,
    _["RootPsi"] = RootPsiInst, 
    _["StemPsi"] = StemPsiInst,
    _["LeafPsi"] = LeafPsiInst,
    _["StemSympPsi"] = StemSympPsiInst,
    _["LeafSympPsi"] = LeafSympPsiInst,
    _["StemPLC"] = StemPLC, 
    _["LeafPLC"] = LeafPLC, 
    _["StemRWC"] = StemRWCInst,
    _["LeafRWC"] = LeafRWCInst,
    _["StemSympRWC"] = StemSympRWCInst,
    _["LeafSympRWC"] = LeafSympRWCInst,
    _["PWB"] = PWBinst);
  return(PlantsInst);
}
List copyLeavesInstOutput(List LeavesInstComm, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  
  NumericMatrix LAI = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["LAI"]), numCohorts, ntimesteps);
  LAI.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Vmax298 = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Vmax298"]), numCohorts, ntimesteps);
  Vmax298.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Jmax298 = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Jmax298"]), numCohorts, ntimesteps);
  Jmax298.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix SWR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Abs_SWR"]), numCohorts, ntimesteps);
  SWR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix PAR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Abs_PAR"]), numCohorts, ntimesteps);
  PAR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix LWR = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Net_LWR"]), numCohorts, ntimesteps);
  LWR.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix An = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["An"]), numCohorts, ntimesteps);
  An.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Ag = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Ag"]), numCohorts, ntimesteps);
  Ag.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Ci = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Ci"]), numCohorts, ntimesteps);
  Ci.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix E = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["E"]), numCohorts, ntimesteps);
  E.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix GSW = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Gsw"]), numCohorts, ntimesteps);
  GSW.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix VPD = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["VPD"]), numCohorts, ntimesteps);
  VPD.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Temp = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Temp"]), numCohorts, ntimesteps);
  Temp.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  NumericMatrix Psi = copyNumericMatrix(as<NumericMatrix>(LeavesInstComm["Psi"]), numCohorts, ntimesteps);
  Psi.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));

  
  List LeavesInst = List::create(
    _["LAI"] = LAI,
    _["Vmax298"] = Vmax298,
    _["Jmax298"] = Jmax298,
    _["Abs_SWR"]=SWR,
    _["Abs_PAR"]=PAR,
    _["Net_LWR"] = LWR,
    _["Ag"] = Ag,
    _["An"] = An,
    _["Ci"] = Ci,
    _["E"] = E,
    _["Gsw"] = GSW,
    _["VPD"] = VPD,
    _["Temp"] = Temp,
    _["Psi"] = Psi);
  
  return(LeavesInst);
}
List copyEnergyBalanceOutput(List EnergyBalanceComm, List x) {
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  
  DataFrame Tinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["Temperature"]), ntimesteps);
  DataFrame CEBinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["CanopyEnergyBalance"]), ntimesteps);
  DataFrame SEBinst = copyDataFrame(as<DataFrame>(EnergyBalanceComm["SoilEnergyBalance"]), ntimesteps);
  NumericMatrix Tcan_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["TemperatureLayers"]), ntimesteps, ncanlayers);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix VPcan_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["VaporPressureLayers"]), ntimesteps, ncanlayers);
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  NumericMatrix Tsoil_mat= copyNumericMatrix(as<NumericMatrix>(EnergyBalanceComm["SoilTemperature"]), ntimesteps, nlayers);
  Tsoil_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,nlayers));
  List EnergyBalance = List::create(_["Temperature"]=Tinst, 
                                    _["SoilTemperature"] = Tsoil_mat,
                                    _["CanopyEnergyBalance"] = CEBinst, 
                                    _["SoilEnergyBalance"] = SEBinst,
                                    _["TemperatureLayers"] = Tcan_mat, 
                                    _["VaporPressureLayers"] = VPcan_mat);
  return(EnergyBalance);
}
List copyAdvancedTranspirationOutput(List atc, List x) {
  List control = x["control"];
  int ntimesteps = control["ndailysteps"];
  String transpirationMode = control["transpirationMode"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = above.nrow();
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  
  NumericMatrix SoilWaterExtractComm = atc["Extraction"];
  NumericMatrix SoilWaterExtract = copyNumericMatrix(SoilWaterExtractComm, numCohorts, nlayers);
  SoilWaterExtract.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  List ExtractionPools(numCohorts);
  List ExtractionPoolsComm = atc["ExtractionPools"];
  if(plantWaterPools) {
    for(int c=0;c<numCohorts;c++) {
      NumericMatrix ExtractionPoolsCohComm = ExtractionPoolsComm[c];
      NumericMatrix ExtractionPoolsCohComm_c = copyNumericMatrix(ExtractionPoolsCohComm, numCohorts, nlayers); // this is final extraction of each cohort from each layer
      ExtractionPoolsCohComm_c.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
      ExtractionPools[c] = ExtractionPoolsCohComm_c;
    }
    ExtractionPools.attr("names") = above.attr("row.names");
  }
  
  NumericMatrix soilLayerExtractInstComm = atc["ExtractionInst"];
  NumericMatrix soilLayerExtractInst = copyNumericMatrix(soilLayerExtractInstComm, nlayers, ntimesteps);
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  
  NumericMatrix minPsiRhizoComm = atc["RhizoPsi"];
  NumericMatrix minPsiRhizo = copyNumericMatrix(minPsiRhizoComm, numCohorts, nlayers);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  
  NumericVector StandComm = atc["Stand"];
  NumericVector Stand = clone(StandComm);
  
  List EnergyBalance = copyEnergyBalanceOutput(as<List>(atc["EnergyBalance"]), x);

  DataFrame Plants = copyDataFrame(as<DataFrame>(atc["Plants"]), numCohorts);
  Plants.attr("row.names") = above.attr("row.names");
  
  DataFrame Sunlit = copyDataFrame(as<DataFrame>(atc["SunlitLeaves"]), numCohorts);
  Sunlit.attr("row.names") = above.attr("row.names");
  DataFrame Shade = copyDataFrame(as<DataFrame>(atc["ShadeLeaves"]), numCohorts);
  Shade.attr("row.names") = above.attr("row.names");
  
  List PlantsInst = copyPlantsInstOutput(as<List>(atc["PlantsInst"]), x);
  
  List SunlitInst = copyLeavesInstOutput(as<List>(atc["SunlitLeavesInst"]), x);
  List ShadeInst = copyLeavesInstOutput(as<List>(atc["ShadeLeavesInst"]), x);
  
  List lwrExtinctionListComm = atc["LWRExtinction"];
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyCommunicationLongWaveRadiation(as<List>(lwrExtinctionListComm[n]), ncanlayers);
  }

  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["EnergyBalance"] = EnergyBalance,
                        _["Extraction"] = SoilWaterExtract,
                        _["ExtractionPools"] = ExtractionPools,
                        _["RhizoPsi"] = minPsiRhizo,
                        _["Stand"] = Stand,
                        _["Plants"] = Plants,
                        _["SunlitLeaves"] = Sunlit,
                        _["ShadeLeaves"] = Shade,
                        _["ExtractionInst"] = soilLayerExtractInst,
                        _["RadiationInputInst"] = clone(as<DataFrame>(atc["RadiationInputInst"])), 
                        _["PlantsInst"] = PlantsInst,
                        _["SunlitLeavesInst"] = SunlitInst,
                        _["ShadeLeavesInst"] = ShadeInst,
                        _["LightExtinction"] =  clone(as<List>(atc["LightExtinction"])), 
                        _["LWRExtinction"] = lwrExtinctionList,
                        _["CanopyTurbulence"] = copyDataFrame(as<DataFrame>(atc["CanopyTurbulence"]), ncanlayers)); //To be replaced
  
  List supply = copyList(as<List>(atc["SupplyFunctions"]), numCohorts);
  supply.attr("names") = above.attr("row.names");
  List outPhotoSunlit = copyList(as<List>(atc["PhotoSunlitFunctions"]), numCohorts);
  outPhotoSunlit.attr("names") = above.attr("row.names");
  List outPhotoShade = copyList(as<List>(atc["PhotoShadeFunctions"]), numCohorts);
  outPhotoShade.attr("names") = above.attr("row.names");
  List outPMSunlit = copyList(as<List>(atc["PMSunlitFunctions"]), numCohorts);
  outPMSunlit.attr("names") = above.attr("row.names");
  List outPMShade = copyList(as<List>(atc["PMShadeFunctions"]), numCohorts);
  outPMShade.attr("names") = above.attr("row.names");
  l.push_back(supply, "SupplyFunctions");
  l.push_back(outPhotoSunlit, "PhotoSunlitFunctions");
  l.push_back(outPhotoShade, "PhotoShadeFunctions");
  l.push_back(outPMSunlit, "PMSunlitFunctions");
  l.push_back(outPMShade, "PMShadeFunctions");
  
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}


List copyBasicSPWBOutput(List boc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  
  NumericVector topo = clone(as<NumericVector>(boc["topography"]));
  NumericVector meteovec_bas = clone(as<NumericVector>(boc["weather"]));
  NumericVector WaterBalance = clone(as<NumericVector>(boc["WaterBalance"]));
  
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_bas,
                        _["WaterBalance"] = WaterBalance);

  if(control["soilResults"]) {
    DataFrame Soil = copyDataFrame(as<DataFrame>(boc["Soil"]), nlayers);
    l.push_back(Soil, "Soil");
  }  
  if(control["standResults"]) {
    NumericVector Stand = clone(as<NumericVector>(boc["Stand"]));
    l.push_back(Stand, "Stand");
  }
  if(control["plantResults"]) {
    DataFrame PlantsComm = Rcpp::as<Rcpp::DataFrame>(boc["Plants"]);
    DataFrame Plants = copyDataFrame(PlantsComm, numCohorts);
    Plants.attr("row.names") = cohorts.attr("row.names");
    l.push_back(Plants, "Plants");
  }
  if(control["fireHazardResults"]) {
    l.push_back(clone(as<NumericVector>(boc["FireHazard"])), "FireHazard");
  }
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

List copyBasicGROWTHOutput(List boc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  
  List spwbOut = copyBasicSPWBOutput(boc, x);

  NumericVector standCB = clone(as<NumericVector>(boc["CarbonBalance"]));
    
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["CarbonBalance"] = standCB);
  if(control["soilResults"]) {
    l.push_back(spwbOut["Soil"], "Soil");
  }
  if(control["standResults"]) {
    l.push_back(spwbOut["Stand"], "Stand");
  }
  if(control["plantResults"]) {
    l.push_back(spwbOut["Plants"], "Plants");
  }
  if(control["labileCarbonBalanceResults"]) {
    DataFrame labileCarbonBalance = copyDataFrame(as<DataFrame>(boc["LabileCarbonBalance"]), numCohorts);
    labileCarbonBalance.attr("row.names") = above.attr("row.names");
    l.push_back(labileCarbonBalance, "LabileCarbonBalance");
  }
  DataFrame plantBiomassBalance = copyDataFrame(as<DataFrame>(boc["PlantBiomassBalance"]), numCohorts);
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  l.push_back(plantBiomassBalance, "PlantBiomassBalance");
  
  if(control["plantStructureResults"]) {
    DataFrame plantStructure = copyDataFrame(as<DataFrame>(boc["PlantStructure"]), numCohorts);
    plantStructure.attr("row.names") = above.attr("row.names");
    l.push_back(plantStructure, "PlantStructure");
  }
  
  if(control["growthMortalityResults"]) {
    DataFrame growthMortality  = copyDataFrame(as<DataFrame>(boc["GrowthMortality"]), numCohorts);
    growthMortality.attr("row.names") = above.attr("row.names");
    l.push_back(growthMortality, "GrowthMortality");
  }
  
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

List copyAdvancedSPWBOutput(List aoc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  int ncanlayers = canopyParams.nrow(); //Number of canopy layers
  int nlayers = soil.nrow();
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  
  NumericVector topo = clone(as<NumericVector>(aoc["topography"]));
  NumericVector meteovec_adv = clone(as<NumericVector>(aoc["weather"]));
  
  NumericVector WaterBalance = clone(as<NumericVector>(aoc["WaterBalance"]));
  List EnergyBalance = copyEnergyBalanceOutput(as<List>(aoc["EnergyBalance"]), x);
  List l = List::create(_["cohorts"] = clone(cohorts),
                        _["topography"] = topo,
                        _["weather"] = meteovec_adv,
                        _["WaterBalance"] = WaterBalance,
                        _["EnergyBalance"] = EnergyBalance);

  if(control["soilResults"]){
    DataFrame Soil = copyDataFrame(as<DataFrame>(aoc["Soil"]), nlayers);
    l.push_back(Soil, "Soil");
  }
  if(control["standResults"]){
    NumericVector Stand = clone(as<NumericVector>(aoc["Stand"]));
    l.push_back(Stand, "Stand");
  }
  if(control["plantResults"]){
    DataFrame Plants = copyDataFrame(as<DataFrame>(aoc["Plants"]), numCohorts);
    Plants.attr("row.names") = cohorts.attr("row.names");
    l.push_back(Plants, "Plants"); 
  }
  if(control["leafResults"]){
    DataFrame Sunlit = copyDataFrame(as<DataFrame>(aoc["SunlitLeaves"]), numCohorts);
    Sunlit.attr("row.names") = above.attr("row.names");
    DataFrame Shade = copyDataFrame(as<DataFrame>(aoc["ShadeLeaves"]), numCohorts);
    Shade.attr("row.names") = above.attr("row.names");
    l.push_back(Sunlit, "SunlitLeaves");
    l.push_back(Shade, "ShadeLeaves");
  }

  NumericMatrix minPsiRhizoComm = aoc["RhizoPsi"];
  NumericMatrix minPsiRhizo = copyNumericMatrix(minPsiRhizoComm, numCohorts, nlayers);
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  l.push_back(minPsiRhizo, "RhizoPsi");

  
  if(control["subdailyResults"]) {
    NumericMatrix soilLayerExtractInstComm = aoc["ExtractionInst"];
    NumericMatrix soilLayerExtractInst = copyNumericMatrix(soilLayerExtractInstComm, nlayers, ntimesteps);
    soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
    l.push_back(soilLayerExtractInst, "ExtractionInst");
    List PlantsInst = copyPlantsInstOutput(as<List>(aoc["PlantsInst"]), x);
    l.push_back(PlantsInst, "PlantsInst");
    l.push_back(clone(as<DataFrame>(aoc["RadiationInputInst"])), "RadiationInputInst");
    List SunlitInst = copyLeavesInstOutput(as<List>(aoc["SunlitLeavesInst"]), x);
    List ShadeInst = copyLeavesInstOutput(as<List>(aoc["ShadeLeavesInst"]), x);
    l.push_back(SunlitInst, "SunlitLeavesInst");
    l.push_back(ShadeInst, "ShadeLeavesInst");
  }
  List lwrExtinctionListComm = aoc["LWRExtinction"];
  List lwrExtinctionList(ntimesteps);
  for(int n=0;n<ntimesteps;n++) {
    lwrExtinctionList[n] = copyCommunicationLongWaveRadiation(as<List>(lwrExtinctionListComm[n]), ncanlayers);
  }
  l.push_back(clone(as<List>(aoc["LightExtinction"])), "LightExtinction");
  l.push_back(lwrExtinctionList, "LWRExtinction");
  l.push_back(copyDataFrame(as<DataFrame>(aoc["CanopyTurbulence"]), ncanlayers), "CanopyTurbulence");
  
  if(control["fireHazardResults"]) l.push_back(aoc["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("spwb_day","list");
  return(l);
}

List copyAdvancedGROWTHOutput(List aoc, List x) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  int numCohorts = cohorts.nrow();
  int ntimesteps = control["ndailysteps"];
  bool subdailyCarbonBalance = control["subdailyCarbonBalance"];

  List spwbOut = copyAdvancedSPWBOutput(aoc, x);

  NumericVector standCB = clone(as<NumericVector>(aoc["CarbonBalance"]));
  
  List l = List::create(_["cohorts"] = spwbOut["cohorts"],
                        _["topography"] = spwbOut["topography"],
                        _["weather"] = spwbOut["weather"],
                        _["WaterBalance"] = spwbOut["WaterBalance"], 
                        _["EnergyBalance"] = spwbOut["EnergyBalance"],
                        _["CarbonBalance"] = standCB);
  
  if(control["soilResults"]) {
    l.push_back(spwbOut["Soil"], "Soil");
  }
  if(control["standResults"]){
    l.push_back(spwbOut["Stand"], "Stand");
  }
  if(control["plantResults"]){
    l.push_back(spwbOut["Plants"], "Plants");
  }
  if(control["labileCarbonBalanceResults"]) {
    DataFrame labileCarbonBalance = copyDataFrame(as<DataFrame>(aoc["LabileCarbonBalance"]), numCohorts);
    labileCarbonBalance.attr("row.names") = above.attr("row.names");
    l.push_back(labileCarbonBalance, "LabileCarbonBalance");
  }
  DataFrame plantBiomassBalance = copyDataFrame(as<DataFrame>(aoc["PlantBiomassBalance"]), numCohorts);
  plantBiomassBalance.attr("row.names") = above.attr("row.names");
  l.push_back(plantBiomassBalance, "PlantBiomassBalance");
  
  if(control["plantStructureResults"]) {
    DataFrame plantStructure = copyDataFrame(as<DataFrame>(aoc["PlantStructure"]), numCohorts);
    plantStructure.attr("row.names") = above.attr("row.names");
    l.push_back(plantStructure, "PlantStructure");
  }
  if(control["growthMortalityResults"]) {
    DataFrame growthMortality  = copyDataFrame(as<DataFrame>(aoc["GrowthMortality"]), numCohorts);
    growthMortality.attr("row.names") = above.attr("row.names");
    l.push_back(growthMortality, "GrowthMortality");
  }
  
  l.push_back(spwbOut["RhizoPsi"], "RhizoPsi");
  
  if(control["leafResults"]){
    l.push_back(spwbOut["SunlitLeaves"], "SunlitLeaves");
    l.push_back(spwbOut["ShadeLeaves"], "ShadeLeaves");
  }
  if(control["subdailyResults"]) {
    l.push_back(spwbOut["ExtractionInst"], "ExtractionInst");
    l.push_back(spwbOut["PlantsInst"], "PlantsInst");
    l.push_back(spwbOut["RadiationInputInst"], "RadiationInputInst");
    l.push_back(spwbOut["SunlitLeavesInst"], "SunlitLeavesInst");
    l.push_back(spwbOut["ShadeLeavesInst"], "ShadeLeavesInst");
  }
  if(subdailyCarbonBalance) {
    List labileCBInstComm = aoc["LabileCarbonBalanceInst"];
    
    NumericMatrix LabileCarbonBalanceInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["LabileCarbonBalance"]), numCohorts, ntimesteps);  
    NumericMatrix GrossPhotosynthesisInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["GrossPhotosynthesis"]), numCohorts, ntimesteps);  
    NumericMatrix MaintenanceRespirationInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["MaintenanceRespiration"]), numCohorts, ntimesteps);    
    NumericMatrix GrowthCostsInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["GrowthCosts"]), numCohorts, ntimesteps);    
    NumericMatrix RootExudationInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["RootExudation"]), numCohorts, ntimesteps);    
    NumericMatrix PlantSugarLeafInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarLeaf"]), numCohorts, ntimesteps);  
    NumericMatrix PlantStarchLeafInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["StarchLeaf"]), numCohorts, ntimesteps);  
    NumericMatrix PlantSugarSapwoodInst= copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarSapwood"]), numCohorts, ntimesteps);  
    NumericMatrix PlantStarchSapwoodInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["StarchSapwood"]), numCohorts, ntimesteps);  
    NumericMatrix PlantSugarTransportInst = copyNumericMatrix(as<NumericMatrix>(labileCBInstComm["SugarTransport"]), numCohorts, ntimesteps);  
    
    GrossPhotosynthesisInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    MaintenanceRespirationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    LabileCarbonBalanceInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    GrowthCostsInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantStarchLeafInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantStarchSapwoodInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    PlantSugarTransportInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    RootExudationInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
    
    List labileCBInst = List::create(
      _["GrossPhotosynthesis"] = GrossPhotosynthesisInst,
      _["MaintenanceRespiration"] = MaintenanceRespirationInst,
      _["GrowthCosts"] = GrowthCostsInst,
      _["RootExudation"] = RootExudationInst,
      _["LabileCarbonBalance"] = LabileCarbonBalanceInst,
      _["SugarLeaf"] = PlantSugarLeafInst,
      _["StarchLeaf"] = PlantStarchLeafInst,
      _["SugarSapwood"] = PlantSugarSapwoodInst,
      _["StarchSapwood"] = PlantStarchSapwoodInst,
      _["SugarTransport"] = PlantSugarTransportInst
    );
    l.push_back(labileCBInst,"LabileCarbonBalanceInst");
  }
  l.push_back(spwbOut["LightExtinction"], "LightExtinction");
  l.push_back(spwbOut["CanopyTurbulence"], "CanopyTurbulence");
  if(control["fireHazardResults"]) l.push_back(spwbOut["FireHazard"], "FireHazard");
  l.attr("class") = CharacterVector::create("growth_day","list");
  return(l);
}

List copySPWBOutput(List internalCommunication, List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List modelOutput;
  if(transpirationMode =="Granier") {
    List modelOutputComm = internalCommunication["basicSPWBOutput"];
    modelOutput = copyBasicSPWBOutput(modelOutputComm, x);
  } else {
    List modelOutputComm = internalCommunication["advancedSPWBOutput"];
    modelOutput = copyAdvancedSPWBOutput(modelOutputComm, x);
  }
  return(modelOutput);
}
List copyGROWTHOutput(List internalCommunication, List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List modelOutput;
  if(transpirationMode =="Granier") {
    List modelOutputComm = internalCommunication["basicGROWTHOutput"];
    modelOutput = copyBasicGROWTHOutput(modelOutputComm, x);
  } else {
    List modelOutputComm = internalCommunication["advancedGROWTHOutput"];
    modelOutput = copyAdvancedGROWTHOutput(modelOutputComm, x);
  }
  return(modelOutput);
}

//' Internal communication
//'
//' Functions for internal communication. Not to be called by users.
//' 
//' @param internalCommunication List for internal communication.
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' @param model String for model, either "spwb" or "growth".
//' 
//' @name communication
//' @keywords internal
// [[Rcpp::export(".copy_model_output")]]
List copyModelOutput(List internalCommunication, List x, String model) {
  List out;
  if(model=="spwb") {
    out= copySPWBOutput(internalCommunication, x);
  } else {
    out= copyGROWTHOutput(internalCommunication, x);
  }
  return(out);
}


