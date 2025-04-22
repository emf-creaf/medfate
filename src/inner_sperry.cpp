#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "photosynthesis.h"
#include "biophysicsutils.h"
#include "hydraulics.h"
#include "soil.h"
#include "tissuemoisture.h"
using namespace Rcpp;

const double eps_xylem = 1e3; // xylem elastic modulus (1 GPa = 1000 MPa)

List initSperryNetwork(int c,
                       DataFrame internalWater, DataFrame paramsTranspiration, DataFrame paramsWaterStorage,
                       NumericVector VCroot_kmax, NumericVector VGrhizo_kmax,
                       NumericVector psiSoil, NumericVector VG_n, NumericVector VG_alpha,
                       List control,
                       double sapFluidityDay = 1.0) {
  
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_kmax"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_kmax"]);
  NumericVector VCleafapo_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleafapo_kmax"]);
  NumericVector kleaf_symp = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["kleaf_symp"]);
  
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCstem_d"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCleaf_d"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector VCroot_c = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_c"]);
  NumericVector VCroot_d = Rcpp::as<Rcpp::NumericVector>(paramsTranspiration["VCroot_d"]);
  
  double leaf_plc = LeafPLCVEC[c];
  if(!control["leafCavitationEffects"]) leaf_plc = 0.0;
  List HN = List::create(_["numericParams"] = control["numericParams"], 
                         _["stemCavitationEffects"] = control["stemCavitationEffects"],
                         _["leafCavitationEffects"] = control["leafCavitationEffects"],
                         _["psisoil"] = psiSoil,
                         _["krhizomax"] = VGrhizo_kmax,_["nsoil"] = VG_n,_["alphasoil"] = VG_alpha,
                         _["krootmax"] = sapFluidityDay*VCroot_kmax, _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                         _["kstemmax"] = sapFluidityDay*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                         _["kleafmax"] = sapFluidityDay*VCleaf_kmax[c], 
                         _["kleafapomax"] = sapFluidityDay*VCleafapo_kmax[c], 
                         _["kleafsymp"] = sapFluidityDay*kleaf_symp[c], 
                         _["leafc"] = VCleaf_c[c], _["leafd"] = VCleaf_d[c],
                         _["PLCstem"] = StemPLCVEC[c],_["PLCleaf"] = leaf_plc);
  
  return(HN);
}

//' @rdname hydraulics_supplyfunctions
//' @keywords internal
// [[Rcpp::export("hydraulics_initSperryNetworks")]]
List initSperryNetworks(List x) {
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  
  List control = Rcpp::as<Rcpp::List>(x["control"]);
  
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  DataFrame paramsTranspiration = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  NumericVector psiSoil = psi(soil, "VG");
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  
  int numCohorts = internalWater.nrow();
  List networks(numCohorts);
  for(int c = 0;c<numCohorts;c++) {
    networks[c] = initSperryNetwork(c, 
                                    internalWater, paramsTranspiration, paramsWaterStorage,
                                    VCroot_kmax(c,_), VGrhizo_kmax(c,_),
                                    psiSoil, VG_n, VG_alpha, 
                                    control, 1.0);
  }
  networks.attr("names") = above.attr("row.names");
  return(networks);
}

List profitMaximization2(List supplyFunction, int initialPos,
                         double Catm, double Patm, double Tair, double vpa, double u, 
                         double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                         double leafWidth, double refLeafArea,
                         double Gswmin, double Gswmax) {
  
  NumericVector E = supplyFunction["E"];
  NumericVector dEdP = supplyFunction["dEdP"];
  NumericVector leafPsi = supplyFunction["psiLeaf"];
  int nsteps = E.size();
  
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  
  int valid = 1;
  for(int i=0;i<nsteps-1;i++) {
    if((i>0) && (E[i] > E[i-1])) valid++; 
    mindEdp = std::min(mindEdp, dEdP[i]);
    maxdEdp = std::max(maxdEdp, dEdP[i]);
  }
  nsteps = valid;
  // Rcout<< valid<< " "<< maxdEdp << " "<<mindEdp<< " "<< mindEdp/maxdEdp<<"\n";
  // initial pos cannot be over the valid steps
  initialPos = std::min(initialPos, nsteps-1);
  
  //Evaluate photosynthesis and profit at Agmax
  NumericVector photoAgMax = leafPhotosynthesisOneFunction2(E[nsteps-1], leafPsi[nsteps - 1], Catm, Patm, Tair, vpa, u, 
                                                            SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                            leafWidth, refLeafArea, false);
  double Agmax = photoAgMax["GrossPhotosynthesis"];
  double profitAgMax = (1.0 - ((maxdEdp - dEdP[nsteps-1])/(maxdEdp - mindEdp))); 
  
  
  //Photosynthesis and profit maximization at current value of initialPos
  NumericVector photoInitial = leafPhotosynthesisOneFunction2(E[initialPos], leafPsi[initialPos], Catm, Patm, Tair, vpa, u, 
                                                              SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                              leafWidth, refLeafArea, false);
  double AgInitial = photoInitial["GrossPhotosynthesis"];
  double profitInitial = (AgInitial/Agmax) - ((maxdEdp - dEdP[initialPos])/(maxdEdp - mindEdp)); 
  if(Agmax ==0.0) profitInitial = - ((maxdEdp - dEdP[initialPos])/(maxdEdp - mindEdp)); //To avoid 0/0 division
  
  NumericVector photoPrev, photoNext, photoCurrent;
  int iPos = initialPos;
  double profitPrev, profitNext, profitCurrent;
  
  photoCurrent = photoInitial;
  profitCurrent = profitInitial;
  
  int prevStep = 0;
  bool cont = true;
  // Rcout<< " initialPos " << initialPos;
  while(cont) {
    // Rcout<<".";
    if((iPos>0) && (prevStep <= 0)) {
      photoPrev = leafPhotosynthesisOneFunction2(E[iPos-1], leafPsi[iPos-1], Catm, Patm, Tair, vpa, u, 
                                                 SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                 leafWidth, refLeafArea, false);
      double AgPrev = photoPrev["GrossPhotosynthesis"];
      profitPrev = (AgPrev/Agmax) - ((maxdEdp - dEdP[iPos-1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitPrev = - ((maxdEdp - dEdP[iPos-1])/(maxdEdp - mindEdp)); 
    } else {
      photoPrev = photoCurrent;
      profitPrev = profitCurrent;
    }
    if((iPos < (nsteps-1)) && (prevStep >= 0)) {
      photoNext = leafPhotosynthesisOneFunction2(E[iPos+1], leafPsi[iPos+1], Catm, Patm, Tair, vpa, u, 
                                                 SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                 leafWidth, refLeafArea, false);
      double AgNext = photoNext["GrossPhotosynthesis"];
      profitNext = (AgNext/Agmax) - ((maxdEdp - dEdP[iPos+1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitNext = - ((maxdEdp - dEdP[iPos+1])/(maxdEdp - mindEdp)); 
    } else {
      photoNext = photoAgMax;
      profitNext = profitAgMax;
    }
    bool selDecr = ((profitPrev >= profitCurrent) && (photoPrev["Gsw"] >= Gswmin)) || (photoCurrent["Gsw"]>Gswmax);
    selDecr = selDecr && (prevStep<=0) && (iPos>0);
    bool selIncr = ((profitNext > profitCurrent) && (photoNext["Gsw"] <= Gswmax)) || (photoCurrent["Gsw"]<Gswmin);
    selIncr = selIncr && (prevStep>=0) && (iPos<(nsteps-1));
    if(selDecr) {
      profitCurrent = profitPrev;
      photoCurrent = photoPrev;
      iPos = iPos-1;
      prevStep = -1;
      // Rcout <<"-";
    } else if(selIncr) {
      profitCurrent = profitNext;
      photoCurrent = photoNext;
      iPos = iPos+1;
      // Rcout <<"+";
      prevStep = 1;
    } else {
      cont = false;
    }
  }
  // Rcout <<"\n";
  // Rcout<< " finalPos " << iPos<< " final profit "<< profitCurrent <<"\n"; 
  
  return(List::create(Named("photosynthesisFunction") = photoCurrent,
                      Named("Profit") = profitCurrent,
                      Named("iMaxProfit")=iPos));
}


//' Stomatal regulation
//' 
//' Set of high-level functions used in the calculation of stomatal conductance and transpiration. 
//' Function \code{transp_profitMaximization} calculates gain and cost functions, 
//' as well as profit maximization from supply and photosynthesis input functions. 
//' Function \code{transp_stomatalRegulationPlot} produces a plot with the cohort supply functions against water potential 
//' and a plot with the cohort photosynthesis functions against water potential, both with the maximum profit values indicated.
//' 
//' @param supplyFunction Water supply function (see \code{\link{hydraulics_supplyFunctionNetwork}}).
//' @param photosynthesisFunction Function returned by \code{photo_photosynthesisFunction()}.
//' @param Gswmin,Gswmax Minimum and maximum stomatal conductance to water vapour (mol·m-2·s-1).
//' 
//' @return
//' Function \code{transp_profitMaximization} returns a list with the following elements:
//' \itemize{
//'   \item{\code{Cost}: Cost function \[0-1\].}
//'   \item{\code{Gain}: Gain function \[0-1\].}
//'   \item{\code{Profit}: Profit function \[0-1\].}
//'   \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
//' }
//' 
//' @references
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{transp_transpirationSperry}}, \code{\link{hydraulics_supplyFunctionNetwork}}, \code{\link{biophysics_leafTemperature}}, \code{\link{photo_photosynthesis}}, \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Define soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' 
//' #Initialize control parameters
//' control <- defaultControl(transpirationMode="Sperry")
//' 
//' #Initialize input
//' x2 <- spwbInput(exampleforest,examplesoil, SpParamsMED, control)
//' 
//' # Stomatal VPD curve and chosen value for the 12th time step at day 100
//' transp_stomatalRegulationPlot(x2, examplemeteo, day=100, timestep = 12,
//'                               latitude = 41.82592, elevation = 100, type="VPD")
//'  
//' @name transp_stomatalregulation
// [[Rcpp::export("transp_profitMaximization")]]
List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax) {
  NumericVector supplyE = supplyFunction["E"];
  NumericVector supplydEdp = supplyFunction["dEdP"];
  NumericVector Ag = photosynthesisFunction["GrossPhotosynthesis"];
  NumericVector leafTemp = photosynthesisFunction["LeafTemperature"];
  NumericVector leafVPD = photosynthesisFunction["LeafVPD"];
  NumericVector Gsw = photosynthesisFunction["Gsw"];
  int nsteps = supplydEdp.size();
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  double Agmax = 0.0;
  //Find valid limits according to stomatal conductance
  int ini = 0, fin = nsteps-1;
  
  for(int i=ini;i<fin;i++) {
    mindEdp = std::min(mindEdp, supplydEdp[i]);
    maxdEdp = std::max(maxdEdp, supplydEdp[i]);
    Agmax = std::max(Agmax, Ag[i]);
  }
  
  //Evaluate profit for valid steps
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  for(int i=ini;i<fin;i++) {
    gain[i] = Ag[i]/Agmax;
    cost[i] = (maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp); 
    profit[i] = gain[i]-cost[i];
  }
  
  while((Gsw[ini]<=Gswmin) && (ini<fin)) ini++;
  while((Gsw[fin]>=Gswmax) && (fin>ini)) fin--; 
  
  //Ensure that ini <=fin
  ini = std::min(ini, fin);
  fin = std::max(ini,fin);
  
  int imaxprofit=ini;
  double maxprofit=profit[ini];
  if(fin>ini) {
    for(int i=ini+1;i<=fin;i++){
      if((profit[i]>maxprofit)) {
        maxprofit = profit[i];
        imaxprofit = i;
      }
    }
  }
  // Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
  if((Gsw[imaxprofit] > Gswmax) && (imaxprofit>ini)) {
    Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
    for(int i=0;i<Gsw.size();i++) {
      Rcout<< i << " Gsw "<< Gsw[i] << " supplyE "<< supplyE[i] << " leafT "<< leafTemp[i]<< " leafVPD "<< leafVPD[i]  << "\n";
    }
    stop("Gsw > Gswmax");
  }
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}


void copyRhizoPsi(int c, int iPM, 
                  NumericMatrix RhizoPsi, NumericMatrix RhizoPsiMAT,
                  LogicalMatrix layerConnected, 
                  List RHOP, List layerConnectedPools,
                  NumericVector VCroot_c, NumericVector VCroot_d,  
                  bool plantWaterPools) {
  int nlayers = layerConnected.ncol();
  int numCohorts = layerConnected.nrow();
  
  if(!plantWaterPools) {
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected(c,l)) {
        RhizoPsiMAT(c,l) = RhizoPsi(iPM,cl);
        cl++;
      } 
    }
  } else {
    NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
    LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
    NumericVector rplv(numCohorts,NA_REAL);
    NumericVector vplv(numCohorts,NA_REAL);
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      int clj = 0;
      for(int j=0;j<numCohorts;j++) {
        if(layerConnectedCoh(j,l)) {
          rplv[clj] = RhizoPsi(iPM,cl);
          vplv[clj] = RHOPcoh(j,l);
          cl++;
          clj++;
        }
      }
      NumericVector pv(clj,NA_REAL);
      NumericVector vv(clj,NA_REAL);
      for(int j=0;j<clj;j++) {
        pv[j] = rplv[j];
        vv[j] = vplv[j];
        // Rcout<< c << " "<< l <<" "<<j<< " "<< pv[j]<< " "<< vv[j]<<"\n";
      }
      RhizoPsiMAT(c,l) = averagePsi(pv,vv,VCroot_c[c], VCroot_d[c]);
      // Rcout<< c << " "<<l<< " "<< RhizoPsiMAT(c,l)<<"\n";
    }
  }
}


void innerSperry(List x, List input, List output, int n, double tstep, 
                 bool verbose = false, int stepFunctions = NA_INTEGER) {
  
  // Extract control variables
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  bool leafCavitationEffects = control["leafCavitationEffects"];
  bool stemCavitationEffects = control["stemCavitationEffects"];
  String stemCavitationRecovery = control["stemCavitationRecovery"];
  String leafCavitationRecovery = control["leafCavitationRecovery"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  bool sunlitShade = control["sunlitShade"];

  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  int numCohorts = cohorts.nrow();

  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  
  DataFrame soil = Rcpp::as<Rcpp::DataFrame>(x["soil"]);
  NumericVector Ws = soil["W"]; //Access to soil state variable
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  // NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  // NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  // NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  int nlayers = Tsoil.length();
  
  // Extract parameters
  // Rcout<<"params\n";
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  NumericVector zlow = canopyParams["zlow"];
  NumericVector zmid = canopyParams["zmid"];
  NumericVector zup = canopyParams["zup"];
  NumericVector Tair = canopyParams["Tair"];
  NumericVector VPair = canopyParams["VPair"];
  NumericVector Cair = canopyParams["Cair"];
  
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCroot_kmax_sum = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCroot_kmax"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  
  //Extract internal variables
  // Rcout<<"internal\n";
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector LeafPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPLC"]);
  NumericVector StemPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPsi"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  
  //Water pools
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  
  //Extract output to be filled
  // Rcout<<"output\n";
  List outPhotoSunlit = output["PhotoSunlitFunctions"];
  List outPhotoShade = output["PhotoShadeFunctions"];
  List outPMSunlit = output["PMSunlitFunctions"];
  List outPMShade = output["PMShadeFunctions"];
  
  NumericMatrix SoilWaterExtract = Rcpp::as<Rcpp::NumericMatrix>(output["Extraction"]);
  List ExtractionPools = Rcpp::as<Rcpp::List>(output["ExtractionPools"]);
  NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(output["ExtractionInst"]);

  NumericMatrix minPsiRhizo = Rcpp::as<Rcpp::NumericMatrix>(output["RhizoPsi"]);

  // Rcout<<"Plants\n";
  List Plants = output["Plants"];
  NumericVector PWB = Plants["WaterBalance"];
  NumericVector Eplant = Plants["Transpiration"];
  NumericVector Agplant = Plants["GrossPhotosynthesis"];
  NumericVector Anplant = Plants["NetPhotosynthesis"];
  NumericVector minLeafPsi = Plants["LeafPsiMin"];
  NumericVector maxLeafPsi = Plants["LeafPsiMax"];
  NumericVector minStemPsi = Plants["StemPsi"];
  NumericVector minRootPsi = Plants["RootPsi"];
  
  
  // Rcout<<"PlantsInst\n";
  List PlantsInst = output["PlantsInst"];
  NumericMatrix Einst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["E"]);
  NumericMatrix Aginst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["Ag"]);
  NumericMatrix Aninst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["An"]);
  NumericMatrix dEdPInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["dEdP"]);
  NumericMatrix PWBinst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["PWB"]);
  NumericMatrix StemSympRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympRWC"]);
  NumericMatrix LeafSympRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympRWC"]);
  NumericMatrix StemRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemRWC"]);
  NumericMatrix LeafRWCInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafRWC"]);
  NumericMatrix StemPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemPsi"]);
  NumericMatrix LeafPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafPsi"]);
  NumericMatrix RootPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["RootPsi"]);
  NumericMatrix StemPLC = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemPLC"]);
  NumericMatrix LeafPLC = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafPLC"]);
  NumericMatrix StemSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["StemSympPsi"]);
  NumericMatrix LeafSympPsiInst = Rcpp::as<Rcpp::NumericMatrix>(PlantsInst["LeafSympPsi"]);
  
  List ShadeInst = output["ShadeLeavesInst"];
  NumericMatrix LAI_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["LAI"]);
  NumericMatrix Vmax298_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Vmax298"]);
  NumericMatrix Jmax298_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Jmax298"]);
  NumericMatrix SWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Abs_SWR"]);
  NumericMatrix PAR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Abs_PAR"]);
  NumericMatrix LWR_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Net_LWR"]);
  NumericMatrix Ag_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Ag"]);
  NumericMatrix An_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["An"]);
  NumericMatrix E_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["E"]);
  NumericMatrix VPD_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["VPD"]);
  NumericMatrix Psi_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Psi"]);
  NumericMatrix Temp_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Temp"]);
  NumericMatrix GSW_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Gsw"]);
  NumericMatrix Ci_SH = Rcpp::as<Rcpp::NumericMatrix>(ShadeInst["Ci"]);
  
  List SunlitInst = output["SunlitLeavesInst"];
  NumericMatrix LAI_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["LAI"]);
  NumericMatrix Vmax298_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Vmax298"]);
  NumericMatrix Jmax298_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Jmax298"]);
  NumericMatrix SWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Abs_SWR"]);
  NumericMatrix PAR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Abs_PAR"]);
  NumericMatrix LWR_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Net_LWR"]);
  NumericMatrix Ag_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ag"]);
  NumericMatrix An_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["An"]);
  NumericMatrix E_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["E"]);
  NumericMatrix VPD_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["VPD"]);
  NumericMatrix Psi_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Psi"]);
  NumericMatrix Temp_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Temp"]);
  NumericMatrix GSW_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Gsw"]);
  NumericMatrix Ci_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ci"]);
  
  //Extract input  
  // Rcout<<"input\n";
  NumericVector zWind = input["zWind"];
  double Patm = input["Patm"];
  double f_dry = input["f_dry"];
  IntegerVector iLayerCohort = input["iLayerCohort"];
  IntegerVector iLayerSunlit = input["iLayerSunlit"];
  IntegerVector iLayerShade = input["iLayerShade"];
  IntegerVector iPMSunlit = input["iPMSunlit"];
  IntegerVector iPMShade = input["iPMShade"];
  IntegerVector nlayerscon = input["nlayerscon"];
  LogicalMatrix layerConnected = input["layerConnected"];
  List layerConnectedPools = input["layerConnectedPools"];
  List supply = input["supply"];
  
  for(int c=0;c<numCohorts;c++) { //Plant cohort loop
    if((LAIphe[c]>0.0) && (LeafPLCVEC[c] < 0.999)) { //Process transpiration and photosynthesis only if there are some leaves

      NumericVector fittedE, dEdP;
      NumericVector LeafPsi, psiRootCrown;
      
      //Retrieve supply functions
      List sFunctionAbove = supply[c];
      List sFunctionBelow = supply[c];
      
      //Retrieve transpiration, LeafPsi and dEdP vectors
      fittedE = sFunctionAbove["E"];
      dEdP = sFunctionAbove["dEdP"];
      LeafPsi = sFunctionAbove["psiLeaf"];
      
      //Get info from sFunctionAbove
      psiRootCrown = sFunctionAbove["psiRootCrown"];
      
      if(fittedE.size()>0) {
        double TPhase_gmin = control["TPhase_gmin"]; 
        double Q10_1_gmin = control["Q10_1_gmin"]; 
        double Q10_2_gmin = control["Q10_2_gmin"];
        //Here canopy air temperature is used instead of Tleaf, because leaf energy balance has not been closed
        double gmin_SL = gmin(Tair[iLayerSunlit[c]], Gswmin[c], TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        double gmin_SH = gmin(Tair[iLayerShade[c]], Gswmin[c], TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        // Rcout<< gmin_SL<< " " << gmin_SH << " " << Gswmax[c]<<"\n";
        //Photosynthesis function for sunlit and shade leaves
        List PMSunlit = profitMaximization2(sFunctionAbove, iPMSunlit[c], 
                                            Cair[iLayerSunlit[c]], Patm,
                                            Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], zWind[iLayerSunlit[c]], 
                                            SWR_SL(c,n), LWR_SL(c,n), irradianceToPhotonFlux(PAR_SL(c,n)), 
                                            Vmax298_SL(c,n), Jmax298_SL(c,n), leafWidth[c], LAI_SL(c,n),
                                            gmin_SL, Gswmax[c]);
        NumericVector photoSunlit = PMSunlit["photosynthesisFunction"];
        iPMSunlit[c] = PMSunlit["iMaxProfit"];
        List PMShade = profitMaximization2(sFunctionAbove, iPMShade[c], 
                                           Cair[iLayerShade[c]], Patm,
                                           Tair[iLayerShade[c]], VPair[iLayerShade[c]], zWind[iLayerShade[c]], 
                                           SWR_SH(c,n), LWR_SH(c,n), irradianceToPhotonFlux(PAR_SH(c,n)), 
                                           Vmax298_SH(c,n), Jmax298_SH(c,n), leafWidth[c], LAI_SH(c,n),
                                           gmin_SH, Gswmax[c]);    
        if(!sunlitShade) PMShade = PMSunlit;
        NumericVector photoShade = PMShade["photosynthesisFunction"];
        iPMShade[c] = PMShade["iMaxProfit"];
        
        
        
        //Store?
        if(!IntegerVector::is_na(stepFunctions)) {
          if(n==stepFunctions) {
            outPhotoSunlit[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerSunlit[c]], Patm,
                                                            Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], 
                                                            zWind[iLayerSunlit[c]], 
                                                            SWR_SL(c,n), LWR_SL(c,n), 
                                                            irradianceToPhotonFlux(PAR_SL(c,n)), 
                                                            Vmax298_SL(c,n), Jmax298_SL(c,n), 
                                                            leafWidth[c], LAI_SL(c,n));
            outPhotoShade[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerShade[c]], Patm,
                                                           Tair[iLayerShade[c]], VPair[iLayerShade[c]], 
                                                           zWind[iLayerShade[c]], 
                                                           SWR_SH(c,n), LWR_SH(c,n), 
                                                           irradianceToPhotonFlux(PAR_SH(c,n)),
                                                           Vmax298_SH(c,n), Jmax298_SH(c,n), 
                                                           leafWidth[c], LAI_SH(c,n));
            outPMSunlit[c] = PMSunlit;
            outPMShade[c] = PMShade;
            if(!sunlitShade) {
              outPMShade[c] = outPMSunlit[c];
              outPhotoShade[c] = outPhotoSunlit[c];
            }
          }
        }
        // Rcout<<iPMSunlit[c]<<" "<<iPMShade[c] <<" "<<GwSunlit[iPMSunlit[c]]<<" "<<GwShade[iPMShade[c]]<<" "<<fittedE[iPMSunlit[c]]<<" "<<fittedE[iPMShade[c]]<<"\n";
        //Get leaf status
        E_SH(c,n) = fittedE[iPMShade[c]];
        E_SL(c,n) = fittedE[iPMSunlit[c]];
        Psi_SH(c,n) = LeafPsi[iPMShade[c]];
        Psi_SL(c,n) = LeafPsi[iPMSunlit[c]];
        An_SH(c,n) = photoShade["NetPhotosynthesis"];
        An_SL(c,n) = photoSunlit["NetPhotosynthesis"];
        Ag_SH(c,n) = photoShade["GrossPhotosynthesis"];
        Ag_SL(c,n) = photoSunlit["GrossPhotosynthesis"];
        Ci_SH(c,n) = photoShade["Ci"];
        Ci_SL(c,n) = photoSunlit["Ci"];
        GSW_SH(c,n)= photoShade["Gsw"];
        GSW_SL(c,n)= photoSunlit["Gsw"];
        VPD_SH(c,n)= photoShade["LeafVPD"];
        VPD_SL(c,n)= photoSunlit["LeafVPD"];
        Temp_SH(c,n)= photoShade["LeafTemperature"];
        Temp_SL(c,n)= photoSunlit["LeafTemperature"];
        
        //Scale photosynthesis
        double Agsum = (Ag_SL(c,n)*LAI_SL(c,n) + Ag_SH(c,n)*LAI_SH(c,n));
        double Ansum = (An_SL(c,n)*LAI_SL(c,n) + An_SH(c,n)*LAI_SH(c,n));
        Aginst(c,n) = (1e-6)*12.01017*Agsum*tstep;
        Aninst(c,n) = (1e-6)*12.01017*Ansum*tstep;
        
        //Average flow from sunlit and shade leaves
        double Eaverage = (E_SL(c,n)*LAI_SL(c,n) + E_SH(c,n)*LAI_SH(c,n))/(LAI_SL(c,n) + LAI_SH(c,n));
        
        
        //Find iPM for  flow corresponding to the  average flow
        double absDiff = 99999999.9;
        int iPM = -1;
        for(int k=0;k<fittedE.size();k++){ //Only check up to the size of fittedE
          double adk = std::abs(fittedE[k]-Eaverage);
          if(adk<absDiff) {
            absDiff = adk;
            iPM = k;
          }
        }
        if(iPM==-1) {
          Rcout<<"\n iPM -1! Eaverage="<< Eaverage << " fittedE.size= "<< fittedE.size()<<" iPMSunlit[c]="<< iPMSunlit[c]<< " fittedE[iPMSunlit[c]]="<<fittedE[iPMSunlit[c]]<<" iPMShade[c]="<<iPMShade[c]<<" fittedE[iPMShade[c]]="<<fittedE[iPMShade[c]]<<"\n";
          stop("");
        }
        
        //Store instantaneous total conductance
        dEdPInst(c,n) = dEdP[iPM];
        
        //Store instantaneous flow and leaf water potential
        EinstVEC[c] = fittedE[iPM]*f_dry;
        LeafPsiVEC[c] = LeafPsi[iPM];
        RootCrownPsiVEC[c] = psiRootCrown[iPM]; 
        
        //Scale from instantaneous flow to water volume in the time step
        Einst(c,n) = EinstVEC[c]*0.001*0.01802*LAIphe[c]*tstep; 
        
        NumericVector Esoilcn(nlayerscon[c],0.0);
        NumericVector ElayersVEC(nlayerscon[c],0.0);
        
        
        //Get info from sFunctionBelow
        NumericMatrix ERhizo = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["ERhizo"]);
        NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["psiRhizo"]);
        
        //Store steady state stem and rootcrown and root surface water potential values
        NumericVector newStemPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionAbove["psiStem"]);
        StemPsiVEC[c] = newStemPsi[iPM]; 
        for(int lc=0;lc<nlayerscon[c];lc++) {
          ElayersVEC[lc] = ERhizo(iPM,lc)*tstep*f_dry; //Scale according to the time step
        }
        
        //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
        copyRhizoPsi(c,iPM, 
                     RhizoPsi, RhizoPsiMAT,
                     layerConnected, 
                     RHOP, layerConnectedPools,
                     VCroot_c, VCroot_d,  
                     plantWaterPools);
        
        StemSympPsiVEC[c] = StemPsiVEC[c]; //Stem symplastic compartment coupled with apoplastic compartment
        LeafSympPsiVEC[c] = LeafPsiVEC[c]; //Leaf symplastic compartment coupled with apoplastic compartment
        
        // Update the leaf and stem PLC
        if(stemCavitationEffects) {
          if(stemCavitationRecovery!="total") {
            StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(StemPsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
          } else { //Immediate refilling
            StemPLCVEC[c] = 1.0 - xylemConductance(StemPsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
          }
        }
        if(leafCavitationEffects) {
          if(leafCavitationRecovery!="total") {
            LeafPLCVEC[c] = std::max(LeafPLCVEC[c], 1.0 - xylemConductance(LeafPsiVEC[c], 1.0, VCleaf_c[c], VCleaf_d[c])); 
          } else { //Immediate refilling
            LeafPLCVEC[c] = 1.0 - xylemConductance(LeafPsiVEC[c], 1.0, VCleaf_c[c], VCleaf_d[c]); 
          }
        }
        
        //Scale soil water extracted from leaf to cohort level
        for(int lc=0;lc<nlayerscon[c];lc++) {
          Esoilcn[lc] = ElayersVEC[lc]*0.001*0.01802*LAIphe[c]; //Scale from flow to water volume in the time step
        }
        
        //Balance between extraction and transpiration
        PWBinst(c,n) = sum(Esoilcn) - Einst(c,n);
        
        //Add step transpiration to daily plant cohort transpiration
        Eplant[c] += Einst(c,n);
        Anplant[c] += Aninst(c,n);
        Agplant[c] += Aginst(c,n);
        //Add PWB
        PWB[c] += PWBinst(c,n); 
        
        
        
        //Copy transpiration and from connected layers to transpiration from soil layers
        //And update soil water content (soil water potential will not be updated until next day!)
        if(!plantWaterPools) {
          int cl = 0;
          for(int l=0;l<nlayers;l++) {
            if(layerConnected(c,l)) {
              SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
              soilLayerExtractInst(l,n) += Esoilcn[cl];
              cl++;
            } 
          }
        } else {
          NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
          LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
          int cl = 0;
          for(int j = 0;j<numCohorts;j++) {
            NumericMatrix ExtractionPoolsCoh = Rcpp::as<Rcpp::NumericMatrix>(ExtractionPools[j]);
            for(int l=0;l<nlayers;l++) {
              if(layerConnectedCoh(j,l)) {
                SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
                soilLayerExtractInst(l,n) += Esoilcn[cl];
                ExtractionPoolsCoh(c,l) += Esoilcn[cl];
                cl++;
              }
            }
          }
        }
      } else {
        if(verbose) Rcout<<"NS!";
        Psi_SH(c,n) = NA_REAL;
        Psi_SL(c,n) = NA_REAL;
        GSW_SH(c,n)= NA_REAL;
        GSW_SL(c,n)= NA_REAL;
        VPD_SH(c,n)= NA_REAL;
        VPD_SL(c,n)= NA_REAL;
        Temp_SH(c,n)= NA_REAL;
        Temp_SL(c,n)= NA_REAL;
      }        
    } else if(LAIlive[c]>0.0) { //Cohorts with living individuals but no LAI should be in equilibrium with soil (i.e. no transpiration)
      
      List sFunctionBelow = supply[c];
      NumericVector  psiRootCrown = sFunctionBelow["psiRootCrown"];
      RootCrownPsiVEC[c] = psiRootCrown[0];
      NumericVector  psiStem1 = sFunctionBelow["psiStem"];
      StemPsiVEC[c] = psiStem1[0];
      StemSympPsiVEC[c] = psiStem1[0];
      NumericVector LeafPsi = sFunctionBelow["psiLeaf"];
      LeafPsiVEC[c] = LeafPsi[0];
      LeafSympPsiVEC[c] = LeafPsi[0];
    }
  } //End of cohort loop
}

