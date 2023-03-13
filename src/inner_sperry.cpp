#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "photosynthesis.h"
using namespace Rcpp;

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
//'   \item{\code{Cost}: Cost function [0-1].}
//'   \item{\code{Gain}: Gain function [0-1].}
//'   \item{\code{Profit}: Profit function [0-1].}
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
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize soil with default soil params (4 layers)
//' examplesoil = soil(defaultSoilParams(4))
//' 
//' #Initialize control parameters
//' control = defaultControl(transpirationMode="Sperry")
//' 
//' #Initialize input
//' x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
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
  NumericVector supplyKterm = supplyFunction["kterm"];
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


