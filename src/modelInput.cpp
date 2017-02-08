#include <Rcpp.h>
#include "swb.h"
#include "root.h"
#include "forestutils.h"
#include "hydraulics.h"

using namespace Rcpp;

/***
 * Parameter checking
 */
// [[Rcpp::export(".checkSpeciesParameters")]]
void checkSpeciesParameters(DataFrame SpParams, CharacterVector params) {
  NumericVector values;
  String s;
  for(int i =0;i<params.size();i++){
    s = params[i];
    if(!SpParams.containsElementNamed(params[i])) {
      Rcout << params[i]<<"\n";
      stop("Parameter missing in species params");
    }
  }
}

/**
 *  Prepare Soil Water Balance input
 */
// [[Rcpp::export("swbInput")]]
List swbInput(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  IntegerVector SP = above["SP"];
  NumericVector LAI = above["LAI"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Simple") & (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Simple' or 'Sperry')");
  double fracTotalTreeResistance = control["fracTotalTreeResistance"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector k(numCohorts);
  NumericVector g(numCohorts);
  NumericVector Sgdd(numCohorts);
  for(int c=0;c<numCohorts;c++){
    k[c]=kSP[SP[c]];
    g[c]=gSP[SP[c]];
    Sgdd[c]=SgddSP[SP[c]];
  }
  List abovedf = DataFrame::create(_["SP"]=SP, _["H"]=H, _["CR"]=CR, _["LAI"]=LAI);
  List df = List::create(_["verbose"] =control["verbose"],_["TranspirationMode"] =transpirationMode, 
                         _["above"] = abovedf, 
                         _["V"] = V,
                         _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  if(transpirationMode=="Simple") {
    NumericVector WUESP = SpParams["WUE"];
    NumericVector Psi_ExtractSP = SpParams["Psi_Extract"];
    NumericVector Al2AsSP = SpParams["Al2As"];
    NumericVector Psi_Extract(numCohorts);
    NumericVector WUE(numCohorts);
    for(int c=0;c<numCohorts;c++){
      Psi_Extract[c]=Psi_ExtractSP[SP[c]];
      WUE[c]=WUESP[SP[c]];
    }
    df["Psi_Extract"]=Psi_Extract; 
    df["WUE"]=WUE; 
  } else if(transpirationMode =="Sperry"){
    NumericVector Al2AsSP = SpParams["Al2As"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector Gwmax(numCohorts);
    NumericVector VCstem_kmax(numCohorts);
    NumericVector VCstem_c(numCohorts), VCstem_d(numCohorts);
    NumericVector VCroot_c(numCohorts), VCroot_d(numCohorts);
    NumericVector Vmax298(numCohorts), Jmax298(numCohorts);
    NumericVector dVec = soil["dVec"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VG_n = soil["VG_n"];
    int nlayers = dVec.size();
    NumericMatrix VCroot_kmax(numCohorts, nlayers); 
    NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
    NumericVector Vc;
    for(int c=0;c<numCohorts;c++){
      Vc = V(c,_);
      //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmaxSP[SP[c]], Al2AsSP[SP[c]],H[c]); 
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      Gwmax[c] = GwmaxSP[SP[c]];
      double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      Vmax298[c] =Vmax298SP[SP[c]];
      Jmax298[c] = 1.67*Vmax298[c];
      for(int l=0;l<nlayers;l++) {
        // Rcout<<Vc[l]<<" ";
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                                                        VCstem_kmax[c], VCstem_c[c], VCstem_d[c]);
      }
      // Rcout<<"\n";
    }
    df["Gwmax"]=Gwmax;
    df["VGrhizo_kmax"] = VGrhizo_kmax;
    df["VCroot_kmax"] = VCroot_kmax;
    df["VCroot_c"]=VCroot_c;
    df["VCroot_d"]=VCroot_d;
    df["VCstem_kmax"]=VCstem_kmax;
    df["VCstem_c"]=VCstem_c;
    df["VCstem_d"]=VCstem_d;
    df["Vmax298"]=Vmax298;
    df["Jmax298"]=Jmax298;
  } 
  df["Transpiration"] = NumericVector(numCohorts, 0.0);
  df["Photosynthesis"] = NumericVector(numCohorts, 0.0);
  df.attr("class") = CharacterVector::create("swbInput","list");
  // df.attr("row.names") = seq(1,numCohorts);
  return(df);
}

// [[Rcpp::export("forest2swbInput")]]
List forest2swbInput(List x, List soil, DataFrame SpParams, List control, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  int nlayers = d.size();
  NumericMatrix V(ntree+nshrub,nlayers);
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ = shrubData["Z"];  
  NumericVector Vi;
  for(int i=0;i<ntree;i++) {
    Vi = ldrRS_one(treeZ50[i], treeZ95[i],d);
    V(i,_) = Vi;
  }
  for(int i=0;i<nshrub;i++) {
    Vi = conicRS_one(shrubZ[i],d);
    V(ntree+i,_) = Vi;
  }
  DataFrame above = forest2aboveground(x, SpParams, gdd);
  return(swbInput(above,  V, soil, SpParams, control));
}



/**
 *  Prepare Forest growth input
 */
// [[Rcpp::export("growthInput")]]
List growthInput(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  IntegerVector SP = above["SP"];
  NumericVector LAI = above["LAI"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Simple") & (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Simple' or 'Sperry')");
  double fracTotalTreeResistance = control["fracTotalTreeResistance"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector k(numCohorts);
  NumericVector g(numCohorts);
  NumericVector Sgdd(numCohorts);
  for(int c=0;c<numCohorts;c++){
    k[c]=kSP[SP[c]];
    g[c]=gSP[SP[c]];
    Sgdd[c]=SgddSP[SP[c]];
  }
  List df = List::create(_["verbose"] =control["verbose"],_["TranspirationMode"] =transpirationMode, 
                         _["above"] = above, 
                         _["V"] = V,
                         _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  if(transpirationMode=="Simple") {
    NumericVector Psi_ExtractSP = SpParams["Psi_Extract"];
    NumericVector Psi_Extract(numCohorts);
    for(int c=0;c<numCohorts;c++){
      Psi_Extract[c]=Psi_ExtractSP[SP[c]];
    }
    df["Psi_Extract"]=Psi_Extract; 
  } else if(transpirationMode =="Sperry"){
    NumericVector Al2AsSP = SpParams["Al2As"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector Gwmax(numCohorts);
    NumericVector xylem_kmax(numCohorts), Al2As(numCohorts);
    NumericVector VCstem_kmax(numCohorts);
    NumericVector VCstem_c(numCohorts), VCstem_d(numCohorts);
    NumericVector VCroot_c(numCohorts), VCroot_d(numCohorts);
    NumericVector Vmax298(numCohorts), Jmax298(numCohorts);
    NumericVector dVec = soil["dVec"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VG_n = soil["VG_n"];
    int nlayers = dVec.size();
    NumericMatrix VCroot_kmax(numCohorts, nlayers); 
    NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
    NumericVector Vc;
    for(int c=0;c<numCohorts;c++){
      Vc = V(c,_);
      xylem_kmax[c] = xylem_kmaxSP[SP[c]];
      Al2As[c] = Al2AsSP[SP[c]];
      //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmax[c], Al2As[c],H[c]); 
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      Gwmax[c] = GwmaxSP[SP[c]];
      double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      Vmax298[c] =Vmax298SP[SP[c]];
      Jmax298[c] = 1.67*Vmax298[c];
      for(int l=0;l<nlayers;l++) {
        // Rcout<<Vc[l]<<" ";
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                                                        VCstem_kmax[c], VCstem_c[c], VCstem_d[c]);
      }
      // Rcout<<"\n";
    }
    df["Al2As"] = Al2As;
    df["Gwmax"]=Gwmax;
    df["VGrhizo_kmax"] = VGrhizo_kmax;
    df["VCroot_kmax"] = VCroot_kmax;
    df["VCroot_c"]=VCroot_c;
    df["VCroot_d"]=VCroot_d;
    df["xylem_kmax"] = xylem_kmax;
    df["VCstem_kmax"]=VCstem_kmax;
    df["VCstem_c"]=VCstem_c;
    df["VCstem_d"]=VCstem_d;
    df["Vmax298"]=Vmax298;
    df["Jmax298"]=Jmax298;
  } 
  df["Transpiration"] = NumericVector(numCohorts, 0.0);
  df["Photosynthesis"] = NumericVector(numCohorts, 0.0);
  df.attr("class") = CharacterVector::create("growthInput","list");
  // df.attr("row.names") = seq(1,numCohorts);
  return(df);
}
// [[Rcpp::export("forest2growthInput")]]
List forest2growthInput(List x, List soil, DataFrame SpParams, List control, double gdd = NA_REAL) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  int nlayers = d.size();
  NumericMatrix V(ntree+nshrub,nlayers);
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ = shrubData["Z"];  
  NumericVector Vi;
  for(int i=0;i<ntree;i++) {
    Vi = ldrRS_one(treeZ50[i], treeZ95[i],d);
    V(i,_) = Vi;
  }
  for(int i=0;i<nshrub;i++) {
    Vi = conicRS_one(shrubZ[i],d);
    V(ntree+i,_) = Vi;
  }
  DataFrame above = forest2aboveground(x, SpParams, gdd);
  return(growthInput(above,  V, soil, SpParams, control));
}
