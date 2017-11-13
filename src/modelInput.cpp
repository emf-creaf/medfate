#include <Rcpp.h>
#include "swb.h"
#include "growth.h"
#include "root.h"
#include "forestutils.h"
#include "hydraulics.h"
#include "stdlib.h"

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
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Simple") & (transpirationMode!="Complex")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Simple' or 'Complex')");
  double fracTotalTreeResistance = control["fracTotalTreeResistance"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  String canopyMode = control["canopyMode"];
  
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector k(numCohorts), g(numCohorts), Sgdd(numCohorts);
  for(int c=0;c<numCohorts;c++){
    k[c]=kSP[SP[c]];
    g[c]=gSP[SP[c]];
    Sgdd[c]=SgddSP[SP[c]];
  }
  
  //Soil layer names
  CharacterVector slnames(V.ncol());
  for(int i=0;i<V.ncol();i++) slnames[i] = i+1;
  V.attr("dimnames") = List::create(above.attr("row.names"), slnames);
  
  //Cohort description
  CharacterVector nameSP = SpParams["Name"];
  CharacterVector nsp(numCohorts);
  for(int i=0;i<numCohorts;i++) nsp[i] = nameSP[SP[i]];
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  //Above 
  DataFrame plantsdf = DataFrame::create(_["H"]=H, _["CR"]=CR, 
                                         _["LAI_live"]=LAI_live, _["LAI_expanded"] = LAI_expanded, _["LAI_dead"] = LAI_dead);
  plantsdf.attr("row.names") = above.attr("row.names");
  
  //Base params
  DataFrame paramsBasedf = DataFrame::create(_["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  paramsBasedf.attr("row.names") = above.attr("row.names");
  
 
  List input;
  if(transpirationMode=="Simple") {
    NumericVector WUESP = SpParams["WUE"];
    NumericVector WUE(numCohorts);
    NumericVector Psi_ExtractSP = SpParams["Psi_Extract"];
    NumericVector Psi_Extract(numCohorts);
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector pRootDisc(numCohorts);
    for(int c=0;c<numCohorts;c++){
      Psi_Extract[c]=Psi_ExtractSP[SP[c]];
      WUE[c]=WUESP[SP[c]];
      pRootDisc[c]=pRootDiscSP[SP[c]];
    }
    DataFrame paramsTranspdf = DataFrame::create(_["Psi_Extract"]=Psi_Extract,_["WUE"] = WUE,  _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    List below = List::create(_["V"] = V);
    List numericParams = control["numericParams"];
    List paramsControl = List::create(_["verbose"] =control["verbose"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["cavitationRefill"] = control["cavitationRefill"]);
    List paramsCanopy = List::create(_["gdd"] = 0);
    input = List::create(_["control"] =paramsControl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsTransp"] = paramsTranspdf);
    
  } else if(transpirationMode =="Complex"){
    NumericVector Al2AsSP = SpParams["Al2As"];
    NumericVector GwminSP = SpParams["Gwmin"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector Gwmax(numCohorts), Gwmin(numCohorts);
    NumericVector xylem_kmax(numCohorts), Al2As(numCohorts);
    NumericVector VCstem_kmax(numCohorts);
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector pRootDisc(numCohorts);
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
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmax[c], Al2As[c],H[c],control["taper"]); 
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      pRootDisc[c]=pRootDiscSP[SP[c]];
      Gwmin[c] = GwminSP[SP[c]];
      Gwmax[c] = GwmaxSP[SP[c]];
      double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      Vmax298[c] =Vmax298SP[SP[c]];
      Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c])); //Walker et al 2014
      for(int l=0;l<nlayers;l++) {
        // Rcout<<Vc[l]<<" ";
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                                                        VCstem_kmax[c], VCstem_c[c], VCstem_d[c]);
      }
      // Rcout<<"\n";
    }
    VGrhizo_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    VCroot_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    
    DataFrame paramsTranspdf = DataFrame::create(
      _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax, _["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,_["xylem_kmax"] = xylem_kmax,
      _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    List below = List::create(_["V"] = V,
                              _["VGrhizo_kmax"] = VGrhizo_kmax,
                              _["VCroot_kmax"] = VCroot_kmax);
    List numericParams = control["numericParams"];
    List paramsControl = List::create(_["verbose"] =control["verbose"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["canopyMode"] =canopyMode, 
                                      _["verticalLayerSize"] = control["verticalLayerSize"],
                                      _["hydraulicCostFunction"] = control["hydraulicCostFunction"],
                                      _["cavitationRefill"] = control["cavitationRefill"],
                                      _["taper"] = control["taper"],
                                      _["ndailysteps"] = control["ndailysteps"],
                                      _["numericParams"] = clone(numericParams));
    List paramsCanopy = List::create(_["gdd"] = 0,_["Temp"] = NA_REAL);
    input = List::create(_["control"] = paramsControl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsTransp"] = paramsTranspdf);
  }
  NumericVector tvec =  NumericVector(numCohorts, 0.0);
  tvec.attr("names") = above.attr("row.names");
  input["Transpiration"] = tvec;
  NumericVector pvec =  NumericVector(numCohorts, 0.0);
  pvec.attr("names") = above.attr("row.names");
  input["Photosynthesis"] = pvec;
  NumericVector cvec =  NumericVector(numCohorts, 0.0);
  cvec.attr("names") = above.attr("row.names");
  input["ProportionCavitated"] = cvec;
  // input["WindSpeed"] = NumericVector(numCohorts, 0.0);
  // input["PAR"] = NumericVector(numCohorts, 0.0);
  // input["AbsorbedSWR"] = NumericVector(numCohorts, 0.0);
  input.attr("class") = CharacterVector::create("swbInput","list");
  // df.attr("row.names") = seq(1,numCohorts);
  return(input);
}

// [[Rcpp::export("forest2swbInput")]]
List forest2swbInput(List x, List soil, DataFrame SpParams, List control) {
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
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(swbInput(above,  V, soil, SpParams, control));
}



/**
 *  Prepare Forest growth input
 */
// [[Rcpp::export("growthInput")]]
List growthInput(DataFrame above, NumericVector Z, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  IntegerVector SP = above["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector N = above["N"];
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Simple") & (transpirationMode!="Complex")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Simple' or 'Complex')");
  String storagePool = control["storagePool"];
  if((storagePool!="none") & (storagePool!="one")& (storagePool!="two")) stop("Wrong storage pool ('storagePool' should be 'none', 'one' or 'two')");
  String canopyMode = control["canopyMode"];
  
  
  double fracTotalTreeResistance = control["fracTotalTreeResistance"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  
  NumericVector Al2AsSP = SpParams["Al2As"];
  NumericVector SLASP = SpParams["SLA"];
  NumericVector WoodCSP = SpParams["WoodC"];
  NumericVector WoodDensSP = SpParams["WoodDens"];
  NumericVector RGRmaxSP = SpParams["RGRmax"];
  NumericVector CstoragepmaxSP = SpParams["Cstoragepmax"];
  
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector k(numCohorts), g(numCohorts), Sgdd(numCohorts);
  NumericVector SLA(numCohorts), Al2As(numCohorts), WoodC(numCohorts), WoodDens(numCohorts);
  NumericVector Cstoragepmax(numCohorts), RGRmax(numCohorts);
  
  NumericVector HmaxSP = SpParams["Hmax"];
  NumericVector ZmaxSP = SpParams["Zmax"];
  NumericVector r635SP = SpParams["r635"];
  NumericVector AashSP = SpParams["a_ash"];
  NumericVector AbshSP = SpParams["a_bsh"];
  NumericVector BbshSP = SpParams["b_bsh"];
  NumericVector AcrSP = SpParams["a_cr"];
  NumericVector B1crSP = SpParams["b_1cr"];
  NumericVector B2crSP = SpParams["b_2cr"];
  NumericVector B3crSP = SpParams["b_3cr"];
  NumericVector C1crSP = SpParams["c_1cr"];
  NumericVector C2crSP = SpParams["c_2cr"];
  NumericVector AcwSP = SpParams["a_cw"];
  NumericVector BcwSP = SpParams["b_cw"];
  NumericVector fHDminSP = SpParams["fHDmin"];
  NumericVector fHDmaxSP = SpParams["fHDmax"];
  
  NumericVector Hmax(numCohorts), Zmax(numCohorts);
  NumericVector Aash(numCohorts), Absh(numCohorts), Bbsh(numCohorts), r635(numCohorts);
  NumericVector fHDmin(numCohorts), fHDmax(numCohorts);
  NumericVector Acr(numCohorts), B1cr(numCohorts), B2cr(numCohorts), B3cr(numCohorts), C1cr(numCohorts), C2cr(numCohorts);
  NumericVector Acw(numCohorts), Bcw(numCohorts);
  for(int c=0;c<numCohorts;c++){
    k[c]=kSP[SP[c]];
    g[c]=gSP[SP[c]];
    Sgdd[c]=SgddSP[SP[c]];
    SLA[c]=SLASP[SP[c]];
    Al2As[c]=Al2AsSP[SP[c]];
    WoodDens[c] = WoodDensSP[SP[c]];
    WoodC[c] = WoodCSP[SP[c]];
    Cstoragepmax[c] = std::max(0.05,CstoragepmaxSP[SP[c]]); //Minimum 5%
    RGRmax[c] = RGRmaxSP[SP[c]];

    //Allometries
    Hmax[c] = HmaxSP[SP[c]];
    Zmax[c] = ZmaxSP[SP[c]];
    
    Aash[c] = AashSP[SP[c]];
    Absh[c] = AbshSP[SP[c]];
    Bbsh[c] = BbshSP[SP[c]];
    r635[c] = r635SP[SP[c]];
    Acr[c] = AcrSP[SP[c]];
    B1cr[c] = B1crSP[SP[c]];
    B2cr[c] = B2crSP[SP[c]];
    B3cr[c] = B3crSP[SP[c]];
    C1cr[c] = C1crSP[SP[c]];
    C2cr[c] = C2crSP[SP[c]];
    Acw[c] = AcwSP[SP[c]];
    Bcw[c] = BcwSP[SP[c]];
    fHDmax[c] = fHDmaxSP[SP[c]];
    fHDmin[c] = fHDminSP[SP[c]];
  }
  NumericVector SA(numCohorts);
  for(int c=0;c<numCohorts;c++){
    SA[c] = 10000.0*(LAI_live[c]/(N[c]/10000.0))/Al2AsSP[SP[c]];//Individual SA in cm2/m2
  }
  
  //Soil layer names
  CharacterVector slnames(V.ncol());
  for(int i=0;i<V.ncol();i++) slnames[i] = i+1;
  V.attr("dimnames") = List::create(above.attr("row.names"), slnames);
  Z.attr("names") = above.attr("row.names");
  
  //Cohort description
  CharacterVector nameSP = SpParams["Name"];
  CharacterVector nsp(numCohorts);
  for(int i=0;i<numCohorts;i++) nsp[i] = nameSP[SP[i]];
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  DataFrame plantsdf, paramsGrowthdf;
  if(storagePool=="one") {
    NumericVector fastCstorage(numCohorts);
    for(int c=0;c<numCohorts;c++){
      NumericVector compartments = carbonCompartments(SA[c], LAI_expanded[c], H[c], 
                                                      Z[c], N[c], SLA[c], WoodDens[c], WoodC[c]);
      fastCstorage[c] = 0.5*Cstoragepmax[c]*(compartments[0]+compartments[1]+compartments[2]);//Pool at 50%
    }
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                   _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                   _["SA"] = SA, _["fastCstorage"] = fastCstorage);
    paramsGrowthdf = DataFrame::create(_["SLA"] = SLA, _["Al2As"] = Al2As,
                                       _["WoodDens"] = WoodDens, _["WoodC"] = WoodC,
                                       _["Cstoragepmax"] = Cstoragepmax, _["RGRmax"] = RGRmax);
  } else if(storagePool=="two") {
    NumericVector slowCstorage(numCohorts), fastCstorage(numCohorts);
    for(int c=0;c<numCohorts;c++){
      NumericVector compartments = carbonCompartments(SA[c], LAI_expanded[c], H[c], 
                                                      Z[c], N[c], SLA[c], WoodDens[c], WoodC[c]);
      slowCstorage[c] = 0.5*(Cstoragepmax[c]-0.05)*(compartments[0]+compartments[1]+compartments[2]); //Slow pool at 50%
      fastCstorage[c] = 0.5*0.05*(compartments[0]+compartments[1]+compartments[2]); //Fast pool at 50%
    }
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                 _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                   _["SA"] = SA, _["fastCstorage"] = fastCstorage, _["slowCstorage"] = slowCstorage);
    paramsGrowthdf = DataFrame::create(_["SLA"] = SLA, _["Al2As"] = Al2As,
                                       _["WoodDens"] = WoodDens, _["WoodC"] = WoodC,
                                       _["Cstoragepmax"] = Cstoragepmax, _["RGRmax"] = RGRmax);
  } else {
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                 _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                 _["SA"] = SA);
    paramsGrowthdf = DataFrame::create(_["SLA"] = SLA, _["Al2As"] = Al2As,
                                                 _["WoodDens"] = WoodDens, _["WoodC"] = WoodC,
                                                 _["RGRmax"] = RGRmax);
  }
  plantsdf.attr("row.names") = above.attr("row.names");
  paramsGrowthdf.attr("row.names") = above.attr("row.names");
  
  //Base params
  DataFrame paramsBasedf = DataFrame::create(_["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  paramsBasedf.attr("row.names") = above.attr("row.names");
  
  //Allometries
  DataFrame paramsAllometriesdf = DataFrame::create(_["Hmax"] = Hmax,
                                                    _["Zmax"] = Zmax,
                                                    _["Aash"] = Aash, _["Absh"] = Absh, _["Bbsh"] = Bbsh,
                                                    _["r635"] = r635,
                                                    _["Acr"] = Acr, _["B1cr"] = B1cr, _["B2cr"] = B2cr, _["B3cr"] = B3cr,
                                                    _["C1cr"] = C1cr, _["C2cr"] = C2cr, 
                                                    _["Acw"] = Acw, _["Bcw"] = Bcw,
                                                    _["fHDmin"] = fHDmin,_["fHDmax"] = fHDmax);
  paramsAllometriesdf.attr("row.names") = above.attr("row.names");
  
  List input;
  if(transpirationMode=="Simple") {
    NumericVector WUESP = SpParams["WUE"];
    NumericVector WUE(numCohorts);
    NumericVector Psi_ExtractSP = SpParams["Psi_Extract"];
    NumericVector Psi_Extract(numCohorts);
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector pRootDisc(numCohorts);
    for(int c=0;c<numCohorts;c++){
      Psi_Extract[c]=Psi_ExtractSP[SP[c]];
      WUE[c]=WUESP[SP[c]];
      pRootDisc[c]=pRootDiscSP[SP[c]];
    }
    
    DataFrame paramsTranspdf = DataFrame::create(_["Psi_Extract"]=Psi_Extract,_["WUE"] = WUE, _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    
    List below = List::create( _["Z"]=Z,_["V"] = V);
    List paramsCanopy = List::create(_["gdd"] = 0);
    List numericParams = control["numericParams"];
    List paramsControl = List::create(_["verbose"] =control["verbose"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["cavitationRefill"] = control["cavitationRefill"],
                                      _["storagePool"] = storagePool);
    input = List::create(_["control"] = paramsControl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsTransp"] = paramsTranspdf,
                         _["paramsGrowth"]= paramsGrowthdf,
                         _["paramsAllometries"] = paramsAllometriesdf);
    
  } else if(transpirationMode =="Complex"){
    NumericVector GwminSP = SpParams["Gwmin"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector Gwmin(numCohorts), Gwmax(numCohorts);
    NumericVector xylem_kmax(numCohorts), Al2As(numCohorts);
    NumericVector VCstem_kmax(numCohorts);
    NumericVector VCstem_c(numCohorts), VCstem_d(numCohorts);
    NumericVector VCroot_c(numCohorts), VCroot_d(numCohorts);
    NumericVector Vmax298(numCohorts), Jmax298(numCohorts);
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector pRootDisc(numCohorts);
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
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmax[c], Al2As[c],H[c], control["taper"]); 
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      Gwmin[c] = GwminSP[SP[c]];
      Gwmax[c] = GwmaxSP[SP[c]];
      pRootDisc[c]=pRootDiscSP[SP[c]];
      double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      Vmax298[c] =Vmax298SP[SP[c]];
      Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c]));//Walker et al 2014
      for(int l=0;l<nlayers;l++) {
        // Rcout<<Vc[l]<<" ";
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                                                        VCstem_kmax[c], VCstem_c[c], VCstem_d[c]);
      }
      // Rcout<<"\n";
    }
    VGrhizo_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    VCroot_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    
    DataFrame paramsTranspdf = DataFrame::create(
      _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax, _["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,_["xylem_kmax"] = xylem_kmax,
      _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    
    List below = List::create( _["Z"]=Z,_["V"] = V,
                              _["VGrhizo_kmax"] = VGrhizo_kmax,
                              _["VCroot_kmax"] = VCroot_kmax);
    List numericParams = control["numericParams"];
    List paramsControl = List::create(_["verbose"] =control["verbose"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["canopyMode"] = canopyMode,
                                      _["verticalLayerSize"] = control["verticalLayerSize"],
                                      _["hydraulicCostFunction"] = control["hydraulicCostFunction"],
                                      _["taper"] = control["taper"],
                                      _["numericParams"] = clone(numericParams),
                                      _["ndailysteps"] = control["ndailysteps"], 
                                      _["cavitationRefill"] = control["cavitationRefill"],
                                      _["storagePool"] = storagePool);
    List paramsCanopy = List::create(_["gdd"] = 0,_["Temp"] = NA_REAL);
    input = List::create(_["control"] =paramsControl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                   _["below"] = below,
                   _["paramsBase"] = paramsBasedf,
                   _["paramsTransp"] = paramsTranspdf,
                   _["paramsGrowth"]= paramsGrowthdf,
                   _["paramsAllometries"] = paramsAllometriesdf);
  } 
  NumericVector tvec =  NumericVector(numCohorts, 0.0);
  tvec.attr("names") = above.attr("row.names");
  input["Transpiration"] = tvec;
  NumericVector pvec =  NumericVector(numCohorts, 0.0);
  pvec.attr("names") = above.attr("row.names");
  input["Photosynthesis"] = pvec;
  NumericVector cvec =  NumericVector(numCohorts, 0.0);
  cvec.attr("names") = above.attr("row.names");
  input["ProportionCavitated"] = cvec;
  // input["WindSpeed"] = NumericVector(numCohorts, 0.0);
  // input["PAR"] = NumericVector(numCohorts, 0.0);
  // input["AbsorbedSWR"] = NumericVector(numCohorts, 0.0);
  input.attr("class") = CharacterVector::create("growthInput","list");
  // df.attr("row.names") = seq(1,numCohorts);
  return(input);
}
// [[Rcpp::export("forest2growthInput")]]
List forest2growthInput(List x, List soil, DataFrame SpParams, List control) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  int nlayers = d.size();
  NumericVector Z(ntree+nshrub); //Rooting depth in cm
  NumericMatrix V(ntree+nshrub,nlayers);
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ = shrubData["Z"];  
  NumericVector Vi;
  for(int i=0;i<ntree;i++) {
    Vi = ldrRS_one(treeZ50[i], treeZ95[i],d);
    V(i,_) = Vi;
    Z[i] = treeZ95[i]/10.0;
  }
  for(int i=0;i<nshrub;i++) {
    Vi = conicRS_one(shrubZ[i],d);
    V(ntree+i,_) = Vi;
    Z[ntree+i] = shrubZ[i]/10.0; 
  }
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(growthInput(above,  Z, V, soil, SpParams, control));
}
