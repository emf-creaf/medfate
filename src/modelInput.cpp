#include <Rcpp.h>
#include "spwb.h"
#include "growth.h"
#include "root.h"
#include "forestutils.h"
#include "tissuemoisture.h"
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
// [[Rcpp::export("spwbInput")]]
List spwbInput(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  
  IntegerVector SP = above["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Simple") & (transpirationMode!="Complex")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Simple' or 'Complex')");

  int numStemSegments = control["nStemSegments"];
  double kver = control["ksymver"];
  
  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") & (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];

  NumericVector albedoSP = SpParams["albedo"];
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector albedo(numCohorts),k(numCohorts), g(numCohorts), Sgdd(numCohorts);
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
    k[c]=kSP[SP[c]];
    albedo[c]=albedoSP[SP[c]];
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
  DataFrame paramsBasedf = DataFrame::create(_["albedo"] = albedo, _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
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

    List paramsCanopy = List::create(_["gdd"] = 0);
    input = List::create(_["control"] = clone(control),
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsTransp"] = paramsTranspdf);
    
    NumericVector tvec =  NumericVector(numCohorts, 0.0);
    tvec.attr("names") = above.attr("row.names");
    input["Transpiration"] = tvec;
    NumericVector pvec =  NumericVector(numCohorts, 0.0);
    pvec.attr("names") = above.attr("row.names");
    input["Photosynthesis"] = pvec;
    input["PLC"] = NumericVector(numCohorts, 0.0);
  } else if(transpirationMode =="Complex"){
    NumericMatrix RWCstemmat =  NumericMatrix(numCohorts, numStemSegments);
    std::fill(RWCstemmat.begin(), RWCstemmat.end(), 1.0);
    RWCstemmat.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numStemSegments));
    NumericMatrix PLCmat =  NumericMatrix(numCohorts, numStemSegments);
    std::fill(PLCmat.begin(), PLCmat.end(), 0.0);
    PLCmat.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numStemSegments));
    NumericMatrix psiStemmat =  NumericMatrix(numCohorts, numStemSegments);
    std::fill(psiStemmat.begin(), psiStemmat.end(), 0.0);
    psiStemmat.attr("dimnames") = List::create(above.attr("row.names"), seq(1,numStemSegments));
    NumericVector Einst = NumericVector(numCohorts, 0.0);
    Einst.attr("names") = above.attr("row.names");
    NumericVector psiRoot = NumericVector(numCohorts, 0.0);
    psiRoot.attr("names") = above.attr("row.names");
    NumericVector psiLeaf = NumericVector(numCohorts, 0.0);
    psiLeaf.attr("names") = above.attr("row.names");
    NumericVector rwcsleaf = NumericVector(numCohorts, 1.0);
    rwcsleaf.attr("names") = above.attr("row.names");
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      warning("Soil pedotransfer functions set to Van Genuchten ('VG').");
    }
    CharacterVector GroupSP = SpParams["Group"];
    NumericVector leafwidthSP = SpParams["LeafWidth"];
    NumericVector HmedSP = SpParams["Hmed"]; //To correct conductivity
    NumericVector Al2AsSP = SpParams["Al2As"];
    NumericVector SLASP = SpParams["SLA"];
    NumericVector r635SP = SpParams["r635"];
    NumericVector LeafDensSP = SpParams["LeafDensity"];
    NumericVector WoodDensSP = SpParams["WoodDens"];
    NumericVector GwminSP = SpParams["Gwmin"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector VCleaf_kmaxSP = SpParams["VCleaf_kmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCleaf_cSP = SpParams["VCleaf_c"];
    NumericVector VCleaf_dSP = SpParams["VCleaf_d"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector LeafPI0SP = SpParams["LeafPI0"];
    NumericVector LeafEPSSP = SpParams["LeafEPS"];
    NumericVector LeafAFSP = SpParams["LeafAF"];
    NumericVector StemPI0SP = SpParams["StemPI0"];
    NumericVector StemEPSSP = SpParams["StemEPS"];
    NumericVector StemAFSP = SpParams["StemAF"];
    NumericVector rootxylem_kmaxSP = SpParams["rootxylem_kmax"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector Gwmax(numCohorts), Gwmin(numCohorts);
    NumericVector VCleaf_kmax(numCohorts), xylem_kmax(numCohorts), rootxylem_kmax(numCohorts), Al2As(numCohorts);
    NumericVector VCroottot_kmax(numCohorts, 0.0), VCstem_kmax(numCohorts),leafwidth(numCohorts);
    NumericVector pRootDisc(numCohorts);
    NumericVector VCleaf_c(numCohorts), VCleaf_d(numCohorts);
    NumericVector VCstem_c(numCohorts), VCstem_d(numCohorts);
    NumericVector VCroot_c(numCohorts), VCroot_d(numCohorts);
    NumericVector StemPI0(numCohorts), StemEPS(numCohorts), StemAF(numCohorts);
    NumericVector LeafPI0(numCohorts), LeafEPS(numCohorts), LeafAF(numCohorts);
    NumericVector ksymver(numCohorts);
    NumericVector Vsapwood(numCohorts), Vleaf(numCohorts);
    NumericVector SLA(numCohorts), LeafDens(numCohorts), WoodDens(numCohorts);
    NumericVector Vmax298(numCohorts), Jmax298(numCohorts);
    NumericVector r635(numCohorts);
    NumericVector dVec = soil["dVec"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VG_n = soil["VG_n"];
    int nlayers = dVec.size();
    NumericMatrix VCroot_kmax(numCohorts, nlayers); 
    NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
    NumericVector Vc;
    for(int c=0;c<numCohorts;c++){
      Vc = V(c,_);
      leafwidth[c] = leafwidthSP[SP[c]];
      xylem_kmax[c] = xylem_kmaxSP[SP[c]];
      rootxylem_kmax[c] = rootxylem_kmaxSP[SP[c]];
      if(NumericVector::is_na(rootxylem_kmax[c])) rootxylem_kmax[c] = xylem_kmax[c];
      Al2As[c] = Al2AsSP[SP[c]];
      SLA[c] = SLASP[SP[c]];
      WoodDens[c] = WoodDensSP[SP[c]];
      r635[c] = r635SP[SP[c]];
      LeafDens[c] = LeafDensSP[SP[c]];
      StemPI0[c] = StemPI0SP[SP[c]];
      StemEPS[c] = StemEPSSP[SP[c]];
      //From: Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N.M., Galbraith, D.R., Baker, T.R., Rowland, L., Fisher, R.A., Binks, O.J., Sevanto, S.A., Xu, C., Jansen, S., Choat, B., Mencuccini, M., McDowell, N.G., & Meir, P. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0: 1–60.
      if(NumericVector::is_na(StemPI0[c])) StemPI0[c] = 0.52 - 4.16*WoodDens[c]; 
      if(NumericVector::is_na(StemEPS[c])) StemEPS[c] = sqrt(1.02*exp(8.5*WoodDens[c])-2.89); 
      StemAF[c] = StemAFSP[SP[c]];
      LeafPI0[c] = LeafPI0SP[SP[c]];
      LeafEPS[c] = LeafEPSSP[SP[c]];
      LeafAF[c] = LeafAFSP[SP[c]];
      //Calculate stem and leaf capacity per leaf area (in m3·m-2)
      Vsapwood[c] = stemWaterCapacity(Al2As[c], H[c], WoodDens[c]); 
      Vleaf[c] = leafWaterCapacity(SLA[c], LeafDens[c]); 
      //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmax[c], HmedSP[SP[c]], Al2As[c],H[c], (GroupSP[SP[c]]=="Angiosperm"),control["taper"]); 
      ksymver[c] = kver/(H[c]/100.0);
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      //Default vulnerability curve parameters if missing
      if(NumericVector::is_na(VCroot_c[c])) VCroot_c[c] = VCstem_c[c];
      if(NumericVector::is_na(VCroot_d[c])) VCroot_d[c] = VCstem_d[c]/2.0;
      VCleaf_kmax[c] = VCleaf_kmaxSP[SP[c]];
      //Sack, L., & Holbrook, N.M. 2006. Leaf Hydraulics. Annual Review of Plant Biology 57: 361–381.
      if(NumericVector::is_na(VCleaf_kmax[c])) { 
        if(GroupSP[SP[c]]=="Angiosperm") {
          VCleaf_kmax[c] = 8.0;
        } else {
          VCleaf_kmax[c] = 6.0;
        }
      } 
      VCleaf_c[c]=VCleaf_cSP[SP[c]];
      VCleaf_d[c]=VCleaf_dSP[SP[c]];
      //Default vulnerability curve parameters if missing
      if(NumericVector::is_na(VCleaf_c[c])) VCleaf_c[c] = VCstem_c[c];
      if(NumericVector::is_na(VCleaf_d[c])) VCleaf_d[c] = VCstem_d[c]/1.5;
      pRootDisc[c]=pRootDiscSP[SP[c]];
      Gwmin[c] = GwminSP[SP[c]];
      //Duursma RA, Blackman CJ, Lopéz R, et al (2018) On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls. New Phytol. doi: 10.1111/nph.15395
      if(NumericVector::is_na(Gwmin[c])) Gwmin[c] = 0.0049;
      //Mencuccini M (2003) The ecological significance of long-distance water transport : short-term regulation , long-term acclimation and the hydraulic costs of stature across plant life forms. Plant Cell Environ 26:163–182
      Gwmax[c] = GwmaxSP[SP[c]];
      if(NumericVector::is_na(Gwmax[c])) Gwmax[c] = 0.12115*pow(VCleaf_kmax[c], 0.633);
      // double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      double VCroot_kmaxc = maximumRootHydraulicConductance(rootxylem_kmax[c],Al2As[c], Vc, dVec);
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      VCroottot_kmax[c] = sum(VCroot_kmax(c,_));
      Vmax298[c] =Vmax298SP[SP[c]];
      //Walker AP, Beckerman AP, Gu L, et al (2014) The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: A meta-analysis and modeling study. Ecol Evol 4:3218–3235. doi: 10.1002/ece3.1173
      Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c])); 
      for(int l=0;l<nlayers;l++) {
        
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                     VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                     VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c]);
        // Rcout<<VGrhizo_kmax(c,l)<<" ";
      }
      
      //Initialize levels of cavitation (to avoid initial water release)
      // double E = VCstem_kmax[c];
      // double psiUp = 0.0, psiDown = 0.0;
      // for(int i=0;i<PLCmat.ncol();i++) {
      //   psiDown = E2psiXylem(E,psiUp, (VCstem_kmax[c]*((double) PLCmat.ncol())), VCstem_c[c], VCstem_d[c]);
      //   PLCmat(c,i) = 1.0 - xylemConductance(psiDown, 1.0, VCstem_c[c], VCstem_d[c]); 
      //   psiUp = psiDown;
      // }
      // // psiDown = E2psiXylem(E,psiUp, VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c]);
      // // psiLeaf[c] = psiDown;
      // Rcout<<"\n";
    }
    VGrhizo_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    VCroot_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    
    DataFrame paramsAnatomydf = DataFrame::create(
       _["Al2As"] = Al2As, _["SLA"] = SLA, _["LeafWidth"] = leafwidth, 
       _["LeafDens"] = LeafDens, _["WoodDens"] = WoodDens, _["r635"] = r635
    );
    paramsAnatomydf.attr("row.names") = above.attr("row.names");
    DataFrame paramsWaterStoragedf = DataFrame::create(
      _["LeafPI0"] = LeafPI0, _["LeafEPS"] = LeafEPS, _["LeafAF"] = LeafAF, _["Vleaf"] = Vleaf,
      _["StemPI0"] = StemPI0, _["StemEPS"] = StemEPS, _["StemAF"] = StemAF, _["Vsapwood"] = Vsapwood);
    paramsWaterStoragedf.attr("row.names") = above.attr("row.names");
    DataFrame paramsTranspdf = DataFrame::create(
      _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax,_["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298, _["xylem_Kmax"] = xylem_kmax, _["root_Kmax"] = rootxylem_kmax,
      _["VCleaf_kmax"]=VCleaf_kmax,_["VCleaf_c"]=VCleaf_c,_["VCleaf_d"]=VCleaf_d,
      _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, 
      _["VCroot_kmax"] = VCroottot_kmax ,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,
      _["ksymver"] = ksymver,
      _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    List below = List::create(_["V"] = V,
                              _["VGrhizo_kmax"] = VGrhizo_kmax,
                              _["VCroot_kmax"] = VCroot_kmax);
    List paramsCanopy = List::create(_["gdd"] = 0,_["Temp"] = NA_REAL);
    input = List::create(_["control"] = clone(control),
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsTransp"] = paramsTranspdf,
                         _["paramsWaterStorage"] = paramsWaterStoragedf);
    
    NumericVector tvec =  NumericVector(numCohorts, 0.0);
    tvec.attr("names") = above.attr("row.names");
    input["Transpiration"] = tvec;
    NumericVector pvec =  NumericVector(numCohorts, 0.0);
    pvec.attr("names") = above.attr("row.names");
    input["Photosynthesis"] = pvec;
    input["PLCstem"] = PLCmat;
    input["RWCsympstem"] = RWCstemmat;
    input["RWCsympleaf"] = rwcsleaf;
    input["Einst"] = Einst;
    input["psiRoot"] = psiRoot;
    input["psiStem"] = psiStemmat;
    input["psiLeaf"] = psiLeaf;
  }

  // input["WindSpeed"] = NumericVector(numCohorts, 0.0);
  // input["PAR"] = NumericVector(numCohorts, 0.0);
  // input["AbsorbedSWR"] = NumericVector(numCohorts, 0.0);
  input.attr("class") = CharacterVector::create("spwbInput","list");
  // df.attr("row.names") = seq(1,numCohorts);
  return(input);
}

// [[Rcpp::export("forest2spwbInput")]]
List forest2spwbInput(List x, List soil, DataFrame SpParams, List control) {
  NumericMatrix V = forest2belowground(x,soil, SpParams);
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(spwbInput(above,  V, soil, SpParams, control));
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

  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") & (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  
  String storagePool = control["storagePool"];
  if((storagePool!="none") & (storagePool!="one")& (storagePool!="two")) stop("Wrong storage pool ('storagePool' should be 'none', 'one' or 'two')");

  
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  
  NumericVector Al2AsSP = SpParams["Al2As"];
  NumericVector SLASP = SpParams["SLA"];
  NumericVector WoodCSP = SpParams["WoodC"];
  NumericVector WoodDensSP = SpParams["WoodDens"];
  NumericVector RGRmaxSP = SpParams["RGRmax"];
  NumericVector CstoragepmaxSP = SpParams["Cstoragepmax"];
  
  NumericVector albedoSP = SpParams["albedo"];
  NumericVector kSP = SpParams["k"];
  NumericVector gSP = SpParams["g"];
  NumericVector SgddSP = SpParams["Sgdd"];
  int numCohorts = SP.size();
  NumericVector albedo(numCohorts), k(numCohorts), g(numCohorts), Sgdd(numCohorts);
  NumericVector SLA(numCohorts), Al2As(numCohorts),  WoodC(numCohorts), WoodDens(numCohorts);
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
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
    albedo[c]=albedoSP[SP[c]];
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
  DataFrame paramsBasedf = DataFrame::create(_["albedo"] = albedo, _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
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
                                      _["soilFunctions"] =soilFunctions, 
                                      _["snowpack"] = control["snowpack"],
                                      _["drainage"] = control["drainage"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["cavitationRefill"] = control["cavitationRefill"],
                                      _["defaultWindSpeed"] = control["defaultWindSpeed"],
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
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      warning("Soil pedotransfer functions set to Van Genuchten ('VG').");
    }
    CharacterVector GroupSP = SpParams["Group"];
    NumericVector HmedSP = SpParams["Hmed"]; //To correct conductivity
    NumericVector leafwidthSP = SpParams["LeafWidth"];
    NumericVector GwminSP = SpParams["Gwmin"];
    NumericVector GwmaxSP = SpParams["Gwmax"];
    NumericVector VCleaf_kmaxSP = SpParams["VCleaf_kmax"];
    NumericVector xylem_kmaxSP = SpParams["xylem_kmax"];
    NumericVector VCleaf_cSP = SpParams["VCleaf_c"];
    NumericVector VCleaf_dSP = SpParams["VCleaf_d"];
    NumericVector VCstem_cSP = SpParams["VCstem_c"];
    NumericVector VCstem_dSP = SpParams["VCstem_d"];
    NumericVector rootxylem_kmaxSP = SpParams["rootxylem_kmax"];
    NumericVector VCroot_cSP = SpParams["VCroot_c"];
    NumericVector VCroot_dSP = SpParams["VCroot_d"];
    NumericVector Vmax298SP = SpParams["Vmax298"];
    NumericVector Gwmin(numCohorts), Gwmax(numCohorts);
    NumericVector VCleaf_kmax(numCohorts), xylem_kmax(numCohorts), rootxylem_kmax(numCohorts),Al2As(numCohorts);
    NumericVector VCroottot_kmax(numCohorts, 0.0), VCstem_kmax(numCohorts),leafwidth(numCohorts);
    NumericVector VCstem_c(numCohorts), VCstem_d(numCohorts);
    NumericVector VCroot_c(numCohorts), VCroot_d(numCohorts);
    NumericVector Vmax298(numCohorts), Jmax298(numCohorts);
    NumericVector pRootDiscSP = SpParams["pRootDisc"];
    NumericVector pRootDisc(numCohorts);
    NumericVector dVec = soil["dVec"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VG_n = soil["VG_n"];
    int nlayers = dVec.size();
    NumericVector VCleaf_c(numCohorts), VCleaf_d(numCohorts);
    NumericMatrix VCroot_kmax(numCohorts, nlayers); 
    NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
    NumericVector Vc;
    for(int c=0;c<numCohorts;c++){
      Vc = V(c,_);
      xylem_kmax[c] = xylem_kmaxSP[SP[c]];
      rootxylem_kmax[c] = rootxylem_kmaxSP[SP[c]];
      if(NumericVector::is_na(rootxylem_kmax[c])) rootxylem_kmax[c] = xylem_kmax[c];
      leafwidth[c]=leafwidthSP[SP[c]];
      Al2As[c] = Al2AsSP[SP[c]];
      //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
      VCstem_kmax[c]=maximumStemHydraulicConductance(xylem_kmax[c], HmedSP[SP[c]], Al2As[c],H[c], GroupSP[SP[c]]=="Angiosperm", control["taper"]); 
      VCstem_c[c]=VCstem_cSP[SP[c]];
      VCstem_d[c]=VCstem_dSP[SP[c]];
      VCroot_c[c]=VCroot_cSP[SP[c]];
      VCroot_d[c]=VCroot_dSP[SP[c]];
      //Default vulnerability curve parameters if missing
      if(NumericVector::is_na(VCroot_c[c])) VCroot_c[c] = VCstem_c[c];
      if(NumericVector::is_na(VCroot_d[c])) VCroot_d[c] = VCstem_d[c]/2.0;
      VCleaf_kmax[c] = VCleaf_kmaxSP[SP[c]];
      if(NumericVector::is_na(VCleaf_kmax[c])) { //Sack, L., & Holbrook, N.M. 2006. Leaf Hydraulics. Annual Review of Plant Biology 57: 361–381.
        if(GroupSP[SP[c]]=="Angiosperm") {
          VCleaf_kmax[c] = 8.0;
        } else {
          VCleaf_kmax[c] = 5.0;
        }
      } 
      VCleaf_c[c]=VCleaf_cSP[SP[c]];
      VCleaf_d[c]=VCleaf_dSP[SP[c]];
      //Default vulnerability curve parameters if missing
      if(NumericVector::is_na(VCleaf_c[c])) VCleaf_c[c] = VCstem_c[c];
      if(NumericVector::is_na(VCleaf_d[c])) VCleaf_d[c] = VCstem_d[c]/1.5;
      Gwmin[c] = GwminSP[SP[c]];
      Gwmax[c] = GwmaxSP[SP[c]];
      pRootDisc[c]=pRootDiscSP[SP[c]];
      // double VCroot_kmaxc = 1.0/((1.0/(VCstem_kmax[c]*fracTotalTreeResistance))-(1.0/VCstem_kmax[c]));
      double VCroot_kmaxc = maximumRootHydraulicConductance(rootxylem_kmax[c],Al2As[c], Vc, dVec);
      VCroot_kmax(c,_) = VCroot_kmaxc*xylemConductanceProportions(Vc,dVec);
      VCroottot_kmax[c] = sum(VCroot_kmax(c,_));
      Vmax298[c] =Vmax298SP[SP[c]];
      Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c]));//Walker et al 2014
      for(int l=0;l<nlayers;l++) {
        // Rcout<<Vc[l]<<" ";
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroot_kmaxc, VCroot_c[c], VCroot_d[c],
                     VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                     VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c]);
      }
      // Rcout<<"\n";
    }
    VGrhizo_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    VCroot_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
    
    DataFrame paramsTranspdf = DataFrame::create(
        _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax, _["LeafWidth"] = leafwidth, _["Vmax298"]=Vmax298,
        _["Jmax298"]=Jmax298,_["xylem_Kmax"] = xylem_kmax, _["root_Kmax"] = rootxylem_kmax,_["Al2As"] = Al2As,  
        _["VCleaf_kmax"]=VCleaf_kmax,_["VCleaf_c"]=VCleaf_c,_["VCleaf_d"]=VCleaf_d,
        _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, 
        _["VCroot_kmax"] = VCroottot_kmax ,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,
        _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    
    List below = List::create( _["Z"]=Z,_["V"] = V,
                              _["VGrhizo_kmax"] = VGrhizo_kmax,
                              _["VCroot_kmax"] = VCroot_kmax);
    List numericParams = control["numericParams"];
    List paramsControl = List::create(_["verbose"] =control["verbose"],
                                      _["transpirationMode"] =transpirationMode, 
                                      _["soilFunctions"] =soilFunctions, 
                                      _["snowpack"] = control["snowpack"],
                                      _["drainage"] = control["drainage"],
                                      _["Catm"] = control["Catm"],                                      
                                      _["averageFracRhizosphereResistance"] = control["averageFracRhizosphereResistance"],
                                      _["verticalLayerSize"] = control["verticalLayerSize"],
                                      _["hydraulicCostFunction"] = control["hydraulicCostFunction"],
                                      _["taper"] = control["taper"],
                                      _["numericParams"] = clone(numericParams),
                                      _["ndailysteps"] = control["ndailysteps"], 
                                      _["cavitationRefill"] = control["cavitationRefill"],
                                      _["thermalCapacityLAI"] = control["thermalCapacityLAI"],
                                      _["defaultWindSpeed"] = control["defaultWindSpeed"],
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
  
  NumericVector Z(ntree+nshrub); //Rooting depth in cm
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ = shrubData["Z"];  
  NumericMatrix V = forest2belowground(x,soil, SpParams);
  for(int i=0;i<ntree;i++) {
    Z[i] = treeZ95[i]/10.0;
  }
  for(int i=0;i<nshrub;i++) {
    Z[ntree+i] = shrubZ[i]/10.0; 
  }
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(growthInput(above,  Z, V, soil, SpParams, control));
}
