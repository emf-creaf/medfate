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

DataFrame paramsAnatomy(DataFrame above, DataFrame SpParams) {
  IntegerVector SP = above["SP"];
  int numCohorts = SP.size();
  NumericVector Al2As = cohortNumericParameter(SP, SpParams, "Al2As");
  NumericVector SLA = cohortNumericParameter(SP, SpParams, "SLA");
  NumericVector LeafDensity = cohortNumericParameter(SP, SpParams, "LeafDensity");
  NumericVector WoodDensity = cohortNumericParameter(SP, SpParams, "WoodDensity");
  NumericVector r635 = cohortNumericParameter(SP, SpParams, "r635");
  NumericVector leafwidth = cohortNumericParameter(SP, SpParams, "LeafWidth");
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(Al2As[c])) Al2As[c] = 2500.0; // = 4 cm2·m-2
  }
  DataFrame paramsAnatomydf = DataFrame::create(
    _["Al2As"] = Al2As, _["SLA"] = SLA, _["LeafWidth"] = leafwidth, 
    _["LeafDensity"] = LeafDensity, _["WoodDensity"] = WoodDensity, _["r635"] = r635
  );
  paramsAnatomydf.attr("row.names") = above.attr("row.names");
  return(paramsAnatomydf);
}

DataFrame paramsWaterStorage(DataFrame above, DataFrame SpParams,
                             DataFrame paramsAnatomydf) {
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  int numCohorts = SP.size();
  
  NumericVector LeafPI0 = cohortNumericParameter(SP, SpParams, "LeafPI0");
  NumericVector LeafEPS = cohortNumericParameter(SP, SpParams, "LeafEPS");
  NumericVector LeafAF = cohortNumericParameter(SP, SpParams, "LeafAF");
  NumericVector StemPI0 = cohortNumericParameter(SP, SpParams, "StemPI0");
  NumericVector StemEPS = cohortNumericParameter(SP, SpParams, "StemEPS");
  NumericVector StemAF = cohortNumericParameter(SP, SpParams, "StemAF");
  
  NumericVector Vsapwood(numCohorts), Vleaf(numCohorts);
  
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector LeafDensity = paramsAnatomydf["LeafDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  
  for(int c=0;c<numCohorts;c++){
    //From: Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N.M., Galbraith, D.R., Baker, T.R., Rowland, L., Fisher, R.A., Binks, O.J., Sevanto, S.A., Xu, C., Jansen, S., Choat, B., Mencuccini, M., McDowell, N.G., & Meir, P. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0: 1–60.
    if(NumericVector::is_na(StemPI0[c])) StemPI0[c] = 0.52 - 4.16*WoodDensity[c]; 
    if(NumericVector::is_na(StemEPS[c])) StemEPS[c] = sqrt(1.02*exp(8.5*WoodDensity[c])-2.89); 
    if(NumericVector::is_na(StemAF[c])) StemAF[c] = 0.8;
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(LeafPI0[c])) LeafPI0[c] = -2.0; //Average values for Mediterranean climate species
    if(NumericVector::is_na(LeafEPS[c])) LeafEPS[c] = 17.0;
    if(NumericVector::is_na(LeafAF[c])) LeafAF[c] = 0.29;
    
    //Calculate stem and leaf capacity per leaf area (in m3·m-2)
    Vsapwood[c] = stemWaterCapacity(Al2As[c], H[c], WoodDensity[c]); 
    Vleaf[c] = leafWaterCapacity(SLA[c], LeafDensity[c]); 
  }
  DataFrame paramsWaterStoragedf = DataFrame::create(
    _["LeafPI0"] = LeafPI0, _["LeafEPS"] = LeafEPS, _["LeafAF"] = LeafAF, _["Vleaf"] = Vleaf,
      _["StemPI0"] = StemPI0, _["StemEPS"] = StemEPS, _["StemAF"] = StemAF, _["Vsapwood"] = Vsapwood);
  paramsWaterStoragedf.attr("row.names") = above.attr("row.names");
  
  return(paramsWaterStoragedf);
}

DataFrame paramsTranspiration(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, 
                              DataFrame paramsAnatomydf, List control) {
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  int numCohorts = SP.size();
  NumericVector Vc;
  double fracRootResistance = control["fracRootResistance"];
  
  NumericVector dVec = soil["dVec"];
  
  CharacterVector Group = cohortCharacterParameter(SP, SpParams, "Group");
  CharacterVector Order = cohortCharacterParameter(SP, SpParams, "Order");
  CharacterVector TreeType = cohortCharacterParameter(SP, SpParams, "TreeType");
  
  NumericVector Hmed = cohortNumericParameter(SP, SpParams, "Hmed"); //To correct conductivity
  NumericVector Gwmin = cohortNumericParameter(SP, SpParams, "Gwmin");
  NumericVector Gwmax = cohortNumericParameter(SP, SpParams, "Gwmax");
  NumericVector VCleaf_kmax = cohortNumericParameter(SP, SpParams, "VCleaf_kmax");
  NumericVector Kmax_stemxylem = cohortNumericParameter(SP, SpParams, "Kmax_stemxylem");
  NumericVector VCleaf_c = cohortNumericParameter(SP, SpParams, "VCleaf_c");
  NumericVector VCleaf_d = cohortNumericParameter(SP, SpParams, "VCleaf_d");
  NumericVector VCstem_c = cohortNumericParameter(SP, SpParams, "VCstem_c");
  NumericVector VCstem_d = cohortNumericParameter(SP, SpParams, "VCstem_d");
  NumericVector Kmax_rootxylem = cohortNumericParameter(SP, SpParams, "Kmax_rootxylem");
  NumericVector VCroot_c = cohortNumericParameter(SP, SpParams, "VCroot_c");
  NumericVector VCroot_d = cohortNumericParameter(SP, SpParams, "VCroot_d");
  NumericVector SLA = cohortNumericParameter(SP, SpParams, "SLA");
  NumericVector Narea = cohortNumericParameter(SP, SpParams, "Narea");
  NumericVector Vmax298 = cohortNumericParameter(SP, SpParams, "Vmax298");
  NumericVector Jmax298 = cohortNumericParameter(SP, SpParams, "Jmax298");
  NumericVector pRootDisc = cohortNumericParameter(SP, SpParams, "pRootDisc");
  
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  
  NumericVector VCstem_kmax(numCohorts);
  NumericVector VCroottot_kmax(numCohorts, 0.0);
  NumericVector Plant_kmax(numCohorts, 0.0);
  
  
  for(int c=0;c<numCohorts;c++){
    Vc = V(c,_);
    
    if(NumericVector::is_na(Kmax_stemxylem[c])) {
      // From: Maherali H, Pockman W, Jackson R (2004) Adaptive variation in the vulnerability of woody plants to xylem cavitation. Ecology 85:2184–2199
      if(Group[c]=="Angiosperm") {
        if(TreeType[c]=="Shrub") {
          Kmax_stemxylem[c] = 1.55; //Angiosperm deciduous shrub
        } else if(TreeType[c]=="Deciduous") {
          Kmax_stemxylem[c] = 1.58; //Angiosperm winter-deciduous tree
        } else { 
          Kmax_stemxylem[c] = 2.43; //Angiosperm evergreen tree
        }
      } else {
        if(TreeType[c]=="Shrub") {
          Kmax_stemxylem[c] = 0.24; //Gymnosperm shrub
        } else {
          Kmax_stemxylem[c] = 0.48; //Gymnosperm tree
        }
      }
    }
    //Oliveras I, Martínez-Vilalta J, Jimenez-Ortiz T, et al (2003) Hydraulic architecture of Pinus halepensis, P . pinea and Tetraclinis articulata in a dune ecosystem of Eastern Spain. Plant Ecol 131–141
    if(NumericVector::is_na(Kmax_rootxylem[c])) Kmax_rootxylem[c] = 4.0*Kmax_stemxylem[c];
    
    //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
    VCstem_kmax[c]=maximumStemHydraulicConductance(Kmax_stemxylem[c], Hmed[c], Al2As[c],H[c], (Group[c]=="Angiosperm"),control["taper"]); 
    
    //Xylem vulnerability curve
    if(NumericVector::is_na(VCstem_d[c]) | NumericVector::is_na(VCstem_c[c])) {
      double psi50 = NA_REAL;
      // From: Maherali H, Pockman W, Jackson R (2004) Adaptive variation in the vulnerability of woody plants to xylem cavitation. Ecology 85:2184–2199
      if(Group[c]=="Angiosperm") {
        if(TreeType[c]=="Shrub") {
          psi50 = -5.09; //Angiosperm evergreen shrub
        } else if(TreeType[c]=="Deciduous") {
          psi50 = -2.34; //Angiosperm winter-deciduous tree
        } else { 
          psi50 = -1.51; //Angiosperm evergreen tree
        }
      } else {
        if(TreeType[c]=="Shrub") {
          psi50 = -8.95; //Gymnosperm shrub
        } else {
          psi50 = -4.17; //Gymnosperm tree
        }
      }
      double psi88 = 1.2593*psi50 - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50, psi88);
      if(NumericVector::is_na(VCstem_c[c])) VCstem_c[c] = par["c"];
      if(NumericVector::is_na(VCstem_d[c])) VCstem_d[c] = par["d"];
    }
    
    //Default vulnerability curve parameters if missing
    if(NumericVector::is_na(VCroot_d[c]) | NumericVector::is_na(VCroot_c[c])) {
      double psi50stem = VCstem_d[c]*pow(0.6931472,1.0/VCstem_c[c]);
      double psi50root = 0.742*psi50stem + 0.4892; //Regression using data from Bartlett et al. 2016
      double psi88root = 1.2593*psi50root - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50root, psi88root);
      if(NumericVector::is_na(VCroot_c[c])) VCroot_c[c] = par["c"];
      if(NumericVector::is_na(VCroot_d[c])) VCroot_d[c] = par["d"];
    }
    //Sack, L., & Holbrook, N.M. 2006. Leaf Hydraulics. Annual Review of Plant Biology 57: 361–381.
    if(NumericVector::is_na(VCleaf_kmax[c])) { 
      if(Group[c]=="Angiosperm") {
        VCleaf_kmax[c] = 8.0;
      } else {
        VCleaf_kmax[c] = 6.0;
      }
    } 
    //Default vulnerability curve parameters if missing
    if(NumericVector::is_na(VCleaf_c[c])) VCleaf_c[c] = VCstem_c[c];
    if(NumericVector::is_na(VCleaf_d[c])) VCleaf_d[c] = VCstem_d[c]/1.5;
    //Duursma RA, Blackman CJ, Lopéz R, et al (2018) On the minimum leaf conductance: its role in models of plant water use, and ecological and environmental controls. New Phytol. doi: 10.1111/nph.15395
    if(NumericVector::is_na(Gwmin[c])) {
      if(Order[c]=="Pinales") {
        Gwmin[c] = 0.003;
      } else if(Order[c]=="Araucariales") {
        Gwmin[c] = 0.003;
      } else if(Order[c]=="Ericales") {
        Gwmin[c] = 0.004;
      } else if(Order[c]=="Fagales") {
        Gwmin[c] = 0.0045;
      } else if(Order[c]=="Rosales") {
        Gwmin[c] = 0.0045;
      } else if(Order[c]=="Cupressales") {
        Gwmin[c] = 0.0045;
      } else if(Order[c]=="Lamiales") {
        Gwmin[c] = 0.0055;
      } else if(Order[c]=="Fabales") {
        Gwmin[c] = 0.0065;
      } else if(Order[c]=="Myrtales") {
        Gwmin[c] = 0.0075;
      } else if(Order[c]=="Poales") {
        Gwmin[c] = 0.0110;
      } else {
        Gwmin[c] = 0.0049;
      }
    }
    //Mencuccini M (2003) The ecological significance of long-distance water transport : short-term regulation , long-term acclimation and the hydraulic costs of stature across plant life forms. Plant Cell Environ 26:163–182
    if(NumericVector::is_na(Gwmax[c])) Gwmax[c] = 0.12115*pow(VCleaf_kmax[c], 0.633);
    
    double VCroot_kmaxc = NA_REAL;
    if(NumericVector::is_na(fracRootResistance)) {
      VCroot_kmaxc = maximumRootHydraulicConductance(Kmax_rootxylem[c],Al2As[c], Vc, dVec);
    } else {
      double rstem = (1.0/VCstem_kmax[c]);
      double rleaf = (1.0/VCleaf_kmax[c]);
      double rtot = (rstem+rleaf)/(1.0 - fracRootResistance);
      VCroot_kmaxc = 1.0/(rtot - rstem - rleaf);
    }
    VCroottot_kmax[c] = VCroot_kmaxc;
    
    if(NumericVector::is_na(Vmax298[c]))  {
      if(!NumericVector::is_na(SLA[c]) & !NumericVector::is_na(Narea[c]))  {
        //Walker AP, Beckerman AP, Gu L, et al (2014) The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: A meta-analysis and modeling study. Ecol Evol 4:3218–3235. doi: 10.1002/ece3.1173
        double lnN = log(Narea[c]);
        double lnSLA = log(SLA[c]/1000.0); //SLA in m2*g-1
        Vmax298[c] = exp(1.993 + 2.555*lnN - 0.372*lnSLA + 0.422*lnN*lnSLA);
      } else {
        Vmax298[c] = 100.0; //Default value of Vmax298 = 100.0
      }
    }
    //Walker AP, Beckerman AP, Gu L, et al (2014) The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: A meta-analysis and modeling study. Ecol Evol 4:3218–3235. doi: 10.1002/ece3.1173
    if(NumericVector::is_na(Jmax298[c])) Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c])); 
    
    //Plant kmax
    Plant_kmax[c] = 1.0/((1.0/VCleaf_kmax[c])+(1.0/VCstem_kmax[c])+(1.0/VCroottot_kmax[c]));
  }
  
  DataFrame paramsTranspdf = DataFrame::create(
    _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax,_["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298, _["Kmax_stemxylem"] = Kmax_stemxylem, _["Kmax_rootxylem"] = Kmax_rootxylem,
        _["VCleaf_kmax"]=VCleaf_kmax,_["VCleaf_c"]=VCleaf_c,_["VCleaf_d"]=VCleaf_d,
        _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, 
        _["VCroot_kmax"] = VCroottot_kmax ,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,
        _["pRootDisc"] = pRootDisc,
        _["Plant_kmax"] = Plant_kmax);
  paramsTranspdf.attr("row.names") = above.attr("row.names");
  return(paramsTranspdf);
}

List paramsBelow(DataFrame above, NumericMatrix V, List soil, 
                 DataFrame paramsTranspirationdf, List control) {
  NumericVector dVec = soil["dVec"];
  NumericVector VG_alpha = soil["VG_alpha"];
  NumericVector VG_n = soil["VG_n"];
  int nlayers = dVec.size();
  
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  
  //Soil layer names
  CharacterVector slnames(V.ncol());
  for(int i=0;i<V.ncol();i++) slnames[i] = i+1;
  V.attr("dimnames") = List::create(above.attr("row.names"), slnames);
  
  NumericVector VCroottot_kmax = paramsTranspirationdf["VCroot_kmax"];
  NumericVector VCroot_c = paramsTranspirationdf["VCroot_c"];
  NumericVector VCroot_d = paramsTranspirationdf["VCroot_d"];
  NumericVector VCstem_kmax = paramsTranspirationdf["VCstem_kmax"];
  NumericVector VCstem_c = paramsTranspirationdf["VCstem_c"];
  NumericVector VCstem_d = paramsTranspirationdf["VCstem_d"];
  NumericVector VCleaf_kmax = paramsTranspirationdf["VCleaf_kmax"];
  NumericVector VCleaf_c = paramsTranspirationdf["VCleaf_c"];
  NumericVector VCleaf_d = paramsTranspirationdf["VCleaf_d"];
  
  int numCohorts = VCroottot_kmax.length();
  
  NumericMatrix VCroot_kmax(numCohorts, nlayers); 
  NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
  NumericVector Vc;
  for(int c=0;c<numCohorts;c++){
    Vc = V(c,_);
    
    for(int l=0;l<nlayers;l++) {
      VCroot_kmax(c,_) = VCroottot_kmax[c]*xylemConductanceProportions(Vc,dVec);
      VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                   VCroot_kmax[c], VCroot_c[c], VCroot_d[c],
                   VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                   VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c]);
    }
  }
  VGrhizo_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
  VCroot_kmax.attr("dimnames") = List::create(above.attr("row.names"), slnames);
  
  
  List below = List::create(_["V"] = V,
                            _["VGrhizo_kmax"] = VGrhizo_kmax,
                            _["VCroot_kmax"] = VCroot_kmax);
  return(below);
}

List paramsBelowZ(DataFrame above, NumericMatrix V, NumericVector Z, List soil, 
                 DataFrame paramsTranspirationdf, List control) {
  
  List belowTemp = paramsBelow(above, V, soil, paramsTranspirationdf, control);

  Z.attr("names") = above.attr("row.names");
  
  List below = List::create(_["V"] = V, _["Z"] = Z,
                            _["VGrhizo_kmax"] = belowTemp["VGrhizo_kmax"],
                            _["VCroot_kmax"] = belowTemp["VCroot_kmax"]);
  
  return(below);
}

DataFrame paramsGrowth(DataFrame above, DataFrame SpParams) {
  IntegerVector SP = above["SP"];
  int numCohorts = SP.size();
  
  NumericVector WoodC = cohortNumericParameter(SP, SpParams, "WoodC");
  NumericVector RGRmax = cohortNumericParameter(SP, SpParams, "RGRmax");
  NumericVector Cstoragepmax = cohortNumericParameter(SP, SpParams, "Cstoragepmax");
  

  for(int c=0;c<numCohorts;c++){
    Cstoragepmax[c] = std::max(0.05,Cstoragepmax[c]); //Minimum 5%
  }
  
  DataFrame paramsGrowthdf = DataFrame::create(_["WoodC"] = WoodC, _["Cstoragepmax"] = Cstoragepmax, _["RGRmax"] = RGRmax);
  paramsGrowthdf.attr("row.names") = above.attr("row.names");
  return(paramsGrowthdf);
}


DataFrame paramsAllometries(DataFrame above, DataFrame SpParams) {
  IntegerVector SP = above["SP"];
  
  NumericVector Hmax = cohortNumericParameter(SP, SpParams, "Hmax");
  NumericVector Zmax = cohortNumericParameter(SP, SpParams, "Zmax");
  NumericVector r635 = cohortNumericParameter(SP, SpParams, "r635");
  NumericVector Aash = cohortNumericParameter(SP, SpParams, "a_ash");
  NumericVector Absh = cohortNumericParameter(SP, SpParams, "a_bsh");
  NumericVector Bbsh = cohortNumericParameter(SP, SpParams, "b_bsh");
  NumericVector Acr = cohortNumericParameter(SP, SpParams, "a_cr");
  NumericVector B1cr = cohortNumericParameter(SP, SpParams, "b_1cr");
  NumericVector B2cr = cohortNumericParameter(SP, SpParams, "b_2cr");
  NumericVector B3cr = cohortNumericParameter(SP, SpParams, "b_3cr");
  NumericVector C1cr = cohortNumericParameter(SP, SpParams, "c_1cr");
  NumericVector C2cr = cohortNumericParameter(SP, SpParams, "c_2cr");
  NumericVector Acw = cohortNumericParameter(SP, SpParams, "a_cw");
  NumericVector Bcw = cohortNumericParameter(SP, SpParams, "b_cw");
  NumericVector fHDmin = cohortNumericParameter(SP, SpParams, "fHDmin");
  NumericVector fHDmax = cohortNumericParameter(SP, SpParams, "fHDmax");

  DataFrame paramsAllometriesdf = DataFrame::create(_["Hmax"] = Hmax,
                                                    _["Zmax"] = Zmax,
                                                    _["Aash"] = Aash, _["Absh"] = Absh, _["Bbsh"] = Bbsh,
                                                    _["r635"] = r635,
                                                    _["Acr"] = Acr, _["B1cr"] = B1cr, _["B2cr"] = B2cr, _["B3cr"] = B3cr,
                                                    _["C1cr"] = C1cr, _["C2cr"] = C2cr, 
                                                    _["Acw"] = Acw, _["Bcw"] = Bcw,
                                                    _["fHDmin"] = fHDmin,_["fHDmax"] = fHDmax);
  paramsAllometriesdf.attr("row.names") = above.attr("row.names");
  return(paramsAllometriesdf);
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
  if((transpirationMode!="Granier") & (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Granier' or 'Sperry')");

  
  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") & (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  

  NumericVector W = soil["W"];
  int nlayers = W.length();
  NumericVector albedo = cohortNumericParameter(SP, SpParams, "albedo");
  NumericVector k = cohortNumericParameter(SP, SpParams, "k");
  NumericVector g = cohortNumericParameter(SP, SpParams, "g");
  NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  

  //Cohort description
  CharacterVector nsp = cohortCharacterParameter(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  //Above 
  DataFrame plantsdf = DataFrame::create(_["H"]=H, _["CR"]=CR, 
                                         _["LAI_live"]=LAI_live, 
                                         _["LAI_expanded"] = LAI_expanded, 
                                         _["LAI_dead"] = LAI_dead);
  plantsdf.attr("row.names") = above.attr("row.names");
  
  //Base params
  DataFrame paramsBasedf = DataFrame::create(_["albedo"] = albedo, _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  paramsBasedf.attr("row.names") = above.attr("row.names");
  
 

  List input;
  if(transpirationMode=="Granier") {
    NumericVector WUE = cohortNumericParameter(SP, SpParams, "WUE");
    NumericVector Psi_Extract = cohortNumericParameter(SP, SpParams, "Psi_Extract");
    NumericVector pRootDisc = cohortNumericParameter(SP, SpParams, "pRootDisc");
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
  } else if(transpirationMode =="Sperry"){
    int numStemSegments = control["nStemSegments"];
    bool capacitance = control["capacitance"];
    
    DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams);
    DataFrame paramsWaterStoragedf = paramsWaterStorage(above, SpParams, paramsAnatomydf);
    DataFrame paramsTranspirationdf = paramsTranspiration(above, V, soil, SpParams,
                                                          paramsAnatomydf, control);
    List below = paramsBelow(above, V, soil, 
                             paramsTranspirationdf, control);
    
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
    NumericMatrix psiRhizo =  NumericMatrix(numCohorts, nlayers);
    psiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
    std::fill(psiRhizo.begin(), psiRhizo.end(), 0.0);
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
    List paramsCanopy = List::create(_["gdd"] = 0,_["Temp"] = NA_REAL);
    List ctl = clone(control);
    if(capacitance) {
      ctl["hydraulicCostFunction"] = 2;
      warning("Hydraulic cost function set to '2'.");
    }
    input = List::create(_["control"] = ctl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsTransp"] = paramsTranspirationdf,
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
    input["psiRhizo"] = psiRhizo;
    input["psiRoot"] = psiRoot;
    input["psiStem"] = psiStemmat;
    input["psiLeaf"] = psiLeaf;
  }

  input.attr("class") = CharacterVector::create("spwbInput","list");
  return(input);
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
  if((transpirationMode!="Granier") & (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Granier' or 'Sperry')");

  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") & (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  
  String storagePool = control["storagePool"];
  if((storagePool!="none") & (storagePool!="one")& (storagePool!="two")) stop("Wrong storage pool ('storagePool' should be 'none', 'one' or 'two')");

  
  NumericVector W = soil["W"];
  int nlayers = W.length();
  
  DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams);
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector SLA = paramsAnatomydf["Al2As"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  
  DataFrame paramsGrowthdf = paramsGrowth(above, SpParams);
  NumericVector WoodC = paramsGrowthdf["WoodC"];
  NumericVector Cstoragepmax = paramsGrowthdf["Cstoragepmax"];
  
  
  DataFrame paramsAllometriesdf = paramsAllometries(above, SpParams);
  
  NumericVector albedo = cohortNumericParameter(SP, SpParams, "albedo");
  NumericVector k = cohortNumericParameter(SP, SpParams, "k");
  NumericVector g = cohortNumericParameter(SP, SpParams, "g");
  NumericVector Sgdd = cohortNumericParameter(SP, SpParams, "Sgdd");
  
  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  
  
  NumericVector SA(numCohorts);
  for(int c=0;c<numCohorts;c++){
    SA[c] = 10000.0*(LAI_live[c]/(N[c]/10000.0))/Al2As[c];//Individual SA in cm2/m2
  }
  

  //Cohort description
  CharacterVector nsp = cohortCharacterParameter(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  DataFrame plantsdf;
  if(storagePool=="one") {
    NumericVector fastCstorage(numCohorts);
    for(int c=0;c<numCohorts;c++){
      NumericVector compartments = carbonCompartments(SA[c], LAI_expanded[c], H[c], 
                                                      Z[c], N[c], SLA[c], WoodDensity[c], WoodC[c]);
      fastCstorage[c] = 0.5*Cstoragepmax[c]*(compartments[0]+compartments[1]+compartments[2]);//Pool at 50%
    }
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                   _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                   _["SA"] = SA, _["fastCstorage"] = fastCstorage);
  } else if(storagePool=="two") {
    NumericVector slowCstorage(numCohorts), fastCstorage(numCohorts);
    for(int c=0;c<numCohorts;c++){
      NumericVector compartments = carbonCompartments(SA[c], LAI_expanded[c], H[c], 
                                                      Z[c], N[c], SLA[c], WoodDensity[c], WoodC[c]);
      slowCstorage[c] = 0.5*(Cstoragepmax[c]-0.05)*(compartments[0]+compartments[1]+compartments[2]); //Slow pool at 50%
      fastCstorage[c] = 0.5*0.05*(compartments[0]+compartments[1]+compartments[2]); //Fast pool at 50%
    }
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                 _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                   _["SA"] = SA, _["fastCstorage"] = fastCstorage, _["slowCstorage"] = slowCstorage);
  } else {
    plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                                 _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                                 _["SA"] = SA);
  }
  plantsdf.attr("row.names") = above.attr("row.names");
  
  
  //Base params
  DataFrame paramsBasedf = DataFrame::create(_["albedo"] = albedo, _["k"] = k, _["g"] = g, _["Sgdd"] = Sgdd);
  paramsBasedf.attr("row.names") = above.attr("row.names");
  

  List input;
  if(transpirationMode=="Granier") {
    NumericVector WUE = cohortNumericParameter(SP, SpParams, "WUE");
    NumericVector Psi_Extract = cohortNumericParameter(SP, SpParams, "Psi_Extract");
    NumericVector pRootDisc = cohortNumericParameter(SP, SpParams, "pRootDisc");
    
    DataFrame paramsTranspdf = DataFrame::create(_["Psi_Extract"]=Psi_Extract,_["WUE"] = WUE, _["pRootDisc"] = pRootDisc);
    paramsTranspdf.attr("row.names") = above.attr("row.names");
    
    List below = List::create( _["Z"]=Z,_["V"] = V);
    List paramsCanopy = List::create(_["gdd"] = 0);
    input = List::create(_["control"] = clone(control),
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsTransp"] = paramsTranspdf,
                         _["paramsGrowth"]= paramsGrowthdf,
                         _["paramsAllometries"] = paramsAllometriesdf);
    
    NumericVector tvec =  NumericVector(numCohorts, 0.0);
    tvec.attr("names") = above.attr("row.names");
    input["Transpiration"] = tvec;
    NumericVector pvec =  NumericVector(numCohorts, 0.0);
    pvec.attr("names") = above.attr("row.names");
    input["Photosynthesis"] = pvec;
    NumericVector cvec =  NumericVector(numCohorts, 0.0);
    cvec.attr("names") = above.attr("row.names");
    input["PLC"] = cvec;

  } else if(transpirationMode =="Sperry"){
    int numStemSegments = control["nStemSegments"];
    bool capacitance = control["capacitance"];

    DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams);
    DataFrame paramsWaterStoragedf = paramsWaterStorage(above, SpParams, paramsAnatomydf);
    DataFrame paramsTranspirationdf = paramsTranspiration(above, V, soil, SpParams,
                                                          paramsAnatomydf, control);
    List below = paramsBelowZ(above, V, Z, soil, 
                             paramsTranspirationdf, control);
    
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
    NumericMatrix psiRhizo =  NumericMatrix(numCohorts, nlayers);
    psiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
    NumericVector rwcsleaf = NumericVector(numCohorts, 1.0);
    rwcsleaf.attr("names") = above.attr("row.names");
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      warning("Soil pedotransfer functions set to Van Genuchten ('VG').");
    }
    List paramsCanopy = List::create(_["gdd"] = 0,_["Temp"] = NA_REAL);
    List ctl = clone(control);
    if(capacitance) {
      ctl["hydraulicCostFunction"] = 2;
      warning("Hydraulic cost function set to '2'.");
    }
    input = List::create(_["control"] = ctl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsBase"] = paramsBasedf,
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsTransp"] = paramsTranspirationdf,
                         _["paramsWaterStorage"] = paramsWaterStoragedf,
                         _["paramsGrowth"]= paramsGrowthdf,
                         _["paramsAllometries"] = paramsAllometriesdf);
    
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
    input["psiRhizo"] = psiRhizo;
    input["psiRoot"] = psiRoot;
    input["psiStem"] = psiStemmat;
    input["psiLeaf"] = psiLeaf;

  } 
  
  input.attr("class") = CharacterVector::create("growthInput","list");
  return(input);
}


// [[Rcpp::export("forest2spwbInput")]]
List forest2spwbInput(List x, List soil, DataFrame SpParams, List control) {
  NumericMatrix V = forest2belowground(x,soil, SpParams);
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(spwbInput(above,  V, soil, SpParams, control));
}


// [[Rcpp::export("forest2growthInput")]]
List forest2growthInput(List x, List soil, DataFrame SpParams, List control) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector d = soil["dVec"];
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  NumericVector treeZ50 = treeData["Z50"];
  NumericVector treeZ95 = treeData["Z95"];
  NumericVector shrubZ50 = shrubData["Z50"];  
  NumericVector shrubZ95 = shrubData["Z95"];  
  NumericMatrix V = forest2belowground(x,soil, SpParams);
  NumericVector Z(ntree+nshrub); //Rooting depth in cm
  for(int i=0;i<ntree;i++) {
    Z[i] = treeZ95[i]/10.0;
  }
  for(int i=0;i<nshrub;i++) {
    Z[ntree+i] = shrubZ95[i]/10.0; 
  }
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(growthInput(above,  Z, V, soil, SpParams, control));
}
