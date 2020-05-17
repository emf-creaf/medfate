#include <Rcpp.h>
#include "spwb.h"
#include "growth.h"
#include "carbon.h"
#include "root.h"
#include "woodformation.h"
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

DataFrame paramsPhenology(DataFrame above, DataFrame SpParams) {
  IntegerVector SP = above["SP"];
  int numCohorts = SP.size();
  CharacterVector phenoType = cohortCharacterParameter(SP, SpParams, "PhenologyType");
  NumericVector leafDuration  = cohortNumericParameter(SP, SpParams, "LeafDuration");
  NumericVector Sgdd  = cohortNumericParameter(SP, SpParams, "Sgdd");
  // NumericVector TbaseGdd = cohortNumericParameter(SP, SpParams, "TbaseGdd");
  // NumericVector Ssen = cohortNumericParameter(SP, SpParams, "Ssen");
  // NumericVector PstartSen = cohortNumericParameter(SP, SpParams, "PstartSen");
  // NumericVector TbaseSen = cohortNumericParameter(SP, SpParams, "TbaseSen");
  NumericVector Tbgdd(numCohorts, 0.0);
  NumericVector Ssen(numCohorts, 0.0);
  NumericVector Psen(numCohorts, 0.0);
  NumericVector Tbsen(numCohorts, 0.0);
  for(int j=0; j<numCohorts;j++) {
    if(phenoType[j] == "winter-deciduous" || phenoType[j] == "winter-semideciduous") { //Delpierre et al 2009
      Tbgdd[j]= 5.0;
      Ssen[j] = 8268.0;
      Psen[j] = 12.5;
      Tbsen[j] = 28.5;
    } else if(phenoType[j] == "oneflush-evergreen") {
      Tbgdd[j]= 5.0;
      Ssen[j] = 8268.0;
      Psen[j] = 12.5;
      Tbsen[j] = 8.5;
    }
  }
  DataFrame paramsPhenologydf = DataFrame::create(
    _["PhenologyType"] = phenoType,
    _["LeafDuration"] = leafDuration,
    _["Sgdd"] = Sgdd, _["Tbgdd"] = Tbgdd, 
    _["Ssen"] = Ssen, _["Psen"] = Psen, _["Tbsen"] = Tbsen 
  );
  paramsPhenologydf.attr("row.names") = above.attr("row.names");
  return(paramsPhenologydf);
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
  NumericVector Hmed = cohortNumericParameter(SP, SpParams, "Hmed"); //To correct conductivity
    
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(Al2As[c])) Al2As[c] = 2500.0; // = 4 cm2·m-2
  }
  DataFrame paramsAnatomydf = DataFrame::create(
    _["Hmed"] = Hmed,
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

DataFrame paramsTranspirationGranier(DataFrame above,  DataFrame SpParams) {
  IntegerVector SP = above["SP"];
  NumericVector WUE = cohortNumericParameter(SP, SpParams, "WUE");
  NumericVector Psi_Extract = cohortNumericParameter(SP, SpParams, "Psi_Extract");
  NumericVector Psi_Critic = cohortNumericParameter(SP, SpParams, "Psi_Critic");
  NumericVector pRootDisc = cohortNumericParameter(SP, SpParams, "pRootDisc");
  DataFrame paramsTranspirationdf = DataFrame::create(_["Psi_Extract"]=Psi_Extract,
                                                      _["Psi_Critic"] = Psi_Critic,
                                                      _["WUE"] = WUE, _["pRootDisc"] = pRootDisc);
  paramsTranspirationdf.attr("row.names") = above.attr("row.names");
  return(paramsTranspirationdf);
}
DataFrame paramsTranspirationSperry(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, 
                              DataFrame paramsAnatomydf, List control) {
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  int numCohorts = SP.size();
  NumericVector Vc;
  double fracRootResistance = control["fracRootResistance"];
  double fracLeafResistance = control["fracLeafResistance"];
  String transpirationMode = control["transpirationMode"];
  
  NumericVector dVec = soil["dVec"];
  
  CharacterVector Group = cohortCharacterParameter(SP, SpParams, "Group");
  CharacterVector Order = cohortCharacterParameter(SP, SpParams, "Order");
  CharacterVector GrowthForm = cohortCharacterParameter(SP, SpParams, "GrowthForm");
  CharacterVector phenoType = cohortCharacterParameter(SP, SpParams, "PhenologyType");
  
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
        if((GrowthForm[c]=="Shrub") && (phenoType[c] == "winter-deciduous")) {
          Kmax_stemxylem[c] = 1.55; //Angiosperm deciduous shrub
        } else if((GrowthForm[c]=="Tree" || GrowthForm[c]=="Tree/Shrub") && (phenoType[c] == "winter-deciduous")) {
          Kmax_stemxylem[c] = 1.58; //Angiosperm winter-deciduous tree
        } else { 
          Kmax_stemxylem[c] = 2.43; //Angiosperm evergreen tree
        }
      } else {
        if(GrowthForm[c]=="Shrub") {
          Kmax_stemxylem[c] = 0.24; //Gymnosperm shrub
        } else {
          Kmax_stemxylem[c] = 0.48; //Gymnosperm tree
        }
      }
    }
    //Oliveras I, Martínez-Vilalta J, Jimenez-Ortiz T, et al (2003) Hydraulic architecture of Pinus halepensis, P . pinea and Tetraclinis articulata in a dune ecosystem of Eastern Spain. Plant Ecol 131–141
    if(NumericVector::is_na(Kmax_rootxylem[c])) Kmax_rootxylem[c] = 4.0*Kmax_stemxylem[c];
    
    //Calculate stem maximum conductance (in mmol·m-2·s-1·MPa-1)
    VCstem_kmax[c]=maximumStemHydraulicConductance(Kmax_stemxylem[c], Hmed[c], Al2As[c],H[c],control["taper"]); 
    
    //Xylem vulnerability curve
    if(NumericVector::is_na(VCstem_d[c]) | NumericVector::is_na(VCstem_c[c])) {
      double psi50 = NA_REAL;
      // From: Maherali H, Pockman W, Jackson R (2004) Adaptive variation in the vulnerability of woody plants to xylem cavitation. Ecology 85:2184–2199
      if(Group[c]=="Angiosperm") {
        if((GrowthForm[c]=="Shrub") && (phenoType[c] != "winter-deciduous")) {
          psi50 = -5.09; //Angiosperm evergreen shrub
        } else if((GrowthForm[c]!="Shrub") && (phenoType[c] == "winter-deciduous")) {
          psi50 = -2.34; //Angiosperm winter-deciduous tree
        } else { 
          psi50 = -1.51; //Angiosperm evergreen tree
        }
      } else {
        if(GrowthForm[c]=="Shrub") {
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
      if(NumericVector::is_na(fracLeafResistance)) {
        if(Group[c]=="Angiosperm") {
          VCleaf_kmax[c] = 8.0;
        } else {
          VCleaf_kmax[c] = 6.0;
        }
      } else {
        double rstem = (1.0/VCstem_kmax[c]);
        double rtot = rstem/(1.0-fracRootResistance - fracLeafResistance);
        VCleaf_kmax[c] = 1.0/(rtot*fracLeafResistance);
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
  DataFrame paramsTranspirationdf = DataFrame::create(
    _["Gwmin"]=Gwmin, _["Gwmax"]=Gwmax,_["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298, _["Kmax_stemxylem"] = Kmax_stemxylem, _["Kmax_rootxylem"] = Kmax_rootxylem,
        _["VCleaf_kmax"]=VCleaf_kmax,_["VCleaf_c"]=VCleaf_c,_["VCleaf_d"]=VCleaf_d,
        _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, 
        _["VCroot_kmax"] = VCroottot_kmax ,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,
        _["Plant_kmax"] = Plant_kmax);
  paramsTranspirationdf.attr("row.names") = above.attr("row.names");
  return(paramsTranspirationdf);
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
  NumericVector fHDmin = cohortNumericParameter(SP, SpParams, "fHDmin");
  NumericVector fHDmax = cohortNumericParameter(SP, SpParams, "fHDmax");
  

  DataFrame paramsGrowthdf = DataFrame::create(_["WoodC"] = WoodC, 
                                               _["RGRmax"] = RGRmax,
                                               _["fHDmin"] = fHDmin,
                                               _["fHDmax"] = fHDmax);
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

  DataFrame paramsAllometriesdf = DataFrame::create(_["Hmax"] = Hmax,
                                                    _["Zmax"] = Zmax,
                                                    _["Aash"] = Aash, _["Absh"] = Absh, _["Bbsh"] = Bbsh,
                                                    _["r635"] = r635,
                                                    _["Acr"] = Acr, _["B1cr"] = B1cr, _["B2cr"] = B2cr, _["B3cr"] = B3cr,
                                                    _["C1cr"] = C1cr, _["C2cr"] = C2cr, 
                                                    _["Acw"] = Acw, _["Bcw"] = Bcw);
  paramsAllometriesdf.attr("row.names") = above.attr("row.names");
  return(paramsAllometriesdf);
}

DataFrame internalPhenologyDataFrame(DataFrame above) {
  int numCohorts = above.nrow();
  NumericVector phi(numCohorts,0.0);
  NumericVector gdd(numCohorts,0.0);
  NumericVector sen(numCohorts,0.0);
  LogicalVector budFormation(numCohorts, false);
  LogicalVector leafUnfolding(numCohorts, false);
  LogicalVector leafSenescence(numCohorts, false);
  LogicalVector leafDormancy(numCohorts, false);
  
  DataFrame df = DataFrame::create(Named("gdd") = gdd,
                                   Named("sen") = sen,
                                   Named("budFormation") = budFormation,
                                   Named("leafUnfolding") = leafUnfolding,
                                   Named("leafSenescence") = leafSenescence,
                                   Named("leafDormancy") = leafDormancy,
                                   Named("phi") = phi);
  df.attr("row.names") = above.attr("row.names");
  return(df);
}
DataFrame internalCarbonDataFrame(DataFrame above, 
                                  List below,
                                  DataFrame paramsAnatomydf,
                                  DataFrame paramsWaterStoragedf,
                                  DataFrame paramsGrowthdf,
                                  List control) {
  int numCohorts = above.nrow();

  double nonSugarConc = control["nonSugarConc"];
  String transpirationMode = control["transpirationMode"];
  
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector LeafDensity = paramsAnatomydf["LeafDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector LeafPI0 = paramsWaterStoragedf["LeafPI0"];
  NumericVector StemPI0 = paramsWaterStoragedf["StemPI0"];
  
  NumericVector WoodC = paramsGrowthdf["WoodC"];

  NumericVector Z = below["Z"];
  IntegerVector SP = above["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector N = above["N"];
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  NumericVector SA = above["SA"];

  // NumericVector sugar(numCohorts,0.0);
  // NumericVector starch(numCohorts,0.0);
  NumericVector sugarLeaf(numCohorts,0.0);
  NumericVector starchLeaf(numCohorts,0.0);
  NumericVector sugarSapwood(numCohorts,0.0);
  NumericVector starchSapwood(numCohorts,0.0);
  // NumericVector longtermStorage(numCohorts,0.0);
  for(int c=0;c<numCohorts;c++){
    double lvol = leafStorageVolume(LAI_expanded[c],  N[c], SLA[c], LeafDensity[c]);
    double svol = sapwoodStorageVolume(SA[c], H[c],Z[c],WoodDensity[c], 0.5);
    
    // 50% in starch storage
    starchLeaf[c] = (0.5/(lvol*glucoseMolarMass))*leafStarchCapacity(LAI_expanded[c], N[c], SLA[c], LeafDensity[c]);
    starchSapwood[c] = (0.5/(svol*glucoseMolarMass))*sapwoodStarchCapacity(SA[c], H[c], Z[c], WoodDensity[c], 0.5);
    // starch[c] = starchLeaf[c]+starchSapwood[c];
    
    //Sugar storage from PI0
    double lconc = sugarConcentration(LeafPI0[c],20.0, nonSugarConc);
    sugarLeaf[c] = lconc;
    double sconc = sugarConcentration(StemPI0[c],20.0, nonSugarConc);
    sugarSapwood[c] = sconc;
    // sugar[c] = sugarLeaf[c] + sugarSapwood[c];
  }
  DataFrame df;
  // if(transpirationMode=="Granier"){
  //   df = DataFrame::create(Named("sugar") = sugar,
  //                          Named("starch") = starch);
  // } else {
    df = DataFrame::create(Named("sugarLeaf") = sugarLeaf,
                           Named("starchLeaf") = starchLeaf,
                           Named("sugarSapwood") = sugarSapwood,
                           Named("starchSapwood") = starchSapwood);
  // }
  df.attr("row.names") = above.attr("row.names");
  return(df);
}  


DataFrame internalAllocationDataFrame(DataFrame above, 
                                  DataFrame paramsAnatomydf,
                                  DataFrame paramsTranspirationdf,
                                  List control) {
  int numCohorts = above.nrow();

  NumericVector allocationTarget(numCohorts,0.0);
  NumericVector leafAreaTarget(numCohorts,0.0);
  
  String transpirationMode = control["transpirationMode"];
  NumericVector SA = above["SA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  if(transpirationMode=="Granier") {
    for(int c=0;c<numCohorts;c++){
      leafAreaTarget[c] = Al2As[c]*(SA[c]/10000.0);
      allocationTarget[c] = Al2As[c];
    }
  } else {
    String allocationStrategy = control["allocationStrategy"];
    NumericVector Plant_kmax = paramsTranspirationdf["Plant_kmax"];
    // NumericVector longtermStorage(numCohorts,0.0);
    for(int c=0;c<numCohorts;c++){
      leafAreaTarget[c] = Al2As[c]*(SA[c]/10000.0);
      if(allocationStrategy=="Plant_kmax") {
        allocationTarget[c] = Plant_kmax[c];
      } else if(allocationStrategy=="Al2As") {
        allocationTarget[c] = Al2As[c];
      }
    }
  }
  
  DataFrame df = DataFrame::create(Named("allocationTarget") = allocationTarget,
                                   Named("leafAreaTarget") = leafAreaTarget);
  df.attr("row.names") = above.attr("row.names");
  return(df);
}  


DataFrame internalWaterDataFrame(DataFrame above, String transpirationMode) {
  int numCohorts = above.nrow();
  DataFrame df;
  if(transpirationMode=="Granier") {
    df = DataFrame::create(Named("PlantPsi") = NumericVector(numCohorts, 0.0),
                           Named("StemPLC") = NumericVector(numCohorts, 0.0));
  } else {
    df = DataFrame::create(Named("Einst") = NumericVector(numCohorts, 0.0),
                           Named("RootCrownPsi") = NumericVector(numCohorts, 0.0),
                           Named("Stem1Psi") = NumericVector(numCohorts, 0.0),
                           Named("Stem2Psi") = NumericVector(numCohorts, 0.0),
                           Named("LeafPsi") = NumericVector(numCohorts, 0.0),
                           Named("StemSympPsi") = NumericVector(numCohorts, 0.0),
                           Named("LeafSympPsi") = NumericVector(numCohorts, 0.0),
                           Named("StemPLC") = NumericVector(numCohorts, 0.0),
                           Named("NSPL") = NumericVector(numCohorts, 1.0));
  }
  df.attr("row.names") = above.attr("row.names");
  return(df);
}

/**
 *  Prepare Soil Water Balance input
 */
// [[Rcpp::export("spwbInput")]]
List spwbInput(DataFrame above, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  
  int nlayers = V.ncol();
  
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
  

  NumericVector alphaSWR = cohortNumericParameter(SP, SpParams, "alphaSWR");
  NumericVector gammaSWR = cohortNumericParameter(SP, SpParams, "gammaSWR");
  NumericVector kPAR = cohortNumericParameter(SP, SpParams, "kPAR");
  NumericVector g = cohortNumericParameter(SP, SpParams, "g");
  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  
  StringVector Status(numCohorts, "alive");
  

  //Cohort description
  CharacterVector nsp = cohortCharacterParameter(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  //Above 
  DataFrame plantsdf = DataFrame::create(_["H"]=H, _["CR"]=CR, 
                                         _["LAI_live"]=LAI_live, 
                                         _["LAI_expanded"] = LAI_expanded, 
                                         _["LAI_dead"] = LAI_dead,
                                         _["Status"] = Status);
  plantsdf.attr("row.names") = above.attr("row.names");
  
 
  NumericVector Wsoil = soil["W"];
  NumericMatrix Wpool = NumericMatrix(numCohorts, nlayers);
  Wpool.attr("dimnames") = V.attr("dimnames");
  for(int c=0;c<numCohorts;c++){
    for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l]; //Init from soil state
  }
  
  
  List input;
  if(transpirationMode=="Granier") {
    
    //Base params
    DataFrame paramsInterceptiondf = DataFrame::create(_["kPAR"] = kPAR, 
                                               _["g"] = g);
    paramsInterceptiondf.attr("row.names") = above.attr("row.names");
  
    List below = List::create(_["V"] = V,
                              _["Wpool"] = Wpool);
    List paramsCanopy = List::create(_["gdd"] = 0.0);
    input = List::create(_["control"] = clone(control),
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsPhenology"] = paramsPhenology(above, SpParams),
                         _["paramsInterception"] = paramsInterceptiondf,
                         _["paramsTranspiration"] = paramsTranspirationGranier(above,SpParams),
                         _["internalPhenology"] = internalPhenologyDataFrame(above),
                         _["internalWater"] = internalWaterDataFrame(above, transpirationMode));
  } else if(transpirationMode =="Sperry"){
    
    //Base params
    DataFrame paramsInterceptiondf = DataFrame::create(_["alphaSWR"] = alphaSWR,
                                               _["gammaSWR"] = gammaSWR, 
                                               _["kPAR"] = kPAR, 
                                               _["g"] = g);
    paramsInterceptiondf.attr("row.names") = above.attr("row.names");
    
    DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams);
    DataFrame paramsTranspirationdf = paramsTranspirationSperry(above, V, soil, SpParams,
                                                          paramsAnatomydf, control);
    DataFrame paramsWaterStoragedf = paramsWaterStorage(above, SpParams, paramsAnatomydf);
    NumericMatrix RhizoPsi =  NumericMatrix(numCohorts, nlayers);
    RhizoPsi.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
    std::fill(RhizoPsi.begin(), RhizoPsi.end(), 0.0);
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      warning("Soil pedotransfer functions set to Van Genuchten ('VG').");
    }
    List below = paramsBelow(above, V, soil, 
                             paramsTranspirationdf, control);
    below["Wpool"] = Wpool;
    below["RhizoPsi"] = RhizoPsi;
    
    List paramsCanopy = List::create(_["gdd"] = 0.0,_["Temp"] = NA_REAL);
    List ctl = clone(control);
    input = List::create(_["control"] = ctl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsPhenology"] = paramsPhenology(above, SpParams),
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsInterception"] = paramsInterceptiondf,
                         _["paramsTranspiration"] = paramsTranspirationdf,
                         _["paramsWaterStorage"] = paramsWaterStoragedf,
                         _["internalPhenology"] = internalPhenologyDataFrame(above),
                         _["internalWater"] = internalWaterDataFrame(above, transpirationMode));
  }

  input.attr("class") = CharacterVector::create("spwbInput","list");
  return(input);
}



/**
 *  Prepare Forest growth input
 */
// [[Rcpp::export("growthInput")]]
List growthInput(DataFrame above, NumericVector Z, NumericMatrix V, List soil, DataFrame SpParams, List control) {
  
  int nlayers = V.ncol();
  
  IntegerVector SP = above["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector N = above["N"];
  NumericVector DBH = above["DBH"];
  NumericVector Cover = above["Cover"];
  NumericVector H = above["H"];
  NumericVector CR = above["CR"];
  
  control["cavitationRefill"] = "growth";
  
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Granier") & (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Granier' or 'Sperry')");

  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") & (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  

  
  DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams);
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  
  DataFrame paramsGrowthdf = paramsGrowth(above, SpParams);

  DataFrame paramsWaterStoragedf = paramsWaterStorage(above, SpParams, paramsAnatomydf);
  
  DataFrame paramsAllometriesdf = paramsAllometries(above, SpParams);
  
  NumericVector alphaSWR = cohortNumericParameter(SP, SpParams, "alphaSWR");
  NumericVector gammaSWR = cohortNumericParameter(SP, SpParams, "gammaSWR");
  NumericVector kPAR = cohortNumericParameter(SP, SpParams, "kPAR");
  NumericVector g = cohortNumericParameter(SP, SpParams, "g");

  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  
  
  NumericVector SA(numCohorts);
  StringVector Status(numCohorts, "alive");
  for(int c=0;c<numCohorts;c++){
    SA[c] = 10000.0*(LAI_live[c]/(N[c]/10000.0))/Al2As[c];//Individual SA in cm2
  }
  
  NumericVector Wsoil = soil["W"];
  NumericMatrix Wpool = NumericMatrix(numCohorts, nlayers);
  Wpool.attr("dimnames") = V.attr("dimnames");
  for(int c=0;c<numCohorts;c++){
    for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l]; //Init from soil state
  }
  
  //Cohort description
  CharacterVector nsp = cohortCharacterParameter(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  DataFrame plantsdf = DataFrame::create(_["SP"]=SP, _["N"]=N,_["DBH"]=DBH, _["Cover"] = Cover, _["H"]=H, _["CR"]=CR,
                               _["LAI_live"]=LAI_live, _["LAI_expanded"]=LAI_expanded, _["LAI_dead"] = LAI_dead,  
                               _["SA"] = SA, _["Status"] = Status);
  plantsdf.attr("row.names") = above.attr("row.names");
  

  List ringList(numCohorts);
  for(int i=0;i<numCohorts;i++) ringList[i] = initialize_ring();
  ringList.attr("names") = above.attr("row.names");
  
  List input;
  if(transpirationMode=="Granier") {
    //Base params
    DataFrame paramsInterceptiondf = DataFrame::create(_["kPAR"] = kPAR, 
                                               _["g"] = g);
    paramsInterceptiondf.attr("row.names") = above.attr("row.names");
    
    DataFrame paramsTranspirationdf = paramsTranspirationGranier(above, SpParams);
    
    List below = List::create( _["Z"]=Z,_["V"] = V,_["Wpool"] = Wpool);
    List paramsCanopy = List::create(_["gdd"] = 0.0);
    input = List::create(_["control"] = clone(control),
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsPhenology"] = paramsPhenology(above, SpParams),
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsInterception"] = paramsInterceptiondf,
                         _["paramsTranspiration"] = paramsTranspirationdf,
                         _["paramsWaterStorage"] = paramsWaterStoragedf,
                         _["paramsGrowth"]= paramsGrowthdf,
                         _["paramsAllometries"] = paramsAllometriesdf,
                         _["internalPhenology"] = internalPhenologyDataFrame(above),
                         _["internalWater"] = internalWaterDataFrame(above, transpirationMode),
                         _["internalCarbon"] = internalCarbonDataFrame(plantsdf, below, 
                                                         paramsAnatomydf, 
                                                         paramsWaterStoragedf,
                                                         paramsGrowthdf, control),
                        _["internalAllocation"] = internalAllocationDataFrame(plantsdf,
                                                            paramsAnatomydf,
                                                            paramsTranspirationdf, control),
                        _["rings"] = ringList);
  } else if(transpirationMode =="Sperry"){
    
    //Base params
    DataFrame paramsInterceptiondf = DataFrame::create(_["alphaSWR"] = alphaSWR,
                                               _["gammaSWR"] = gammaSWR, 
                                               _["kPAR"] = kPAR, 
                                               _["g"] = g);
    paramsInterceptiondf.attr("row.names") = above.attr("row.names");
    
    DataFrame paramsTranspirationdf = paramsTranspirationSperry(above, V, soil, SpParams,
                                                          paramsAnatomydf, control);
    List below = paramsBelowZ(above, V, Z, soil, 
                              paramsTranspirationdf, control);
    
    NumericMatrix RhizoPsi =  NumericMatrix(numCohorts, nlayers);
    RhizoPsi.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
    std::fill(RhizoPsi.begin(), RhizoPsi.end(), 0.0);
    below["RhizoPsi"] = RhizoPsi;
    below["Wpool"] = Wpool;
    
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      warning("Soil pedotransfer functions set to Van Genuchten ('VG').");
    }
    List paramsCanopy = List::create(_["gdd"] = 0.0,_["Temp"] = NA_REAL);
    List ctl = clone(control);
    input = List::create(_["control"] = ctl,
                         _["canopy"] = paramsCanopy,
                         _["cohorts"] = cohortDescdf,
                         _["above"] = plantsdf,
                         _["below"] = below,
                         _["paramsPhenology"] = paramsPhenology(above, SpParams),
                         _["paramsAnatomy"] = paramsAnatomydf,
                         _["paramsInterception"] = paramsInterceptiondf,
                         _["paramsTranspiration"] = paramsTranspirationdf,
                         _["paramsWaterStorage"] = paramsWaterStoragedf,
                         _["paramsGrowth"]= paramsGrowthdf,
                         _["paramsAllometries"] = paramsAllometriesdf,
                         _["internalPhenology"] = internalPhenologyDataFrame(above),
                         _["internalWater"] = internalWaterDataFrame(above, transpirationMode),
                         _["internalCarbon"] = internalCarbonDataFrame(plantsdf, below,
                                                         paramsAnatomydf, 
                                                         paramsWaterStoragedf,
                                                         paramsGrowthdf, control),
                         _["internalAllocation"] = internalAllocationDataFrame(plantsdf,
                                                         paramsAnatomydf,
                                                         paramsTranspirationdf, control),
                         _["rings"] = ringList);
    
  } 
  
  input.attr("class") = CharacterVector::create("growthInput","list");
  return(input);
}


// [[Rcpp::export("forest2spwbInput")]]
List forest2spwbInput(List x, List soil, DataFrame SpParams, List control, String mode = "MED") {
  NumericMatrix V = forest2belowground(x,soil, SpParams);
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL, mode);
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

// [[Rcpp::export("resetInputs")]]
void resetInputs(List x, List soil) {
  List can = x["canopy"];
  NumericVector Wsoil = soil["W"];
  NumericVector Temp = soil["Temp"];
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  List below = x["below"];
  NumericMatrix Wpool = below["Wpool"];
  int nlayers = Wsoil.size();
  int numCohorts = Wpool.nrow();
  
  can["gdd"] = 0.0;
  can["Temp"] = NA_REAL;
  for(int i=0;i<nlayers;i++) {
    Wsoil[i] = 1.0; //Defaults to soil at field capacity
    Temp[i] = NA_REAL;
  }
  for(int c=0;c<numCohorts;c++) {
    for(int l=0;l<nlayers;l++) {
      Wpool(c,l) = 1.0;
    }
  }
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  
  if(transpirationMode=="Sperry") {
    NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(below["RhizoPsi"]);
    NumericVector RootCrownPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
    NumericVector Stem1Psi = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
    NumericVector Stem2Psi = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem2Psi"]);
    NumericVector StemSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
    NumericVector LeafSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
    NumericVector LeafPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
    NumericVector NSPL = Rcpp::as<Rcpp::NumericVector>(internalWater["NSPL"]);
    for(int i=0;i<LeafPsi.size();i++) {
      Einst[i] = 0.0;
      RootCrownPsi[i] = 0.0;
      Stem1Psi[i] = 0.0;
      Stem2Psi[i] = 0.0;
      LeafPsi[i] = 0.0;
      LeafSympPsi[i] = 0.0;
      StemSympPsi[i] = 0.0;
      StemPLC[i] = 0.0;
      NSPL[i] = 1.0;
      for(int j=0;j<RhizoPsi.ncol();j++) RhizoPsi(i,j) = 0.0;
    }
  } else {
    NumericVector PlantPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
    for(int i=0;i<StemPLC.length();i++) {
      PlantPsi[i] = 0.0;
      StemPLC[i] = 0.0;
    }
  }
}
