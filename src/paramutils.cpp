#define STRICT_R_HEADERS
#include <numeric>
#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "forestutils.h"
using namespace Rcpp;

int findRowIndex(int sp, DataFrame SpParams) {
  IntegerVector spIndexSP = SpParams["SpIndex"];
  for(int i=0;i<spIndexSP.length();i++) if(spIndexSP[i]==sp) return(i);
  Rcerr << sp << " not found!\n";
  stop("Species code not found in SpParams");
  return(NA_INTEGER);
}

// [[Rcpp::export(".checkSpeciesParameters")]]
void checkSpeciesParameters(DataFrame SpParams, CharacterVector params) {
  NumericVector values;
  String s;
  for(int i =0;i<params.size();i++){
    s = params[i];
    if(!SpParams.containsElementNamed(params[i])) {
      Rcout << params[i]<<"\n";
      stop("Parameter column missing in species params");
    }
  }
}

NumericVector speciesNumericParameter(IntegerVector SP, DataFrame SpParams, String parName){
  NumericVector par(SP.size(), NA_REAL);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    NumericVector parSP = Rcpp::as<Rcpp::NumericVector>(SpParams[parName]);
    for(int i=0;i<SP.size();i++) {
      int iSP = findRowIndex(SP[i], SpParams);
      par[i] = parSP[iSP];
    }
  } else {
    Rcerr << "Variable '" << parName.get_cstring() << "' was not found in SpParams!\n";
  }
  return(par);
}

// [[Rcpp::export("species_characterParameter")]]
CharacterVector speciesCharacterParameter(IntegerVector SP, DataFrame SpParams, String parName){
  CharacterVector par(SP.size(), NA_STRING);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    CharacterVector parSP = SpParams[parName];
    for(int i=0;i<SP.size();i++) {
      int iSP = findRowIndex(SP[i], SpParams);
      par[i] = parSP[iSP];
    }
  } else {
    Rcerr << "Variable '" << parName.get_cstring() << "' was not found in SpParams!\n";
  }
  return(par);
}

NumericVector cohortNumericParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector tSP = treeData["Species"];
  IntegerVector shSP = shrubData["Species"];
  NumericVector par(tSP.size()+shSP.size());
  NumericVector parTrees = speciesNumericParameter(tSP, SpParams, parName);
  NumericVector parShrubs = speciesNumericParameter(shSP, SpParams, parName);
  for(int i=0;i<tSP.size();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shSP.size();i++) {
    par[i+tSP.size()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x);
  return(par);
}

// [[Rcpp::export("plant_characterParameter")]]
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector tSP = treeData["Species"];
  IntegerVector shSP = shrubData["Species"];
  CharacterVector par(tSP.size()+shSP.size());
  CharacterVector parTrees = speciesCharacterParameter(tSP, SpParams, parName);
  CharacterVector parShrubs = speciesCharacterParameter(shSP, SpParams, parName);
  for(int i=0;i<tSP.size();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shSP.size();i++) {
    par[i+tSP.size()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x);
  return(par);
}

/** Parameter retrieval with imputation */
NumericVector kPARWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  NumericVector kPAR = speciesNumericParameter(SP, SpParams, "kPAR");
  for(int j=0;j<kPAR.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.55;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.45;
    } else if(leafShape[j]=="Needle" | leafShape[j]=="Scale"){
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.50;
    }
  }
  return(kPAR);
}
NumericVector gammaSWRWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  NumericVector gammaSWR = speciesNumericParameter(SP, SpParams, "gammaSWR");
  for(int j=0;j<gammaSWR.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.18;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.15;
    } else if(leafShape[j]=="Needle" | leafShape[j]=="Scale"){
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.14;
    }
  }
  return(gammaSWR);
}
NumericVector alphaSWRWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector alphaSWR = speciesNumericParameter(SP, SpParams, "alphaSWR");
  for(int j=0;j<alphaSWR.size();j++) {
    if(NumericVector::is_na(alphaSWR[j])) alphaSWR[j] = 0.7;
  }
  return(alphaSWR);
}
NumericVector gWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  NumericVector g = speciesNumericParameter(SP, SpParams, "g");
  for(int j=0;j<g.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(g[j])) g[j] = 0.5;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(g[j])) g[j] = 0.8;
    } else if(leafShape[j]=="Needle" | leafShape[j]=="Scale"){
      if(NumericVector::is_na(g[j])) g[j] = 1.0;
    }
  }
  return(g);
}
NumericVector fineFoliarRatioWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector ffr = speciesNumericParameter(SP, SpParams, "r635");
  for(int i=0;i<ffr.size();i++) {
    if(NumericVector::is_na(ffr[i])) {
      if(leafShape[i]=="Scale") {
        ffr[i] = 1.482607;
      } else if(leafShape[i]=="Spines") {
        ffr[i] = NA_REAL;
      } else if(leafShape[i]=="Linear" ) {
        ffr[i] = 3.260730;
      } else if(leafShape[i]=="Needle" ) {
        ffr[i] = 1.715895;
      } else { //Broad
        if(leafSize[i]=="Small") {
          ffr[i] = 3.025709;
        } else if(leafSize[i] == "Medium") {
          ffr[i] = 2.358575;
        } else { //large
          ffr[i] = 2.277993;
        }
      }
    }
  }
  return(ffr);
}
NumericVector specificLeafAreaWithImputation(IntegerVector SP, DataFrame SpParams){
  CharacterVector LeafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  CharacterVector LeafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  NumericVector SLA = speciesNumericParameter(SP, SpParams, "SLA");
  for(int c=0;c<SLA.size();c++){
    if(NumericVector::is_na(SLA[c])) {
      if(!CharacterVector::is_na(LeafShape[c]) && !CharacterVector::is_na(LeafSize[c])) {
        if(LeafShape[c]=="Linear") {
          if(LeafSize[c]=="Small") {
            SLA[c] = 13.189;
          } else if(LeafSize[c]=="Medium") {
            SLA[c] = 4.144;
          } else if(LeafSize[c]=="Large") {
            SLA[c]= 5.522;
          }
        } else if(LeafShape[c]=="Broad") {
          if(LeafSize[c]=="Small") {
            SLA[c] = 9.540;
          } else if(LeafSize[c]=="Medium") {
            SLA[c] = 11.499;
          } else if(LeafSize[c]=="Large") {
            SLA[c]= 16.039;
          }
        } else if(LeafShape[c]=="Needle") { 
          SLA[c] = 9.024;
        } else if(LeafShape[c]=="Scale") { 
          SLA[c] = 4.544;
        }
      } else {
        SLA[c] = 4.0;
      }
    } 
  }
  return(SLA);
}
NumericVector ligninPercentWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector cohLigninPercent = speciesNumericParameter(SP, SpParams, "LigninPercent");
  for(int i=0;i<cohLigninPercent.size();i++) {
    if(NumericVector::is_na(cohLigninPercent[i])) {
      if(leafShape[i]=="Scale") {
        cohLigninPercent[i] = 14.55;
      } else if(leafShape[i]=="Spines") {
        cohLigninPercent[i] = 14.55;
      } else if(leafShape[i]=="Linear" | leafShape[i]=="Needle") {
        if(leafSize[i]=="Small") {
          cohLigninPercent[i] = 18.55;
        } else if(leafSize[i] == "Medium") {
          cohLigninPercent[i] = 24.52;
        } else { 
          cohLigninPercent[i] = 24.52;
        }
      } else { //Broad
        if(leafSize[i]=="Small") {
          cohLigninPercent[i] = 22.32;
        } else if(leafSize[i] == "Medium") {
          cohLigninPercent[i] = 20.21;
        } else { 
          cohLigninPercent[i] = 15.50;
        }
      }
    }
  }
  return(cohLigninPercent);
}
NumericVector surfaceToAreaRatioWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector cohSAV = speciesNumericParameter(SP, SpParams, "SAV");
  for(int i=0;i<cohSAV.size();i++) {
    if(NumericVector::is_na(cohSAV[i])) {
      if(leafShape[i]=="Scale") {
        cohSAV[i] = 1120.0;
      } else if(leafShape[i]=="Spines") {
        cohSAV[i] = 6750.0;
      } else if(leafShape[i]=="Linear" | leafShape[i]=="Needle") {
        if(leafSize[i]=="Small") {
          cohSAV[i] = 3620.0;
        } else if(leafSize[i] == "Medium") {
          cohSAV[i] = 4758.0;
        } else { 
          cohSAV[i] = 3620.0;
        }
      } else { //Broad
        if(leafSize[i]=="Small") {
          cohSAV[i] = 4386.0;
        } else if(leafSize[i] == "Medium") {
          cohSAV[i] = 4039.0;
        } else { 
          cohSAV[i] = 5740.0;
        }
      }
    }
  }
  return(cohSAV);
}
NumericVector heatContentWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector cohHeatContent = speciesNumericParameter(SP, SpParams, "HeatContent");
  for(int i=0;i<cohHeatContent.size();i++) {
    if(NumericVector::is_na(cohHeatContent[i])) {
      if(leafShape[i]=="Scale") {
        cohHeatContent[i] = 20504.0;
      } else if(leafShape[i]=="Spines") {
        cohHeatContent[i] = 20433.0;
      } else if(leafShape[i]=="Linear" | leafShape[i]=="Needle") {
        if(leafSize[i]=="Small") {
          cohHeatContent[i] = 21888.0;
        } else if(leafSize[i] == "Medium") {
          cohHeatContent[i] = 21182.0;
        } else { 
          cohHeatContent[i] = 18250.0;
        }
      } else { //Broad
        if(leafSize[i]=="Small") {
          cohHeatContent[i] = 20062.0;
        } else if(leafSize[i] == "Medium") {
          cohHeatContent[i] = 19825.0;
        } else { 
          cohHeatContent[i] = 19740.0;
        }
      }
    }
  }
  return(cohHeatContent);
}
NumericVector leafWidthWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector leafwidth = speciesNumericParameter(SP, SpParams, "LeafWidth");
  for(int c=0;c<leafwidth.size();c++) {
    if(NumericVector::is_na(leafwidth[c])) {
      if(leafShape[c]=="Linear") {
        leafwidth[c]= 0.6393182;
      } else if(leafShape[c]=="Needle") {
        leafwidth[c]= 0.3792844;
      } else if(leafShape[c]=="Broad") {
        if(leafSize[c]=="Small") {
          leafwidth[c] = 0.6439761;
        } else if(leafSize[c]=="Medium") {
          leafwidth[c] = 3.0537686;
        } else if(leafSize[c]=="Large") {
          leafwidth[c]= 6.8980354;
        }
      } else if(leafShape[c]=="Scale") { 
        leafwidth[c] = 0.1007839;
      }
    }
  }
  return(leafwidth);
}
NumericVector Al2AsWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameter(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameter(SP, SpParams, "LeafSize");
  NumericVector Al2As = speciesNumericParameter(SP, SpParams, "Al2As");
  for(int c=0;c<Al2As.size();c++) {
    if(NumericVector::is_na(Al2As[c])) {
      if(leafShape[c]=="Linear") {
        Al2As[c]= 2156.0;
      } else if(leafShape[c]=="Needle") {
        Al2As[c]= 2751.7;
      } else if(leafShape[c]=="Broad") {
        if(leafSize[c]=="Small") {
          Al2As[c] = 2284.9;
        } else if(leafSize[c]=="Medium") {
          Al2As[c] = 2446.1;
        } else if(leafSize[c]=="Large") {
          Al2As[c]= 4768.7;
        }
      } else if(leafShape[c]=="Scale") { 
        Al2As[c] = 1696.6;
      }
    }
  }
  return(Al2As);
}
NumericVector woodDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector woodDensity = speciesNumericParameter(SP, SpParams, "WoodDensity");
  for(int c=0;c<woodDensity.size();c++) {
    if(NumericVector::is_na(woodDensity[c])) {
      woodDensity[c] = 0.652;
    }
  }
  return(woodDensity);
}
NumericVector leafDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafDensity = speciesNumericParameter(SP, SpParams, "LeafDensity");
  for(int c=0;c<leafDensity.size();c++) {
    if(NumericVector::is_na(leafDensity[c])) {
      leafDensity[c] = 0.7;
    }
  }
  return(leafDensity);
}
NumericVector fineRootDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector fineRootDensity = speciesNumericParameter(SP, SpParams, "FineRootDensity");
  for(int c=0;c<fineRootDensity.size();c++) {
    if(NumericVector::is_na(fineRootDensity[c])) {
      fineRootDensity[c] = 0.165;
    }
  }
  return(fineRootDensity);
}
NumericVector specificRootLengthWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector specificRootLength = speciesNumericParameter(SP, SpParams, "SRL");
  for(int c=0;c<specificRootLength.size();c++) {
    if(NumericVector::is_na(specificRootLength[c])) {
      specificRootLength[c] = 3870.0;
    }
  }
  return(specificRootLength);
}
NumericVector rootLengthDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector rootLengthDensity = speciesNumericParameter(SP, SpParams, "RLD");
  for(int c=0;c<rootLengthDensity.size();c++) {
    if(NumericVector::is_na(rootLengthDensity[c])) {
      rootLengthDensity[c] = 10.0;
    }
  }
  return(rootLengthDensity);
}
NumericVector conduit2sapwoodWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector Group = speciesCharacterParameter(SP, SpParams, "Group");
  NumericVector conduit2sapwood = speciesNumericParameter(SP, SpParams, "conduit2sapwood");
  for(int c=0;c<conduit2sapwood.size();c++) {
    if(NumericVector::is_na(conduit2sapwood[c])) {
      if(Group[c]=="Angiosperm") conduit2sapwood[c] = 0.70; //20-40% parenchyma in angiosperms.
      else conduit2sapwood[c] = 0.925; //5-10% parenchyma in gymnosperms (https://link.springer.com/chapter/10.1007/978-3-319-15783-2_8)
    }
  }
  return(conduit2sapwood);
}
NumericVector stemPI0WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WoodDensity = woodDensityWithImputation(SP, SpParams);
  NumericVector StemPI0 = speciesNumericParameter(SP, SpParams, "StemPI0");
  for(int c=0;c<StemPI0.size();c++) {
    //From: Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N.M., Galbraith, D.R., Baker, T.R., Rowland, L., Fisher, R.A., Binks, O.J., Sevanto, S.A., Xu, C., Jansen, S., Choat, B., Mencuccini, M., McDowell, N.G., & Meir, P. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0: 1–60.
    if(NumericVector::is_na(StemPI0[c])) {
      StemPI0[c] = 0.52 - 4.16*WoodDensity[c];
    }
  }
  return(StemPI0);
}
NumericVector stemEPSWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WoodDensity = woodDensityWithImputation(SP, SpParams);
  NumericVector StemEPS = speciesNumericParameter(SP, SpParams, "StemEPS");
  for(int c=0;c<StemEPS.size();c++) {
    //From: Christoffersen, B.O., Gloor, M., Fauset, S., Fyllas, N.M., Galbraith, D.R., Baker, T.R., Rowland, L., Fisher, R.A., Binks, O.J., Sevanto, S.A., Xu, C., Jansen, S., Choat, B., Mencuccini, M., McDowell, N.G., & Meir, P. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0: 1–60.
    if(NumericVector::is_na(StemEPS[c])) {
      StemEPS[c] = sqrt(1.02*exp(8.5*WoodDensity[c])-2.89); 
    }
  }
  return(StemEPS);
}
NumericVector stemAFWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector conduit2sapwood = conduit2sapwoodWithImputation(SP, SpParams);
  NumericVector stemAF = speciesNumericParameter(SP, SpParams, "StemAF");
  for(int c=0;c<stemAF.size();c++) {
    if(NumericVector::is_na(stemAF[c])) {
      stemAF[c] = conduit2sapwood[c]; 
    }
  }
  return(stemAF);
}
NumericVector leafPI0WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafPI0 = speciesNumericParameter(SP, SpParams, "LeafPI0");
  for(int c=0;c<leafPI0.size();c++) {
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(leafPI0[c])) {
      leafPI0[c] = -2.0;//Average for Mediterranean climate species
    }
  }
  return(leafPI0);
}
NumericVector leafEPSWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafEPS = speciesNumericParameter(SP, SpParams, "LeafEPS");
  for(int c=0;c<leafEPS.size();c++) {
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(leafEPS[c])) {
      leafEPS[c] = 17.0;//Average for Mediterranean climate species
    }
  }
  return(leafEPS);
}
NumericVector leafAFWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafAF = speciesNumericParameter(SP, SpParams, "LeafAF");
  for(int c=0;c<leafAF.size();c++) {
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(leafAF[c])) {
      leafAF[c] = 0.29; //Average for Mediterranean climate species
    }
  }
  return(leafAF);
}

/** Allometric coefficient retrieval with imputation */
NumericVector shrubAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName) {
  NumericVector coef = speciesNumericParameter(SP,SpParams, parName);
  CharacterVector lifeForm = speciesCharacterParameter(SP, SpParams, "LifeForm");
  NumericVector Hmax = speciesNumericParameter(SP, SpParams, "Hmax");
  for(int i=0;i<coef.size();i++) { // From De Caceres et al. 2019
    if(NumericVector::is_na(coef[i])) {
      if(parName=="a_ash") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 24.5888;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 1.0083;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 5.8458;
        else coef[i] = 24.5888; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="b_ash") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 1.1662;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 1.8700;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 1.4944;
        else coef[i] = 1.1662; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="a_bsh") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 0.7963;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 0.7900;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 0.3596;
        else coef[i] = 0.7963; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="b_bsh") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 0.3762;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 0.6942;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 0.7138;
        else coef[i] = 0.3762; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="a_btsh") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 1.9189;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 1.2694;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 0.7856;
        else coef[i] = 1.9189; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="b_btsh") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 0.6873;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 0.7610;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 0.8101;
        else coef[i] = 0.6873; //Cryptophytes and hemicryptophytes like chamaephytes
      }
      else if(parName=="cr") {
        if(lifeForm[i]=="Chamaephyte") coef[i] = 0.8076;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]<300)) coef[i] = 0.663;
        else if(lifeForm[i]=="Phanerophyte" && (Hmax[i]>=300)) coef[i] = 0.719;
        else coef[i] = 0.95; 
      }
    }
  }
  return(coef);
}

/** Allometric coefficient retrieval with imputation */
NumericVector treeAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName) {
  NumericVector coef = speciesNumericParameter(SP,SpParams, parName);
  CharacterVector group = speciesCharacterParameter(SP, SpParams, "Group");
  for(int i=0;i<coef.size();i++) { // From De Caceres et al. 2019
    if(NumericVector::is_na(coef[i])) {
      if(parName=="a_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = 0.1989;
        else coef[i] = 0.0709;
      }
      else if(parName=="b_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = 1.3805;
        else coef[i] = 1.5120;
      }
      else if(parName=="c_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = -0.0363;
        else coef[i] = -0.0183;
      }
      else if(parName=="d_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = 0.066;
        else coef[i] = 0.0057;
      }
      else if(parName=="a_cw") {
        if(group[i]=="Gymnosperm") coef[i] = 0.747;
        else coef[i] = 0.839;
      }
      else if(parName=="b_cw") {
        if(group[i]=="Gymnosperm") coef[i] = 0.672;
        else coef[i] = 0.735;
      }
      else if(parName=="a_cr") {
        if(group[i]=="Gymnosperm") coef[i] = 1.995;
        else coef[i] = 1.506;
      }
      else if(parName=="b_1cr") {
        if(group[i]=="Gymnosperm") coef[i] = -0.649;
        else coef[i] = -0.706;
      }
      else if(parName=="b_2cr") {
        if(group[i]=="Gymnosperm") coef[i] = -0.020;
        else coef[i] = -0.078;
      }
      else if(parName=="b_3cr") {
        if(group[i]=="Gymnosperm") coef[i] = -0.00012;
        else coef[i] = 0.00018;
      }
      else if(parName=="c_1cr") {
        if(group[i]=="Gymnosperm") coef[i] = -0.004;
        else coef[i] = -0.007;
      }
      else if(parName=="c_2cr") {
        if(group[i]=="Gymnosperm") coef[i] = -0.159;
        else coef[i] = 0.0;
      }
    }
  }
  return(coef);
}





// [[Rcpp::export("species_parameter")]]
NumericVector speciesNumericParameterWithImputation(IntegerVector SP, DataFrame SpParams, String parName, bool fillMissing = true){
  if(fillMissing) {
    if(parName == "kPAR") return(kPARWithImputation(SP,SpParams));
    else if(parName == "gammaSWR") return(gammaSWRWithImputation(SP,SpParams));
    else if(parName == "alphaSWR") return(alphaSWRWithImputation(SP,SpParams));
    else if(parName == "g") return(gWithImputation(SP,SpParams));
    else if(parName == "r635") return(fineFoliarRatioWithImputation(SP,SpParams));
    else if(parName == "SLA") return(specificLeafAreaWithImputation(SP,SpParams));
    else if(parName == "LigninPercent") return(ligninPercentWithImputation(SP, SpParams));
    else if(parName == "SAV") return(surfaceToAreaRatioWithImputation(SP, SpParams));
    else if(parName == "HeatContent") return(heatContentWithImputation(SP, SpParams));
    else if(parName == "LeafWidth") return(leafWidthWithImputation(SP, SpParams));
    else if(parName == "Al2As") return(Al2AsWithImputation(SP, SpParams));
    else if(parName == "WoodDensity") return(woodDensityWithImputation(SP, SpParams));
    else if(parName == "LeafDensity") return(leafDensityWithImputation(SP, SpParams));
    else if(parName == "FineRootDensity") return(fineRootDensityWithImputation(SP, SpParams));
    else if(parName == "SRL") return(specificRootLengthWithImputation(SP, SpParams));
    else if(parName == "RLD") return(rootLengthDensityWithImputation(SP, SpParams));
    else if(parName == "conduit2sapwood") return(conduit2sapwoodWithImputation(SP, SpParams));
    else if(parName == "StemPI0") return(stemPI0WithImputation(SP, SpParams));
    else if(parName == "StemEPS") return(stemEPSWithImputation(SP, SpParams));
    else if(parName == "StemAF") return(stemAFWithImputation(SP, SpParams));
    else if(parName == "LeafPI0") return(leafPI0WithImputation(SP, SpParams));
    else if(parName == "LeafEPS") return(leafEPSWithImputation(SP, SpParams));
    else if(parName == "LeafAF") return(leafAFWithImputation(SP, SpParams));
  }
  return(speciesNumericParameter(SP, SpParams,parName));
}


// [[Rcpp::export("plant_parameter")]]
NumericVector cohortNumericParameterWithImputation(List x, DataFrame SpParams, String parName, bool fillMissing = true){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  IntegerVector tSP = treeData["Species"];
  IntegerVector shSP = shrubData["Species"];
  NumericVector par(tSP.size()+shSP.size());
  NumericVector parTrees = speciesNumericParameterWithImputation(tSP, SpParams, parName, fillMissing);
  NumericVector parShrubs = speciesNumericParameterWithImputation(shSP, SpParams, parName, fillMissing);
  for(int i=0;i<tSP.size();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shSP.size();i++) {
    par[i+tSP.size()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x);
  return(par);
}
