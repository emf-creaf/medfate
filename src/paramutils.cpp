#define STRICT_R_HEADERS
#include <numeric>
#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "forestutils.h"
#include "hydraulics.h"
#include "tissuemoisture.h"
using namespace Rcpp;

int findSpParamsRowByName(String spname, DataFrame SpParams) {
  CharacterVector spNameSP = SpParams["Name"];
  for(int i=0;i<spNameSP.length();i++) if(spNameSP[i]==spname) return(i);
  String s = "Species name '";
  s += spname;
  s +="' not found in SpParams";
  stop(s);
  return(NA_INTEGER);
}

int findSpParamsRowBySpIndex(int sp, DataFrame SpParams) {
  IntegerVector spIndexSP = SpParams["SpIndex"];
  for(int i=0;i<spIndexSP.length();i++) if(spIndexSP[i]==sp) return(i);
  String s = "Species index '";
  s += sp;
  s +="' not found in SpParams";
  stop(s);
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

IntegerVector speciesIndex(CharacterVector species, DataFrame SpParams){
  IntegerVector spIndex(species.size(), NA_INTEGER);
  IntegerVector spIndexSP = Rcpp::as<Rcpp::IntegerVector>(SpParams["SpIndex"]);
  for(int i=0;i<species.size();i++) {
    int iSP = findSpParamsRowByName(species[i], SpParams);
    spIndex[i] = spIndexSP[iSP];
  }
  return(spIndex);
}

// [[Rcpp::export(".speciesNumericParameterFromSpIndex")]]
NumericVector speciesNumericParameterFromIndex(IntegerVector SP, DataFrame SpParams, String parName){
  NumericVector par(SP.size(), NA_REAL);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    NumericVector parSP = Rcpp::as<Rcpp::NumericVector>(SpParams[parName]);
    for(int i=0;i<SP.size();i++) {
      int iSP = findSpParamsRowBySpIndex(SP[i], SpParams);
      par[i] = parSP[iSP];
    }
  } else {
    Rcerr << "Variable '" << parName.get_cstring() << "' was not found in SpParams!\n";
  }
  return(par);
}

NumericVector speciesNumericParameter(CharacterVector species, DataFrame SpParams, String parName){
  NumericVector par(species.size(), NA_REAL);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    NumericVector parSP = Rcpp::as<Rcpp::NumericVector>(SpParams[parName]);
    for(int i=0;i<species.size();i++) {
      int iSP = findSpParamsRowByName(species[i], SpParams);
      par[i] = parSP[iSP];
    }
  } else {
    Rcerr << "Variable '" << parName.get_cstring() << "' was not found in SpParams!\n";
  }
  return(par);
}


// [[Rcpp::export(".speciesCharacterParameterFromSpIndex")]]
CharacterVector speciesCharacterParameterFromIndex(IntegerVector SP, DataFrame SpParams, String parName){
  CharacterVector par(SP.size(), NA_STRING);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    CharacterVector parSP = SpParams[parName];
    for(int i=0;i<SP.size();i++) {
      int iSP = findSpParamsRowBySpIndex(SP[i], SpParams);
      par[i] = parSP[iSP];
    }
  } else {
    Rcerr << "Variable '" << parName.get_cstring() << "' was not found in SpParams!\n";
  }
  return(par);
}

//' @rdname species_values
// [[Rcpp::export("species_characterParameter")]]
 CharacterVector speciesCharacterParameter(CharacterVector species, DataFrame SpParams, String parName){
  CharacterVector par(species.size(), NA_STRING);
  if(SpParams.containsElementNamed(parName.get_cstring())) {
    CharacterVector parSP = SpParams[parName];
    for(int i=0;i<species.size();i++) {
      int iSP = findSpParamsRowByName(species[i], SpParams);
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
  NumericVector par(treeData.nrow() + shrubData.nrow());
  NumericVector parTrees, parShrubs;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    IntegerVector tSP = treeData["Species"];
    parTrees = speciesNumericParameterFromIndex(tSP, SpParams, parName);
  } else {
    CharacterVector tspecies = treeData["Species"];
    parTrees = speciesNumericParameter(tspecies, SpParams, parName);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    IntegerVector shSP = shrubData["Species"];
    parShrubs = speciesNumericParameterFromIndex(shSP, SpParams, parName);
  } else {
    CharacterVector sspecies = shrubData["Species"];
    parShrubs = speciesNumericParameter(sspecies, SpParams, parName);
  }
  for(int i=0;i<treeData.nrow();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shrubData.nrow();i++) {
    par[i + treeData.nrow()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x, SpParams);
  return(par);
}

//' @rdname plant_values
// [[Rcpp::export("plant_characterParameter")]]
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  CharacterVector par(treeData.nrow()+shrubData.nrow());
  CharacterVector parTrees, parShrubs;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    IntegerVector tSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
    parTrees = speciesCharacterParameterFromIndex(tSP, SpParams, parName);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    parTrees = speciesCharacterParameter(tspecies, SpParams, parName);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    IntegerVector shSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);
    parShrubs = speciesCharacterParameterFromIndex(shSP, SpParams, parName);
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    parShrubs = speciesCharacterParameter(sspecies, SpParams, parName);
  }
  for(int i=0;i<treeData.nrow();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shrubData.nrow();i++) {
    par[i + treeData.nrow()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x, SpParams);
  return(par);
}

/** Parameter retrieval with imputation */
NumericVector kPARWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  NumericVector kPAR = speciesNumericParameterFromIndex(SP, SpParams, "kPAR");
  for(int j=0;j<kPAR.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.55;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.45;
    } else if((leafShape[j]=="Needle") || (leafShape[j]=="Scale")){
      if(NumericVector::is_na(kPAR[j])) kPAR[j] = 0.50;
    }
  }
  return(kPAR);
}
NumericVector gammaSWRWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  NumericVector gammaSWR = speciesNumericParameterFromIndex(SP, SpParams, "gammaSWR");
  for(int j=0;j<gammaSWR.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.18;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.15;
    } else if((leafShape[j]=="Needle") || (leafShape[j]=="Scale")){
      if(NumericVector::is_na(gammaSWR[j])) gammaSWR[j] = 0.14;
    }
  }
  return(gammaSWR);
}
NumericVector alphaSWRWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector alphaSWR = speciesNumericParameterFromIndex(SP, SpParams, "alphaSWR");
  for(int j=0;j<alphaSWR.size();j++) {
    if(NumericVector::is_na(alphaSWR[j])) alphaSWR[j] = 0.7;
  }
  return(alphaSWR);
}
NumericVector gWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  NumericVector g = speciesNumericParameterFromIndex(SP, SpParams, "g");
  for(int j=0;j<g.size();j++) {
    if(leafShape[j] == "Broad") {
      if(NumericVector::is_na(g[j])) g[j] = 0.5;
    } else if(leafShape[j]=="Linear") {
      if(NumericVector::is_na(g[j])) g[j] = 0.8;
    } else if((leafShape[j]=="Needle") || (leafShape[j]=="Scale")){
      if(NumericVector::is_na(g[j])) g[j] = 1.0;
    }
  }
  return(g);
}
NumericVector fineFoliarRatioWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector ffr = speciesNumericParameterFromIndex(SP, SpParams, "r635");
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
  CharacterVector LeafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  CharacterVector LeafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  NumericVector SLA = speciesNumericParameterFromIndex(SP, SpParams, "SLA");
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
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector cohLigninPercent = speciesNumericParameterFromIndex(SP, SpParams, "LigninPercent");
  for(int i=0;i<cohLigninPercent.size();i++) {
    if(NumericVector::is_na(cohLigninPercent[i])) {
      if(leafShape[i]=="Scale") {
        cohLigninPercent[i] = 14.55;
      } else if(leafShape[i]=="Spines") {
        cohLigninPercent[i] = 14.55;
      } else if((leafShape[i]=="Linear") || (leafShape[i]=="Needle")) {
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
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector cohSAV = speciesNumericParameterFromIndex(SP, SpParams, "SAV");
  for(int i=0;i<cohSAV.size();i++) {
    if(NumericVector::is_na(cohSAV[i])) {
      if(leafShape[i]=="Scale") {
        cohSAV[i] = 1120.0;
      } else if(leafShape[i]=="Spines") {
        cohSAV[i] = 6750.0;
      } else if((leafShape[i]=="Linear") || (leafShape[i]=="Needle")) {
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
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector cohHeatContent = speciesNumericParameterFromIndex(SP, SpParams, "HeatContent");
  for(int i=0;i<cohHeatContent.size();i++) {
    if(NumericVector::is_na(cohHeatContent[i])) {
      if(leafShape[i]=="Scale") {
        cohHeatContent[i] = 20504.0;
      } else if(leafShape[i]=="Spines") {
        cohHeatContent[i] = 20433.0;
      } else if((leafShape[i]=="Linear") || (leafShape[i]=="Needle")) {
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
NumericVector proportionDeadWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector pDead = speciesNumericParameterFromIndex(SP, SpParams, "pDead");
  for(int c=0;c<pDead.size();c++) {
    if(NumericVector::is_na(pDead[c])) {
      pDead[c] = 0.05;
    }
  }
  return(pDead);
}
NumericVector leafWidthWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector leafwidth = speciesNumericParameterFromIndex(SP, SpParams, "LeafWidth");
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
NumericVector Ar2AlWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Ar2Al = speciesNumericParameterFromIndex(SP, SpParams, "Ar2Al");
  for(int c=0;c<Ar2Al.size();c++) {
    if(NumericVector::is_na(Ar2Al[c])) {
      Ar2Al[c] = 1.0;
    }
  }
  return(Ar2Al);
}
NumericVector Al2AsWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector Al2As = speciesNumericParameterFromIndex(SP, SpParams, "Al2As");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_Al2As = TFM["Al2As"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Al2As.size();c++) {
    if(NumericVector::is_na(Al2As[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Al2As[c] = fam_Al2As[i];
        }
      }
    }
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
  NumericVector woodDensity = speciesNumericParameterFromIndex(SP, SpParams, "WoodDensity");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_wood_density = TFM["WoodDensity"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<woodDensity.size();c++) {
    if(NumericVector::is_na(woodDensity[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          woodDensity[c] = fam_wood_density[i];
        }
      }
    }
    if(NumericVector::is_na(woodDensity[c])) {
      woodDensity[c] = 0.652;
    }
  }
  return(woodDensity);
}
NumericVector leafDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafDensity = speciesNumericParameterFromIndex(SP, SpParams, "LeafDensity");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_leaf_density = TFM["LeafDensity"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<leafDensity.size();c++) {
    if(NumericVector::is_na(leafDensity[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          leafDensity[c] = fam_leaf_density[i];
        }
      }
    }
    if(NumericVector::is_na(leafDensity[c])) {
      leafDensity[c] = 0.7;
    }
  }
  return(leafDensity);
}
NumericVector fineRootDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector fineRootDensity = speciesNumericParameterFromIndex(SP, SpParams, "FineRootDensity");
  for(int c=0;c<fineRootDensity.size();c++) {
    if(NumericVector::is_na(fineRootDensity[c])) {
      fineRootDensity[c] = 0.165;
    }
  }
  return(fineRootDensity);
}
NumericVector specificRootLengthWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector specificRootLength = speciesNumericParameterFromIndex(SP, SpParams, "SRL");
  for(int c=0;c<specificRootLength.size();c++) {
    if(NumericVector::is_na(specificRootLength[c])) {
      specificRootLength[c] = 3870.0;
    }
  }
  return(specificRootLength);
}
NumericVector rootLengthDensityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector rootLengthDensity = speciesNumericParameterFromIndex(SP, SpParams, "RLD");
  for(int c=0;c<rootLengthDensity.size();c++) {
    if(NumericVector::is_na(rootLengthDensity[c])) {
      rootLengthDensity[c] = 10.0;
    }
  }
  return(rootLengthDensity);
}
NumericVector conduit2sapwoodWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector Group = speciesCharacterParameterFromIndex(SP, SpParams, "Group");
  NumericVector conduit2sapwood = speciesNumericParameterFromIndex(SP, SpParams, "conduit2sapwood");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_conduit2sapwood = TFM["conduit2sapwood"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<conduit2sapwood.size();c++) {
    if(NumericVector::is_na(conduit2sapwood[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          conduit2sapwood[c] = fam_conduit2sapwood[i];
        }
      }
    }
    if(NumericVector::is_na(conduit2sapwood[c])) {
      if(Group[c]=="Angiosperm") conduit2sapwood[c] = 0.70; //20-40% parenchyma in angiosperms.
      else conduit2sapwood[c] = 0.925; //5-10% parenchyma in gymnosperms (https://link.springer.com/chapter/10.1007/978-3-319-15783-2_8)
    }
  }
  return(conduit2sapwood);
}
NumericVector maxFMCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector maxFMC = speciesNumericParameterFromIndex(SP, SpParams, "maxFMC");
  for(int c=0;c<maxFMC.size();c++) {
    if(NumericVector::is_na(maxFMC[c])) {
      maxFMC[c] = 120.0; //To be improved
    }
  }
  return(maxFMC);
}
NumericVector minFMCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector minFMC = speciesNumericParameterFromIndex(SP, SpParams, "minFMC");
  for(int c=0;c<minFMC.size();c++) {
    if(NumericVector::is_na(minFMC[c])) {
      minFMC[c] = 80.0; //To be improved
    }
  }
  return(minFMC);
}
NumericVector stemPI0WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WoodDensity = woodDensityWithImputation(SP, SpParams);
  NumericVector StemPI0 = speciesNumericParameterFromIndex(SP, SpParams, "StemPI0");
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
  NumericVector StemEPS = speciesNumericParameterFromIndex(SP, SpParams, "StemEPS");
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
  NumericVector stemAF = speciesNumericParameterFromIndex(SP, SpParams, "StemAF");
  for(int c=0;c<stemAF.size();c++) {
    if(NumericVector::is_na(stemAF[c])) {
      stemAF[c] = conduit2sapwood[c]; 
    }
  }
  return(stemAF);
}
NumericVector leafPI0WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafPI0 = speciesNumericParameterFromIndex(SP, SpParams, "LeafPI0");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_leaf_PI0 = TFM["LeafPI0"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<leafPI0.size();c++) {
    if(NumericVector::is_na(leafPI0[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          leafPI0[c] = fam_leaf_PI0[i];
        }
      }
    }
    if(NumericVector::is_na(leafPI0[c])) {
      //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
      leafPI0[c] = -2.0;//Average for Mediterranean climate species
    }
  }
  return(leafPI0);
}
NumericVector leafEPSWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafEPS = speciesNumericParameterFromIndex(SP, SpParams, "LeafEPS");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_leaf_EPS = TFM["LeafEPS"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<leafEPS.size();c++) {
    if(NumericVector::is_na(leafEPS[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          leafEPS[c] = fam_leaf_EPS[i];
        }
      }
    }
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(leafEPS[c])) {
      leafEPS[c] = 17.0;//Average for Mediterranean climate species
    }
  }
  return(leafEPS);
}
NumericVector leafAFWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafAF = speciesNumericParameterFromIndex(SP, SpParams, "LeafAF");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_leaf_AF = TFM["LeafAF"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<leafAF.size();c++) {
    if(NumericVector::is_na(leafAF[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          leafAF[c] = fam_leaf_AF[i];
        }
      }
    }
    //From: Bartlett MK, Scoffoni C, Sack L (2012) The determinants of leaf turgor loss point and prediction of drought tolerance of species and biomes: a global meta-analysis. Ecol Lett 15:393–405. doi: 10.1111/j.1461-0248.2012.01751.x
    if(NumericVector::is_na(leafAF[c])) {
      leafAF[c] = 0.29; //Average for Mediterranean climate species
    }
  }
  return(leafAF);
}

NumericVector TmaxLAIWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Tmax_LAI = speciesNumericParameterFromIndex(SP, SpParams, "Tmax_LAI");
  for(int c=0;c<Tmax_LAI.size();c++) {
    if(NumericVector::is_na(Tmax_LAI[c])) {
      Tmax_LAI[c] = 0.134;//Granier coefficient for LAI
    }
  }
  return(Tmax_LAI);
}
NumericVector TmaxLAIsqWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Tmax_LAIsq = speciesNumericParameterFromIndex(SP, SpParams, "Tmax_LAIsq");
  for(int c=0;c<Tmax_LAIsq.size();c++) {
    if(NumericVector::is_na(Tmax_LAIsq[c])) {
      Tmax_LAIsq[c] = -0.006; //Granier coefficient for LAI^2
    }
  }
  return(Tmax_LAIsq);
}
NumericVector WUEWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector leafShape = speciesCharacterParameterFromIndex(SP, SpParams, "LeafShape");
  CharacterVector leafSize = speciesCharacterParameterFromIndex(SP, SpParams, "LeafSize");
  NumericVector WUE = speciesNumericParameterFromIndex(SP, SpParams, "WUE");
  for(int c=0;c<WUE.size();c++) {
    if(NumericVector::is_na(WUE[c])) {
      WUE[c] = 7.9; 
    }
  }
  //Access internal data frame "trait_family_means"
  // Environment pkg = Environment::namespace_env("medfate");
  // DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  // CharacterVector fams = TFM.attr("row.names");
  // NumericVector fam_WUE = TFM["WUE"];
  // CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  // for(int c=0;c<WUE.size();c++) {
  //   if(NumericVector::is_na(WUE[c])) {
  //     for(int i=0;i<fams.size();i++) {
  //       if(fams[i]==family[c]) {
  //         WUE[c] = fam_WUE[i];
  //       }
  //     }
  //   }
  //   if(NumericVector::is_na(WUE[c])) {
  //     if(leafShape[c]=="Linear") {
  //       WUE[c]= 3.707131;
  //     } else if(leafShape[c]=="Needle") {
  //       WUE[c]= 3.707131;
  //     } else if(leafShape[c]=="Broad") {
  //       if(leafSize[c]=="Small") {
  //         WUE[c] = 4.289629;
  //       } else if(leafSize[c]=="Medium") {
  //         WUE[c] = 3.982086;
  //       } else if(leafSize[c]=="Large") {
  //         WUE[c]= 3.027647;
  //       }
  //     } else if(leafShape[c]=="Scale") { 
  //       WUE[c] = 1.665034;
  //     }
  //   }
  // }
  return(WUE);
}
NumericVector WUEPARWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WUE_par = speciesNumericParameterFromIndex(SP, SpParams, "WUE_par");
  for(int c=0;c<WUE_par.size();c++) {
    if(NumericVector::is_na(WUE_par[c])) {
      WUE_par[c] = 0.3643; //default value
    }
  }
  return(WUE_par);
}
NumericVector WUECO2WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WUE_co2 = speciesNumericParameterFromIndex(SP, SpParams, "WUE_co2");
  for(int c=0;c<WUE_co2.size();c++) {
    if(NumericVector::is_na(WUE_co2[c])) {
      WUE_co2[c] = 0.002757;
    }
  }
  return(WUE_co2);
}
NumericVector WUEVPDWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WUE_vpd = speciesNumericParameterFromIndex(SP, SpParams, "WUE_vpd");
  for(int c=0;c<WUE_vpd.size();c++) {
    if(NumericVector::is_na(WUE_vpd[c])) {
      WUE_vpd[c] = -0.4636;
    }
  }
  return(WUE_vpd);
}
NumericVector psi50Imputation(NumericVector psi50, IntegerVector SP, DataFrame SpParams) {
  CharacterVector Group = speciesCharacterParameterFromIndex(SP, SpParams, "Group");
  CharacterVector GrowthForm = speciesCharacterParameterFromIndex(SP, SpParams, "GrowthForm");
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_P50 = TFM["P50"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<psi50.size();c++) {
    if(NumericVector::is_na(psi50[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          psi50[c] = fam_P50[i];
        }
      }
    }
    if(NumericVector::is_na(psi50[c])) {
      // From: Maherali H, Pockman W, Jackson R (2004) Adaptive variation in the vulnerability of woody plants to xylem cavitation. Ecology 85:2184–2199
      if(Group[c]=="Angiosperm") {
        if((GrowthForm[c]=="Shrub") && (phenoType[c] != "winter-deciduous") && (phenoType[c] != "winter-semideciduous")) {
          psi50[c] = -5.09; //Angiosperm evergreen shrub
        } else if((GrowthForm[c]!="Shrub") && ((phenoType[c] == "winter-deciduous") || (phenoType[c] == "winter-semideciduous"))) {
          psi50[c] = -2.34; //Angiosperm winter-deciduous tree
        } else { 
          psi50[c] = -1.51; //Angiosperm evergreen tree
        }
      } else {
        if(GrowthForm[c]=="Shrub") {
          psi50[c] = -8.95; //Gymnosperm shrub
        } else {
          psi50[c] = -4.17; //Gymnosperm tree
        }
      }
    }
  }
  return(psi50);
}

NumericVector expExtractWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Exp_Extract = speciesNumericParameterFromIndex(SP, SpParams, "Exp_Extract");
  for(int c=0;c<Exp_Extract.size();c++) {
    if(NumericVector::is_na(Exp_Extract[c])) {
      Exp_Extract[c] = 1.3;
    }
  }
  return(Exp_Extract);
}

NumericVector psiExtractWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector leafPI0 = leafPI0WithImputation(SP, SpParams);
  NumericVector leafEPS = leafEPSWithImputation(SP, SpParams);
  NumericVector Exp_Extract = expExtractWithImputation(SP, SpParams);
  NumericVector Psi_Extract = speciesNumericParameterFromIndex(SP, SpParams, "Psi_Extract");
  for(int c=0;c<Psi_Extract.size();c++) {
    if(NumericVector::is_na(Psi_Extract[c])) {
      double corr = pow(log(0.5)/log(0.10), 1.0/Exp_Extract[c]); //Assumes TLP corresponds to 10% stomatal conductance (Martin-StPaul 2017 Ecol. Lett)
      Psi_Extract[c] = corr*turgorLossPoint(leafPI0[c], leafEPS[c]);
    }
  }
  return(Psi_Extract);
}
NumericVector GswmaxWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Gswmax = speciesNumericParameterFromIndex(SP, SpParams, "Gswmax");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_gswmax = TFM["Gswmax"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  
  for(int c=0;c<Gswmax.size();c++) {
    if(NumericVector::is_na(Gswmax[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Gswmax[c] = fam_gswmax[i];
        }
      }
    }
    if(NumericVector::is_na(Gswmax[c])) Gswmax[c] = 0.200;
  }
  return(Gswmax);
}
NumericVector GswminWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Gswmin = speciesNumericParameterFromIndex(SP, SpParams, "Gswmin");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_gswmin = TFM["Gswmin"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Gswmin.size();c++) {
    if(NumericVector::is_na(Gswmin[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Gswmin[c] = fam_gswmin[i];
        }
      }
    }
    if(NumericVector::is_na(Gswmin[c])) Gswmin[c] = 0.0049;
  }
  return(Gswmin);
}

NumericVector KmaxStemXylemWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Kmax_stemxylem = speciesNumericParameterFromIndex(SP, SpParams, "Kmax_stemxylem");
  CharacterVector Group = speciesCharacterParameterFromIndex(SP, SpParams, "Group");
  CharacterVector GrowthForm = speciesCharacterParameterFromIndex(SP, SpParams, "GrowthForm");
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_Ks = TFM["Kmax_stemxylem"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Kmax_stemxylem.size();c++) {
    if(NumericVector::is_na(Kmax_stemxylem[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Kmax_stemxylem[c] = fam_Ks[i];
        }
      }
    }
    if(NumericVector::is_na(Kmax_stemxylem[c])) {
      // From: Maherali H, Pockman W, Jackson R (2004) Adaptive variation in the vulnerability of woody plants to xylem cavitation. Ecology 85:2184–2199
      if(Group[c]=="Angiosperm") {
        if((GrowthForm[c]=="Shrub") && ((phenoType[c] == "winter-deciduous") || (phenoType[c] == "winter-semideciduous"))) {
          Kmax_stemxylem[c] = 1.55; //Angiosperm deciduous shrub
        } else if(((GrowthForm[c]=="Tree") || (GrowthForm[c]=="Tree/Shrub")) && ((phenoType[c] == "winter-deciduous") || (phenoType[c] == "winter-semideciduous"))) {
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
  }
  return(Kmax_stemxylem);
}
NumericVector KmaxRootXylemWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Kmax_stemxylem = KmaxStemXylemWithImputation(SP, SpParams);
  NumericVector Kmax_rootxylem = speciesNumericParameterFromIndex(SP, SpParams, "Kmax_rootxylem");
  for(int c=0;c<Kmax_rootxylem.size();c++) {
    //Oliveras I, Martínez-Vilalta J, Jimenez-Ortiz T, et al (2003) Hydraulic architecture of Pinus halepensis, P . pinea and Tetraclinis articulata in a dune ecosystem of Eastern Spain. Plant Ecol 131–141
    if(NumericVector::is_na(Kmax_rootxylem[c])) Kmax_rootxylem[c] = 4.0*Kmax_stemxylem[c];
  }
  return(Kmax_rootxylem);
}
NumericVector VCleafkmaxWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCleaf_kmax = speciesNumericParameterFromIndex(SP, SpParams, "VCleaf_kmax");
  NumericVector Gswmax = GswmaxWithImputation(SP, SpParams);
  for(int c=0;c<VCleaf_kmax.size();c++) {
    //Franks, P. J. (2006). Higher rates of leaf gas exchange are associated with higher leaf hydrodynamic pressure gradients. Plant, Cell and Environment, 29(4), 584–592. https://doi.org/10.1111/j.1365-3040.2005.01434.x
    //Original coefficients were c=0.02 and m = 1.4 but did not match graphical representation
    if(NumericVector::is_na(VCleaf_kmax[c])) {
      VCleaf_kmax[c] = pow(Gswmax[c]/0.015,1.0/1.3);
    }
  }
  return(VCleaf_kmax);
}
NumericVector NleafWithImputation(IntegerVector SP, DataFrame SpParams){
  NumericVector Nleaf = speciesNumericParameterFromIndex(SP, SpParams, "Nleaf");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_Nleaf = TFM["Nleaf"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Nleaf.size();c++) {
    if(NumericVector::is_na(Nleaf[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Nleaf[c] = fam_Nleaf[i];
        }
      }
    }
    if(NumericVector::is_na(Nleaf[c])) Nleaf[c] = 20.088;
  }
  return(Nleaf);
}
NumericVector NsapwoodWithImputation(IntegerVector SP, DataFrame SpParams){
  NumericVector Nsapwood = speciesNumericParameterFromIndex(SP, SpParams, "Nsapwood");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_Nsapwood = TFM["Nsapwood"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Nsapwood.size();c++) {
    if(NumericVector::is_na(Nsapwood[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Nsapwood[c] = fam_Nsapwood[i];
        }
      }
    }
    if(NumericVector::is_na(Nsapwood[c])) Nsapwood[c] = 3.9791;
  }
  return(Nsapwood);
}
NumericVector NfinerootWithImputation(IntegerVector SP, DataFrame SpParams){
  NumericVector Nfineroot = speciesNumericParameterFromIndex(SP, SpParams, "Nfineroot");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_Nfineroot = TFM["Nfineroot"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<Nfineroot.size();c++) {
    if(NumericVector::is_na(Nfineroot[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          Nfineroot[c] = fam_Nfineroot[i];
        }
      }
    }
    if(NumericVector::is_na(Nfineroot[c])) Nfineroot[c] = 12.207;
  }
  return(Nfineroot);
}

NumericVector Vmax298WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector SLA = specificLeafAreaWithImputation(SP, SpParams);
  NumericVector Nleaf = NleafWithImputation(SP, SpParams);
  NumericVector Vmax298 = speciesNumericParameterFromIndex(SP, SpParams, "Vmax298");
  for(int c=0;c<Vmax298.size();c++) {
    if(NumericVector::is_na(Vmax298[c]))  {
      if(!NumericVector::is_na(SLA[c]) && !NumericVector::is_na(Nleaf[c]))  {
        //Walker AP, Beckerman AP, Gu L, et al (2014) The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: A meta-analysis and modeling study. Ecol Evol 4:3218–3235. doi: 10.1002/ece3.1173
        double lnN = log(Nleaf[c]/SLA[c]); // Narea (g · m-2) = Nmass/SLA = (mgN/gdry)*(1 gN/1000 mg N)*(1000 g dry/ 1 kg dry)* (kg dry / m2 )
        double lnSLA = log(SLA[c]/1000.0); //SLA in m2*g-1
        Vmax298[c] = exp(1.993 + 2.555*lnN - 0.372*lnSLA + 0.422*lnN*lnSLA);
      } else {
        Vmax298[c] = 100.0; //Default value of Vmax298 = 100.0
      }
    }
  }
  return(Vmax298);
}
NumericVector Jmax298WithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector Vmax298 = Vmax298WithImputation(SP, SpParams);
  NumericVector Jmax298 = speciesNumericParameterFromIndex(SP, SpParams, "Jmax298");
  for(int c=0;c<Jmax298.size();c++) {
    //Walker AP, Beckerman AP, Gu L, et al (2014) The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: A meta-analysis and modeling study. Ecol Evol 4:3218–3235. doi: 10.1002/ece3.1173
    if(NumericVector::is_na(Jmax298[c])) Jmax298[c] = exp(1.197 + 0.847*log(Vmax298[c])); 
  }
  return(Jmax298);
}
NumericVector VCstemDWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCstem_d = speciesNumericParameterFromIndex(SP, SpParams, "VCstem_d");
  NumericVector psi50(VCstem_d.size(), NA_REAL);
  NumericVector psi50Imp = psi50Imputation(psi50, SP, SpParams);
  for(int c=0;c<VCstem_d.size();c++) {
    if(NumericVector::is_na(VCstem_d[c])) {
      double psi88 = 1.2593*psi50Imp[c] - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50Imp[c], psi88);
      VCstem_d[c] = par["d"];
    }
  }
  return(VCstem_d);
}
NumericVector VCstemCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCstem_c = speciesNumericParameterFromIndex(SP, SpParams, "VCstem_c");
  NumericVector psi50(VCstem_c.size(), NA_REAL);
  NumericVector psi50Imp = psi50Imputation(psi50, SP, SpParams);
  for(int c=0;c<VCstem_c.size();c++) {
    if(NumericVector::is_na(VCstem_c[c])) {
      double psi88 = 1.2593*psi50Imp[c] - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50Imp[c], psi88);
      VCstem_c[c] = par["c"];
    }
  }
  return(VCstem_c);
}
NumericVector VCleafDWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCleaf_d = speciesNumericParameterFromIndex(SP, SpParams, "VCleaf_d");
  NumericVector leafPI0 = leafPI0WithImputation(SP, SpParams);
  NumericVector leafEPS = leafEPSWithImputation(SP, SpParams);
  for(int c=0;c<VCleaf_d.size();c++) {
    if(NumericVector::is_na(VCleaf_d[c])) {
      double leaf_tlp = turgorLossPoint(leafPI0[c], leafEPS[c]);
      //From Bartlett,et al (2016). The correlations and sequence of plant stomatal, hydraulic, and wilting responses to drought. Proceedings of the National Academy of Sciences of the United States of America, 113(46), 13098–13103. https://doi.org/10.1073/pnas.1604088113
      double leaf_P50 = std::min(0.0, 0.9944*leaf_tlp + 0.2486);
      double leaf_P88 = 1.2593*leaf_P50 - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(leaf_P50, leaf_P88);
      VCleaf_d[c] = par["d"];
    }
  }
  return(VCleaf_d);
}
NumericVector VCleafCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCleaf_c = speciesNumericParameterFromIndex(SP, SpParams, "VCleaf_c");
  NumericVector leafPI0 = leafPI0WithImputation(SP, SpParams);
  NumericVector leafEPS = leafEPSWithImputation(SP, SpParams);
  for(int c=0;c<VCleaf_c.size();c++) {
    if(NumericVector::is_na(VCleaf_c[c])) {
      double leaf_tlp = turgorLossPoint(leafPI0[c], leafEPS[c]);
      //From Bartlett,et al (2016). The correlations and sequence of plant stomatal, hydraulic, and wilting responses to drought. Proceedings of the National Academy of Sciences of the United States of America, 113(46), 13098–13103. https://doi.org/10.1073/pnas.1604088113
      double leaf_P50 = std::min(0.0, 0.9944*leaf_tlp + 0.2486);
      double leaf_P88 = 1.2593*leaf_P50 - 1.4264; //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(leaf_P50, leaf_P88);
      VCleaf_c[c] = par["c"];
    }
  }
  return(VCleaf_c);
}
NumericVector VCrootDWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCroot_d = speciesNumericParameterFromIndex(SP, SpParams, "VCroot_d");
  NumericVector VCstem_d = VCstemDWithImputation(SP, SpParams);
  NumericVector VCstem_c = VCstemCWithImputation(SP, SpParams);
  for(int c=0;c<VCroot_d.size();c++) {
    if(NumericVector::is_na(VCroot_d[c])) {
      double psi50stem = std::min(-0.25, VCstem_d[c]*pow(0.6931472,1.0/VCstem_c[c]));
      double psi50root = std::min(-0.25, 0.742*psi50stem + 0.4892); //Regression using data from Bartlett et al. 2016
      double psi88root = std::min(-0.5, 1.2593*psi50root - 1.4264); //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50root, psi88root);
      VCroot_d[c] = par["d"];
      VCroot_d[c] = std::min(-0.25, VCroot_d[c]);
      // Rcout<< c<< " d: "<< VCroot_d[c]<<"\n";
    }
  }
  return(VCroot_d);
}
NumericVector VCrootCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector VCroot_c = speciesNumericParameterFromIndex(SP, SpParams, "VCroot_c");
  NumericVector VCstem_d = VCstemDWithImputation(SP, SpParams);
  NumericVector VCstem_c = VCstemCWithImputation(SP, SpParams);
  for(int c=0;c<VCroot_c.size();c++) {
    if(NumericVector::is_na(VCroot_c[c])) {
      double psi50stem = std::min(-0.25, VCstem_d[c]*pow(0.6931472,1.0/VCstem_c[c]));
      double psi50root = std::min(-0.25, 0.742*psi50stem + 0.4892); //Regression using data from Bartlett et al. 2016
      double psi88root = std::min(-0.5, 1.2593*psi50root - 1.4264); //Regression using data from Choat et al. 2012
      NumericVector par = psi2Weibull(psi50root, psi88root);
      VCroot_c[c] = par["c"];
      // Rcout<< c<< " c: "<< VCroot_c[c]<<"\n";
    }
  }
  return(VCroot_c);
}
NumericVector LeafRespirationRateWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector RERleaf = speciesNumericParameterFromIndex(SP, SpParams, "RERleaf");
  NumericVector Nleaf = NleafWithImputation(SP, SpParams);
  for(int c=0;c<RERleaf.size();c++) {
    if(NumericVector::is_na(RERleaf[c])) {
      //Reich, P. B., M. G. Tjoelker, K. S. Pregitzer, I. J. Wright, J. Oleksyn, and J. L. Machado. 2008. Scaling of respiration to nitrogen in leaves, stems and roots of higher land plants. Ecology Letters 11:793–801.
      double Nleaf_mmol_g = Nleaf[c]/14.0;
      double RER_nmolCO2_g_s = pow(10.0, 0.691 + 1.639*log10(Nleaf_mmol_g)); //nmol CO2·g-1·s-1
      RERleaf[c] = (24.0*3600.0)*(RER_nmolCO2_g_s/6.0)*(1e-9)*180.156; // nmol CO2·g-1·s-1 to g gluc·g-1·d-1
    }
  }
  return(RERleaf);
}
NumericVector SapwoodRespirationRateWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector RERsapwood = speciesNumericParameterFromIndex(SP, SpParams, "RERsapwood");
  NumericVector Nsapwood = NsapwoodWithImputation(SP, SpParams);
  for(int c=0;c<RERsapwood.size();c++) {
    if(NumericVector::is_na(RERsapwood[c])) {
      //Reich, P. B., M. G. Tjoelker, K. S. Pregitzer, I. J. Wright, J. Oleksyn, and J. L. Machado. 2008. Scaling of respiration to nitrogen in leaves, stems and roots of higher land plants. Ecology Letters 11:793–801.
      // double Nsapwood_mmol_g = Nsapwood[c]/14.0;
      // double RER_nmolCO2_g_s = pow(10.0, 1.024 + 1.344*log10(Nsapwood_mmol_g)); //nmol CO2·g-1·s-1
      // RERsapwood[c] = 24.0*3600.0*(RER_nmolCO2_g_s/6.0)*(1e-9)*180.156; // nmol CO2·g-1·s-1 to g gluc·g-1·d-1
      // ESTIMATES ARE TOO HIGH
      RERsapwood[c] = 5.15e-05;
    }
  }
  return(RERsapwood);
}
NumericVector SapwoodSenescenceRateWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector SRsapwood = speciesNumericParameterFromIndex(SP, SpParams, "SRsapwood");
  for(int c=0;c<SRsapwood.size();c++) {
    if(NumericVector::is_na(SRsapwood[c])) {
      SRsapwood[c] = 0.000135;
    }
  }
  return(SRsapwood);
}

NumericVector FinerootRespirationRateWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector RERfineroot = speciesNumericParameterFromIndex(SP, SpParams, "RERfineroot");
  NumericVector Nfineroot = NsapwoodWithImputation(SP, SpParams);
  for(int c=0;c<RERfineroot.size();c++) {
    if(NumericVector::is_na(RERfineroot[c])) {
      //Reich, P. B., M. G. Tjoelker, K. S. Pregitzer, I. J. Wright, J. Oleksyn, and J. L. Machado. 2008. Scaling of respiration to nitrogen in leaves, stems and roots of higher land plants. Ecology Letters 11:793–801.
      double Nfineroot_mmol_g = Nfineroot[c]/14.0;
      double RER_nmolCO2_g_s = pow(10.0, 0.980 + 1.352*log10(Nfineroot_mmol_g)); //nmol CO2·g-1·s-1
      RERfineroot[c] = 24.0*3600.0*(RER_nmolCO2_g_s/6.0)*(1e-9)*180.156; // nmol CO2·g-1·s-1 to g gluc·g-1·d-1
    }
  }
  return(RERfineroot);
}

NumericVector WoodCWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector WoodC = speciesNumericParameterFromIndex(SP, SpParams, "WoodC");
  //Access internal data frame "trait_family_means"
  Environment pkg = Environment::namespace_env("medfate");
  DataFrame TFM = Rcpp::as<Rcpp::DataFrame>(pkg["trait_family_means"]);
  CharacterVector fams = TFM.attr("row.names");
  NumericVector fam_WoodC = TFM["WoodC"];
  CharacterVector family = speciesCharacterParameterFromIndex(SP, SpParams, "Family");
  for(int c=0;c<WoodC.size();c++) {
    if(NumericVector::is_na(WoodC[c])) {
      for(int i=0;i<fams.size();i++) {
        if(fams[i]==family[c]) {
          WoodC[c] = fam_WoodC[i];
        }
      }
    }
    if(NumericVector::is_na(WoodC[c])) {
      WoodC[c] = 0.5; 
    }
  }
  return(WoodC);
}
NumericVector leafDurationWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector leafDuration = speciesNumericParameterFromIndex(SP, SpParams, "LeafDuration");
  for(int c=0;c<leafDuration.size();c++) {
    if(NumericVector::is_na(leafDuration[c])) {
      if((phenoType[c]=="winter-deciduous") || (phenoType[c]=="winter-semideciduous")) leafDuration[c] = 1.0; 
      else leafDuration[c] = 2.42; //Average value
      // From Falster et al. (2018) 1/LeafDuration = b1*(LMA/LMA_0)^(-b2), with LMA_0 = 0.1978791, b1 = 0.4565655, b2 = 1.71
    }
  }
  return(leafDuration);
}
NumericVector t0gddWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector t0gdd = speciesNumericParameterFromIndex(SP, SpParams, "t0gdd");
  for(int c=0;c<t0gdd.size();c++) {
    if(NumericVector::is_na(t0gdd[c])) {
      t0gdd[c] = 50.0; //Default
    }
  }
  return(t0gdd);
}
NumericVector SgddWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector Sgdd = speciesNumericParameterFromIndex(SP, SpParams, "Sgdd");
  for(int c=0;c<Sgdd.size();c++) {
    if(NumericVector::is_na(Sgdd[c])) {
      Sgdd[c] = 200.0; //Default
    }
  }
  return(Sgdd);
}
NumericVector TbgddWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector Tbgdd = speciesNumericParameterFromIndex(SP, SpParams, "Tbgdd");
  for(int c=0;c<Tbgdd.size();c++) {
    if(NumericVector::is_na(Tbgdd[c])) {
      Tbgdd[c] = 0.0; //Default
    }
  }
  return(Tbgdd);
}
NumericVector SsenWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector Ssen = speciesNumericParameterFromIndex(SP, SpParams, "Ssen");
  for(int c=0;c<Ssen.size();c++) {
    if(NumericVector::is_na(Ssen[c])) {//Delpierre et al 2009
      Ssen[c] = 8268.0; //Default (broadleaved deciduous)
    }
  }
  return(Ssen);
}
NumericVector PhsenWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector Phsen = speciesNumericParameterFromIndex(SP, SpParams, "Phsen");
  for(int c=0;c<Phsen.size();c++) {
    if(NumericVector::is_na(Phsen[c])) {//Delpierre et al 2009
      Phsen[c] = 12.5; //Default (broadleaved deciduous)
    }
  }
  return(Phsen);
}
NumericVector TbsenWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector Tbsen = speciesNumericParameterFromIndex(SP, SpParams, "Tbsen");
  for(int c=0;c<Tbsen.size();c++) {
    if(NumericVector::is_na(Tbsen[c])) {//Delpierre et al 2009
      Tbsen[c] = 28.5; //Default (broadleaved deciduous)
    }
  }
  return(Tbsen);
}
NumericVector xsenWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector xsen = speciesNumericParameterFromIndex(SP, SpParams, "xsen");
  for(int c=0;c<xsen.size();c++) {
    if(NumericVector::is_na(xsen[c])) {//Delpierre et al 2009
      xsen[c] = 2.0; //Default
    }
  }
  return(xsen);
}
NumericVector ysenWithImputation(IntegerVector SP, DataFrame SpParams) {
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  NumericVector ysen = speciesNumericParameterFromIndex(SP, SpParams, "ysen");
  for(int c=0;c<ysen.size();c++) {
    if(NumericVector::is_na(ysen[c])) {//Delpierre et al 2009
      ysen[c] = 2.0; //Default
    }
  }
  return(ysen);
}
NumericVector seedMassWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector seedMass = speciesNumericParameterFromIndex(SP, SpParams, "SeedMass");
  for(int c=0;c<seedMass.size();c++) {
    if(NumericVector::is_na(seedMass[c])) {
      seedMass[c] = 50.0; //Default 50 mg
    }
  }
  return(seedMass);
}
NumericVector seedLongevityWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector seedLongevity = speciesNumericParameterFromIndex(SP, SpParams, "SeedLongevity");
  for(int c=0;c<seedLongevity.size();c++) {
    if(NumericVector::is_na(seedLongevity[c])) {
      seedLongevity[c] = 2.0; //Default 2 years
    }
  }
  return(seedLongevity);
}
NumericVector dispersalDistanceWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector dispersalDistance = speciesNumericParameterFromIndex(SP, SpParams, "DispersalDistance");
  for(int c=0;c<dispersalDistance.size();c++) {
    if(NumericVector::is_na(dispersalDistance[c])) {
      dispersalDistance[c] = 50.0; //Default 50 m
      //TO BE DONE: IMPLEMENT RULES in Bullock et al. 2019
    }
  }
  return(dispersalDistance);
}
NumericVector dispersalShapeWithImputation(IntegerVector SP, DataFrame SpParams) {
  NumericVector dispersalShape = speciesNumericParameterFromIndex(SP, SpParams, "DispersalShape");
  for(int c=0;c<dispersalShape.size();c++) {
    if(NumericVector::is_na(dispersalShape[c])) {
      dispersalShape[c] = 2.0; //Default 2 (bell shaped)
      //TO BE DONE: IMPLEMENT RULES in Bullock et al. 2019
    }
  }
  return(dispersalShape);
}

/** Allometric coefficient retrieval with imputation */
NumericVector shrubAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName) {
  NumericVector coef = speciesNumericParameterFromIndex(SP,SpParams, parName);
  CharacterVector GrowthForm = speciesCharacterParameterFromIndex(SP, SpParams, "GrowthForm");
  CharacterVector lifeForm = speciesCharacterParameterFromIndex(SP, SpParams, "LifeForm");
  NumericVector Hmax = speciesNumericParameterFromIndex(SP, SpParams, "Hmax");
  for(int i=0;i<coef.size();i++) { // From De Caceres et al. 2019
    if(NumericVector::is_na(coef[i]) && (GrowthForm[i] != "Tree")) {
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
      else if(parName=="BTsh") {
        coef[i] = 2.0;
      }
    }
  }
  return(coef);
}

/** Allometric coefficient retrieval with imputation */
NumericVector treeAllometricCoefficientWithImputation(IntegerVector SP, DataFrame SpParams, String parName) {
  NumericVector coef = speciesNumericParameterFromIndex(SP,SpParams, parName);
  CharacterVector GrowthForm = speciesCharacterParameterFromIndex(SP, SpParams, "GrowthForm");
  CharacterVector group = speciesCharacterParameterFromIndex(SP, SpParams, "Group");
  for(int i=0;i<coef.size();i++) { // From De Caceres et al. 2019
    //Limit to group values for 'c' 
    if(parName=="c_fbt") {
      if(group[i]=="Gymnosperm") coef[i] = std::min(coef[i],-0.0066);
      else coef[i] = std::min(coef[i],-0.0147);
    }
    if(parName=="d_fbt") {
      coef[i] = 0.0;
    }
    if(NumericVector::is_na(coef[i]) && (GrowthForm[i]!="Shrub")) {
      if(parName=="a_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = 0.0527;
        else coef[i] = 0.1300;
      }
      else if(parName=="b_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = 1.5782;
        else coef[i] = 1.2285;
      }
      else if(parName=="c_fbt") {
        if(group[i]=="Gymnosperm") coef[i] = -0.0066;
        else coef[i] = -0.0147;
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
      else if(parName=="a_bt") {
        coef[i] = 0.886;
      }
      else if(parName=="b_bt") {
        coef[i] = 0.969;
      }
    }
  }
  return(coef);
}


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
    else if(parName == "pDead") return(proportionDeadWithImputation(SP,SpParams));
    else if(parName == "LeafWidth") return(leafWidthWithImputation(SP, SpParams));
    else if(parName == "Ar2Al") return(Ar2AlWithImputation(SP, SpParams));
    else if(parName == "Al2As") return(Al2AsWithImputation(SP, SpParams));
    else if(parName == "WoodDensity") return(woodDensityWithImputation(SP, SpParams));
    else if(parName == "LeafDensity") return(leafDensityWithImputation(SP, SpParams));
    else if(parName == "FineRootDensity") return(fineRootDensityWithImputation(SP, SpParams));
    else if(parName == "SRL") return(specificRootLengthWithImputation(SP, SpParams));
    else if(parName == "RLD") return(rootLengthDensityWithImputation(SP, SpParams));
    else if(parName == "conduit2sapwood") return(conduit2sapwoodWithImputation(SP, SpParams));
    else if(parName == "maxFMC") return(maxFMCWithImputation(SP, SpParams));
    else if(parName == "minFMC") return(minFMCWithImputation(SP, SpParams));
    else if(parName == "StemPI0") return(stemPI0WithImputation(SP, SpParams));
    else if(parName == "StemEPS") return(stemEPSWithImputation(SP, SpParams));
    else if(parName == "StemAF") return(stemAFWithImputation(SP, SpParams));
    else if(parName == "LeafPI0") return(leafPI0WithImputation(SP, SpParams));
    else if(parName == "LeafEPS") return(leafEPSWithImputation(SP, SpParams));
    else if(parName == "LeafAF") return(leafAFWithImputation(SP, SpParams));
    else if(parName == "Tmax_LAI") return(TmaxLAIWithImputation(SP, SpParams));
    else if(parName == "Tmax_LAIsq") return(TmaxLAIsqWithImputation(SP, SpParams));
    else if(parName == "WUE") return(WUEWithImputation(SP, SpParams));
    else if(parName == "WUE_par") return(WUEPARWithImputation(SP, SpParams));
    else if(parName == "WUE_co2") return(WUECO2WithImputation(SP, SpParams));
    else if(parName == "WUE_vpd") return(WUEVPDWithImputation(SP, SpParams));
    else if(parName == "Psi_Extract") return(psiExtractWithImputation(SP, SpParams));
    else if(parName == "Exp_Extract") return(expExtractWithImputation(SP, SpParams));
    else if(parName == "Kmax_stemxylem") return(KmaxStemXylemWithImputation(SP, SpParams));
    else if(parName == "Kmax_rootxylem") return(KmaxRootXylemWithImputation(SP, SpParams));
    else if(parName == "VCleaf_kmax") return(VCleafkmaxWithImputation(SP, SpParams));
    else if(parName == "Gswmax") return(GswmaxWithImputation(SP, SpParams));
    else if(parName == "Gswmin") return(GswminWithImputation(SP, SpParams));
    else if(parName == "Nleaf") return(NleafWithImputation(SP, SpParams));
    else if(parName == "Nsapwood") return(NsapwoodWithImputation(SP, SpParams));
    else if(parName == "Nfineroot") return(NfinerootWithImputation(SP, SpParams));
    else if(parName == "RERleaf") return(LeafRespirationRateWithImputation(SP, SpParams));
    else if(parName == "RERsapwood") return(SapwoodRespirationRateWithImputation(SP, SpParams));
    else if(parName == "RERfineroot") return(FinerootRespirationRateWithImputation(SP, SpParams));
    else if(parName == "SRsapwood") return(SapwoodSenescenceRateWithImputation(SP, SpParams));
    else if(parName == "Vmax298") return(Vmax298WithImputation(SP, SpParams));
    else if(parName == "Jmax298") return(Jmax298WithImputation(SP, SpParams));
    else if(parName == "VCstem_c") return(VCstemCWithImputation(SP, SpParams));
    else if(parName == "VCstem_d") return(VCstemDWithImputation(SP, SpParams));
    else if(parName == "VCleaf_c") return(VCleafCWithImputation(SP, SpParams));
    else if(parName == "VCleaf_d") return(VCleafDWithImputation(SP, SpParams));
    else if(parName == "VCroot_c") return(VCrootCWithImputation(SP, SpParams));
    else if(parName == "VCroot_d") return(VCrootDWithImputation(SP, SpParams));
    else if(parName == "WoodC") return(WoodCWithImputation(SP, SpParams));
    else if(parName == "LeafDuration") return(leafDurationWithImputation(SP, SpParams));
    else if(parName == "t0gdd") return(t0gddWithImputation(SP, SpParams));
    else if(parName == "Sgdd") return(SgddWithImputation(SP, SpParams));
    else if(parName == "Tbgdd") return(TbgddWithImputation(SP, SpParams));
    else if(parName == "Ssen") return(SsenWithImputation(SP, SpParams));
    else if(parName == "Phsen") return(PhsenWithImputation(SP, SpParams));
    else if(parName == "Tbsen") return(TbsenWithImputation(SP, SpParams));
    else if(parName == "xsen") return(xsenWithImputation(SP, SpParams));
    else if(parName == "ysen") return(ysenWithImputation(SP, SpParams));
    else if(parName == "SeedMass") return(seedMassWithImputation(SP, SpParams));
    else if(parName == "SeedLongevity") return(seedLongevityWithImputation(SP, SpParams));
    else if(parName == "DispersalDistance") return(dispersalDistanceWithImputation(SP, SpParams));
    else if(parName == "DispersalShape") return(dispersalShapeWithImputation(SP, SpParams));
    else if((parName == "a_fbt") || (parName == "b_fbt") || (parName == "c_fbt") || (parName == "d_fbt")) return(treeAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_cw") || (parName == "b_cw")) return(treeAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_cr") || (parName == "b_1cr") || (parName == "b_2cr") || (parName == "b_3cr")) return(treeAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "c_1cr") || (parName == "c_2cr")) return(treeAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_bt") || (parName == "b_bt")) return(treeAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_ash") || (parName == "b_ash")) return(shrubAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_bsh") || (parName == "b_bsh")) return(shrubAllometricCoefficientWithImputation(SP, SpParams, parName));
    else if((parName == "a_btsh") || (parName == "b_btsh") || (parName == "cr") || (parName == "BTsh")) return(shrubAllometricCoefficientWithImputation(SP, SpParams, parName));
  }
  return(speciesNumericParameterFromIndex(SP, SpParams,parName));
}

//' @rdname species_values
// [[Rcpp::export("species_parameter")]]
NumericVector speciesNumericParameterWithImputation(CharacterVector species, DataFrame SpParams, String parName, bool fillMissing = true){
  if(fillMissing) {
    IntegerVector SP = speciesIndex(species, SpParams);
    return(speciesNumericParameterWithImputation(SP, SpParams, parName, fillMissing));
  }
  return(speciesNumericParameter(species, SpParams,parName));
}

//' @rdname plant_values
// [[Rcpp::export("plant_parameter")]]
NumericVector cohortNumericParameterWithImputation(List x, DataFrame SpParams, String parName, bool fillMissing = true){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  NumericVector par(treeData.nrow() + shrubData.nrow());
  NumericVector parTrees, parShrubs;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    IntegerVector tSP = treeData["Species"];
    parTrees = speciesNumericParameterWithImputation(tSP, SpParams, parName, fillMissing);
  } else {
    CharacterVector tspecies = treeData["Species"];
    parTrees = speciesNumericParameterWithImputation(tspecies, SpParams, parName, fillMissing);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    IntegerVector shSP = shrubData["Species"];
    parShrubs = speciesNumericParameterWithImputation(shSP, SpParams, parName, fillMissing);
  } else {
    CharacterVector sspecies = shrubData["Species"];
    parShrubs = speciesNumericParameterWithImputation(sspecies, SpParams, parName, fillMissing);
  }
  for(int i=0;i<treeData.nrow();i++) {
    par[i] = parTrees[i];
  }
  for(int i=0;i<shrubData.nrow();i++) {
    par[i + treeData.nrow()] = parShrubs[i];
  }
  par.attr("names") = cohortIDs(x, SpParams);
  return(par);
}
