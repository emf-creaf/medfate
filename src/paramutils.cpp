#define STRICT_R_HEADERS
#include <numeric>
#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
using namespace Rcpp;

int findRowIndex(int sp, DataFrame SpParams) {
  IntegerVector spIndexSP = SpParams["SpIndex"];
  for(int i=0;i<spIndexSP.length();i++) if(spIndexSP[i]==sp) return(i);
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

// [[Rcpp::export("species_parameter")]]
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

/** Parameter retrieval with imputation */
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

/** Parameter retrieval with imputation */
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
