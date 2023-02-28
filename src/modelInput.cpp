#include <Rcpp.h>
#include "spwb.h"
#include "growth.h"
#include "carbon.h"
#include "root.h"
#include "soil.h"
#include "woodformation.h"
#include "forestutils.h"
#include "paramutils.h"
#include "tissuemoisture.h"
#include "hydraulics.h"
#include "stdlib.h"

using namespace Rcpp;


DataFrame paramsPhenology(DataFrame above, DataFrame SpParams, bool fillMissingSpParams) {
  IntegerVector SP = above["SP"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_dead = above["LAI_dead"];
  int numCohorts = SP.size();
  
  NumericVector leafDuration  = speciesNumericParameterWithImputation(SP, SpParams, "LeafDuration", fillMissingSpParams);
  NumericVector t0gdd  = speciesNumericParameterWithImputation(SP, SpParams, "t0gdd", fillMissingSpParams);
  NumericVector Sgdd  = speciesNumericParameterWithImputation(SP, SpParams, "Sgdd", fillMissingSpParams);
  NumericVector Tbgdd = speciesNumericParameterWithImputation(SP, SpParams, "Tbgdd", fillMissingSpParams);
  NumericVector Ssen = speciesNumericParameterWithImputation(SP, SpParams, "Ssen", fillMissingSpParams);
  NumericVector Phsen = speciesNumericParameterWithImputation(SP, SpParams, "Phsen", fillMissingSpParams);
  NumericVector Tbsen  = speciesNumericParameterWithImputation(SP, SpParams, "Tbsen", fillMissingSpParams);
  NumericVector xsen  = speciesNumericParameterWithImputation(SP, SpParams, "xsen", fillMissingSpParams);
  NumericVector ysen  = speciesNumericParameterWithImputation(SP, SpParams, "ysen", fillMissingSpParams);
  
  CharacterVector phenoType = speciesCharacterParameterFromIndex(SP, SpParams, "PhenologyType");
  for(int j=0; j<numCohorts;j++) {
    if(phenoType[j] == "winter-deciduous" || phenoType[j] == "winter-semideciduous") { 
      LAI_expanded[j] = 0.0; //Set initial LAI to zero, assuming simulations start at Jan 1st
      if(phenoType[j] == "winter-semideciduous") LAI_dead[j] = LAI_live[j];
    }
    if(phenoType[j]=="oneflush-evergreen") {
      //Do not allow flushing all leaves at once (i.e. limit leaf duration to 1.25 yrs)
      leafDuration[j] = std::max(leafDuration[j], 1.25);
    }
  } 
  DataFrame paramsPhenologydf = DataFrame::create(
    _["PhenologyType"] = phenoType,
    _["LeafDuration"] = leafDuration,
    _["t0gdd"] = t0gdd,_["Sgdd"] = Sgdd, _["Tbgdd"] = Tbgdd, 
    _["Ssen"] = Ssen, _["Phsen"] = Phsen, _["Tbsen"] = Tbsen, _["xsen"] = xsen, _["ysen"] = ysen 
  );
  paramsPhenologydf.attr("row.names") = above.attr("row.names");
  return(paramsPhenologydf);
}

DataFrame paramsInterception(DataFrame above, DataFrame SpParams, List control) {
  IntegerVector SP = above["SP"];

  String transpirationMode = control["transpirationMode"];
  bool fillMissingSpParams = control["fillMissingSpParams"];
  
  NumericVector alphaSWR = speciesNumericParameterWithImputation(SP, SpParams, "alphaSWR", fillMissingSpParams);
  NumericVector gammaSWR = speciesNumericParameterWithImputation(SP, SpParams, "gammaSWR", fillMissingSpParams);
  NumericVector kPAR = speciesNumericParameterWithImputation(SP, SpParams, "kPAR", fillMissingSpParams);
  NumericVector g = speciesNumericParameterWithImputation(SP, SpParams, "g", fillMissingSpParams);
  DataFrame paramsInterceptiondf;
  if(transpirationMode=="Granier") {
    paramsInterceptiondf = DataFrame::create(_["kPAR"] = kPAR, 
                                             _["g"] = g);
  } else {
    paramsInterceptiondf = DataFrame::create(_["alphaSWR"] = alphaSWR,
                                             _["gammaSWR"] = gammaSWR, 
                                             _["kPAR"] = kPAR, 
                                             _["g"] = g);
  }
  paramsInterceptiondf.attr("row.names") = above.attr("row.names");
  return(paramsInterceptiondf);
}


DataFrame paramsAnatomy(DataFrame above, DataFrame SpParams, bool fillMissingSpParams, 
                        String model = "spwb", String transpirationMode = "Granier") {
  IntegerVector SP = above["SP"];

  NumericVector Hmax = speciesNumericParameterFromIndex(SP, SpParams, "Hmax");
  NumericVector Hmed = speciesNumericParameterFromIndex(SP, SpParams, "Hmed"); //To correct conductivity
  
  NumericVector Al2As = speciesNumericParameterWithImputation(SP, SpParams, "Al2As", fillMissingSpParams);
  NumericVector Ar2Al = speciesNumericParameterWithImputation(SP, SpParams, "Ar2Al", fillMissingSpParams);
  NumericVector SLA = speciesNumericParameterWithImputation(SP, SpParams, "SLA", fillMissingSpParams);
  NumericVector r635 = speciesNumericParameterWithImputation(SP, SpParams, "r635", fillMissingSpParams);
  NumericVector leafwidth = speciesNumericParameterWithImputation(SP, SpParams, "LeafWidth", fillMissingSpParams);
  NumericVector WoodDensity = speciesNumericParameterWithImputation(SP, SpParams, "WoodDensity", fillMissingSpParams);
  NumericVector LeafDensity = speciesNumericParameterWithImputation(SP, SpParams, "LeafDensity", fillMissingSpParams);
  NumericVector FineRootDensity = speciesNumericParameterWithImputation(SP, SpParams, "FineRootDensity", fillMissingSpParams);
  NumericVector SRL = speciesNumericParameterWithImputation(SP, SpParams, "SRL", fillMissingSpParams);  
  NumericVector RLD = speciesNumericParameterWithImputation(SP, SpParams, "RLD", fillMissingSpParams);  
  NumericVector conduit2sapwood = speciesNumericParameterWithImputation(SP, SpParams, "conduit2sapwood", fillMissingSpParams);

  DataFrame paramsAnatomydf;
  if(model=="spwb") {
    if(transpirationMode=="Sperry") {
      paramsAnatomydf = DataFrame::create(
        _["Hmed"] = Hmed,
        _["Al2As"] = Al2As, _["SLA"] = SLA, _["LeafWidth"] = leafwidth, 
        _["LeafDensity"] = LeafDensity, _["WoodDensity"] = WoodDensity, _["FineRootDensity"] = LeafDensity, 
          _["conduit2sapwood"] = conduit2sapwood,
          _["SRL"] = SRL, _["RLD"] = RLD,  
          _["r635"] = r635);
    } else {
      paramsAnatomydf = DataFrame::create(
          _["Al2As"] = Al2As, _["Ar2Al"] = Ar2Al, _["SLA"] = SLA,
          _["LeafDensity"] = LeafDensity, _["WoodDensity"] = WoodDensity, _["FineRootDensity"] = LeafDensity, 
          _["SRL"] = SRL, _["RLD"] = RLD,  
          _["r635"] = r635);
    }
  } else if(model=="growth") {
    if(transpirationMode=="Sperry") {
      paramsAnatomydf = DataFrame::create(
        _["Hmax"] = Hmax,_["Hmed"] = Hmed,
        _["Al2As"] = Al2As, _["SLA"] = SLA, _["LeafWidth"] = leafwidth, 
        _["LeafDensity"] = LeafDensity, _["WoodDensity"] = WoodDensity, _["FineRootDensity"] = LeafDensity, 
          _["conduit2sapwood"] = conduit2sapwood,
          _["SRL"] = SRL, _["RLD"] = RLD,  
          _["r635"] = r635);
    } else {
      paramsAnatomydf = DataFrame::create(
        _["Hmax"] = Hmax,_["Hmed"] = Hmed,
        _["Al2As"] = Al2As, _["Ar2Al"] = Ar2Al, _["SLA"] = SLA, _["LeafWidth"] = leafwidth, 
          _["LeafDensity"] = LeafDensity, _["WoodDensity"] = WoodDensity, _["FineRootDensity"] = LeafDensity, 
            _["conduit2sapwood"] = conduit2sapwood,
            _["SRL"] = SRL, _["RLD"] = RLD,  
            _["r635"] = r635);
    }
  }
  paramsAnatomydf.attr("row.names") = above.attr("row.names");
  return(paramsAnatomydf);
}

DataFrame paramsWaterStorage(DataFrame above, List belowLayers,
                             DataFrame SpParams,
                             DataFrame paramsAnatomydf, bool fillMissingSpParams) {
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  int numCohorts = SP.size();
  
  NumericMatrix V = belowLayers["V"];
  NumericMatrix L = belowLayers["L"];
  
  NumericVector StemPI0 = speciesNumericParameterWithImputation(SP, SpParams, "StemPI0", fillMissingSpParams);
  NumericVector StemEPS = speciesNumericParameterWithImputation(SP, SpParams, "StemEPS", fillMissingSpParams);
  NumericVector StemAF = speciesNumericParameterWithImputation(SP, SpParams, "StemAF", fillMissingSpParams);
  NumericVector LeafPI0 = speciesNumericParameterWithImputation(SP, SpParams, "LeafPI0", fillMissingSpParams);
  NumericVector LeafEPS = speciesNumericParameterWithImputation(SP, SpParams, "LeafEPS", fillMissingSpParams);
  NumericVector LeafAF = speciesNumericParameterWithImputation(SP, SpParams, "LeafAF", fillMissingSpParams);
  
  NumericVector maxFMC = speciesNumericParameterWithImputation(SP, SpParams, "maxFMC", fillMissingSpParams);
  
  NumericVector Vsapwood(numCohorts), Vleaf(numCohorts);
  
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector LeafDensity = paramsAnatomydf["LeafDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];

  //Calculate stem and leaf capacity per leaf area (in l·m-2)
  for(int c=0;c<numCohorts;c++){
    Vsapwood[c] = sapwoodWaterCapacity(Al2As[c], H[c], V, L, WoodDensity[c]); 
    Vleaf[c] = leafWaterCapacity(SLA[c], LeafDensity[c]); 
  }
  DataFrame paramsWaterStoragedf = DataFrame::create(_["maxFMC"] = maxFMC,
      _["LeafPI0"] = LeafPI0, _["LeafEPS"] = LeafEPS, _["LeafAF"] = LeafAF, _["Vleaf"] = Vleaf,
      _["StemPI0"] = StemPI0, _["StemEPS"] = StemEPS, _["StemAF"] = StemAF, _["Vsapwood"] = Vsapwood);
  paramsWaterStoragedf.attr("row.names") = above.attr("row.names");
  
  return(paramsWaterStoragedf);
}

DataFrame paramsTranspirationGranier(DataFrame above,  DataFrame SpParams, bool fillMissingSpParams) {
  IntegerVector SP = above["SP"];
  
  NumericVector Tmax_LAI = speciesNumericParameterWithImputation(SP, SpParams, "Tmax_LAI", true);
  NumericVector Tmax_LAIsq = speciesNumericParameterWithImputation(SP, SpParams, "Tmax_LAIsq", true);
  NumericVector WUE = speciesNumericParameterWithImputation(SP, SpParams, "WUE", fillMissingSpParams);
  NumericVector WUE_par = speciesNumericParameterWithImputation(SP, SpParams, "WUE_par", true);
  NumericVector WUE_co2 = speciesNumericParameterWithImputation(SP, SpParams, "WUE_co2", true);
  NumericVector WUE_vpd = speciesNumericParameterWithImputation(SP, SpParams, "WUE_vpd", true);
  NumericVector Psi_Extract = speciesNumericParameterWithImputation(SP, SpParams, "Psi_Extract", fillMissingSpParams);
  NumericVector Exp_Extract = speciesNumericParameterWithImputation(SP, SpParams, "Exp_Extract", fillMissingSpParams);
  NumericVector VCstem_c = speciesNumericParameterWithImputation(SP, SpParams, "VCstem_c", fillMissingSpParams);
  NumericVector VCstem_d = speciesNumericParameterWithImputation(SP, SpParams, "VCstem_d", fillMissingSpParams);
  NumericVector Gswmin = speciesNumericParameterWithImputation(SP, SpParams, "Gswmin", fillMissingSpParams);
  
  DataFrame paramsTranspirationdf = DataFrame::create(_["Gswmin"] = Gswmin,
                                                      _["Tmax_LAI"] = Tmax_LAI,
                                                      _["Tmax_LAIsq"] = Tmax_LAIsq,
                                                      _["Psi_Extract"]=Psi_Extract,
                                                      _["Exp_Extract"]=Exp_Extract,
                                                      _["VCstem_c"] = VCstem_c,
                                                      _["VCstem_d"] = VCstem_d,
                                                      _["WUE"] = WUE, 
                                                      _["WUE_par"] = WUE_par, 
                                                      _["WUE_co2"] = WUE_co2,
                                                      _["WUE_vpd"] = WUE_vpd);
  paramsTranspirationdf.attr("row.names") = above.attr("row.names");
  return(paramsTranspirationdf);
}
DataFrame paramsTranspirationSperry(DataFrame above, List soil, DataFrame SpParams, 
                              DataFrame paramsAnatomydf, List control) {
  IntegerVector SP = above["SP"];
  NumericVector H = above["H"];
  int numCohorts = SP.size();
  
  double maximumStemConductance = control["maximumStemConductance"];
  double fracRootResistance = control["fracRootResistance"];
  double fracLeafResistance = control["fracLeafResistance"];
  String transpirationMode = control["transpirationMode"];
  
  bool fillMissingSpParams = control["fillMissingSpParams"];
  
  NumericVector dVec = soil["dVec"];
  
  NumericVector Vmax298 = speciesNumericParameterWithImputation(SP, SpParams, "Vmax298", fillMissingSpParams);
  NumericVector Jmax298 = speciesNumericParameterWithImputation(SP, SpParams, "Jmax298", fillMissingSpParams);
  NumericVector VCleaf_kmax = speciesNumericParameterWithImputation(SP, SpParams, "VCleaf_kmax", fillMissingSpParams);
  NumericVector Gswmax = speciesNumericParameterWithImputation(SP, SpParams, "Gswmax", fillMissingSpParams);
  NumericVector Gswmin = speciesNumericParameterWithImputation(SP, SpParams, "Gswmin", fillMissingSpParams);
  NumericVector Kmax_stemxylem = speciesNumericParameterWithImputation(SP, SpParams, "Kmax_stemxylem", fillMissingSpParams);
  NumericVector Kmax_rootxylem = speciesNumericParameterWithImputation(SP, SpParams, "Kmax_rootxylem", fillMissingSpParams);
  NumericVector VCstem_c = speciesNumericParameterWithImputation(SP, SpParams, "VCstem_c", fillMissingSpParams);
  NumericVector VCstem_d = speciesNumericParameterWithImputation(SP, SpParams, "VCstem_d", fillMissingSpParams);
  NumericVector VCleaf_c = speciesNumericParameterWithImputation(SP, SpParams, "VCleaf_c", fillMissingSpParams);
  NumericVector VCleaf_d = speciesNumericParameterWithImputation(SP, SpParams, "VCleaf_d", fillMissingSpParams);
  NumericVector VCroot_c = speciesNumericParameterWithImputation(SP, SpParams, "VCroot_c", fillMissingSpParams);
  NumericVector VCroot_d = speciesNumericParameterWithImputation(SP, SpParams, "VCroot_d", fillMissingSpParams);
  
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector Hmed =  paramsAnatomydf["Hmed"];
  
  NumericVector VCstem_kmax(numCohorts);
  NumericVector VCroottot_kmax(numCohorts, 0.0);
  NumericVector VGrhizotot_kmax(numCohorts, 0.0);
  NumericVector Plant_kmax(numCohorts, 0.0);

  // Scaled conductance parameters parameters
  for(int c=0;c<numCohorts;c++){
    //Stem maximum conductance (in mmol·m-2·s-1·MPa-1)
    VCstem_kmax[c]=maximumStemHydraulicConductance(Kmax_stemxylem[c], Hmed[c], Al2As[c],H[c],control["taper"]); 
    VCstem_kmax[c]=std::min(VCstem_kmax[c], maximumStemConductance);
    
    //Root maximum conductance
    double rstem = (1.0/VCstem_kmax[c]);
    double rleaf = (1.0/VCleaf_kmax[c]);
    double rtot = (rstem+rleaf)/(1.0 - fracRootResistance);
    double VCroot_kmaxc = 1.0/(rtot - rstem - rleaf);
    VCroottot_kmax[c] = VCroot_kmaxc;
    
    //Leaf maximum conductance
    if(!NumericVector::is_na(fracLeafResistance)) {
      double rstem = (1.0/VCstem_kmax[c]);
      double rtot = rstem/(1.0-fracRootResistance - fracLeafResistance);
      VCleaf_kmax[c] = 1.0/(rtot*fracLeafResistance);
    }
    //Plant kmax
    Plant_kmax[c] = 1.0/((1.0/VCleaf_kmax[c])+(1.0/VCstem_kmax[c])+(1.0/VCroottot_kmax[c]));
  }
  
  DataFrame paramsTranspirationdf = DataFrame::create(
    _["Gswmin"]=Gswmin, _["Gswmax"]=Gswmax,_["Vmax298"]=Vmax298,
      _["Jmax298"]=Jmax298, _["Kmax_stemxylem"] = Kmax_stemxylem, _["Kmax_rootxylem"] = Kmax_rootxylem,
        _["VCleaf_kmax"]=VCleaf_kmax,_["VCleaf_c"]=VCleaf_c,_["VCleaf_d"]=VCleaf_d,
        _["VCstem_kmax"]=VCstem_kmax,_["VCstem_c"]=VCstem_c,_["VCstem_d"]=VCstem_d, 
        _["VCroot_kmax"] = VCroottot_kmax ,_["VCroot_c"]=VCroot_c,_["VCroot_d"]=VCroot_d,
        _["VGrhizo_kmax"] = VGrhizotot_kmax,
        _["Plant_kmax"] = Plant_kmax);
  paramsTranspirationdf.attr("row.names") = above.attr("row.names");
  return(paramsTranspirationdf);
}

// [[Rcpp::export(".paramsBelow")]]
List paramsBelow(DataFrame above, NumericVector Z50, NumericVector Z95, List soil, 
                 DataFrame paramsAnatomydf, DataFrame paramsTranspirationdf, List control) {

  NumericVector dVec = soil["dVec"];
  // NumericVector bd = soil["bd"];
  NumericVector rfc = soil["rfc"];
  NumericVector VG_alpha = soil["VG_alpha"];
  NumericVector VG_n = soil["VG_n"];
  int nlayers = dVec.size();

  NumericVector LAI_live = above["LAI_live"];
  NumericVector N = above["N"];
  int numCohorts = N.size();
  
  
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  String transpirationMode = control["transpirationMode"];
  double averageFracRhizosphereResistance = control["averageFracRhizosphereResistance"];
  

  NumericMatrix V = ldrDistribution(Z50, Z95, dVec);
  V.attr("dimnames") = List::create(above.attr("row.names"), layerNames(nlayers));
  // CharacterVector slnames(V.ncol());
  // for(int i=0;i<V.ncol();i++) slnames[i] = i+1;
  
  NumericVector Wsoil = soil["W"];
  NumericMatrix Wpool = NumericMatrix(numCohorts, nlayers);
  Wpool.attr("dimnames") = V.attr("dimnames");
  
  NumericMatrix L = NumericMatrix(numCohorts, nlayers);
  L.attr("dimnames") = V.attr("dimnames");
  
  for(int c=0;c<numCohorts;c++){
    for(int l=0;l<nlayers;l++) Wpool(c,l) = Wsoil[l]; //Init from soil state
  }
  
  NumericVector poolProportions = LAI_live/sum(LAI_live);
  
  List belowLayers;
  DataFrame belowdf;
  
  if(transpirationMode == "Granier") {
    NumericVector Ar2Al = paramsAnatomydf["Ar2Al"];
    NumericVector FineRootDensity = paramsAnatomydf["FineRootDensity"];
    NumericVector SRL = paramsAnatomydf["SRL"];
    NumericVector RLD = paramsAnatomydf["RLD"];
    
    NumericVector CRSV(numCohorts);
    NumericVector FRB(numCohorts, NA_REAL);
    for(int c=0;c<numCohorts;c++){
      L(c,_) = coarseRootLengths(V(c,_), dVec, 0.5); //Arbitrary ratio (to revise some day)
      CRSV[c] = coarseRootSoilVolume(V(c,_), dVec, 0.5);
      //Assume fine root biomass is half leaf structural biomass
      double LA = leafArea(LAI_live[c], N[c]); //m2 leaf area per individual
      double fineRootArea = Ar2Al[c]*LA;//fine root area in m2
      FRB[c] = fineRootArea/(specificRootSurfaceArea(SRL[c], FineRootDensity[c])*1e-4);
    }
    if(rhizosphereOverlap!="total") {
      belowdf = DataFrame::create(_["Z50"] = Z50,
                                  _["Z95"] = Z95,
                                  _["fineRootBiomass"] = FRB,
                                  _["coarseRootSoilVolume"] = CRSV,
                                  _["poolProportions"] = poolProportions);
      List RHOP;
      if(rhizosphereOverlap=="none") RHOP = nonoverlapHorizontalProportions(V);
      else RHOP = horizontalProportions(poolProportions, CRSV, N, V, dVec, rfc);
      belowLayers = List::create(_["V"] = V,
                                 _["L"] = L,
                                 _["Wpool"] = Wpool,
                                 _["RHOP"] = RHOP);
    } else {
      belowdf = DataFrame::create(_["Z50"] = Z50,
                                  _["Z95"] = Z95,
                                  _["fineRootBiomass"] = FRB,
                                  _["coarseRootSoilVolume"] = CRSV);
      belowLayers = List::create(_["V"] = V,
                                 _["L"] = L,
                                 _["Wpool"] = Wpool);
    }
  } else {
    NumericVector Al2As = paramsAnatomydf["Al2As"];
    NumericVector FineRootDensity = paramsAnatomydf["FineRootDensity"];
    NumericVector SRL = paramsAnatomydf["SRL"];
    NumericVector RLD = paramsAnatomydf["RLD"];
    
    NumericVector Kmax_stemxylem = paramsTranspirationdf["Kmax_stemxylem"];
    NumericVector VCroottot_kmax = paramsTranspirationdf["VCroot_kmax"];
    NumericVector VCroot_c = paramsTranspirationdf["VCroot_c"];
    NumericVector VCroot_d = paramsTranspirationdf["VCroot_d"];
    NumericVector VCstem_kmax = paramsTranspirationdf["VCstem_kmax"];
    NumericVector VCstem_c = paramsTranspirationdf["VCstem_c"];
    NumericVector VCstem_d = paramsTranspirationdf["VCstem_d"];
    NumericVector VCleaf_kmax = paramsTranspirationdf["VCleaf_kmax"];
    NumericVector VCleaf_c = paramsTranspirationdf["VCleaf_c"];
    NumericVector VCleaf_d = paramsTranspirationdf["VCleaf_d"];
    NumericVector VGrhizotot_kmax = paramsTranspirationdf["VGrhizo_kmax"];
    
    
    NumericMatrix RhizoPsi =  NumericMatrix(numCohorts, nlayers);
    RhizoPsi.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
    std::fill(RhizoPsi.begin(), RhizoPsi.end(), -0.033);

    
    NumericVector FRB(numCohorts), CRSV(numCohorts),FRAI(numCohorts);
    NumericVector Ksat = soil["Ksat"];
    for(int c=0;c<numCohorts;c++)  {
      //We use Kmax_stemxylem instead of Kmax_rootxylem because of reliability
      CRSV[c] = coarseRootSoilVolumeFromConductance(Kmax_stemxylem[c], VCroottot_kmax[c], Al2As[c],
                                                    V(c,_), dVec, rfc);
    }
    
    
    NumericMatrix VCroot_kmax(numCohorts, nlayers); 
    NumericMatrix VGrhizo_kmax(numCohorts, nlayers);
    VGrhizo_kmax.attr("dimnames") = V.attr("dimnames");
    VCroot_kmax.attr("dimnames") = V.attr("dimnames");
    NumericVector Vc;
    for(int c=0;c<numCohorts;c++){
      Vc = V(c,_);
      L(c,_) = coarseRootLengthsFromVolume(CRSV[c], V(c,_), dVec, rfc); 
      NumericVector xp = rootxylemConductanceProportions(L(c,_), V(c,_));
      for(int l=0;l<nlayers;l++) {
        VCroot_kmax(c,_) = VCroottot_kmax[c]*xp;
        VGrhizo_kmax(c,l) = V(c,l)*findRhizosphereMaximumConductance(averageFracRhizosphereResistance*100.0, VG_n[l], VG_alpha[l],
                     VCroottot_kmax[c], VCroot_c[c], VCroot_d[c],
                     VCstem_kmax[c], VCstem_c[c], VCstem_d[c],
                     VCleaf_kmax[c], VCleaf_c[c], VCleaf_d[c]);
        VGrhizotot_kmax[c] += VGrhizo_kmax(c,l); 
      }
      FRB[c] = fineRootBiomassPerIndividual(Ksat, VGrhizo_kmax(c,_), LAI_live[c], N[c], 
                                            SRL[c], FineRootDensity[c], RLD[c]);
    }
    belowLayers = List::create(_["V"] = V,
                               _["L"] = L,
                               _["VGrhizo_kmax"] = VGrhizo_kmax,
                               _["VCroot_kmax"] = VCroot_kmax,
                               _["Wpool"] = Wpool,
                               _["RhizoPsi"] = RhizoPsi);
    if(rhizosphereOverlap!="total") {
      belowdf = DataFrame::create(_["Z50"]=Z50,
                                  _["Z95"]=Z95,
                                  _["fineRootBiomass"] = FRB,
                                  _["coarseRootSoilVolume"] = CRSV,
                                  _["poolProportions"] = poolProportions);
      List RHOP;
      if(rhizosphereOverlap=="none") RHOP = nonoverlapHorizontalProportions(V);
      else RHOP = horizontalProportions(poolProportions, CRSV, N, V, dVec, rfc);
      belowLayers["RHOP"] = RHOP;
    } else {
      belowdf = DataFrame::create(_["Z50"]=Z50,
                                  _["Z95"]=Z95,
                                  _["fineRootBiomass"] = FRB,
                                  _["coarseRootSoilVolume"] = CRSV);
    }
  } 
  belowdf.attr("row.names") = above.attr("row.names");
  
  List below = List::create(_["below"] = belowdf,
                            _["belowLayers"] = belowLayers);
  return(below);
}

DataFrame paramsGrowth(DataFrame above, DataFrame SpParams, List control) {
  
  bool fillMissingSpParams = control["fillMissingSpParams"];
  
  IntegerVector SP = above["SP"];
  NumericVector DBH = above["DBH"];
  int numCohorts = SP.size();

  NumericVector WoodC = speciesNumericParameterWithImputation(SP, SpParams, "WoodC", fillMissingSpParams);
  
  NumericVector RERleaf = speciesNumericParameterWithImputation(SP, SpParams, "RERleaf", fillMissingSpParams);
  NumericVector RERsapwood = speciesNumericParameterWithImputation(SP, SpParams, "RERsapwood", fillMissingSpParams);
  NumericVector RERfineroot = speciesNumericParameterWithImputation(SP, SpParams, "RERfineroot", fillMissingSpParams);
  NumericVector SRsapwood = speciesNumericParameterWithImputation(SP, SpParams, "SRsapwood", fillMissingSpParams);
  
  
  NumericVector CCleaf = speciesNumericParameterFromIndex(SP, SpParams, "CCleaf");
  NumericVector CCsapwood = speciesNumericParameterFromIndex(SP, SpParams, "CCsapwood");
  NumericVector CCfineroot = speciesNumericParameterFromIndex(SP, SpParams, "CCfineroot");
  NumericVector RGRleafmax = speciesNumericParameterFromIndex(SP, SpParams, "RGRleafmax");
  NumericVector RGRsapwoodmax = speciesNumericParameterFromIndex(SP, SpParams, "RGRsapwoodmax");
  NumericVector RGRcambiummax = speciesNumericParameterFromIndex(SP, SpParams, "RGRcambiummax");
  NumericVector RGRfinerootmax = speciesNumericParameterFromIndex(SP, SpParams, "RGRfinerootmax");
  NumericVector SRfineroot = speciesNumericParameterFromIndex(SP, SpParams, "SRfineroot");
  NumericVector fHDmin = speciesNumericParameterFromIndex(SP, SpParams, "fHDmin");
  NumericVector fHDmax = speciesNumericParameterFromIndex(SP, SpParams, "fHDmax");

  double minimumRelativeStarchForGrowth_default = control["minimumRelativeStarchForGrowth"];
  NumericVector RSSG = speciesNumericParameterFromIndex(SP, SpParams, "RSSG");
  
  NumericVector MortalityBaselineRate = speciesNumericParameterFromIndex(SP, SpParams, "MortalityBaselineRate");
  
  List maximumRelativeGrowthRates = control["maximumRelativeGrowthRates"];
  double RGRleafmax_default = maximumRelativeGrowthRates["leaf"];
  double RGRsapwoodmax_default = maximumRelativeGrowthRates["sapwood"];
  double RGRcambiummax_default = maximumRelativeGrowthRates["cambium"];
  double RGRfinerootmax_default = maximumRelativeGrowthRates["fineroot"];
  
  
  List constructionCosts = control["constructionCosts"];
  double CCleaf_default = constructionCosts["leaf"];
  double CCsapwood_default = constructionCosts["sapwood"];
  double CCfineroot_default = constructionCosts["fineroot"];
  
  List senescenceRates = control["senescenceRates"];
  double SRsapwood_default = senescenceRates["sapwood"];
  double SRfineroot_default = senescenceRates["fineroot"];
  
  double mortalityBaselineRate_default = control["mortalityBaselineRate"];
  
  
  if(fillMissingSpParams) {
    for(int c=0;c<numCohorts;c++){
      if(NumericVector::is_na(CCleaf[c])) CCleaf[c] = CCleaf_default;
      if(NumericVector::is_na(CCsapwood[c])) CCsapwood[c] = CCsapwood_default;
      if(NumericVector::is_na(CCfineroot[c])) CCfineroot[c] = CCfineroot_default;
      if(NumericVector::is_na(RGRleafmax[c])) RGRleafmax[c] = RGRleafmax_default;
      if(NumericVector::is_na(RGRcambiummax[c]) && !NumericVector::is_na(DBH[c])) RGRcambiummax[c] = RGRcambiummax_default;
      if(NumericVector::is_na(RGRsapwoodmax[c]) && NumericVector::is_na(DBH[c])) RGRsapwoodmax[c] = RGRsapwoodmax_default;
      if(NumericVector::is_na(RGRfinerootmax[c])) RGRfinerootmax[c] = RGRfinerootmax_default;
      if(NumericVector::is_na(SRsapwood[c])) SRsapwood[c] = SRsapwood_default;
      if(NumericVector::is_na(SRfineroot[c])) SRfineroot[c] = SRfineroot_default;
      if(NumericVector::is_na(RSSG[c])) RSSG[c] = minimumRelativeStarchForGrowth_default;
      if(NumericVector::is_na(MortalityBaselineRate[c])) MortalityBaselineRate[c] = mortalityBaselineRate_default;
    }
  }
  
  DataFrame paramsGrowthdf = DataFrame::create(_["RERleaf"] = RERleaf,
                                               _["RERsapwood"] = RERsapwood,
                                               _["RERfineroot"] = RERfineroot,
                                               _["CCleaf"] = CCleaf,
                                               _["CCsapwood"] = CCsapwood,
                                               _["CCfineroot"] = CCfineroot,
                                               _["RGRleafmax"] = RGRleafmax,
                                               _["RGRsapwoodmax"] = RGRsapwoodmax,
                                               _["RGRcambiummax"] = RGRcambiummax,
                                               _["RGRfinerootmax"] = RGRfinerootmax,
                                               _["SRsapwood"] = SRsapwood,
                                               _["SRfineroot"] = SRfineroot,
                                               _["RSSG"] = RSSG,
                                               _["fHDmin"] = fHDmin,
                                               _["fHDmax"] = fHDmax,
                                               _["WoodC"] = WoodC,
                                               _["MortalityBaselineRate"] = MortalityBaselineRate);
  paramsGrowthdf.attr("row.names") = above.attr("row.names");
  return(paramsGrowthdf);
}


DataFrame paramsAllometries(DataFrame above, DataFrame SpParams, bool fillMissingSpParams) {
  IntegerVector SP = above["SP"];
  
  NumericVector Afbt = speciesNumericParameterWithImputation(SP, SpParams, "a_fbt",fillMissingSpParams);
  NumericVector Bfbt = speciesNumericParameterWithImputation(SP, SpParams, "b_fbt",fillMissingSpParams);
  NumericVector Cfbt = speciesNumericParameterWithImputation(SP, SpParams, "c_fbt",fillMissingSpParams);
  NumericVector Aash = speciesNumericParameterWithImputation(SP, SpParams, "a_ash",fillMissingSpParams);
  NumericVector Bash = speciesNumericParameterWithImputation(SP, SpParams, "b_ash",fillMissingSpParams);
  NumericVector Absh = speciesNumericParameterWithImputation(SP, SpParams, "a_bsh",fillMissingSpParams);
  NumericVector Bbsh = speciesNumericParameterWithImputation(SP, SpParams, "b_bsh",fillMissingSpParams);
  NumericVector Acr = speciesNumericParameterWithImputation(SP, SpParams, "a_cr",fillMissingSpParams);
  NumericVector B1cr = speciesNumericParameterWithImputation(SP, SpParams, "b_1cr",fillMissingSpParams);
  NumericVector B2cr = speciesNumericParameterWithImputation(SP, SpParams, "b_2cr",fillMissingSpParams);
  NumericVector B3cr = speciesNumericParameterWithImputation(SP, SpParams, "b_3cr",fillMissingSpParams);
  NumericVector C1cr = speciesNumericParameterWithImputation(SP, SpParams, "c_1cr",fillMissingSpParams);
  NumericVector C2cr = speciesNumericParameterWithImputation(SP, SpParams, "c_2cr",fillMissingSpParams);
  NumericVector Acw = speciesNumericParameterWithImputation(SP, SpParams, "a_cw",fillMissingSpParams);
  NumericVector Bcw = speciesNumericParameterWithImputation(SP, SpParams, "b_cw",fillMissingSpParams);

  DataFrame paramsAllometriesdf = DataFrame::create(_["Afbt"] = Afbt, _["Bfbt"] = Bfbt, _["Cfbt"] = Cfbt,
                                                    _["Aash"] = Aash, _["Bash"] = Bash, _["Absh"] = Absh, _["Bbsh"] = Bbsh,
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
                                  DataFrame belowdf,
                                  List belowLayers,
                                  DataFrame paramsAnatomydf,
                                  DataFrame paramsWaterStoragedf,
                                  DataFrame paramsGrowthdf,
                                  List control) {
  int numCohorts = above.nrow();

  double nonSugarConcentration = control["nonSugarConcentration"];
  String transpirationMode = control["transpirationMode"];
  
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector LeafDensity = paramsAnatomydf["LeafDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector conduit2sapwood = paramsAnatomydf["conduit2sapwood"];
  NumericVector LeafPI0 = paramsWaterStoragedf["LeafPI0"];
  NumericVector StemPI0 = paramsWaterStoragedf["StemPI0"];
  
  NumericVector WoodC = paramsGrowthdf["WoodC"];

  NumericMatrix V = belowLayers["V"];
  NumericMatrix L = belowLayers["L"];
  
  NumericVector Z95 = belowdf["Z95"];

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
    double svol = sapwoodStorageVolume(SA[c], H[c], L(c,_), V(c,_),WoodDensity[c], conduit2sapwood[c]);
    
    // 70% in starch storage for sapwood and 1% in leaves
    if(LAI_expanded[c]>0.0) starchLeaf[c] = (0.01/(lvol))*leafStarchCapacity(LAI_expanded[c], N[c], SLA[c], LeafDensity[c]);
    starchSapwood[c] = (0.7/(svol))*sapwoodStarchCapacity(SA[c], H[c], L, V(c,_), WoodDensity[c], conduit2sapwood[c]);
    // starch[c] = starchLeaf[c]+starchSapwood[c];
    
    //Sugar storage from PI0
    if(LAI_expanded[c]>0.0) {
      double lconc = sugarConcentration(LeafPI0[c],20.0, nonSugarConcentration);
      sugarLeaf[c] = lconc;
    }
    double sconc = sugarConcentration(StemPI0[c],20.0, nonSugarConcentration);
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

DataFrame internalMortalityDataFrame(DataFrame above) {
  int numCohorts = above.nrow();
  NumericVector N_dead(numCohorts, 0.0);
  NumericVector Cover_dead(numCohorts, 0.0);
  DataFrame df = DataFrame::create(Named("N_dead") = N_dead,
                                   Named("Cover_dead") = Cover_dead);
  df.attr("row.names") = above.attr("row.names");
  return(df);
}
DataFrame internalAllocationDataFrame(DataFrame above, 
                                      DataFrame belowdf, 
                                      DataFrame paramsAnatomydf,
                                      DataFrame paramsTranspirationdf,
                                      List control) {
  int numCohorts = above.nrow();

  NumericVector allocationTarget(numCohorts,0.0);
  NumericVector leafAreaTarget(numCohorts,0.0);
  NumericVector sapwoodAreaTarget(numCohorts,0.0);
  NumericVector fineRootBiomassTarget(numCohorts, 0.0);
  
  
  String transpirationMode = control["transpirationMode"];
  NumericVector SA = above["SA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  NumericVector fineRootBiomass = belowdf["fineRootBiomass"];
  DataFrame df;
  if(transpirationMode=="Granier") {
    for(int c=0;c<numCohorts;c++){
      leafAreaTarget[c] = Al2As[c]*(SA[c]/10000.0);
      sapwoodAreaTarget[c] = SA[c];
      allocationTarget[c] = Al2As[c];
      fineRootBiomassTarget[c] = fineRootBiomass[c];
    }
  } else {
    String allocationStrategy = control["allocationStrategy"];
    NumericVector Plant_kmax = paramsTranspirationdf["Plant_kmax"];
    NumericVector VGrhizo_kmax = paramsTranspirationdf["VGrhizo_kmax"];
    // NumericVector longtermStorage(numCohorts,0.0);
    for(int c=0;c<numCohorts;c++){
      leafAreaTarget[c] = Al2As[c]*(SA[c]/10000.0);
      if(allocationStrategy=="Plant_kmax") {
        allocationTarget[c] = Plant_kmax[c];
      } else if(allocationStrategy=="Al2As") {
        allocationTarget[c] = Al2As[c];
      }
      sapwoodAreaTarget[c] = SA[c];
      fineRootBiomassTarget[c] = fineRootBiomass[c];
    }
  }
  df = DataFrame::create(Named("allocationTarget") = allocationTarget,
                         Named("leafAreaTarget") = leafAreaTarget,
                         Named("sapwoodAreaTarget") = sapwoodAreaTarget,
                         Named("fineRootBiomassTarget") = fineRootBiomassTarget);
  df.attr("row.names") = above.attr("row.names");
  return(df);
}  


DataFrame internalWaterDataFrame(DataFrame above, String transpirationMode) {
  int numCohorts = above.nrow();
  DataFrame df;
  if(transpirationMode=="Granier") {
    df = DataFrame::create(Named("PlantPsi") = NumericVector(numCohorts, -0.033),
                           Named("StemPLC") = NumericVector(numCohorts, 0.0));
  } else {
    df = DataFrame::create(Named("Einst") = NumericVector(numCohorts, 0.0),
                           Named("RootCrownPsi") = NumericVector(numCohorts, -0.033),
                           Named("Stem1Psi") = NumericVector(numCohorts, -0.033),
                           Named("Stem2Psi") = NumericVector(numCohorts, -0.033),
                           Named("LeafPsi") = NumericVector(numCohorts, -0.033),
                           Named("StemSympPsi") = NumericVector(numCohorts, -0.033),
                           Named("LeafSympPsi") = NumericVector(numCohorts, -0.033),
                           Named("StemPLC") = NumericVector(numCohorts, 0.0));
  }
  df.attr("row.names") = above.attr("row.names");
  return(df);
}

DataFrame paramsCanopy(DataFrame above, List control) {
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector H = above["H"];
  int numCohorts = H.size();
  //Determine number of vertical layers
  double verticalLayerSize = control["verticalLayerSize"];
  double boundaryLayerSize = control["boundaryLayerSize"];
  double canopyHeight = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if((canopyHeight<H[c]) && ((LAI_live[c]+LAI_dead[c])>0.0)) canopyHeight = H[c];
  }
  int nz = ceil((canopyHeight+boundaryLayerSize)/verticalLayerSize); //Number of vertical layers (adding 2 m above to match wind measurement height)
  NumericVector zlow(nz,0.0);
  NumericVector zmid(nz, verticalLayerSize/2.0);
  NumericVector zup(nz, verticalLayerSize);
  for(int i=1;i<nz;i++) {
    zlow[i] = zlow[i-1] + verticalLayerSize;
    zmid[i] = zmid[i-1] + verticalLayerSize;
    zup[i] = zup[i-1] + verticalLayerSize;
  }
  DataFrame paramsCanopy = DataFrame::create(_["zlow"] = zlow,
                                             _["zmid"] = zmid,
                                             _["zup"] = zup,
                                             _["Tair"] = NumericVector(nz, NA_REAL),
                                             _["Cair"] = NumericVector(nz, NA_REAL),
                                             _["VPair"] = NumericVector(nz, NA_REAL));
  return(paramsCanopy);
}

//' Input for simulation models
//'
//' Functions \code{forest2spwbInput} and \code{forest2growthInput} take an object of class \code{\link{forest}} 
//' and calculate input data for functions \code{\link{spwb}}, \code{\link{pwb}} and \code{\link{growth}}, respectively. 
//' Functions \code{spwbInput} and \code{growthInput} do the same but starting from different input data. 
//' Function \code{forest2aboveground} calculates aboveground variables that may be used in \code{spwbInput} and \code{growthInput} functions. 
//' Function \code{forest2belowground} calculates belowground variables such as fine root distribution.
//' 
//' @param x An object of class \code{\link{forest}}.
//' @param SpParams A data frame with species parameters (see \code{\link{SpParamsDefinition}} and \code{\link{SpParamsMED}}).
//' @param gdd Growth degree days to account for leaf phenology effects (in Celsius). This should be left \code{NA} in most applications.
//' @param mode Calculation mode, either "MED" or "US".
//' @param soil An object of class \code{\link{soil}}.
//' @param control A list with default control parameters (see \code{\link{defaultControl}}).
//' @param above A data frame with aboveground plant information (see the return value of \code{forest2aboveground} below). In the case of \code{spwbInput} the variables should include \code{SP}, \code{N}, \code{LAI_live}, \code{LAI_dead}, \code{H} and \code{CR}. In the case of \code{growthInput} variables should include \code{DBH} and \code{Cover}.
//' @param Z50,Z95 Numeric vectors with cohort depths (in mm) corresponding to 50\% and 95\% of fine roots.
//' 
//' @details
//' Functions \code{forest2spwbInput} and \code{forest2abovegroundInput} extracts height and species identity from plant cohorts of \code{x}, 
//' and calculate leaf area index and crown ratio. Function \code{forest2spwbInput} also calculates the distribution of fine roots 
//' across soil. Both \code{forest2spwbInput} and \code{spwbInput} find parameter values for each plant cohort 
//' according to the parameters of its species as specified in \code{SpParams}. If \code{control$transpirationMode = "Sperry"} 
//' the functions also estimate the maximum conductance of rhizosphere, root xylem and stem xylem elements.
//' 
//' @return 
//' Function \code{forest2aboveground()} returns a data frame with the following columns (rows are identified as specified by function \code{\link{plant_ID}}):
//' \itemize{
//'   \item{\code{SP}: Species identity (an integer) (first species is 0).}
//'   \item{\code{N}: Cohort density (ind/ha) (see function \code{\link{plant_density}}).}
//'   \item{\code{DBH}: Tree diameter at breast height (cm).}
//'   \item{\code{H}: Plant total height (cm).}
//'   \item{\code{CR}: Crown ratio (crown length to total height) (between 0 and 1).}
//'   \item{\code{LAI_live}: Live leaf area index (m2/m2) (one-side leaf area relative to plot area), includes leaves in winter dormant buds.}
//'   \item{\code{LAI_expanded}: Leaf area index of expanded leaves (m2/m2) (one-side leaf area relative to plot area).}
//'   \item{\code{LAI_dead}: Dead leaf area index (m2/m2) (one-side leaf area relative to plot area).}
//' }
//' 
//' Functions \code{forest2spwbInput()} and \code{spwbInput()} return a list of class \code{spwbInput} with the following elements (rows of data frames are identified as specified by function \code{\link{plant_ID}}):
//'   \itemize{
//'     \item{\code{control}: List with control parameters (see \code{\link{defaultControl}}).}
//'     \item{\code{canopy}: A list of stand-level state variables.}
//'     \item{\code{cohorts}: A data frame with cohort information, with columns \code{SP} and \code{Name}.}
//'     \item{\code{above}: A data frame with columns  \code{H}, \code{CR} and \code{LAI} (see function \code{forest2aboveground}).}
//'     \item{\code{below}: A data frame with columns \code{Z50}, \code{Z95}.  If \code{control$transpirationMode = "Sperry"} additional columns are \code{fineRootBiomass} and \code{coarseRootSoilVolume}.}
//'     \item{\code{belowLayers}: A list. If \code{control$transpirationMode = "Granier"} it contains elements: 
//'       \itemize{
//'         \item{\code{V}: A matrix with the proportion of fine roots of each cohort (in rows) in each soil layer (in columns).}
//'         \item{\code{L}: A matrix with the length of coarse roots of each cohort (in rows) in each soil layer (in columns).}
//'         \item{\code{Wpool}: A matrix with the soil moisture relative to field capacity around the rhizosphere of each cohort (in rows) in each soil layer (in columns).}
//'       }
//'       If \code{control$transpirationMode = "Sperry"} there are the following additional elements:
//'       \itemize{
//'         \item{\code{VGrhizo_kmax}: A matrix with maximum rhizosphere conductance values of each cohort (in rows) in each soil layer (in columns).}
//'         \item{\code{VGroot_kmax}: A matrix with maximum root xylem conductance values of each cohort (in rows) in each soil layer (in columns).}
//'         \item{\code{RhizoPsi}: A matrix with the water potential around the rhizosphere of each cohort (in rows) in each soil layer (in columns).}
//'       }
//'     }
//'     \item{\code{paramsPhenology}: A data frame with leaf phenology parameters:
//'       \itemize{
//'         \item{\code{PhenologyType}: Leaf phenology type.}
//'         \item{\code{LeafDuration}: Leaf duration (in years).}
//'         \item{\code{Sgdd}: Degree days needed for leaf budburst (for winter decideous species).}
//'         \item{\code{Tbgdd}: Base temperature for the calculation of degree days to leaf budburst.}
//'         \item{\code{Ssen}: Degree days corresponding to leaf senescence.}
//'         \item{\code{Phsen}: Photoperiod corresponding to start counting senescence degree-days.}
//'         \item{\code{Tbsen}: Base temperature for the calculation of degree days to leaf senescence.}
//'       }
//'     }
//'     \item{\code{paramsAnatomy}: A data frame with plant anatomy parameters for each cohort:
//'       \itemize{
//'         \item{\code{Hmax}: Maximum plant height (cm).}
//'         \item{\code{Hmed}: Median plant height (cm).}
//'         \item{\code{Al2As}: Leaf area to sapwood area ratio (in m2·m-2).}
//'         \item{\code{Ar2Al}: Fine root area to leaf area ratio (in m2·m-2).}
//'         \item{\code{SLA}: Specific leaf area (mm2/mg = m2/kg).}
//'         \item{\code{LeafWidth}: Leaf width (in cm).}
//'         \item{\code{LeafDensity}: Density of leaf tissue (dry weight over volume).}
//'         \item{\code{WoodDensity}: Density of wood tissue (dry weight over volume).}
//'         \item{\code{FineRootDensity}: Density of fine root tissue (dry weight over volume).}
//'         \item{\code{SRL}: Specific Root length (cm·g-1).}
//'         \item{\code{RLD}: Root length density (cm·cm-3).}
//'         \item{\code{r635}: Ratio between the weight of leaves plus branches and the weight of leaves alone for branches of 6.35 mm.}
//'       }
//'     }
//'     \item{\code{paramsInterception}: A data frame with rain interception and light extinction parameters for each cohort:
//'       \itemize{
//'         \item{\code{kPAR}: PAR extinction coefficient.}
//'         \item{\code{g}: Canopy water retention capacity per LAI unit (mm/LAI).}
//'       }
//'     If \code{control$transpirationMode = "Sperry"} additional columns are:
//'       \itemize{
//'         \item{\code{gammaSWR}: Reflectance (albedo) coefficient for SWR .}
//'         \item{\code{alphaSWR}: Absorbance coefficient for SWR .}
//'       }
//'     }
//'     \item{\code{paramsTranspiration}: A data frame with parameters for transpiration and photosynthesis. If \code{control$transpirationMode = "Granier"}, columns are:
//'       \itemize{
//'         \item{\code{Gswmin}: Minimum stomatal conductance to water vapor (in mol H2O·m-2·s-1).}
//'         \item{\code{Tmax_LAI}: Coefficient relating LAI with the ratio of maximum transpiration over potential evapotranspiration.}
//'         \item{\code{Tmax_LAIsq}: Coefficient relating squared LAI with the ratio of maximum transpiration over potential evapotranspiration.}
//'         \item{\code{Psi_Extract}: Water potential corresponding to 50\% relative transpiration (in MPa).}
//'         \item{\code{Exp_Extract}: Parameter of the Weibull function regulating transpiration reduction.}
//'         \item{\code{VCstem_c}, \code{VCstem_d}: Parameters of the stem xylem vulnerability curve.}
//'         \item{\code{WUE}: Daily water use efficiency (gross photosynthesis over transpiration) under no light, water or CO2 limitations and VPD = 1kPa (g C/mm water).}
//'         \item{\code{WUE_par}: Coefficient regulating the influence of \% PAR on gross photosynthesis.}
//'         \item{\code{WUE_par}: Coefficient regulating the influence of atmospheric CO2 concentration on gross photosynthesis.}
//'         \item{\code{WUE_par}: Coefficient regulating the influence of vapor pressure deficit (VPD) on gross photosynthesis.}
//'       }
//'     If \code{control$transpirationMode = "Sperry"} columns are:
//'       \itemize{
//'         \item{\code{Gswmin}: Minimum stomatal conductance to water vapor (in mol H2O·m-2·s-1).}
//'         \item{\code{Gswmax}: Maximum stomatal conductance to water vapor (in mol H2O·m-2·s-1).}
//'         \item{\code{Vmax298}: Maximum Rubisco carboxilation rate at 25ºC (in micromol CO2·s-1·m-2).}
//'         \item{\code{Jmax298}: Maximum rate of electron transport at 25ºC (in micromol photons·s-1·m-2).}
//'         \item{\code{Kmax_stemxylem}: Sapwood-specific hydraulic conductivity of stem xylem (in kg H2O·s-1·m-2).}
//'         \item{\code{Kmax_rootxylem}: Sapwood-specific hydraulic conductivity of root xylem (in kg H2O·s-1·m-2).}
//'         \item{\code{VCleaf_kmax}: Maximum leaf hydraulic conductance.}
//'         \item{\code{VCleaf_c}, \code{VCleaf_d}: Parameters of the leaf vulnerability curve.}
//'         \item{\code{VCstem_kmax}: Maximum stem xylem conductance.}
//'         \item{\code{VCstem_c}, \code{VCstem_d}: Parameters of the stem xylem vulnerability curve.}
//'         \item{\code{VCroot_c}, \code{VCroot_d}: Parameters of the root xylem vulnerability curve.}
//'         \item{\code{Plant_kmax}: Maximum whole-plant conductance.}
//'       }
//'     }
//'     \item{\code{paramsWaterStorage}: A data frame with plant water storage parameters for each cohort:
//'       \itemize{
//'         \item{\code{LeafPI0}: Osmotic potential at full turgor of leaves (MPa).}
//'         \item{\code{LeafEPS}: Modulus of elasticity (capacity of the cell wall to resist changes in volume in response to changes in turgor) of leaves (MPa).}
//'         \item{\code{LeafAF}: Apoplastic fraction (proportion of water outside the living cells) in leaves.}
//'         \item{\code{Vleaf}: Storage water capacity in leaves, per leaf area (L/m2).}
//'         \item{\code{StemPI0}: Osmotic potential at full turgor of symplastic xylem tissue (MPa).}
//'         \item{\code{StemEPS}: Modulus of elasticity (capacity of the cell wall to resist changes in volume in response to changes in turgor) of symplastic xylem tissue (Mpa).}
//'         \item{\code{StemAF}: Apoplastic fraction (proportion of water outside the living cells) in stem xylem.}
//'         \item{\code{Vstem}: Storage water capacity in sapwood, per leaf area (L/m2).}
//'       }
//'     }
//'     \item{\code{internalPhenology} and \code{internalWater}: data frames to store internal state variables.}
//'   }
//' Functions \code{forest2growthInput} and \code{growthInput} return a list of class \code{growthInput} with the same elements as \code{spwbInput}, but with additional information. 
//' \itemize{
//' \item{Element \code{above} includes the following additional columns:
//'     \itemize{
//'       \item{\code{LA_live}: Live leaf area per individual (m2/ind).}
//'       \item{\code{LA_dead}: Dead leaf area per individual (m2/ind).}
//'       \item{\code{SA}: Live sapwood area per individual (cm2/ind).} 
//'   }
//'   }
//'   \item{\code{paramsGrowth}: A data frame with growth parameters for each cohort:
//'     \itemize{
//'       \item{\code{RERleaf}: Maintenance respiration rates (at 20ºC) for leaves (in g gluc·g dry-1·day-1).}
//'       \item{\code{RERsapwood}: Maintenance respiration rates (at 20ºC) for sapwood (in g gluc·g dry-1·day-1).}
//'       \item{\code{RERfineroot}: Maintenance respiration rates (at 20ºC) for fine roots (in g gluc·g dry-1·day-1).}
//'       \item{\code{CCleaf}: Leaf construction costs (in g gluc·g dry-1).}
//'       \item{\code{CCsapwood}: Sapwood construction costs (in g gluc·g dry-1).}
//'       \item{\code{CCfineroot}: Fine root construction costs (in g gluc·g dry-1).}
//'       \item{\code{RGRleafmax}: Maximum leaf relative growth rate (in m2·cm-2·day-1).}
//'       \item{\code{RGRsapwoodmax}: Maximum sapwood relative growth rate (in cm2·cm-2·day-1).}
//'       \item{\code{RGRfinerootmax}: Maximum fine root relative growth rate (in g dry·g dry-1·day-1).}
//'       \item{\code{SRsapwood}: Sapwood daily senescence rate (in day-1).}
//'       \item{\code{SRfineroot}: Fine root daily senescence rate (in day-1).}
//'       \item{\code{RSSG}: Minimum relative starch for sapwood growth (proportion).}
//'       \item{\code{fHDmin}: Minimum value of the height-to-diameter ratio (dimensionless).}
//'       \item{\code{fHDmax}: Maximum value of the height-to-diameter ratio (dimensionless).}
//'       \item{\code{WoodC}: Wood carbon content per dry weight (g C /g dry).}
//'       \item{\code{MortalityBaselineRate}: Deterministic proportion or probability specifying the baseline reduction of cohort's density occurring in a year.}
//'     }
//'   }
//'   \item{\code{paramsAllometry}: A data frame with allometric parameters for each cohort:
//'     \itemize{
//'       \item{\code{Aash}: Regression coefficient relating the square of shrub height with shrub area.}
//'       \item{\code{Absh}, \code{Bbsh}: Allometric coefficients relating phytovolume with dry weight of shrub individuals.}
//'       \item{\code{Acr}, \code{B1cr}, \code{B2cr}, \code{B3cr}, \code{C1cr}, \code{C2cr}: Regression coefficients used to calculate crown ratio of trees.}
//'       \item{\code{Acw}, \code{Bcw}: Regression coefficients used to calculated crown width of trees.}
//'     }
//'   }
//'   \item {\code{internalAllocation}: A data frame with internal allocation variables for each cohort:
//'     \itemize{
//'       \item{\code{allocationTarget}: Value of the allocation target variable.}
//'       \item{\code{leafAreaTarget}: Target leaf area (m2) per individual.}
//'       \item{\code{fineRootBiomassTarget}: Target fine root biomass (g dry) per individual (only if \code{transpirationMode = "Sperry"}).}
//'     }
//'   }
//'   \item{\code{internalCarbon} and \code{internalRings}: data structures to store other internal state variables.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{resetInputs}}, \code{\link{spwb}}, \code{\link{soil}},  
//' \code{\link{forest}}, \code{\link{SpParamsMED}}, \code{\link{defaultSoilParams}}, \code{\link{plant_ID}}
//' 
//' @examples
//' #Load example plot plant data
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' # Aboveground parameters
//' above = forest2aboveground(exampleforestMED, SpParamsMED)
//' above
//' 
//' # Initialize soil with default soil params
//' examplesoil = soil(defaultSoilParams())
//' 
//' # Rooting depths
//' Z50 = c(exampleforestMED$treeData$Z50, exampleforestMED$shrubData$Z50)
//' Z95 = c(exampleforestMED$treeData$Z95, exampleforestMED$shrubData$Z95)
//' 
//' # Initialize control parameters
//' control = defaultControl("Granier")
//' 
//' # Prepare spwb input
//' spwbInput(above, Z50, Z95, examplesoil,SpParamsMED, control)
//' 
//' # When starting from an object of class 'forest' the whole process
//' # can be simplified:
//' forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)
//'                 
//' # Prepare input for Sperry transpiration mode
//' control = defaultControl("Sperry")
//' forest2spwbInput(exampleforestMED,examplesoil,SpParamsMED, control)
//' 
//' @name modelInput
// [[Rcpp::export("spwbInput")]]
List spwbInput(DataFrame above, NumericVector Z50, NumericVector Z95, List soil, DataFrame SpParams, List control) {
  
  
  IntegerVector SP = above["SP"];
  NumericVector LAI_live = above["LAI_live"];
  NumericVector LAI_expanded = above["LAI_expanded"];
  NumericVector LAI_dead = above["LAI_dead"];
  NumericVector H = above["H"];
  NumericVector DBH = above["DBH"];
  NumericVector N = above["N"];
  NumericVector CR = above["CR"];
  
  String transpirationMode = control["transpirationMode"];
  if((transpirationMode!="Granier") && (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Granier' or 'Sperry')");

  bool fillMissingSpParams = control["fillMissingSpParams"];
  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") && (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  
  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  

  //Cohort description
  CharacterVector nsp = speciesCharacterParameterFromIndex(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  //Above 
  DataFrame plantsdf = DataFrame::create(_["H"]=H, _["CR"]=CR, _["N"] = N,
                                         _["LAI_live"]=LAI_live, 
                                         _["LAI_expanded"] = LAI_expanded, 
                                         _["LAI_dead"] = LAI_dead);
  plantsdf.attr("row.names") = above.attr("row.names");
  
  DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams, fillMissingSpParams, "spwb", transpirationMode);
  
  DataFrame paramsTranspirationdf;
  if(transpirationMode=="Granier") {
    paramsTranspirationdf = paramsTranspirationGranier(above,SpParams, fillMissingSpParams);
  } else {
    paramsTranspirationdf = paramsTranspirationSperry(above, soil, SpParams, paramsAnatomydf, control);
  }

  List below = paramsBelow(above, Z50, Z95, soil, 
                           paramsAnatomydf, paramsTranspirationdf, control);
  List belowLayers = below["belowLayers"];
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(below["below"]);
  
  DataFrame paramsWaterStoragedf = paramsWaterStorage(above, belowLayers, SpParams, paramsAnatomydf, fillMissingSpParams);
  
  DataFrame paramsCanopydf;
  List ctl = clone(control);
  if(transpirationMode=="Granier") {
    paramsCanopydf = List::create();
  } else if(transpirationMode =="Sperry"){
    paramsCanopydf = paramsCanopy(above, control);
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      ctl["soilFunctions"] = soilFunctions;
      Rcerr<<"Soil pedotransfer functions set to Van Genuchten ('VG').\n";
    }
  }
  List input = List::create(_["control"] = ctl,
                            _["soil"] = clone(soil),
                            _["canopy"] = paramsCanopydf,
                            _["cohorts"] = cohortDescdf,
                            _["above"] = plantsdf,
                            _["below"] = belowdf,
                            _["belowLayers"] = belowLayers,
                            _["paramsPhenology"] = paramsPhenology(above, SpParams, fillMissingSpParams),
                            _["paramsAnatomy"] = paramsAnatomydf,
                            _["paramsInterception"] = paramsInterception(above, SpParams, control),
                            _["paramsTranspiration"] = paramsTranspirationdf,
                            _["paramsWaterStorage"] = paramsWaterStoragedf,
                            _["internalPhenology"] = internalPhenologyDataFrame(above),
                            _["internalWater"] = internalWaterDataFrame(above, transpirationMode));
  
  input.attr("class") = CharacterVector::create("spwbInput","list");
  return(input);
}



//' @rdname modelInput
// [[Rcpp::export("growthInput")]]
List growthInput(DataFrame above, NumericVector Z50, NumericVector Z95, List soil, DataFrame SpParams, List control) {

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
  if((transpirationMode!="Granier") && (transpirationMode!="Sperry")) stop("Wrong Transpiration mode ('transpirationMode' should be either 'Granier' or 'Sperry')");

  bool fillMissingSpParams = control["fillMissingSpParams"];
  
  String soilFunctions = control["soilFunctions"]; 
  if((soilFunctions!="SX") && (soilFunctions!="VG")) stop("Wrong soil functions ('soilFunctions' should be either 'SX' or 'VG')");
  
  
  DataFrame paramsAnatomydf = paramsAnatomy(above, SpParams, fillMissingSpParams, "growth",transpirationMode);
  NumericVector WoodDensity = paramsAnatomydf["WoodDensity"];
  NumericVector SLA = paramsAnatomydf["SLA"];
  NumericVector Al2As = paramsAnatomydf["Al2As"];
  
  DataFrame paramsGrowthdf = paramsGrowth(above, SpParams, control);
  
  DataFrame paramsAllometriesdf = paramsAllometries(above, SpParams, fillMissingSpParams);
  
  
  DataFrame paramsTranspirationdf;
  if(transpirationMode=="Granier") {
    paramsTranspirationdf = paramsTranspirationGranier(above,SpParams, fillMissingSpParams);
  } else {
    paramsTranspirationdf = paramsTranspirationSperry(above, soil, SpParams, paramsAnatomydf, control);
  }


  int numCohorts = SP.size();
  for(int c=0;c<numCohorts;c++){
    if(NumericVector::is_na(CR[c])) CR[c] = 0.5; //PATCH TO AVOID MISSING VALUES!!!!
  }
  
  
  NumericVector SA(numCohorts);
  for(int c=0;c<numCohorts;c++){
    SA[c] = 10000.0*(LAI_live[c]/(N[c]/10000.0))/Al2As[c];//Individual SA in cm2
  }
  

  
  //Cohort description
  CharacterVector nsp = speciesCharacterParameterFromIndex(SP, SpParams, "Name");
  DataFrame cohortDescdf = DataFrame::create(_["SP"] = SP, _["Name"] = nsp);
  cohortDescdf.attr("row.names") = above.attr("row.names");
  
  DataFrame plantsdf = DataFrame::create(_["SP"]=SP, 
                                         _["N"]=N,
                                         _["DBH"]=DBH, 
                                         _["Cover"] = Cover, 
                                         _["H"]=H, 
                                         _["CR"]=CR,
                                         _["SA"] = SA, 
                                         _["LAI_live"]=LAI_live, 
                                         _["LAI_expanded"]=LAI_expanded, 
                                         _["LAI_dead"] = LAI_dead);
  plantsdf.attr("row.names") = above.attr("row.names");
  

  // List ringList(numCohorts);
  // for(int i=0;i<numCohorts;i++) ringList[i] = initialize_ring();
  // ringList.attr("names") = above.attr("row.names");
  
  List below = paramsBelow(above, Z50, Z95, soil, 
                           paramsAnatomydf, paramsTranspirationdf, control);
  List belowLayers = below["belowLayers"];
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(below["below"]);
  
  DataFrame paramsWaterStoragedf = paramsWaterStorage(above, belowLayers, SpParams, paramsAnatomydf, fillMissingSpParams);
  
  
  DataFrame paramsCanopydf;
  List ctl = clone(control);
  if(transpirationMode=="Granier") {
    paramsCanopydf = List::create();
  } else if(transpirationMode =="Sperry"){
    paramsCanopydf = paramsCanopy(above, control);
    if(soilFunctions=="SX") {
      soilFunctions = "VG"; 
      ctl["soilFunctions"] = soilFunctions;
      Rcerr<<"Soil pedotransfer functions set to Van Genuchten ('VG').\n";
    }
  } 
  List input = List::create(_["control"] = ctl,
                       _["soil"] = clone(soil),
                       _["canopy"] = paramsCanopydf,
                       _["cohorts"] = cohortDescdf,
                       _["above"] = plantsdf,
                       _["below"] = belowdf,
                       _["belowLayers"] = belowLayers,
                       _["paramsPhenology"] = paramsPhenology(above, SpParams, fillMissingSpParams),
                       _["paramsAnatomy"] = paramsAnatomydf,
                       _["paramsInterception"] = paramsInterception(above, SpParams, control),
                       _["paramsTranspiration"] = paramsTranspirationdf,
                       _["paramsWaterStorage"] = paramsWaterStoragedf,
                       _["paramsGrowth"]= paramsGrowthdf,
                       _["paramsAllometries"] = paramsAllometriesdf,
                       _["internalPhenology"] = internalPhenologyDataFrame(above),
                       _["internalWater"] = internalWaterDataFrame(above, transpirationMode),
                       _["internalCarbon"] = internalCarbonDataFrame(plantsdf, belowdf, belowLayers,
                                                       paramsAnatomydf, 
                                                       paramsWaterStoragedf,
                                                       paramsGrowthdf, control),
                       _["internalAllocation"] = internalAllocationDataFrame(plantsdf, belowdf,
                                                      paramsAnatomydf,
                                                      paramsTranspirationdf, control),
                       _["internalMortality"] = internalMortalityDataFrame(plantsdf));
  
  input.attr("class") = CharacterVector::create("growthInput","list");
  return(input);
}

// [[Rcpp::export(".cloneInput")]]
List cloneInput(List input) {
  return(clone(input));
}

List rootDistributionComplete(List x, DataFrame SpParams, bool fillMissingRootParams){
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  NumericVector Z95(ntree+nshrub), Z50(ntree+nshrub);

  NumericVector treeZ95 = treeData["Z95"];
  NumericVector treeZ50 = treeData["Z50"];

  IntegerVector treeSP, shrubSP;
  if((TYPEOF(treeData["Species"]) == INTSXP) || (TYPEOF(treeData["Species"]) == REALSXP)) {
    treeSP = Rcpp::as<Rcpp::IntegerVector>(treeData["Species"]);
  } else {
    CharacterVector tspecies = Rcpp::as<Rcpp::CharacterVector>(treeData["Species"]);
    treeSP = speciesIndex(tspecies, SpParams);
  }
  if((TYPEOF(shrubData["Species"]) == INTSXP) || (TYPEOF(shrubData["Species"]) == REALSXP)) {
    shrubSP = Rcpp::as<Rcpp::IntegerVector>(shrubData["Species"]);  
  } else {
    CharacterVector sspecies = Rcpp::as<Rcpp::CharacterVector>(shrubData["Species"]);
    shrubSP = speciesIndex(sspecies, SpParams);
  }
  
  NumericVector treeSPZ50 = speciesNumericParameterFromIndex(treeSP, SpParams, "Z50");
  NumericVector treeSPZ95 = speciesNumericParameterFromIndex(treeSP, SpParams, "Z95");
  for(int i=0;i<ntree;i++) {
    Z50[i] = treeZ50[i];
    Z95[i] = treeZ95[i];
    if(fillMissingRootParams) {
      if(NumericVector::is_na(Z50[i])) Z50[i] = treeSPZ50[i];
      if(NumericVector::is_na(Z95[i])) Z95[i] = treeSPZ95[i];
      if(NumericVector::is_na(Z50[i]) && !NumericVector::is_na(Z95[i])) Z50[i] = exp(log(Z95[i])/1.4);
    }
  }
  NumericVector shrubZ95 = shrubData["Z95"];  
  NumericVector shrubZ50 = shrubData["Z50"];  
  NumericVector shrubSPZ50 = speciesNumericParameterFromIndex(shrubSP, SpParams, "Z50");
  NumericVector shrubSPZ95 = speciesNumericParameterFromIndex(shrubSP, SpParams, "Z95");
  for(int i=0;i<nshrub;i++) {
    Z50[ntree+i] = shrubZ50[i]; 
    Z95[ntree+i] = shrubZ95[i]; 
    if(fillMissingRootParams) {
      if(NumericVector::is_na(Z50[ntree+i])) Z50[ntree+i] = shrubSPZ50[i];
      if(NumericVector::is_na(Z95[ntree+i])) Z95[ntree+i] = shrubSPZ95[i];
      if(NumericVector::is_na(Z50[ntree+i]) && !NumericVector::is_na(Z95[ntree+i])) Z50[ntree+i] = exp(log(Z95[ntree+i])/1.4);
    }
  }

  return(List::create(_["Z50"] = Z50, _["Z95"] = Z95));  
}

//' @rdname modelInput
// [[Rcpp::export("forest2spwbInput")]]
List forest2spwbInput(List x, List soil, DataFrame SpParams, List control, String mode = "MED") {
  List rdc = rootDistributionComplete(x, SpParams, control["fillMissingRootParams"]);
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL, mode);
  return(spwbInput(above, rdc["Z50"], rdc["Z95"], soil, SpParams, control));
}


//' @rdname modelInput
// [[Rcpp::export("forest2growthInput")]]
List forest2growthInput(List x, List soil, DataFrame SpParams, List control) {
  List rdc = rootDistributionComplete(x, SpParams, control["fillMissingRootParams"]);
  DataFrame above = forest2aboveground(x, SpParams, NA_REAL);
  return(growthInput(above,  rdc["Z50"], rdc["Z95"], soil, SpParams, control));
}

//' Reset simulation inputs
//' 
//' Function \code{resetInputs()} allows resetting state variables in \code{x} to their defaults.
//' 
//' @param x An object of class \code{\link{spwbInput}} or \code{\link{growthInput}}.
//' 
//' @return Does not return any value. Instead, it modifies input object \code{x}.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso \code{\link{spwbInput}}, \code{\link{growthInput}}, \code{\link{spwb}}
//' 
// [[Rcpp::export("resetInputs")]]
void resetInputs(List x) {
  List control = x["control"];
  List soil = x["soil"];
  String transpirationMode = control["transpirationMode"];
  //Reset of canopy layer state variables 
  if(transpirationMode=="Sperry") {
    DataFrame can = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
    NumericVector Tair = can["Tair"];
    NumericVector Cair = can["Cair"];
    NumericVector VPair = can["VPair"];
    int ncanlayers = can.nrow();
    for(int i=0;i<ncanlayers;i++) {
      Tair[i] = NA_REAL;
      Cair[i] = NA_REAL;
      VPair[i] = NA_REAL;
    }
  }
  //Reset of soil state variables
  NumericVector Wsoil = soil["W"];
  NumericVector Temp = soil["Temp"];
  List belowLayers = x["belowLayers"];
  NumericMatrix Wpool = belowLayers["Wpool"];
  int nlayers = Wsoil.size();
  for(int i=0;i<nlayers;i++) {
    Wsoil[i] = 1.0; //Defaults to soil at field capacity
    Temp[i] = NA_REAL;
  }
  //Reset of cohort-based state variables
  int numCohorts = Wpool.nrow();
  for(int c=0;c<numCohorts;c++) {
    for(int l=0;l<nlayers;l++) {
      Wpool(c,l) = 1.0;
    }
  }
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  
  if(transpirationMode=="Sperry") {
    NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
    NumericVector RootCrownPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
    NumericVector Stem1Psi = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
    NumericVector Stem2Psi = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem2Psi"]);
    NumericVector StemSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
    NumericVector LeafSympPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
    NumericVector LeafPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
    NumericVector Einst = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);
    for(int i=0;i<LeafPsi.size();i++) {
      Einst[i] = 0.0;
      RootCrownPsi[i] = -0.033;
      Stem1Psi[i] = -0.033;
      Stem2Psi[i] = -0.033;
      LeafPsi[i] = -0.033;
      LeafSympPsi[i] = -0.033;
      StemSympPsi[i] = -0.033;
      StemPLC[i] = 0.0;
      for(int j=0;j<RhizoPsi.ncol();j++) RhizoPsi(i,j) = -0.033;
    }
  } else {
    NumericVector PlantPsi = Rcpp::as<Rcpp::NumericVector>(internalWater["PlantPsi"]);
    for(int i=0;i<StemPLC.length();i++) {
      PlantPsi[i] = -0.033;
      StemPLC[i] = 0.0;
    }
  }
}

void updatePlantKmax(List x) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  
  if(transpirationMode=="Sperry") {
    DataFrame paramsTranspirationdf =  Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    NumericVector Plant_kmax = paramsTranspirationdf["Plant_kmax"];
    NumericVector VCleaf_kmax = paramsTranspirationdf["VCleaf_kmax"];
    NumericVector VCstem_kmax = paramsTranspirationdf["VCstem_kmax"];
    NumericVector VCroot_kmax = paramsTranspirationdf["VCroot_kmax"];
    int numCohorts = Plant_kmax.size();
    for(int i=0;i<numCohorts;i++) {
      Plant_kmax[i] = 1.0/((1.0/VCleaf_kmax[i])+(1.0/VCstem_kmax[i])+(1.0/VCroot_kmax[i]));
    }
  }
}
void updateBelowgroundConductances(List x) {
  List soil = x["soil"];
  NumericVector dVec = soil["dVec"];
  NumericVector rfc = soil["dVec"];
  List belowLayers = x["belowLayers"];
  NumericMatrix V = belowLayers["V"];
  NumericMatrix L = belowLayers["L"];
  int numCohorts = V.nrow();
  int nlayers = V.ncol();
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  if(transpirationMode=="Sperry") {
    DataFrame paramsTranspirationdf = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    NumericVector VCroot_kmax = belowLayers["VCroot_kmax"];
    NumericVector VGrhizo_kmax = belowLayers["VGrhizo_kmax"];
    NumericVector VCroottot_kmax = paramsTranspirationdf["VCroot_kmax"];
    NumericVector VGrhizotot_kmax = paramsTranspirationdf["VGrhizo_kmax"];
    for(int c=0;c<numCohorts;c++) {
      NumericVector xp = rootxylemConductanceProportions(L(c,_), V(c,_));
      for(int l=0;l<nlayers;l++)  {
        VCroot_kmax(c,l) = VCroottot_kmax[c]*xp[l]; 
        VGrhizo_kmax(c,l) = VGrhizo_kmax[c]*V(c,l);
      }
    }
  }
}
void updateFineRootDistribution(List x) {
  List soil = x["soil"];
  NumericVector dVec = soil["dVec"];
  DataFrame belowdf =  Rcpp::as<Rcpp::DataFrame>(x["below"]);
  NumericVector Z50 = belowdf["Z50"];
  NumericVector Z95 = belowdf["Z95"];
  List belowLayers = x["belowLayers"];
  NumericMatrix V = belowLayers["V"];
  int numCohorts = V.nrow();
  int nlayers = V.ncol();
  for(int c=0;c<numCohorts;c++) {
    NumericVector PC = ldrRS_one(Z50[c], Z95[c], dVec);
    for(int l=0;l<nlayers;l++) V(c,l) = PC[l]; 
  }
  updateBelowgroundConductances(x);
}
// [[Rcpp::export(".updateBelow")]]
void updateBelow(List x) {
  List control = x["control"];
  List soil = x["soil"];

  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  
  DataFrame paramsAnatomydf = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  DataFrame paramsTranspirationdf = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Z50 = belowdf["Z50"];
  NumericVector Z95 = belowdf["Z95"];
  List newBelowList = paramsBelow(above, Z50, Z95, soil, 
                               paramsAnatomydf, paramsTranspirationdf, control);
  x["below"] = newBelowList["below"];
  x["belowLayers"] = newBelowList["belowLayers"];
}

double getInputParamValue(List x, String paramType, String paramName, int cohort) {
  DataFrame paramdf = Rcpp::as<Rcpp::DataFrame>(x[paramType]);
  NumericVector param = paramdf[paramName];
  return(param[cohort]);
}

void modifyMessage(String paramName, String cohortName, double newValue) {
  Rcerr<< "[Message] Modifying parameter "<< paramName.get_cstring()<< "' of cohort '" << cohortName.get_cstring() << "' to value " << newValue <<".\n";
}
void multiplyMessage(String paramName, String cohortName, double factor) {
  Rcerr<< "[Message] Multiplying parameter "<< paramName.get_cstring()<< "' of cohort '" << cohortName.get_cstring() << "' by factor " << factor <<".\n";
}

void modifyInputParamSingle(List x, String paramType, String paramName, int cohort, double newValue) {
  DataFrame paramdf = Rcpp::as<Rcpp::DataFrame>(x[paramType]);
  NumericVector param = paramdf[paramName];
  param[cohort] = newValue;
}
void multiplyInputParamSingle(List x, String paramType, String paramName, int cohort, double f) {
  DataFrame paramdf = Rcpp::as<Rcpp::DataFrame>(x[paramType]);
  NumericVector param = paramdf[paramName];
  param[cohort] = param[cohort]*f;
}

// [[Rcpp::export(".multiplyInputParam")]]
void multiplyInputParam(List x, String paramType, String paramName, 
                        int cohort, double f,
                        bool message) {
  List control = x["control"];
  String transpirationMode = control["transpirationMode"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  CharacterVector cohNames = cohorts.attr("row.names");
  if(paramName=="Z50/Z95") {
    multiplyInputParamSingle(x, "below", "Z50", cohort, f);
    multiplyInputParamSingle(x, "below", "Z95", cohort, f);
    if(message) Rcerr<< "[Message] Updating fine root distribution for cohort " << cohNames[cohort] <<".\n";
    updateFineRootDistribution(x);
  } else  if(paramName=="WaterStorage") {
    multiplyInputParamSingle(x, "paramsWaterStorage", "Vsapwood", cohort, f);
    multiplyInputParamSingle(x, "paramsWaterStorage", "Vleaf", cohort, f);
  } else if(paramName=="Plant_kmax") {
    multiplyInputParamSingle(x, "paramsTranspiration", "Plant_kmax", cohort, f);
    if(message) multiplyMessage("VCleaf_kmax", cohNames[cohort], f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCleaf_kmax", cohort, f);
    if(message) multiplyMessage("VCstem_kmax", cohNames[cohort], f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCstem_kmax", cohort, f);
    if(message) multiplyMessage("VCroot_kmax", cohNames[cohort], f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCroot_kmax", cohort, f);
    if(message) Rcerr<< "[Message] Updating below-ground conductances for cohort " << cohNames[cohort] <<".\n";
    updateBelowgroundConductances(x);
  } else if(paramName=="LAI_live") {
    multiplyInputParamSingle(x, "above", "LAI_live", cohort, f);
    if(message) multiplyMessage("LAI_expanded", cohNames[cohort], f);
    multiplyInputParamSingle(x, "above", "LAI_expanded", cohort, f);
  } else if(paramName=="c") {
    multiplyInputParamSingle(x, "paramsTranspiration", "VCleaf_c", cohort, f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCstem_c", cohort, f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCroot_c", cohort, f);
  } else if(paramName=="d") {
    multiplyInputParamSingle(x, "paramsTranspiration", "VCleaf_d", cohort, f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCstem_d", cohort, f);
    multiplyInputParamSingle(x, "paramsTranspiration", "VCroot_d", cohort, f);
  }else if(paramName=="Al2As") {
    multiplyInputParamSingle(x, "paramsAnatomy", "Al2As", cohort, f);
    if(message) multiplyMessage("Vsapwood", cohNames[cohort], f);
    multiplyInputParamSingle(x, "paramsWaterStorage", "Vsapwood", cohort, 1.0/f);
    if(transpirationMode=="Sperry") {
      if(message) multiplyMessage("VCstem_kmax", cohNames[cohort], f);
      multiplyInputParamSingle(x, "paramsTranspiration", "VCstem_kmax", cohort, 1.0/f);
      if(message) multiplyMessage("VCroot_kmax", cohNames[cohort], f);
      multiplyInputParamSingle(x, "paramsTranspiration", "VCroot_kmax", cohort, 1.0/f);
    }
    if(x.containsElementNamed("internalAllocation")) {
      DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
      DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
      DataFrame paramsTranspirationdf = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
      DataFrame paramsAnatomydf = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
      if(message) Rcerr<< "[Message] Rebuilding allocation targets for cohort " << cohNames[cohort] <<".\n";
      x["internalAllocation"]  = internalAllocationDataFrame(above, 
                                           belowdf, 
                                           paramsAnatomydf,
                                           paramsTranspirationdf,
                                           control);
    }
  } else if(paramName=="Vmax298/Jmax298") {
    multiplyInputParamSingle(x, "paramsTranspiration", "Vmax298", cohort, f);
    multiplyInputParamSingle(x, "paramsTranspiration", "Jmax298", cohort, f);
  } else {
    multiplyInputParamSingle(x, paramType, paramName, cohort, f);
  }
  if(transpirationMode=="Sperry") {
    if(message) Rcerr<< "[Message] Recalculating plant maximum conductances.\n";
    updatePlantKmax(x);
  }
  if(message) Rcerr<< "[Message] Updating below-ground parameters.\n";
  updateBelow(x);
}


// [[Rcpp::export(".modifyInputParam")]]
void modifyInputParam(List x, String paramType, String paramName, 
                      int cohort, double newValue,
                      bool message) {
  List control = x["control"];
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  CharacterVector cohNames = cohorts.attr("row.names");
  
  String transpirationMode = control["transpirationMode"];
  if(paramName=="LAI_live") {
    double old = getInputParamValue(x, "above", "LAI_live", cohort);
    double f = newValue/old;
    if(message) modifyMessage("LAI_live", cohNames[cohort], newValue);
    modifyInputParamSingle(x, "above", "LAI_live", cohort, newValue);
    if(message) multiplyMessage("LAI_expanded", cohNames[cohort], f);
    multiplyInputParamSingle(x, "above", "LAI_expanded", cohort, f);
    if(above.containsElementNamed("SA")) {
      if(message) multiplyMessage("SA", cohNames[cohort], f);
      multiplyInputParamSingle(x, "above", "SA", cohort, f);
    }
  } else if(paramName=="SLA") {
    double old = getInputParamValue(x, "paramsAnatomy", "SLA", cohort);
    double f = newValue/old;
    if(message) modifyMessage("SLA", cohNames[cohort], newValue);
    modifyInputParamSingle(x, "paramsAnatomy", "SLA", cohort, newValue);
    if(message) multiplyMessage("LAI_live", cohNames[cohort], f);
    multiplyInputParamSingle(x, "above", "LAI_live", cohort, f);
    if(message) multiplyMessage("LAI_expanded", cohNames[cohort], f);
    multiplyInputParamSingle(x, "above", "LAI_expanded", cohort, f);
    if(above.containsElementNamed("SA")) {
      if(message) multiplyMessage("SA", cohNames[cohort], f);
      multiplyInputParamSingle(x, "above", "SA", cohort, f);
    }
    if(x.containsElementNamed("paramsWaterStorage")) {
      if(message) multiplyMessage("Vleaf", cohNames[cohort], 1.0/f);
      multiplyInputParamSingle(x, "paramsWaterStorage", "Vleaf", cohort, 1.0/f);
    }
  } else if(paramName=="Al2As") {
    double old = getInputParamValue(x, "paramsAnatomy", "Al2As", cohort);
    double f = newValue/old;
    if(message) modifyMessage("Al2As", cohNames[cohort], newValue);
    modifyInputParamSingle(x, "paramsAnatomy", "Al2As", cohort, newValue);
    if(x.containsElementNamed("paramsWaterStorage")) {
      if(message) multiplyMessage("Vsapwood", cohNames[cohort], 1.0/f);
      multiplyInputParamSingle(x, "paramsWaterStorage", "Vsapwood", cohort, 1.0/f);
    }
    if(above.containsElementNamed("SA")) {
      if(message) multiplyMessage("SA", cohNames[cohort], 1.0/f);
      multiplyInputParamSingle(x, "above", "SA", cohort, 1.0/f);
    }
    if(transpirationMode=="Sperry") {
      if(message) multiplyMessage("VCstem_kmax", cohNames[cohort], 1.0/f);
      multiplyInputParamSingle(x, "paramsTranspiration", "VCstem_kmax", cohort, 1.0/f);
      if(message) multiplyMessage("VCroot_kmax", cohNames[cohort], 1.0/f);
      multiplyInputParamSingle(x, "paramsTranspiration", "VCroot_kmax", cohort, 1.0/f);
    }
  } else {
    if(message) modifyMessage(paramName, cohNames[cohort], newValue);
    modifyInputParamSingle(x, paramType, paramName, cohort, newValue);
  }
  if(transpirationMode=="Sperry") {
    if(message) Rcerr<< "[Message] Recalculating plant maximum conductances.\n";
    updatePlantKmax(x);
  }
  if(message) Rcerr<< "[Message] Updating below-ground parameters.\n";
  updateBelow(x);
  if(x.containsElementNamed("internalAllocation")) {
    DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
    DataFrame paramsTranspirationdf = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    DataFrame paramsAnatomydf = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
    if(message) Rcerr<< "[Message] Rebuilding allocation targets.\n";
    x["internalAllocation"]  = internalAllocationDataFrame(above, 
                                         belowdf, 
                                         paramsAnatomydf,
                                         paramsTranspirationdf,
                                         control);
  }
  if(x.containsElementNamed("internalCarbon")) {
    DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
    List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
    DataFrame paramsGrowthdf = Rcpp::as<Rcpp::DataFrame>(x["paramsGrowth"]);
    DataFrame paramsWaterStoragedf = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
    DataFrame paramsTranspirationdf = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
    DataFrame paramsAnatomydf = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
    if(message) Rcerr<< "[Message] Rebuilding internal carbon.\n";
    x["internalCarbon"]  = internalCarbonDataFrame(above, 
                                     belowdf,
                                     belowLayers,
                                     paramsAnatomydf,
                                     paramsWaterStoragedf,
                                     paramsGrowthdf,
                                     control);
  }
}
