#include <Rcpp.h>
#include "fuelmoisture.h"
#include "communication_structures.h"
using namespace Rcpp;


//' Low-level decomposition functions
//' 
//' Functions related to litter and soil carbon decomposition processes
//' 
//' @param AET Actual evapotranspiration (mm)
//' @param lignin Lignin percent
//' 
//' @details
//' Function \code{decomposition_moistureEffect} follows Kelly et al. (2000) 
//' Function \code{decomposition_snagFallProbability} follows Vanderwell et al. (2006) 
//' 
//' @return Functions \code{decomposition_moistureEffect}, \code{decomposition_pHEffect} and \code{decomposition_temperatureEffect} return
//' a scalar value representing a factor that should modify a decomposition rate. Function \code{decomposition_annualLitterDecompositionRate} 
//' directly returns a scalar value with the annual decomposition rate (yr-1). Function \code{decomposition_litterMetabolicFraction} returns
//' a scalar with the fraction of litter that corresponds to metabolic carbon.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' 
//' @seealso \code{\link{decomposition_DAYCENT}}
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' Vanderwel et al. (2006) Snag dynamics in partially harvested and unmanaged northern hardwood forests. Canadian Journal of Forest Research 36: 2769-2779.
//' 
//' Meentemeyer (1978)
//' 
//' Kelly et al (2000)
//' 
//' @keywords internal
//' @name decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_annualLitterDecompositionRate")]]
double annualLitterDecompositionRate(double AET, double lignin) {
  double ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin;
  return(ki);
}

//' @param DBH Diameter at breast height
//' @param decayClass Decay class, from 1 to 5
//' @param durabilityEffect Effect of wood durability
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_snagFallProbability")]]
double snagFallProbability(double DBH, int decayClass, double durabilityEffect = 0.0) {
  double gamma_dur[5];
  gamma_dur[0] = 0.0;
  gamma_dur[1] = 0.177;
  gamma_dur[2] = 0.542;
  gamma_dur[3] = 0.702;
  gamma_dur[4] = 0.528;
  
  double lnDBH = log(DBH);
  double lp  = 5.691 + durabilityEffect + gamma_dur[decayClass - 1] - 3.777*lnDBH + 0.531*pow(lnDBH,2.0);
  double p5 = exp(lp)/(1.0 + exp(lp));
  double p1 = 1.0 - pow(1.0 - p5, 1.0/5.0);
  return(p1);
}


//' @param ligninPercent lignin content (% of dry)
//' @param Nmass  nitrogen content (mg N / g dry)
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_litterMetabolicFraction")]]
double litterMetabolicFraction(double ligninPercent, double Nmass) {
  double fnit = Nmass/1000.0; //to g N/g dry
  double flig = ligninPercent/100.0; //to g lignin /g dry
  double rlig2n = flig/fnit;
  double fmet = 0.85 - 0.013 * rlig2n;
  return(fmet);
}

// [[Rcpp::export(".decomposition_addLeafTwigLitter")]]
void addLeafTwigLitter(String species_litter, double leaf_litter, double twig_litter,
                       DataFrame litter, 
                       DataFrame paramsLitterDecomposition,
                       NumericVector SOC) {
  
  NumericVector structural_litter_leaves = litter["Leaves"];
  NumericVector structural_litter_twigs = litter["Twigs"];
  NumericVector Nleaf = paramsLitterDecomposition["Nleaf"];
  NumericVector LeafLignin = paramsLitterDecomposition["LeafLignin"];
  
  CharacterVector Species = litter["Species"];
  int row = -1;
  for(int j=0;j<Species.length();j++) if(Species[j]==species_litter) row = j;
  if(row !=-1) {
    double fmet = litterMetabolicFraction(LeafLignin[row], Nleaf[row]);
    //Distribute between metabolic and structural
    SOC["SurfaceMetabolic"] = SOC["SurfaceMetabolic"] + leaf_litter*fmet;
    structural_litter_leaves[row] += leaf_litter*(1.0 - fmet);
    structural_litter_twigs[row] += twig_litter;
  }
}

// [[Rcpp::export(".decomposition_addSmallBranchLitter")]]
void addSmallBranchLitter(String species_litter, double smallbranch_litter, 
                          DataFrame litter) {
  NumericVector structural_litter_smallbranches = litter["SmallBranches"];
  // All small branch goes to structural litter compartment
  CharacterVector Species = litter["Species"];
  int row = -1;
  for(int j=0;j<Species.length();j++) if(Species[j]==species_litter) row = j;
  if(row !=-1) {
    structural_litter_smallbranches[row] += smallbranch_litter;
  }
}

// [[Rcpp::export(".decomposition_addLargeWoodLitter")]]
void addLargeWoodLitter(String species_litter, double largewood_litter, 
                        DataFrame litter) {
  NumericVector structural_litter_largewood = litter["LargeWood"];
  // All small branch goes to structural litter compartment
  CharacterVector Species = litter["Species"];
  int row = -1;
  for(int j=0;j<Species.length();j++) if(Species[j]==species_litter) row = j;
  if(row !=-1) {
    structural_litter_largewood[row] += largewood_litter;
  }
}

// [[Rcpp::export(".decomposition_addCoarseRootLitter")]]
void addCoarseRootLitter(String species_litter, double coarsewood_litter, 
                        DataFrame litter) {
  NumericVector structural_litter_coarseroots = litter["CoarseRoots"];
  // All small branch goes to structural litter compartment
  CharacterVector Species = litter["Species"];
  int row = -1;
  for(int j=0;j<Species.length();j++) if(Species[j]==species_litter) row = j;
  if(row !=-1) {
    structural_litter_coarseroots[row] += coarsewood_litter;
  }
}

// [[Rcpp::export(".decomposition_addFineRootLitter")]]
void addFineRootLitter(String species_litter, double fineroot_litter, 
                       DataFrame litter, 
                       DataFrame paramsLitterDecomposition,
                       NumericVector SOC) {
  NumericVector structural_litter_fineroots = litter["FineRoots"];
  NumericVector Nfineroot = paramsLitterDecomposition["Nfineroot"];

  CharacterVector Species = litter["Species"];
  int row = -1;
  for(int j=0;j<Species.length();j++) if(Species[j]==species_litter) row = j;
  if(row !=-1) {
    double fmet = litterMetabolicFraction(34.9, Nfineroot[row]); //Fine root lignin fraction 0.349
    //Distribute between metabolic and structural
    SOC["SoilMetabolic"] = SOC["SoilMetabolic"] + fineroot_litter*fmet;
    structural_litter_fineroots[row] += fineroot_litter*(1.0 - fmet);
  }
}
  


double pHEffect(double x, double a, double b, double c, double d) {
  double pi = 3.141592;
  double pHeff = b + (c / pi) * atan(d * (x - a) * pi);
  pHeff = std::max(std::min(pHeff, 1.0), 0.0);
  return(pHeff);
}

//' @param x Soil water pH (0-14)
//' @param pool String indicating the decomposition pool
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_pHEffect")]]
double pHEffect(double x, String pool) {
  double a, b, c, d;
  if(pool=="SurfaceMetabolic") {
    a = 4.8; b=0.5; c=1.14; d = 0.7;
  } else if(pool=="SoilMetabolic") {
      a = 4.8; b=0.5; c=1.14; d = 0.7;
  } else if(pool=="Leaves") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="FineRoots") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="Twigs") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="SmallBranches") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="LargeWood") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="CoarseRoots") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="SurfaceActive") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="SoilActive") {
    a = 4.8; b= 0.5; c = 1.14; d = 0.7;
  } else if(pool=="SurfaceSlow") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="SoilSlow") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="SoilPassive") {
    a = 3.0; b= 0.5; c = 1.1; d = 0.7;
  } else {
    stop("Wrong carbon pool");
  }
  return(pHEffect(x, a, b, c, d));
}

//' @param sand,clay Soil texture values in percent volume.
//' @param soilMoisture Soil moisture content, relative to saturation.
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_moistureEffect")]]
double moistureEffect(double sand, double clay, double soilMoisture) {
  
  bool fineTexture = (sand < 20.0); //TO BE REVISED
  
  double a = 0.6;
  double b = 1.27;
  double c = 0.0012;
  double d = 2.84;
  if(!fineTexture) {
    a = 0.55;
    b = 1.7;
    c = -0.007;
    d = 3.22;
  }
  double e = d*(b - a)/(a - c);
  double f1 = pow((soilMoisture - b)/(a - b), e);
  double f2 = pow((soilMoisture - c)/(a - c), d);
  return(f1*f2);
}


//' @param soilTemperature Soil temperature (in Celsius).
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_temperatureEffect")]]
double temperatureEffect(double soilTemperature) {
  double pi = 3.141592;
  double tEff = 0.56  + (1.46/pi) * atan(0.031*pi*(soilTemperature - 15.7)); //15.7 Celsius = 288.85 K
  return(tEff);
}

void updateBaseRates(List commDecomp,
                     NumericVector baseAnnualRates, double annualTurnoverRate) {
  NumericVector K = commDecomp["K"];

  K[DECOMPCOM_SURFACE_METABOLIC] = baseAnnualRates["SurfaceMetabolic"]/365.25;
  K[DECOMPCOM_SOIL_METABOLIC] = baseAnnualRates["SoilMetabolic"]/365.25;
  K[DECOMPCOM_SURFACE_ACTIVE] = baseAnnualRates["SurfaceActive"]/365.25;
  K[DECOMPCOM_SOIL_ACTIVE] = baseAnnualRates["SoilActive"]/365.25;
  K[DECOMPCOM_SURFACE_SLOW] = baseAnnualRates["SurfaceActive"]/365.25;
  K[DECOMPCOM_SOIL_SLOW] = baseAnnualRates["SoilActive"]/365.25;
  K[DECOMPCOM_SOIL_PASSIVE] = baseAnnualRates["SoilPassive"]/365.25;
  
  commDecomp["Kmix"] = annualTurnoverRate/365.25;
}
  
// Environmental scalar for each carbon pool adjusts base
// rate for soil temperature and soil moisture scalars
// (cdi) and additionally pH, lignin, texture, anaerobic,
// and cultivation
 //'  @param soilPH  soil pH
 //'  @param soilO2 effect of soil anaerobic conditions on decomposition (0-1)
 //'  @param sand,clay percent sand, clay
 //'  @param strlig lignin fraction: (1) surface and (2) soil structural litter (g lignin/g biomass)
 //'  @param cwdlig lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
 //'  @param cultfac effect of cultivation on decomposition (1:SOM1, 2:SOM2, 3:SOM3, 4:structural)
 //' 
 //' Updates
 //'    K_s21       ! rate constant: total loss from SOM2(surface), 1/sec
 //'    xi          ! environmental scalar
 void updateDecompositionRateScalars(List commDecomp, 
                                     double sand, double clay,
                                     double soilTemperature, double soilMoisture, double soilPH, 
                                     double soilO2, double cultfac) {
   
   NumericVector xi = commDecomp["xi"];
   NumericVector K = commDecomp["K"];
   double Kmix = commDecomp["Kmix"];
   
   double pHeff, textureEff;
   
   double tempEff = temperatureEffect(soilTemperature);
   double moistEff = moistureEffect(sand, clay, soilMoisture);
   double cdi = tempEff*moistEff;
   
   //  metabolic litter (surface)
   pHeff = pHEffect(soilPH, "SurfaceMetabolic");
   xi[DECOMPCOM_SURFACE_METABOLIC] = cdi * pHeff;
   //  metabolic litter (soil)
   pHeff = pHEffect(soilPH, "SoilMetabolic");
   xi[DECOMPCOM_SOIL_METABOLIC] = cdi * pHeff * soilO2;
   //  active soil organic matter: SOM1 (surface)
   pHeff = pHEffect(soilPH, "SurfaceActive");
   xi[DECOMPCOM_SURFACE_ACTIVE] = cdi * pHeff;
   //  active soil organic matter: SOM1 (SOIL)
   pHeff = pHEffect(soilPH, "SoilActive");
   textureEff = 0.25 + 0.75 * (sand/100.0);
   xi[DECOMPCOM_SOIL_ACTIVE] = cdi * pHeff * soilO2 * textureEff * cultfac;
   //  slow soil organic matter: SOM2 (surface)
   // som2(surface) -> som1(surface)
   pHeff = pHEffect(soilPH, "SurfaceSlow");
   double K_s21_to_s11 = K[DECOMPCOM_SURFACE_SLOW] * pHeff;
   // som2(surface) -> som2(soil): mixing
   double K_s21_to_s22 = Kmix;
   // total loss from som2(surface)
   double K_s21 = K_s21_to_s11 + K_s21_to_s22;
   // effective environmental scalar
   xi[DECOMPCOM_SURFACE_SLOW] = cdi * (K_s21 / K[DECOMPCOM_SURFACE_SLOW]);
   // slow soil organic matter: SOM2 (soil)
   pHeff = pHEffect(soilPH, "SoilSlow");
   xi[DECOMPCOM_SOIL_SLOW] = cdi * pHeff * soilO2 * cultfac;
   //  passive soil organic matter: SOM3
   pHeff = pHEffect(soilPH, "SoilPassive");
   xi[DECOMPCOM_SOIL_PASSIVE] = cdi * pHeff * soilO2 * cultfac;
   
   commDecomp["K_s21"] = K_s21;
 }


void updateCarbonTransferMatrices(List commDecomp, 
                                   double sand, double clay, double soilO2) {
   
   int npool = 7;
   NumericMatrix A = commDecomp["A"];
   NumericMatrix respf = commDecomp["respf"];
   NumericMatrix pathf = commDecomp["pathf"];
   double Kmix = commDecomp["Kmix"];
   double K_s21 = commDecomp["K_s21"];
   
   // anaerobic factor
   double fanerb = 1.0 + 5.0 * (1.0 - soilO2);
   
   // Updates pathf: fractional carbon flow from pool j to pool i
   pathf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_METABOLIC) = 1.0;
   pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_METABOLIC) = 1.0;
   pathf(DECOMPCOM_SURFACE_SLOW,DECOMPCOM_SURFACE_ACTIVE) = 1.0;
   pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE) = (0.003 + 0.032 * clay/100) * fanerb;
   pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE) = 1.0 - pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE);
   pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW) = Kmix / K_s21;
   pathf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_SLOW) = 1.0 - pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW);
   pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW) = (0.003 + 0.009 * clay/100) * fanerb;
   pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW) = 1.0 - pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW);
   pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_PASSIVE) = 1.0;
   
   // Updates respf: fractional respiration loss for carbon flow from pool j to pool i
   respf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_METABOLIC) = 0.55;
   respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_METABOLIC) = 0.55;
   respf(DECOMPCOM_SURFACE_SLOW,DECOMPCOM_SURFACE_ACTIVE) = 0.60;
   respf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE) = 0.0;
   respf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE) = (0.17 + 0.68 * sand/100.0) / pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE);
   respf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW) = 0.0;
   respf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_SLOW) = 0.55;
   respf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW) = 0.0;
   respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW) = 0.55 / pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW);
   respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_PASSIVE) = 0.55 * soilO2;
   
   // Update carbon transfer matrix: fractional carbon flow from pool j that enters pool i
   for(int i=0;i<npool;i++) {
     A(i,i) = -1.0; 
     for(int j=0;j<npool;j++) {
       if(j!=i) {
         A(i,j) = pathf(i,j) * (1.0 - respf(i,j));
       }
     }
   }
 }

void DAYCENTsnagsInner(NumericVector snagDecompositionOutput, 
                       DataFrame snags,
                       NumericVector baseAnnualRates,
                       double airTemperature, double airRelativeHumidity, 
                       double tstep = 1.0) {
  
  NumericVector structural_smallbranches = snags["SmallBranches"];
  NumericVector structural_largewood = snags["LargeWood"];
  int numCohorts = snags.nrow();
  
  // Reset output
  snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_ACTIVE] = 0.0;
  snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_SLOW] = 0.0;
  snagDecompositionOutput[SNAGDECOMPCOM_FLUX_RESPIRATION] = 0.0;
  
  
  //Temperature effect
  double tempEff = temperatureEffect(airTemperature);
  //Moisture effect for fine texture
  double fuelMaximumMoisture = 100.0; //% dry Depends on wood density
  double fuelMoisture = EMCSimard(airTemperature, airRelativeHumidity); //Ranges between 0 and 30% dry
  double moistEff = moistureEffect(0, 80, fuelMoisture/fuelMaximumMoisture);
  
  double k, flig, loss;
  
  // STRUCTURAL small branches
  flig = 0.25; //Lignin fraction for small branches
  for(int i=0;i<numCohorts;i++) {
    k = (baseAnnualRates["SmallBranches"]/365.25)*tempEff*moistEff*exp(-3.0*flig);
    loss = structural_smallbranches[i]*k*tstep;
    snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    snagDecompositionOutput[SNAGDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_smallbranches[i] -= loss;
  }
  
  // STRUCTURAL large wood
  flig = 0.25; //Lignin fraction for large wood
  for(int i=0;i<numCohorts;i++) {
    k = (baseAnnualRates["LargeWood"]/365.25)*tempEff*moistEff*exp(-3.0*flig);
    loss = structural_largewood[i]*k*tstep;
    snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    snagDecompositionOutput[SNAGDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_largewood[i] -= loss;
  }
}

//' @rdname decomposition_DAYCENT
//' @keywords internal
// [[Rcpp::export("decomposition_DAYCENTsnags")]]
NumericVector DAYCENTsnags(DataFrame snags, 
                           NumericVector baseAnnualRates,
                           double airTemperature, double airRelativeHumidity, 
                           double tstep = 1.0) {
  NumericVector snagDecompositionOutput = communicationSnagDecomposition();
  DAYCENTsnagsInner(snagDecompositionOutput, 
                    snags,  
                    baseAnnualRates,
                    airTemperature, airRelativeHumidity, 
                    tstep);
  return(snagDecompositionOutput);
}

void DAYCENTlitterInner(NumericVector litterDecompositionOutput, 
                        DataFrame litter, DataFrame paramsLitterDecomposition,
                        NumericVector baseAnnualRates,
                        double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                        double soilO2 = 1.0, double cultfac = 1.0,
                        double tstep = 1.0) {
  
  
  NumericVector structural_leaves = litter["Leaves"];
  NumericVector structural_twigs = litter["Twigs"];
  NumericVector structural_smallbranches = litter["SmallBranches"];
  NumericVector structural_largewood = litter["LargeWood"];
  NumericVector structural_coarseroots = litter["CoarseRoots"];
  NumericVector structural_fineroots = litter["FineRoots"];
  
  NumericVector LeafLignin = paramsLitterDecomposition["LeafLignin"];
  int numCohorts = litter.nrow();

  // Reset output
  litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] = 0.0;
  litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] = 0.0;
  litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_ACTIVE] = 0.0;
  litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_SLOW] = 0.0;
  litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] = 0.0;

  //Combined effect of temperature and moisture
  double pHeff;
  

  double tempEff = temperatureEffect(soilTemperature);
  double moistEff = moistureEffect(sand, clay, soilMoisture);

  double k, flig, loss;
  
  // STRUCTURAL leaves
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "Leaves");
    flig = LeafLignin[i]/100.0;
    k = (baseAnnualRates["Leaves"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_leaves[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_leaves[i] -= loss;
  }
  
  // STRUCTURAL twigs
  flig = 0.25; //Lignin fraction for small branches
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "Twigs");
    k = (baseAnnualRates["Twigs"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_smallbranches[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_smallbranches[i] -= loss;
  }
  
  // STRUCTURAL small branches
  flig = 0.25; //Lignin fraction for small branches
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "SmallBranches");
    k = (baseAnnualRates["SmallBranches"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_smallbranches[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_smallbranches[i] -= loss;
  }

  // STRUCTURAL large wood
  flig = 0.25; //Lignin fraction for large wood
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "LargeWood");
    k = (baseAnnualRates["LargeWood"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_largewood[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_largewood[i] -= loss;
  }

  // STRUCTURAL coarse root
  flig = 0.25; //Lignin fraction for coarse root
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "CoarseRoots");
    k = (baseAnnualRates["CoarseRoots"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_coarseroots[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.55);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.55);
    structural_coarseroots[i] -= loss;
  }
  
  // STRUCTURAL fineroots
  //Fine root lignin fraction 0.349
  flig = 0.349;
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "FineRoots");
    k = (baseAnnualRates["FineRoots"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig)*soilO2*cultfac;
    loss = structural_fineroots[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_fineroots[i] -= loss;
  }
}


//' @rdname decomposition_DAYCENT
//' @keywords internal
// [[Rcpp::export("decomposition_DAYCENTlitter")]]
NumericVector DAYCENTlitter(DataFrame litter, DataFrame paramsLitterDecomposition,
                            NumericVector baseAnnualRates,
                            double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                            double soilO2 = 1.0, double cultfac = 1.0,
                            double tstep = 1.0) {
  NumericVector litterDecompositionOutput = communicationLitterDecomposition();
  DAYCENTlitterInner(litterDecompositionOutput, 
                     litter, paramsLitterDecomposition, 
                     baseAnnualRates,
                     sand, clay, soilTemperature, soilMoisture, soilPH,
                     soilO2, cultfac,
                     tstep);
  return(litterDecompositionOutput);
}

double DAYCENTInner(List commDecomp,
                    DataFrame snags, DataFrame litter, NumericVector SOC,
                    DataFrame paramsLitterDecomposition,
                    NumericVector baseAnnualRates, double annualTurnoverRate,
                    double airTemperature, double airRelativeHumidity, 
                    double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                    double soilO2 = 1.0, double cultfac = 1.0,
                    double tstep = 1.0) {
  
  int npool = 7;
  NumericVector snagDecompositionOutput = commDecomp["sdo"];
  DAYCENTsnagsInner(snagDecompositionOutput, 
                    snags, 
                    baseAnnualRates,
                    airTemperature, airRelativeHumidity,
                    tstep);
  
  NumericVector litterDecompositionOutput = commDecomp["ldo"];
  DAYCENTlitterInner(litterDecompositionOutput, 
                     litter, paramsLitterDecomposition, 
                     baseAnnualRates,
                     sand, clay, soilTemperature, soilMoisture, soilPH,
                     tstep);
  
  updateBaseRates(commDecomp, baseAnnualRates, annualTurnoverRate);
  updateDecompositionRateScalars(commDecomp, 
                                 sand, clay,
                                 soilTemperature, soilMoisture, soilPH, 
                                 soilO2, cultfac);
  updateCarbonTransferMatrices(commDecomp, 
                               sand, clay, soilO2);

  NumericMatrix respf = commDecomp["respf"];
  NumericMatrix pathf = commDecomp["pathf"];
  NumericVector xi = commDecomp["xi"];
  NumericVector K = commDecomp["K"];
  
  NumericVector dC(npool, 0.0);
  dC[DECOMPCOM_SURFACE_ACTIVE] = snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_ACTIVE] + litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE];
  dC[DECOMPCOM_SURFACE_SLOW] = snagDecompositionOutput[SNAGDECOMPCOM_TRANSFER_SURFACE_SLOW] + litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW];
  dC[DECOMPCOM_SOIL_ACTIVE] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_ACTIVE];
  dC[DECOMPCOM_SOIL_SLOW] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_SLOW];
  for(int i=0;i<npool;i++) {
    // carbon transfer from pool j to pool i
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        dC[i] = dC[i] + (1.0 - respf(i,j)) * pathf(i,j) * xi[j] * K[j] * SOC[j]*tstep;
      }
    }
    //  carbon loss from pool i
    dC[i] = dC[i] - xi[i] * K[i] * SOC[i]*tstep;
  }
  //   heterotrophic respiration
  double RH = snagDecompositionOutput[SNAGDECOMPCOM_FLUX_RESPIRATION] + litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION];
  for(int i=0;i<npool;i++) {
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        RH = RH + respf(i,j) * pathf(i,j) * xi[j] * K[j] * SOC[j] * tstep;
      }
    }
  }
  // update pools
  for(int i=0;i<npool;i++) {
    SOC[i] += dC[i];
  }
  //return respiration
  return(RH);
}

//' DAYCENT decomposition
//' 
//' Functions implementing a modification of the DAYCENT carbon decomposition model, inspired by the description given in Chapter 18 of Bonan (2019).
//' Functions \code{decompositionDAYCENTsnags} and \code{decompositionDAYCENTlitter} conduct snag and litter decomposition, respectively, 
//' whereas function \code{decomposition_DAYCENT} performs the whole model for carbon decomposition.
//' 
//' @param snags A data frame with dead standing (snag) cohort information (see \code{\link{growthInput}}).
//' @param litter A data frame with aboveground and belowground structural carbon pools corresponding to plant cohorts, in g C/m2  (see \code{\link{growthInput}}).
//' @param SOC A named numeric vector with metabolic, active, slow and passive carbon pools for surface and soil, in g C/m2  (see \code{\link{growthInput}}).
//' @param paramsLitterDecomposition A data frame of species-specific litter decomposition parameters (see \code{\link{growthInput}}).
//' @param baseAnnualRates A named vector of annual decomposition rates, in yr-1 (see \code{\link{defaultControl}}).
//' @param annualTurnoverRate Annual turnover rate, in yr-1  (see \code{\link{defaultControl}}).
//' @param sand,clay Soil texture (sand and sand) in percent volume (%). 
//' @param airTemperature Mean daily air temperature (in Celsius).
//' @param airRelativeHumidity Mean daily relative humidity (%).
//' @param soilTemperature Soil temperature (in Celsius).
//' @param soilMoisture Soil moisture content, relative to saturation (0-1).
//' @param soilPH Soil pH (0-14).
//' @param soilO2 Soil oxygen factor (0-1).
//' @param cultfac Cultivation factor (0-1).
//' @param tstep Time step in days. By default, one day. For annual time steps, use \code{tstep = 365.25}.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @details Each call to functions \code{decomposition_DAYCENTlitter} or \code{decomposition_DAYCENTsnags} conducts one time step of the snag or litter dynamics, respectively. 
//' Each call to function \code{decomposition_DAYCENT} conducts one time step of the whole DAYCENT 
//' model and returns the heterotrophic respiration for that day. 
//' 
//' \emph{IMPORTANT NOTE}: Decomposition functions modify the input data (i.e. \code{snags}, \code{litter} and/or \code{SOC}) according to decomposition rates and carbon transfer rates. When used as part of \code{\link{growth}} simulations,
//' soil physical and chemical parameters correspond to the uppermost soil layer.
//' 
//' @returns 
//' Function \code{decomposition_DAYCENTsnags} returns a vector containing transfer carbon flows to SOC pools and heterotrophic respiration from snag decomposition.
//' Function \code{decomposition_DAYCENTlitter} returns a vector containing transfer carbon flows to SOC pools and heterotrophic respiration from litter decomposition. 
//' Function \code{decomposition_DAYCENT} returns scalar value with heterotrophic respiration (snags + litter + soil), in g C/m2.
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' @seealso \code{\link{decomposition_temperatureEffect}}, \code{\link{growthInput}}, \code{\link{growth}}
//' 
//' @keywords internal
// [[Rcpp::export("decomposition_DAYCENT")]]
double DAYCENT(DataFrame snags, DataFrame litter, NumericVector SOC,
               DataFrame paramsLitterDecomposition,
               NumericVector baseAnnualRates, double annualTurnoverRate,
               double airTemperature, double airRelativeHumidity, 
               double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
               double soilO2 = 1.0, double cultfac = 1.0,
               double tstep = 1.0) {
  List commDecomp = communicationDecomposition();
  return(DAYCENTInner(commDecomp,
                      snags, litter, SOC,
                      paramsLitterDecomposition,
                      baseAnnualRates, annualTurnoverRate,
                      airTemperature, airRelativeHumidity,  
                      sand, clay, soilTemperature, soilMoisture, soilPH, 
                      soilO2, cultfac,
                      tstep));
}
