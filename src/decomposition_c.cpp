#include "fuelmoisture.h"
#include "communication_structures.h"
#include "modelInput_c.h"
#include "decomposition_c.h"
#include "control_c.h"
#include "fuelmoisture_c.h"

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
double annualLitterDecompositionRate_c(double AET, double lignin) {
  double ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin;
  return(ki);
}

//' @param DBH Diameter at breast height
//' @param decayClass Decay class, from 1 to 5
//' @param durabilityEffect Effect of wood durability
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_snagFallProbability")]]
double snagFallProbability_c(double DBH, int decayClass, double durabilityEffect = 0.0) {
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
double litterMetabolicFraction_c(double ligninPercent, double Nmass) {
  double fnit = Nmass/1000.0; //to g N/g dry
  double flig = ligninPercent/100.0; //to g lignin /g dry
  double rlig2n = flig/fnit;
  double fmet = 0.85 - 0.013 * rlig2n;
  return(fmet);
}


double pHEffect_c(double x, double a, double b, double c, double d) {
  double pi = 3.141592;
  double pHeff = b + (c / pi) * atan(d * (x - a) * pi);
  pHeff = std::max(std::min(pHeff, 1.0), 0.0);
  return(pHeff);
}

double pHEffect_c(double x, const std::string& pool) {
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
  return(pHEffect_c(x, a, b, c, d));
}


//' @param sand,clay Soil texture values in percent volume.
//' @param soilMoisture Soil moisture content, relative to saturation.
//' 
//' @rdname decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_moistureEffect")]]
double moistureEffect_c(double sand, double clay, double soilMoisture) {
  
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
double temperatureEffect_c(double soilTemperature) {
  double pi = 3.141592;
  double tEff = 0.56  + (1.46/pi) * atan(0.031*pi*(soilTemperature - 15.7)); //15.7 Celsius = 288.85 K
  return(tEff);
}


int getRowLitter_c(std::string& species_litter, InternalLitter& litter) {
  int row = -1;
  int nlitter = litter.Species.size();
  for(int j=0;j<nlitter;j++) if(litter.Species[j]==species_litter) row = j;
  return(row);
}
void addLeafTwigLitter_c(std::string& species_litter, double leaf_litter, double twig_litter,
                         InternalLitter& litter, 
                         LitterDecompositionParams& paramsLitterDecomposition,
                         InternalSOC& SOC) {
  int row = getRowLitter_c(species_litter, litter);
  if(row !=-1) {
    double fmet = litterMetabolicFraction_c(paramsLitterDecomposition.LeafLignin[row], paramsLitterDecomposition.Nleaf[row]);
    //Distribute between metabolic and structural
    SOC.SurfaceMetabolic += leaf_litter*fmet;
    litter.Leaves[row] += leaf_litter*(1.0 - fmet);
    litter.Twigs[row] += twig_litter;
  }
}
void addSmallBranchLitter_c(std::string& species_litter, double smallbranch_litter, 
                            InternalLitter& litter) {
  int row = getRowLitter_c(species_litter, litter);
  if(row !=-1) {
    litter.SmallBranches[row] += smallbranch_litter;
  }
}
void addLargeWoodLitter_c(std::string& species_litter, double largewood_litter, 
                          InternalLitter& litter) {
  // All small branch goes to structural litter compartment
  int row = getRowLitter_c(species_litter, litter);
  if(row !=-1) {
    litter.LargeWood[row] += largewood_litter;
  }
}

void addCoarseRootLitter_c(std::string& species_litter, double coarsewood_litter, 
                           InternalLitter& litter) {
  // All small branch goes to structural litter compartment
  int row = getRowLitter_c(species_litter, litter);
  if(row !=-1) {
    litter.CoarseRoots[row] += coarsewood_litter;
  }
}

void addFineRootLitter_c(std::string& species_litter, double fineroot_litter, 
                         InternalLitter& litter, 
                         LitterDecompositionParams& paramsLitterDecomposition,
                         InternalSOC& SOC) {
  int row = getRowLitter_c(species_litter, litter);
  if(row !=-1) {
    double fmet = litterMetabolicFraction_c(34.9, paramsLitterDecomposition.Nfineroot[row]); //Fine root lignin fraction 0.349
    //Distribute between metabolic and structural
    SOC.SoilMetabolic += fineroot_litter*fmet;
    litter.FineRoots[row] += fineroot_litter*(1.0 - fmet);
  }
}

void updateBaseRates_c(Decomposition_COMM& DECcomm,
                       DecompositionAnnualBaseRates& baseAnnualRates, double annualTurnoverRate) {

  DECcomm.K[DECOMPCOM_SURFACE_METABOLIC] = baseAnnualRates.SurfaceMetabolic/365.25;
  DECcomm.K[DECOMPCOM_SOIL_METABOLIC] = baseAnnualRates.SoilMetabolic/365.25;
  DECcomm.K[DECOMPCOM_SURFACE_ACTIVE] = baseAnnualRates.SurfaceActive/365.25;
  DECcomm.K[DECOMPCOM_SOIL_ACTIVE] = baseAnnualRates.SoilActive/365.25;
  DECcomm.K[DECOMPCOM_SURFACE_SLOW] = baseAnnualRates.SurfaceSlow/365.25;
  DECcomm.K[DECOMPCOM_SOIL_SLOW] = baseAnnualRates.SoilSlow/365.25;
  DECcomm.K[DECOMPCOM_SOIL_PASSIVE] = baseAnnualRates.SoilPassive/365.25;
  
  DECcomm.Kmix = annualTurnoverRate/365.25;
}

// Environmental scalar for each carbon pool adjusts base
// rate for soil temperature and soil moisture scalars
// (cdi) and additionally pH, lignin, texture, anaerobic,
// and cultivation
//  @param soilPH  soil pH
//  @param soilO2 effect of soil anaerobic conditions on decomposition (0-1)
//  @param sand,clay percent sand, clay
//  @param strlig lignin fraction: (1) surface and (2) soil structural litter (g lignin/g biomass)
//  @param cwdlig lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
//  @param cultfac effect of cultivation on decomposition (1:SOM1, 2:SOM2, 3:SOM3, 4:structural)
// 
// Updates
//    K_s21       ! rate constant: total loss from SOM2(surface), 1/sec
//    xi          ! environmental scalar
void updateDecompositionRateScalars_c(Decomposition_COMM& DECcomm, 
                                       double sand, double clay,
                                       double soilTemperature, double soilMoisture, double soilPH, 
                                       double soilO2, double cultfac) {
   double pHeff, textureEff;
   
   double tempEff = temperatureEffect_c(soilTemperature);
   double moistEff = moistureEffect_c(sand, clay, soilMoisture);
   double cdi = tempEff*moistEff;
   
   //  metabolic litter (surface)
   pHeff = pHEffect_c(soilPH, "SurfaceMetabolic");
   DECcomm.xi[DECOMPCOM_SURFACE_METABOLIC] = cdi * pHeff;
   //  metabolic litter (soil)
   pHeff = pHEffect_c(soilPH, "SoilMetabolic");
   DECcomm.xi[DECOMPCOM_SOIL_METABOLIC] = cdi * pHeff * soilO2;
   //  active soil organic matter: SOM1 (surface)
   pHeff = pHEffect_c(soilPH, "SurfaceActive");
   DECcomm.xi[DECOMPCOM_SURFACE_ACTIVE] = cdi * pHeff;
   //  active soil organic matter: SOM1 (SOIL)
   pHeff = pHEffect_c(soilPH, "SoilActive");
   textureEff = 0.25 + 0.75 * (sand/100.0);
   DECcomm.xi[DECOMPCOM_SOIL_ACTIVE] = cdi * pHeff * soilO2 * textureEff * cultfac;
   //  slow soil organic matter: SOM2 (surface)
   // som2(surface) -> som1(surface)
   pHeff = pHEffect_c(soilPH, "SurfaceSlow");
   double K_s21_to_s11 = DECcomm.K[DECOMPCOM_SURFACE_SLOW] * pHeff;
   // som2(surface) -> som2(soil): mixing
   double K_s21_to_s22 = DECcomm.Kmix;
   // total loss from som2(surface)
   DECcomm.K_s21 = K_s21_to_s11 + K_s21_to_s22;
   // effective environmental scalar
   DECcomm.xi[DECOMPCOM_SURFACE_SLOW] = cdi * (DECcomm.K_s21 / DECcomm.K[DECOMPCOM_SURFACE_SLOW]);
   // slow soil organic matter: SOM2 (soil)
   pHeff = pHEffect_c(soilPH, "SoilSlow");
   DECcomm.xi[DECOMPCOM_SOIL_SLOW] = cdi * pHeff * soilO2 * cultfac;
   //  passive soil organic matter: SOM3
   pHeff = pHEffect_c(soilPH, "SoilPassive");
   DECcomm.xi[DECOMPCOM_SOIL_PASSIVE] = cdi * pHeff * soilO2 * cultfac;
 }

void updateCarbonTransferMatrices_c(Decomposition_COMM& DECcomm, 
                                    double sand, double clay, double soilO2) {
  
  int npool = DECcomm.K.size();

  // anaerobic factor
  double fanerb = 1.0 + 5.0 * (1.0 - soilO2);
  
  // Updates pathf: fractional carbon flow from pool j to pool i
  DECcomm.pathf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_METABOLIC) = 1.0;
  DECcomm.pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_METABOLIC) = 1.0;
  DECcomm.pathf(DECOMPCOM_SURFACE_SLOW,DECOMPCOM_SURFACE_ACTIVE) = 1.0;
  DECcomm.pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE) = (0.003 + 0.032 * clay/100.0) * fanerb;
  DECcomm.pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE) = 1.0 - DECcomm.pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE);
  DECcomm.pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW) = DECcomm.Kmix / DECcomm.K_s21;
  DECcomm.pathf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_SLOW) = 1.0 - DECcomm.pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW);
  DECcomm.pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW) = (0.003 + 0.009 * clay/100.0) * fanerb;
  DECcomm.pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW) = 1.0 - DECcomm.pathf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW);
  DECcomm.pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_PASSIVE) = 1.0;
  
  // Updates respf: fractional respiration loss for carbon flow from pool j to pool i
  DECcomm.respf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_METABOLIC) = 0.55;
  DECcomm.respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_METABOLIC) = 0.55;
  DECcomm.respf(DECOMPCOM_SURFACE_SLOW,DECOMPCOM_SURFACE_ACTIVE) = 0.60;
  DECcomm.respf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_ACTIVE) = 0.0;
  DECcomm.respf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE) = (0.17 + 0.68 * sand/100.0) / DECcomm.pathf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SOIL_ACTIVE);
  DECcomm.respf(DECOMPCOM_SOIL_SLOW,DECOMPCOM_SURFACE_SLOW) = 0.0;
  DECcomm.respf(DECOMPCOM_SURFACE_ACTIVE,DECOMPCOM_SURFACE_SLOW) = 0.55;
  DECcomm.respf(DECOMPCOM_SOIL_PASSIVE,DECOMPCOM_SOIL_SLOW) = 0.0;
  DECcomm.respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW) = 0.55 / DECcomm.pathf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_SLOW);
  DECcomm.respf(DECOMPCOM_SOIL_ACTIVE,DECOMPCOM_SOIL_PASSIVE) = 0.55 * soilO2;
  
  // Update carbon transfer matrix: fractional carbon flow from pool j that enters pool i
  for(int i=0;i<npool;i++) {
    DECcomm.A(i,i) = -1.0; 
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        DECcomm.A(i,j) = DECcomm.pathf(i,j) * (1.0 - DECcomm.respf(i,j));
      }
    }
  }
}


void DAYCENTsnagsInner_c(SnagDecomposition_COMM& sdo,
                         InternalSnags& snags, LitterDecompositionParams& paramsLitterDecomposition,
                         DecompositionAnnualBaseRates& baseAnnualRates,
                         double airTemperature, double airRelativeHumidity,
                         double tstep) {
  int numSnagCohorts = snags.DBH.size();

  // Reset output
  sdo.transfer_surface_active = 0.0;
  sdo.transfer_surface_slow = 0.0;
  sdo.surface_flux_respiration = 0.0;

  //Temperature effect
  double tempEff = temperatureEffect_c(airTemperature);
  //Moisture effect for fine texture
  double fuelMaximumMoisture = 100.0; //% dry Depends on wood density
  double fuelMoisture = EMCSimard_c(airTemperature, airRelativeHumidity); //Ranges between 0 and 30% dry
  double moistEff = moistureEffect_c(0, 80, fuelMoisture/fuelMaximumMoisture);

  double k, flig, loss;

  // STRUCTURAL small branches
  for(int i=0;i<numSnagCohorts;i++) {
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.SmallBranches/365.25)*tempEff*moistEff*exp(-3.0*flig);
    loss = snags.SmallBranches[i]*k*tstep;
    sdo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    sdo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    sdo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    snags.SmallBranches[i] -= loss;
  }

  // STRUCTURAL large wood
  for(int i=0;i<numSnagCohorts;i++) {
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.LargeWood/365.25)*tempEff*moistEff*exp(-3.0*flig);
    loss = snags.LargeWood[i]*k*tstep;
    sdo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    sdo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    sdo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    snags.LargeWood[i] -= loss;
  }
}


void DAYCENTlitterInner_c(LitterDecomposition_COMM& ldo, 
                          InternalLitter& litter, LitterDecompositionParams& paramsLitterDecomposition,
                          DecompositionAnnualBaseRates& baseAnnualRates,
                          double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                          double soilO2, double cultfac,
                          double tstep) {
  int numLitterCohorts = litter.CoarseRoots.size();
  
  // Reset output
  ldo.transfer_surface_active = 0.0;
  ldo.transfer_surface_slow = 0.0;
  ldo.transfer_soil_active = 0.0;
  ldo.transfer_soil_slow = 0.0;
  ldo.surface_flux_respiration = 0.0;
  ldo.soil_flux_respiration = 0.0;
  
  //Combined effect of temperature and moisture
  double pHeff;
  double tempEff = temperatureEffect_c(soilTemperature);
  double moistEff = moistureEffect_c(sand, clay, soilMoisture);
  
  double k, flig, loss;
  
  // STRUCTURAL leaves
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "Leaves");
    flig = paramsLitterDecomposition.LeafLignin[i]/100.0;
    k = (baseAnnualRates.Leaves/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = litter.Leaves[i]*k*tstep;
    ldo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    ldo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    ldo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    litter.Leaves[i] -= loss;
    // if(i==0) Rcout<< structural_leaves[i] << " T" <<soilTemperature <<  " Teff"<< tempEff << " M"<< moistEff << " pH"<< pHeff << " K" << k << " Loss: "<< loss << "\n"; 
  }
  
  // STRUCTURAL twigs
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "Twigs");
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.Twigs/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = litter.Twigs[i]*k*tstep;
    ldo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    ldo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    ldo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    litter.Twigs[i] -= loss;
  }
  
  // STRUCTURAL small branches
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "SmallBranches");
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.SmallBranches/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = litter.SmallBranches[i]*k*tstep;
    ldo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    ldo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    ldo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    litter.SmallBranches[i] -= loss;
  }
  
  // STRUCTURAL large wood
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "LargeWood");
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.LargeWood/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = litter.LargeWood[i]*k*tstep;
    ldo.transfer_surface_active += loss*(1.0 - flig)*(1.0 - 0.45);
    ldo.transfer_surface_slow += loss*flig*(1.0 - 0.30);
    ldo.surface_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    litter.LargeWood[i] -= loss;
  }
  
  // STRUCTURAL coarse root
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "CoarseRoots");
    flig = paramsLitterDecomposition.WoodLignin[i]/100.0; //Lignin fraction for wood
    k = (baseAnnualRates.CoarseRoots/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = litter.CoarseRoots[i]*k*tstep;
    ldo.transfer_soil_active += loss*(1.0 - flig)*(1.0 - 0.55);
    ldo.transfer_soil_slow += loss*flig*(1.0 - 0.30);
    ldo.soil_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.55);
    litter.CoarseRoots[i] -= loss;
  }
  
  // STRUCTURAL fineroots
  for(int i=0;i<numLitterCohorts;i++) {
    pHeff = pHEffect_c(soilPH, "FineRoots");
    flig = paramsLitterDecomposition.FineRootLignin[i]/100.0; //Lignin fraction for fine roots
    k = (baseAnnualRates.FineRoots/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig)*soilO2*cultfac;
    loss = litter.FineRoots[i]*k*tstep;
    ldo.transfer_soil_active += loss*(1.0 - flig)*(1.0 - 0.45);
    ldo.transfer_soil_slow += loss*flig*(1.0 - 0.30);
    ldo.soil_flux_respiration += loss*(flig*0.30 + (1.0-flig)*0.45);
    litter.FineRoots[i] -= loss;
  }
}

double DAYCENTInner_c(Decomposition_COMM& DECcomm,
                      InternalSnags& snags, InternalLitter& litter, InternalSOC& SOC,
                      LitterDecompositionParams& paramsLitterDecomposition,
                      DecompositionAnnualBaseRates& baseAnnualRates, double annualTurnoverRate,
                      double airTemperature, double airRelativeHumidity, 
                      double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                      double soilO2, double cultfac,
                      double tstep) {
  
  int npool = DECcomm.K.size();
  double RH = 0.0;
  
  // Rcpp::Rcout << "n pools " << npool <<"\n";
  std::vector<double> SOC_v(npool);
  SOC_v[DECOMPCOM_SURFACE_METABOLIC] = SOC.SurfaceMetabolic;
  SOC_v[DECOMPCOM_SOIL_METABOLIC] = SOC.SoilMetabolic;
  SOC_v[DECOMPCOM_SURFACE_ACTIVE] = SOC.SurfaceActive;
  SOC_v[DECOMPCOM_SOIL_ACTIVE] = SOC.SoilActive;
  SOC_v[DECOMPCOM_SURFACE_SLOW] = SOC.SurfaceSlow;
  SOC_v[DECOMPCOM_SOIL_SLOW] = SOC.SoilSlow;
  SOC_v[DECOMPCOM_SOIL_PASSIVE] = SOC.SoilPassive;

  // Rcpp::Rcout << "Snags-";
  DAYCENTsnagsInner_c(DECcomm.sdo,
                      snags, paramsLitterDecomposition,
                      baseAnnualRates,
                      airTemperature, airRelativeHumidity,
                      tstep);

  // Rcpp::Rcout << "Litter-";
  DAYCENTlitterInner_c(DECcomm.ldo,
                       litter, paramsLitterDecomposition,
                       baseAnnualRates,
                       sand, clay, soilTemperature, soilMoisture, soilPH,
                       tstep);

  // Rcpp::Rcout << "Update rates-";
  updateBaseRates_c(DECcomm, baseAnnualRates, annualTurnoverRate);
  updateDecompositionRateScalars_c(DECcomm,
                                   sand, clay,
                                   soilTemperature, soilMoisture, soilPH,
                                   soilO2, cultfac);
  updateCarbonTransferMatrices_c(DECcomm,
                                 sand, clay, soilO2);

  std::vector<double> dC(npool, 0.0);
  dC[DECOMPCOM_SURFACE_ACTIVE] = DECcomm.sdo.transfer_surface_active + DECcomm.ldo.transfer_surface_active;
  dC[DECOMPCOM_SURFACE_SLOW] = DECcomm.sdo.transfer_surface_slow + DECcomm.ldo.transfer_surface_slow;
  dC[DECOMPCOM_SOIL_ACTIVE] = DECcomm.ldo.transfer_soil_active;
  dC[DECOMPCOM_SOIL_SLOW] = DECcomm.ldo.transfer_soil_slow;
  for(int i=0;i<npool;i++) {
    // carbon transfer from pool j to pool i
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        dC[i] = dC[i] + (1.0 - DECcomm.respf(i,j)) * DECcomm.pathf(i,j) * DECcomm.xi[j] * DECcomm.K[j] * SOC_v[j]*tstep;
      }
    }
    //  carbon loss from pool i
    dC[i] = dC[i] - DECcomm.xi[i] * DECcomm.K[i] * SOC_v[i]*tstep;
  }
  //   heterotrophic respiration
  RH += DECcomm.sdo.surface_flux_respiration + DECcomm.ldo.surface_flux_respiration + DECcomm.ldo.soil_flux_respiration;
  for(int i=0;i<npool;i++) {
    DECcomm.flux_respiration_pools[i] = 0.0;
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        DECcomm.flux_respiration_pools[i] += DECcomm.respf(i,j) * DECcomm.pathf(i,j) * DECcomm.xi[j] * DECcomm.K[j] * SOC_v[j] * tstep;
      }
    }
    RH = RH + DECcomm.flux_respiration_pools[i];
  }
  // update pools
  for(int i=0;i<npool;i++) {
    SOC_v[i] += dC[i];
  }
  SOC.SurfaceMetabolic = SOC_v[DECOMPCOM_SURFACE_METABOLIC];
  SOC.SoilMetabolic = SOC_v[DECOMPCOM_SOIL_METABOLIC];
  SOC.SurfaceActive = SOC_v[DECOMPCOM_SURFACE_ACTIVE];
  SOC.SoilActive = SOC_v[DECOMPCOM_SOIL_ACTIVE];
  SOC.SurfaceSlow = SOC_v[DECOMPCOM_SURFACE_SLOW];
  SOC.SoilSlow = SOC_v[DECOMPCOM_SOIL_SLOW];
  SOC.SoilPassive = SOC_v[DECOMPCOM_SOIL_PASSIVE];

  //return respiration
  return(RH);
}