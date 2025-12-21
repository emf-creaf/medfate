#include <Rcpp.h>
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
//' 
//' @return Functions \code{decomposition_moistureEffect}, \code{decomposition_pHEffect} and \code{decomposition_temperatureEffect} return
//' a scalar value representing a factor that should modify a decomposition rate. Function \code{decomposition_annualLitterDecompositionRate} 
//' directly returns a scalar value with the annual decomposition rate (yr-1). Function \code{decomposition_litterMetabolicFraction} returns
//' a scalar with the fraction of litter that corresponds to metabolic carbon.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' @seealso \code{\link{decomposition_DAYCENT}}
//' 
//' @references 
//' Meentemeyer (1978)
//' Kelly et al (2000)
//' 
//' @keywords internal
//' @name decomposition_annualLitterDecompositionRate
// [[Rcpp::export("decomposition_annualLitterDecompositionRate")]]
double annualLitterDecompositionRate(double AET, double lignin) {
  double ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin;
  return(ki);
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
  if(pool=="surface/metabolic") {
    a = 4.8; b=0.5; c=1.14; d = 0.7;
  } else if(pool=="soil/metabolic") {
      a = 4.8; b=0.5; c=1.14; d = 0.7;
  } else if(pool=="surface/structural") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="soil/structural") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="cwd/smallbranch") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="cwd/largewood") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="cwd/coarseroot") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="surface/active") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="soil/active") {
    a = 4.8; b= 0.5; c = 1.14; d = 0.7;
  } else if(pool=="surface/slow") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="soil/slow") {
    a = 4.0; b= 0.5; c = 1.1; d = 0.7;
  } else if(pool=="soil/passive") {
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

  K[DECOMPCOM_SURFACE_METABOLIC] = baseAnnualRates["surface/metabolic"]/365.25;
  K[DECOMPCOM_SOIL_METABOLIC] = baseAnnualRates["soil/metabolic"]/365.25;
  K[DECOMPCOM_SURFACE_ACTIVE] = baseAnnualRates["surface/active"]/365.25;
  K[DECOMPCOM_SOIL_ACTIVE] = baseAnnualRates["soil/active"]/365.25;
  K[DECOMPCOM_SURFACE_SLOW] = baseAnnualRates["surface/active"]/365.25;
  K[DECOMPCOM_SOIL_SLOW] = baseAnnualRates["soil/active"]/365.25;
  K[DECOMPCOM_SOIL_PASSIVE] = baseAnnualRates["soil/passive"]/365.25;
  
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
   pHeff = pHEffect(soilPH, "surface/metabolic");
   xi[DECOMPCOM_SURFACE_METABOLIC] = cdi * pHeff;
   //  metabolic litter (soil)
   pHeff = pHEffect(soilPH, "soil/metabolic");
   xi[DECOMPCOM_SOIL_METABOLIC] = cdi * pHeff * soilO2;
   //  active soil organic matter: SOM1 (surface)
   pHeff = pHEffect(soilPH, "surface/active");
   xi[DECOMPCOM_SURFACE_ACTIVE] = cdi * pHeff;
   //  active soil organic matter: SOM1 (SOIL)
   pHeff = pHEffect(soilPH, "soil/active");
   textureEff = 0.25 + 0.75 * (sand/100.0);
   xi[DECOMPCOM_SOIL_ACTIVE] = cdi * pHeff * soilO2 * textureEff * cultfac;
   //  slow soil organic matter: SOM2 (surface)
   // som2(surface) -> som1(surface)
   pHeff = pHEffect(soilPH, "surface/slow");
   double K_s21_to_s11 = K[DECOMPCOM_SURFACE_SLOW] * pHeff;
   // som2(surface) -> som2(soil): mixing
   double K_s21_to_s22 = Kmix;
   // total loss from som2(surface)
   double K_s21 = K_s21_to_s11 + K_s21_to_s22;
   // effective environmental scalar
   xi[DECOMPCOM_SURFACE_SLOW] = cdi * (K_s21 / K[DECOMPCOM_SURFACE_SLOW]);
   // slow soil organic matter: SOM2 (soil)
   pHeff = pHEffect(soilPH, "soil/slow");
   xi[DECOMPCOM_SOIL_SLOW] = cdi * pHeff * soilO2 * cultfac;
   //  passive soil organic matter: SOM3
   pHeff = pHEffect(soilPH, "soil/passive");
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


void DAYCENTlitterInner(NumericVector litterDecompositionOutput, 
                        DataFrame structuralLitter, DataFrame paramsDecomposition,
                        NumericVector baseAnnualRates,
                        double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                        double soilO2 = 1.0, double cultfac = 1.0,
                        double tstep = 1.0) {
  
  
  NumericVector structural_leaves = structuralLitter["leaves"];
  NumericVector structural_smallbranches = structuralLitter["smallbranches"];
  NumericVector structural_fineroots = structuralLitter["fineroots"];
  NumericVector structural_largewood = structuralLitter["largewood"];
  NumericVector structural_coarseroots = structuralLitter["coarseroots"];
  
  NumericVector LeafLignin = paramsDecomposition["LeafLignin"];
  int numCohorts = structuralLitter.nrow();

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
    pHeff = pHEffect(soilPH, "surface/structural");
    flig = LeafLignin[i]/100.0;
    k = (baseAnnualRates["surface/structural"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_leaves[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_leaves[i] -= loss;
  }

  // STRUCTURAL small branches
  flig = 0.25; //Lignin fraction for small branches
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "cwd/smallbranch");
    k = (baseAnnualRates["cwd/smallbranch"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_smallbranches[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_smallbranches[i] -= loss;
  }

  // STRUCTURAL large wood
  flig = 0.25; //Lignin fraction for large wood
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "cwd/largewood");
    k = (baseAnnualRates["cwd/largewood"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
    loss = structural_largewood[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_largewood[i] -= loss;
  }

  // STRUCTURAL coarse root
  flig = 0.25; //Lignin fraction for coarse root
  for(int i=0;i<numCohorts;i++) {
    pHeff = pHEffect(soilPH, "cwd/coarseroot");
    k = (baseAnnualRates["cwd/coarseroot"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig);
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
    pHeff = pHEffect(soilPH, "soil/structural");
    k = (baseAnnualRates["soil/structural"]/365.25)*tempEff*moistEff*pHeff*exp(-3.0*flig)*soilO2*cultfac;
    loss = structural_fineroots[i]*k*tstep;
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_ACTIVE] += loss*(1.0 - flig)*(1.0 - 0.45);
    litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_SLOW] += loss*flig*(1.0 - 0.30);
    litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION] += loss*(flig*0.30 + (1.0-flig)*0.45);
    structural_fineroots[i] -= loss;
  }
}

// [[Rcpp::export("decomposition_DAYCENTlitter")]]
NumericVector DAYCENTlitter(DataFrame structuralLitter, DataFrame paramsDecomposition,
                            NumericVector baseAnnualRates,
                            double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                            double soilO2 = 1.0, double cultfac = 1.0,
                            double tstep = 1.0) {
  NumericVector litterDecompositionOutput = communicationLitterDecomposition();
  DAYCENTlitterInner(litterDecompositionOutput, 
                     structuralLitter, paramsDecomposition, 
                     baseAnnualRates,
                     sand, clay, soilTemperature, soilMoisture, soilPH,
                     soilO2, cultfac,
                     tstep);
  return(litterDecompositionOutput);
}

double DAYCENTInner(List commDecomp,
                    DataFrame structuralLitter, NumericVector CENTURYPools,
                    DataFrame paramsDecomposition,
                    NumericVector baseAnnualRates, double annualTurnoverRate,
                    double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
                    double soilO2 = 1.0, double cultfac = 1.0,
                    double tstep = 1.0) {
  
  int npool = 7;
  NumericVector litterDecompositionOutput = commDecomp["ldo"];
  DAYCENTlitterInner(litterDecompositionOutput, 
                     structuralLitter, paramsDecomposition, 
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
  dC[DECOMPCOM_SURFACE_ACTIVE] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_ACTIVE];
  dC[DECOMPCOM_SOIL_ACTIVE] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_ACTIVE];
  dC[DECOMPCOM_SURFACE_SLOW] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SURFACE_SLOW];
  dC[DECOMPCOM_SOIL_SLOW] = litterDecompositionOutput[LITDECOMPCOM_TRANSFER_SOIL_SLOW];
  for(int i=0;i<npool;i++) {
    // carbon transfer from pool j to pool i
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        dC[i] = dC[i] + (1.0 - respf(i,j)) * pathf(i,j) * xi[j] * K[j] * CENTURYPools[j]*tstep;
      }
    }
    //  carbon loss from pool i
    dC[i] = dC[i] - xi[i] * K[i] * CENTURYPools[i]*tstep;
  }
  //   heterotrophic respiration
  double RH = litterDecompositionOutput[LITDECOMPCOM_FLUX_RESPIRATION];
  for(int i=0;i<npool;i++) {
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        RH = RH + respf(i,j) * pathf(i,j) * xi[j] * K[j] * CENTURYPools[j] * tstep;
      }
    }
  }
  // update pools
  for(int i=0;i<npool;i++) {
    CENTURYPools[i] += dC[i];
  }
  //return respiration
  return(RH);
}

//' DAYCENT decomposition
//' 
//' This function implements the DAYCENT carbon decomposition model, following the description in Bonan (2019) 
//' 
//' @param structuralLitter A data frame with structural carbon pools corresponding to plant cohorts, in g C/m2  (see \code{\link{growthInput}}).
//' @param CENTURYPools A named numeric vector with metabolic, active, slow and passive carbon pools for surface and soil, in g C/m2  (see \code{\link{growthInput}}).
//' @param paramsDecomposition A data frame of species-specific decomposition parameters (see \code{\link{growthInput}}).
//' @param baseAnnualRates A named vector of annual decomposition rates, in yr-1 (see \code{\link{defaultControl}}).
//' @param annualTurnoverRate Annual turnover rate, in yr-1  (see \code{\link{defaultControl}}).
//' @param sand,clay Soil texture (sand and sand) in percent volume (%). 
//' @param soilTemperature Soil temperature (in Celsius).
//' @param soilMoisture Soil moisture content, relative to saturation (0-1).
//' @param soilPH Soil pH (0-14).
//' @param soilO2 Soil oxygen factor (0-1).
//' @param cultfac Cultivation factor (0-1).
//' @param tstep Time step in days. By default, one day. For annual time steps, use \code{tstep = 365.25}.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @details Each call to function \code{decomposition_DAYCENT} conducts one time step of the DAYCENT
//' model and returns the heterotrophic respiration for that day. The function modifies input data \code{structuralLitter}
//' and \code{CENTURYPools} according to decomposition rates and carbon transfer rates. When used as part of \code{\link{growth}} simulations,
//' soil physical and chemical characteristics correspond to the uppermost soil layer.
//' 
//' @returns A scalar value with heterotrophic respiration, in g C/m2
//' 
//' @references
//' Bonan, G. (2019). Climate change and terrestrial ecosystem modeling. Cambridge University Press, Cambridge, UK.
//' 
//' @seealso \code{\link{decomposition_temperatureEffect}}, \code{\link{growthInput}}, \code{\link{growth}}
//' 
// [[Rcpp::export("decomposition_DAYCENT")]]
double DAYCENT(DataFrame structuralLitter, NumericVector CENTURYPools,
               DataFrame paramsDecomposition,
               NumericVector baseAnnualRates, double annualTurnoverRate,
               double sand, double clay, double soilTemperature, double soilMoisture, double soilPH, 
               double soilO2 = 1.0, double cultfac = 1.0,
               double tstep = 1.0) {
  List commDecomp = communicationDecomposition();
  return(DAYCENTInner(commDecomp,
                      structuralLitter, CENTURYPools,
                      paramsDecomposition,
                      baseAnnualRates, annualTurnoverRate,
                      sand, clay, soilTemperature, soilMoisture, soilPH, 
                      soilO2, cultfac,
                      tstep));
}
