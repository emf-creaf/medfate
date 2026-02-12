#include <RcppArmadillo.h>
#include "tissuemoisture_c.h"
#include "carbon_c.h"
#include "forestutils_c.h"
#include "biophysicsutils_c.h"
#include "meteoland/utils_c.hpp"

const double Rn = 0.008314472; // The perfect gas constant MPa·l/K·mol = kJ/K·mol

//Sains Malaysiana 44(7)(2015): 973–977
//Modeling of Sago Starch Hydrolysis Using Glucoamylase
// double starchHydrolysis(double starchConc) {
//   double ks_hydrolysis_starch_vol = 2.324; // = g·L-1
//   double vmax_hydrolysis_starch_vol = 0.424; // = g·L-1·min-1
//   starchConc = starchConc*starchMolarMass; // from mol·L-1 to g·L-1
//   double v = (vmax_hydrolysis_starch_vol*starchConc)/(ks_hydrolysis_starch_vol + starchConc);
//   v = v/(starchMolarMass*60.0); // from g·L-1·min to mol·L-1·s-1
//   return(v); //return mol·L-1·s-1
// }

double sugarStarchDynamics_c(double sugarConc, double starchConc,
                             double kmsyn, double vmaxsyn, double khyd, double eqSugarConc) {
  double sugarM = sugarConc*(18.0/1000.0); //from mol·L-1 to mol·mol-1
  double starchM = starchConc*(18.0/1000.0); //from mol·L-1 to mol·mol-1
  double STsyn = vmaxsyn*(sugarM/(kmsyn + sugarM));
  double SThyd = khyd*starchM;
  double dSdt;
  if(sugarConc > eqSugarConc) {
    dSdt = STsyn*(1000.0/18.0);
  } else { //Downregulate starch synthesis
    dSdt = -SThyd*(1000.0/18.0);
  }
  return(dSdt/(3600.0*24.0)); //return mol·l-1·s-1
}

//' Carbon-related functions
//' 
//' Low-level functions used in the calculation of carbon balance.
//' 
//' @param LAI Leaf area index.
//' @param N Density (ind·ha-1).
//' @param DBH Tree diameter (cm).
//' @param SLA Specific leaf area (mm2/mg = m2/kg).
//' @param leafDensity  Density of leaf tissue (dry weight over volume).
//' @param SA Sapwood area (cm2).
//' @param H Plant height (cm).
//' @param L Coarse root length (mm) for each soil layer.
//' @param V Proportion of fine roots in each soil layer.
//' @param woodDensity Wood density (dry weight over volume).
//' @param conduit2sapwood Proportion of sapwood corresponding to conducive elements (vessels or tracheids) as opposed to parenchymatic tissue.
//' @param osmoticWP Osmotic water potential (MPa).
//' @param temp Temperature (degrees Celsius).
//' @param nonSugarConc Concentration of inorganic solutes (mol/l).
//' @param sugarConc Concentration of soluble sugars (mol/l).
//' @param starchConc Concentration of starch (mol/l)
//' @param eqSugarConc Equilibrium concentration of soluble sugars (mol/l).
//' 
//' @return Values returned for each function are:
//' \itemize{
//'   \item{\code{carbon_leafStarchCapacity}: Capacity of storing starch in the leaf compartment (mol gluc/ind.).}
//'   \item{\code{carbon_leafStructuralBiomass}: Leaf structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sapwoodStarchCapacity}: Capacity of storing starch in the sapwood compartment (mol gluc/ind.).}
//'   \item{\code{carbon_sapwoodStructuralBiomass}: Sapwood structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sapwoodStructuralLivingBiomass}: Living sapwood (parenchyma) structural biomass (g dry/ind.)}
//'   \item{\code{carbon_sugarConcentration}: Sugar concentration (mol gluc/l)}
//'   \item{\code{carbon_osmoticWaterPotential}: Osmotic component of water potential (MPa)}
//'   \item{\code{carbon_relativeSapViscosity}: Relative viscosity of sapwood with respect to pure water (according to Forst et al. (2002)).}
//'   \item{\code{carbon_sugarStarchDynamicsLeaf}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
//'   \item{\code{carbon_sugarStarchDynamicsStem}: Rate of conversion from sugar to starch in leaf (mol gluc/l/s).}
//'   \item{\code{carbon_carbonCompartments}: A data frame with the size of compartments for each plant cohort, in the specified units.}
//' }
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references
//' Forst P, Wermer F, Delgado A (2002). On the pressure dependence of the viscosity of aqueous sugar solutions. Rheol Acta 41: 369–374 DOI 10.1007/s00397-002-0238-y
//' 
//' @seealso \code{\link{growth}}
//' 
//' @examples
//' #Load example plot plant data
//' data(exampleforest)
//' #Default species parameterization
//' data(SpParamsMED)
//' #Initialize control parameters
//' control <- defaultControl("Granier")
//' #Initialize soil with default soil params (4 layers)
//' examplesoil <- defaultSoilParams(4)
//' #Initialize model input
//' x1 <- growthInput(exampleforest, examplesoil, SpParamsMED, control)
//' 
//' # Estimate carbon compartments
//' carbon_carbonCompartments(x1)
//' 
//' @name carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarStarchDynamicsLeaf")]]
double sugarStarchDynamicsLeaf_c(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics_c(sugarConc, starchConc, 0.1, 0.3, 1, eqSugarConc));
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarStarchDynamicsStem")]]
double sugarStarchDynamicsStem_c(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics_c(sugarConc, starchConc, 0.1, 0.15, 0.4, eqSugarConc));
}
double sugarStarchDynamicsRoot_c(double sugarConc, double starchConc, double eqSugarConc) {
  return(sugarStarchDynamics_c(sugarConc, starchConc, 0.1, 0.6, 0.4, eqSugarConc));
}

/**
 * Van 't Hoff equation
 *  conc - mol/l 
 *  temp - deg C
 *  wp - MPa
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_osmoticWaterPotential")]]
double osmoticWaterPotential_c(double sugarConc, double temp, double nonSugarConc) {
  sugarConc = std::max(0.0, sugarConc); // To avoid positive osmotic water potential
  return(- (sugarConc + nonSugarConc)*Rn*(temp + 273.15));
}
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_sugarConcentration")]]
double sugarConcentration_c(double osmoticWP, double temp, double nonSugarConc) {
  return(- osmoticWP/(Rn*(temp + 273.15)) - nonSugarConc);
}


/**
* 
*  sugarConc - sugar concentration (mol/l)
*  temp - temperature (degrees C)
*/
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_relativeSapViscosity")]]
double relativeSapViscosity_c(double sugarConc, double temp) {
  double x = sugarConc*glucoseMolarMass/1e3; //from mol/l to g*cm-3
  double Tkelvin = temp + 273.15;
  double q0a = 1.12; //g*cm-3
  double q1 = -0.248;
  double Ea = 2.61; //kJ*mol-1 energy of activation
  double va = x/(q0a*exp(-1.0*Ea/(Rn*Tkelvin)));
  double relVisc = exp(va/(1.0 + q1*va)); // relative viscosity
  double relWat = waterDynamicViscosity_c(temp); 
  return(relWat*relVisc);
}

/**
 *  Turgor (MPa)
 *  conc - mol/l 
 *  temp - deg C
 *  psi - water potential (MPa)
 */
double turgor_c(double psi, double sugarConc, double temp, double nonSugarConc) {
  return(std::max(0.0, psi-osmoticWaterPotential_c(sugarConc,temp, nonSugarConc)));
}

/**
 * Leaf area in m2 · ind-1
 */
double leafArea_c(double LAI, double N) {
  if(N==0.0) return(0.0);
  return(10000.0*LAI/N);
}
/**
 * leaf volume in l
 */
double leafStorageVolume_c(double LAI, double N, double SLA, double leafDensity) {
  return(leafArea_c(LAI,N)*leafWaterCapacity_c(SLA, leafDensity)); 
}

/*
 * Leaf structural biomass in g dw · ind-1
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_leafStructuralBiomass")]]
double leafStructuralBiomass_c(double LAI, double N, double SLA) {
  return(1000.0*leafArea_c(LAI,N)/SLA);  
}

/*
 * Twig structural biomass in g dw · ind-1
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_twigStructuralBiomass")]]
double twigStructuralBiomass_c(double LAI, double N, double SLA, double r635) {
  return(leafStructuralBiomass_c(LAI, N, SLA)*(r635 - 1.0));  
}
/*
 * Leaf starch storage capacity in mol · ind-1
 * Up to 10% of leaf cell volume
 */
//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_leafStarchCapacity")]]
double leafStarchCapacity_c(double LAI, double N, double SLA, double leafDensity) {
  return(0.1*1000.0*leafStorageVolume_c(LAI,N,SLA,leafDensity)*starchDensity/starchMolarMass);
}

/**
 * aboveground wood volume in l = dm3
 * 
 * area - cm2
 * H - cm
 */
double abovegroundWoodVolume_c(double area, double H) {
  return(0.001*area*H);
}

double belowgroundWoodVolume_c(double area, const std::vector<double>& L, const std::vector<double>& V) {
  int nlayers = V.size();
  double vBelow = 0.0;
  for(int i=0;i<nlayers;i++) {
    vBelow += 0.001*area*V[i]*(L[i]/10.0);
  }
  return(vBelow);
}

double abovegroundSapwoodVolume_c(double SA, double H) {
  return(abovegroundWoodVolume_c(SA, H));
}
double belowgroundSapwoodVolume_c(double SA, const std::vector<double>& L, const std::vector<double>& V) {
  return(belowgroundWoodVolume_c(SA, L, V));
}
/**
 * sapwood volume in l = dm3
 * 
 * SA - cm2
 * H - cm
 * Z - mm
 */
double sapwoodVolume_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V) {
  double vAbove = abovegroundSapwoodVolume_c(SA, H);
  double vBelow = belowgroundSapwoodVolume_c(SA, L, V);
  return(vAbove+vBelow);
}

double abovegroundHeartwoodVolume_c(double DBH, double SA, double H) {
  double section_cm2 = pow(DBH/2.0, 2.0)*M_PI;
  double vAbove = abovegroundWoodVolume_c(section_cm2, H);
  double vAboveSap = abovegroundSapwoodVolume_c(SA, H);
  return(std::max(0.0,vAbove- vAboveSap));
}
double belowgroundHeartwoodVolume_c(double DBH, double SA, const std::vector<double>& L, const std::vector<double>& V) {
  double section_cm2 = pow(DBH/2.0, 2.0)*M_PI;
  double vBelow = belowgroundWoodVolume_c(section_cm2, L, V);
  double vBelowSap = belowgroundSapwoodVolume_c(SA, L, V);
  return(std::max(0.0,vBelow - vBelowSap));
}

/**
 * heartwood volume in l = dm3
 * 
 * DBH - cm
 * H - cm
 * Z - mm
 */
double heartwoodVolume(double DBH, double SA, double H, const std::vector<double>& L, const std::vector<double>& V) {
  double vAbove = abovegroundHeartwoodVolume_c(DBH, SA, H);
  double vBelow = belowgroundHeartwoodVolume_c(DBH, SA, L, V);
  return(vAbove+vBelow);
}

/**
 * sapwood storage volume in l
 *  
 *  SA - cm2
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
double sapwoodStorageVolume_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                            double woodDensity, double conduit2sapwood) { 
  double woodPorosity = (1.0- (woodDensity/1.54));
  return((1.0 - conduit2sapwood)*sapwoodVolume_c(SA,H,L,V)*woodPorosity);
}




//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_abovegroundSapwoodStructuralBiomass")]]
double abovegroundSapwoodStructuralBiomass_c(double SA, double H, double woodDensity) {
   return(1000.0*abovegroundSapwoodVolume_c(SA, H)*woodDensity);
}

double belowgroundSapwoodStructuralBiomass_c(double SA, const std::vector<double>& L, const std::vector<double>& V, 
                                             double woodDensity) {
  return(1000.0*belowgroundSapwoodVolume_c(SA, L, V)*woodDensity);
}

/**
 * sapwood structural biomass in g dw / ind
 * 
 *  SA - cm2·ind-1
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
double sapwoodStructuralBiomass_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                                double woodDensity) {
  double B_above = abovegroundSapwoodStructuralBiomass_c(SA,H,woodDensity);
  double B_below = belowgroundSapwoodStructuralBiomass_c(SA,L, V,woodDensity);
  return(B_above + B_below);
}

//' @rdname carbon
//' @keywords internal
// [[Rcpp::export("carbon_abovegroundHeartwoodStructuralBiomass")]]
double abovegroundHeartwoodStructuralBiomass_c(double DBH, double SA, double H, double woodDensity) {
  return(1000.0*abovegroundHeartwoodVolume_c(DBH, SA, H)*woodDensity);
}

double belowgroundHeartwoodStructuralBiomass_c(double DBH, double SA, const std::vector<double>& L, const std::vector<double>& V, 
                                             double woodDensity) {
  return(1000.0*belowgroundHeartwoodVolume_c(DBH, SA, L, V)*woodDensity);
}

double heartwoodStructuralBiomass_c(double DBH, double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                                    double woodDensity) {
  double B_above = abovegroundHeartwoodStructuralBiomass_c(DBH, SA,H,woodDensity);
  double B_below = belowgroundHeartwoodStructuralBiomass_c(DBH, SA,L, V,woodDensity);
  return(B_above + B_below);
}


/*
 *  Sapwood starch storage capacity in mol · ind-1
 *  Up to 50% of volume of non-conductive cells
 *  
 *  SA - cm2
 *  H - cm
 *  Z - mm
 *  woodDensity - g/cm3
 */
double sapwoodStarchCapacity_c(double SA, double H, const std::vector<double>& L, const std::vector<double>& V, 
                               double woodDensity, double conduit2sapwood) {
  return(0.5*1000.0*sapwoodStorageVolume_c(SA,H,L,V,woodDensity,conduit2sapwood)*starchDensity/starchMolarMass);
}



void fillCarbonCompartments_c(CarbonCompartments& cc, ModelInput& x, std::string& biomassUnits) {
  
  if((biomassUnits!="g_m2") && (biomassUnits !="g_ind") && (biomassUnits !="gC_m2")) {
    throw medfate::MedfateInternalError("Wrong biomass units"); 
  }
  
  //Aboveground parameters
  std::vector<double>& H = x.above.H;
  std::vector<double>& N = x.above.N;
  std::vector<double>& LAI_expanded = x.above.LAI_expanded;
  std::vector<double>& LAI_dead = x.above.LAI_dead;
  std::vector<double>& SA = x.above.SA;
  std::vector<double>& DBH = x.above.DBH;
  
  std::vector<std::string> ctype = cohortType_c(x.cohorts.CohortCode);

  arma::mat& V = x.belowLayers.V;
  arma::mat& L = x.belowLayers.L;
   
  std::vector<double>& SLA = x.paramsAnatomy.SLA;
  std::vector<double>& r635 = x.paramsAnatomy.r635;
  std::vector<double>& WoodDensity = x.paramsAnatomy.WoodDensity;
  std::vector<double>& LeafDensity = x.paramsAnatomy.LeafDensity;
  std::vector<double>& conduit2sapwood = x.paramsAnatomy.conduit2sapwood;
  std::vector<double>& WoodC = x.paramsGrowth.WoodC;
  
  int numCohorts = x.cohorts.CohortCode.size();
  int nlayers = x.soil.getNlayers();
  for(int j=0;j<numCohorts;j++){
    std::vector<double> Vj(nlayers);
    std::vector<double> Lj(nlayers);
    for(int l=0;l<nlayers;l++) {
      Vj[l] = V(j,l);
      Lj[l] = L(j,l);
    }
    cc.FineRootBiomass[j] = x.below.fineRootBiomass[j];
    cc.LeafStructuralBiomass[j] = leafStructuralBiomass_c(LAI_expanded[j],N[j],x.paramsAnatomy.SLA[j]);
    cc.TwigStructuralBiomass[j] = twigStructuralBiomass_c(LAI_expanded[j],N[j],SLA[j], x.paramsAnatomy.r635[j]);
    cc.DeadLeafStructuralBiomass[j] = leafStructuralBiomass_c(LAI_dead[j],N[j],SLA[j]);
    cc.DeadTwigStructuralBiomass[j] = twigStructuralBiomass_c(LAI_dead[j],N[j],SLA[j], r635[j]);
    cc.TwigLivingStructuralBiomass[j] = cc.TwigStructuralBiomass[j]*(1.0 - x.internalWater.StemPLC[j])*(1.0-conduit2sapwood[j]);
    cc.SapwoodStructuralBiomass[j] = sapwoodStructuralBiomass_c(SA[j], H[j], Lj,Vj, WoodDensity[j]);
    cc.SapwoodLivingStructuralBiomass[j] = cc.SapwoodStructuralBiomass[j]*(1.0 - x.internalWater.StemPLC[j])*(1.0-conduit2sapwood[j]);
    cc.AbovegroundWoodBiomass[j] = abovegroundSapwoodStructuralBiomass_c(SA[j], H[j], WoodDensity[j]);
    cc.BelowgroundWoodBiomass[j] = belowgroundSapwoodStructuralBiomass_c(SA[j], Lj, Vj, WoodDensity[j]);
    if(ctype[j] =="tree") {
      cc.HeartwoodStructuralBiomass[j] = heartwoodStructuralBiomass_c(DBH[j], SA[j], H[j], Lj, Vj, WoodDensity[j]);
      cc.AbovegroundWoodBiomass[j] += abovegroundHeartwoodStructuralBiomass_c(DBH[j], SA[j], H[j], WoodDensity[j]);
      cc.BelowgroundWoodBiomass[j] += belowgroundHeartwoodStructuralBiomass_c(DBH[j], SA[j], Lj,Vj, WoodDensity[j]);
    } else {
      cc.HeartwoodStructuralBiomass[j] = 0.0;
    }
    cc.LeafStorageVolume[j] = leafStorageVolume_c(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    cc.SapwoodStorageVolume[j] = sapwoodStorageVolume_c(SA[j], H[j], Lj, Vj, WoodDensity[j], conduit2sapwood[j]);
    cc.LeafStarchCapacity[j] = leafStarchCapacity_c(LAI_expanded[j],  N[j], SLA[j], LeafDensity[j]);
    cc.SapwoodStarchCapacity[j] =  sapwoodStarchCapacity_c(SA[j], H[j], Lj, Vj, WoodDensity[j], conduit2sapwood[j]);
    cc.LeafStarchMaximumConcentration[j] = cc.LeafStarchCapacity[j]/cc.LeafStorageVolume[j];
    cc.SapwoodStarchMaximumConcentration[j] = cc.SapwoodStarchCapacity[j]/cc.SapwoodStorageVolume[j];
    if(cc.LeafStorageVolume[j]==0.0) cc.LeafStarchMaximumConcentration[j] = 0.0;
    if(cc.SapwoodStorageVolume[j]==0.0) cc.SapwoodStarchMaximumConcentration[j] = 0.0;
    
    double labileMassLeaf = (x.internalCarbon.sugarLeaf[j]+x.internalCarbon.starchLeaf[j])*(glucoseMolarMass*cc.LeafStorageVolume[j]);
    double labileMassSapwood = (x.internalCarbon.sugarSapwood[j]+x.internalCarbon.starchSapwood[j])*(glucoseMolarMass*cc.SapwoodStorageVolume[j]);
    cc.LabileBiomass[j] = labileMassSapwood+labileMassLeaf;
    
    if(biomassUnits=="g_m2") {
      double f = N[j]/(10000.0);
      cc.LeafStructuralBiomass[j] = cc.LeafStructuralBiomass[j]*f;
      cc.TwigStructuralBiomass[j] = cc.TwigStructuralBiomass[j]*f;
      cc.TwigLivingStructuralBiomass[j] = cc.TwigLivingStructuralBiomass[j]*f;
      cc.DeadLeafStructuralBiomass[j] = cc.DeadLeafStructuralBiomass[j]*f;
      cc.DeadTwigStructuralBiomass[j] = cc.DeadTwigStructuralBiomass[j]*f;
      cc.SapwoodStructuralBiomass[j] = cc.SapwoodStructuralBiomass[j]*f;
      cc.SapwoodLivingStructuralBiomass[j] = cc.SapwoodLivingStructuralBiomass[j]*f;
      cc.HeartwoodStructuralBiomass[j] = cc.HeartwoodStructuralBiomass[j]*f;
      cc.AbovegroundWoodBiomass[j] = cc.AbovegroundWoodBiomass[j]*f;
      cc.BelowgroundWoodBiomass[j] = cc.BelowgroundWoodBiomass[j]*f;
      cc.FineRootBiomass[j] = cc.FineRootBiomass[j]*f;
      cc.LabileBiomass[j] = cc.LabileBiomass[j]*f;
    } else if(biomassUnits=="gC_m2") {
      double f = N[j]/(10000.0);
      cc.LeafStructuralBiomass[j] = cc.LeafStructuralBiomass[j]*f*leafCperDry;
      cc.TwigStructuralBiomass[j] = cc.TwigStructuralBiomass[j]*f*WoodC[j];
      cc.TwigLivingStructuralBiomass[j] = cc.TwigLivingStructuralBiomass[j]*f*WoodC[j];
      cc.DeadLeafStructuralBiomass[j] = cc.DeadLeafStructuralBiomass[j]*f*leafCperDry;
      cc.DeadTwigStructuralBiomass[j] = cc.DeadTwigStructuralBiomass[j]*f*WoodC[j];
      cc.SapwoodStructuralBiomass[j] = cc.SapwoodStructuralBiomass[j]*f*WoodC[j];
      cc.SapwoodLivingStructuralBiomass[j] = cc.SapwoodLivingStructuralBiomass[j]*f*WoodC[j];
      cc.HeartwoodStructuralBiomass[j] = cc.HeartwoodStructuralBiomass[j]*f*WoodC[j];
      cc.AbovegroundWoodBiomass[j] = cc.AbovegroundWoodBiomass[j]*f*WoodC[j];
      cc.BelowgroundWoodBiomass[j] = cc.BelowgroundWoodBiomass[j]*f*WoodC[j];
      cc.FineRootBiomass[j] = cc.FineRootBiomass[j]*f*rootCperDry;
      cc.LabileBiomass[j] = cc.LabileBiomass[j]*f*(6.0*carbonMolarMass/glucoseMolarMass);
    }
    
    cc.StructuralBiomass[j] = cc.LeafStructuralBiomass[j]  + cc.TwigStructuralBiomass[j]   + cc.SapwoodStructuralBiomass[j]  + cc.FineRootBiomass[j];
    cc.DeadBiomass[j] = cc.DeadLeafStructuralBiomass[j] + cc.DeadTwigStructuralBiomass[j] + cc.HeartwoodStructuralBiomass[j];
    cc.TotalLivingBiomass[j] = cc.LeafStructuralBiomass[j] + cc.TwigLivingStructuralBiomass[j] + cc.SapwoodLivingStructuralBiomass[j] + cc.FineRootBiomass[j] + cc.LabileBiomass[j];
    cc.TotalBiomass[j] = cc.StructuralBiomass[j] + cc.DeadBiomass[j] + cc.LabileBiomass[j];
  }
  // if(biomassUnits=="gC_m2"){
  //   Rcout << " B leaves " << sum(leafStructBiomass) << " fineroots "<< sum(fineRootBiomass) <<" sapwood " << sum(sapwoodStructBiomass) << " heartwood "<< sum(heartwoodStructBiomass)<<" struct " << sum(structuralBiomass)<< " dead " << sum(deadBiomass) << " labile " << sum(labileBiomass) << " total " << sum(totalBiomass) <<"\n";
  // }
}

Rcpp::DataFrame copyCarbonCompartments_c(const CarbonCompartments& cc, ModelInput& x) {
  Rcpp::DataFrame outDF = Rcpp::DataFrame::create(
    Rcpp::Named("LeafStorageVolume") = Rcpp::wrap(cc.LeafStorageVolume),
    Rcpp::Named("SapwoodStorageVolume") = Rcpp::wrap(cc.SapwoodStorageVolume),
    Rcpp::Named("LeafStarchMaximumConcentration") = Rcpp::wrap(cc.LeafStarchMaximumConcentration),
    Rcpp::Named("SapwoodStarchMaximumConcentration") = Rcpp::wrap(cc.SapwoodStarchMaximumConcentration),
    Rcpp::Named("LeafStarchCapacity") = Rcpp::wrap(cc.LeafStarchCapacity),
    Rcpp::Named("SapwoodStarchCapacity") = Rcpp::wrap(cc.SapwoodStarchCapacity),
    Rcpp::Named("LeafStructuralBiomass") = Rcpp::wrap(cc.LeafStructuralBiomass),
    Rcpp::Named("TwigStructuralBiomass") = Rcpp::wrap(cc.TwigStructuralBiomass),
    Rcpp::Named("TwigLivingStructuralBiomass") = Rcpp::wrap(cc.TwigLivingStructuralBiomass),
    Rcpp::Named("DeadLeafStructuralBiomass") = Rcpp::wrap(cc.DeadLeafStructuralBiomass),
    Rcpp::Named("DeadTwigStructuralBiomass") = Rcpp::wrap(cc.DeadTwigStructuralBiomass),
    Rcpp::Named("SapwoodStructuralBiomass") = Rcpp::wrap(cc.SapwoodStructuralBiomass),
    Rcpp::Named("SapwoodLivingStructuralBiomass") = Rcpp::wrap(cc.SapwoodLivingStructuralBiomass),
    Rcpp::Named("HeartwoodStructuralBiomass") = Rcpp::wrap(cc.HeartwoodStructuralBiomass),
    Rcpp::Named("AbovegroundWoodBiomass") = Rcpp::wrap(cc.AbovegroundWoodBiomass),
    Rcpp::Named("BelowgroundWoodBiomass") = Rcpp::wrap(cc.BelowgroundWoodBiomass),
    Rcpp::Named("FineRootBiomass") = Rcpp::wrap(cc.FineRootBiomass),
    Rcpp::Named("StructuralBiomass") = Rcpp::wrap(cc.StructuralBiomass),
    Rcpp::Named("DeadBiomass") = Rcpp::wrap(cc.DeadBiomass),
    Rcpp::Named("LabileBiomass") = Rcpp::wrap(cc.LabileBiomass),
    Rcpp::Named("TotalLivingBiomass") = Rcpp::wrap(cc.TotalLivingBiomass),
    Rcpp::Named("TotalBiomass") = Rcpp::wrap(cc.TotalBiomass));
  outDF.attr("row.names") = x.cohorts.CohortCode;
  return(outDF);
}