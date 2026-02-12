#include <RcppArmadillo.h>
#include "tissuemoisture_c.h"
#include "carbon_c.h"
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
    Rcpp::Named("SapwoodStructuralBiomass") = Rcpp::wrap(cc.SapwoodStructuralBiomass),
    Rcpp::Named("SapwoodLivingStructuralBiomass") = Rcpp::wrap(cc.SapwoodLivingStructuralBiomass),
    Rcpp::Named("HeartwoodStructuralBiomass") = Rcpp::wrap(cc.HeartwoodStructuralBiomass),
    Rcpp::Named("AbovegroundWoodBiomass") = Rcpp::wrap(cc.AbovegroundWoodBiomass),
    Rcpp::Named("BelowgroundWoodBiomass") = Rcpp::wrap(cc.BelowgroundWoodBiomass),
    Rcpp::Named("FineRootBiomass") = Rcpp::wrap(cc.FineRootBiomass),
    Rcpp::Named("StructuralBiomass") = Rcpp::wrap(cc.StructuralBiomass),
    Rcpp::Named("LabileBiomass") = Rcpp::wrap(cc.LabileBiomass),
    Rcpp::Named("TotalLivingBiomass") = Rcpp::wrap(cc.TotalLivingBiomass),
    Rcpp::Named("TotalBiomass") = Rcpp::wrap(cc.TotalBiomass));
  outDF.attr("row.names") = x.cohorts.CohortCode;
  return(outDF);
}