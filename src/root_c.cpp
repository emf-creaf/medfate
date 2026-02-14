#include <RcppArmadillo.h>
#include <vector>
#include "root_c.h"

void ldrRS_one_c(std::vector<double>& ldr, 
                 double Z50, double Z95, double Z100, const std::vector<double>&  d){
  int nlayer = d.size();
  double* dCum = new double[nlayer];
  double c = 2.94/log(Z50/Z95);
  ldr[0] = 1.0/(1.0+pow(d[0]/Z50,c));
  //Cumulate d
  dCum[0] = d[0];
  for(int i=1;i<nlayer;i++) dCum[i] = d[i]+dCum[i-1];
  //LDR equation
  for(int i=1;i<nlayer;i++){
    ldr[i] = 1.0/(1.0+pow(dCum[i]/Z50,c)) -1.0/(1.0+pow(dCum[i-1]/Z50,c));
  }
  //Truncate distribution if cumulative depth of the previous layer is larger than maximum rooting depth
  if(!std::isnan(Z100)) {
    for(int i=1;i<nlayer;i++){
      if(dCum[i-1]>Z100) ldr[i] = 0.0;
    }
  }
  //Rescale proportions so that they sum 1
  double Vtot = std::accumulate(ldr.begin(), ldr.end(), 0.0);
  for(int i=0;i<nlayer; i++) {
    ldr[i] = ldr[i]/Vtot;
  }
  delete[] dCum;
}


/**
 * Fine root radius in cm
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootRadius")]]
double fineRootRadius_c(double specificRootLength, double rootTissueDensity) {
  return(sqrt(1.0/(M_PI*specificRootLength*rootTissueDensity)));
}

/**
 *  specificRootSurfaceArea (SRSA; cm2/g) as function of: 
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root tissue density (RTD; g/cm3) e.g. 0.165 g/cm3
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_specificRootSurfaceArea")]]
double specificRootSurfaceArea_c(double specificRootLength, double rootTissueDensity) {
  return(2.0*sqrt(M_PI*specificRootLength/rootTissueDensity));
}

/**
 *   Estimates soil volume (m3) occupied with fine roots
 *    . fine root biomass (g dry)
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root length density (RLD; cm/cm3) e.g. 10 cm/cm3 = 0.1 mm/mm3
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootSoilVolume")]]
double fineRootSoilVolume_c(double fineRootBiomass, double specificRootLength, double rootLengthDensity) {
  return(fineRootBiomass*(specificRootLength/rootLengthDensity)*1e-6);
}


double frv_c(double vol, double B, const std::vector<double>& v, const std::vector<double>& ax, const std::vector<double>& ra) {
  int numLayers = ax.size();
  double s = 0.0;
  double li = 0.0;
  for(int i=0;i<numLayers;i++) {
    li = ax[i]+sqrt(vol)*ra[i];
    s +=(v[i]/li); //No taper effect
    // s +=(V[i]/(li*taperFactorSavage(li*100.0))); //TODO: Improve usage of Savage taper factor for roots
  }
  return(B*s - 1.0);
}


/**
 *   Estimates soil volume (m3) occupied with coarse roots
 *    . sapwood area (cm2)
 *    . rooting depth (cm)
 */
double coarseRootSoilVolumeFromConductance_c(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                                             const std::vector<double>& v, const std::vector<double>& d, const std::vector<double>& rfc) {
  int numLayers = v.size();
  std::vector<double> ra(numLayers, 0.0);
  std::vector<double> ax(numLayers, 0.0);
  for(int j=0;j<numLayers;j++) {
    ra[j] = sqrt(v[j]/((d[j]/1000.0)*M_PI*(1.0 - (rfc[j]/100.0))));
    if(j==0) ax[j] = (d[j]/1000.0);
    else ax[j] = ax[j-1]+(d[j]/1000.0);
    // Rcout<<j<<" "<<ax[j]<<" "<<ra[j]<<"\n";
  }
  double B = (1000.0/0.018)*Kmax_rootxylem/(VCroot_kmax*Al2As);
  // Rcout<<" B: " << B<<"\n";
  double step = 1.0;
  double fTol = 0.005;
  double vol = 0.0;
  double f = frv_c(vol, B, v, ax, ra);
  int nsteps = 0;
  int maxnsteps = 200;
  while((std::abs(f)>fTol) && (nsteps < maxnsteps)) {
    // Rcout<<vol<<"\n";
    if((f > 0.0)) {
      vol += step; 
    } else {
      vol -= step;
      step = step/2.0;
    }
    f = frv_c(vol,B,v, ax,ra);
    nsteps++;
  }
  if(nsteps==maxnsteps) throw medfate::MedfateInternalError("Maximum number of steps reached in coarse root volume estimation");
  // for(int j=0;j<numLayers;j++) {
  // Rcout<<j<<" "<<ax[j]<<" "<<sqrt(vol)*ra[j]<<" "<<((d[j]/1000.0)*M_PI*pow(sqrt(vol)*ra[j],2.0))<<"\n";
  // }
  return(std::max(0.25,vol));
}


/**
 *  Root lengths
 * 
 * Calculates the sum of radial and vertical root lengths.
 * 
 * Sperry, J. S., Y. Wang, B. T. Wolfe, D. S. Mackay, W. R. L. Anderegg, N. G. Mcdowell, and W. T. Pockman. 2016. 
 * Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. 
 * New Phytologist 212:577–589.
 * 
 * Returs: coarse root length in mm (same units as d)
 * 
 */
std::vector<double> coarseRootLengthsFromVolume_c(double VolInd, const std::vector<double>& v, const std::vector<double>& d, const std::vector<double>& rfc) {
  int nlayers = v.size();
  std::vector<double> rl(nlayers), vl(nlayers), tl(nlayers);
  for(int j=0;j<nlayers;j++) {
    if(j==0) vl[j] = d[j];
    else vl[j] = vl[j-1]+d[j];
    rl[j] = 1000.0*sqrt((VolInd*v[j])/((d[j]/1000.0)*M_PI*(1.0 - (rfc[j]/100.0))));
    // Rcout<<vl[j]<<" "<< rl[j]<<"\n";
    tl[j] = vl[j] + rl[j];
  }
  return(tl);
}



RadialAxialLengths coarseRootRadialAxialLengths_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio) {
  int nlayers = v.size();
  RadialAxialLengths radax(nlayers);
  double maxRootDepth = 0.0;
  
  //Vertical lengths
  std::vector<double> zini(nlayers);
  for(int i=0;i<nlayers;i++) {
    if(i==0) {
      zini[i] = 0.0;
    } else {
      zini[i] = zini[i-1]+ d[i-1];
    }
    if(v[i]>0.0) {
      radax.axial[i] = zini[i]+ d[i]/2.0;
      maxRootDepth +=d[i];
    } else {
      radax.axial[i] = 0.0;
    }
    // Rcout<<vl[i]<<" ";
  }
  // Rcout<<"\n";
  int nlayerseff = nlayers;
  for(int i=(nlayers-1);i>=0;i--) if(radax.axial[i]>0.0) nlayerseff = i;
  
  //Radial lengths
  std::vector<double> r(nlayers, 0.0);
  double maxr = 0.0;
  for(int i=0;i<nlayerseff;i++) {
    r[i] = sqrt(v[i]/(d[i]*M_PI));
    maxr = std::max(r[i],maxr); 
  }
  // Rcout<<maxr<<"\n";
  for(int i=0;i<nlayerseff;i++) {
    radax.radial[i] = maxRootDepth*depthWidthRatio*(r[i]/maxr);
    // Rcout<<rl[i]<<" ";
  }
  return(radax);
}

std::vector<double> coarseRootLengths_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio) {
  int nlayers = v.size();
  RadialAxialLengths radax = coarseRootRadialAxialLengths_c(v, d, depthWidthRatio);
  std::vector<double> l(nlayers, 0.0);
  for(int i=0;i<nlayers;i++) {
    l[i]= (radax.radial[i]+ radax.axial[i]);
  }
  return(l);
}


double coarseRootSoilVolume_c(const std::vector<double>& v, const std::vector<double>& d, double depthWidthRatio) {
  int nlayers = v.size();
  RadialAxialLengths radax = coarseRootRadialAxialLengths_c(v, d, depthWidthRatio);
  //Weights
  double volInd = 0.0;
  for(int i=0;i<nlayers;i++) {
    volInd += 1e-9*(std::pow(radax.radial[i],2.0)*M_PI)*d[i];
  }
  return(volInd);
}


//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootHalfDistance")]]
double fineRootHalfDistance_c(double rootLengthDensity) {
  return(1.0/sqrt(M_PI*rootLengthDensity));
}

/**
 * Derivation of fine root length per ground area (m·m-2) 
 * from the conductance factor for cylindrical flow geometry of the rhizosphere
 * 
 */
double fineRootLengthPerArea_c(double Ksoil, double krhizo, double lai,
                               double radius, double rootLengthDensity) {
  double Xi = krhizo*lai/Ksoil;
  double rmax = fineRootHalfDistance_c(rootLengthDensity);
  return(log(pow(rmax,2.0)/pow(radius,2.0))*Xi/(4.0*M_PI));
}

double fineRootBiomassPerIndividual_c(const std::vector<double>& Ksoil, const std::vector<double>& krhizo, double lai, double N,
                                      double specificRootLength, double rootTissueDensity,  
                                      double rootLengthDensity) {
  double r = fineRootRadius_c(specificRootLength, rootTissueDensity); //cm
  int numLayers = Ksoil.size();
  double frb = 0.0;
  for(int l=0;l<numLayers;l++) {
    double LA = fineRootLengthPerArea_c(Ksoil[l], krhizo[l], lai, r, rootLengthDensity);//m·m-2 
    // Rcout<<l<<" "<<LA<<"\n";
    frb += (10000.0*LA)/(N*0.01*specificRootLength);
  }
  return(frb);
}

/**
 * Derivation of maximum rhizosphere conductance (mmol·m-2·s·MPa-1)
 * from the conductance factor for cylindrical flow geometry of the rhizosphere
 * 
 */
double fineRootMaximumConductance_c(double Ksoil, double fineRootLengthPerArea, double lai,
                                  double radius, double rootLengthDensity) {
  double rmax = fineRootHalfDistance_c(rootLengthDensity);
  return((fineRootLengthPerArea*Ksoil*4.0*M_PI)/(lai*log(pow(rmax,2.0)/pow(radius,2.0))));
}

//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_rhizosphereMaximumConductance")]]
std::vector<double> rhizosphereMaximumConductance_c(const std::vector<double>& Ksoil, const std::vector<double>& fineRootBiomass, double lai, double N,
                                                    double specificRootLength, double rootTissueDensity,  
                                                    double rootLengthDensity) {
   double r = fineRootRadius_c(specificRootLength, rootTissueDensity); //cm
   int numLayers = Ksoil.size();
   std::vector<double> krhizo(numLayers, 0.0);
   for(int l=0;l<numLayers;l++) {
     double FRLA = 1e-6*fineRootBiomass[l]*N*specificRootLength;
     krhizo[l] = fineRootMaximumConductance_c(Ksoil[l], FRLA, lai, r, rootLengthDensity);
   }
   return(krhizo);
 }
