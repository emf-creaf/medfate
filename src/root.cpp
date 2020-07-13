#include <Rcpp.h>
#include <numeric>
#include "hydraulics.h"
using namespace Rcpp;
using namespace std;
using std::exp;
using std::log;
using std::sqrt;



/**
 *  Root distribution
 */
NumericVector ldrRS_one(double Z50, double Z95, NumericVector d){
  int nlayer = d.size();
  NumericVector dCum = clone(d);
  NumericVector Vd(nlayer);
  double c = 2.94/log(Z50/Z95);
  double Vtot = 0.0;
  Vd[0] = 1.0/(1.0+pow(d[0]/Z50,c));
  Vtot = Vd[0];
  for(int i=1;i<nlayer;i++) dCum[i] = dCum[i]+dCum[i-1];
  for(int i=1;i<nlayer;i++){
    Vd[i] = 1.0/(1.0+pow(dCum[i]/Z50,c)) -1.0/(1.0+pow(dCum[i-1]/Z50,c));
    Vtot +=Vd[i];
  }
  //Rescale proportions so that they sum 1
  for(int i=0;i<nlayer; i++) {
    Vd[i] = Vd[i]/Vtot;
  }
  return(Vd);
}
NumericVector conicRS_one(double Zcone, NumericVector d){
  int nlayers  = d.size();
  NumericVector Zd(nlayers,0.0);
  NumericVector Vd(nlayers,0.0);
  
  //Proportion of root length on each layer
  double sumZd=0.0;
  for(int l=0;l<nlayers;l++) {
     if(l==0) Zd[l] = std::min(d[l],Zcone)/Zcone;
     else if(l==(nlayers-1)) Zd[l] = 1.0-sumZd;
     else {
       Zd[l] = std::max(std::min(d[l]/Zcone,1.0-sumZd),0.0);  
     }
     sumZd = sumZd+Zd[l];
  }

  //Proportion of fine roots in each layer (assuming soil depth = maximum root depth)
  for(int l1=0;l1<nlayers;l1++) {
    if(l1==0) {// Vd[0] = 1.0-pow(Zd[1]+Zd[2],3.0);
      sumZd = 0.0;
      for(int l2=1;l2<nlayers;l2++) sumZd += Zd[l2];
      // Rcout<<"\n "<<l1<<" S "<<sumZd<<"\n";
      Vd[l1] = 1.0 - pow(sumZd,3.0);
    } else if(l1==(nlayers-1)) {
      // Vd[2] = pow(1.0-Zd[0]-Zd[1],3.0);
      sumZd = 0.0;
      for(int l2=0;l2<l1;l2++) {
        sumZd += Zd[l2]; 
        // Rcout<<"\n "<<l2<<" S "<<sumZd<<"\n";
      }
      Vd[l1] = pow(1.0-sumZd, 3.0);
    } else {
      // Vd[1] = pow(1.0-Zd[0],3.0)-pow(Zd[2],3.0);
      double sumZd1 = 0.0, sumZd2 = 0.0;
      for(int l2=0;l2<l1;l2++) sumZd1 += Zd[l2];
      for(int l2=l1+1;l2<nlayers;l2++) sumZd2 += Zd[l2];
      Vd[l1] = pow(1.0-sumZd1,3.0)-pow(sumZd2,3.0);
    }
    // Rcout<<"("<<Vd[l1]<<", "<<Zd[l1]<<")\n";
  }
  // Rcout<<"\n";
  //Remove non-existing part
  double Zsoil = 0.0;
  for(int l=0;l<nlayers;l++) Zsoil +=d[l];
  double Zreal = std::max((Zcone-Zsoil)/Zcone,0.0);  
  double Vrem = pow(Zreal,3.0);
  double Vtot = 0.0;
  for(int i=(nlayers-1);i>=0; i--) {
    Vd[i] = Vd[i]-Vrem;
    Vrem = (-1)*std::min(Vd[i],0.0);
    Vd[i] = std::max(Vd[i],0.0);
    Vtot +=Vd[i];
  }
  
  //Rescale proportions so that they sum 1
  for(int l=0;l<nlayers; l++) {
    Vd[l] = Vd[l]/Vtot;
  }
  return(Vd);
}

// [[Rcpp::export("root_conicDistribution")]]
NumericMatrix conicDistribution(NumericVector Zcone, NumericVector d) {
  int numCohorts = Zcone.size();
  NumericMatrix P(numCohorts,d.size());
  for(int c=0;c<numCohorts;c++){
    NumericVector PC = conicRS_one(Zcone[c],d);
    for(int i=0;i<PC.size();i++) P(c,i) = PC[i];
  }
  return(P);
}

/**
 * Root distribution according to the LDR model
 * 
 * Schenk, H., Jackson, R., 2002. The global biogeography of roots. Ecol. Monogr. 72, 311–328.
 */
// [[Rcpp::export("root_ldrDistribution")]]
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector d) {
  int numCohorts = Z50.size();
  NumericMatrix P(numCohorts,d.size());
  NumericVector PC;
  for(int c=0;c<numCohorts;c++){
    PC = ldrRS_one(Z50[c], Z95[c],d);
    for(int i=0;i<d.size();i++) P(c,i) = PC[i];
  }
  return(P);
}
NumericMatrix ldrDistribution(NumericVector treeZ50, NumericVector shrubZ50, 
                              NumericVector treeZ95, NumericVector shrubZ95, NumericVector d) {
  int ntree = treeZ50.size();
  int nshrub = shrubZ50.size();
  int nlayers = d.size();
  NumericMatrix V(ntree+nshrub,nlayers);
  for(int i=0;i<ntree;i++) {
    V(i,_) = ldrRS_one(treeZ50[i], treeZ95[i],d);
  }
  for(int i=0;i<nshrub;i++) {
    V(ntree+i,_) = ldrRS_one(shrubZ50[i],shrubZ95[i],d);
  }
  return(V);
}

// [[Rcpp::export(".rootDistribution")]]
NumericMatrix rootDistribution(NumericVector z, List x) {
  DataFrame treeData = Rcpp::as<Rcpp::DataFrame>(x["treeData"]);
  DataFrame shrubData = Rcpp::as<Rcpp::DataFrame>(x["shrubData"]);
  int ntree = treeData.nrows();
  int nshrub = shrubData.nrows();
  
  NumericVector treeZ50 = Rcpp::as<Rcpp::NumericVector>(treeData["Z50"]);
  NumericVector treeZ95 = Rcpp::as<Rcpp::NumericVector>(treeData["Z95"]);
  NumericVector shrubZ50 =  Rcpp::as<Rcpp::NumericVector>( shrubData["Z50"]);  
  NumericVector shrubZ95 =  Rcpp::as<Rcpp::NumericVector>( shrubData["Z95"]);  
  NumericMatrix rdtree = ldrDistribution(treeZ50, treeZ95, z);
  NumericMatrix rdshrub = ldrDistribution(shrubZ50, shrubZ95, z);
  NumericMatrix rd = NumericMatrix(ntree+nshrub, z.length());
  for(int i=0;i<ntree;i++) {
    rd(i,_) = rdtree(i,_);
  }
  for(int i=0;i<nshrub;i++) {
    rd(ntree+i,_) = rdshrub(i,_);
  }
  return(rd);
}

// DataFrame rootSpatialDimensions(double rootVolumeIndividual, NumericVector v, NumericVector d, NumericVector bulkDensity) {
//   int numLayers = v.size();
//   NumericVector lvol(numLayers,0.0);
//   NumericVector larea(numLayers,0.0);
//   NumericVector laxial(numLayers,0.0);
//   NumericVector lradial(numLayers,0.0);
//   for(int i=0;i<numLayers;i++) {
//     lvol[i] = rootVolumeIndividual*v[i]; //m3
//     laxial[i] = d[i]/1000.0; //mm to m
//     lradial[i] = sqrt(lvol[i]/(laxial[i]*PI*(1.0 - (bulkDensity[i]/2.65))));
//     larea[i] = PI*pow(lradial[i],2.0);
//   }
//   DataFrame df = DataFrame::create(_["Volume"] = lvol,
//                                    _["GroundArea"] = larea,
//                                    _["AxialLength"]  = laxial,
//                                    _["RadialLength"] = lradial);
//   return(df);
// }

/**
 * Calculates root ground area for each layer of the root system of 
 * an individual with a given total fine root volume (in m3) and fine root proportions
 */
// [[Rcpp::export("root_individualRootedGroundArea")]]
NumericMatrix individualRootedGroundArea(NumericVector VolInd, NumericMatrix V, NumericVector d, 
                                         NumericVector rfc) {
  int numCohorts = V.nrow();
  int numLayers = V.ncol();
  NumericMatrix larea(numCohorts,numLayers);
  for(int i=0;i<numCohorts;i++) {
    for(int j=0;j<numLayers;j++) {
      double lvol = VolInd[i]*V(i,j); //m3
      double laxial = d[j]/1000.0; //mm to m
      double lradial = sqrt(lvol/(laxial*PI*(1.0 - (rfc[j]/100.0))));
      larea(i,j) = PI*pow(lradial,2.0); //m2
    }
  }
  larea.attr("dimnames") = V.attr("dimnames");
  return(larea);
}



/**
 *  specificRootSurfaceArea (SRSA; cm2/g) as function of: 
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root tissue density (RTD; g/cm3) e.g. 0.165 g/cm3
 */
// [[Rcpp::export("root_specificRootSurfaceArea")]]
double specificRootSurfaceArea(double specificRootLength, double rootTissueDensity) {
  return(2.0*sqrt(PI*specificRootLength/rootTissueDensity));
}
/**
 * Fine root radius in cm
 */
// [[Rcpp::export("root_fineRootRadius")]]
double fineRootRadius(double specificRootLength, double rootTissueDensity) {
  return(sqrt(1.0/(PI*specificRootLength*rootTissueDensity)));
}
// [[Rcpp::export("root_fineRootHalfDistance")]]
double fineRootHalfDistance(double rootLengthDensity) {
  return(1.0/sqrt(PI*rootLengthDensity));
}
/**
 * Derivation of fine root length per ground area (m·m-2) 
 * from the conductance factor for cylindrical flow geometry of the rhizosphere
 * 
 */
double fineRootLengthPerArea(double Ksoil, double krhizo, double lai,
                      double radius, double rootLengthDensity) {
  double Xi = krhizo*lai/Ksoil;
  double rmax = fineRootHalfDistance(rootLengthDensity);
  return(log(pow(rmax,2.0)/pow(radius,2.0))*Xi/(4.0*PI));
}
/**
 * Derivation of maximum rhizosphere conductance (mmol·m-2·s·MPa-1)
 * from the conductance factor for cylindrical flow geometry of the rhizosphere
 * 
 */
double fineRootMaximumConductance(double Ksoil, double fineRootLengthPerArea, double lai,
                                   double radius, double rootLengthDensity) {
  double rmax = fineRootHalfDistance(rootLengthDensity);
  return((fineRootLengthPerArea*Ksoil*4.0*PI)/(lai*log(pow(rmax,2.0)/pow(radius,2.0))));
}

// [[Rcpp::export("root_fineRootAreaIndex")]]
double fineRootAreaIndex(NumericVector Ksoil, NumericVector krhizo, double lai,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  double r = fineRootRadius(specificRootLength, rootTissueDensity); //cm
  int numLayers = Ksoil.size();
  double frai = 0.0;
  for(int l=0;l<numLayers;l++) {
    double LA = fineRootLengthPerArea(Ksoil[l], krhizo[l], lai, r, rootLengthDensity);//m·m-2 
    // Rcout<<l<<" "<<LA<<"\n";
    frai += LA*2.0*PI*(r/100.0);
  }
  return(frai);
}
/**
 * Fine root biomass in g dry · ind-1
 */
// [[Rcpp::export("root_fineRootBiomass")]]
double fineRootBiomassPerIndividual(NumericVector Ksoil, NumericVector krhizo, double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  double r = fineRootRadius(specificRootLength, rootTissueDensity); //cm
  int numLayers = Ksoil.size();
  double frb = 0.0;
  for(int l=0;l<numLayers;l++) {
    double LA = fineRootLengthPerArea(Ksoil[l], krhizo[l], lai, r, rootLengthDensity);//m·m-2 
    // Rcout<<l<<" "<<LA<<"\n";
    frb += (10000.0*LA)/(N*0.01*specificRootLength);
  }
  return(frb);
}

// [[Rcpp::export("root_rhizosphereMaximumConductance")]]
NumericVector rhizosphereMaximumConductance(NumericVector Ksoil, NumericVector fineRootBiomass, double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  double r = fineRootRadius(specificRootLength, rootTissueDensity); //cm
  int numLayers = Ksoil.size();
  NumericVector krhizo(numLayers, 0.0);
  for(int l=0;l<numLayers;l++) {
    double FRLA = 1e-6*fineRootBiomass[l]*N*specificRootLength;
    krhizo[l] = fineRootMaximumConductance(Ksoil[l], FRLA, lai, r, rootLengthDensity);
  }
  return(krhizo);
}


/**
 *   Estimates soil volume (m3) occupied with fine roots
 *    . fine root biomass (g dry)
 *    . specific root length (SRL; cm/g) e.g. 3870 cm/g
 *    . root length density (RLD; cm/cm3) e.g. 10 cm/cm3 = 0.1 mm/mm3
 */
// [[Rcpp::export("root_fineRootSoilVolume")]]
double fineRootSoilVolume(double fineRootBiomass, double specificRootLength, double rootLengthDensity) {
  return(fineRootBiomass*(specificRootLength/rootLengthDensity)*1e-6);
}

double frv(double vol, double B, NumericVector v, NumericVector ax, NumericVector ra) {
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
// [[Rcpp::export("root_coarseRootSoilVolumeFromConductance")]]
double coarseRootSoilVolumeFromConductance(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                                   NumericVector v, NumericVector d, NumericVector rfc) {
  int numLayers = v.size();
  NumericVector ra(numLayers, 0.0);
  NumericVector ax(numLayers, 0.0);
  for(int j=0;j<numLayers;j++) {
    ra[j] = sqrt(v[j]/((d[j]/1000.0)*PI*(1.0 - (rfc[j]/100.0))));
    if(j==0) ax[j] = (d[j]/1000.0);
    else ax[j] = ax[j-1]+(d[j]/1000.0);
    // Rcout<<j<<" "<<ax[j]<<" "<<ra[j]<<"\n";
  }
  double B = (1000.0/0.018)*Kmax_rootxylem/(VCroot_kmax*Al2As);
  // Rcout<<" B: " << B<<"\n";
  double step = 1.0;
  double fTol = 0.005;
  double vol = 0.0;
  double f = frv(vol, B, v, ax, ra);
  while(std::abs(f)>fTol) {
    if(f>0.0) {
      vol += step; 
    } else {
      vol -= step;
      step = step/2.0;
    }
    f = frv(vol,B,v, ax,ra);
  }
  // for(int j=0;j<numLayers;j++) {
    // Rcout<<j<<" "<<ax[j]<<" "<<sqrt(vol)*ra[j]<<" "<<((d[j]/1000.0)*PI*pow(sqrt(vol)*ra[j],2.0))<<"\n";
  // }
  return(vol);
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
// [[Rcpp::export("root_coarseRootLengthsFromVolume")]]
NumericVector coarseRootLengthsFromVolume(double VolInd, NumericVector v, NumericVector d, NumericVector rfc) {
  int nlayers = v.size();
  NumericVector rl(nlayers), vl(nlayers), tl(nlayers);
  for(int j=0;j<nlayers;j++) {
    if(j==0) vl[j] = d[j];
    else vl[j] = vl[j-1]+d[j];
    rl[j] = 1000.0*sqrt((VolInd*v[j])/((d[j]/1000.0)*PI*(1.0 - (rfc[j]/100.0))));
    // Rcout<<vl[j]<<" "<< rl[j]<<"\n";
    tl[j] = vl[j] + rl[j];
  }
  return(tl);
}

List coarseRootRadialAxialLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
  int nlayers = v.size();
  double maxRootDepth = 0.0;
  //Vertical lengths
  NumericVector vl(nlayers), zini(nlayers);
  for(int i=0;i<nlayers;i++) {
    if(i==0) {
      zini[i] = 0.0;
    } else {
      zini[i] = zini[i-1]+ d[i-1];
    }
    if(v[i]>0.0) {
      vl[i] = zini[i]+ d[i]/2.0;
      maxRootDepth +=d[i];
    } else {
      vl[i] = NA_REAL;
    }
    // Rcout<<vl[i]<<" ";
  }
  // Rcout<<"\n";
  int nlayerseff = nlayers;
  for(int i=(nlayers-1);i>=0;i--) if(NumericVector::is_na(vl[i])) nlayerseff = i;
  
  //Radial lengths
  NumericVector r(nlayerseff), rl(nlayerseff);
  double maxr = 0.0;
  for(int i=0;i<nlayerseff;i++) {
    r[i] = sqrt(v[i]/(d[i]*PI));
    maxr = std::max(r[i],maxr); 
  }
  // Rcout<<maxr<<"\n";
  for(int i=0;i<nlayerseff;i++) {
    rl[i] = maxRootDepth*depthWidthRatio*(r[i]/maxr);
    // Rcout<<rl[i]<<" ";
  }
  return(List::create(_["radial"] = rl, _["axial"] = vl));
}

// [[Rcpp::export("root_coarseRootLengths")]]
NumericVector coarseRootLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
  List radax = coarseRootRadialAxialLengths(v, d, depthWidthRatio);
  NumericVector rl = radax["radial"];
  NumericVector vl = radax["axial"];
  int nlayers = rl.size();
  NumericVector l(nlayers, 0.0);
  for(int i=0;i<nlayers;i++) {
    l[i]= (rl[i]+vl[i]);
  }
  return(l);
}

// [[Rcpp::export("root_coarseRootSoilVolume")]]
double coarseRootSoilVolume(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
  List radax = coarseRootRadialAxialLengths(v, d, depthWidthRatio);
  NumericVector rl = radax["radial"];
  int nlayers = rl.size();
  //Weights
  double volInd = 0.0;
  for(int i=0;i<nlayers;i++) {
    volInd += 1e-9*(std::pow(rl[i],2.0)*PI)*d[i];
  }
  return(volInd);
}

// List horizontalProportionsBasic(NumericVector poolProportions, NumericMatrix V, 
//                                 double LAIcell, double poolOverlapFactor) {
//   int numCohorts = V.nrow();
//   int numlayers = V.ncol();
//   double ropmax = (1.0 - exp(-(poolOverlapFactor*LAIcell)));
//   List l(numCohorts);
//   for(int coh=0;coh<numCohorts;coh++) {
//     NumericMatrix RHOP(numCohorts,numlayers);
//     double Vmax = 0.0;
//     for(int l=0;l<numlayers;l++) {
//       Vmax = std::max(Vmax, V(coh,l));
//     }
//     for(int l=0;l<numlayers;l++) {
//       double s = 0.0;
//       for(int c=0;c<numCohorts;c++) {
//         if(c!=coh) {
//           RHOP(c,l) = poolProportions[c]*ropmax*sqrt(V(coh,l)*V(c,l))/Vmax;
//           s +=RHOP(c,l);
//         } 
//       }
//       RHOP(coh,l) = 1.0 - s;
//     }
//     RHOP.attr("dimnames") = V.attr("dimnames");
//     l[coh] = RHOP;
//   }
//   l.attr("names") = rownames(V);
//   return(l);
// }

// [[Rcpp::export("root_horizontalProportions")]]
List horizontalProportions(NumericVector poolProportions, NumericVector VolInd, NumericVector N, NumericMatrix V, 
                           NumericVector d, NumericVector rfc) {
  
  int numCohorts = V.nrow();
  int numlayers = V.ncol();
  List l(numCohorts);
  NumericVector poolAreaInd(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    poolAreaInd[c] = 10000.0*poolProportions[c]/N[c]; //area of the pool per individual of the cohort
  }
  
  NumericMatrix iga = individualRootedGroundArea(VolInd,V,d,rfc);

  for(int coh=0;coh<numCohorts;coh++) {
    NumericMatrix RHOP(numCohorts,numlayers);
    for(int l=0;l<numlayers;l++) {
      // Rcout<<coh<< " "<< l<< " "<<iga(coh,l)<<" "<< poolAreaInd[coh]<<"\n";
      RHOP(coh,l) = std::min(poolAreaInd[coh],iga(coh,l))/iga(coh,l);
      if(iga(coh,l)>poolAreaInd[coh]) {
        double dif = iga(coh,l) - poolAreaInd[coh];
        // Rcout<<dif<<"\n";
        for(int c2=0; c2<numCohorts;c2++) { //divide the remainder among all cohorts (including itself) depending on their proportions
          RHOP(c2,l) += poolProportions[c2]*dif/iga(coh,l);
        }
      }
    }
    RHOP.attr("dimnames") = V.attr("dimnames");
    l[coh] = RHOP;
  }
  l.attr("names") = rownames(V);
  return(l);
}

