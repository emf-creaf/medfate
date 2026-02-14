#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include "hydraulics.h"
#include "root_c.h"
using namespace Rcpp;
using namespace std;
using std::exp;
using std::log;
using std::sqrt;



/**
 *  Root distribution
 */
NumericVector ldrRS_one(double Z50, double Z95, double Z100, NumericVector d){
  std::vector<double> d_c = Rcpp::as<std::vector<double>>(d);
  std::vector<double> lrd_c(d.size());
  ldrRS_one_c(lrd_c, Z50, Z95, Z100, d_c);
  NumericVector Vd = Rcpp::wrap(lrd_c);
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

//' Root functions
//' 
//' Functions to calculate properties of fine/coarse roots within the soil, given root system parameters and soil layer definition.
//' 
//' @param Z50 A vector of depths (in mm) corresponding to 50% of roots.
//' @param Z95 A vector of depths (in mm) corresponding to 95% of roots.
//' @param Z100 A vector of depths (in mm) corresponding to 100% of roots.
//' @param Zcone A vector of depths (in mm) corresponding to the root cone tip.
//' @param d The width (in mm) corresponding to each soil layer.
//' @param v Vector of proportions of fine roots in each soil layer.
//' @param depthWidthRatio Ratio between radius of the soil layer with the largest radius and maximum rooting depth.
//' @param rfc Percentage of rock fragment content (volume basis) for each layer.
//' @param Kmax_rootxylem Sapwood-specific hydraulic conductivity of root xylem (in kg H2O·s-1·m-1·MPa-1).
//' @param VCroot_kmax Root xylem maximum conductance per leaf area (mmol·m-2·s-1·MPa-1). 
//' @param Al2As Leaf area to sapwood area ratio (in m2·m-2).
//' @param specificRootLength Specific fine root length (length of fine roots over weight).
//' @param rootTissueDensity Fine root tissue density (weight over volume at turgidity).
//' @param Ksoil Soil saturated conductivity (mmol·m-1·s-1·MPa-1).
//' @param krhizo Rhizosphere maximum conductance per leaf area (mmol·m-2·s-1·MPa-1).
//' @param lai Leaf area index.
//' @param rootLengthDensity Fine root length density (length of fine roots over soil volume; cm/cm3)
//' @param fineRootBiomass Biomass of fine roots (g).
//' @param V Matrix of proportions of fine roots (cohorts x soil layers).
//' @param VolInd Volume of soil (in m3) occupied by coarse roots per individual. 
//' @param N Density of individuals per hectare.
//' @param poolProportions Division of the stand area among plant cohorts (proportions).
//' 
//' @details
//' \itemize{
//'   \item{\code{root_conicDistribution()} assumes a (vertical) conic distribution of fine roots, whereas \code{root_ldrDistribution()} distributes fine roots according to the linear dose response model of Schenck & Jackson (2002). Return a matrix of fine root proportions in each layer with as many rows as elements in \code{Z} (or \code{Z50}) and as many columns as soil layers.}
//'   \item{\code{root_coarseRootLengths()} and \code{root_coarseRootLengthsFromVolume()} estimate the length of coarse roots (mm) for each soil layer, including axial and radial lengths.}
//'   \item{\code{root_coarseRootSoilVolume} estimates the soil volume (m3) occupied by coarse roots of an individual.}
//'   \item{\code{root_coarseRootSoilVolumeFromConductance} estimates the soil volume (m3) occupied by coarse roots of an individual from root xylem conductance.}
//'   \item{\code{root_fineRootHalfDistance()} calculates the half distance (cm) between neighbouring fine roots.}
//'   \item{\code{root_fineRootRadius()} calculates the radius of fine roots (cm).}
//'   \item{\code{root_fineRootAreaIndex()} estimates the fine root area index for a given soil conductivity and maximum rhizosphere conductance.}
//'   \item{\code{root_fineRootBiomass()} estimates the biomass of fine roots (g dry/individual) for a given soil conductivity and maximum rhizosphere conductance.}
//'   \item{\code{root_rhizosphereMaximumConductance()} is the inverse of the preceeding function, i.e. it estimates rhizosphere conductance from soil conductivity and fine root biomass.}
//'   \item{\code{root_fineRootSoilVolume()} calculates the soil volume (m3) occupied with fine roots.}
//'   \item{\code{root_specificRootSurfaceArea()} returns the specific fine root area (cm2/g).}
//'   \item{\code{root_individualRootedGroundArea()} calculates the area (m2) covered by roots of an individual, for each soil layer.}
//'   \item{\code{root_horizontalProportions()} calculates the (horizontal) proportion of roots of each cohort in the water pool corresponding to itself and that of other cohorts, for each soil layer. Returns a list (with as many elements as cohorts) with each element being a matrix.}
//'   }
//' 
//' @return See details.
//' 
//' @references
//' Schenk, H., Jackson, R., 2002. The global biogeography of roots. Ecol. Monogr. 72, 311–328.
//' 
//' Sperry, J. S., Y. Wang, B. T. Wolfe, D. S. Mackay, W. R. L. Anderegg, N. G. Mcdowell, and W. T. Pockman. 2016. Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. New Phytologist 212, 577–589.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//'  \code{\link{spwb}},  \code{\link{spwbInput}}, \code{\link{soil}}
//'
//' @examples
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' ntree <- nrow(exampleforest$treeData)
//' 
//' #Initialize soil with default soil params
//' s <- defaultSoilParams(4)
//' 
//' #Calculate conic root system for trees
//' V1 <- root_conicDistribution(Z=rep(2000,ntree), s$widths)            
//' print(V1)
//'      
//' #Calculate LDR root system for trees (Schenck & Jackson 2002)
//' V2 <- root_ldrDistribution(Z50 = rep(200,ntree), 
//'                            Z95 = rep(1000,ntree),
//'                            Z100 = rep(NA, ntree), s$widths)
//' print(V2)     
//' 
//' @name root
//' @keywords internal
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
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_ldrDistribution")]]
NumericMatrix ldrDistribution(NumericVector Z50, NumericVector Z95, NumericVector Z100, NumericVector d) {
  int numCohorts = Z50.size();
  NumericMatrix P(numCohorts,d.size());
  NumericVector PC;
  for(int c=0;c<numCohorts;c++){
    PC = ldrRS_one(Z50[c], Z95[c], Z100[c], d);
    for(int i=0;i<d.size();i++) P(c,i) = PC[i];
  }
  return(P);
}
NumericMatrix ldrDistribution(NumericVector treeZ50, NumericVector shrubZ50, NumericVector herbZ50, 
                              NumericVector treeZ95, NumericVector shrubZ95, NumericVector herbZ95, 
                              NumericVector treeZ100, NumericVector shrubZ100, NumericVector herbZ100, 
                              NumericVector d) {
  int ntree = treeZ50.size();
  int nshrub = shrubZ50.size();
  int nherb = herbZ50.size();
  int nlayers = d.size();
  NumericMatrix V(ntree+nshrub+nherb,nlayers);
  for(int i=0;i<ntree;i++) {
    V(i,_) = ldrRS_one(treeZ50[i], treeZ95[i], treeZ100[i], d);
  }
  for(int i=0;i<nshrub;i++) {
    V(ntree+i,_) = ldrRS_one(shrubZ50[i], shrubZ95[i], shrubZ100[i], d);
  }
  for(int i=0;i<nherb;i++) {
    V(ntree+nshrub+i,_) = ldrRS_one(herbZ50[i], herbZ95[i], herbZ100[i], d);
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
  NumericVector treeZ100(ntree, NA_REAL);
  if(treeData.containsElementNamed("Z100")) treeZ100 = Rcpp::as<Rcpp::NumericVector>(treeData["Z100"]);
  NumericVector shrubZ50 =  Rcpp::as<Rcpp::NumericVector>( shrubData["Z50"]);  
  NumericVector shrubZ95 =  Rcpp::as<Rcpp::NumericVector>( shrubData["Z95"]);
  NumericVector shrubZ100(nshrub, NA_REAL);
  if(shrubData.containsElementNamed("Z100")) shrubZ100 = Rcpp::as<Rcpp::NumericVector>(shrubData["Z100"]);

  NumericMatrix rdtree = ldrDistribution(treeZ50, treeZ95, treeZ100, z);
  NumericMatrix rdshrub = ldrDistribution(shrubZ50, shrubZ95, shrubZ100, z);
  
  int nherb = 0;
  NumericMatrix rdherb;
  if(x.containsElementNamed("herbData")) {
    DataFrame herbData = Rcpp::as<Rcpp::DataFrame>(x["herbData"]);
    nherb = herbData.nrows();
    NumericVector herbZ50 = herbData["Z50"];
    NumericVector herbZ95 = herbData["Z95"];
    NumericVector herbZ100(nherb, NA_REAL);
    if(herbData.containsElementNamed("Z100")) herbZ100 = Rcpp::as<Rcpp::NumericVector>(herbData["Z100"]);
    rdherb = ldrDistribution(herbZ50, herbZ95, herbZ100, z);
  }
  
  NumericMatrix rd = NumericMatrix(ntree+nshrub+nherb, z.length());
  for(int i=0;i<ntree;i++) {
    rd(i,_) = rdtree(i,_);
  }
  for(int i=0;i<nshrub;i++) {
    rd(ntree+i,_) = rdshrub(i,_);
  }
  for(int i=0;i<nherb;i++) {
    rd(ntree+nshrub+i,_) = rdherb(i,_);
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
//     lradial[i] = sqrt(lvol[i]/(laxial[i]*M_PI*(1.0 - (bulkDensity[i]/2.65))));
//     larea[i] = M_PI*pow(lradial[i],2.0);
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
//' @rdname root
//' @keywords internal
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
      double lradial = sqrt(lvol/(laxial*M_PI*(1.0 - (rfc[j]/100.0))));
      larea(i,j) = M_PI*pow(lradial,2.0); //m2
    }
  }
  larea.attr("dimnames") = V.attr("dimnames");
  return(larea);
}









//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootAreaIndex")]]
double fineRootAreaIndex(NumericVector Ksoil, NumericVector krhizo, double lai,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  double r = fineRootRadius_c(specificRootLength, rootTissueDensity); //cm
  int numLayers = Ksoil.size();
  double frai = 0.0;
  for(int l=0;l<numLayers;l++) {
    double LA = fineRootLengthPerArea_c(Ksoil[l], krhizo[l], lai, r, rootLengthDensity);//m·m-2 
    // Rcout<<l<<" "<<LA<<"\n";
    frai += LA*2.0*M_PI*(r/100.0);
  }
  return(frai);
}
/**
 * Fine root biomass in g dry · ind-1
 */
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_fineRootBiomass")]]
double fineRootBiomassPerIndividual(NumericVector Ksoil, NumericVector krhizo, double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  return(fineRootBiomassPerIndividual_c(as<std::vector<double>>(Ksoil),
                                        as<std::vector<double>>(krhizo),
                                        lai, N, specificRootLength, rootTissueDensity, rootLengthDensity));
}

//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_rhizosphereMaximumConductance")]]
NumericVector rhizosphereMaximumConductance(NumericVector Ksoil, NumericVector fineRootBiomass, double lai, double N,
                                    double specificRootLength, double rootTissueDensity,  
                                    double rootLengthDensity) {
  std::vector<double> krhizo = rhizosphereMaximumConductance_c(as<std::vector<double>>(Ksoil),
                                                         as<std::vector<double>>(fineRootBiomass),
                                                         lai, N,
                                                         specificRootLength, rootTissueDensity,
                                                         rootLengthDensity);
  return(Rcpp::wrap(krhizo));
}



//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_coarseRootSoilVolumeFromConductance")]]
double coarseRootSoilVolumeFromConductance(double Kmax_rootxylem, double VCroot_kmax, double Al2As,
                                           NumericVector v, NumericVector d, NumericVector rfc) {
  return(coarseRootSoilVolumeFromConductance_c(Kmax_rootxylem, VCroot_kmax, Al2As,
                                               as<std::vector<double>>(v),
                                               as<std::vector<double>>(d),
                                               as<std::vector<double>>(rfc)));
}


//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_coarseRootLengthsFromVolume")]]
NumericVector coarseRootLengthsFromVolume(double VolInd, NumericVector v, NumericVector d, NumericVector rfc) {
  std::vector<double> tl = coarseRootLengthsFromVolume_c(VolInd, 
                                                         as<std::vector<double>>(v),
                                                         as<std::vector<double>>(d),
                                                         as<std::vector<double>>(rfc));
  return(Rcpp::wrap(tl));
}

// List coarseRootRadialAxialLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
//   int nlayers = v.size();
//   double maxRootDepth = 0.0;
//   //Vertical lengths
//   NumericVector vl(nlayers), zini(nlayers);
//   for(int i=0;i<nlayers;i++) {
//     if(i==0) {
//       zini[i] = 0.0;
//     } else {
//       zini[i] = zini[i-1]+ d[i-1];
//     }
//     if(v[i]>0.0) {
//       vl[i] = zini[i]+ d[i]/2.0;
//       maxRootDepth +=d[i];
//     } else {
//       vl[i] = 0.0;
//     }
//     // Rcout<<vl[i]<<" ";
//   }
//   // Rcout<<"\n";
//   int nlayerseff = nlayers;
//   for(int i=(nlayers-1);i>=0;i--) if(vl[i]>0.0) nlayerseff = i;
//   
//   //Radial lengths
//   NumericVector r(nlayers, 0.0), rl(nlayers, 0.0);
//   double maxr = 0.0;
//   for(int i=0;i<nlayerseff;i++) {
//     r[i] = sqrt(v[i]/(d[i]*M_PI));
//     maxr = std::max(r[i],maxr); 
//   }
//   // Rcout<<maxr<<"\n";
//   for(int i=0;i<nlayerseff;i++) {
//     rl[i] = maxRootDepth*depthWidthRatio*(r[i]/maxr);
//     // Rcout<<rl[i]<<" ";
//   }
//   return(List::create(_["radial"] = rl, _["axial"] = vl));
// }

//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_coarseRootLengths")]]
NumericVector coarseRootLengths(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
  NumericVector l = Rcpp::wrap(coarseRootLengths_c(as<std::vector<double>>(v),
                                                   as<std::vector<double>>(d),
                                                   depthWidthRatio));
  return(l);
}

//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_coarseRootSoilVolume")]]
double coarseRootSoilVolume(NumericVector v, NumericVector d, double depthWidthRatio = 1.0) {
  return(coarseRootSoilVolume_c(as<std::vector<double>>(v),
                                as<std::vector<double>>(d),
                                depthWidthRatio));
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

List nonoverlapHorizontalProportions(NumericMatrix V) {
  int numCohorts = V.nrow();
  int numlayers = V.ncol();
  List l(numCohorts);
  for(int coh=0;coh<numCohorts;coh++) {
    NumericMatrix RHOP(numCohorts,numlayers);
    std::fill(RHOP.begin(), RHOP.end(), 0.0);
    for(int l=0;l<numlayers;l++) RHOP(coh,l) = 1.0; 
    RHOP.attr("dimnames") = V.attr("dimnames");
    l[coh] = RHOP;
  }
  l.attr("names") = rownames(V);
  return(l);
}
List equaloverlapHorizontalProportions(NumericVector poolProportions, NumericMatrix V) {
  int numCohorts = V.nrow();
  int numlayers = V.ncol();
  List l(numCohorts);
  for(int coh=0;coh<numCohorts;coh++) {
    NumericMatrix RHOP(numCohorts,numlayers);
    for(int l=0;l<numlayers;l++) RHOP(_,l) = poolProportions;
    RHOP.attr("dimnames") = V.attr("dimnames");
    l[coh] = RHOP;
  }
  l.attr("names") = rownames(V);
  return(l);
}
//' @rdname root
//' @keywords internal
// [[Rcpp::export("root_horizontalProportions")]]
List horizontalProportions(NumericVector poolProportions, NumericVector VolInd, NumericVector N, NumericMatrix V, 
                           NumericVector d, NumericVector rfc) {
  
  int numCohorts = V.nrow();
  int numlayers = V.ncol();
  List l(numCohorts);
  NumericVector poolAreaInd(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    poolAreaInd[c] = 10000.0*poolProportions[c]/N[c]; //area of the pool per individual of the cohort
    // Rcout<<poolAreaInd[c]<<"\n";
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

