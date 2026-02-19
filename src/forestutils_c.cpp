#include "RcppArmadillo.h"
#include "medfate.h"
#include "forestutils_c.h"
#include "incgamma_c.h"
#include "medfate.h"
#include "strings.h"
#include "vector"

double leafAreaProportion_c(double z1, double z2, double zmin, double zmax) {
  double mu = (zmax+zmin)/2.0;
  double sd15 = (zmax-zmin)/2.0;
  double sd = sd15/1.5;
  z1 = std::max(z1, zmin);
  z2 = std::max(z2, zmin);
  z1 = std::min(z1, zmax);
  z2 = std::min(z2, zmax);
  double x1 = (z1-mu)/sd;
  double x2 = (z2-mu)/sd;
  double p1 = 0.5*(1.0+errorfunction_c(x1/sqrt(2.0), false, false));
  double p2 = 0.5*(1.0+errorfunction_c(x2/sqrt(2.0), false, false));
  double v = (p2-p1)/0.8663856; //truncated to -1.5 to 1.5
  return(v);
}

void updateLAIdistributionVectors_c(arma::mat& LAIdist, 
                                    const std::vector<double>& z, 
                                    const std::vector<double>& LAI, 
                                    const std::vector<double>& H, 
                                    const std::vector<double>& CR) {
  int ncanlayers = LAIdist.n_rows;
  int numCohorts = LAIdist.n_cols;
  for(int ci=0;ci<numCohorts;ci++) {
    double cbh = H[ci]*(1.0-CR[ci]);
    for(int hi=0;hi<ncanlayers;hi++) {
      if(z[hi]<= H[ci]) {
        LAIdist(hi,ci) = LAI[ci]*leafAreaProportion_c(z[hi],z[hi+1], cbh,H[ci]);
      } else {
        LAIdist(hi,ci) = 0.0;
      }
    }
  }
}

std::vector<std::string> cohortType_c(const std::vector<std::string> IDs) {
  std::vector<std::string> types(IDs.size());
  int numCohorts = IDs.size();
  for(int i=0;i<numCohorts;i++) {
    if(IDs[i].substr(0,1) =="T") {
      types[i] = "tree";
    } else if(IDs[i].substr(0,1) =="S") {
      types[i] = "shrub";
    } else if(IDs[i].substr(0,1) =="H") {
      types[i] = "herb";
    } else {
      throw medfate::MedfateInternalError("Wrong ID start");
    }
  }
  return(types);
}

std::vector<double> treeBasalArea_c(const std::vector<double>& N, const std::vector<double>& dbh) {
  int ncoh = N.size(); //N is density of individuals (ind/ha) in the cell
  std::vector<double> BA(ncoh, medfate::NA_DOUBLE); 
  for(int i=0;i<ncoh;i++) {
    if(!std::isnan(dbh[i])) BA[i] = N[i]*3.141593*pow(dbh[i]/200,2.0); //Basal area in m2/ha
  }
  return(BA);
}

std::vector<double> largerTreeBasalArea_c(const std::vector<double>& N, const std::vector<double>& dbh, 
                                          double self_include_prop) {
  int ncoh = N.size();
  std::vector<double> BA = treeBasalArea_c(N, dbh); 
  std::vector<double> ltBA(ncoh, medfate::NA_DOUBLE);
  for(int i=0;i<ncoh;i++) {
    if(!std::isnan(BA[i])) {
      ltBA[i] = 0.0;
      for(int j=0;j<ncoh;j++) {
        if(i==j) ltBA[i] += (BA[j]*self_include_prop); //add half of its own basal area
        else if(dbh[j]>dbh[i]) ltBA[i] += BA[j];
      }
    }
  }
  return(ltBA);
}

double crownCompetitionFactorAllometric_c(const std::vector<double>& N, const std::vector<double>& dbh, 
                                          const std::vector<double>& Acw, const std::vector<double>& Bcw) {
  int ntree = N.size();
  double ccf = 0.0;
  for(int i=0;i<ntree;i++) {
    if(!std::isnan(dbh[i])) {
      double cw = Acw[i]*pow(dbh[i], Bcw[i]);
      ccf = ccf + (N[i]*M_PI*pow(cw/2.0,2.0)/100.0);
    }
  }
  return(ccf);
}

std::vector<double> treeCrownRatioAllometric_c(std::vector<double>& N, std::vector<double>& dbh, std::vector<double>& H, 
                                               const std::vector<double>& Acw, const std::vector<double>& Bcw,
                                               const std::vector<double>& Acr, const std::vector<double>& B1cr, 
                                               const std::vector<double>& B2cr, 
                                               const std::vector<double>& B3cr,
                                               const std::vector<double>& C1cr, const std::vector<double>& C2cr) {
  std::vector<double> BAL = largerTreeBasalArea_c(N, dbh);
  double ccf = crownCompetitionFactorAllometric_c(N, dbh, Acw, Bcw);
  // Rcout<<ccf<<"\n";
  int ntree = N.size();
  std::vector<double> treeCR(ntree, medfate::NA_DOUBLE);
  for(int i=0;i<ntree;i++) {
    if(!std::isnan(dbh[i])) {
      double lm = Acr[i]+ B1cr[i]*(H[i]/(100.0*dbh[i]))+B2cr[i]*(H[i]/100.0)+B3cr[i]*pow(dbh[i],2.0)+C1cr[i]*BAL[i]+C2cr[i]*log(ccf);
      treeCR[i] = 1.0/(1.0+ exp(-1.0*lm));
    }
  }
  return(treeCR);
}