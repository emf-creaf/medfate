#include "inner_sperry_c.h"
#include "hydraulics_c.h"
#include "tissuemoisture_c.h"


double ludcmp_c(arma::mat& a, int n, std::vector<int>& indx) {
  const double TINY=1.0e-20;
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  std::vector<double> vv(n);
  double d = 1.0;
  for(i=0;i<n;i++) {
    big = 0.0;
    for(j=0;j<n;j++) if((temp=std::abs(a(i,j)))>big) big=temp;
    if(big==0.0) throw std::range_error("Singular matrix in routine ludcmp");
    vv[i] = 1.0/big; //Save the scaling
  }
  //Loop over columns of Crout's method
  for(j=0;j<n;j++){
    for(i=0;i<j;i++) {
      sum=a(i,j);
      for(k=0;k<i;k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
    }
    big=0.0;
    for(i=j;i<n;i++) {
      sum=a(i,j);
      for(k=0;k<j;k++) sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
      if((dum=vv[i]*std::abs(sum))>=big) {
        big=dum;
        imax=i;
      }
    }
    if(j!=imax) {
      for(k=0;k<n;k++){
        dum=a(imax,k);
        a(imax,k) = a(j,k);
        a(j,k) = dum;
      }
      d=-d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a(j,j)==0.0) a(j,j) = TINY;
    if(j!=n) {
      dum=1.0/(a(j,j));
      for(i=j+1;i<n;i++) a(i,j)*=dum;
    }
  }
  return(d);
}

void lubksb_c(const arma::mat& a, const int n, const std::vector<int>& indx, 
              std::vector<double>& b) {
  int i,ii=-1,ip,j;
  double sum;
  for(i = 0;i<n;i++){
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii>=0) for(j=ii;j<=i-1;j++) sum-=a(i,j)*b[j];
    else if(sum) ii=i;
    b[i] = sum;
  }
  for(i=(n-1);i>=0;i--) {
    sum=b[i];
    for(j=i+1;j<n;j++) sum-=a(i,j)*b[j];
    b[i] = sum/a(i,i);
  }
}

void E2psiBelowground_c(NetworkSteadyState& nss, SperryNetwork& hydraulicNetwork, std::vector<double>& psiIni) {
  
  int nlayers = hydraulicNetwork.psisoil.size();

  int ntrial = hydraulicNetwork.sperryParams.ntrial;
  double psiTol = hydraulicNetwork.sperryParams.psiTol;
  double ETol = hydraulicNetwork.sperryParams.ETol;
  // Rcpp::Rcout<< "parameters: " << ntrial << " " << psiTol << " " << ETol << "\n";
  std::vector<double>& psiSoil = hydraulicNetwork.psisoil; 
  std::vector<double>& krhizomax = hydraulicNetwork.krhizomax; 
  std::vector<double>& nsoil = hydraulicNetwork.nsoil; 
  std::vector<double>& alphasoil = hydraulicNetwork.alphasoil; 
  std::vector<double>& krootmax = hydraulicNetwork.krootmax;
  double rootc = hydraulicNetwork.rootc;
  double rootd  = hydraulicNetwork.rootd; 
  
  //Initialize
  if(((int) psiIni.size())==(nlayers+1)){
    for(int l=0;l<(nlayers+1);l++) {
      nss.x[l] = psiIni[l];
    }
  } else{
    double minPsi = -0.00001;
    for(int l=0;l<nlayers;l++) {
      nss.x[l] =psiSoil[l];
      minPsi = std::min(minPsi, psiSoil[l]);
    }
    nss.x[nlayers] = minPsi;
  }
  // for(int l=0;l<nlayers +1;l++) {
  //   Rcpp::Rcout << nss.x[l]<<" ";
  // }
  // Rcpp::Rcout<<"\n";
  //Newton-Raphson algorithm
  std::vector<double> p(nlayers+1), fvec(nlayers+1);
  std::vector<int> indx(nlayers+1);
  arma::mat fjac(nlayers+1,nlayers+1);
  double Esum = 0.0;
  for(int k=0;k<ntrial;k++) {
    // Rcpp::Rcout<<"trial "<<k<<"\n";
    //Calculate steady-state flow functions
    Esum = 0.0;
    bool stop = false;
    for(int l=0;l<nlayers;l++) {
      nss.ERoot[l] = EXylem_c(nss.x[nlayers], nss.x[l], krootmax[l], rootc, rootd, true, 0.0);
      nss.ERhizo[l] = EVanGenuchten_c(nss.x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l]);
      fvec[l] = nss.ERhizo[l] - nss.ERoot[l];
      // Rcpp::Rcout<<" Erhizo"<<l<<": "<< nss.ERhizo[l]<< " psiSoil"<<l <<": "<< psiSoil[l]<< " krhizo max"<<l<<": "<< krhizomax[l]<< " nsoil "<< nsoil[l] << " alphasoil "<< alphasoil[l] << " kroot max"<<l<<": "<< krootmax[l]<< " Eroot"<<l<<": "<<nss.ERoot[l]<<" fvec: "<<fvec[l]<<"\n";
      Esum +=nss.ERoot[l];
    }
    fvec[nlayers] = Esum-nss.E;
    // Rcpp::Rcout<<"fvec_nlayers: "<<fvec[nlayers]<<"\n";
    //Fill Jacobian
    for(int l1=0;l1<nlayers;l1++) { //funcio
      for(int l2=0;l2<nlayers;l2++) { //derivada
        if(l1==l2) {
          fjac(l1,l2) = -vanGenuchtenConductance_c(nss.x[l2],krhizomax[l2], nsoil[l2], alphasoil[l2])-xylemConductance_c(nss.x[l2], krootmax[l2], rootc, rootd);  
        }
        else fjac(l1,l2) = 0.0;
      }
    }
    fjac(nlayers,nlayers) = 0.0;
    for(int l=0;l<nlayers;l++) { 
      fjac(l,nlayers) = xylemConductance_c(nss.x[nlayers], krootmax[l], rootc, rootd); //funcio l derivada psi_rootcrown
      fjac(nlayers,l) = xylemConductance_c(nss.x[l], krootmax[l], rootc, rootd);//funcio nlayers derivada psi_l
      // funcio nlayers derivada psi_rootcrown
      fjac(nlayers,nlayers) +=-xylemConductance_c(nss.x[nlayers], krootmax[l], rootc, rootd);
    }
    // for(int l1=0;l1<=nlayers;l1++) { //funcio
    //   for(int l2=0;l2<=nlayers;l2++) { //derivada
    //     Rcpp::Rcout<<fjac(l1,l2)<<" ";
    //   }
    //   Rcpp::Rcout<<"\n";
    // }
    //Check function convergence
    double errf = 0.0;
    for(int fi=0;fi<=nlayers;fi++) errf += std::abs(fvec[fi]);
    // Rcpp::Rcout<< "ERRF: "<< errf<<"\n";
    if(errf<=ETol) break;
    //Right-hand side of linear equations
    for(int fi=0;fi<=nlayers;fi++) p[fi] = -fvec[fi];
    //Solve linear equations using LU decomposition
    ludcmp_c(fjac,nlayers+1,indx);
    lubksb_c(fjac,nlayers+1,indx,p);
    //Check root convergence
    double errx = 0.0;
    for(int fi=0;fi<=nlayers;fi++) {
      errx +=std::abs(p[fi]);
      nss.x[fi]+=p[fi];
      nss.x[fi] = std::min(0.0, nss.x[fi]);
      if(nss.x[fi]<-40.0) {
        nss.x[fi] = medfate::NA_DOUBLE;
        stop = true;
      }
      // Rcpp::Rcout<<"("<<nss.x[fi]<<") ";
    }
    // Rcpp::Rcout<<"\n";
    if(errx<=psiTol) break;
    else if(k==(ntrial-1)) { //Last trial and no convergence
      for(int fi=0;fi<=nlayers;fi++) nss.x[fi] = medfate::NA_DOUBLE;
      // Rcout<<"LC";
      stop = true;
    }
    if(stop) break;
  }
  
  
  //Initialize and copy output
  for(int l=0;l<nlayers;l++) {
    nss.psiRhizo[l] = nss.x[l];
  }
  nss.psiRootCrown = nss.x[nlayers];
  
  //Calculate final flows
  Esum = 0.0;
  for(int l=0;l<(nlayers-1);l++) {
    nss.ERhizo[l] = EVanGenuchten_c(nss.x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l]);
    Esum += nss.ERhizo[l];
  }
  nss.ERhizo[nlayers-1] = nss.E - Esum; //Define as difference to match input
} 

void E2psiAboveground_c(NetworkSteadyState& nss, SperryNetwork& hydraulicNetwork) {
  double kstemmax = hydraulicNetwork.kstemmax;
  double stemc = hydraulicNetwork.stemc;
  double stemd = hydraulicNetwork.stemd;
  double kleafmax = hydraulicNetwork.kleafmax;
  double leafc = hydraulicNetwork.leafc;
  double leafd = hydraulicNetwork.leafd;
  double PLCstem = hydraulicNetwork.PLCstem;
  double PLCleaf = hydraulicNetwork.PLCleaf;

  double psiPLCStem =  0.0;
  double psiPLCLeaf= 0.0;
  if(hydraulicNetwork.sperryParams.stemCavitationEffects) psiPLCStem = apoplasticWaterPotential_c(std::max(0.0001, 1.0 - PLCstem), stemc, stemd);
  if(hydraulicNetwork.sperryParams.leafCavitationEffects) psiPLCLeaf = apoplasticWaterPotential_c(std::max(0.0001,1.0 - PLCleaf), leafc, leafd);
  nss.psiStem = E2psiXylem_c(nss.E, nss.psiRootCrown, kstemmax, stemc, stemd, psiPLCStem); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
  nss.psiLeaf = E2psiXylem_c(nss.E, nss.psiStem, kleafmax, leafc, leafd, psiPLCLeaf); 
}

void E2psiNetwork_c(NetworkSteadyState& nss, SperryNetwork& hydraulicNetwork,
                    std::vector<double>& psiIni) {
  E2psiBelowground_c(nss, hydraulicNetwork, psiIni);
  if(!std::isnan(nss.psiRootCrown)) E2psiAboveground_c(nss, hydraulicNetwork);
} 

void fillSupplyFunctionNetwork_c(SperryNetwork&  hydraulicNetwork, double minFlow, double pCrit) {

  int maxNsteps  = hydraulicNetwork.sperryParams.maxNsteps;
  double ETol = hydraulicNetwork.sperryParams.ETol;
  std::vector<double>& psiSoil = hydraulicNetwork.psisoil;
  int nlayers = psiSoil.size();

  NetworkSteadyState sol(nlayers, minFlow);
  SupplyFunction SF(nlayers, maxNsteps);
  std::vector<double> psiIni = {0.0};
  E2psiNetwork_c(sol, hydraulicNetwork, psiIni);
  for(int l=0;l<nlayers;l++) {
    SF.ERhizo(0,l) = sol.ERhizo[l];
    SF.psiRhizo(0,l) = sol.psiRhizo[l];
  }
  SF.psiStem[0] = sol.psiStem;
  SF.psiLeaf[0] = sol.psiLeaf;
  SF.psiRoot[0] = sol.psiRootCrown;
  SF.E[0] = minFlow;
  //Store solution
  std::vector<double> prevX = sol.x;
  
  //Calculate initial slope
  sol.E = minFlow+ETol*2.0;
  E2psiNetwork_c(sol, hydraulicNetwork, prevX);
  double psiLeafI = sol.psiLeaf;
  double maxdEdp = (ETol*2.0)/std::abs(psiLeafI - SF.psiLeaf[0]);
  // Rcpp::Rcout<< " maxdEdp " << maxdEdp <<"\n";

  int nsteps = 1;
  double dE = std::min(0.0005,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // Rcpp::Rcout<< " step: "<< i;
    
    SF.E[i] = SF.E[i-1]+dE;
    sol.E = SF.E[i];
    E2psiNetwork_c(sol, hydraulicNetwork, prevX);
    //Replace stored solution
    prevX = sol.x;
    for(int l=0;l<nlayers;l++) {
      SF.ERhizo(i,l) = sol.ERhizo[l];
      SF.psiRhizo(i,l) = sol.psiRhizo[l];
    }
    SF.psiStem[i] = sol.psiStem;
    SF.psiLeaf[i] = sol.psiLeaf;
    SF.psiRoot[i] = sol.psiRootCrown;

    if(!std::isnan(SF.psiLeaf[i])) {
      // Rcpp::Rcout<<" psiRC: "<< SF.psiRoot[i] << " psileaf: "<< SF.psiLeaf[i] <<"\n";
      if(i==1) {
        SF.dEdp[0] = (SF.E[1]-SF.E[0])/std::abs(SF.psiLeaf[1] - SF.psiLeaf[0]);
      } else {
        double d1 = (SF.E[i-1]-SF.E[i-2])/std::abs(SF.psiLeaf[i-1] - SF.psiLeaf[i-2]);
        double d2 = (SF.E[i]-SF.E[i-1])/std::abs(SF.psiLeaf[i] - SF.psiLeaf[i-1]);
        SF.dEdp[i-1] = (d1+d2)/2.0;
      }
      if(SF.E[i]>0.1) dE = std::min(0.05,SF.dEdp[i-1]*0.05);
      else if(SF.E[i]>0.05) dE = std::min(0.01,SF.dEdp[i-1]*0.05);
      else if(SF.E[i]>0.01) dE = std::min(0.005,SF.dEdp[i-1]*0.05);
      nsteps++;
      if(SF.dEdp[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) SF.dEdp[nsteps-1] = (SF.E[nsteps-1]-SF.E[nsteps-2])/std::abs(SF.psiLeaf[nsteps-1] - SF.psiLeaf[nsteps-2]);
  //Copy values to nsteps
  SupplyFunction def(nlayers, nsteps);
  for(int i=0;i<nsteps;i++) {
    if(std::isnan(SF.E[i])) throw medfate::MedfateInternalError("NA E in supplyFunctionNetwork");
    def.E[i] = SF.E[i];
    def.dEdp[i] = SF.dEdp[i];
    def.psiRoot[i] = SF.psiRoot[i];
    for(int l=0;l<nlayers;l++) {
      def.ERhizo(i,l) = SF.ERhizo(i,l);
      def.psiRhizo(i,l) = SF.psiRhizo(i,l);
    }
    def.psiStem[i] = SF.psiStem[i];
    def.psiLeaf[i] = SF.psiLeaf[i];
  }

  // Assign to hydraulic network supply function
  hydraulicNetwork.supply = def;
}

void initSperryNetwork_inner_c(SperryNetwork& network,
                               int c,
                               const InternalWater& internalWater, 
                               const TranspirationParams& paramsTranspiration, 
                               const WaterStorageParams& paramsWaterStorage,
                               const std::vector<double>& VCroot_kmax, 
                               const std::vector<double>& VGrhizo_kmax,
                               const std::vector<double>& psiSoil, 
                               const std::vector<double>& VG_n, 
                               const std::vector<double>& VG_alpha,
                               const ControlParameters& control,
                               double sapFluidityDay) {
  
  int nlayers = VG_n.size();
  
  network.sperryParams = control.sperry;
  network.psisoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.krhizomax = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.nsoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.alphasoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.krootmax = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  for(int l=0;l<nlayers;l++) {
    network.psisoil[l] = psiSoil[l];
    network.krhizomax[l] = VGrhizo_kmax[l];
    network.krootmax[l] = sapFluidityDay*VCroot_kmax[l];
    network.alphasoil[l] = VG_alpha[l];
    network.nsoil[l] = VG_n[l];
  }
  network.kstemmax = sapFluidityDay*paramsTranspiration.VCstem_kmax[c];
  network.kleafmax = sapFluidityDay*paramsTranspiration.VCleaf_kmax[c];
  network.kleafapomax = sapFluidityDay*paramsTranspiration.VCleafapo_kmax[c];
  network.kleafsymp = sapFluidityDay*paramsTranspiration.kleaf_symp[c];
  network.stemc = paramsTranspiration.VCstem_c[c];
  network.stemd = paramsTranspiration.VCstem_d[c];
  network.leafc = paramsTranspiration.VCleaf_c[c];
  network.leafd = paramsTranspiration.VCleaf_d[c];
  network.rootc = paramsTranspiration.VCroot_c[c];
  network.rootd = paramsTranspiration.VCroot_d[c];
  if(control.sperry.leafCavitationEffects) {
    network.PLCleaf = internalWater.LeafPLC[c];
  } else {
    network.PLCleaf = 0.0;
  }
  network.PLCstem = internalWater.StemPLC[c];
  
  // Calculates supply function
  fillSupplyFunctionNetwork_c(network, 0.0, 0.001); 
}
