#include "inner_sperry_c.h"
#include "hydraulics_c.h"
#include "photosynthesis_c.h"
#include "tissuemoisture_c.h"
#include "transpiration_advanced_c.h"
#include "biophysicsutils_c.h"



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

void E2psiBelowground_c(NetworkSteadyState& nss, SperryNetwork& hydraulicNetwork, const std::vector<double>& psiIni) {
  
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
                    const std::vector<double>& psiIni) {
  E2psiBelowground_c(nss, hydraulicNetwork, psiIni);
  if(!std::isnan(nss.psiRootCrown)) E2psiAboveground_c(nss, hydraulicNetwork);
} 

SupplyFunction buildSupplyFunctionNetwork_c(SperryNetwork& hydraulicNetwork, double minFlow, double pCrit) {

  int maxNsteps  = hydraulicNetwork.sperryParams.maxNsteps;
  double ETol = hydraulicNetwork.sperryParams.ETol;
  std::vector<double>& psiSoil = hydraulicNetwork.psisoil;
  int nlayers = psiSoil.size();
  
  SupplyFunction SF = SupplyFunction(nlayers, maxNsteps);
  
  NetworkSteadyState sol(nlayers, minFlow);
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
  // std::vector<double> prevX = sol.x;
  
  //Calculate initial slope
  sol.E = minFlow+ETol*2.0;
  E2psiNetwork_c(sol, hydraulicNetwork, sol.x);
  double psiLeafI = sol.psiLeaf;
  double maxdEdp = (ETol*2.0)/std::abs(psiLeafI - SF.psiLeaf[0]);
  // Rcpp::Rcout<< " maxdEdp " << maxdEdp <<"\n";

  SF.nsteps = 1;
  double dE = std::min(0.0005,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // Rcpp::Rcout<< " step: "<< i;

    SF.E[i] = SF.E[i-1]+dE;
    sol.E = SF.E[i];
    E2psiNetwork_c(sol, hydraulicNetwork, sol.x);
    
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
        SF.dEdP[0] = (SF.E[1]-SF.E[0])/std::abs(SF.psiLeaf[1] - SF.psiLeaf[0]);
      } else {
        double d1 = (SF.E[i-1]-SF.E[i-2])/std::abs(SF.psiLeaf[i-1] - SF.psiLeaf[i-2]);
        double d2 = (SF.E[i]-SF.E[i-1])/std::abs(SF.psiLeaf[i] - SF.psiLeaf[i-1]);
        SF.dEdP[i-1] = (d1+d2)/2.0;
      }
      if(SF.E[i]>0.1) dE = std::min(0.05,SF.dEdP[i-1]*0.05);
      else if(SF.E[i]>0.05) dE = std::min(0.01,SF.dEdP[i-1]*0.05);
      else if(SF.E[i]>0.01) dE = std::min(0.005,SF.dEdP[i-1]*0.05);
      SF.nsteps++;
      if(SF.dEdP[i-1]<(pCrit*maxdEdp)) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(SF.nsteps>1) SF.dEdP[SF.nsteps-1] = (SF.E[SF.nsteps-1]-SF.E[SF.nsteps-2])/std::abs(SF.psiLeaf[SF.nsteps-1] - SF.psiLeaf[SF.nsteps-2]);
  
  return(SF);
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
  
  // Rcpp::Rcout<< "parameters: " << network.sperryParams.ntrial << " " << network.sperryParams.psiTol << " " << network.sperryParams.ETol << "\n";
  // Rcpp::Rcout<< " Stem PLC " << network.PLCstem << " Leaf PLC " << network.PLCleaf<< "\n";
  // Rcpp::Rcout<< " kstemmax " << network.kstemmax << " kleafmax " << network.kleafmax<< " kleafapomax " << network.kleafapomax<< " kleafsymp " << network.kleafsymp<< "\n";
  // Rcpp::Rcout<< " root c " << network.rootc << " root d " << network.rootd<< "\n";
  // Rcpp::Rcout<< " stem c " << network.stemc << " stem d " << network.stemd<< "\n";
  // Rcpp::Rcout<< " leaf c " << network.leafc << " leaf d " << network.leafd<< "\n";
  // for(int l=0;l<nlayers; l++) {
  //   Rcpp::Rcout<< " l "<< l << " psisoil " << network.psisoil[l] << 
  //     " krootmax " << network.krootmax[l] << " krhizomax "<< network.krhizomax[l] <<
  //       " alpha " << network.alphasoil[l] << " n "<< network.nsoil[l] << "\n";
  // }
}


void profitMaximization2_c(ProfitMaximization& PM, 
                           SupplyFunction& supply, int initialPos,
                           double Catm, double Patm, double Tair, double vpa, double u, 
                           double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                           double leafWidth, double refLeafArea,
                           double Gswmin, double Gswmax) {
  
  int nsteps = supply.nsteps;
  
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  
  int valid = 1;
  for(int i=0;i<nsteps-1;i++) {
    if((i>0) && (supply.E[i] > supply.E[i-1])) valid++; 
    mindEdp = std::min(mindEdp, supply.dEdP[i]);
    maxdEdp = std::max(maxdEdp, supply.dEdP[i]);
  }
  nsteps = valid;
  // Rcpp::Rcout<< " valid " << valid<< " maxdEdp "<< maxdEdp << " "<<mindEdp<< " mindEdp "<< mindEdp/maxdEdp<<"\n";
  // initial pos cannot be over the valid steps
  initialPos = std::min(initialPos, nsteps-1);
  
  //Evaluate photosynthesis and profit at Agmax
  PhotoFunction photoAgMax;
  leafPhotosynthesisOneFunction2_c(photoAgMax, supply.E[nsteps-1], supply.psiLeaf[nsteps - 1], Catm, Patm, Tair, vpa, u, 
                                   SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                   leafWidth, refLeafArea);
  double Agmax = photoAgMax.GrossPhotosynthesis;
  double profitAgMax = (1.0 - ((maxdEdp - supply.dEdP[nsteps-1])/(maxdEdp - mindEdp))); 
  
  //Photosynthesis and profit maximization at current value of initialPos
  PhotoFunction photoInitial;
  leafPhotosynthesisOneFunction2_c(photoInitial, supply.E[initialPos], supply.psiLeaf[initialPos], Catm, Patm, Tair, vpa, u, 
                                   SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                   leafWidth, refLeafArea);
  double AgInitial = photoInitial.GrossPhotosynthesis;
  double profitInitial = (AgInitial/Agmax) - ((maxdEdp - supply.dEdP[initialPos])/(maxdEdp - mindEdp)); 
  if(Agmax ==0.0) profitInitial = - ((maxdEdp - supply.dEdP[initialPos])/(maxdEdp - mindEdp)); //To avoid 0/0 division
  
  PhotoFunction photoPrev, photoNext, photoCurrent;
  int iPos = initialPos;
  double profitPrev, profitNext, profitCurrent;
  
  photoCurrent = photoInitial;
  profitCurrent = profitInitial;
  
  int prevStep = 0;
  bool cont = true;
  // Rcpp::Rcout<< " initialPos " << initialPos;
  while(cont) {
    // Rcpp::Rcout<<".";
    if((iPos>0) && (prevStep <= 0)) {
      leafPhotosynthesisOneFunction2_c(photoPrev, supply.E[iPos-1], supply.psiLeaf[iPos-1], Catm, Patm, Tair, vpa, u, 
                                       SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                       leafWidth, refLeafArea);
      double AgPrev = photoPrev.GrossPhotosynthesis;
      profitPrev = (AgPrev/Agmax) - ((maxdEdp - supply.dEdP[iPos-1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitPrev = - ((maxdEdp - supply.dEdP[iPos-1])/(maxdEdp - mindEdp)); 
    } else {
      photoPrev = photoCurrent;
      profitPrev = profitCurrent;
    }
    if((iPos < (nsteps-1)) && (prevStep >= 0)) {
      leafPhotosynthesisOneFunction2_c(photoNext, supply.E[iPos+1], supply.psiLeaf[iPos+1], Catm, Patm, Tair, vpa, u, 
                                       SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                       leafWidth, refLeafArea);
      double AgNext = photoNext.GrossPhotosynthesis;
      profitNext = (AgNext/Agmax) - ((maxdEdp - supply.dEdP[iPos+1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitNext = - ((maxdEdp - supply.dEdP[iPos+1])/(maxdEdp - mindEdp)); 
    } else {
      photoNext = photoAgMax;
      profitNext = profitAgMax;
    }
    bool selDecr = ((profitPrev >= profitCurrent) && (photoPrev.Gsw >= Gswmin)) || (photoCurrent.Gsw>Gswmax);
    selDecr = selDecr && (prevStep<=0) && (iPos>0);
    bool selIncr = ((profitNext > profitCurrent) && (photoNext.Gsw <= Gswmax)) || (photoCurrent.Gsw<Gswmin);
    selIncr = selIncr && (prevStep>=0) && (iPos<(nsteps-1));
    if(selDecr) {
      profitCurrent = profitPrev;
      photoCurrent = photoPrev;
      iPos = iPos-1;
      prevStep = -1;
      // Rcpp::Rcout <<"-";
    } else if(selIncr) {
      profitCurrent = profitNext;
      photoCurrent = photoNext;
      iPos = iPos+1;
      // Rcpp::Rcout <<"+";
      prevStep = 1;
    } else {
      cont = false;
    }
  }
  // Rcpp::Rcout <<"\n";
  // Rcpp::Rcout<< " finalPos " << iPos<< " final profit "<< profitCurrent <<"\n";
  PM.photosynthesisFunction = photoCurrent;
  PM.Profit = profitCurrent;
  PM.iMaxProfit = iPos;
}

void copyRhizoPsi_c(int c, int iPM, 
                    arma::mat& RhizoPsi, 
                    arma::mat& RhizoPsiMAT,
                    const arma::Mat<uint8_t>& layerConnected, 
                    const std::vector<arma::mat>& RHOP, 
                    const std::vector<arma::Mat<uint8_t>>& layerConnectedPools,
                    const std::vector<double>& VCroot_c, 
                    const std::vector<double>& VCroot_d,  
                    bool plantWaterPools) {
  int nlayers = layerConnected.n_cols;
  int numCohorts = layerConnected.n_rows;
  
  if(!plantWaterPools) {
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected(c,l)) {
        RhizoPsiMAT(c,l) = RhizoPsi(iPM,cl);
        cl++;
      } 
    }
  } else {
    const arma::mat& RHOPcoh = RHOP[c];
    const arma::Mat<uint8_t>& layerConnectedCoh = layerConnectedPools[c];
    std::vector<double> rplv(numCohorts,medfate::NA_DOUBLE);
    std::vector<double> vplv(numCohorts,medfate::NA_DOUBLE);
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      int clj = 0;
      for(int j=0;j<numCohorts;j++) {
        if(layerConnectedCoh(j,l)) {
          rplv[clj] = RhizoPsi(iPM,cl);
          vplv[clj] = RHOPcoh(j,l);
          cl++;
          clj++;
        }
      }
      std::vector<double> pv(clj,medfate::NA_DOUBLE);
      std::vector<double> vv(clj,medfate::NA_DOUBLE);
      for(int j=0;j<clj;j++) {
        pv[j] = rplv[j];
        vv[j] = vplv[j];
        // Rcout<< c << " "<< l <<" "<<j<< " "<< pv[j]<< " "<< vv[j]<<"\n";
      }
      RhizoPsiMAT(c,l) = averagePsi_c(pv,vv,VCroot_c[c], VCroot_d[c]);
      // Rcout<< c << " "<<l<< " "<< RhizoPsiMAT(c,l)<<"\n";
    }
  }
}
void innerSperry_c(ModelInput& x,
                   std::vector<SperryNetwork>& networks, 
                   InnerTranspirationInput_COMM& input, 
                   AdvancedTranspiration_RESULT& output, int n, double tstep, 
                   int stepFunctions) {
  
  int numCohorts = x.cohorts.CohortCode.size();
  int nlayers = x.soil.getNlayers();
  
  // Extract control variables
  bool plantWaterPools = (x.control.rhizosphereOverlap!="total");

  // Extract parameters
  // //Water pools
  // DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  // List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  // NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  // List RHOP;
  // NumericVector poolProportions(numCohorts);
  // if(plantWaterPools) {
  //   RHOP = belowLayers["RHOP"];
  //   poolProportions = belowdf["poolProportions"];
  // }
  // NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  // 
  // //Extract output to be filled
  // // Rcout<<"output\n";
  // List outPhotoSunlit = output["PhotoSunlitFunctions"];
  // List outPhotoShade = output["PhotoShadeFunctions"];
  // List outPMSunlit = output["PMSunlitFunctions"];
  // List outPMShade = output["PMShadeFunctions"];
  // 
  // NumericMatrix SoilWaterExtract = Rcpp::as<Rcpp::NumericMatrix>(output["Extraction"]);
  // List ExtractionPools = Rcpp::as<Rcpp::List>(output["ExtractionPools"]);
  // NumericMatrix soilLayerExtractInst = Rcpp::as<Rcpp::NumericMatrix>(output["ExtractionInst"]);
  // 
  // NumericMatrix minPsiRhizo = Rcpp::as<Rcpp::NumericMatrix>(output["RhizoPsi"]);
  // 

  arma::mat& LAI_SH = output.shade_inst.LAI;
  arma::mat& Vmax298_SH = output.shade_inst.Vmax298;
  arma::mat& Jmax298_SH = output.shade_inst.Jmax298;
  arma::mat& SWR_SH = output.shade_inst.Abs_SWR;
  arma::mat& PAR_SH = output.shade_inst.Abs_PAR;
  arma::mat& LWR_SH = output.shade_inst.Net_LWR;

  arma::mat& LAI_SL = output.sunlit_inst.LAI;
  arma::mat& Vmax298_SL = output.sunlit_inst.Vmax298;
  arma::mat& Jmax298_SL = output.sunlit_inst.Jmax298;
  arma::mat& SWR_SL = output.sunlit_inst.Abs_SWR;
  arma::mat& PAR_SL = output.sunlit_inst.Abs_PAR;
  arma::mat& LWR_SL = output.sunlit_inst.Net_LWR;
  // NumericMatrix Ag_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ag"]);
  // NumericMatrix An_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["An"]);
  // NumericMatrix E_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["E"]);
  // NumericMatrix VPD_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["VPD"]);
  // NumericMatrix Psi_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Psi"]);
  // NumericMatrix Temp_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Temp"]);
  // NumericMatrix GSW_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Gsw"]);
  // NumericMatrix Ci_SL = Rcpp::as<Rcpp::NumericMatrix>(SunlitInst["Ci"]);
  // 
  // //Extract input  
  // std::vector<int>& iLayerCohort = input.iLayerCohort;
  std::vector<int>& iLayerSunlit = input.iLayerSunlit;
  std::vector<int>& iLayerShade = input.iLayerShade;
  // IntegerVector iPMSunlit = input["iPMSunlit"];
  // IntegerVector iPMShade = input["iPMShade"];
  // std::vector<int>& nlayerscon = input.nlayerscon;
  // LogicalMatrix layerConnected = input["layerConnected"];
  // List layerConnectedPools = input["layerConnectedPools"];
  
  for(int c=0;c<numCohorts;c++) { //Plant cohort loop
    SupplyFunction& supply = *networks[c].supply;
    if((x.above.LAI_expanded[c]>0.0) && (x.internalWater.LeafPLC[c] < 0.999)) { //Process transpiration and photosynthesis only if there are some leaves
      if(supply.E.size()>0) {
        double TPhase_gmin = x.control.advancedWB.TPhase_gmin;
        double Q10_1_gmin = x.control.advancedWB.Q10_1_gmin;
        double Q10_2_gmin = x.control.advancedWB.Q10_2_gmin;
        //Here canopy air temperature is used instead of Tleaf, because leaf energy balance has not been closed
        double gmin_SL = gmin_c(x.canopy.Tair[iLayerSunlit[c]], x.paramsTranspiration.Gswmin[c], TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        double gmin_SH = gmin_c(x.canopy.Tair[iLayerShade[c]], x.paramsTranspiration.Gswmin[c], TPhase_gmin, Q10_1_gmin, Q10_2_gmin);
        // Rcpp::Rcout<< " gmin_SL "<< gmin_SL<< " gmin_SH " << gmin_SH << " " << x.paramsTranspiration.Gswmax[c]<<"\n";
        //Photosynthesis function for sunlit and shade leaves
        ProfitMaximization PMSunlit;
        profitMaximization2_c(PMSunlit, supply, input.iPMSunlit[c],
                              x.canopy.Cair[iLayerSunlit[c]], input.Patm,
                              x.canopy.Tair[iLayerSunlit[c]], x.canopy.VPair[iLayerSunlit[c]], input.zWind[iLayerSunlit[c]],
                              SWR_SL(c,n), LWR_SL(c,n), irradianceToPhotonFlux_c(PAR_SL(c,n), defaultLambda),
                              Vmax298_SL(c,n), Jmax298_SL(c,n), x.paramsAnatomy.LeafWidth[c], LAI_SL(c,n),
                              gmin_SL, x.paramsTranspiration.Gswmax[c]);
        PhotoFunction& photoSunlit = PMSunlit.photosynthesisFunction;
        input.iPMSunlit[c] = PMSunlit.iMaxProfit;
        // Rcpp::Rcout << c << " / " << n << " iPMSunlit " << input.iPMSunlit[c] <<"\n";
        ProfitMaximization PMShade;
        profitMaximization2_c(PMShade, supply, input.iPMShade[c], 
                              x.canopy.Cair[iLayerShade[c]], input.Patm,
                              x.canopy.Tair[iLayerShade[c]], x.canopy.VPair[iLayerShade[c]], input.zWind[iLayerShade[c]],
                              SWR_SH(c,n), LWR_SH(c,n), irradianceToPhotonFlux_c(PAR_SH(c,n), defaultLambda),
                              Vmax298_SH(c,n), Jmax298_SH(c,n), x.paramsAnatomy.LeafWidth[c], LAI_SH(c,n),
                              gmin_SH, x.paramsTranspiration.Gswmax[c]);
        if(!x.control.advancedWB.sunlitShade) PMShade = PMSunlit;
        PhotoFunction& photoShade = PMShade.photosynthesisFunction;
        input.iPMShade[c] = PMShade.iMaxProfit;

  //       //Store?
  //       if(!IntegerVector::is_na(stepFunctions)) {
  //         if(n==stepFunctions) {
  //           outPhotoSunlit[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerSunlit[c]], Patm,
  //                                                           Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], 
  //                                                                                       zWind[iLayerSunlit[c]], 
  //                                                                                            SWR_SL(c,n), LWR_SL(c,n), 
  //                                                                                            irradianceToPhotonFlux_c(PAR_SL(c,n), defaultLambda), 
  //                                                                                            Vmax298_SL(c,n), Jmax298_SL(c,n), 
  //                                                                                            leafWidth[c], LAI_SL(c,n));
  //           outPhotoShade[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerShade[c]], Patm,
  //                                                          Tair[iLayerShade[c]], VPair[iLayerShade[c]], 
  //                                                                                     zWind[iLayerShade[c]], 
  //                                                                                          SWR_SH(c,n), LWR_SH(c,n), 
  //                                                                                          irradianceToPhotonFlux_c(PAR_SH(c,n), defaultLambda),
  //                                                                                          Vmax298_SH(c,n), Jmax298_SH(c,n), 
  //                                                                                          leafWidth[c], LAI_SH(c,n));
  //           outPMSunlit[c] = PMSunlit;
  //           outPMShade[c] = PMShade;
  //           if(!sunlitShade) {
  //             outPMShade[c] = outPMSunlit[c];
  //             outPhotoShade[c] = outPhotoSunlit[c];
  //           }
  //         }
  //       }
  //       // Rcout<<iPMSunlit[c]<<" "<<iPMShade[c] <<" "<<GwSunlit[iPMSunlit[c]]<<" "<<GwShade[iPMShade[c]]<<" "<<fittedE[iPMSunlit[c]]<<" "<<fittedE[iPMShade[c]]<<"\n";
        //Get leaf status
        output.shade_inst.E(c,n) = supply.E[input.iPMShade[c]];
        output.sunlit_inst.E(c,n) = supply.E[input.iPMSunlit[c]];
        output.shade_inst.Psi(c,n) = supply.psiLeaf[input.iPMShade[c]];
        output.sunlit_inst.Psi(c,n) = supply.psiLeaf[input.iPMSunlit[c]];
        output.shade_inst.An(c,n) = photoShade.NetPhotosynthesis;
        output.sunlit_inst.An(c,n) = photoSunlit.NetPhotosynthesis;
        output.shade_inst.Ag(c,n) = photoShade.GrossPhotosynthesis;
        output.sunlit_inst.Ag(c,n) = photoSunlit.GrossPhotosynthesis;
        output.shade_inst.Ci(c,n) = photoShade.Ci;
        output.sunlit_inst.Ci(c,n) = photoSunlit.Ci;
        output.shade_inst.Gsw(c,n)= photoShade.Gsw;
        output.sunlit_inst.Gsw(c,n)= photoSunlit.Gsw;
        output.shade_inst.VPD(c,n)= photoShade.LeafVPD;
        output.sunlit_inst.VPD(c,n)= photoSunlit.LeafVPD;
        output.shade_inst.Temp(c,n)= photoShade.LeafTemperature;
        output.sunlit_inst.Temp(c,n)= photoSunlit.LeafTemperature;

        //Scale photosynthesis
        double Agsum = (output.sunlit_inst.Ag(c,n)*LAI_SL(c,n) + output.shade_inst.Ag(c,n)*LAI_SH(c,n));
        double Ansum = (output.sunlit_inst.An(c,n)*LAI_SL(c,n) + output.shade_inst.An(c,n)*LAI_SH(c,n));
        output.plants_inst.Ag(c,n) = (1e-6)*12.01017*Agsum*tstep;
        output.plants_inst.An(c,n) = (1e-6)*12.01017*Ansum*tstep;

        //Average flow from sunlit and shade leaves
        double Eaverage = (output.sunlit_inst.E(c,n)*LAI_SL(c,n) + output.shade_inst.E(c,n)*LAI_SH(c,n))/(LAI_SL(c,n) + LAI_SH(c,n));

        //Find iPM for  flow corresponding to the  average flow
        double absDiff = 99999999.9;
        int iPM = -1;
        for(int k=0;k< (int) supply.E.size();k++){ //Only check up to the size of fittedE
          double adk = std::abs(supply.E[k]-Eaverage);
          // Rcpp::Rcout << supply.E[k] << "\n";
          if(adk<absDiff) {
            absDiff = adk;
            iPM = k;
          }
        }
        // Rcpp::Rcout<< c << " / " << n << " iPM " << iPM <<"\n";
        if(iPM==-1) {
          throw medfate::MedfateInternalError("iPM -1 (could not determine regulation)");
        }

        //Store instantaneous total conductance
        output.plants_inst.dEdP(c,n) = supply.dEdP[iPM];

        //Store instantaneous flow and leaf water potential
        x.internalWater.Einst[c] = supply.E[iPM]*input.f_dry;
        x.internalWater.LeafPsi[c] = supply.psiLeaf[iPM];
        x.internalWater.RootCrownPsi[c] = supply.psiRoot[iPM];

        //Scale from instantaneous flow to water volume in the time step
        output.plants_inst.E(c,n) = x.internalWater.Einst[c]*0.001*0.01802*x.above.LAI_expanded[c]*tstep;
        // Rcpp::Rcout<< c << " / " << n << " E: "<< output.plants_inst.E(c,n) << "\n";
        std::vector<double> Esoilcn(input.nlayerscon[c],0.0);
        std::vector<double> ElayersVEC(input.nlayerscon[c],0.0);


        //Get info from sFunctionBelow
        // NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["psiRhizo"]);
 
        //Store steady state stem and rootcrown and root surface water potential values
        x.internalWater.StemPsi[c] = supply.psiStem[iPM];
        for(int lc=0;lc<input.nlayerscon[c];lc++) {
          ElayersVEC[lc] = supply.ERhizo(iPM,lc)*tstep*input.f_dry; //Scale according to the time step
        }

        //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
        copyRhizoPsi_c(c,iPM,
                       supply.psiRhizo, x.belowLayers.RhizoPsi,
                       input.layerConnected,
                       x.belowLayers.RHOP, input.layerConnectedPools,
                       x.paramsTranspiration.VCroot_c, x.paramsTranspiration.VCroot_d,
                       plantWaterPools);
        
        x.internalWater.StemSympPsi[c] = x.internalWater.StemPsi[c]; //Stem symplastic compartment coupled with apoplastic compartment
        x.internalWater.LeafSympPsi[c] = x.internalWater.LeafPsi[c]; //Leaf symplastic compartment coupled with apoplastic compartment

        // Update the leaf and stem PLC
        if(x.control.sperry.stemCavitationEffects) {
          if(x.control.commonWB.stemCavitationRecovery!="total") {
            x.internalWater.StemPLC[c] = std::max(x.internalWater.StemPLC[c], 1.0 - xylemConductance_c(x.internalWater.StemPsi[c], 1.0, x.paramsTranspiration.VCstem_c[c], x.paramsTranspiration.VCstem_d[c]));
          } else { //Immediate refilling
            x.internalWater.StemPLC[c] = 1.0 - xylemConductance_c(x.internalWater.StemPsi[c], 1.0, x.paramsTranspiration.VCstem_c[c], x.paramsTranspiration.VCstem_d[c]);
          }
        }
        if(x.control.sperry.leafCavitationEffects) {
          if(x.control.commonWB.leafCavitationRecovery!="total") {
            x.internalWater.LeafPLC[c] = std::max(x.internalWater.LeafPLC[c], 1.0 - xylemConductance_c(x.internalWater.LeafPsi[c], 1.0, x.paramsTranspiration.VCleaf_c[c], x.paramsTranspiration.VCleaf_d[c]));
          } else { //Immediate refilling
            x.internalWater.LeafPLC[c] = 1.0 - xylemConductance_c(x.internalWater.LeafPsi[c], 1.0, x.paramsTranspiration.VCleaf_c[c], x.paramsTranspiration.VCleaf_d[c]);
          }
        }

        //Scale soil water extracted from leaf to cohort level
        for(int lc=0;lc<input.nlayerscon[c];lc++) {
          Esoilcn[lc] = ElayersVEC[lc]*0.001*0.01802*x.above.LAI_expanded[c]; //Scale from flow to water volume in the time step
        }

        //Balance between extraction and transpiration
        output.plants_inst.PWB(c,n) = std::accumulate(Esoilcn.begin(), Esoilcn.end(), 0.0) - output.plants_inst.E(c,n);

        //Add step transpiration to daily plant cohort transpiration
        output.plants.Transpiration[c] += output.plants_inst.E(c,n);
        output.plants.NetPhotosynthesis[c] += output.plants_inst.An(c,n);
        output.plants.GrossPhotosynthesis[c] += output.plants_inst.Ag(c,n);
        //Add PWB
        output.plants.WaterBalance[c] += output.plants_inst.PWB(c,n);

        //Copy transpiration and from connected layers to transpiration from soil layers
        //And update soil water content (soil water potential will not be updated until next day!)
        if(!plantWaterPools) {
          int cl = 0;
          for(int l=0;l<nlayers;l++) {
            if(input.layerConnected(c,l)) {
              output.extraction(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
              output.extractionInst(l,n) += Esoilcn[cl];
              cl++;
            }
          }
        } else {
          arma::Mat<uint8_t>&  layerConnectedCoh = input.layerConnectedPools[c];
          int cl = 0;
          for(int j = 0;j<numCohorts;j++) {
            arma::mat& ExtractionPoolsCoh = output.extractionPools[j];
            for(int l=0;l<nlayers;l++) {
              if(layerConnectedCoh(j,l)) {
                output.extraction(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
                output.extractionInst(l,n) += Esoilcn[cl];
                ExtractionPoolsCoh(c,l) += Esoilcn[cl];
                cl++;
              }
            }
          }
        }
      } else {
        output.shade_inst.Psi(c,n) = medfate::NA_DOUBLE;
        output.sunlit_inst.Psi(c,n) = medfate::NA_DOUBLE;
        output.shade_inst.Gsw(c,n)= medfate::NA_DOUBLE;
        output.sunlit_inst.Gsw(c,n)= medfate::NA_DOUBLE;
        output.shade_inst.VPD(c,n)= medfate::NA_DOUBLE;
        output.sunlit_inst.VPD(c,n)= medfate::NA_DOUBLE;
        output.shade_inst.Temp(c,n)= medfate::NA_DOUBLE;
        output.sunlit_inst.Temp(c,n)= medfate::NA_DOUBLE;
      }
    } else if(x.above.LAI_expanded[c]>0.0) { //Cohorts with living individuals but no LAI should be in equilibrium with soil (i.e. no transpiration)
      x.internalWater.RootCrownPsi[c] = supply.psiRoot[0];
      x.internalWater.StemPsi[c] = supply.psiStem[0];
      x.internalWater.StemSympPsi[c] = supply.psiStem[0];
      x.internalWater.LeafPsi[c] = supply.psiLeaf[0];
      x.internalWater.LeafSympPsi[c] = supply.psiLeaf[0];
    }
  } //End of cohort loop
}
