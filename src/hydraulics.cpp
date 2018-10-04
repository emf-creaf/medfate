#include "Rcpp.h"
#include "root.h"
#include "tissuemoisture.h"
#include <math.h>


using namespace Rcpp;
using namespace std;

double const maxPsi = -0.000001;
double const cmhead2MPa = 0.00009804139; //Constant to transform cm head to MPa

/**
 * Whole-plant conductance function (simple water balance model)
 */
// [[Rcpp::export("hydraulics.psi2K")]]
double Psi2K(double psi, double Psi_extract, double ws = 3.0) {
  return(exp(-0.6931472*pow(std::abs(psi/Psi_extract),ws)));
}
NumericVector Psi2K(double psi, NumericVector Psi_extract, double ws = 3.0) {
  int n = Psi_extract.size();
  NumericVector k(n);
  for(int i=0; i<n; i++) {
    k[i] = Psi2K(psi,Psi_extract[i],ws);
  }
  return k;
}

/**
 * Inverse of the whole-plant conductance function. Used to obtain the 'average' soil water
 * potential perceived by each plant cohort.
 */
// [[Rcpp::export("hydraulics.K2Psi")]]
double K2Psi(double K, double Psi_extract, double ws = 3.0) {
  double psi = Psi_extract*pow(log(K)/(-0.6931472),1.0/ws);
  if(psi>0.0) psi = -psi; //Usually psi_extr is a positive number
  return psi;
}
NumericVector K2Psi(NumericVector K, NumericVector Psi_extract, double ws = 3.0) {
  int n = Psi_extract.size();
  NumericVector psi(n);
  for(int i=0; i<n; i++) {
    psi[i] = K2Psi(K[i], Psi_extract[i], ws);
  }
  return psi;
}

// [[Rcpp::export("hydraulics.averagePsi")]]
double averagePsi(NumericVector psi, NumericVector v, double c, double d) {
  int nlayers = psi.size();
  NumericVector K(nlayers);
  for(int l=0;l<nlayers;l++) K[l]= exp(-0.6931472*pow(std::abs(psi[l]/d),c));
  double psires =  d*pow(log(sum(K*v))/(-0.6931472),1.0/c);
  return(psires);
}

// [[Rcpp::export(".gammds")]]
double gammds ( double x, double p)
  
  //****************************************************************************80
  //
  //  Purpose:
  //
  //    GAMMDS computes the incomplete Gamma integral.
  //
  //  Discussion:
  //
  //    The parameters must be positive.  An infinite series is used.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license. 
  //
  //  Modified:
  //
  //    22 January 2008
  //
  //  Author:
  //
  //    Original FORTRAN77 version by Chi Leung Lau.
  //    C++ version by John Burkardt.
  //
  //  Reference:
  //
  //    Chi Leung Lau,
  //    Algorithm AS 147:
  //    A Simple Series for the Incomplete Gamma Integral,
  //    Applied Statistics,
  //    Volume 29, Number 1, 1980, pages 113-114.
  //
  //  Parameters:
  //
  //    Input, double X, P, the arguments of the incomplete
  //    Gamma integral.  X and P must be greater than 0.
  //
  //    Output, double GAMMDS, the value of the incomplete
  //    Gamma integral.
  //
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  double uflo = 1.0E-37;
  double value;
  //
  //  Check the input.
  //
  if ( x <= 0.0 )
  {
    stop("x <= 0.0 in gammds");
    value = 0.0;
    return value;
  }
  
  if ( p <= 0.0 ) 
  {
    stop("p <= 0.0 in gammds");
    value = 0.0;
    return value;
  }
  //
  //  LGAMMA is the natural logarithm of the gamma function.
  //
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;
  
  if ( arg < log ( uflo ) )
  {
    // stop("underflow during the computation in gammds");
    value = NA_REAL;
    return value;
  }
  
  f = exp ( arg );
  
  if ( f == 0.0 )
  {
    // stop("underflow during the computation in gammds");
    value = NA_REAL;
    return value;
  }
  
  //
  //  Series begins.
  //
  c = 1.0;
  value = 1.0;
  a = p;
  
  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;
    
    if ( c <= e * value )
    {
      break;
    }
  }
  
  value = value * f;
  
  return value;
}

// [[Rcpp::export("hydraulics.xylemConductance")]]
double xylemConductance(double psi, double kxylemmax, double c, double d) {
  if(psi>0.0) {
    Rcout<< psi<<"\n";
    stop("psi has to be negative"); 
  } else if(psi==0.0) return(kxylemmax);
  return(kxylemmax*exp(-pow(psi/d,c)));
}

// [[Rcpp::export("hydraulics.xylemPsi")]]
double xylemPsi(double kxylem, double kxylemmax, double c, double d) {
  return(d*pow(-log(kxylem/kxylemmax),1.0/c));
}

// [[Rcpp::export("hydraulics.psiCrit")]]
double psiCrit(double c, double d) {
  return(d * pow(-log(0.01), 1.0/c));
}

/**
 * Van genuchten-mualem conductance equation (m = 1 - 1/n; l = 0.5)
 */
// [[Rcpp::export("hydraulics.vanGenuchtenConductance")]]
double vanGenuchtenConductance(double psi, double krhizomax, double n, double alpha) {
  double v = 1.0/(pow(alpha*std::abs(psi),n)+1.0);
  return(krhizomax*pow(v,(n-1.0)/(2.0*n))*pow(pow((1.0-v),(n-1.0)/n)-1.0,2.0));
}

// [[Rcpp::export(".Egamma")]]
double Egamma(double psi, double kxylemmax, double c, double d, double psiCav = 0.0) {
  if(psi>0.0) stop("psi has to be negative");
  else if(psi==0.0) return(0.0);
  double h = 1.0/c;
  double z = pow(psi/d,c);
  double g = exp(lgamma(h))*gammds(z,h); //Upper incomplete gamma, without the normalizing factor
  double E = kxylemmax*(-d/c)*g;
  if(psiCav<0.0) { //Decrease E from 0 to psiCav (avoid recursiveness!)
    if(psiCav < psi) {
      E = xylemConductance(psiCav,kxylemmax,c,d)*(-psi); //square integral
    } else {
      double Epsimin = kxylemmax*(-d/c)*exp(lgamma(h))*gammds(pow(psiCav/d,c),h);
      E = E - Epsimin + xylemConductance(psiCav,kxylemmax,c,d)*(-psiCav); //Remove part of the integral corresponding to psimin and add square integral
    }
  }
  return(E);
}


/*
 * Integral of the xylem vulnerability curve
 */
// [[Rcpp::export("hydraulics.EXylem")]]
double EXylem(double psiPlant, double psiUpstream, double kxylemmax, double c, double d, bool allowNegativeFlux = true, double psiCav = 0.0) {
  if((psiPlant > psiUpstream) & !allowNegativeFlux) ::Rf_error("Downstream potential larger (less negative) than upstream potential");
  return(Egamma(psiPlant, kxylemmax, c, d, psiCav)-Egamma(psiUpstream, kxylemmax, c,d, psiCav));
}


// [[Rcpp::export("hydraulics.EVanGenuchten")]]
double EVanGenuchten(double psiRhizo, double psiSoil, double krhizomax, double n, double alpha, double psiStep = -0.001, double psiTol = 0.0001, bool allowNegativeFlux = true) {
  if((psiRhizo>psiSoil) & !allowNegativeFlux) ::Rf_error("Downstream potential larger (less negative) than upstream potential");
  bool reverse = false;
  if(psiRhizo>psiSoil) reverse = true;
  if(reverse) {
    double tmp = psiSoil;
    psiSoil = psiRhizo;
    psiRhizo = tmp;
  }
  double psi = psiSoil;
  double vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
  double E = 0.0, vgPrev = vg;
  psiStep = std::max(psiStep, (psiRhizo-psiSoil)/10.0); //Check that step is not too large
  do {
    psi = psi + psiStep;
    if(psi>psiRhizo) {
      vgPrev = vg;
      vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
      E += ((vg+vgPrev)/2.0)*std::abs(psiStep);
    } else {
      psi = psi - psiStep; //retrocedeix
      psiStep = psiStep/2.0; //canvia pas
    }
  } while (std::abs(psi-psiRhizo)>psiTol);
  if(reverse) E = -E;
  return(E);
}

// [[Rcpp::export("hydraulics.Ecrit")]]
double Ecrit(double psiUpstream, double kxylemmax, double c, double d) {
  return(EXylem(psiCrit(c,d), psiUpstream, kxylemmax, c, d));
}




// [[Rcpp::export("hydraulics.E2psiXylem")]]
double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0, 
                   double psiStep = -0.01, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiUpstream);
  double psi = psiUpstream;
  double k = xylemConductance(std::min(psi, psiCav),kxylemmax, c, d);
  double Eg = 0.0;
  double psiPrev = psi;
  double kprev = k;
  while(Eg<E) {
    psiPrev = psi;
    kprev = k;
    psi = psi + psiStep;
    k = xylemConductance(std::min(psi, psiCav),kxylemmax, c, d);
    Eg = Eg + (-1.0*psiStep)*((kprev+k)/2.0);
    if(psi<psiMax) return(NA_REAL);
    if(NumericVector::is_na(Eg)) return(NA_REAL);
  }
  return((psiPrev+psi)/2.0);
}
// double E2psiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiCav = 0.0,
//                   double psiStep = -0.01, double psiMax = -10.0) {
//   if(E<0.0) stop("E has to be positive");
//   if(E==0) return(psiUpstream);
//   double psi = psiUpstream;
//   double Eg = 0.0;
//   double psiPrev = psi;
//   while(Eg<E) {
//     psiPrev = psi;
//     psi = psi + psiStep;
//     Eg = EXylem(psi,psiUpstream, kxylemmax, c, d, true, psiCav=psiCav);
//     if(psi<psiMax) return(psiMax);
//     if(NumericVector::is_na(Eg)) return(NA_REAL);
//   }
//   return((psiPrev+psi)/2.0);
// }






// [[Rcpp::export("hydraulics.E2psiVanGenuchten")]]
double E2psiVanGenuchten(double E, double psiSoil, double krhizomax, double n, double alpha, double psiStep = -0.01, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psi = psiSoil;
  double psiPrev = psi;
  double vgPrev = vanGenuchtenConductance(psi, krhizomax, n, alpha);
  double vg = vgPrev;
  double Eg = 0.0;
  while(Eg<E) {
    psiPrev = psi;
    vgPrev = vg;
    psi = psi + psiStep;
    vg = vanGenuchtenConductance(psi, krhizomax, n, alpha);
    Eg = Eg + ((vg+vgPrev)/2.0)*std::abs(psiStep);
    if(psi<psiMax) return(NA_REAL);
  }
  return(psiPrev);
}


// [[Rcpp::export("hydraulics.E2psiTwoElements")]]
double E2psiTwoElements(double E, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double psiCav = 0.0,
                        double psiStep = -0.001, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psiRoot = E2psiVanGenuchten(E, psiSoil, krhizomax, n, alpha, psiStep, psiMax);
  if(NumericVector::is_na(psiRoot)) return(NA_REAL);
  return(E2psiXylem(E, psiRoot, kxylemmax, c, d, psiCav, psiStep, psiMax));
}


double ludcmp(NumericMatrix a, int n, IntegerVector indx) {
    const double TINY=1.0e-20;
    int i,imax=0,j,k;
    double big,dum,sum,temp;
    NumericVector vv(n);
    double d = 1.0;
    for(i=0;i<n;i++) {
      big = 0.0;
      for(j=0;j<n;j++) if((temp=std::abs(a(i,j)))>big) big=temp;
      if(big==0.0) ::Rf_error("Singular matrix in routine ludcmp");
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
void lubksb(NumericMatrix a, int n, IntegerVector indx, NumericVector b) {
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
// [[Rcpp::export("hydraulics.E2psiRootSystem")]]
List E2psiRootSystem(double E, NumericVector psiSoil, 
                  NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                  NumericVector krootmax, double rootc, double rootd, 
                  NumericVector psiIni = NumericVector::create(0),
                  double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  int nlayers = psiSoil.length();
  //Initialize
  NumericVector x(nlayers+1);
  if(psiIni.size()==(nlayers+1)){
    for(int l=0;l<(nlayers+1);l++) {
      x[l] = psiIni[l];
    }
  } else{
    double minPsi = -0.00001;
    for(int l=0;l<nlayers;l++) {
      x[l] =psiSoil[l];
      minPsi = std::min(minPsi, psiSoil[l]);
      // Rcout<<"("<<x[l]<<") ";
    }
    x[nlayers] = minPsi;
    // Rcout<<"("<<x[nlayers]<<")\n";
    
  }
  
  //Flow across root xylem and rhizosphere elements
  NumericVector Eroot(nlayers), Erhizo(nlayers);
  
  //Newton-Raphson algorithm
  NumericVector p(nlayers+1), fvec(nlayers+1);
  IntegerVector indx(nlayers+1);
  NumericMatrix fjac(nlayers+1,nlayers+1);
  for(int k=0;k<ntrial;k++) {
    // Rcout<<"trial "<<k<<"\n";
    //Calculate steady-state flow functions
    double Esum = 0.0;
    bool stop = false;
    for(int l=0;l<nlayers;l++) {
      Eroot[l] = EXylem(x[nlayers], x[l], krootmax[l], rootc, rootd);
      // Rcout<<"("<<Eroot[l]<<"\n";
      Erhizo[l] = EVanGenuchten(x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l], psiStep, psiTol/1000.0);
      fvec[l] = Erhizo[l] - Eroot[l];
      // Rcout<<" Erhizo"<<l<<": "<< Erhizo[l]<<" Eroot"<<l<<": "<<Eroot[l]<<" fvec: "<<fvec[l]<<"\n";
      // Rcout<<"der psi_l "<<d_psi_l<<"der Eroot psi_l "<<d_Eroot_psi_l<<"  der psi_root "<<d_psi_root<<"\n";
      Esum +=Eroot[l];
    }
    fvec[nlayers] = Esum-E;
    // Rcout<<"fvec_nlayers: "<<fvec[nlayers]<<"\n";
    //Fill Jacobian
    for(int l1=0;l1<nlayers;l1++) { //funcio
      for(int l2=0;l2<nlayers;l2++) { //derivada
        if(l1==l2) fjac(l1,l2) = -vanGenuchtenConductance(x[l2],krhizomax[l2], nsoil[l2], alphasoil[l2])-xylemConductance(x[l2], krootmax[l2], rootc, rootd);
        else fjac(l1,l2) = 0.0;
      }
    }
    fjac(nlayers,nlayers) = 0.0;
    for(int l=0;l<nlayers;l++) { 
      fjac(l,nlayers) = xylemConductance(x[nlayers], krootmax[l], rootc, rootd); //funcio l derivada psi_rootcrown
      fjac(nlayers,l) = xylemConductance(x[l], krootmax[l], rootc, rootd);//funcio nlayers derivada psi_l
      // funcio nlayers derivada psi_rootcrown
      fjac(nlayers,nlayers) +=-xylemConductance(x[nlayers], krootmax[l], rootc, rootd);
    }
    // for(int l1=0;l1<=nlayers;l1++) { //funcio
    //   for(int l2=0;l2<=nlayers;l2++) { //derivada
    //     Rcout<<fjac(l1,l2)<<" ";
    //   }
    //   Rcout<<"\n";
    // }
    //Check function convergence
    double errf = 0.0;
    for(int fi=0;fi<=nlayers;fi++) errf += std::abs(fvec[fi]);
    if(errf<=ETol) break;
    //Right-hand side of linear equations
    for(int fi=0;fi<=nlayers;fi++) p[fi] = -fvec[fi];
    //Solve linear equations using LU decomposition
    ludcmp(fjac,nlayers+1,indx);
    lubksb(fjac,nlayers+1,indx,p);
    //Check root convergence
    double errx = 0.0;
    for(int fi=0;fi<=nlayers;fi++) {
      errx +=std::abs(p[fi]);
      x[fi]+=p[fi];
      x[fi] = std::min(0.0, x[fi]);
      if(x[fi]<psiMax) {
        x[fi] = NA_REAL;
        stop = true;
      }
      // Rcout<<"("<<x[fi]<<") ";
    }
    // Rcout<<"\n";
    if(errx<=psiTol) break;
    else if(k==(ntrial-1)) { //Last trial and no convergence
      for(int fi=0;fi<=nlayers;fi++) x[fi] = NA_REAL;
      Rcout<<"LC";
      stop = true;
    }
    if(stop) break;
  }
  
  //Initialize and copy output
  NumericVector psi(nlayers+1);
  for(int l=0;l<=nlayers;l++) {
    psi[l] = x[l];
    // Rcout<<psi[l]<<" ";
  }
  for(int l=0;l<(nlayers-1);l++) {
    Erhizo[l] = EVanGenuchten(x[l], psiSoil[l], krhizomax[l], nsoil[l], alphasoil[l], psiStep, psiTol/1000.0);
    // Rcout<<Erhizo[l]<<" ";
  }
  return(List::create(Named("Psi") = psi, Named("E")=Erhizo));
} 


// [[Rcpp::export("hydraulics.E2psiNetwork")]]
List E2psiNetwork(double E, NumericVector psiSoil, 
                  NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                  NumericVector krootmax, double rootc, double rootd, 
                  double kstemmax, double stemc, double stemd,
                  double kleafmax, double leafc, double leafd,
                  NumericVector psiIni = NumericVector::create(0),
                  double psiCav = 0.0, //Minimum plant potential (cavitation)
                  double psiStep = -0.001, double psiMax = -10.0, int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  
  List E2psiRS = E2psiRootSystem(E, psiSoil, 
                                 krhizomax, nsoil, alphasoil,
                                 krootmax, rootc, rootd, 
                                 psiIni,
                                 psiStep, psiMax, ntrial, psiTol, ETol);

  NumericVector x = E2psiRS["Psi"];  
  NumericVector Erhizo = E2psiRS["E"];  
  int nlayers = psiSoil.length();

  //Copy output
  NumericVector psi(nlayers+3);
  for(int l=0;l<=nlayers;l++) {
    psi[l] = x[l];
  }
  if(!NumericVector::is_na(x[nlayers])) {
    psi[nlayers+1] = E2psiXylem(E, x[nlayers], kstemmax, stemc, stemd, psiCav, psiStep, psiMax); //Apliquem la fatiga per cavitacio a la caiguda de potencial a la tija 
    psi[nlayers+2] = E2psiXylem(E, psi[nlayers+1], kleafmax, leafc, leafd, 0.0, psiStep, psiMax); 
  } 
  else {
    psi[nlayers+1] = NA_REAL;
    psi[nlayers+2] = NA_REAL;
  } 
  double kterm = xylemConductance(psi[nlayers+2], kleafmax, leafc, leafd);
  return(List::create(Named("Psi") = psi, Named("E")=Erhizo, Named("kterm") = kterm));
} 

/*
 * E -  Flow rate from roots (mmol·s-1·m-2)
 * psiRootCrown - Root crown water potential (MPa)
 * PLC - Proportion of loss conductance [0-1]
 * RWC - Relative water content [0-1]
 * kxylemmax - Whole-stem leaf-specific conductance in mmol·s-1·MPa-1·m-2
 * c, d - parameters of the xylem VC
 * pi0, epsilon - parameters of symplastic xylem tissue 
 * klat - lateral conductance between symplastic and apoplastic tissue (mmol·s-1·MPa-1)
 * ksto - whole-stem vertical symplastic conductance (mmol·s-1·MPa-1)
 * Vmax - Maximum stem volum per leaf area in m ( = m3·m-2)
 * tstep - Time step (s)
 */
// [[Rcpp::export("hydraulics.E2psiXylemCapacitance")]]
List E2psiXylemCapacitance(double E, double psiRootCrown, 
                           NumericVector PLC, NumericVector RWCstorage, 
                           double kxylemmax, double c, double d, 
                           double Vmax, double fapo, double pi0, double epsilon,
                           double klat, double ksto,
                           double tstep = 3600, 
                           double psiStep = -0.0001, double psiMax = -10.0) {
  int n = PLC.size();
  NumericVector newPLC(n, NA_REAL),newRWCstorage(n, NA_REAL), newPsiStem(n, NA_REAL);
  NumericVector Eout(n, NA_REAL),Flat(n, NA_REAL), Fver1(n, NA_REAL), Fver2(n, NA_REAL), V(n, NA_REAL);
  double kxsegmax = kxylemmax*((double) n);
  double kstoseg = ksto*((double) n);
  double Vsegmax = Vmax/((double) n);
  double Vprev;
  double psiPLC, psiStorage;
  //Initial values from root
  double Ein = E;
  double Einc = Ein;
  double psiUp = psiRootCrown;
  double psiStorageUp = psiRootCrown;
  double m3tommol = 55555556.0;
  for(int i=0;i<n;i++) {
    //Transform PLC and RWC into water potential of storage compartments
    psiPLC = d*pow(-1.0*log(1.0-PLC[i]),1.0/c);
    psiStorage = -pi0 -epsilon*(1.0 - RWCstorage[i]);
    
    //Lateral flow
    Flat[i] =  klat*(psiStorage-psiPLC);
    
    //Symplastic flow from below
    Fver1[i] = kstoseg*(psiStorageUp-psiStorage);
    
    if(i==0) {
      if(Ein<Fver1[i]) Fver1[i] = Ein; //Do not substract more than input flow
      Ein = Ein - Fver1[i]; //Substract or add flow from/to root input
      Einc = Ein;
    }
    
    //Store psi value before adding lateral flow
    newPsiStem[i] = E2psiXylem(Ein, psiUp, kxsegmax, c,d, psiPLC, psiStep, psiMax);
    
    //Upwards symplastic flow
    Fver2[i] = kstoseg*(psiStorage-newPsiStem[i]);
    
    //Increase in flow due to new cavitation
    Vprev = Vsegmax*fapo*(1.0-PLC[i]);
    V[i] = std::min(Vprev, Vsegmax*fapo*exp(-pow(newPsiStem[i]/d,c))); //Only allow decreases in volume (i.e. refilling cannot occur unless there is lateral flow)
    Eout[i] = Ein - (m3tommol/tstep)*(V[i]-Vprev);
    
    if(i==(n-1)) {
      Eout[i] = Eout[i] - Fver2[i];
    }
    
    //Update variables for next segment
    psiUp = newPsiStem[i];
    Ein = Eout[i];
    psiStorageUp = psiStorage;
  }
  //Update compartments
  for(int i=0;i<n;i++) {
    //Apoplastic compartment
    V[i] = V[i] + (tstep/m3tommol)*Flat[i];
    if(V[i]>Vsegmax*fapo) { //Correct if flow is too high
      double Vdif = Vsegmax*fapo - V[i];
      V[i] = Vsegmax*fapo;
      Flat[i] = Flat[i]-(Vdif*(m3tommol/tstep));
    }
    newPLC[i] = 1.0- (V[i]/(Vsegmax*fapo));
    //Symplastic compartment
    newRWCstorage[i] = (Vsegmax*(1.0-fapo)*RWCstorage[i] + (tstep/m3tommol)*(Fver1[i] - Fver2[i] - Flat[i]))/(Vsegmax*(1.0-fapo));
  }
  
  return(List::create( _["Einc"] = Einc, _["Eout"] = Eout,_["newPsiStem"] = newPsiStem, _["newPLC"] = newPLC, _["newRWCstorage"] = newRWCstorage));
}


/*
 * E -  Flow rate from roots (mmol·s-1·m-2)
 * psiRootCrown - Root crown water potential (MPa)
 * PLC - Proportion of loss conductance [0-1]
 * RWC - Relative water content [0-1]
 * kstemmax - Whole-stem leaf-specific conductance in mmol·s-1·MPa-1·m-2
 * stemc, stemd - parameters of the stem VC
 * kleafmax - Whole-stem leaf-specific conductance in mmol·s-1·MPa-1·m-2
 * leafc, leafd - parameters of the leaf VC
 * pi0, epsilon - parameters of symplastic xylem tissue 
 * klat - lateral conductance between symplastic and apoplastic tissue (mmol·s-1·MPa-1)
 * ksto - whole-stem vertical symplastic conductance (mmol·s-1·MPa-1)
 * Vmax - Maximum stem volum per leaf area in m ( = m3·m-2)
 * tstep - Time step (s)
 */
// [[Rcpp::export("hydraulics.E2psiAboveground")]]
List E2psiAboveGround(double E, double psiRootCrown,                           
                      NumericVector PLC, NumericVector RWCstorage, 
                      double kstemmax, double stemc, double stemd,
                      double kleafmax, double leafc, double leafd,
                      double Vmax, double fapo, double pi0, double epsilon,
                      double klat, double ksto,
                      double tstep = 3600, 
                      double psiStep = -0.0001, double psiMax = -10.0) {
  
  List E2psiCap = E2psiXylemCapacitance(E, psiRootCrown, 
                                        PLC, RWCstorage, 
                                        kstemmax, stemc,  stemd, 
                                        Vmax, fapo, pi0, epsilon,
                                        klat, ksto,
                                        tstep, 
                                        psiStep, psiMax);
  
  NumericVector newPsiStem = E2psiCap["newPsiStem"];
  NumericVector Eout = E2psiCap["Eout"];
  int n = newPsiStem.size();
  double psiLeaf = E2psiXylem(Eout[n-1], newPsiStem[n-1], kleafmax, leafc, leafd, 0.0, psiStep, psiMax); 
  double kterm = xylemConductance(psiLeaf, kleafmax, leafc, leafd);
  return(List::create( _["Einc"] = E2psiCap["Einc"], 
                       _["Eout"] = Eout,
                       _["E"] = Eout[n-1],
                       _["psiLeaf"] = psiLeaf,
                       _["newPsiStem"] = newPsiStem, 
                       _["newPLC"] = E2psiCap["newPLC"], 
                       _["newRWCstorage"] = E2psiCap["newRWCstorage"],
                       _["kterm"] = kterm));
}

// [[Rcpp::export("hydraulics.supplyFunctionOneXylem")]]
List supplyFunctionOneXylem(NumericVector psiSoil, NumericVector v,
                            double kstemmax, double stemc, double stemd, double psiCav = 0.0,
                            int maxNsteps=200, double psiStep = -0.001, double psiMax = -10.0, double dE=0.01) {
  int nlayers = psiSoil.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyElayers(maxNsteps,nlayers);
  NumericVector supplyPsi(maxNsteps);
  
  supplyE[0] = 0;
  for(int l=0;l<nlayers;l++) supplyElayers[l] = 0.0;
  supplyPsi[0] = averagePsi(psiSoil, v, stemc, stemd);
  NumericVector Psilayers(nlayers);
  //Calculate initial slope
  for(int l=0;l<nlayers;l++) {
    Psilayers[l] = E2psiXylem(dE, psiSoil[l], 
                              kstemmax, stemc,stemd, psiCav, 
                              psiStep, psiMax);
    // Rcout<<Psilayers[l]<<" ";
  }
  Rcout<<"\n";
  double psiI = averagePsi(Psilayers, v, stemc, stemd);
  double maxdEdp = dE/std::abs(psiI-supplyPsi[0]);
  // Rcout<<maxdEdp<<"\n";
  
  int nsteps = 1;
  dE = maxdEdp*0.1;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    // Rcout<<supplyE[i]<<" ";
    for(int l=0;l<nlayers;l++) {
      Psilayers[l] = E2psiXylem(supplyE[i], psiSoil[l],
                                kstemmax, stemc,stemd, psiCav, 
                                psiStep, psiMax);
      // Rcout<<Psilayers[l]<<" ";
    }
    // Rcout<<"\n";
    supplyPsi[i] = averagePsi(Psilayers, v, stemc, stemd);
    for(int l=0;l<nlayers;l++) {
      supplyElayers(i,l) = supplyE[i]*v[l];
    }
    
    if(!NumericVector::is_na(supplyPsi[i])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsi[1]-supplyPsi[0]);
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsi[i-1]-supplyPsi[i-2]);
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsi[i]-supplyPsi[i-1]);
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      dE = supplydEdp[i-1]*0.1;
      nsteps++;
      if(supplydEdp[i-1]<0.01*maxdEdp) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsi[nsteps-1]-supplyPsi[nsteps-2]);
  //Copy values tp nsteps
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyElayersDef(nsteps,nlayers);
  NumericVector supplyPsiPlant(nsteps);
  for(int i=0;i<nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyPsiPlant[i] = supplyPsi[i];
    for(int l=0;l<nlayers;l++) {
      supplyElayersDef(i,l) = supplyElayers(i,l);
    }
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("Elayers") = supplyElayersDef,
                      Named("PsiPlant")=supplyPsiPlant,
                      Named("dEdP")=supplydEdpDef));
  
}

// [[Rcpp::export("hydraulics.supplyFunctionTwoElements")]]
List supplyFunctionTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, 
                               double psiCav = 0.0, 
                               double dE = 0.1, double psiMax = -10.0) {
  dE = std::min(dE,Emax/5.0);
  int maxNsteps = round(Emax/dE)+1;
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyFittedE(maxNsteps);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiPlant(maxNsteps);
  double Eg1 = 0.0;
  double Eg2 = 0.0;
  double psiStep1 = -0.1;
  double psiStep2 = -0.1;
  double psiRoot = psiSoil;
  double psiPlant = psiSoil;
  double vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
  double vg = 0.0;
  double wPrev = xylemConductance(std::min(psiCav,psiPlant), kxylemmax, c, d); //conductance can decrease if psiCav < psiPlant
  double w = 0.0;
  double incr = 0.0;
  supplyPsiRoot[0] = psiSoil;
  supplyPsiPlant[0] = psiSoil;
  supplyE[0] = 0.0;
  double psiPrec = -0.000001;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    psiStep1 = -0.01;
    psiRoot = supplyPsiRoot[i-1];
    vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
    while((psiStep1<psiPrec) & (psiRoot>psiMax))  {
      vg = vanGenuchtenConductance(psiRoot+psiStep1, krhizomax, n, alpha);
      incr = ((vg+vgPrev)/2.0)*std::abs(psiStep1);
      if((Eg1+incr)>supplyE[i]) {
        psiStep1 = psiStep1*0.5;
      } else {
        psiRoot = psiRoot + psiStep1;
        Eg1 = Eg1+incr;
        vgPrev = vg;
      }
    }
    supplyPsiRoot[i] = psiRoot;
    if(supplyPsiRoot[i]<psiMax) supplyPsiRoot[i] = psiMax;
    
    psiStep2 = -0.01;
    Eg2 = 0.0;
    psiPlant = psiRoot;
    wPrev = xylemConductance(std::min(psiCav,psiPlant), kxylemmax, c, d);
    while((psiStep2<psiPrec) & (psiPlant>psiMax))  {
      w = xylemConductance(std::min(psiCav,psiPlant+psiStep2), kxylemmax, c, d);
      incr = ((w+wPrev)/2.0)*std::abs(psiStep2);
      if((Eg2+incr)>supplyE[i]) {
        psiStep2 = psiStep2*0.5;
      } else {
        psiPlant = psiPlant + psiStep2;
        Eg2 = Eg2+incr;
        wPrev = w;
      }
    }
    supplyPsiPlant[i] = psiPlant;
    if(supplyPsiPlant[i]<psiMax) supplyPsiPlant[i] = psiMax;
    supplyFittedE[i] = std::max(supplyFittedE[i-1], Eg2); //Ensure non-decreasing function
    // supplyFittedE[i] = Eg2;
    if(i==1) {
      supplydEdp[0] = (supplyFittedE[1]-supplyFittedE[0])/(supplyPsiPlant[0]-supplyPsiPlant[1]); 
      if((supplyPsiPlant[0]-supplyPsiPlant[1])==0.0) supplydEdp[0] = 0.0;
    }
    else if(i>1) {
      supplydEdp[i-1] = 0.5*(supplyFittedE[i-1]-supplyFittedE[i-2])/(supplyPsiPlant[i-2]-supplyPsiPlant[i-1])+0.5*(supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-2]-supplyPsiPlant[i-1])==0.0) supplydEdp[i-1] = 0.0;
      else if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i-1] = 0.0;
    }
    if(i==(maxNsteps-1)) {
      supplydEdp[i] = (supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i] = 0.0;
    }
  }
  return(List::create(Named("E") = supplyE,
                      Named("FittedE") = supplyFittedE,
                      Named("PsiRoot")=supplyPsiRoot, 
                      Named("PsiPlant")=supplyPsiPlant,
                      Named("dEdP")=supplydEdp));
}

// [[Rcpp::export("hydraulics.supplyFunctionThreeElements")]]
List supplyFunctionThreeElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double kleafmax, 
                                 double n, double alpha, 
                                 double stemc, double stemd, 
                                 double leafc, double leafd,
                                 double psiCav = 0.0, 
                                 double dE = 0.1, double psiMax = -10.0) {
  dE = std::min(dE,Emax/5.0);
  int maxNsteps = round(Emax/dE)+1;
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyFittedE(maxNsteps);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiStem(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  double Eg1 = 0.0;
  double Eg2 = 0.0;
  double Eg3 = 0.0;
  double psiStep1 = -0.1;
  double psiStep2 = -0.1;
  double psiStep3 = -0.1;
  double psiRoot = psiSoil;
  double psiStem = psiSoil;
  double psiLeaf = psiSoil;
  double vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
  double vg = 0.0;
  //conductance can decrease if psiCav < psiStem/psiLeaf
  double wPrevStem = xylemConductance(std::min(psiCav,psiStem), kxylemmax, stemc, stemd); 
  double wPrevLeaf = xylemConductance(std::min(psiCav,psiLeaf), kleafmax, leafc, leafd); 
  double w = 0.0;
  double incr = 0.0;
  supplyPsiRoot[0] = psiSoil;
  supplyPsiStem[0] = psiSoil;
  supplyPsiLeaf[0] = psiSoil;
  supplyE[0] = 0.0;
  double psiPrec = -0.000001;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    psiStep1 = -0.01;
    
    // Root
    psiRoot = supplyPsiRoot[i-1];
    vgPrev = vanGenuchtenConductance(psiRoot, krhizomax, n, alpha);
    while((psiStep1<psiPrec) & (psiRoot>psiMax))  {
      vg = vanGenuchtenConductance(psiRoot+psiStep1, krhizomax, n, alpha);
      incr = ((vg+vgPrev)/2.0)*std::abs(psiStep1);
      if((Eg1+incr)>supplyE[i]) {
        psiStep1 = psiStep1*0.5;
      } else {
        psiRoot = psiRoot + psiStep1;
        Eg1 = Eg1+incr;
        vgPrev = vg;
      }
    }
    supplyPsiRoot[i] = psiRoot;
    if(supplyPsiRoot[i]<psiMax) supplyPsiRoot[i] = psiMax;
    
    //Stem
    psiStep2 = -0.01;
    Eg2 = 0.0;
    psiStem = psiRoot;
    wPrevStem = xylemConductance(std::min(psiCav,psiStem), kxylemmax, stemc, stemd);
    while((psiStep2<psiPrec) & (psiStem>psiMax))  {
      w = xylemConductance(std::min(psiCav,psiStem+psiStep2), kxylemmax, stemc, stemd);
      incr = ((w+wPrevStem)/2.0)*std::abs(psiStep2);
      if((Eg2+incr)>supplyE[i]) {
        psiStep2 = psiStep2*0.5;
      } else {
        psiStem = psiStem + psiStep2;
        Eg2 = Eg2+incr;
        wPrevStem = w;
      }
    }
    supplyPsiStem[i] = psiStem;
    if(supplyPsiStem[i]<psiMax) supplyPsiStem[i] = psiMax;
    
    //Leaf
    psiStep3 = -0.01;
    Eg3 = 0.0;
    psiLeaf = psiStem;
    wPrevLeaf = xylemConductance(psiLeaf, kleafmax,  leafc, leafd);
    while((psiStep3<psiPrec) & (psiLeaf>psiMax))  {
      w = xylemConductance(psiLeaf+psiStep3, kleafmax,  leafc, leafd);
      incr = ((w+wPrevLeaf)/2.0)*std::abs(psiStep3);
      if((Eg3+incr)>supplyE[i]) {
        psiStep3 = psiStep3*0.5;
      } else {
        psiLeaf = psiLeaf + psiStep3;
        Eg3 = Eg3+incr;
        wPrevLeaf = w;
      }
    }
    supplyPsiLeaf[i] = psiLeaf;
    if(supplyPsiLeaf[i]<psiMax) supplyPsiLeaf[i] = psiMax;
    
    
    //Ensure non-decreasing function
    supplyFittedE[i] = std::max(supplyFittedE[i-1], Eg3); 
    
    
    // supplyFittedE[i] = Eg2;
    if(i==1) {
      supplydEdp[0] = (supplyFittedE[1]-supplyFittedE[0])/(supplyPsiLeaf[0]-supplyPsiLeaf[1]); 
      if((supplyPsiLeaf[0]-supplyPsiLeaf[1])==0.0) supplydEdp[0] = 0.0;
    }
    else if(i>1) {
      supplydEdp[i-1] = 0.5*(supplyFittedE[i-1]-supplyFittedE[i-2])/(supplyPsiLeaf[i-2]-supplyPsiLeaf[i-1])+0.5*(supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiLeaf[i-1]-supplyPsiLeaf[i]); 
      if((supplyPsiLeaf[i-2]-supplyPsiLeaf[i-1])==0.0) supplydEdp[i-1] = 0.0;
      else if((supplyPsiLeaf[i-1]-supplyPsiLeaf[i])==0.0) supplydEdp[i-1] = 0.0;
    }
    if(i==(maxNsteps-1)) {
      supplydEdp[i] = (supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiLeaf[i-1]-supplyPsiLeaf[i]); 
      if((supplyPsiLeaf[i-1]-supplyPsiLeaf[i])==0.0) supplydEdp[i] = 0.0;
    }
  }
  return(List::create(Named("E") = supplyE,
                      Named("FittedE") = supplyFittedE,
                      Named("PsiRoot")=supplyPsiRoot, 
                      Named("PsiStem")=supplyPsiStem,
                      Named("PsiLeaf")=supplyPsiLeaf,
                      Named("dEdP")=supplydEdp));
}

// [[Rcpp::export("hydraulics.supplyFunctionRootSystem")]]
List supplyFunctionRootSystem(NumericVector psiSoil, 
                           NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                           NumericVector krootmax, double rootc, double rootd, 
                           double minFlow = 0.0, int maxNsteps=400, double psiStep = -0.001, double psiMax = -10.0, 
                           int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  int nlayers = psiSoil.size();
  int nnodes = nlayers+1; //Rhizosphere+rootcollar
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyElayers(maxNsteps,nlayers);
  NumericMatrix supplyPsi(maxNsteps,nnodes);

  List sol = E2psiRootSystem(minFlow, psiSoil, 
                          krhizomax,  nsoil,  alphasoil,
                          krootmax,  rootc,  rootd, 
                          NumericVector::create(0),
                          psiStep, psiMax, ntrial,psiTol, ETol);
  NumericVector solE = sol["E"];
  NumericVector solPsi = sol["Psi"];
  for(int l=0;l<nlayers;l++) supplyElayers(0,l) = solE[l];
  for(int n=0;n<nnodes;n++) supplyPsi(0,n) = solPsi[n];
  supplyE[0] = minFlow;
  
  //Calculate initial slope
  List solI = E2psiRootSystem(minFlow+ETol*2.0, psiSoil, 
                           krhizomax,  nsoil,  alphasoil,
                           krootmax,  rootc,  rootd, 
                           solPsi,
                           psiStep, psiMax, ntrial,psiTol, ETol);
  NumericVector solPsiI = solI["Psi"];
  double maxdEdp = (ETol*2.0)/std::abs(solPsiI[nnodes-1]-solPsi[nnodes-1]);
  
  int nsteps = 1;
  double dE = std::min(0.05,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE; 
    sol = E2psiRootSystem(supplyE[i], psiSoil, 
                       krhizomax,  nsoil,  alphasoil,
                       krootmax,  rootc,  rootd, 
                       solPsi,
                       psiStep, psiMax, ntrial,psiTol, ETol);
    solE = sol["E"];
    solPsi = sol["Psi"];
    // Rcout<<supplyE[i]<<" ";
    for(int l=0;l<nlayers;l++) {
      supplyElayers(i,l) = solE[l];
      // Rcout<<supplyElayers(i,l)<<" ";
    }
    for(int n=0;n<nnodes;n++) {
      supplyPsi(i,n) = solPsi[n];
      // Rcout<<supplyPsi(i,n)<<" ";
    }
    // Rcout<<"\n";
    
    if(!NumericVector::is_na(solPsi[nnodes-1])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsi(1,nnodes-1)-supplyPsi(0,nnodes-1));
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsi(i-1,nnodes-1)-supplyPsi(i-2,nnodes-1));
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsi(i,nnodes-1)-supplyPsi(i-1,nnodes-1));
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      dE = std::min(0.05,supplydEdp[i-1]*0.05);
      nsteps++; 
      if(supplydEdp[i-1]<0.01*maxdEdp) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsi(nsteps-1,nnodes-1)-supplyPsi(nsteps-2,nnodes-1));
  //Copy values tp nsteps
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyElayersDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizo(nsteps,nlayers);
  NumericVector supplyPsiStem(nsteps);
  NumericVector supplyPsiLeaf(nsteps);
  NumericVector supplyPsiRoot(nsteps);
  for(int i=0;i<nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyPsiRoot[i] = supplyPsi(i,nnodes-1);
    for(int l=0;l<nlayers;l++) {
      supplyElayersDef(i,l) = supplyElayers(i,l);
      supplyPsiRhizo(i,l) = supplyPsi(i,l);
    }
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("Elayers") = supplyElayersDef,
                      Named("PsiRhizo")=supplyPsiRhizo, 
                      Named("PsiRoot")=supplyPsiRoot,
                      Named("dEdP")=supplydEdpDef));
  
}


// [[Rcpp::export("hydraulics.supplyFunctionNetwork")]]
List supplyFunctionNetwork(NumericVector psiSoil, 
                           NumericVector krhizomax, NumericVector nsoil, NumericVector alphasoil,
                           NumericVector krootmax, double rootc, double rootd, 
                           double kstemmax, double stemc, double stemd,
                           double kleafmax, double leafc, double leafd,
                           double psiCav = 0.0,
                           double minFlow = 0.0, int maxNsteps=400, double psiStep = -0.001, double psiMax = -10.0, 
                           int ntrial = 10, double psiTol = 0.0001, double ETol = 0.0001) {
  int nlayers = psiSoil.size();
  int nnodes = nlayers+3; //Rhizosphere+rootcollar+stem+leaf
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericMatrix supplyElayers(maxNsteps,nlayers);
  NumericMatrix supplyPsi(maxNsteps,nnodes);
  NumericMatrix supplyKterm(maxNsteps);

  List sol = E2psiNetwork(minFlow, psiSoil, 
                          krhizomax,  nsoil,  alphasoil,
                          krootmax,  rootc,  rootd, 
                          kstemmax,  stemc,  stemd,
                          kleafmax,  leafc,  leafd,
                          NumericVector::create(0),
                          psiCav,
                          psiStep, psiMax, ntrial,psiTol, ETol);
  NumericVector solE = sol["E"];
  NumericVector solPsi = sol["Psi"];
  for(int l=0;l<nlayers;l++) supplyElayers(0,l) = solE[l];
  for(int n=0;n<nnodes;n++) supplyPsi(0,n) = solPsi[n];
  supplyKterm[0] = sol["kterm"];
  supplyE[0] = minFlow;
  
  //Calculate initial slope
  List solI = E2psiNetwork(minFlow+ETol*2.0, psiSoil, 
                           krhizomax,  nsoil,  alphasoil,
                           krootmax,  rootc,  rootd, 
                           kstemmax,  stemc,  stemd,
                           kleafmax,  leafc,  leafd,
                           solPsi,
                           psiCav,
                           psiStep, psiMax, ntrial,psiTol, ETol);
  NumericVector solPsiI = solI["Psi"];
  double maxdEdp = (ETol*2.0)/std::abs(solPsiI[nnodes-1]-solPsi[nnodes-1]);
  
  int nsteps = 1;
  double dE = std::min(0.05,maxdEdp*0.05);
  for(int i=1;i<maxNsteps;i++) {
    // if(i==3) stop("kk");
    supplyE[i] = supplyE[i-1]+dE; 
    sol = E2psiNetwork(supplyE[i], psiSoil, 
                            krhizomax,  nsoil,  alphasoil,
                            krootmax,  rootc,  rootd, 
                            kstemmax,  stemc,  stemd,
                            kleafmax,  leafc,  leafd,
                            solPsi,
                            psiCav,
                            psiStep, psiMax, ntrial,psiTol, ETol);
    solE = sol["E"];
    solPsi = sol["Psi"];
    supplyKterm[i] = sol["kterm"];
    // Rcout<<supplyE[i]<<" ";
    for(int l=0;l<nlayers;l++) {
      supplyElayers(i,l) = solE[l];
      // Rcout<<supplyElayers(i,l)<<" ";
    }
    for(int n=0;n<nnodes;n++) {
      supplyPsi(i,n) = solPsi[n];
      // Rcout<<supplyPsi(i,n)<<" ";
    }
    // Rcout<<"\n";
    
    if(!NumericVector::is_na(solPsi[nnodes-1])) {
      if(i==1) {
        supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsi(1,nnodes-1)-supplyPsi(0,nnodes-1));
      } else {
        double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsi(i-1,nnodes-1)-supplyPsi(i-2,nnodes-1));
        double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsi(i,nnodes-1)-supplyPsi(i-1,nnodes-1));
        supplydEdp[i-1] = (d1+d2)/2.0;
      }
      dE = std::min(0.05,supplydEdp[i-1]*0.05);
      nsteps++; 
      if(supplydEdp[i-1]<0.01*maxdEdp) break;
    } else {
      break;
    }
  }
  //Calculate last dEdP
  if(nsteps>1) supplydEdp[nsteps-1] = (supplyE[nsteps-1]-supplyE[nsteps-2])/std::abs(supplyPsi(nsteps-1,nnodes-1)-supplyPsi(nsteps-2,nnodes-1));
  //Copy values tp nsteps
  NumericVector supplyKtermDef(nsteps);
  NumericVector supplyEDef(nsteps);
  NumericVector supplydEdpDef(nsteps);
  NumericMatrix supplyElayersDef(nsteps,nlayers);
  NumericMatrix supplyPsiRhizo(nsteps,nlayers);
  NumericVector supplyPsiStem(nsteps);
  NumericVector supplyPsiLeaf(nsteps);
  NumericVector supplyPsiRoot(nsteps);
  for(int i=0;i<nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiRoot[i] = supplyPsi(i,nnodes-3);
    supplyPsiStem[i] = supplyPsi(i,nnodes-2);
    supplyPsiLeaf[i] = supplyPsi(i,nnodes-1);
    for(int l=0;l<nlayers;l++) {
      supplyElayersDef(i,l) = supplyElayers(i,l);
      supplyPsiRhizo(i,l) = supplyPsi(i,l);
    }
  }
  return(List::create(Named("E") = supplyEDef,
                      Named("Elayers") = supplyElayersDef,
                      Named("PsiRhizo")=supplyPsiRhizo, 
                      Named("PsiRoot")=supplyPsiRoot,
                      Named("PsiStem")=supplyPsiStem,
                      Named("PsiLeaf")=supplyPsiLeaf,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}


// [[Rcpp::export("hydraulics.supplyFunctionAboveGround")]]
List supplyFunctionAboveground(NumericVector Erootcrown, NumericVector psiRootcrown, 
                               NumericVector PLC, NumericVector RWCstorage, 
                               double kstemmax, double stemc, double stemd,
                               double kleafmax, double leafc, double leafd,
                               double Vmax, double fapo, double pi0, double epsilon,
                               double klat, double ksto,
                               double tstep = 3600, 
                               double psiStep = -0.001, double psiMax = -10.0) {
  int nnodes = PLC.size(); // stem nodes + leaf
  int maxNsteps = Erootcrown.size();
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyPsiLeaf(maxNsteps);
  NumericVector supplyKterm(maxNsteps);
  NumericMatrix supplyPsiStem(maxNsteps,nnodes);
  NumericMatrix supplyPLC(maxNsteps,nnodes);
  NumericMatrix supplyRWCstorage(maxNsteps,nnodes);
  

 
  int Nsteps = 0;
  for(int i=0;i<maxNsteps;i++) {
    List sol = E2psiAboveGround(Erootcrown[i], psiRootcrown[i],                           
                                PLC, RWCstorage, 
                                kstemmax, stemc, stemd,
                                kleafmax, leafc, leafd,
                                Vmax, fapo, pi0, epsilon,
                                klat, ksto,
                                tstep, 
                                psiStep, psiMax);
    NumericVector solNewPsiStem = sol["newPsiStem"];
    NumericVector solNewPLC = sol["newPLC"];
    NumericVector solNewRWCstorage = sol["newRWCstorage"];
    for(int n=0;n<nnodes;n++) {
      supplyPsiStem(i,n) = solNewPsiStem[n]; 
      supplyPLC(i,n) = solNewPLC[n]; 
      supplyRWCstorage(i,n) = solNewRWCstorage[n]; 
    }
    supplyPsiLeaf[i] = sol["psiLeaf"];
    if(NumericVector::is_na(supplyPsiLeaf[i])) {
      break; 
    } else {
      Nsteps = Nsteps + 1;
    }
    supplyKterm[i] = sol["kterm"];
    supplyE[i] = sol["E"];
    
    if(i==1) {
      supplydEdp[0] = (supplyE[1]-supplyE[0])/std::abs(supplyPsiLeaf[1]-supplyPsiLeaf[0]);
    } else {
      double d1 = (supplyE[i-1]-supplyE[i-2])/std::abs(supplyPsiLeaf[i-1]-supplyPsiLeaf[i-2]);
      double d2 = (supplyE[i]-supplyE[i-1])/std::abs(supplyPsiLeaf[i]-supplyPsiLeaf[i-1]);
      supplydEdp[i-1] = (d1+d2)/2.0;
    }
    
  }
  //Calculate last dEdP
  if(Nsteps>1) supplydEdp[Nsteps-1] = (supplyE[Nsteps-1]-supplyE[Nsteps-2])/std::abs(supplyPsiLeaf[Nsteps-1]-supplyPsiLeaf[Nsteps-2]);
  
  //Copy values tp nsteps
  NumericVector supplyKtermDef(Nsteps);
  NumericVector supplyEDef(Nsteps);
  NumericVector supplydEdpDef(Nsteps);
  NumericVector supplyPsiLeafDef(Nsteps);
  NumericMatrix supplyPsiStemDef(Nsteps,nnodes);
  NumericMatrix supplyPLCDef(Nsteps,nnodes);
  NumericMatrix supplyRWCstorageDef(Nsteps,nnodes);
  
  for(int i=0;i<Nsteps;i++) {
    supplyEDef[i] = supplyE[i];
    supplydEdpDef[i] = supplydEdp[i];
    supplyKtermDef[i] = supplyKterm[i];
    supplyPsiLeafDef[i] = supplyPsiLeaf[i];
    for(int n=0;n<nnodes;n++) {
      supplyPsiStemDef(i,n) = supplyPsiStem(i,n); 
      supplyPLCDef(i,n) = supplyPLC(i,n); 
      supplyRWCstorageDef(i,n) = supplyRWCstorage(i,n); 
    }
  }
  
  return(List::create(Named("E") = supplyEDef,
                      Named("PLC")=supplyPLCDef,
                      Named("RWC")=supplyRWCstorageDef,
                      Named("PsiStem")=supplyPsiStemDef,
                      Named("PsiLeaf")=supplyPsiLeafDef,
                      Named("dEdP")=supplydEdpDef,
                      Named("kterm") = supplyKtermDef));
  
}


// [[Rcpp::export("hydraulics.regulatedPsiXylem")]]
NumericVector regulatedPsiXylem(double E, double psiUpstream, double kxylemmax, double c, double d, double psiStep = -0.01) {
  //If Ein > Ecrit then set Ein to Ecrit
  double psiUnregulated = E2psiXylem(E, psiUpstream, kxylemmax, c, d, 0.0, psiStep);
  double Ec = Ecrit(psiUpstream, kxylemmax,c,d);
  double Ein = E;
  if(Ein > Ec) {
    Ein = Ec;
    psiUnregulated = psiCrit(c,d);
  }
  double deltaPsiUnregulated = psiUnregulated - psiUpstream;
  double kp = xylemConductance(psiUpstream, kxylemmax, c, d);
  double deltaPsiRegulated = deltaPsiUnregulated*(xylemConductance(psiUnregulated, kxylemmax, c, d)/kp);
  //replace by maximum if found for lower psi values
  // Rcout <<"Initial "<<psiUnregulated << " "<< deltaPsiRegulated <<"\n";
  for(double psi = psiUpstream; psi > psiUnregulated; psi +=psiStep) {
    double deltaPsi = (psi-psiUpstream)*(xylemConductance(psi, kxylemmax, c,d)/kp);
    // Rcout <<psi << " "<< deltaPsi<< " "<< deltaPsiRegulated <<"\n";
    if(NumericVector::is_na(deltaPsiRegulated)) deltaPsiRegulated = deltaPsi;
    else if(deltaPsi < deltaPsiRegulated) deltaPsiRegulated = deltaPsi;
  }
  double psiRegulated = psiUpstream + deltaPsiRegulated;
  double Efin = EXylem(psiRegulated, psiUpstream, kxylemmax, c, d);
  double relativeConductance1 = Efin/Ein;
  double relativeConductance2 = Efin/E;
  return(NumericVector::create(psiUnregulated, psiRegulated, Ein, Efin, relativeConductance1, relativeConductance2));
}

// [[Rcpp::export("hydraulics.regulatedPsiTwoElements")]]
NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double krhizomax, double kxylemmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0) {
  List s = supplyFunctionTwoElements(Emax, psiSoil, krhizomax, kxylemmax, n, alpha, c, d, 0.0, dE,psiMax);
  NumericVector supplyPsi = s["PsiPlant"];
  NumericVector Efitted = s["FittedE"];
  NumericVector dEdP = s["dEdP"];
  int maxNsteps = Efitted.size();
  double deltaPsiRegulated = 0.0;
  double deltaPsiRegulatedi=0.0;
  double dEdP0 = dEdP[0];
  for(int i=1;i<maxNsteps;i++) {
    if(supplyPsi[i]>psiMax) {
      deltaPsiRegulatedi = (supplyPsi[i] - psiSoil)*std::min(1.0, dEdP[i]/dEdP0);
      // Rcout<<supplydEdp <<" "<<deltaPsiRegulatedi<<"\n";
      if(deltaPsiRegulatedi < deltaPsiRegulated) {
        deltaPsiRegulated = deltaPsiRegulatedi;
      }
    }
  }
  //Regulated potential
  double psiRegulated = psiSoil + deltaPsiRegulated;
  //Find transpiration corresponding to regulated potential
  double ERegulated = 0.0, dEdPRegulated = 0.0;
  for(int i=1;i<maxNsteps;i++) {
    if((supplyPsi[i-1] >= psiRegulated) & (supplyPsi[i]<psiRegulated)) {
      ERegulated = Efitted[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + Efitted[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      ERegulated = std::min(ERegulated, Emax);
      psiRegulated = supplyPsi[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + supplyPsi[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      dEdPRegulated = dEdP[i]*std::abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + dEdP[i-1]*std::abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      if((supplyPsi[i-1]-supplyPsi[i])==0.0) dEdPRegulated = 0.0;
      // Rcout<<dEdP[i]<< " "<<dEdP[i-1]<< " "<<dEdPRegulated<<"\n";
      break;
    }
  }
  return(NumericVector::create(supplyPsi[maxNsteps-1], psiRegulated, Efitted[maxNsteps-1], ERegulated, dEdPRegulated));
}



/*
 * Parametrization of rhizosphere conductance
 */

double rhizosphereResistancePercent(double psiSoil, 
                                    double krhizomax, double n, double alpha,
                                    double krootmax, double rootc, double rootd,
                                    double kstemmax, double stemc, double stemd,
                                    double kleafmax, double leafc, double leafd) {
  double krhizo = vanGenuchtenConductance(psiSoil, krhizomax, n, alpha);
  double kroot = xylemConductance(psiSoil, krootmax, rootc, rootd);
  double kstem = xylemConductance(psiSoil, kstemmax, stemc, stemd);
  double kleaf = xylemConductance(psiSoil, kleafmax, leafc, leafd);
  return(100.0*(1.0/krhizo)/((1.0/kroot)+(1.0/kstem)+(1.0/kleaf)+(1.0/krhizo)));
}

// [[Rcpp::export("hydraulics.averageRhizosphereResistancePercent")]]
double averageRhizosphereResistancePercent(double krhizomax, double n, double alpha,
                                           double krootmax, double rootc, double rootd,
                                           double kstemmax, double stemc, double stemd, 
                                           double kleafmax, double leafc, double leafd,
                                           double psiStep = -0.01){
  double psiC = psiCrit(stemc, stemd);
  double cnt = 0.0;
  double sum = 0.0;
  for(double psi=0.0; psi>psiC;psi += psiStep) {
    sum +=rhizosphereResistancePercent(psi, krhizomax, n,alpha,krootmax, rootc, rootd,
                                       kstemmax, stemc,stemd,
                                       kleafmax, leafc,leafd);
    cnt+=1.0;
  }
  return(sum/cnt);
}

// [[Rcpp::export("hydraulics.findRhizosphereMaximumConductance")]]
double findRhizosphereMaximumConductance(double averageResistancePercent, double n, double alpha,
                                         double krootmax, double rootc, double rootd,
                                         double kstemmax, double stemc, double stemd,
                                         double kleafmax, double leafc, double leafd) {
  double step = 1.0;
  double fTol = 0.1;
  double krhizomaxlog = 0.0;
  double f = averageRhizosphereResistancePercent(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
                                                 kstemmax, stemc,stemd,
                                                 kleafmax, leafc,leafd);
  while(std::abs(f-averageResistancePercent)>fTol) {
    if(f>averageResistancePercent) {
      krhizomaxlog += step; 
    } else {
      krhizomaxlog -= step;
      step = step/2.0;
    }
    f = averageRhizosphereResistancePercent(exp(krhizomaxlog), n,alpha,krootmax, rootc, rootd,
                                            kstemmax, stemc,stemd,
                                            kleafmax, leafc,leafd);
  }
  return(exp(krhizomaxlog));
}



/**
 * BIOMECHANICS
 * 
 * Savage, V. M., L. P. Bentley, B. J. Enquist, J. S. Sperry, D. D. Smith, P. B. Reich, and E. I. von Allmen. 2010. Hydraulic trade-offs and space filling enable better predictions of vascular structure and function in plants. Proceedings of the National Academy of Sciences of the United States of America 107:22722–7.
 * Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0:1–60.
 * 
 * height - Tree height in cm
 */

// [[Rcpp::export("hydraulics.taperFactorSavage")]]
double taperFactorSavage(double height) {
  double b_p0 = 1.32, b_p13 = 1.85; //normalizing constants (p = 1/3)
  double a_p0 = 7.20E-13, a_p13 = 6.67E-13;
  double n_ext = 2.0; //Number of daughter branches per parent
  double N = ((3.0*log(1.0-(height/4.0)*(1.0-pow(n_ext, 1.0/3.0))))/log(n_ext))-1.0;
  double K_0 = a_p0*pow(pow(n_ext, N/2.0),b_p0);
  double K_13 = a_p13*pow(pow(n_ext, N/2.0),b_p13);
  return(K_13/K_0);
}

/**
 *  Returns the terminal conduit radius (in micras)
 *  
 *  height - plant height in cm
 */
// [[Rcpp::export("hydraulics.terminalConduitRadius")]]
double terminalConduitRadius(double height) {
  double dh  = pow(10,1.257 +  0.24*log10(height/100.0));//Olson, M.E., Anfodillo, T., Rosell, J.A., Petit, G., Crivellaro, A., Isnard, S., León-Gómez, C., Alvarado-Cárdenas, L.O., & Castorena, M. 2014. Universal hydraulics of the flowering plants: Vessel diameter scales with stem length across angiosperm lineages, habits and climates. Ecology Letters 17: 988–997.
  return(dh/2.0);
}


// [[Rcpp::export("hydraulics.referenceConductivityHeightFactor")]]
double referenceConductivityHeightFactor(double refheight, double height) {
  double rhref  = terminalConduitRadius(refheight);
  double rh  = terminalConduitRadius(height);
  double df = pow(rh/rhref,2.0);
  return(df);
}


/**
 * Calculate maximum leaf-specific stem hydraulic conductance (in mmol·m-2·s-1·MPa-1)
 * 
 * xylemConductivity - Sapwood-specific conductivity of stem xylem (in kg·m-1·s-1·MPa-1), 
 *                     assumed to be measured at distal twigs
 * refheight - Reference plant height (on which xylem conductivity was measured)
 * Al2As - Leaf area to sapwood area ratio (in m2·m-2)
 * height - plant height (in cm)
 * taper - boolean to apply taper
 */
// [[Rcpp::export("hydraulics.maximumStemHydraulicConductance")]]
double maximumStemHydraulicConductance(double xylemConductivity, double refheight, double Al2As, double height, 
                                       bool angiosperm = true, bool taper = false) {
  
  
  // Christoffersen, B. O., M. Gloor, S. Fauset, N. M. Fyllas, D. R. Galbraith, T. R. Baker, L. Rowland, R. A. Fisher, O. J. Binks, S. A. Sevanto, C. Xu, S. Jansen, B. Choat, M. Mencuccini, N. G. McDowell, and P. Meir. 2016. Linking hydraulic traits to tropical forest function in a size-structured and trait-driven model (TFS v.1-Hydro). Geoscientific Model Development Discussions 0:1–60.
  double kmax = 0.0;
  if(!taper) {
    double xylemConductivityCorrected = xylemConductivity*referenceConductivityHeightFactor(refheight, height);
    kmax =   (1000.0/0.018)*(xylemConductivityCorrected/Al2As)*(100.0/height);
  } else {
    double petioleConductivity = xylemConductivity*referenceConductivityHeightFactor(refheight, 100.0);
    // Correct reference conductivity in relation to the reference plant height in which it was measured
    kmax =   (1000.0/0.018)*(petioleConductivity/Al2As)*(100.0/height)*(taperFactorSavage(height)/(taperFactorSavage(100.0)));
  } 
  return(kmax); 
}

/**
 * Calculate maximum leaf-specific root hydraulic conductance (in mmol·m-2·s-1·MPa-1)
 * 
 * xylemConductivity - Sapwood-specific conductivity of root xylem (in kg·m-1·s-1·MPa-1)
 * Al2As - Leaf area to sapwood area ratio (in m2·m-2)
 * v - proportion of fine roots in each soil layer
 * widths - soil layer depths (in mm)
 */
// [[Rcpp::export("hydraulics.maximumRootHydraulicConductance")]]
double maximumRootHydraulicConductance(double xylemConductivity, double Al2As, NumericVector v, 
                                       NumericVector widths, double depthWidthRatio = 1.0){
  NumericVector rl = rootLengths(v,widths, depthWidthRatio);
  NumericVector w = xylemConductanceProportions(v,widths, depthWidthRatio);
  int nlayers = v.length();
  double kmax = 0.0;
  for(int i=0;i<nlayers;i++) kmax = kmax + w[i]*(1000.0/0.018)*(xylemConductivity/((rl[i]/1000.0)*Al2As));
  return(kmax); 
}