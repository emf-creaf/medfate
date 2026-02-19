#include <math.h>
#include "incbeta_c.h"
#include "medfate.h"


// [[Rcpp::export(".gammln")]]
double gammln_c(double xx) {
  static double cof[6] = {76.18009172947146,-86.50532032941677,
                          24.01409824083091, -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  double x = xx - 1.0;
  double tmp = x + 5.5;
  tmp = tmp - (x + 0.5) * std::log(tmp);
  double ser = 1.000000000190015;
  for ( j = 0; j <= 5; j++ ) {
    x = x + 1.0;
    ser = ser + cof[j]/x;
  }
  return (-tmp + std::log(2.5066282746310005*ser));
}

// Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz’s method
// [[Rcpp::export(".betacf")]]
double betacf_c(double a, double b, double x) {
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  qab = a + b; //These q’s will be used in factors that occur in the coefficients
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0; //First step of Lentz’s method.
  d = 1.0 - qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for(m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; //One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h = h*d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; //Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d = 1.0/d;
    del = d*c;
    h = h*del;
    if (fabs(del-1.0) < EPS) break; //Are we done?
  }
  if (m > MAXIT) throw medfate::MedfateInternalError("a or b too big, or MAXIT too small in betacf");
  return h;
}

// Returns the regularized incomplete beta function Ix(a, b).
// [[Rcpp::export(".incbeta")]]
double incbeta_c(double a, double b, double x) {
  double bt;
  if (x < 0.0 || x > 1.0) throw medfate::MedfateInternalError("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) {
    bt=0.0; 
  } else {//Factors in front of the continued fraction.
    bt=exp(gammln_c(a+b)-gammln_c(a)-gammln_c(b)+a*log(x)+b*log(1.0-x));
  }
  if (x < (a + 1.0)/(a + b + 2.0)){ //Use continued fraction directly.
    return bt*betacf_c(a,b,x)/a;
  } else { //Use continued fraction after making the symmetry transformation.
    return 1.0-bt*betacf_c(b,a,1.0-x)/b;
  }
}