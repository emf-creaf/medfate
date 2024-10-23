#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

const double dwarf=0.0000001;
const double giant=999999999.9;
const double explow = -300.0;
const double machtol = 1.0E-15;
const double sqrttwopi=2.5066282746310005024;
const double lnsqrttwopi=0.9189385332046727418;
const double twopi=6.2831853071795864769;
const double oneoversqrtpi=0.5641895835477562869;
const double epss=1.0E-15;

// FUNCTION alfa(x)
// USE Someconstants
// IMPLICIT NONE
// REAL(r8) :: x, alfa, lnx
// lnx=log(x)
// IF (x>0.25) THEN
// alfa=x+0.25_r8
// ELSEIF (x>=dwarf) THEN
// alfa=-0.6931_r8/lnx
// ELSE
// alfa=-0.6931_r8/log(dwarf)
// ENDIF
// END FUNCTION alfa
double alfa(double x){
  double lnx = log(x);
  double alfa;
  if(x>0.25) {
    alfa = x + 0.25;
  } else if (x>=dwarf) {
    alfa = -0.6931/lnx;
  } else {
    alfa = -0.6931/log(dwarf);
  }
  return(alfa);
}

// FUNCTION exmin1(x,eps)
//   USE Someconstants  
//   IMPLICIT NONE
//   ! computes (exp(x)-1)/x 
//   REAL(r8) :: exmin1, x, eps
//   REAL(r8) :: e, t, y
//   IF (x==0) THEN
//   y=1.0_r8
// ELSEIF ((x<-0.69).OR.(x > 0.4)) THEN
//   y=(exp(x)-1.0_r8)/x
//   ELSE
//   t=x/2.0_r8;
// y=exp(t)*sinh(t,eps)/t
//   ENDIF
//   exmin1= y
//   END FUNCTION exmin1
double exmin1(double x) {
  double y;
  if(x==0.0) {
    y=1.0;
  } else if ((x<-0.69) || (x > 0.4)) {
    y=(exp(x)-1.0)/x;
  } {
    double t = x/2.0;
    y = exp(t)*sinh(t)/t;
  }
  return(y);
}

// FUNCTION exmin1minx(x,eps)
//   USE Someconstants    
//   IMPLICIT NONE
//   !{computes (exp(x)-1-x)/(0.5*x*x) }
// REAL(r8) :: exmin1minx, x, eps
//   REAL(r8) :: t, t2, y
//   IF (x==0) THEN
//   y=1.0_r8
// ELSEIF (abs(x)>0.9) THEN
//   y=(exp(x)-1-x)/(x*x/2.0_r8)
//   ELSE
//   t=sinh(x/2.0_r8,eps);
// t2=t*t;
// y=(2*t2+(2*t*sqrt(1.0_r8+t2)-x))/(x*x/2.0_r8)
//   ENDIF
//   exmin1minx=y
//   END FUNCTION exmin1minx
double exmin1minx(double x) {
  double y;
  if(x==0.0) {
    y = 1.0;
  } else if(std::abs(x)>0.9) {
    y = (exp(x)-1.0-x)/(x*x/2.0);
  } else {
    double t=sinh(x/2.0);
    double t2=t*t;
    y=(2.0*t2+(2.0*t*sqrt(1.0+t2)-x))/(x*x/2.0);
  }
  return(y);
}  
// FUNCTION logoneplusx(x)
//   USE Someconstants  
//   IMPLICIT NONE
//   !{x >-1; computes ln(1+x) with good}
// !{relative precision when |x| is small}
// REAL(r8) :: x, logoneplusx
//   REAL(r8) :: y0, r, s
//   y0=log(1.0_r8+x);
// IF ((-0.2928 < x).AND.(x < 0.4142)) THEN
//   s=y0*exmin1(y0, machtol);
// r=(s-x)/(s+1.0_r8);
// y0=y0-r*(6-r)/(6-4*r)
//   ENDIF
//   logoneplusx= y0
//   END FUNCTION logoneplusx
double logoneplusx(double x) {
  double y0=log(1.0+x);
  if((-0.2928 < x) && (x < 0.4142)){
    double s = y0*exmin1(y0);
    double r = (s-x)/(s+1.0);
    y0=y0-r*(6.0-r)/(6.0-4.0*r);
  }
  return(y0);
}

// FUNCTION lnec(x)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: lnec, ln1, x,  y0, z, e2, r, s
//   !{x>-1; lnec:=ln1:=ln(1+x)-x}
// z=logoneplusx(x);
// y0=z-x;
// e2=exmin1minx(z,machtol);
// s=e2*z*z/2;
// r=(s+y0)/(s+1+z);
// ln1=y0-r*(6-r)/(6-4*r);
// lnec=ln1
//   ! lnec := x + ln1
//   END FUNCTION lnec
double lnec(double x) {
  double z = logoneplusx(x);
  double y0 = z - x;
  double e2 = exmin1minx(z);
  double s = e2*z*z/2.0;
  double r = (s+y0)/(s+1.0+z);
  double ln1 = y0 - r*(6.0-r)/(6.0-4.0*r);
  return(ln1);
}

//   FUNCTION chepolsum(n,x,a)
//   USE Someconstants
//   IMPLICIT NONE
//   INTEGER :: n
//   REAL(r8) :: chepolsum, x, a(0:n)
//   REAL(r8) :: h, r, s, tx
//   INTEGER :: k
//   !{a[0]/2+a[1]T1(x)+...a[n]Tn(x); series of Chebychev polynomials}
// IF (n==0) THEN
//   chepolsum=a(0)/2.0_r8
// ELSEIF (n==1) THEN
//   chepolsum=a(0)/2.0_r8+a(1)*x
//   ELSE
//   tx=x+x;
// r=a(n);
// h=a(n-1)+r*tx;
// DO k=n-2,1,-1 
// s=r;
// r=h;
// h=a(k)+r*tx-s
//   ENDDO
//   chepolsum=a(0)/2.0_r8-r+h*x
//   ENDIF
//   END FUNCTION chepolsum
double chepolsum(double x, NumericVector a) {
  int n = a.size() - 1;
  if(n==0) {
    return(a[0]/2.0);
  } else if (n==1) {
    return(a[0]/2.0+a[1]*x);
  } else {
    double tx=x+x;
    double r=a[n];
    double h=a[n-1]+r*tx;
    for(int k=n-2;k>0;k--) {
      double s = r;
      r = h;
      h = a[k]+r*tx-s;
    }
    return(a[0]/2.0-r+h*x);
  }
}
// RECURSIVE FUNCTION auxgam(x) RESULT(auxgamm)
//   USE Someconstants
//   IMPLICIT NONE    
//   !{function g in 1/gamma(x+1)=1+x*(x-1)*g(x), -1<=x<=1}
// REAL(r8) :: auxgamm, x
//   REAL(r8) :: t, dr(0:17)
//   IF (x<0) THEN
//   auxgamm=-(1.0_r8+(1+x)*(1+x)*auxgam(1+x))/(1.0_r8-x)
//   ELSE
//   dr(0)= -1.013609258009865776949_r8;
// dr(1)= 0.784903531024782283535e-1_r8;
// dr(2)= 0.67588668743258315530e-2_r8;
// dr(3)= -0.12790434869623468120e-2_r8;
// dr(4)= 0.462939838642739585e-4_r8;
// dr(5)= 0.43381681744740352e-5_r8;
// dr(6)= -0.5326872422618006e-6_r8;
// dr(7)= 0.172233457410539e-7_r8;
// dr(8)= 0.8300542107118e-9_r8;
// dr(9)= -0.10553994239968e-9_r8;
// dr(10)= 0.39415842851e-11_r8;
// dr(11)= 0.362068537e-13_r8;
// dr(12)= -0.107440229e-13_r8;
// dr(13)= 0.5000413e-15_r8;
// dr(14)= -0.62452e-17_r8;
// dr(15)= -0.5185e-18_r8;
// dr(16)= 0.347e-19_r8;
// dr(17)= -0.9e-21_r8;
// t=2*x-1.0_r8;
// auxgamm=chepolsum(17,t,dr);
// ENDIF
//   END FUNCTION auxgam
double auxgam(double x) {
  NumericVector dr(18);
  double auxgamm;
  if(x<0.0) {
    auxgamm = -(1.0+(1.0+x)*(1.0+x)*auxgam(1.0+x))/(1.0-x);
  } else {
    dr[0]= -1.013609258009865776949;
    dr[1]= 0.784903531024782283535e-1;
    dr[2]= 0.67588668743258315530e-2;
    dr[3]= -0.12790434869623468120e-2;
    dr[4]= 0.462939838642739585e-4;
    dr[5]= 0.43381681744740352e-5;
    dr[6]= -0.5326872422618006e-6;
    dr[7]= 0.172233457410539e-7;
    dr[8]= 0.8300542107118e-9;
    dr[9]= -0.10553994239968e-9;
    dr[10]= 0.39415842851e-11;
    dr[11]= 0.362068537e-13;
    dr[12]= -0.107440229e-13;
    dr[13]= 0.5000413e-15;
    dr[14]= -0.62452e-17;
    dr[15]= -0.5185e-18;
    dr[16]= 0.347e-19;
    dr[17]= -0.9e-21;
    double t=2*x-1.0;
    auxgamm=chepolsum(t,dr);
  }
  return(auxgamm);
}

// FUNCTION lngam1(x)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: x, lngam1   
//   ! {ln(gamma(1+x)), -1<=x<=1}
// lngam1=-logoneplusx(x*(x-1.0_r8)*auxgam(x))
//   END FUNCTION lngam1
double lngam1(double x) {
  return(-logoneplusx(x*(x-1.0)*auxgam(x)));
}
  

  
// FUNCTION stirling(x)            
//   !{Stirling series, function corresponding with}
// !{asymptotic series for log(gamma(x))}
// !{that is:  1/(12x)-1/(360x**3)...; x>= 3}
// USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: stirling, x, a(0:17), c(0:6), z
//   IF (x<dwarf) THEN
//   stirling=giant
//   ELSEIF (x<1) THEN
//   stirling= lngam1(x)-(x+0.5)*log(x)+x-lnsqrttwopi
//   ELSEIF (x<2) THEN
//   stirling=lngam1(x-1)-(x-0.5)*log(x)+x-lnsqrttwopi
//   ELSEIF (x<3) THEN
//   stirling=lngam1(x-2)-(x-0.5)*log(x)+x-lnsqrttwopi+log(x-1)
//   ELSEIF (x<12) THEN
//   a(0)=1.996379051590076518221_r8;
// a(1)=-0.17971032528832887213e-2_r8;
// a(2)=0.131292857963846713e-4_r8;
// a(3)=-0.2340875228178749e-6_r8;
// a(4)=0.72291210671127e-8_r8;
// a(5)=-0.3280997607821e-9_r8;
// a(6)=0.198750709010e-10_r8;
// a(7)=-0.15092141830e-11_r8;
// a(8)=0.1375340084e-12_r8;
// a(9)=-0.145728923e-13_r8;
// a(10)=0.17532367e-14_r8;
// a(11)=-0.2351465e-15_r8;
// a(12)=0.346551e-16_r8;
// a(13)=-0.55471e-17_r8;
// a(14)=0.9548e-18_r8;
// a(15)=-0.1748e-18_r8;
// a(16)=0.332e-19_r8;
// a(17)=-0.58e-20_r8;
// z=18.0_r8/(x*x)-1.0_r8;
// stirling=chepolsum(17,z,a)/(12.0_r8*x);
// ELSE
//   z=1.0_r8/(x*x);
// IF (x<1000) THEN
//   c(0)=0.25721014990011306473e-1_r8;
// c(1)=0.82475966166999631057e-1_r8;
// c(2)=-0.25328157302663562668e-2_r8;
// c(3)=0.60992926669463371e-3_r8;
// c(4)=-0.33543297638406e-3_r8;
// c(5)=0.250505279903e-3_r8;
// c(6)=0.30865217988013567769_r8;
// stirling=((((((c(5)*z+c(4))*z+c(3))*z+c(2))*z+c(1))*z+c(0))/(c(6)+z)/x)
//   ELSE
//   stirling=(((-z/1680.0_r8+1.0_r8/1260.0_r8)*z-1.0_r8/360.0_r8)*z+1.0_r8/12.0_r8)/x
//   ENDIF
//   ENDIF 
//   END FUNCTION stirling
double  stirling(double x) {
  double stirling, z;
  NumericVector a(18);
  NumericVector c(7);
  if(x<dwarf) {
    stirling = giant; 
  } else if(x<1.0) {
    stirling= lngam1(x)-(x+0.5)*log(x)+x-lnsqrttwopi;
  }
  else if(x<2.0) {
    stirling=lngam1(x-1.0)-(x-0.5)*log(x)+x-lnsqrttwopi;
  }
  else if(x<3.0) {
    stirling=lngam1(x-2.0)-(x-0.5)*log(x)+x-lnsqrttwopi+log(x-1.0);
  }
  else if(x<12.0) {
    a[0]=1.996379051590076518221;
    a[1]=-0.17971032528832887213e-2;
    a[2]=0.131292857963846713e-4;
    a[3]=-0.2340875228178749e-6;
    a[4]=0.72291210671127e-8;
    a[5]=-0.3280997607821e-9;
    a[6]=0.198750709010e-10;
    a[7]=-0.15092141830e-11;
    a[8]=0.1375340084e-12;
    a[9]=-0.145728923e-13;
    a[10]=0.17532367e-14;
    a[11]=-0.2351465e-15;
    a[12]=0.346551e-16;
    a[13]=-0.55471e-17;
    a[14]=0.9548e-18;
    a[15]=-0.1748e-18;
    a[16]=0.332e-19;
    a[17]=-0.58e-20;
    z=18.0/(x*x)-1.0;
    stirling=chepolsum(z,a)/(12.0*x);
  } else {
    z=1.0/(x*x);
    if(x<1000.0) {
      c[0]=0.25721014990011306473e-1;
      c[1]=0.82475966166999631057e-1;
      c[2]=-0.25328157302663562668e-2;
      c[3]=0.60992926669463371e-3;
      c[4]=-0.33543297638406e-3;
      c[5]=0.250505279903e-3;
      c[6]=0.30865217988013567769;
      stirling=((((((c[5]*z+c[4])*z+c[3])*z+c[2])*z+c[1])*z+c[0])/(c[6]+z)/x);
    } else {
      stirling=(((-z/1680.0+1.0/1260.0)*z-1.0/360.0)*z+1.0/12.0)/x;
    }
  }
  return(stirling);
}

// FUNCTION gamstar(x) 
//   USE Someconstants
//   !  {gamstar(x)=exp(stirling(x)), x>0; or }
// !  {gamma(x)/(exp(-x+(x-0.5)*ln(x))/sqrt(2pi)}
// IMPLICIT NONE
//   REAL(r8) :: gamstar, x
//   IF (x>=3) THEN
//   gamstar=exp(stirling(x))
//   ELSEIF (x>0) THEN
//   gamstar=gamma(x)/(exp(-x+(x-0.5)*log(x))*sqrttwopi)
//   ELSE
//   gamstar=giant
//   ENDIF
//   END FUNCTION gamstar
double gamstar(double x) {
  double gamstar;
  if(x>=3.0) {
    gamstar = exp(stirling(x));
  } else if(x>0.0) {
    gamstar= tgamma(x)/(exp(-x+(x-0.5)*log(x))*sqrttwopi);
  } else {
    gamstar = giant;
  }
  return(gamstar);
}
// FUNCTION  dompart(a,x,qt)
//   ! dompart is approx. of  x^a * exp(-x) / gamma(a+1)   
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: dompart, a, x
//   REAL(r8) :: lnx, c, dp, la, mu, r 
//   LOGICAL :: qt
//   lnx=log(x)
//   IF (a<=1) THEN                     
//      r=-x+a*lnx
//   ELSE
//    IF (x==a) THEN
//      r=0
//    ELSE
//      la=x/a
//      r=a*(1.0_r8-la+log(la))
//    ENDIF
//    r=r-0.5_r8*log(6.2832_r8*a)
//   ENDIF
//   IF (r<explow) THEN
//   dp=0.0_r8
// ELSE
//   dp=exp(r)
//   ENDIF
//   IF (qt) THEN
//   dompart=dp
//   ELSE
//   IF ((a<3).OR.(x<0.2_r8)) THEN
//   dompart=exp(a*lnx-x)/gamma(a+1.0_r8)
//   ELSE
//   mu=(x-a)/a;
// c=lnec(mu);
// IF ((a*c)>log(giant)) THEN
//   dompart=-100
// ELSE
//   dompart=exp(a*c)/(sqrt(a*2*pi)*gamstar(a))
//   ENDIF
//   ENDIF
//   ENDIF
//   END FUNCTION dompart
double dompart(double a, double x, bool qt) {
  double lnx = log(x);
  double r, la, dp, dompart;
  if(a<=1) {
    r = -x + (a*lnx);
  } else {
    if(x==a) {
      r = 0.0;
    } else {
      la = x/a;
      r = a*(1.0 - la + log(la));
    }
    r = r -0.5*log(6.2832*a);
  }
  if(r<explow)  {
    dp = 0.0;
  } else {
    dp = exp(r);
  }
  if(qt) {
    dompart = dp;
  } else {
    if((a<3.0) || (x<0.2)) {
      dompart = exp(a*lnx - x)/tgamma(a+1.0); //using in-built tgamma instead of fortran gamma
    } else {
      double mu = (x-a)/a;
      double c = lnec(mu);
      if((a*c)>log(giant)) {
        dompart = -100.0;
      } else {
        dompart = exp(a*c)/(sqrt(a*2.0*M_PI)*gamstar(a));
      }
    }
  }
  return(dompart);
}
 
 
 // FUNCTION fractio(x,n,r,s)
 //  IMPLICIT NONE
 //  INTEGER n,k
 //  REAL(r8) :: x, fractio, r(0:8), s(0:8), a, b
 //  a=r(n); b=1
 //  DO k=n-1,0,-1 
 //  a=a*x+r(k); b=b*x+s(k) 
 //    ENDDO
 //    fractio=a/b
 //    END FUNCTION fractio
double fractio(double x, int n, double r[9], double s[9]) {
  double a=r[n];
  double b=1.0;
  for(int k=n-1;k>=0;k--) {
    a=a*x+r[k]; 
    b=b*x+s[k];
  } 
  return(a/b);
}
  
//  RECURSIVE FUNCTION errorfunction (x, erfcc, expo) RESULT(errfu)
//   ! coefficients are from Cody (1969), Math. Comp., 23, 631-637
// USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: x, y, z, r(0:8), s(0:8), errfu
//   LOGICAL erfcc, expo
//   IF (erfcc) THEN
//     IF (x < -6.5) THEN
//        y= 2.0_r8 
//     ELSEIF (x < 0) THEN
//        y= 2.0_r8 - errorfunction(-x, .true., .false.) 
//     ELSEIF (x == 0) THEN
//        y= 1.0_r8 
//     ELSEIF (x < 0.5) THEN
//       IF (expo) THEN
//        y=exp(x*x)
//       ELSE
//         y=1.0_r8
//       ENDIF
//       y=y*(1.0_r8 - errorfunction(x, .false., .false.))
//     ELSEIF (x < 4) THEN
//       IF (expo) THEN
//         y= 1.0_r8 
//       ELSE
//        y= exp(-x*x)
//       ENDIF
//   r(0)= 1.230339354797997253e3_r8
// r(1)= 2.051078377826071465e3_r8
// r(2)= 1.712047612634070583e3_r8
// r(3)= 8.819522212417690904e2_r8
// r(4)= 2.986351381974001311e2_r8
// r(5)= 6.611919063714162948e1_r8
// r(6)= 8.883149794388375941_r8
// r(7)= 5.641884969886700892e-1_r8
// r(8)= 2.153115354744038463e-8_r8
// s(0)= 1.230339354803749420e3_r8
// s(1)= 3.439367674143721637e3_r8
// s(2)= 4.362619090143247158e3_r8
// s(3)= 3.290799235733459627e3_r8
// s(4)= 1.621389574566690189e3_r8
// s(5)= 5.371811018620098575e2_r8
// s(6)= 1.176939508913124993e2_r8
// s(7)= 1.574492611070983473e1_r8
// y=y*fractio(x,8,r,s)
//   ELSE
//   z=x*x
//   IF (expo) THEN
//   y=1.0_r8 
// ELSE
//   y= exp(-z)
//   ENDIF
//   z=1.0_r8/z
//   r(0)=6.587491615298378032e-4_r8
// r(1)=1.608378514874227663e-2_r8
// r(2)=1.257817261112292462e-1_r8
// r(3)=3.603448999498044394e-1_r8
// r(4)=3.053266349612323440e-1_r8
// r(5)=1.631538713730209785e-2_r8
// s(0)=2.335204976268691854e-3_r8
// s(1)=6.051834131244131912e-2_r8
// s(2)=5.279051029514284122e-1_r8
// s(3)=1.872952849923460472_r8
// s(4)=2.568520192289822421_r8
// y=y*(oneoversqrtpi-z*fractio(z,5,r,s))/x
//   ENDIF
//   errfu=y
//   ELSE
//   IF (x == 0.0_r8) THEN 
//   y=0
// ELSEIF (abs(x) > 6.5) THEN 
//   y=x/abs(x)
//   ELSEIF (x > 0.5) THEN
//   y=1.0_r8 - errorfunction(x, .true., .false.) 
//   ELSEIF (x < -0.5) THEN
//   y=errorfunction(-x, .true., .false.)-1.0_r8
// ELSE
//   r(0)=3.209377589138469473e3_r8
// r(1)=3.774852376853020208e2_r8
// r(2)=1.138641541510501556e2_r8
// r(3)=3.161123743870565597e0_r8
// r(4)=1.857777061846031527e-1_r8
// s(0)=2.844236833439170622e3_r8
// s(1)=1.282616526077372276e3_r8
// s(2)=2.440246379344441733e2_r8
// s(3)=2.360129095234412093e1_r8
// z=x*x
//   y=x*fractio(z,4,r,s)
//   ENDIF  
//   errfu= y
//   ENDIF        
//   END FUNCTION errorfunction
double errorfunction(double x, bool erfcc, bool expo){
  double y, z, errfu;
  double r[9], s[9];
  if(erfcc) {
    if(x<-6.5) {
      y = 2.0;
    } else if(x<0.0) {
      y= 2.0 - errorfunction(-x, true, false);
    } else if(x==0.0) {
      y = 1.0;
    } else if(x<0.5) {
      if(expo) {
        y=exp(x*x);
      } else {
        y=1.0;  
      }
      y=y*(1.0 - errorfunction(x, false, false));
    } else if(x<4.0) {
      if(expo) {
        y= 1.0;  
      } else {
        y= exp(-x*x);  
      }
      r[0]= 1.230339354797997253e3;
      r[1]= 2.051078377826071465e3;
      r[2]= 1.712047612634070583e3;
      r[3]= 8.819522212417690904e2;
      r[4]= 2.986351381974001311e2;
      r[5]= 6.611919063714162948e1;
      r[6]= 8.883149794388375941;
      r[7]= 5.641884969886700892e-1;
      r[8]= 2.153115354744038463e-8;
      s[0]= 1.230339354803749420e3;
      s[1]= 3.439367674143721637e3;
      s[2]= 4.362619090143247158e3;
      s[3]= 3.290799235733459627e3;
      s[4]= 1.621389574566690189e3;
      s[5]= 5.371811018620098575e2;
      s[6]= 1.176939508913124993e2;
      s[7]= 1.574492611070983473e1;
      y=y*fractio(x,8,r,s);
    } else {
      z=x*x;
      if(expo) {
        y = 1.0;
      } else {
        y= exp(-z);
      }
      z=1.0/z;
      r[0]=6.587491615298378032e-4;
      r[1]=1.608378514874227663e-2;
      r[2]=1.257817261112292462e-1;
      r[3]=3.603448999498044394e-1;
      r[4]=3.053266349612323440e-1;
      r[5]=1.631538713730209785e-2;
      s[0]=2.335204976268691854e-3;
      s[1]=6.051834131244131912e-2;
      s[2]=5.279051029514284122e-1;
      s[3]=1.872952849923460472;
      s[4]=2.568520192289822421;
      y=y*(oneoversqrtpi-z*fractio(z,5,r,s))/x;
    }
    errfu = y;
  } else {
    if(x==0.0) {
      y = 0.0;
    } else if(std::abs(x)>6.5) {
      y = x/std::abs(x);
    } else if(x > 0.5) {
      y=1.0 - errorfunction(x, true, false);
    } else if(x < -0.5) {
      y=errorfunction(-x, true, false)-1.0;
    } else {
      r[0]=3.209377589138469473e3;
      r[1]=3.774852376853020208e2;
      r[2]=1.138641541510501556e2;
      r[3]=3.161123743870565597e0;
      r[4]=1.857777061846031527e-1;
      s[0]=2.844236833439170622e3;
      s[1]=1.282616526077372276e3;
      s[2]=2.440246379344441733e2;
      s[3]=2.360129095234412093e1;
      z=x*x;
      y=x*fractio(z,4,r,s);
    }
    errfu= y;
  }
  return(errfu);
}

// FUNCTION saeta(a,eta)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: a, eta, saeta, y, s, t
//   REAL(r8) :: eps, fm(0:26), bm(0:26)
//   INTEGER  ::  m
//   eps=epss
//   fm(0)=1.0_r8;
// fm(1)=-1.0_r8/3.0_r8;
// fm(2)=1.0_r8/12.0_r8;
// fm(3)=-2.0_r8/135.0_r8;
// fm(4)=1.0_r8/864.0_r8;
// fm(5)=1.0_r8/ 2835.0_r8;
// fm(6)=-139.0_r8/777600.0_r8;
// fm(7)=1.0_r8/25515.0_r8;
// fm(8)=-571.0_r8/261273600.0_r8;
// fm(9)=-281.0_r8/151559100.0_r8;
// fm(10)=8.29671134095308601e-7_r8;
// fm(11)=-1.76659527368260793e-7_r8;
// fm(12)=6.70785354340149857e-9_r8;
// fm(13)=1.02618097842403080e-8_r8;
// fm(14)=-4.38203601845335319e-9_r8;
// fm(15)=9.14769958223679023e-10_r8;
// fm(16)=-2.55141939949462497e-11_r8;
// fm(17)=-5.83077213255042507e-11_r8;
// fm(18)=2.43619480206674162e-11_r8;
// fm(19)=-5.02766928011417559e-12_r8;
// fm(20)=1.10043920319561347e-13_r8;
// fm(21)=3.37176326240098538e-13_r8;
// fm(22)=-1.39238872241816207e-13_r8;
// fm(23)=2.85348938070474432e-14_r8;
// fm(24)=-5.13911183424257258e-16_r8;
// fm(25)=-1.97522882943494428e-15_r8;
// fm(26)= 8.09952115670456133e-16_r8;
// bm(25)=fm(26);
// bm(24)=fm(25);
// DO m=24,1,-1 
// bm(m-1)=fm(m)+(m+1)*bm(m+1)/a;
// ENDDO
//   s=bm(0);
// t=s;
// y=eta;
// m=1;
// DO WHILE ((abs(t/s)>eps).AND.(m<25))
//   t=bm(m)*y;
// s=s+t;
// m=m+1;
// y=y*eta
//   ENDDO 
//   saeta=s/(1.0_r8+bm(1)/a);
// END FUNCTION saeta
double saeta(double a, double eta){
  double saeta, y, s, t, eps;
  int m;
  double fm[27], bm[27];
  eps=epss;
  fm[0]=1.0;
  fm[1]=-1.0/3.0;
  fm[2]=1.0/12.0;
  fm[3]=-2.0/135.0;
  fm[4]=1.0/864.0;
  fm[5]=1.0/ 2835.0;
  fm[6]=-139.0/777600.0;
  fm[7]=1.0/25515.0;
  fm[8]=-571.0/261273600.0;
  fm[9]=-281.0/151559100.0;
  fm[10]=8.29671134095308601e-7;
  fm[11]=-1.76659527368260793e-7;
  fm[12]=6.70785354340149857e-9;
  fm[13]=1.02618097842403080e-8;
  fm[14]=-4.38203601845335319e-9;
  fm[15]=9.14769958223679023e-10;
  fm[16]=-2.55141939949462497e-11;
  fm[17]=-5.83077213255042507e-11;
  fm[18]=2.43619480206674162e-11;
  fm[19]=-5.02766928011417559e-12;
  fm[20]=1.10043920319561347e-13;
  fm[21]=3.37176326240098538e-13;
  fm[22]=-1.39238872241816207e-13;
  fm[23]=2.85348938070474432e-14;
  fm[24]=-5.13911183424257258e-16;
  fm[25]=-1.97522882943494428e-15;
  fm[26]= 8.09952115670456133e-16;
  bm[25]=fm[26];
  bm[24]=fm[25];
  for(m=24;m>=1; m--) {
    bm[m-1]=fm[m]+(((double)m)+1.0)*bm[m+1]/a;
  }
  s=bm[0];
  t=s;
  y=eta;
  m=1;
  while((std::abs(t/s)>eps) && (m<25)){
    t=bm[m]*y;
    s=s+t;
    m=m+1;
    y=y*eta;
  }
  saeta=s/(1.0+bm[1]/a);
  return(saeta);
}

//  FUNCTION qtaylor(a,x,dp)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: qtaylor, a, x, dp
//   REAL(r8) :: eps, lnx, p, q, r, s, t, u, v
//   eps=epss
//   lnx=log(x)
//   IF (dp==0) THEN
//   q=0.0_r8
// ELSE
//   r=a*lnx;
// q=r*exmin1(r,eps);   ! {q = x^a - 1 }
// s=a*(1.0_r8-a)*auxgam(a); ! {s = 1-1/Gamma(1+a) }
// q=(1-s)*q;
// u=s-q;               ! {u = 1 - x^a/Gamma(1+a)}
// p=a*x;
// q=a+1;
// r=a+3;
// t=1.0_r8;
// v=1.0_r8;
// DO WHILE (abs(t / v) > eps)
//   p=p+x;
// q=q+r;
// r=r+2;
// t=-p*t/q;
// v=v+t
//   ENDDO
//   v=a*(1-s)*exp((a+1.0_r8)*lnx)*v/(a+1.0_r8);
// q=u+v
//   ENDIF
//   qtaylor=q
//   END FUNCTION qtaylor
double qtaylor(double a, double x, double dp){
  double eps, lnx, p, q, r, s, t, u, v;
  eps=epss;
  lnx=log(x);
  if(dp==0.0)  {
    q=0.0;
  } else {
    r = a*lnx;
    q=r*exmin1(r);
    s=a*(1.0-a)*auxgam(a); // {s = 1-1/Gamma(1+a) }
    q=(1.0-s)*q;
    u=s-q;  // {u = 1 - x^a/Gamma(1+a)}
    p=a*x;
    q=a+1.0;
    r=a+3.0;
    t=1.0;
    v=1.0;
    while(std::abs(t/v)>eps) {
      p=p+x;
      q=q+r;
      r=r+2.0;
      t=-p*t/q;
      v=v+t;
    }
    v=a*(1.0-s)*exp((a+1.0)*lnx)*v/(a+1.0);
    q=u+v;
  }
  return(q);
}

// FUNCTION ptaylor(a,x,dp)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: ptaylor,a,x,dp
//   REAL(r8) :: eps,p,c,r
//   eps=epss
//   IF (dp==0) THEN
//   p=0.0_r8
// ELSE
//   p=1.0_r8
// c=1.0_r8
// r=a
//   DO WHILE ((c/p)>eps)
//   r=r+1
// c=x*c/r
//   p=p+c
//   ENDDO
//   p=p*dp
//   ENDIF
//   ptaylor=p
//   END FUNCTION ptaylor
double ptaylor(double a, double x, double dp) {
  double eps,p,c,r;
  eps=epss;
  if(dp==0) {
    p = 0.0;
  } else {
    p = 1.0;
    c = 1.0;
    r = a;
    while((c/p)>eps) {
      r=r+1.0;
      c=x*c/r;
      p=p+c;
      
    }
    p=p*dp;
  }
  return(p);
}  
  
//   FUNCTION  pqasymp (a,x,dp,p)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: pqasymp, a, x, dp
//   REAL(r8) :: y, mu, eta, u, v 
//   INTEGER :: s
//   LOGICAL :: p
//   IF (dp==0.0_r8) THEN
//   IF (p) THEN
//   pqasymp=0.0_r8
// ELSE
//   pqasymp=1.0_r8
// ENDIF
//   ELSE
//   IF (p) THEN
//   s=-1
// ELSE
//   s=1
// ENDIF
//   mu=(x-a)/a;
// y=-lnec(mu)
//   IF (y<0) THEN
//   eta=0.0_r8
// ELSE
//   eta=sqrt(2.0_r8*y)
//   ENDIF
//   y=y*a;
// v=sqrt(abs(y));
// IF (mu<0.0_r8) THEN     
//   eta=-eta
//   v=-v
//   ENDIF  
//   u=0.5_r8*errorfunction(s*v,.true.,.false.);
// v=s*exp(-y)*saeta(a,eta)/sqrt(2.0_r8*pi*a);
// pqasymp=u+v
//   ENDIF
//   END FUNCTION pqasymp
double pqasymp(double a, double x, double dp, bool p) {
  double y, mu, eta, u, v;
  double pqasymp;
  int s;
  if(dp==0.0) {
    if(p) pqasymp = 0.0;
    else pqasymp = 1.0;
  } else {
    if(p) {
      s = -1;
    } else {
      s = 1;
    }
    mu=(x-a)/a;
    y=-lnec(mu);
    if(y<0.0) {
      eta = 0.0;
    } else {
      eta=sqrt(2.0*y);
    }
    y=y*a;
    v=sqrt(std::abs(y));
    if(mu<0.0) {
      eta = -eta;
      v = -v;
    }
    u=0.5*errorfunction(s*v,true,false);
    v=s*exp(-y)*saeta(a,eta)/sqrt(2.0*M_PI*a);
    pqasymp=u+v;
  }
  return(pqasymp);
}
  
//   FUNCTION qfraction(a,x,dp)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: qfraction, a, x, dp
//   REAL(r8) :: eps, g, p, q, r, s, t, tau, ro
//   eps=epss
//   IF (dp==0) THEN
//   q=0.0_r8
// ELSE
//   p=0;
// q=(x-1.0_r8-a)*(x+1.0_r8-a);
// r=4*(x+1.0_r8-a);
// s=1.0_r8-a;
// ro=0.0_r8;
// t=1.0_r8;
// g=1.0_r8;
// DO WHILE(abs(t/g)>=eps)
//   p=p+s;
// q=q+r;
// r=r+8;
// s=s+2;
// tau=p*(1.0_r8+ro);
// ro=tau/(q-tau);
// t=ro*t;
// g=g+t
//   ENDDO
//   q=(a/(x+1.0_r8-a))*g*dp; 
// ENDIF
//   qfraction= q
//   END FUNCTION qfraction
double qfraction(double a, double x, double dp){
  double eps, g, p, q, r, s, t, tau, ro;
  eps=epss;
  if(dp==0.0) {
    q = 0.0;
  } else {
    p=0;
    q=(x-1.0-a)*(x+1.0-a);
    r=4.0*(x+1.0-a);
    s=1.0-a;
    ro=0.0;
    t=1.0;
    g=1.0;
    while(std::abs(t/g)>eps) {
      p=p+s;
      q=q+r;
      r=r+8.0;
      s=s+2.0;
      tau=p*(1.0+ro);
      ro=tau/(q-tau);
      t=ro*t;
      g=g+t;
    } 
    q=(a/(x+1.0-a))*g*dp; 
  }
  return(q);
}  
  
  
//   SUBROUTINE incgam(a,x,p,q,ierr)
//   ! -------------------------------------------------------------
//     ! Calculation of the incomplete gamma functions ratios P(a,x)
//     ! and Q(a,x).
//     ! -------------------------------------------------------------
//     ! Inputs:
//     !   a ,    argument of the functions
//     !   x ,    argument of the functions
//     ! Outputs:
//     !   p,     function P(a,x)
//     !   q,     function Q(a,x)  
//     !   ierr , error flag
//     !          ierr=0, computation succesful
//     !          ierr=1, overflow/underflow problems. The function values 
//     !          (P(a,x) and Q(a,x)) are set to zero.
//     ! ----------------------------------------------------------------------
//     ! Authors:
//     !  Amparo Gil    (U. Cantabria, Santander, Spain)
//     !                 e-mail: amparo.gil@unican.es
//     !  Javier Segura (U. Cantabria, Santander, Spain)
//     !                 e-mail: javier.segura@unican.es
//     !  Nico M. Temme (CWI, Amsterdam, The Netherlands)
//     !                 e-mail: nico.temme@cwi.nl
//     ! -------------------------------------------------------------
//     !  References: "Efficient and accurate algorithms for 
//     !  the computation and inversion of the incomplete gamma function ratios",    
//     !  A. Gil, J. Segura and N.M. Temme, submitted to SIAM J Sci Comput
//     ! -------------------------------------------------------------------
//     USE Someconstants    
//       IMPLICIT NONE
//       REAL(r8), INTENT(IN) :: a
//       REAL(r8), INTENT(IN) :: x
//       REAL(r8), INTENT(OUT) :: p
//       REAL(r8), INTENT(OUT) :: q   
//       REAL(r8) :: lnx, dp   
//       INTEGER,  INTENT(OUT) :: ierr
//       ierr=0
//     IF (x<dwarf) THEN
//       lnx=log(dwarf)
//     ELSE
//       lnx=log(x)
//     ENDIF
//     IF (a>alfa(x)) THEN  
//       dp=dompart(a,x,.false.);
//       IF (dp<0) THEN
//        ierr=1
//         p=0; q=0;
//      ELSE
//       IF ((x < 0.3*a).OR.(a<12)) THEN
//         p=ptaylor(a,x,dp)
//       ELSE
//         p=pqasymp(a,x,dp,.true.)
//       ENDIF
//       q=1.0_r8-p
//      ENDIF
//     ELSE
//      IF (a<-dwarf/lnx) THEN
//        q=0.0_r8
//      ELSE
//        IF (x<1.0_r8) THEN
//          dp=dompart(a,x,.true.)
//          IF (dp<0) THEN
//            ierr=1
//            q=0; p=0;
//          ELSE
//           q=qtaylor(a,x,dp)
//           p=1.0_r8-q
//          ENDIF
//       ELSE
//         dp=dompart(a,x,.false.);
//         IF (dp<0) THEN
//          ierr=1
//          p=0; q=0;
//         ELSE
//          IF ((x>2.35_r8*a).OR.(a<12)) THEN
//            q=qfraction(a,x,dp)
//          ELSE
//            print*,'a, dp',a,dp
//            q=pqasymp(a,x,dp,.false.)
//          ENDIF
//          p=1.0_r8-q
//         ENDIF
//       ENDIF
//     ENDIF
//   ENDIF
//   END SUBROUTINE incgam
  
// [[Rcpp::export(".incgam")]]
NumericVector incgam(double a, double x) {
  double lnx, p = NA_REAL, q = NA_REAL;
  double dp;
  if(x<dwarf) {
    lnx = log(dwarf);
  } else {
    lnx = log(x); 
  }
  if(a>alfa(x)) {
    dp = dompart(a,x, false);
    if(dp<0.0) {
      stop("dp < 0");
    } else {
      if ((x < 0.3*a) || (a<12.0)) {
        p=ptaylor(a,x,dp);
      } else {
        p=pqasymp(a,x,dp, true);
      }
      q = 1.0 - p;
    }
  } else {
    if(a< -dwarf/lnx) {
      q = 0.0;
    } else {
      if(x<1.0) {
        dp=dompart(a,x,true);
        if(dp<0.0) {
          stop("dp < 0");
        } else {
          q=qtaylor(a,x,dp);
          p=1.0-q;
        }
      } else {
        dp=dompart(a,x,false);
        if(dp<0.0) {
          stop("dp < 0");
        } else {
          if((x>2.35*a) || (a<12.0)) {
            q = qfraction(a,x,dp);
          } else {
            q=pqasymp(a,x,dp,false);
          }
          p=1.0-q;
        }
      }
    }
  }
  return(NumericVector::create(p,q));
}


// FUNCTION invq(x)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: invq, x, t
//   !  Abramowitx & Stegun 26.2.23; 
// t=sqrt(-2*log(x));
// t=t-(2.515517_r8+t*(0.802853_r8+t*0.010328_r8))/&
//   (1.0_r8+t*(1.432788+t*(0.189269_r8+t*0.001308_r8)))
//   invq=t
//   END FUNCTION invq
double invq(double x)  {
  double t;
  t=sqrt(-2.0*log(x));
  t=t-(2.515517+t*(0.802853+t*0.010328))/(1.0+t*(1.432788+t*(0.189269+t*0.001308)));
  return(t);
}

// RECURSIVE FUNCTION inverfc(x) RESULT(y)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) ::  x, y, y0, y02, h, r, f, fp, c1, c2, c3, c4, c5;
// IF (x > 1) THEN
//   y=-inverfc(2-x)
//   ELSE
//   y0=0.70710678*invq(x/2.0_r8);
// f= erfc(y0)-x;
// f=errorfunction(y0,.true.,.false.)-x;
// y02= y0*y0;
// fp=-2.0_r8/sqrt(pi)*exp(-y02);
// c1=-1.0_r8/fp;
// c2= y0;
// c3=(4*y02+1)/3.0_r8;
// c4=y0*(12*y02+7)/6.0_r8;
// c5=(8*y02+7)*(12*y02+1)/30.0_r8;
// r= f*c1;
// h=r*(1+r*(c2+r*(c3+r*(c4+r*c5))));
// y=y0+h
//   ENDIF
//   END FUNCTION inverfc
double inverfc(double x)  {
  double y, y0, y02, h, r, f, fp, c1, c2, c3, c4, c5;
  if(x > 1.0)  {
    y=-inverfc(2.0-x);
  } else {
    y0=0.70710678*invq(x/2.0);
    // f= erfc(y0)-x;
    f=errorfunction(y0,true,false)-x;
    y02= y0*y0;
    fp=-2.0/sqrt(M_PI)*exp(-y02);
    c1=-1.0/fp;
    c2= y0;
    c3=(4.0*y02+1.0)/3.0;
    c4=y0*(12.0*y02+7.0)/6.0;
    c5=(8.0*y02+7.0)*(12.0*y02+1.0)/30.0;
    r= f*c1;
    h=r*(1.0+r*(c2+r*(c3+r*(c4+r*c5))));
    y=y0+h;
  }
  return(y);
}

// FUNCTION ratfun(x,ak,bk)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: ratfun, x, ak(0:4), bk(0:4), p, q
//   p= ak(0)+x*(ak(1)+x*(ak(2)+x*(ak(3)+x*ak(4))));
// q= bk(0)+x*(bk(1)+x*(bk(2)+x*(bk(3)+x*bk(4))));
// ratfun=p/q
//   END FUNCTION ratfun
double ratfun(double x, double ak[5], double bk[5]){
  double p= ak[0]+x*(ak[1]+x*(ak[2]+x*(ak[3]+x*ak[4])));
  double q= bk[0]+x*(bk[1]+x*(bk[2]+x*(bk[3]+x*bk[4])));
  return(p/q);
}
  
  // FUNCTION lambdaeta(eta)
  // ! lambdaeta is the positive number satisfying
  // ! eta^2/2=lambda-1-ln(lambda)
  // ! with sign(lambda-1)=sign(eta);
  // USE Someconstants
  //   IMPLICIT NONE
  //   REAL(r8) :: eta, lambdaeta, ak(6), q, r, s, L, la
  //   REAL(r8) :: L2, L3, L4, L5
  //   s=eta*eta*0.5_r8
  // IF (eta==0) THEN
  //   la=1
  // ELSEIF (eta < -1) THEN
  //   r=exp(-1-s);
  // ak(1)=1.0_r8;
  // ak(2)=1.0_r8;
  // ak(3)=3.0_r8/2.0_r8;
  // ak(4)=8.0_r8/3.0_r8;
  // ak(5)=125.0_r8/24.0_r8;
  // ak(6)=54.0_r8/5.0_r8;
  // la=r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
  //   ELSEIF (eta<1) THEN
  //   ak(1)= 1.0_r8;
  // ak(2)= 1.0_r8/3.0_r8;
  // ak(3)=1.0_r8/36.0_r8;
  // ak(4)= -1.0_r8/270.0_r8;
  // ak(5)= 1.0_r8/4320.0_r8;
  // ak(6)= 1.0_r8/17010.0_r8;
  // r=eta;
  // la=1+r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
  //   ELSE
  //   r=11+s; L=log(r); la=r+L; r=1.0_r8/r;
  //   L2=L*L
  //     L3=L2*L
  //     L4=L3*L
  //     L5=L4*L
  //     ak(1)= 1;
  //   ak(2)=(2-L)*0.5_r8;
  //   ak(3)=(-9*L+6+2*L2)/6.0_r8;
  //   ak(4)= -(3*L3+36*L-22*L2-12)/12.0_r8;
  //   ak(5)=(60+350*L2-300*L-125*L3+12*L4)/60.0_r8;
  //   ak(6)=-(-120-274*L4+900*L-1700*L2+1125*L3+20*L5)/120.0_r8;
  //   la=la+L*r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
  //     ENDIF
  //     r= 1;
  //   IF (((eta>-3.5).AND.(eta<-0.03)).OR.((eta>0.03).AND.(eta<40))) THEN
  //     r=1;
  //   q=la;
  //   DO WHILE (r > 1.0e-8_r8)
  //     la=q*(s+log(q))/(q-1.0_r8);
  //   r= abs(q/la-1);
  //   q= la
  //     ENDDO
  //     ENDIF
  //     lambdaeta=la
  //     END FUNCTION lambdaeta
  //     
double lambdaeta(double eta) {
  double q, r, s, L, la;
  NumericVector ak(6);
  double L2, L3, L4, L5;
  s=eta*eta*0.5;
  if(eta==0.0) {
    la = 1.0;
  } else if(eta<-1.0) {
    r=exp(-1-s);
    ak[0]=1.0;
    ak[1]=1.0;
    ak[2]=3.0/2.0;
    ak[3]=8.0/3.0;
    ak[4]=125.0/24.0;
    ak[5]=54.0/5.0;
    la=r*(ak[0]+r*(ak[1]+r*(ak[2]+r*(ak[3]+r*(ak[4]+r*ak[5])))));
  } else if(eta < 1.0) {
    ak[0]= 1.0;
    ak[1]= 1.0/3.0;
    ak[2]=1.0/36.0;
    ak[3]= -1.0/270.0;
    ak[4]= 1.0/4320.0;
    ak[5]= 1.0/17010.0;
    r=eta;
    la=1.0+r*(ak[0]+r*(ak[1]+r*(ak[2]+r*(ak[3]+r*(ak[4]+r*ak[5])))));
  } else {
    r=11.0+s; L=log(r); la=r+L; r=1.0/r;
    L2=L*L;
    L3=L2*L;
    L4=L3*L;
    L5=L4*L;
    ak[0]= 1.0;
    ak[1]=(2.0-L)*0.5;
    ak[2]=(-9.0*L+6.0+2.0*L2)/6.0;
    ak[3]= -(3.0*L3+36.0*L-22.0*L2-12.0)/12.0;
    ak[4]=(60.0+350.0*L2-300.0*L-125.0*L3+12.0*L4)/60.0;
    ak[5]=-(-120.0-274.0*L4+900.0*L-1700.0*L2+1125.0*L3+20.0*L5)/120.0;
    la=la+L*r*(ak[0]+r*(ak[1]+r*(ak[2]+r*(ak[3]+r*(ak[4]+r*ak[5])))));
  }
  r= 1.0;
  if(((eta>-3.5) && (eta<-0.03)) || ((eta>0.03) && (eta<40))) {
    r=1.0;
    q=la;
    while(r > 1.0e-8) {
      la = q*(s+log(q))/(q-1.0);
      r = std::abs(q/la-1.0);
      q = la;
    }
  }
  return(la);
}      
  
// FUNCTION eps1(eta)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: eps1, eta, la, ak(0:4), bk(0:4) 
//   IF (abs(eta)<1.0) THEN
//   ak(0)=-3.333333333438e-1_r8;  bk(0)= 1.000000000000e+0_r8;     
//   ak(1)=-2.070740359969e-1_r8;  bk(1)= 7.045554412463e-1_r8;     
//   ak(2)=-5.041806657154e-2_r8;  bk(2)= 2.118190062224e-1_r8;     
//   ak(3)=-4.923635739372e-3_r8;  bk(3)= 3.048648397436e-2_r8;     
//   ak(4)=-4.293658292782e-5_r8;  bk(4)= 1.605037988091e-3_r8;     
//   eps1=ratfun(eta,ak,bk)
//     ELSE
//     la=lambdaeta(eta);
//   eps1=log(eta/(la-1.0_r8))/eta
//     ENDIF
//     END FUNCTION eps1
double eps1(double eta) {
  double eps1, la;
  double ak[5], bk[5];
  if(std::abs(eta)<1.0) {
    ak[0]=-3.333333333438e-1;  bk[0]= 1.000000000000e+0;
    ak[1]=-2.070740359969e-1;  bk[1]= 7.045554412463e-1;
    ak[2]=-5.041806657154e-2;  bk[2]= 2.118190062224e-1;
    ak[3]=-4.923635739372e-3;  bk[3]= 3.048648397436e-2;
    ak[4]=-4.293658292782e-5;  bk[4]= 1.605037988091e-3;
    eps1=ratfun(eta,ak,bk);
  } else {
    la = lambdaeta(eta);
    eps1=log(eta/(la-1.0))/eta;
  }
  return(eps1);
}
     
//     FUNCTION eps2(eta)
//     USE Someconstants
//     IMPLICIT NONE
//     REAL(r8) :: eta, eps2, x, lnmeta, ak(0:4), bk(0:4)
//     IF (eta < -5.0) THEN
//     x=eta*eta;
//   lnmeta=log(-eta)
//     eps2=(12.0_r8-x-6.0*(lnmeta*lnmeta))/(12.0*x*eta)
//     ELSEIF (eta<-2.0) THEN
//     ak(0)=-1.72847633523e-2_r8;  bk(0)=1.00000000000e+0_r8;     
//     ak(1)= -1.59372646475e-2_r8;  bk(1)= 7.64050615669e-1_r8;     
//     ak(2)= -4.64910887221e-3_r8;  bk(2)= 2.97143406325e-1_r8;     
//     ak(3)= -6.06834887760e-4_r8;  bk(3)= 5.79490176079e-2_r8;     
//     ak(4)= -6.14830384279e-6_r8;  bk(4)= 5.74558524851e-3_r8;     
//     eps2= ratfun(eta,ak,bk)
//       ELSEIF (eta < 2.0) THEN
//       ak(0)=-1.72839517431e-2_r8;  bk(0)= 1.00000000000e+0_r8;     
//       ak(1)=-1.46362417966e-2_r8;  bk(1)= 6.90560400696e-1_r8;     
//       ak(2)=-3.57406772616e-3_r8;  bk(2)= 2.49962384741e-1_r8;     
//       ak(3)=-3.91032032692e-4_r8;  bk(3)= 4.43843438769e-2_r8;     
//       ak(4)=2.49634036069e-6_r8;   bk(4)= 4.24073217211e-3_r8;     
//       eps2= ratfun(eta,ak,bk)
//         ELSEIF (eta < 1000.0) THEN
//         ak(0)= 9.99944669480e-1_r8;  bk(0)= 1.00000000000e+0_r8;     
//         ak(1)= 1.04649839762e+2_r8;  bk(1)= 1.04526456943e+2_r8;     
//         ak(2)= 8.57204033806e+2_r8;  bk(2)= 8.23313447808e+2_r8;     
//         ak(3)= 7.31901559577e+2_r8;  bk(3)= 3.11993802124e+3_r8;     
//         ak(4)= 4.55174411671e+1_r8;  bk(4)= 3.97003311219e+3_r8;     
//         x=1.0_r8/eta;
//         eps2=ratfun(x,ak,bk)/(-12.0*eta)
//           ELSE
//           eps2=-1.0/(12.0*eta)
//           ENDIF
//           END FUNCTION
double eps2(double eta) {
  double eps2, x, lnmeta;
  double ak[5], bk[5];
  if(eta < -5.0) {
    x=eta*eta;
    lnmeta=log(-eta);
    eps2=(12.0-x-6.0*(lnmeta*lnmeta))/(12.0*x*eta);
  } else if(eta <-2.0) {
    ak[0]=-1.72847633523e-2;  bk[0]=1.00000000000e+0;
    ak[1]= -1.59372646475e-2;  bk[1]= 7.64050615669e-1;
    ak[2]= -4.64910887221e-3;  bk[2]= 2.97143406325e-1;
    ak[3]= -6.06834887760e-4;  bk[3]= 5.79490176079e-2;
    ak[4]= -6.14830384279e-6;  bk[4]= 5.74558524851e-3;
    eps2= ratfun(eta,ak,bk);
  } else if(eta <2.0) {
    ak[0]=-1.72839517431e-2;  bk[0]= 1.00000000000e+0;
    ak[1]=-1.46362417966e-2;  bk[1]= 6.90560400696e-1;
    ak[2]=-3.57406772616e-3;  bk[2]= 2.49962384741e-1;
    ak[3]=-3.91032032692e-4;  bk[3]= 4.43843438769e-2;
    ak[4]=2.49634036069e-6;   bk[4]= 4.24073217211e-3;
    eps2= ratfun(eta,ak,bk);
  } else if(eta <1000.0) {
    ak[0]= 9.99944669480e-1;  bk[0]= 1.00000000000e+0;
    ak[1]= 1.04649839762e+2;  bk[1]= 1.04526456943e+2;
    ak[2]= 8.57204033806e+2;  bk[2]= 8.23313447808e+2;
    ak[3]= 7.31901559577e+2;  bk[3]= 3.11993802124e+3;
    ak[4]= 4.55174411671e+1;  bk[4]= 3.97003311219e+3;
    x=1.0/eta;
    eps2=ratfun(x,ak,bk)/(-12.0*eta);
  } else {
    eps2=-1.0/(12.0*eta);
  }
  return(eps2);
}

           
//            FUNCTION eps3(eta)
//   USE Someconstants
//   IMPLICIT NONE
//   REAL(r8) :: eta, eps3, eta3, x, y, ak(0:4), bk(0:4)
//   IF (eta <-8.0) THEN
//   x=eta*eta
//   y=log(-eta)/eta;
// eps3=(-30.0+eta*y*(6.0_r8*x*y*y-12.0+x))/(12.0*eta*x*x)
//   ELSEIF (eta <-4.0) THEN
//   ak(0)= 4.95346498136e-2_r8;  bk(0)= 1.00000000000e+0_r8;
//   ak(1)= 2.99521337141e-2_r8;  bk(1)= 7.59803615283e-1_r8;
//   ak(2)= 6.88296911516e-3_r8;  bk(2)= 2.61547111595e-1_r8;
//   ak(3)= 5.12634846317e-4_r8;  bk(3)= 4.64854522477e-2_r8;
//   ak(4)= -2.01411722031e-5_r8; bk(4)= 4.03751193496e-3_r8;
//   eps3=ratfun(eta,ak,bk)/(eta*eta)
//     ELSEIF (eta <-2.0) THEN
//     ak(0)=4.52313583942e-3_r8;  bk(0)= 1.00000000000e+0_r8;
//     ak(1)=1.20744920113e-3_r8;  bk(1)= 9.12203410349e-1_r8;
//     ak(2)=-7.89724156582e-5_r8; bk(2)= 4.05368773071e-1_r8;
//     ak(3)=-5.04476066942e-5_r8; bk(3)= 9.01638932349e-2_r8;
//     ak(4)=-5.35770949796e-6_r8; bk(4)= 9.48935714996e-3_r8;
//     eps3=ratfun(eta,ak,bk)
//       ELSEIF  (eta < 2.0) THEN
//       ak(0)= 4.39937562904e-3_r8;  bk(0)= 1.00000000000e+0_r8;
//       ak(1)= 4.87225670639e-4_r8;  bk(1)= 7.94435257415e-1_r8;
//       ak(2)= -1.28470657374e-4_r8; bk(2)= 3.33094721709e-1_r8;
//       ak(3)= 5.29110969589e-6_r8;  bk(3)= 7.03527806143e-2_r8;
//       ak(4)= 1.57166771750e-7_r8;  bk(4)= 8.06110846078e-3_r8;
//       eps3= ratfun(eta,ak,bk)
//         ELSEIF (eta < 10.0) THEN
//         ak(0)= -1.14811912320e-3_r8;  bk(0)= 1.00000000000e+0_r8;
//         ak(1)= -1.12850923276e-1_r8;  bk(1)= 1.42482206905e+1_r8;
//         ak(2)= 1.51623048511e+0_r8;   bk(2)= 6.97360396285e+1_r8;
//         ak(3)= -2.18472031183e-1_r8;  bk(3)= 2.18938950816e+2_r8;
//         ak(4)= 7.30002451555e-2_r8;   bk(4)= 2.77067027185e+2_r8;
//         x= 1.0_r8/eta;
//         eps3= ratfun(x,ak,bk)/(eta*eta)
//           ELSEIF (eta < 100.0) THEN
//           ak(0)= -1.45727889667e-4_r8;  bk(0)= 1.00000000000e+0_r8;
//           ak(1)= -2.90806748131e-1_r8;  bk(1)= 1.39612587808e+2_r8;
//           ak(2)= -1.33085045450e+1_r8;  bk(2)= 2.18901116348e+3_r8;
//           ak(3)= 1.99722374056e+2_r8;   bk(3)= 7.11524019009e+3_r8;
//           ak(4)= -1.14311378756e+1_r8;  bk(4)= 4.55746081453e+4_r8;
//           x= 1.0_r8/eta;
//           eps3= ratfun(x,ak,bk)/(eta*eta)
//             ELSE
//             eta3=eta*eta*eta
//             eps3=-log(eta)/(12.0*eta3)
//             ENDIF
//             END FUNCTION eps3
double eps3(double eta) {
  double eps3, eta3, x, y;
  double ak[5], bk[5];
  if(eta<-8.0) {
    x=eta*eta;
    y=log(-eta)/eta;
    eps3=(-30.0+eta*y*(6.0*x*y*y-12.0+x))/(12.0*eta*x*x);
  }
  else if(eta<-4.0) {
    ak[0]= 4.95346498136e-2;  bk[0]= 1.00000000000e+0;
    ak[1]= 2.99521337141e-2;  bk[1]= 7.59803615283e-1;
    ak[2]= 6.88296911516e-3;  bk[2]= 2.61547111595e-1;
    ak[3]= 5.12634846317e-4;  bk[3]= 4.64854522477e-2;
    ak[4]= -2.01411722031e-5; bk[4]= 4.03751193496e-3;
    eps3=ratfun(eta,ak,bk)/(eta*eta);
  }
  else if(eta<-2.0) {
    ak[0]=4.52313583942e-3;  bk[0]= 1.00000000000e+0;
    ak[1]=1.20744920113e-3;  bk[1]= 9.12203410349e-1;
    ak[2]=-7.89724156582e-5; bk[2]= 4.05368773071e-1;
    ak[3]=-5.04476066942e-5; bk[3]= 9.01638932349e-2;
    ak[4]=-5.35770949796e-6; bk[4]= 9.48935714996e-3;
    eps3=ratfun(eta,ak,bk);
  }
  else if(eta<2.0) {
    ak[0]= 4.39937562904e-3;  bk[0]= 1.00000000000e+0;
    ak[1]= 4.87225670639e-4;  bk[1]= 7.94435257415e-1;
    ak[2]= -1.28470657374e-4; bk[2]= 3.33094721709e-1;
    ak[3]= 5.29110969589e-6;  bk[3]= 7.03527806143e-2;
    ak[4]= 1.57166771750e-7;  bk[4]= 8.06110846078e-3;
    eps3= ratfun(eta,ak,bk);
  }
  else if(eta < 10.0) {
    ak[0]= -1.14811912320e-3;  bk[0]= 1.00000000000e+0;
    ak[1]= -1.12850923276e-1;  bk[1]= 1.42482206905e+1;
    ak[2]= 1.51623048511e+0;   bk[2]= 6.97360396285e+1;
    ak[3]= -2.18472031183e-1;  bk[3]= 2.18938950816e+2;
    ak[4]= 7.30002451555e-2;   bk[4]= 2.77067027185e+2;
    x= 1.0/eta;
    eps3= ratfun(x,ak,bk)/(eta*eta);
  }
  else if(eta < 100.0) {
    ak[0]= -1.45727889667e-4;  bk[0]= 1.00000000000e+0;
    ak[1]= -2.90806748131e-1;  bk[1]= 1.39612587808e+2;
    ak[2]= -1.33085045450e+1;  bk[2]= 2.18901116348e+3;
    ak[3]= 1.99722374056e+2;   bk[3]= 7.11524019009e+3;
    ak[4]= -1.14311378756e+1;  bk[4]= 4.55746081453e+4;
    x= 1.0/eta;
    eps3= ratfun(x,ak,bk)/(eta*eta);
  }
  else {
    eta3=eta*eta*eta;
    eps3=-log(eta)/(12.0*eta3);
  }
  return(eps3);
}

// SUBROUTINE invincgam(a,p,q,xr,ierr)
//   ! -------------------------------------------------------------
//     ! invincgam computes xr in the equations P(a,xr)=p and Q(a,xr)=q
//     ! with a as a given positive parameter.
//     ! In most cases, we invert the equation with min(p,q)
//     ! -------------------------------------------------------------
//     ! Inputs:
//     !   a ,    argument of the functions
//     !   p,     function P(a,x)
//     !   q,     function Q(a,x)  
//     ! Outputs:
//     !   xr   , soluction of the equations P(a,xr)=p and Q(a,xr)=q
//     !          with a as a given positive parameter.
//     !   xini,  initial value for the Newton
//     !   nit  , number of iterations in the Newton
//     !   ierr , error flag
//     !          ierr=0,  computation succesful
//     !          ierr=-1, overflow problem in the computation of one of the 
//     !                   gamma factors before starting the Newton iteration.
//     !                   The initial approximation to the root is given
//     !                   as output.
//     !          ierr=-2, the number of iterations in the Newton method
//     !                   reached the upper limit N=15. The last value
//     !                   obtained for the root is given as output.
//     ! ------------------------------------------------------------------
//     USE Someconstants    
//       IMPLICIT NONE
//       REAL(r8), INTENT(IN) :: a
//       REAL(r8), INTENT(IN) :: p
//       REAL(r8), INTENT(IN) :: q  
//       REAL(r8), INTENT(OUT) :: xr
//       REAL(r8) :: porq, s, dlnr, logr, r, a2, a3, a4, ap1, ap12, ap13, ap14,&
//         ap2, ap22, x0, ck(1:5), b, eta, L, L2, L3, L4,&
//         b2, b3, x, x2, t, px, qx, y, fp
//       INTEGER,  INTENT(OUT) :: ierr
//       INTEGER :: n, m, ierrf
//       LOGICAL :: pcase
//       ierr=0
//     IF (p<0.5) THEN
//       pcase=.true.
//     porq=p
//       s=-1 
//     ELSE
//       pcase=.false.
//     porq=q 
//       s=1 
//     ENDIF
//       logr=(1.0_r8/a)*(log(p)+loggam(a+1.0_r8))
//       IF (logr <log(0.2*(1+a))) THEN
//       r=exp(logr)
//       m=0
//     a2=a*a
//       a3=a2*a
//       a4=a3*a
//       ap1=a+1.0_r8
//     ap12=(a+1.0_r8)*ap1
//       ap13=(a+1.0_r8)*ap12
//       ap14=ap12*ap12
//       ap2=a+2
//     ap22=ap2*ap2
//       ck(1)= 1.0_r8;
// ck(2)= 1.0_r8/(1.0_r8+a);
// ck(3)=0.5_r8*(3*a+5)/(ap12*(a+2));
// ck(4)= (1.0_r8/3.0_r8)*(31+8*a2+33*a)/(ap13*ap2*(a+3));
// ck(5)= (1.0_r8/24.0_r8)*(2888+1179*a3+125*a4+3971*a2+5661*a)/(ap14*ap22*(a+3)*(a+4)); 
// x0=r*(1+r*(ck(2)+r*(ck(3)+r*(ck(4)+r*ck(5)))))
//   ELSEIF ((q < min(0.02_r8,exp(-1.5*a)/gamma(a))).AND.(a<10)) THEN
//   m=0
// b=1.0_r8-a; 
// b2=b*b;
// b3=b2*b;
// eta=sqrt(-2/a*log(q*gamstar(a)*sqrttwopi/sqrt(a)));
// x0=a*lambdaeta(eta); 
// L=log(x0); 
// IF ((a>0.12).OR.(x0>5)) THEN
//   L2=L*L;
// L3=L2*L;
// L4=L3*L
//   r=1.0_r8/x0;
// ck(1)=L-1;
// ck(2)=(3*b-2*b*L+L2-2*L+2)/2.0_r8;
// ck(3)=(24*b*L-11*b2-24*b-6*L2+12*L-12-9*b*L2+6*b2*L+2*L3)/6.0_r8;
// ck(4)=(-12*b3*L+84*b*L2-114*b2*L+72+36*L2+3*L4-72*L+162*b-168*b*L-12*L3+25*b3&
//   -22*b*L3+36*b2*L2+120*b2)/12.0_r8;
//   x0=x0-L+b*r*(ck(1)+r*(ck(2)+r*(ck(3)+r*ck(4))))
//     ELSE
//     r=1.0_r8/x0;
//   L2=L*L;
//   ck(1)=L-1;
//   x0=x0-L+b*r*ck(1);    
//   ENDIF 
//     ELSEIF (abs(porq-0.5)< 1.0e-5_r8) THEN
//     m=0
//   x0=a-1.0_r8/3.0_r8+(8.0_r8/405.0_r8+184.0_r8/25515.0_r8/a)/a
//     ELSEIF (abs(a-1)<1.0e-4_r8) THEN
//     m=0
//   IF (pcase) THEN
//     x0=-log(1.0_r8-p) 
//     ELSE
//     x0=-log(q) 
//     ENDIF
//     ELSEIF (a<1.0) THEN
//     m=0 
//   IF (pcase) THEN 
//     x0=exp((1.0_r8/a)*(log(porq)+loggam(a+1.0_r8)))
//     ELSE
//     x0=exp((1.0_r8/a)*(log(1.0_r8-porq)+loggam(a+1.0_r8)))
//     ENDIF
//     ELSE 
//     m=1  
//   r=inverfc(2*porq);
//   eta=s*r/sqrt(a*0.5_r8);
//   eta=eta+(eps1(eta)+(eps2(eta)+eps3(eta)/a)/a)/a;
//   x0=a*lambdaeta(eta)
//     ENDIF
//     t=1;
//   x=x0; 
//   n=1;
//   a2=a*a
//     a3=a2*a
//     ! Implementation of the high order Newton-like method
//     DO WHILE ((t>1.0e-15_r8).AND.(n< 15))
//     x=x0;
//   x2=x*x
//     IF (m==0) THEN
//     dlnr=(1.0_r8-a)*log(x)+x+loggam(a);
//   IF (dlnr>log(giant)) THEN
//     n=20
//   ierr=-1
//   ELSE
//     r=exp(dlnr);
//   IF (pcase) THEN 
//     CALL  incgam(a,x,px,qx,ierrf)
//     ck(1)=-r*(px-p);  
//   ELSE
//     CALL  incgam(a,x,px,qx,ierrf)
//     ck(1)=r*(qx-q);
//   ENDIF   
//     ck(2)=(x-a+1.0_r8)/(2.0_r8*x);
//   ck(3)=(2*x2-4*x*a+4*x+2*a2-3*a+1)/(6*x2);
//   r=ck(1);
//   IF (a>0.1) THEN
//     x0=x+r*(1+r*(ck(2)+r*ck(3)));
//   ELSE  
//     IF (a>0.05) THEN
//     x0=x+r*(1+r*(ck(2)));
//   ELSE
//     x0=x+r;
//   ENDIF
//     ENDIF 
//     ENDIF       
//     ELSE
//     y=eta
//     fp=-sqrt(a/twopi)*exp(-0.5*a*y*y)/(gamstar(a));
//   r=-(1/fp)*x
//     IF (pcase) THEN 
//     CALL  incgam(a,x,px,qx,ierrf)
//     ck(1)=-r*(px-p);  
//   ELSE
//     CALL  incgam(a,x,px,qx,ierrf)
//     ck(1)=r*(qx-q);
//   ENDIF  
//     ck(2)=(x-a+1.0_r8)/(2.0_r8*x);
//   ck(3)=(2*x2-4*x*a+4*x+2*a2-3*a+1)/(6*x2);
//   r=ck(1);
//   IF (a>0.1) THEN
//     x0=x+r*(1+r*(ck(2)+r*ck(3)));
//   ELSE  
//     IF (a>0.05) THEN
//     x0=x+r*(1+r*(ck(2)));
//   ELSE
//     x0=x+r;
//   ENDIF
//     ENDIF 
//     ENDIF
//     t=abs(x/x0-1.0_r8);
//   n=n+1; 
//   x=x0
//     ENDDO
//     IF (n==15) ierr=-2
//   xr=x
//     END SUBROUTINE invincgam
// [[Rcpp::export(".invincgam")]]
double invincgam(double a, double p, double q) {
  double porq, s, dlnr, logr, r, a2, a3, a4, ap1, ap12, ap13, ap14;
  double ap2, ap22, x0, b, eta, L, L2, L3, L4;
  double b2, b3, x, x2, t, px, qx, y, fp;
  NumericVector ck(5); //ck(1:5)
  int n, m;
  bool pcase;
    
  // !          ierr=0,  computation succesful
  // !          ierr=-1, overflow problem in the computation of one of the 
  // !                   gamma factors before starting the Newton iteration.
  // !                   The initial approximation to the root is given
  // !                   as output.
  // !          ierr=-2, the number of iterations in the Newton method
  // !                   reached the upper limit N=15. The last value
  // !                   obtained for the root is given as output.
  if(p<0.5) {
    pcase=true;
    porq=p;
    s=-1.0;
  } else {
    pcase=false;
    porq=q; 
    s=1.0;
  }
  logr=(1.0/a)*(log(p)+lgamma(a+1.0)); //loggam in FORTRAN
  if(logr <log(0.2*(1.0+a))) {
    r=exp(logr);
    m=0;
    a2=a*a;
    a3=a2*a;
    a4=a3*a;
    ap1=a+1.0;
    ap12=(a+1.0)*ap1;
    ap13=(a+1.0)*ap12;
    ap14=ap12*ap12;
    ap2=a+2.0;
    ap22=ap2*ap2;
    ck[0]= 1.0;
    ck[1]= 1.0/(1.0+a);
    ck[2]=0.5*(3.0*a+5.0)/(ap12*(a+2.0));
    ck[3]= (1.0/3.0)*(31.0+8.0*a2+33.0*a)/(ap13*ap2*(a+3.0));
    ck[4]= (1.0/24.0)*(2888.0+1179.0*a3+125.0*a4+3971.0*a2+5661.0*a)/(ap14*ap22*(a+3.0)*(a+4.0)); 
    x0=r*(1.0+r*(ck[1]+r*(ck[2]+r*(ck[3]+r*ck[4]))));
  } else if((q < std::min(0.02,exp(-1.5*a)/tgamma(a))) && (a<10.0)) {
    m=0;
    b=1.0-a; 
    b2=b*b;
    b3=b2*b;
    eta=sqrt(-2.0/a*log(q*gamstar(a)*sqrttwopi/sqrt(a)));
    x0=a*lambdaeta(eta); 
    L=log(x0); 
    if((a>0.12) || (x0>5.0))  {
      L2=L*L;
      L3=L2*L;
      L4=L3*L;
      r=1.0/x0;
      ck[0]=L-1.0;
      ck[1]=(3.0*b-2.0*b*L+L2-2.0*L+2.0)/2.0;
      ck[2]=(24.0*b*L-11.0*b2-24.0*b-6.0*L2+12.0*L-12.0-9.0*b*L2+6.0*b2*L+2.0*L3)/6.0;
      ck[3]=(-12.0*b3*L+84.0*b*L2-114.0*b2*L+72.0+36.0*L2+3.0*L4-72.0*L+162.0*b-168.0*b*L-12.0*L3+25.0*b3-22.0*b*L3+36.0*b2*L2+120.0*b2)/12.0;
      x0=x0-L+b*r*(ck[0]+r*(ck[1]+r*(ck[2]+r*ck[3])));
    } else {
      r=1.0/x0;
      L2=L*L;
      ck[0]=L-1.0;
      x0=x0-L+b*r*ck[0];    
    }
  } else if(std::abs(porq-0.5)< 1.0e-5) {
    m=0;
    x0=a-1.0/3.0+(8.0/405.0+184.0/25515.0/a)/a;
  } else if(std::abs(a-1.0)<1.0e-4) {
    m=0;
    if(pcase) {
      x0=-log(1.0-p);   
    } else {
      x0=-log(q);
    }
  } else if(a<1.0) {
    m=0;
    if(pcase) {
      x0=exp((1.0/a)*(log(porq)+lgamma(a+1.0)));  //loggam in fortran
    } else {
      x0=exp((1.0/a)*(log(1.0-porq)+lgamma(a+1.0)));  //loggam in fortran
    }
  } else {
    m=1;
    r=inverfc(2.0*porq);
    eta=s*r/sqrt(a*0.5);
    eta=eta+(eps1(eta)+(eps2(eta)+eps3(eta)/a)/a)/a;
    x0=a*lambdaeta(eta);
  }
  t=1.0;
  x=x0; 
  n=1;
  a2=a*a;
  a3=a2*a;
  //Implementation of the high order Newton-like method
  while((t>1.0e-15) && (n< 15)) {
    x=x0;
    x2=x*x;
    if(m==0) {
      dlnr=(1.0-a)*log(x)+x+lgamma(a); //loggam in fortran
      if(dlnr > log(giant)) {
        n=20;
        warning("overflow problem in the computation of one of the gamma factors before starting the Newton iteration. The initial approximation to the root is given as output.");
      } else {
        r=exp(dlnr);
        if(pcase) {
          NumericVector pq = incgam(a,x);
          px = pq[0];
          qx = pq[1];
          ck[0]=-r*(px-p);
        } else {
          NumericVector pq = incgam(a,x);
          px = pq[0];
          qx = pq[1];
          ck[0]=r*(qx-q);
        }
        ck[1]=(x-a+1.0)/(2.0*x);
        ck[2]=(2.0*x2-4.0*x*a+4.0*x+2.0*a2-3.0*a+1.0)/(6.0*x2);
        r=ck[0];
        if(a>0.1) {
          x0=x+r*(1.0+r*(ck[1]+r*ck[2]));
        } else {
          if(a>0.05) {
            x0=x+r*(1.0+r*(ck[1]));
          } else {
            x0 = x + r;
          }
        }
      }
    }
    else {
      y=eta;
      fp=-sqrt(a/twopi)*exp(-0.5*a*y*y)/(gamstar(a));
      r=-(1.0/fp)*x;
      if(pcase) {
        NumericVector pq = incgam(a,x);
        px = pq[0];
        qx = pq[1];
        ck[0]=-r*(px-p);
      } else {
        NumericVector pq = incgam(a,x);
        px = pq[0];
        qx = pq[1];
        ck[0]=r*(qx-q);
      }
      ck[1]=(x-a+1.0)/(2.0*x);
      ck[2]=(2.0*x2-4.0*x*a+4.0*x+2.0*a2-3.0*a+1.0)/(6.0*x2);
      r=ck[0];
      if(a>0.1) {
        x0=x+r*(1.0+r*(ck[1]+r*ck[2]));
      } else {
        if(a>0.05) {
          x0=x+r*(1.0+r*(ck[1]));
        } else {
          x0=x+r;
        }
      }
    }
    t=std::abs(x/x0-1.0);
    n=n+1; 
    x=x0;
  }
  // if(n==15) {
  //   warning("The number of iterations in the Newton method reached the upper limit N=15."); 
  // }
  return(x);
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
  //    This code is distributed under the GNU-LGPL v.3 license,
  //    described in: 
  //               https://www.gnu.org/licenses/lgpl-3.0.txt
  //
  //  Modified:
  //
  //    22 January 2008
  //
  //  Author:
  //
  //    Original FORTRAN77 version by Chi Leung Lau.
  //    C++ version by John Burkardt available at:
  //           https://people.math.sc.edu/Burkardt/cpp_src/asa032/asa032.html
  //    Modified by Miquel De Caceres
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
    warning("x <= 0.0 in gammds");
    value = 0.0;
    return value;
  }
  
  if ( p <= 0.0 ) 
  {
    warning("p <= 0.0 in gammds");
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