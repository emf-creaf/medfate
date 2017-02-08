#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

const double pi = 3.141592;
const double g = 9.8; //m/s



/**
 * From Richards (1995). Eq. 14
 * phi - angle from ignition point
 * theta - angle of the maximum spread rate (virtual wind direction)
 * a - major axis
 * b - minor axis
 * c - distance from ignition point (i.e. focus) to center of the ellipse
 * 
 */
double ellros(double phi, double theta, double a, double b, double c){
  return(sqrt(pow(a,2.0)*pow(cos(phi-theta),2.0)+pow(b,2.0)*pow(sin(phi-theta),2.0))+ c*cos(phi-theta));
}
 // [[Rcpp::export(".genros")]]
double genros(double phi, double theta, double a1, double a2, double b, double n1, double n2, double c) {
  double pi = 3.141592;
  bool cond = (((pi/2.0)<= phi) &  (phi<= (pi*(3.0/2.0))));
  double a =  cond ? a2 : a1;
  double n = cond ? n2 : n1;
  double dif = phi-theta;
  double first = pow(std::abs(a*cos(dif)), 2.0/(2.0-n));
  double second = pow(std::abs(b*sin(dif)), 2.0/(2.0-n));
//  Rcout<<phi<< " , "<<a1<<", "<<a2<<", "  <<a << ", "<<first <<", "<<second<<"\n";
  return(pow(first+ second ,(2.0-n)/2.0)+c*cos(dif));
}
/**
 * Vector addition from polar coordinates (length, angles in radians)
 * Angles are measured from the y-axis (north)
 */
NumericVector vectorAddition(NumericVector v1, NumericVector v2) {
  //add coordinates
  double x = v1[0]*sin(v1[1])+v2[0]*sin(v2[1]);
  double y = v1[0]*cos(v1[1])+v2[0]*cos(v2[1]);
//  Rcout << x << " "<< y<<"\n";
  return(NumericVector::create(sqrt(pow(x,2.0)+pow(y,2.0)), atan2 (x, y)));
}

int getEllipseIndex(int angle) {
  if(angle==0) return(0);
  else if(angle==45) return(1);
  else if(angle==90) return(2);
  else if(angle==135) return(3);
  else if(angle==180) return(4);
  else if(angle==225) return(5);
  else if(angle==270) return(6);
  else if(angle==315) return(7);
  return(NA_INTEGER);
}

/**
 * From Finney (2004): FARSITE
*/
 // [[Rcpp::export(".ellipseROS")]]
NumericVector ellipseROS(NumericVector phi, double theta, double vws, double ros) {
  //From Alexander (1985)
  double LB = 0.936*exp(0.2566*vws)+0.461*exp(-0.1548*vws)-0.397;
  double k = sqrt(pow(LB,2.0)-1.0);
  double HB = (LB+k)/(LB-k);
  double a = 0.5*(ros+(ros/HB));
  double b = a/LB;
  double c = a - (ros/HB);
//  Rcout<< " LB " << LB <<" HB " << HB <<" a "<< a << " b "<< b << " c "<< c<<"\n";
  int n = phi.size();
  NumericVector R(n);
  double fc = (3.141592/4.0);
  theta = fc*round(theta/fc);
  int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i7 = 0, i8 = 0;
  for(int i=0;i<n;i++) {
    R[i] = ellros(phi[i], theta, a, b, c);
    int t = round((theta-phi[i])/fc);
    while(t<=0) t +=8;
    while(t>8) t -=8;
//    Rcout << theta << ", " << phi[i]<<", "<< (theta-phi[i])<<", "<<t<<"\n";
    if(t==1) i1 = i;
    else if(t==2) i2 = i;
    else if(t==3) i3 = i;
    else if(t==4) i4 = i;
    else if(t==5) i5 = i;
    else if(t==7) i7 = i;
    else if(t==8) i8 = i;
  }
//  Rcout<< "i1 "<< i1 << " i2 "<< i2 << " i3 "<< i3 << " i4 "<< i4;
//  stop("kk");
  double DF = R[i8];
  double CF = R[i2];
  double RnewF = sqrt(pow(CF,2.0) - (pow(CF,4.0)/(pow(CF,2.0)+pow(DF,2.0))));
  RnewF = std::max(RnewF,std::min(CF, DF));
  R[i1] = RnewF;
  R[i7] = RnewF;
  double DB = R[i4];
  double CB = R[i2];
  double RnewB = sqrt(pow(CB,2.0) - (pow(CB,4.0)/(pow(CB,2.0)+pow(DB,2.0))));
  RnewB = std::max(RnewB,std::min(CB, DB));
  R[i3] = RnewB;
  R[i5] = RnewB;
  return(R);
}

 // [[Rcpp::export(".doubleEllipseROS")]]
NumericVector doubleEllipseROS(NumericVector phi, double theta, double vws, double ros) {
  //From Alexander (1985)
  double LB = 0.936*exp(0.2566*vws)+0.461*exp(-0.1548*vws)-0.397;
  
  double n1 = 1.0;
  double n2 = 1.0;
  double aRatio = 1.0/2.0;
  double LB1 = (2.0*LB)/(1.0+(1.0/aRatio));
  double LB2 = LB1/aRatio;
  double k = sqrt(pow(LB2,2.0)-1.0);
  double HB = (LB1+k)/(LB2-k);
  double a1 = (ros+(ros/HB))/(1.0+(1.0/aRatio));
  double a2 = a1/aRatio;
  double b = (a1+a2)/(2.0*LB);
  double c = ros - a1;
  Rcout<< "k "<< k <<" HB " << HB <<" LB " << LB <<" LB1 " << LB1 <<" LB2 " << LB2 <<" a1 "<< a1 << " a2 " << a2 << " b "<< b << " c "<< c<<"\n";
  int n = phi.size();
  NumericVector R(n);
  double fc = (3.141592/4.0);
  theta = fc*round(theta/fc);
  int i1=0, i2=0, i3=0, i4=0, i5=0, i7=0, i8=0;
  for(int i=0;i<n;i++) {
    R[i] = genros(phi[i], theta, a1, a2, b, n1, n2, c);
    int t = round((theta-phi[i])/fc);
    while(t<=0) t +=8;
    while(t>8) t -=8;
//    Rcout << theta << ", " << phi[i]<<", "<< (theta-phi[i])<<", "<<t<<"\n";
    if(t==1) i1 = i;
    else if(t==2) i2 = i;
    else if(t==3) i3 = i;
    else if(t==4) i4 = i;
    else if(t==5) i5 = i;
    else if(t==7) i7 = i;
    else if(t==8) i8 = i;
  }
//  Rcout<< "i1 "<< i1 << " i2 "<< i2 << " i3 "<< i3 << " i4 "<< i4;
//  stop("kk");
  double DF = R[i8];
  double CF = R[i2];
  double RnewF = sqrt(pow(CF,2.0) - (pow(CF,4.0)/(pow(CF,2.0)+pow(DF,2.0))));
  RnewF = std::max(RnewF,std::min(CF, DF));
  R[i1] = RnewF;
  R[i7] = RnewF;
  double DB = R[i4];
  double CB = R[i2];
  double RnewB = sqrt(pow(CB,2.0) - (pow(CB,4.0)/(pow(CB,2.0)+pow(DB,2.0))));
  RnewB = std::max(RnewB,std::min(CB, DB));
  R[i3] = RnewB;
  R[i5] = RnewB;

  return(R);
}




// [[Rcpp::export(".fireBrandFallingHeight")]]
double fireBrandFallingHeight(double initialHeight, double timeFalling, double Dp) {
  double Cd = 1.2; //dimensionless
  double rho_s = 0.3; // density of charred wood cylinder g/cm3
  double rho_a = 0.0012; //density of air g/cm3
  double v0 = sqrt((pi*g*rho_s*Dp)/(2.0*Cd*rho_a));
  double tau = (4.0*Cd*v0)/(0.0064*pi*g);
  double z = initialHeight-v0*tau*((timeFalling/tau)-0.5*pow(timeFalling/tau,2.0));
  return(z);
}



/**
 * Albini's spotting model (1979). Description also in Finney (2004)
 * zF: flame height (m)
 * z0: initial height of the particle
 * Dp: diameter of the particle
 */
// [[Rcpp::export(".totalFirebrandLoftingTime")]]
double totalFirebrandLoftingTime(double z, double z0, double zF, double Dp) {
  double zRatioIni = z0/zF;
  double vRatio = 40.0*sqrt(Dp/zF);
  double t1 = 1.0 - sqrt(zRatioIni)+vRatio*log((1.0-vRatio)/(sqrt(zRatioIni)-vRatio));
  double t2 = 0.2 + vRatio*(1.0+vRatio*log(1.0+(1.0/(1.0-sqrt(Dp/zF)))));
  double r = sqrt((4.563+ z/zF)/5.963);
  double carro = log((1.0-0.8*vRatio)/(1.0-0.8*r*vRatio))-(0.8*vRatio*(r-1.0))-(0.5*pow(0.8*vRatio,2.0)*(pow(r,2.0)-1.0));
  double t3 = (5.963/pow(0.8*vRatio,3.0))*carro;
  return(t1+t2+t3);
}
// [[Rcpp::export(".totalGasFlowPersistenceTime")]]
double totalGasFlowPersistenceTime(double z, double t0, double zF) {
  return(t0+ 1.2 + (5.963/3.0)*(pow((4.563+ z/zF)/5.963,1.5)-1.0));
}
// [[Rcpp::export(".findFireBrandLoftedHeight")]]
double findFireBrandLoftedHeight(double t0, double z0, double zF,double Dp) {
  double z = z0;
  double zRatioFin = z0/zF;
  double diff = totalGasFlowPersistenceTime(z, t0, zF)-totalFirebrandLoftingTime(z ,z0, zF, Dp);
  while(diff>0.0) {
    zRatioFin = zRatioFin+0.1;
    z = zRatioFin*zF;
    diff = totalGasFlowPersistenceTime(z, t0, zF)-totalFirebrandLoftingTime(z ,z0, zF, Dp);
  }
  return(zRatioFin*zF);
}
// [[Rcpp::export(".willBurnWhenHitFloor")]]
bool willBurnWhenHitFloor(double zIni, double Dp) {
  double maxZ = 0.39*pow(10,5.0)*Dp;
  return(maxZ > zIni);
}

double numberOfTreesBurningTogether(double LAIc) {
  return((LAIc/0.05)/(100*3.1415*0.20*0.20)); //Calculates density per 10x10 m squares assuming trees of 15 cm diameter and LAI2BA = 0.05
}
// [[Rcpp::export(".fireBrandBurningTimeFromCanopyStructure")]]
double fireBrandBurningTimeFromCanopyStructure(double LAIc) {
  return(4.0*pow(numberOfTreesBurningTogether(LAIc),-0.2)); //= 4*N^(-1/5)
}
// [[Rcpp::export(".fireBrandFlameHeightFromCanopyStructure")]]
double fireBrandFlameHeightFromCanopyStructure(double crownLength, double LAIc) {
  return(3.0*crownLength*pow(numberOfTreesBurningTogether(LAIc),0.4)); //= 3*CL*N^(2/5)
}

/**
 * From van Wagner (1977) model
 * CBH - crown base height (in m)
 * M - Crown foliar moisture (in percent)
 */
// [[Rcpp::export(".criticalFirelineIntensity")]]
double criticalFirelineIntensity(double CBH, double M) {
  return(pow(0.010*(CBH)*(460.0+25.9*M),1.5));
}


