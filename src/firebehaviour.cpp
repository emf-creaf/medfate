// [[Rcpp::depends(meteoland)]]
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <numeric>
#include <math.h>
#include "forestutils.h"
#include "firebehaviour.h"
#include "firebehaviour_c.h"
#include "fuelstructure.h"
#include "fuelmoisture.h"
#include "windextinction_c.h"
#include <meteoland.h>
using namespace Rcpp;

/**
 * Vector addition from polar coordinates (length, angles in radians)
 * Angles are measured from the y-axis (north)
 */
NumericVector vectorAddition(NumericVector v1, NumericVector v2) {
  //add coordinates
  double x = v1[0]*sin(v1[1])+v2[0]*sin(v2[1]);
  double y = v1[0]*cos(v1[1])+v2[0]*cos(v2[1]);
  //  Rcout << x << " "<< y<<"\n";
  return(NumericVector::create(sqrt(pow(x,2.0)+pow(y,2.0)), atan2(x, y)));
}


/**
 * FCCS
 * 
 *  Default moisture, slope and windspeed values are benchmark conditions used 
 *  to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire
 *  
 *  Strata indices:  0 - Canopy, 1 - Shrub, 2- Herb, 3 - Woody, 4 - Litter
 */
//' Fire behaviour functions
//' 
//' Function \code{fire_FCCS()} implements a modification of the fire behavior models 
//' described for the Fuel Characteristics Classification System (FCCS) in Prichard et al. (2013). 
//' Function \code{fire_Rothermel()} implements Rothermel's (1972) fire behaviour 
//' model (modified from package 'Rothermel' (Giorgio Vacchiano, Davide Ascoli)).
//' 
//' @param FCCSpropsSI A data frame describing the properties of five fuel strata (canopy, shrub, herbs, dead woody and litter) returned by \code{\link{fuel_FCCS}}.
//' @param MliveSI Moisture of live fuels (in percent of dry weight) for canopy, shrub, and herb strata. Live moisture values are drawn from column \code{ActFCM} in \code{FCCSpropsSI} if available (see \code{\link{fuel_FCCS}}). Otherwise, moisture values supplied for \code{MliveSI} are used.
//' @param MdeadSI Moisture of dead fuels (in percent of dry weight) for canopy, shrub, herb, woody and litter strata.
//' @param slope Slope (in degrees).
//' @param windSpeedSI Wind speed (in m/s) at 20 ft (6 m) over vegetation (default 11 m/s = 40 km/h)
//' 
//' @details Default moisture, slope and windspeed values are benchmark conditions 
//' used to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire.
//' 
//' @return Both functions return list with fire behavior variables. 
//' 
//' In the case of \code{fire_FCCS}, the function returns the variables in three blocks (lists \code{SurfaceFire}, \code{CrownFire} and \code{FirePotentials}), and the values are:
//' \itemize{
//'   \item{\code{SurfaceFire$`midflame_WindSpeed [m/s]`}: Midflame wind speed in the surface fire.}
//'   \item{\code{SurfaceFire$phi_wind}: Spread rate modifier due to wind.}
//'   \item{\code{SurfaceFire$phi_slope}: Spread rate modifier due to slope.}
//'   \item{\code{SurfaceFire$`I_R_surf [kJ/m2/min]`}: Intensity of the surface fire reaction.}
//'   \item{\code{SurfaceFire$`I_R_litter [kJ/m2/min]`}: Intensity of the litter fire reaction.}
//'   \item{\code{SurfaceFire$`q_surf [kJ/m2]`}: Heat sink of the surface fire.}
//'   \item{\code{SurfaceFire$`q_litter [kJ/m2]`}: Heat sink of the litter fire.}
//'   \item{\code{SurfaceFire$xi_surf}: Propagating flux ratio of the surface fire.}
//'   \item{\code{SurfaceFire$xi_litter}: Propagating flux ratio of the litter fire.}
//'   \item{\code{SurfaceFire$`ROS_surf [m/min]`}: Spread rate of the surface fire(without accounting for faster spread in the litter layer).}
//'   \item{\code{SurfaceFire$`ROS_litter [m/min]`}: Spread rate of the litter fire.}
//'   \item{\code{SurfaceFire$`ROS_windslopecap [m/min]`}: Maximum surface fire spread rate according to wind speed.}
//'   \item{\code{SurfaceFire$`ROS [m/min]`}: Final spread rate of the surface fire.}
//'   \item{\code{SurfaceFire$`I_b [kW/m]`}: Fireline intensity of the surface fire.}
//'   \item{\code{SurfaceFire$`FL [m]`}: Flame length of the surface fire.}
//'   \item{\code{CrownFire$`I_R_canopy [kJ/m2/min]`}: Intensity of the canopy fire reaction.}
//'   \item{\code{CrownFire$`I_R_crown [kJ/m2/min]`}: Intensity of the crown fire reaction (adding surface and canopy reactions).}
//'   \item{\code{CrownFire$`q_canopy [kJ/m2]`}: Heat sink of the canopy fire.}
//'   \item{\code{CrownFire$`q_crown [kJ/m2]`}: Heat sink of the crown fire (adding surface and canopy heat sinks).}
//'   \item{\code{CrownFire$xi_surf}: Propagating flux ratio of the crown fire.}
//'   \item{\code{CrownFire$`canopy_WindSpeed [m/s]`}: Wind speed in the canopy fire (canopy top wind speed).}
//'   \item{\code{CrownFire$WAF}: Wind speed adjustment factor for crown fires.}
//'   \item{\code{CrownFire$`ROS [m/min]`}: Spread rate of the crown fire.}
//'   \item{\code{CrownFire$Ic_ratio}: Crown initiation ratio.}
//'   \item{\code{CrownFire$`I_b [kW/m]`}: Fireline intensity of the crown fire.}
//'   \item{\code{CrownFire$`FL [m]`}: Flame length of the crown fire.}
//'   \item{\code{FirePotentials$RP}: Surface fire reaction potential (\[0-9\]).}
//'   \item{\code{FirePotentials$SP}: Surface fire spread rate potential (\[0-9\]).}
//'   \item{\code{FirePotentials$FP}: Surface fire flame length potential (\[0-9\]).}
//'   \item{\code{FirePotentials$SFP}: Surface fire potential (\[0-9\]).}
//'   \item{\code{FirePotentials$IC}: Crown initiation potential (\[0-9\]).}
//'   \item{\code{FirePotentials$TC}: Crown-to-crown transmission potential (\[0-9\]).}
//'   \item{\code{FirePotentials$RC}: Crown fire spread rate potential (\[0-9\]).}
//'   \item{\code{FirePotentials$CFC}: Crown fire potential (\[0-9\]).}
//' }
//' 
//' @references
//' Albini, F. A. (1976). Computer-based models of wildland fire behavior: A users' manual. Ogden, UT: US Department of Agriculture, Forest Service, Intermountain Forest and Range Experiment Station.
//' 
//' Rothermel, R. C. 1972. A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Research Paper INT USA.
//' 
//' Prichard, S. J., D. V Sandberg, R. D. Ottmar, E. Eberhardt, A. Andreu, P. Eagle, and K. Swedin. 2013. Classification System Version 3.0: Technical Documentation.
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @note Default moisture, slope and windspeed values are benchmark conditions used to calculate fire potentials (Sandberg et al. 2007) and map vulnerability to fire.
//' 
//' @seealso \code{\link{fuel_FCCS}}
//' 
//' @examples
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Calculate fuel properties according to FCCS
//' fccs <- fuel_FCCS(exampleforest, SpParamsMED)
//' fccs
//'   
//' #Calculate fire behavior according to FCCS
//' fire_FCCS(fccs)
//'   
//'  
//' @name fire_behaviour
// [[Rcpp::export("fire_FCCS")]]
List FCCSbehaviour(DataFrame FCCSpropsSI,
                   NumericVector MliveSI = NumericVector::create(90, 90, 60), 
                   NumericVector MdeadSI = NumericVector::create(6, 6, 6, 6, 6), 
                   double slope = 0.0, double windSpeedSI = 11.0) {
  
  InternalFCCS internalFCCS;
  internalFCCS.w = Rcpp::as< std::vector<double> >(FCCSpropsSI["w"]);
  internalFCCS.cover = Rcpp::as< std::vector<double> >(FCCSpropsSI["cover"]);
  internalFCCS.hbc = Rcpp::as< std::vector<double> >(FCCSpropsSI["hbc"]);
  internalFCCS.htc = Rcpp::as< std::vector<double> >(FCCSpropsSI["htc"]);
  internalFCCS.habc = Rcpp::as< std::vector<double> >(FCCSpropsSI["habc"]);
  internalFCCS.hatc = Rcpp::as< std::vector<double> >(FCCSpropsSI["hatc"]);
  internalFCCS.delta = Rcpp::as< std::vector<double> >(FCCSpropsSI["delta"]);
  internalFCCS.rhob = Rcpp::as< std::vector<double> >(FCCSpropsSI["rhob"]);
  internalFCCS.rhop = Rcpp::as< std::vector<double> >(FCCSpropsSI["rhop"]);
  internalFCCS.PV = Rcpp::as< std::vector<double> >(FCCSpropsSI["PV"]);
  internalFCCS.beta = Rcpp::as< std::vector<double> >(FCCSpropsSI["beta"]);
  internalFCCS.betarel = Rcpp::as< std::vector<double> >(FCCSpropsSI["betarel"]);
  internalFCCS.etabetarel = Rcpp::as< std::vector<double> >(FCCSpropsSI["etabetarel"]);
  internalFCCS.sigma = Rcpp::as< std::vector<double> >(FCCSpropsSI["sigma"]);
  internalFCCS.pDead = Rcpp::as< std::vector<double> >(FCCSpropsSI["pDead"]);
  internalFCCS.FAI = Rcpp::as< std::vector<double> >(FCCSpropsSI["FAI"]);
  internalFCCS.h = Rcpp::as< std::vector<double> >(FCCSpropsSI["h"]);
  internalFCCS.RV = Rcpp::as< std::vector<double> >(FCCSpropsSI["RV"]);
  internalFCCS.MinFMC = Rcpp::as< std::vector<double> >(FCCSpropsSI["MinFMC"]);
  internalFCCS.MaxFMC = Rcpp::as< std::vector<double> >(FCCSpropsSI["MaxFMC"]);
  internalFCCS.ActFMC = Rcpp::as< std::vector<double> >(FCCSpropsSI["ActFMC"]);
  
  std::vector<double> MliveSIvec = Rcpp::as< std::vector<double> >(MliveSI);
  std::vector<double> MdeadSIvec = Rcpp::as< std::vector<double> >(MdeadSI);
  FCCSBehaviour_RESULT res;
  FCCSbehaviour_c(res,
                  internalFCCS,
                  MliveSIvec, 
                  MdeadSIvec, 
                  slope, windSpeedSI);
  return(copyFCCSBehaviour_Result_c(res));
}

/** ROTHERMEL
 *  Recodified from function 'ros' in package 'Rothermel' (Vacchiano & Ascoli)
 *  Fuel classes are: 1-hour, 10-hour, 100-hour, live herbs and live woody
 * 
 *  modeltype: 'S'(tatic) or 'D'(ynamic) 
 *  wSI: vector of fuel load (t/ha) for five fuel classes
 *  sSI: vector of surface-to-volume ratio (m2/m3) for five fuel classes
 *  delta: a value of fuel bed depth (cm)
 *  mx_dead: a value of dead fuel moisture of extinction (percent)
 *  hSI: a vector of heat content (kJ/kg) for five fuel classes
 *  mSI: a vector of percent moisture on a dry weight basis (percent) for five fuel classes
 *  u: a value of windspeed (m/s) at midflame height. 
 *  windDir: wind direction (in degrees from north). North means blowing from north to south
 *  slope: a value of site slope (in degrees). Negative values are possible.
 *  aspect: aspect (in degrees from north)
 * 
 */
//' @rdname fire_behaviour
//' 
//' @param modeltype 'S'(tatic) or 'D'(ynamic)
//' @param wSI A vector of fuel load (t/ha) for five fuel classes.
//' @param sSI A vector of surface-to-volume ratio (m2/m3) for five fuel classes.
//' @param delta A value of fuel bed depth (cm).
//' @param mx_dead A value of dead fuel moisture of extinction (percent).
//' @param hSI A vector of heat content (kJ/kg) for five fuel classes.
//' @param mSI A vector of percent moisture on a dry weight basis (percent) for five fuel classes.
//' @param u A value of windspeed (m/s) at midflame height.
//' @param windDir Wind direction (in degrees from north). North means blowing from north to south.
//' @param aspect Aspect (in degrees from north).
//' 
// [[Rcpp::export("fire_Rothermel")]]
List rothermel(String modeltype, NumericVector wSI, NumericVector sSI, double delta, double mx_dead,
                  NumericVector hSI, NumericVector mSI, double u, double windDir, double slope, double aspect) {
  //Rescale variables to units of rothermel model
  NumericVector m = mSI/100.0; //from percent to proportions
  NumericVector w = wSI*0.02048161; //from t/ha to lbs/ft2
  NumericVector s = sSI/3.281; //from m-1 to ft-1
  NumericVector h = hSI*0.429922614; //from kJ/kg to ?
  delta = delta*0.0328084; //from cm to feet
  mx_dead = mx_dead/100.0;//from percent to proportions
  u=u*196.8504; //from m/s to ft/min
  slope = tan(slope*180.0/M_PI)/100.0; //from degrees to proportions

  // transfer partially cured herbaceous fuels to dead
  if (modeltype=="D") {
    double kt = 0.0;
    if((m[3]>=0.3) && (m[3]<1.2)) kt = (1.20-m[3])/0.9;
    // weighting SAV from transferred cured herbaceous
    double f1=w[0]*s[0]/32.0;
    double f4 = w[3]*kt*s[3]/32.0;
    s[0]=(f1*s[0]+f4*s[3])/(f1+f4);        
    if((f1+f4)==0.0) s[0] = 0.0;
    w[0]=w[0]+w[3]*kt;
    w[3]=w[3]-w[3]*kt;    
  }
  
  //Constants 
  double rhop= 32.0; // = 513*0.0624279606 Scott and Burgan (2005)
  double st=0.0555;
  double se=0.01;
  
  //Area fractions
  NumericVector a = s*w;
  a = a/rhop;
  double a_dead=a[0]+a[1]+a[2];
  double a_live=a[3]+a[4];
  double a_tot=(a_dead+a_live);
  NumericVector f = NumericVector::create(a[0]/a_dead,a[1]/a_dead,a[2]/a_dead,a[3]/a_live,a[4]/a_live);
  if(a_live==0.0) {
    f[3] = 0.0;
    f[4] = 0.0;
  }
  if(a_dead==0.0) {
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
  }
  double f_dead=a_dead/a_tot;
  double f_live=a_live/a_tot;
  if(a_tot==0.0) {
    f_dead = 0.0;
    f_live = 0.0;
  }
  
  //net (weighted) fuel loadings
  NumericVector wn=w*(1.0-st); // Albini 1976
  double wn_dead=f[0]*wn[0]+f[1]*wn[1]+f[2]*wn[2];
  double wn_live=wn[3]+wn[4]; // corrects models w/ 2 live fuel classes  (undocumented)  
  //weighted fuel moisture
  double mf_dead=f[0]*m[0]+f[1]*m[1]+f[2]*m[2];
  double mf_live=f[3]*m[3]+f[4]*m[4];  
  //weighted SAV ratio
  double sigma_dead=f[0]*s[0]+f[1]*s[1]+f[2]*s[2];
  double sigma_live=f[3]*s[3]+f[4]*s[4];
  double sigma_tot=(f_dead*sigma_dead+f_live*sigma_live); //characteristic SAV
  
  //weighted heat content
  double h_dead=f[0]*h[0]+f[1]*h[1]+f[2]*h[2];
  double h_live=f[3]*h[3]+f[4]*h[4];
  //mean packing ratio for fuel complex
  double beta=(1.0/delta)*(w[0]/rhop+w[1]/rhop+w[2]/rhop+w[3]/rhop+w[4]/rhop);
  if(delta==0.0) beta = 0.0;
  
  //live fuel moisture of extinction
  double Wden = (w[3]*exp(-500.0/s[3])+w[4]*exp(-500.0/s[4]));
  double W=(w[0]*exp(-138.0/s[0])+w[1]*exp(-138.0/s[1])+w[2]*exp(-138.0/s[2]))/Wden;
  if(Wden==0.0) W = 0.0;
  double den = (w[0]*exp(-138/s[0])+w[1]*exp(-138.0/s[1])+w[2]*exp(-138.0/s[2]));
  double mfpd=(w[0]*m[0]*exp(-138.0/s[0])+w[1]*m[1]*exp(-138.0/s[1])+w[2]*m[2]*exp(-138.0/s[2]))/den;
  if(den==0.0) mfpd = 0.0;
  double mx_live=std::max(2.9*W*(1.0-mfpd/mx_dead)-0.226, mx_dead);
  if(Wden==0) mx_live = mx_dead;
//  Rcout<<mx_dead<<", "<< mfpd<<", "<< mx_live<<"\n";

  //damping coefficients
  double ns=0.174*pow(se,-0.19);
  if(se==0.0) ns = 0.0;
  double nm_dead=1-2.59*(mf_dead/mx_dead)+5.11*pow(mf_dead/mx_dead,2.0)-3.52*pow(mf_dead/mx_dead,3.0);
  double nm_live=1-2.59*(mf_live/mx_live)+5.11*pow(mf_live/mx_live,2.0)-3.52*pow(mf_live/mx_live,3.0);
  if(mx_dead==0.0) nm_dead = 0.0;
  if(mx_live==0.0) nm_live = 0.0;
  if(mf_dead>mx_dead) nm_dead = 0.0;  
  if(mf_live>mx_live) nm_live = 0.0;  
//  Rcout<<nm_dead<<", "<< nm_live<<"\n";
  
  // Andrews 2013 pag.E
  //optimum packing ratio
  double beta_op=3.348*pow(sigma_tot,-0.8189);
  if(sigma_tot==0.0) beta_op = 0.0;
//  Rcout<<beta_op<<"\n";

  //relative packing ratio
  double rpr=beta/beta_op; 
  if(beta_op==0.0) rpr = 0.0;
//  Rcout<<rpr<<"\n";

  //maximum reaction velocity
  double gamma_max=pow(sigma_tot,1.5)/(495.0+0.0594*pow(sigma_tot,1.5));
  if(sigma_tot==0.0) gamma_max = 0.0;
//  Rcout<<gamma_max<<"\n";
  
  //reaction intensity (in kJ/ft2/min)
  double sum_dead=(wn_dead*h_dead*nm_dead*ns);
  double sum_live=(wn_live*h_live*nm_live*ns);
  //alternate formulation from Albini (1976)
  double A=133.0*pow(sigma_tot,-0.7913);
  if(sigma_tot==0.0) A = 0.0;
  double ir_dead=gamma_max*pow(rpr*exp(1.0-rpr), A)*sum_dead; //*f.dead removed by Frandsen 73
  double ir_live=gamma_max*pow(rpr*exp(1.0-rpr), A)*sum_live; //*f.live removed by Frandsen 73
  if(rpr==0.0) {
    ir_dead=0.0;
    ir_live=0.0;
  }
  double ir=ir_dead+ir_live;
//  Rcout<<sum_dead<<", "<<sum_live<<"\n";
//  Rcout<<ir<<"\n";
  
  //propagating flux ratio
  double xi=pow(192+0.2595*sigma_tot,-1.0)*exp((0.792+0.681*sqrt(sigma_tot))*(beta+0.1));  
//  Rcout<<xi<<"\n";

//wind coefficient  
  double C=7.47*exp(-0.133*pow(sigma_tot,0.55));
  if(sigma_tot==0.0) C = 0.0;
  double B=0.02526*pow(sigma_tot,0.54);
  if(sigma_tot==0.0) B = 0.0;
  double E=0.715*exp(-0.000359*sigma_tot);
  double fw=C*pow(u,B)*pow(rpr,-E);
  if(rpr==0.0) fw = 0.0;
  
  //slope coefficient
  double fs=5.275*pow(beta,-0.3)*pow(slope,2.0);
  if(beta==0.0) fs = 0.0;
  
  //calculate 'effective wind' direction and spread
  NumericVector vS = NumericVector::create(fs, (aspect*M_PI/180.0)+M_PI); //add 180 degrees for pushing fire in the direction opposite of the aspect
  NumericVector vW = NumericVector::create(fw, (windDir*M_PI/180.0)+M_PI); //add 180 degrees for blowing 'in the direction of'
  NumericVector vEff = vectorAddition(vS,vW);
//  Rcout << " fw " << fw << " dir " << (windDir+3.141592)*(180.0/3.141592) << " fs " << fs << " dir " << (slopeDir+3.141592)*(180.0/3.141592) << "res: " << vEff[0] <<","<<vEff[1]*(180.0/3.141592)<<"\n";
  
  //resulting effect and virtual wind speed
  double phi = vEff[0];
  double vws = pow(phi/(C*pow(rpr,-E)), 1.0/B);
  if((B==0.0) || (rpr==0.0)) vws = 0.0;
  vws = vws/196.8504; //from ft/min to m/s
//  Rcout<<vws<<"\n";

  //for heat sink (denominator ROS)
  double rhob=(w[0]+w[1]+w[2]+w[3]+w[4])/delta; 
  if(delta==0.0) rhob = 0.0;

  NumericVector qig=250+1116*m;
  for(int i=0;i<5;i++) if(w[i]==0.0) qig[i]=0.0;
  
  double eps=f_dead*(f[0]*qig[0]*exp(-138.0/s[0])+f[1]*qig[1]*exp(-138.0/s[1])+f[2]*qig[2]*exp(-138.0/s[2]))+f_live*(f[3]*qig[3]*exp(-138.0/s[3])+f[4]*qig[4]*exp(-138.0/s[4]));
//  Rcout<<eps<<"\n";

  //ROS
  double r = (ir*xi*(1+phi))/(rhob*eps);
  if(rhob==0.0) r = 0.0;
  if(eps==0.0) r = 0.0;
//  Rcout<<r<<"\n";  
  r = 0.3048*r; //from feet/min to m/min

  //Limit the rate of spread to be smaller or equal the effective wind speed
  //Andrews, P.L., Cruz, M.G., Rothermel, R.C., 2013. Examination of the wind speed limit function in the Rothermel surface fire spread model. Int. J. Wildl. Fire 22, 959–969. doi:10.1071/WF12122
  r = std::min(r, vws*60.0);
  
  
  double tr = 1259.843/(sigma_tot*3.281); //Residence time (in min) Anderson (1969)
  double Ib = (ir*0.1893)*tr*r;//Fireline intensity (kW/m)
  
  //return values for rothermel function
  List output=List::create(Named("Characteristic dead fuel moisture [%]")=mf_dead*100.0,
                           Named("Characteristic live fuel moisture [%]")=mf_live*100.0,
                           Named("Live fuel moisture of extinction [%]")=mx_live*100.0,
                           Named("Characteristic SAV [m2/m3]")=sigma_tot*3.281,
                           Named("Bulk density [kg/m3]") = rhob*16.0184634,
                           Named("Packing ratio [dimensionless]")=beta,
                           Named("Relative packing ratio [dimensionless]")=rpr,
                           Named("Dead fuel Reaction intensity [kW/m2]") = ir_dead* 0.1893,
                           Named("Live fuel Reaction intensity [kW/m2]") = ir_live* 0.1893,
                           Named("Reaction intensity [kW/m2]") = ir*0.1893,
                           Named("Fireline intensity [kW/m]") = Ib,
                           Named("Wind factor [dimensionless]") = fw,
                           Named("Slope factor [dimensionless]") = fs,
                           Named("Slope-wind vector") = vEff,
                           Named("Virtual wind speed [m/s]") = vws,
                           Named("Heat source [kW/m2]") = (ir*xi*(1+phi))*0.1893,
                           Named("Heat sink [kJ/m3]") = (rhob*eps)*37.25894580781,
                           Named("ROS [m/min]") = r);
  return(output);
}


/** FIRE BEHAVIOUR
 * 
 *  x: An object of class 'forest'
 *  soil: An object of class 'soil'
 *  slope: A value of site slope (in degrees). Negative values are possible.
 *  aspect: aspect (in degrees from north)
 *  meteo: A data frame with meteorological variables for each day.
 *  SpParams: A data frame with species parameters
 *  surfaceToVolumeRatios: A vector of surface-to-volume ratio (m2/m3) for five fuel classes
 *  heatContent: A vector of heat content (kJ/kg) for five fuel classes
 *  deadFuelMoistureExtinction: a value of dead fuel moisture of extinction (percent)
 *  control: A list with default parameter values
 */
// List fb(List x, DataFrame soil, double latitude, double slope, double aspect, DataFrame meteo, DataFrame SpParams, 
//         DataFrame FuelModelParams, List control) {
//   bool verbose = control["verbose"];
//   String liveFMCmode = control["liveFMCmode"];
//   bool useModelForLive = control["useModelForLive"]; 
//   if(verbose) Rcout<<"Initializing";
//   
//   //Prepare swb input
//   DataFrame swbInput = forest2swbInput(x, soil, SpParams, control);
//   NumericVector pBole = swbInput["pBole"];
//   NumericVector H = swbInput["H"];
//   NumericVector Sgdd = swbInput["Sgdd"];
//   IntegerVector SP = swbInput["SP"];
//   NumericVector LAI = swbInput["LAI"];
//   
//   NumericVector cohLoading = cohortFuelLoading(x, SpParams);
// 
//   //Initialize other cohort-based variables
//   int numCohorts = LAI.size();
//   NumericVector LAIphe(numCohorts), cohFMC(numCohorts), cohLoadingPhe(numCohorts);
//   
//   //Meteorological variables
//   IntegerVector DOY = meteo["DOY"]; // day of the year
//   NumericVector Precipitation = meteo["Precipitation"]; // in mm of water (L/m2)
//   NumericVector MeanTemperature = meteo["MeanTemperature"]; // in degrees Celsius
//   NumericVector MinTemperature = meteo["MinTemperature"]; // in degrees Celsius
//   NumericVector MaxTemperature = meteo["MaxTemperature"]; // in degrees Celsius
//   NumericVector RHmean = meteo["MeanRelativeHumidity"]; // in percentage
//   NumericVector RHmin = meteo["MinRelativeHumidity"]; // in percentage
//   NumericVector RHmax = meteo["MaxRelativeHumidity"]; // in percentage
//   NumericVector WS = meteo["WindSpeed"]; // in m/s
//   NumericVector WD = meteo["WindDirection"]; // in radians from north
//   NumericVector Rn = meteo["Radiation"];  //in MJ/m2
//   NumericVector PET = meteo["PET"]; //in mm of water
// 
//   NumericVector GDD = gdd(DOY, MeanTemperature, 5.0);
//   NumericVector ER = er(DOY);
// 
//   int numDays = Precipitation.size();
//   
//   //Output FS variables
//   NumericVector fuelbedHeight(numDays), fine1hLoading(numDays), 
//                 coarse10hLoading(numDays), coarse100hLoading(numDays), 
//                 herbaceousLoading(numDays), woodyLoading(numDays),
//                 canopyBaseHeight(numDays), canopyTopHeight(numDays),
//                 canopyLength(numDays), canopyBulkDensity(numDays),
//                 canopyLAI(numDays);
// 
//   //Output FM variables
//   NumericVector fuelWind(numDays), fine1hFMC(numDays), coarse10hFMC(numDays), 
//                 coarse100hFMC(numDays), canopyFMC(numDays), fuelbedWoodyFMC(numDays),
//                 fuelbedHerbaceousFMC(numDays);
// 
//   //Output FB variables
//   IntegerVector burningType(numDays);
//   NumericVector firelineIntensity(numDays), spreadRate(numDays), surfaceSpreadRate(numDays), 
//                 midflameWind(numDays), virtualWindSpeed(numDays), spreadDirection(numDays), I0(numDays);
//   
//   //Initialize fuel moisture
//   NumericVector fMoisture(5);
//   fMoisture[0] = x["dead1hourMoisture"];
//   fMoisture[1] = x["dead10hourMoisture"];
//   fMoisture[2] = x["dead100hourMoisture"];
//   
//   NumericVector fcafternoonprev; //Stores fuel conditions of afternoon of previous day
//   if(verbose) Rcout<<" - Daily loop ";
//   for(int i=0;i<numDays;i++) {
//     //SOIL WATER BALANCE
//     if(verbose) Rcout<<".";
//     List s = swbDay1(swbInput, soil, GDD[i], PET[i], Precipitation[i], ER[i], 0.0); //No Runon in simulations for a single cell    
//     double Lground = s["Lground"]; //Percent light at the ground level
// 
//     //FUEL STRUCTURE
//     NumericVector Phe = leafDevelopmentStatus(Sgdd, GDD[i]);
//     for(int c=0;c<numCohorts;c++) {
//       if(Sgdd[c]==0.0) Phe[c]=1.0;
//       LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology 
//       cohLoadingPhe[c] = cohLoading[c]*Phe[c]; //LAI modified by phenology 
//     }
//     //Loadings and FMC by cohort
//     cohFMC = cohortFuelMoistureContent(s, swbInput, SpParams);
//     
//     List fs = fuelStructure(x, SpParams, FuelModelParams, GDD[i], 30, 5000, 0.01, useModelForLive);
//     fuelbedHeight[i] = fs["fuelbedHeight [cm]"];
//     NumericVector fBiomass = fs["fuelbedLoading [Mg/ha]"];
//     NumericVector surfaceToVolumeRatios = fs["fuelbedSAV [m2/m3]"];
//     fine1hLoading[i] = fBiomass[0];
//     coarse10hLoading[i] = fBiomass[1];
//     coarse100hLoading[i] = fBiomass[2];
//     herbaceousLoading[i] = fBiomass[3];
//     woodyLoading[i] = fBiomass[4];
//     canopyBaseHeight[i] = fs["canopyBaseHeight [cm]"];
//     canopyTopHeight[i] = fs["canopyTopHeight [cm]"];
//     canopyLength[i] = fs["canopyLength [cm]"];
//     canopyBulkDensity[i] = fs["canopyBulkDensity [kg/m3]"];
//     canopyLAI[i] = fs["canopyLAI [dimensionless]"];
// 
//     //WIND EXTINCTION
//     //Check NA
//     if(NumericVector::is_na(WS[i])) WS[i] = 2.0; //Default value
//     if(NumericVector::is_na(WD[i])) WD[i] = 180.0; //Wind from south
//     
//     //Calculates midflameWind
//     if(NumericVector::is_na(canopyBaseHeight[i])) { //NO canopy over fire
//        midflameWind[i] = unshelteredMidflameWindSpeed_c(WS[i], fuelbedHeight[i]/100.0); 
//        fuelWind[i] = windSpeedAtCanopyHeight_c(WS[i], fuelbedHeight[i]/100.0);
//     } else { //Wind extinction due to crowns
//        midflameWind[i] = windSpeedMassmanExtinction_c(canopyBaseHeight[i]/100.0,WS[i], canopyLAI[i], canopyTopHeight[i]/100.0);
//        fuelWind[i] = midflameWind[i];
//     }    
// 
//     //FUEL MOISTURE
//     double netPrec = s["NetPrec"];
//     double rainDuration = std::min(24.0,Precipitation[i]/6.35);
//     
//     //Number of hours of insolation
//     NumericVector srs = meteoland::radiation_sunRiseSet(latitude*(M_PI/180.0),  slope*(M_PI/180.0), aspect*(M_PI/180.0),DOY[i]);
//     double n = std::max(0.0, (srs[1] - srs[0])*12.0/M_PI); 
//     if(n<=0.0) n = 0.0;
//     
//     //Dead fuel moisture
//     double radSec = Rn[i]*(Lground/100.0)*(1000000.0/(n*3600.0)); //Transform to MJ/m2.day to W/m2 assuming n hours of equal sunlight
//     if(n<=0.0) radSec = 0.0;
// 
//     NumericVector fcmean = fuelConditions(MeanTemperature[i],RHmean[i], radSec, fuelWind[i]);
//     NumericVector fcmorning = fuelConditions(MinTemperature[i],RHmax[i], radSec, fuelWind[i]);
//     NumericVector fcafternoon = fuelConditions(MaxTemperature[i],RHmin[i], radSec, fuelWind[i]);
//     if(i==0) fcafternoonprev = fcafternoon; 
//     fine1hFMC[i] = fine1hday(fMoisture[0], fcmean[0], fcmean[1], fuelWind[i], netPrec);
//     coarse10hFMC[i] = coarse10hday(fMoisture[1], 
//                                    (fcafternoonprev[0] +fcmorning[0])/2.0, (fcafternoonprev[1] + fcmorning[1])/2.0,
//                                    (fcmorning[0]+fcafternoon[0])/2.0, (fcmorning[1]+fcafternoon[1])/2.0, rainDuration);
//     coarse100hFMC[i] = coarse100hday(fMoisture[2], 
//                                      fcmorning[0], fcmorning[1],
//                                      fcafternoon[0], fcafternoon[1], n, rainDuration);
//     //Store for next day (and for fire behaviour)
//     fMoisture[0] = fine1hFMC[i]; 
//     fMoisture[1] =  coarse10hFMC[i]; 
//     fMoisture[2] = coarse100hFMC[i]; 
//     fcafternoonprev = fcafternoon;
//       
//     //Live fuel moisture
//     NumericVector psi = s["psiVec"];
//     if(liveFMCmode == "constant") {
//       fMoisture[3] = 200.0; 
//       fMoisture[4] = 100.0; 
//       canopyFMC[i] = 100.0;
//     } else if(liveFMCmode =="ffmc") {
//       fMoisture[3] = fine1hFMC[i]; 
//       fMoisture[4] = fine1hFMC[i]; 
//       canopyFMC[i] = fine1hFMC[i];
//     } else if(liveFMCmode =="swb") {
//       //Level of physiological activity is defined in terms of growth degree days
//       double act = leafDevelopmentStatus_c(200.0, GDD[i]);
//       //In winter, herbaceous FMC is modelled as fine (1h) dead fuels
//       //When physiologically active, FMC is modelled as function of water content in the topsoil
//       fMoisture[3] = act*std::max(30.0, 200.0+20*(psi[0]/1000.0))+(1.0-act)*fine1hFMC[i]; 
//       
//       fMoisture[4] = fuelbedLiveFuelMoisture(fuelbedHeight[i], cohFMC, cohLoadingPhe, H, pBole); 
//       canopyFMC[i] = canopyLiveFuelMoisture(canopyBaseHeight[i], canopyTopHeight[i], cohFMC, cohLoadingPhe, H, pBole); 
//     } else {
//       stop("Wrong 'liveFMCmode'");
//     }
//     
//     fuelbedHerbaceousFMC[i] = fMoisture[3];
//     fuelbedWoodyFMC[i] = fMoisture[4];
// 
// 
//    CharacterVector FMcode = x["FuelModelCode"];
//    CharacterVector models = FuelModelParams.attr("row.names");
//    NumericVector heatContent(5);
//    double deadFuelMoistureExtinction = 30.0; //Default value
//    for(int j=0;j<models.length();j++) {
//      if(models[j]==FMcode[0]) {
//        heatContent[0] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_1h"])[j];
//        heatContent[1] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_10h"])[j];
//        heatContent[2] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_100h"])[j];
//        heatContent[3] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Herb"])[j];
//        heatContent[4] = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Heat_Live_Woody"])[j];
//        deadFuelMoistureExtinction = Rcpp::as<Rcpp::NumericVector>(FuelModelParams["Mx_dead"])[j];
//      }
//    }
//    
//     //FIRE BEHAVIOUR
//     //1. Calculate surface fire characteristics
//     //1.2 Call to Rothermel model to obtain surface fire spread rate and intensity
//     List rtm = rothermel("D", fBiomass, surfaceToVolumeRatios, fuelbedHeight[i],
//                         deadFuelMoistureExtinction, heatContent, fMoisture, 
//                         midflameWind[i], WD[i], 
//                         slope, aspect);
//     surfaceSpreadRate[i] = rtm["ROS [m/min]"];
//     spreadRate[i] = surfaceSpreadRate[i];
//     firelineIntensity[i] = rtm["Fireline intensity [kW/m]"];
//     virtualWindSpeed[i] = rtm["Virtual wind speed [m/s]"];
//     NumericVector slopeWindVec = rtm["Slope-wind vector"];
//     
//     //1.3 Determine whether crowning will occur  (Van Wagner 1989)
//     if(!NumericVector::is_na(canopyBaseHeight[i])) {
//       I0[i] = criticalFirelineIntensity_c(canopyBaseHeight[i]/100.0, canopyFMC[i]); //critical fireline intensity (MW/m)
//       //Is crown fire (active or passive)?
//       if(firelineIntensity[i] >= I0[i]) { 
//         double RAC = 3.0/canopyBulkDensity[i]; // Threshold for active crown spread rate (m/min)
//         double Ei = 1.0;
//         NumericVector fBiomass10 = NumericVector::create(0.38, 0.092, 0.230,0.0, 0.092); //lb/square feet
//         fBiomass10 = fBiomass10*(1.0/(0.3048*0.3048))*0.45359237*(10.0); // transform to tons/hectare
//         List rtm10 = rothermel("D", fBiomass10, surfaceToVolumeRatios, 30.48,
//                         25, heatContent, fMoisture, WS[i]*0.4, WD[i], 
//                         slope, aspect);
//         double R10 = rtm10["ROS [m/min]"]; //Spread rate with fuel model 10 (Rothermel 1991)
//         double R0 = (I0[i]/firelineIntensity[i])*surfaceSpreadRate[i]; //critical surface spread rate (m/s) associated with critical fireline intensity
//         double ac = -log(0.1)/(0.9*(RAC-R0)); //Scaling coefficient
//         double cFractionBurn = 1.0 - exp(-ac*(surfaceSpreadRate[i]-R0)); //Fraction of canopy burn with crowning
//         double Rcmax = (3.34)*R10*Ei; // Maximum crown spread rate (m/s)
//         double Rcactual = surfaceSpreadRate[i] + cFractionBurn*(Rcmax - surfaceSpreadRate[i]); //Actual crown spread rate (m/s)
//         if(Rcactual >= RAC) {
//           burningType[i] = 4;
//           spreadRate[i] = Rcactual;
//         } else {
//           burningType[i] = 3;
//         }
//         //Modify fireline intensity for crown fire
//         firelineIntensity[i] = 300.0 * ((firelineIntensity[i]/(300.0*surfaceSpreadRate[i])) + cFractionBurn*canopyBulkDensity[i]*(canopyTopHeight[i]-canopyBaseHeight[i])*0.01)*spreadRate[i];        
//       } else {
//         burningType[i] = 2; // Surface fire
//       }
//     } else {
//       burningType[i] = 0; // Shrubland (open) fire
//     }
// 
//     spreadDirection[i] = slopeWindVec[1]*(180.0/3.141592);  
//     if(spreadDirection[i]<0) spreadDirection[i] +=360.0; 
//   }
//   
//   if(verbose) Rcout<<"Building output ...";  
//   Rcpp::DataFrame FS = DataFrame::create(_["fuelbedHeight"] = fuelbedHeight,
//                       _["fine1hLoading"] = fine1hLoading,
//                       _["coarse10hLoading"] = coarse10hLoading,
//                       _["coarse100hLoading"] = coarse100hLoading,
//                       _["herbaceousLoading"] = herbaceousLoading,
//                       _["woodyLoading"] = woodyLoading,
//                       _["canopyBaseHeight"] = canopyBaseHeight,
//                       _["canopyTopHeight"] = canopyTopHeight,
//                       _["canopyLength"] = canopyLength,
//                       _["canopyBulkDensity"] = canopyBulkDensity,
//                       _["canopyLAI"] = canopyLAI);
//   Rcpp::DataFrame FM = DataFrame::create(_["fuelWind"]=fuelWind, 
//                                          _["fine1hFMC"]=fine1hFMC, 
//                                          _["coarse10hFMC"]=coarse10hFMC,
//                                          _["coarse100hFMC"]=coarse100hFMC,
//                                          _["liveHerbaceousFMC"] = fuelbedHerbaceousFMC,
//                                          _["liveWoodyFMC"] = fuelbedWoodyFMC,
//                                          _["canopyFMC"] = canopyFMC);
// 
//   Rcpp::DataFrame FB = DataFrame::create(_["burningType"]=burningType,
//                                          _["midflameWind"] = midflameWind,
//                                          _["virtualWindSpeed"]=virtualWindSpeed,
//                                          _["spreadDirection"]=spreadDirection,
//                                          _["surfaceSpreadRate"]=surfaceSpreadRate,
//                                          _["intensityCrowningThreshold"]=I0,
//                                          _["spreadRate"]=spreadRate,
//                                          _["firelineIntensity"]=firelineIntensity);
//   
//   FS.attr("row.names") = meteo.attr("row.names");
//   FM.attr("row.names") = meteo.attr("row.names");
//   FB.attr("row.names") = meteo.attr("row.names");
//   List l = List::create(Named("FuelStructure")=FS, Named("FuelMoisture")=FM, Named("Behaviour")=FB);
//   l.attr("class") = CharacterVector::create("fb","list");
//   if(verbose) Rcout<<"done.\n";
//   return(l);  
// }
