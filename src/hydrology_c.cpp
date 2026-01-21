#include <cmath>
#include <vector>
#include <meteoland.h>

/*=============================================================================
 * Implementation of hydrology routines using C++ code
 *=============================================================================*/

//' @rdname hydrology_soilEvaporation
//' 
//' @param DEF Water deficit in the (topsoil) layer.
//' @param PETs Potential evapotranspiration at the soil surface.
//' @param Gsoil Gamma parameter (maximum daily evaporation).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_soilEvaporationAmount")]]
double soilEvaporationAmount_c(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t + 1.0)-sqrt(t)), PETs);
  return(Esoil);
}


//Old defaults
//ERconv=0.05, ERsyn = 0.2
//New defaults
//Rconv = 5.6, Rsyn = 1.5

// Rainfall intensity calculation
double rainfallIntensity_c(int month, double prec, std::vector<double> rainfallIntensityPerMonth){
  double Ri_month = rainfallIntensityPerMonth[month - 1];
  double Ri = std::max(prec/24.0,Ri_month);
  return(Ri); 
}


// [[Rcpp::export(".hydrology_interceptionGashDay")]]
double interceptionGashDay_c(double Rainfall, double Cm, double p, double ER=0.05) {
  double I = 0.0;
  double PG = (-Cm/(ER*(1.0-p)))*log(1.0-ER); //Rainfall need to saturate the canopy
  if(Cm==0.0 || p==1.0) PG = 0.0; //Avoid NAs
  if(Rainfall>PG) {
    I = (1-p)*PG + (1-p)*ER*(Rainfall-PG);
  } else {
    I = (1-p)*Rainfall;
  }
  return(I);
}

// [[Rcpp::export(".hydrology_interceptionLiuDay")]]
double interceptionLiuDay_c(double Rainfall, double Cm, double p, double ER=0.05){
  double I = Cm*(1.0 - exp(-1.0*(Rainfall)*((1.0 - p)/Cm)))*(1.0 - (ER/(1.0 - p))) + (ER*Rainfall);
  return(I);
}

void infiltrationRepartition_c(int nlayers, double I, 
                               double* Ivec, double* widths, double* macro, 
                               double a = -0.005, double b = 3.0) {
  double pi = 0.0;
  double z1 = 0.0;
  double p1 = 1.0;
  for(int i=0;i<nlayers;i++) {
    double ai = a*pow(1.0-macro[i],b);
    if(i<(nlayers-1)) {
      pi = p1*(1.0-exp(ai*widths[i]));
    } else {
      pi = p1;
    }
    p1 = p1*exp(ai*widths[i]);
    z1 = z1 + widths[i];
    Ivec[i] = I*pi;
  }
}


//' @rdname hydrology_infiltration
//' 
//' @param input A numeric vector of (daily) water input (in mm of water).
//' @param Ssoil Soil water storage capacity (can be referred to topsoil) (in mm of water).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationBoughton")]]
double infiltrationBoughton_c(double input, double Ssoil) {
  double I = 0;
  if(input>0.2*Ssoil) {
    I = input-(pow(input-0.2*Ssoil,2.0)/(input+0.8*Ssoil));
  } else {
    I = input;
  }
  return(I);
}

double fGreenAmpt_c(double x, double t, double psi_w, double Ksat, double delta_theta) {
  double f = Ksat*t + std::abs(psi_w)*delta_theta*log(1.0 + (x/(std::abs(psi_w)*delta_theta))) - x;
  return(f);
}
double fGreenAmptDer_c(double x, double t, double psi_w, double Ksat, double delta_theta) {
  double fder = (log(1.0 + (x/(std::abs(psi_w)*delta_theta)))/(1.0+ (x/(std::abs(psi_w)*delta_theta))))-1.0;
  return(fder);
}

//' @rdname hydrology_infiltration
//' 
//' @param t Time of the infiltration event
//' @param psi_w Matric potential at the wetting front
//' @param Ksat hydraulic conductivity at saturation
//' @param theta_sat volumetric content at saturation
//' @param theta_dry volumetric content at the dry side of the wetting front
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_infiltrationGreenAmpt")]]
double infitrationGreenAmpt_c(double t, double psi_w, double Ksat, double theta_sat, double theta_dry) {
   double delta_theta = theta_sat - theta_dry;
   double x,x1,e,fx,fx1;
   x1 = 0.0;//initial guess
   e = 0.001; // accuracy in mm
   int cnt = 0;
   int mxiter = 100;
   do {
     x=x1; /*make x equal to the last calculated value of  x1*/
     fx=fGreenAmpt_c(x, t, psi_w, Ksat, delta_theta);            //simplifying f(x)to fx
     fx1=fGreenAmptDer_c(x, t, psi_w, Ksat, delta_theta);            //simplifying fprime(x) to fx1
     x1=x-(fx/fx1);/*calculate x{1} from x, fx and fx1*/ 
     cnt++;
   } while ((std::abs(x1-x)>=e) && (cnt < mxiter));
   return(x);
 }


//' @rdname hydrology_verticalInputs
//' 
//' @param tday Average day temperature (ºC).
//' @param rad Solar radiation (in MJ/m2/day).
//' @param elevation Altitude above sea level (m).
//' 
//' @keywords internal
// [[Rcpp::export("hydrology_snowMelt")]]
double snowMelt_c(double tday, double rad, double LgroundSWR, double elevation) {
  //missing data checks
  if(std::isnan(rad)) Rcpp::stop("Missing radiation data for snow melt!\n");
  if(std::isnan(elevation)) Rcpp::stop("Missing elevation data for snow melt!\n");
  double rho = meteoland::utils_airDensity(tday, meteoland::utils_atmosphericPressure(elevation));
  double ten = (86400.0*tday*rho*1013.86*1e-6/100.0); //ten can be negative if temperature is below zero
  double ren = (rad*(LgroundSWR/100.0))*(0.1); //90% albedo of snow
  double melt = std::max(0.0,(ren+ten)/0.33355); //Do not allow negative melting values
  return(melt);
}

//Imbitition from macropores to micropores following Larsbo et al. 2005, eq. 6-7
//Diffusion equation with gradients in water content as driving force
//Returns m3/m3/s
double microporeImbibitionRate_c(double theta_b, double theta_micro, 
                                 double D_theta_b, double D_theta_micro,
                                 double S_macro) {
  double G_f = 3.0;//geometry factor
  double gamma_w = 0.1; //Scaling factor
  double D_w = ((D_theta_b + D_theta_micro)/2.0)*S_macro;//Effective water diffusivity
  double d = 9.35*1e-3; //Effective diffusion pathlength in m
  double S_w = std::max(0.0, ((G_f*D_w*gamma_w)/(d*d))*(theta_b - theta_micro));
  // Rcout<< "D_w " << D_w << " S_w "<< S_w<<"\n";
  return(S_w);
}


double rootFindingMacropores_c(double S_t, double K_up, double Ksat_ms, double Ksat_b_ms, double kin_exp,
                               double e_macro, double lambda, double dZ_m, double sourceSink_macro_m3s, double tstep, 
                               int Nmax = 100) {
  double a = 0.0;
  //function at 'a'
  double Ka = (Ksat_ms - Ksat_b_ms)*pow(a, kin_exp);
  double f_a = S_t - a + (tstep/(e_macro*lambda*dZ_m))*((K_up - Ka) + sourceSink_macro_m3s);
  double b = a + 1.0;
  //function at 'b'
  double Kb = (Ksat_ms - Ksat_b_ms)*pow(b, kin_exp);
  double f_b = S_t - b + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kb) + sourceSink_macro_m3s);
  while(f_b > 0.0) {
    b = b + 1.0;
    Kb = (Ksat_ms - Ksat_b_ms)*pow(b, kin_exp);
    f_b = S_t - b + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kb) + sourceSink_macro_m3s);
    if(b>10.0) Rcpp::stop("Could not find appropriate bounds for macropore circulation");
  }
  // Rcout<< "Ka "<< Ka <<" f_a "<< f_a<<" Kb "<< Kb <<" f_b "<< f_b<<"\n";
  
  int N = 1;
  double c = NA_REAL;
  double Kc = NA_REAL;
  double f_c = NA_REAL;
  double tol = 1e-7;
  bool has_to_stop = false;
  bool found = false;
  while(!has_to_stop) {
    c = (a + b)/2.0; // new midpoint
    //function at 'c'
    Kc = (Ksat_ms - Ksat_b_ms)*pow(c, kin_exp);
    f_c = S_t - c + (tstep/(e_macro*lambda*dZ_m))*((K_up - Kc) + sourceSink_macro_m3s);
    // Rcout<<N << " std::abs((b - a)/2.0 "<< std::abs((b - a)/2.0) << " c "<< c <<" f_c "<< f_c<<"\n";
    if((f_c == 0) || (std::abs((b - a)/2.0) < tol)) { // solution found
      found = true;
      has_to_stop = true;
    }
    if(((f_c > 0) && (f_a > 0)) || ((f_c < 0) && (f_a < 0)) ) { //new interval
      a = c;
      f_a = f_c;
    } else {
      b = c;
      f_b = f_c;
    }
    N = N + 1; //Increment step counter
    if(N == Nmax) {
      has_to_stop = true;
    }
  }
  if(!found) {
    Rcpp::stop("Not found");
  }
  return(c);
}


 