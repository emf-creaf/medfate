#include <RcppArmadillo.h>
#include "windKatul_c.h"
#include "windKatul.h"
#include <math.h> 
using namespace Rcpp;


void windCanopyTurbulence_inner(DataFrame output, NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u, double windMeasurementHeight = 200, String model = "k-epsilon") {
  int n = zmid.size();
  CanopyTurbulence_RESULT res(n);
  CanopyTurbulenceModel_RESULT comm(n);
  std::vector<double> zmid_vec = as<std::vector<double>>(zmid);
  std::vector<double> LAD_vec = as<std::vector<double>>(LAD);
  windCanopyTurbulence_inner_c(res, comm, 
                               zmid_vec, 
                               LAD_vec, 
                               canopyHeight,
                               u, windMeasurementHeight, 
                               model.get_cstring());
  output["zmid"] = Rcpp::wrap(res.zmid);
  output["u"] = Rcpp::wrap(res.u);
  output["du"] = Rcpp::wrap(res.du);
  output["epsilon"] = Rcpp::wrap(res.epsilon);
  output["k"] = Rcpp::wrap(res.k);
  output["uw"] = Rcpp::wrap(res.uw);
}

/* (C)
 * 
 *  K-epsilon model (equations 1-7 with equation 4b)
 *  K-U model (equations 9 and 10)
 *  
 *  k_epsilon_CSL(z, Cx, h, do,zo)
 *  output [z1, U1, k1, uw1, Lmix1]
 */
//' Models for canopy turbulence
//' 
//' Models for canopy turbulence by Katul et al (2004).
//' 
//' @param zm A numeric vector with height values (m).
//' @param Cx Effective drag = Cd x leaf area density.
//' @param hm Canopy height (m).
//' @param d0 Zero displacement height (m).
//' @param z0 Momentum roughness height (m).
//' @param zmid A numeric vector of mid-point heights (in cm) for canopy layers.
//' @param LAD A numeric vector of leaf area density values (m3/m2).
//' @param canopyHeight Canopy height (in cm).
//' @param u Measured wind speed (m/s).
//' @param windMeasurementHeight Height of wind speed measurement with respect to canopy height (cm).
//' @param model Closure model.
//' 
//' @return 
//' Function \code{wind_canopyTurbulenceModel} returns a data frame of vertical profiles for variables:
//' \itemize{
//'   \item{\code{z1}: Height values.}
//'   \item{\code{U1}: U/u*, where U is mean velocity and u* is friction velocity.}
//'   \item{\code{dU1}: dUdz/u*, where dUdz is mean velocity gradient and u* is friction velocity.}
//'   \item{\code{epsilon1}: epsilon/(u^3/h) where epsilon is the turbulent kinetic dissipation rate, u* is friction velocity and h is canopy height.}
//'   \item{\code{k1}: k/(u*^2), where k is the turbulent kinetic energy and u* is friction velocity.}
//'   \item{\code{uw1}: uw/(u*^2), where uw is the Reynolds stress and u* is friction velocity.}
//'   \item{\code{Lmix1}: Mixing length.}
//' }
//' 
//' Function \code{wind_canopyTurbulence} returns a data frame of vertical profiles for transformed variables:
//'   \itemize{
//'     \item{\code{zmid}: Input mid-point heights (in cm) for canopy layers.}
//'     \item{\code{u}: Wind speed (m/s).}
//'     \item{\code{du}: Mean velocity gradient (1/s).}
//'     \item{\code{epsilon}: Turbulent kinetic dissipation rate.}
//'     \item{\code{k}: Turbulent kinetic energy.}
//'     \item{\code{uw}: Reynolds stress.}
//'   }
//' 
//' @details
//' Implementation in Rcpp of the K-epsilon canopy turbulence models by Katul et al (2004) originally in Matlab code (https://nicholas.duke.edu/people/faculty/katul/k_epsilon_model.htm).
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @references 
//' Katul GG, Mahrt L, Poggi D, Sanz C (2004) One- and two-equation models for canopy turbulence. Boundary-Layer Meteorol 113:81–109. https://doi.org/10.1023/B:BOUN.0000037333.48760.e5
//' 
//' @seealso
//' \code{\link{vprofile_windExtinction}}
//' 
//' @examples
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Load example plot plant data
//' data(exampleforest)
//' 
//' #Canopy height (in m)
//' h= max(exampleforest$treeData$Height/100) 
//' d0 = 0.67*h
//' z0 = 0.08*h
//' 
//' #Height values (cm)
//' z = seq(50,1000, by=50)
//' zm = z/100 # (in m)
//' 
//' # Leaf area density
//' lad = vprofile_leafAreaDensity(exampleforest, SpParamsMED, draw = FALSE,
//'                                z = c(0,z))
//'   
//' # Effective drag
//' Cd = 0.2
//' Cx = Cd*lad
//'   
//' # canopy turbulence model
//' wind_canopyTurbulenceModel(zm, Cx,h,d0,z0)
//' 
//' @name wind
//' @keywords internal
// [[Rcpp::export("wind_canopyTurbulenceModel")]]
DataFrame windCanopyTurbulenceModel(NumericVector zm, NumericVector Cx, double hm, double d0, double z0,
                                        String model = "k-epsilon") {
  CanopyTurbulenceModel_RESULT comm(zm.size());
  std::vector<double> zm_vec = as<std::vector<double>>(zm);
  std::vector<double> Cx_vec = as<std::vector<double>>(Cx);
  windCanopyTurbulenceModel_inner_c(comm, 
                                    zm_vec, 
                                    Cx_vec, 
                                    hm, d0, z0,
                                    model.get_cstring());
  return(copyCanopyTurbulenceModelResult_c(comm));
}
/*
 *   zmid - Vector of mid heights for canopy layers (cm)
 *   LAD - Vector of leaf area density for canopy layers (m2/m3)
 *   canopyHeight - Canopy height (cm)
 *   u - Wind speed over the canopy (m/s)
 *   windMeasurementHeight - Height of wind measurement over canopy
 */
//' @rdname wind
//' @keywords internal
// [[Rcpp::export("wind_canopyTurbulence")]]
DataFrame windCanopyTurbulence(NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u, double windMeasurementHeight = 200, String model = "k-epsilon") {
  
  int n = zmid.size();
  CanopyTurbulence_RESULT res(n);
  CanopyTurbulenceModel_RESULT comm(n);
  std::vector<double> zmid_vec = as<std::vector<double>>(zmid);
  std::vector<double> LAD_vec = as<std::vector<double>>(LAD);
  windCanopyTurbulence_inner_c(res, comm, 
                               zmid_vec, 
                               LAD_vec, 
                               canopyHeight,
                               u, windMeasurementHeight, 
                               model.get_cstring());
  return(copyCanopyTurbulenceResult_c(res));
}
