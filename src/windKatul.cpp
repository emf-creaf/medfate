#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

/* K-Epsilon Models of Katul et al. (2004)
 * Code adapted from Matlab to Rcpp
 * www: https://nicholas.duke.edu/people/faculty/katul/k_epsilon_model.htm
*/

/* (A)
 *   
 * Second order turbulence model closure constants for canopy-layer simulation.
 * Revised by csz: Fri Mar 28 2003
 */

//Rate of MKE lost by drag converted into TKE (in [0, 1])
// (roughly-) Fitted value for canopy-layer velocity profile simulation
const double Bp=1.0;
// Mixing lenght constant (Seginer, 1974; Massman and Weil, 1999)
// alpha=0.09 (Katul and Chang, 1999)
double alphaCNT=0.05;

// ASL Values for sigma_u/u*, sigma_v/u*, sigma_w/u*
// (set upper boundary condition on TKE)
const double AAu=2.3;
const double AAv=2.1;
const double AAw=1.25;
const double Aq=0.5*((AAu*AAu) + (AAv*AAv) + (AAw*AAw));
      
// Determine the Kolmogorov constant from Au, Av, Aw
// (constant for the turbulent viscocity)
const double Cu=1.0/(Aq*Aq);
      
// Von karman constant
const double kv=0.4;
      
// Constants for the dissipation budget (Launder and Spalding, 1974)
const double Ce1=1.44;
const double Ce2=1.92;
      
// Prandtl numbers
// Prandtl number (PrTKE = 1) for TKE hardwired in solver routines.
const double PrTKE=1.0;
// dissipation budget (Detering and Etling, 1974)
const double Pr=(kv*kv)/(sqrt(Cu)*(Ce2-Ce1));
      
// Wake TKE budget coefficients (Sanz, 2003)
const double Cg=pow(2.0/alphaCNT,2.0/3.0);
      
// Shortcircuit in Dissipation (linear in Bp)
const double Bd=sqrt(Cu)*Cg*Bp + 3.0/PrTKE;
      
// Dissipation budget (constant with alpha)
const double Ce4=PrTKE*(2.0/Pr-sqrt(Cu)/6.0*Cg*(Ce2-Ce1));                
      
// Ce5[=1.5]=Ce4 (Green, 1992), not [=0.6]!=Ce4 (Liu et al., 1996).
// *** Proof: see ``A note on k-epsilon modelling...'' (Sanz, 2003)
// *** The '1.5' value above for engineering (e.g. wind tunnel) flows (Cu=0.09)
const double Ce5=Ce4;

/* (B)
 * 
 * tri-diagonal solver needed for the implicit schemes used in solving the non-linear ODEs.
 */
NumericVector thomas(NumericVector aa, NumericVector bb, NumericVector cc, NumericVector dd) {
  int n=bb.length();
  NumericVector bet(n);
  NumericVector gam(n);
  NumericVector q(n);
  bet[0]=bb[0];
  gam[0]=dd[0]/bb[0];
  for(int i=1;i<n;i++){
    bet[i]=bb[i]-(aa[i]*cc[i-1]/bet[i-1]);
    gam[i]=(dd[i]-aa[i]*gam[i-1])/bet[i];
  }
  q[n-1]=gam[n-1];
  for(int i=(n-2); i >=0;i--) {
    q[i]=gam[i]-(cc[i]*q[i+1]/bet[i]);  
  }
  return(q);
}

// Generates a vector with linear spacing
NumericVector linspace(double x1, double x2, int N) {
  NumericVector v(N);
  v[0] = x1;
  double s = (x2 - x1)/((double)(N - 1));
  for(int i=1;i<N;i++) {
    v[i] = v[i-1] + s;
  }
  return(v);
}
// Implementation of R which
IntegerVector which(LogicalVector l) {
  int c = 0;
  for(int i=0;i<l.size();i++) if(l[i]) c++;
  IntegerVector w(c);
  int cnt=0;
  for(int i=0;i<l.size();i++) if(l[i]) {w[cnt] = i;cnt++;}
  return(w);
}

NumericMatrix windCanopyTurbulenceModel_inner(NumericVector zm, NumericVector Cx, double hm, double d0, double z0,
                                              String model = "k-epsilon") {
  if(zm.size() != Cx.size()) stop("Height and effective drag vectors should have the same length!");
  int N=zm.size();
  double zmax=max(zm);
  double dz=zm[1]-zm[0];
  // ------- Define starting conditions for U/u*, k/(u*^2), epsilon/(u*3/h)
  double Ulow=0.0;
  double Uhigh=(1.0/kv)*log((zmax-d0)/z0);
  NumericVector U=linspace(Ulow,Uhigh,N);
  
  double khigh=0.5*((AAu*AAu) + (AAv*AAv) + (AAw*AAw));
  double klow=0.001*khigh;
  NumericVector k=linspace(klow,khigh,N);
  
  double epsilonhigh=(1.0/(kv*(zmax-d0)));
  double epsilonlow=0.001*epsilonhigh;
  NumericVector epsilon=linspace(epsilonlow,epsilonhigh,N);
  
  // Mixing Length model
  double alpha = kv*(hm-d0)/hm;  //fraction of mixing length (Lmixing = alpha x h) inside the canopy up to canopy top
  NumericVector Lmix(N, alpha*hm);
  IntegerVector c1 = which(zm < hm);
  int nn=max(c1);
  for(int i=nn;i<N;i++) {
    Lmix[i] = kv*(zm[i] - d0);
  }
  Lmix[nn]=(Lmix[nn-1]+Lmix[nn+1])/2.0;
  
  double am3=-(Aq*Aq*Aq)*((AAu*AAu)-(AAw*AAw))/((AAw*AAw)-(Aq*Aq)/3.0);
  NumericVector lambda3=am3*Lmix;
  
  NumericVector y(N); 
  NumericVector vt(N); //viscosity
  NumericVector dvt(N); //derivative
  NumericVector dU(N); //derivative
  NumericVector upd(N), dia(N), lod(N);
  NumericVector aa(N), bb(N), cc(N), dd(N);
  NumericVector Su(N),Sk(N),Se1(N),Se2(N); 
  NumericVector a1(N), a2(N), a3(N);
  NumericVector uw(N);
  double eps1=0.1;
  double maxerr=9999.9;
  
  int cnt=0;
  int maxcnt = 100;
  while((maxerr>0.1) && (cnt < maxcnt)) {
    // Viscocity (and derivative) Model
    for(int i=0;i<N;i++) {
      vt[i]=pow(Cu,1.0/4.0)*Lmix[i]*sqrt(std::abs(k[i])); 
    }
    for(int i=1;i<N;i++) {
      y[i] = (vt[i] - vt[i-1])/dz;
    }
    y[0] = y[1];
    for(int i=1;i<N;i++) {
      dvt[i] = y[i];
    }
    //  dU/dz
    for(int i=1;i<N;i++) {
      y[i] = (U[i] - U[i-1])/dz;
    }
    y[0] = y[1];
    for(int i=0;i<N;i++) {
      dU[i] = y[i];
    }
    //  Compute the Reynolds stress
    for(int i=0;i<N;i++) {
      uw[i]= (-1.0)*vt[i]*dU[i];
    }
    uw[N-1]=uw[N-2];
    uw[0]=uw[1];
    if(model=="k-U") { // Estimate dissipation rate in k-U model
      for(int i=0;i<N;i++) {
        epsilon[i] = (3.0/3.0)*(pow(sqrt(2.0*k[i]),3.0)/lambda3[i]);
      }
    }
    //   Set up coefficients for Mean Momentum ODE-------------------------------------------
    for(int i=0;i<N;i++) {
      a1[i]=vt[i];
      a2[i]=dvt[i];
      a3[i]=(-1.0)*Cx[i]*std::abs(U[i]);
    }
    double dx=dz;
    //  ------ Set the elements of the Tri-diagonal Matrix
    for(int i=0;i<N;i++) {
      upd[i]=(a1[i]/(dx*dx)+a2[i]/(2.0*dx));
      dia[i]=(a1[i]*(-2.0)/(dx*dx)+a3[i]);
      lod[i]=(a1[i]/(dx*dx)-a2[i]/(2.0*dx));
      aa[i]=lod[i];
      bb[i]=dia[i];
      cc[i]=upd[i];
      dd[i] = 0.0;
    }
    aa[0]=0.0;
    bb[0]=1.0;
    cc[0]=-1.0;
    dd[0]= 0.0;
    aa[N-1]=0.0;
    bb[N-1]=1.0;
    cc[N-1]=0.0;
    dd[N-1]=Uhigh;
    //   Use the Thomas Algorithm to solve the tridiagonal matrix
    NumericVector Un=thomas(aa,bb,cc,dd);
    //   Use successive relaxations in iterations
    for(int i=0;i<N;i++) {
      U[i]=std::abs(eps1*Un[i]+(1.0-eps1)*U[i]);
    }
    //   Set up coefficients for TKE ODE------------------------------------------------------------
    for(int i=0;i<N;i++) {
      a1[i]=vt[i];
      a2[i]=dvt[i];
      a3[i]=(-1.0)*Bd*std::abs(U[i])*Cx[i];
    }
    dx=dz;
    //  ------ Set the elements of the Tri-diagonal Matrix
    for(int i=0;i<N;i++) {
      upd[i]=(a1[i]/(dx*dx)+a2[i]/(2.0*dx));
      dia[i]=(a1[i]*(-2.0)/(dx*dx)+a3[i]);
      lod[i]=(a1[i]/(dx*dx)-a2[i]/(2.0*dx));
      aa[i]=lod[i];
      bb[i]=dia[i];
      cc[i]=upd[i];
      dd[i] = epsilon[i]-Cx[i]*Bp*pow(std::abs(U[i]),3.0)-vt[i]*(dU[i]*dU[i]);
    }
    aa[0]=0.0;
    bb[0]=1.0;
    cc[0]=-1.0;
    dd[0]= 0.0;
    aa[N-1]=0.0;
    bb[N-1]=1.0;
    cc[N-1]=0.0;
    dd[N-1]=khigh;
    
    //   ------Use the Thomas Algorithm to solve the tridiagonal matrix
    NumericVector Kn=thomas(aa,bb,cc,dd);
    maxerr=max(abs(Kn-k));
    // -----Use successive relaxations in iterations
    for(int i=0;i<N;i++){
      k[i]=std::abs(eps1*Kn[i]+(1.0-eps1)*k[i]);
    }
    
    if(model=="k-epsilon") {
      //   Set up coefficients for dissipation ODE ---------------------------------------------------------------------------------
      for(int i=0;i<N;i++){
        Su[i]=Cx[i]*(U[i]*U[i]); 
        Sk[i]=Cx[i]*(Bp*(U[i]*U[i]*U[i])-Bd*U[i]*k[i]);
        Se1[i]=(Ce4*epsilon[i]/k[i])*Sk[i];
        Se2[i]=epsilon[i]*Cx[i]*((Ce4*Bp*(U[i]*U[i]*U[i])/k[i])-Bd*Ce5*U[i]);
      }
      for(int i=0;i<N;i++){
        a1[i]=vt[i]/Pr;
        a2[i]=dvt[i]/Pr;
        a3[i]=((-1.0)*Ce2*epsilon[i]/k[i]);
      }
      dx=dz;
      //  ------ Set the elements of the Tri-diagonal Matrix
      for(int i=0;i<N;i++) {
        upd[i]=(a1[i]/(dx*dx)+a2[i]/(2.0*dx));
        dia[i]=(a1[i]*(-2.0)/(dx*dx)+a3[i]);
        lod[i]=(a1[i]/(dx*dx)-a2[i]/(2.0*dx));
        aa[i]=lod[i];
        bb[i]=dia[i];
        cc[i]=upd[i];
        dd[i] = -Ce1*Cu*k[i]*(dU[i]*dU[i])-Se2[i];
      }
      aa[0]=0.0;
      bb[0]=1.0;
      cc[0]=-1.0;
      dd[0]= 0.0;
      aa[N-1]=0.0;
      bb[N-1]=1.0;
      cc[N-1]=0.0;
      dd[N-1]=epsilonhigh;
      //   Use the Thomas Algorithm to solve the tridiagonal matrix
      NumericVector epsilonn=thomas(aa,bb,cc,dd);
      //   Use successive relaxations in iterations
      for(int i=0;i<N;i++) {
        epsilon[i]=std::abs(eps1*epsilonn[i]+(1.0-eps1)*epsilon[i]);
      }
    }
    cnt++;
    if(cnt==maxcnt) warning("too many iterations in canopy turbulence model");
  }
  NumericMatrix res(zm.size(), 7);
  res(_,0) = zm;
  res(_,1) = U;
  res(_,2) = dU;
  res(_,3) = epsilon;
  res(_,4) = k;
  res(_,5) = uw/(abs(uw[N-1]));
  res(_,6) = Lmix;  
  res.attr("dimnames") = List::create(seq(1,zm.size()), 
           CharacterVector::create("z1", "U1" ,"dU1", "epsilon1", "k1", "uw1", "Lmix1"));
  return(res);
}

void windCanopyTurbulence_inner(DataFrame output, NumericVector zmid, NumericVector LAD, double canopyHeight,
                                double u, double windMeasurementHeight = 200, String model = "k-epsilon") {
  
  //z - height vector in m
  NumericVector zm(zmid.size());
  for(int i=0;i<zmid.size();i++) zm[i]= zmid[i]/100.0;
  //Effective drag = Cd x leaf area density
  NumericVector Cx(LAD.size());
  for(int i=0;i<LAD.size();i++) Cx[i]= LAD[i]*0.2;
  //hm - canopy height (m)
  double hm = (canopyHeight/100.0);
  //d0 - displacement height (m)
  double d0 = 0.67*hm;
  //z0 - Momentum roughness height (m)
  double z0 = 0.08*hm;
  
  //u_f - Friction velocity
  double u_f = u*kv/log(((hm + windMeasurementHeight/100.0)-d0)/z0);
  NumericMatrix cmout = windCanopyTurbulenceModel_inner(zm, Cx, hm, d0, z0);
  
  //"z1", "U1" ,"dU1", "epsilon1", "k1", "uw1", "Lmix1"
  NumericVector U1 = cmout(_,1); 
  NumericVector dU1 = cmout(_,2);
  NumericVector epsilon1 = cmout(_,3);
  NumericVector k1 = cmout(_,4);
  NumericVector uw1 = cmout(_,5);
  NumericVector out_zmid = output["zmid"];
  NumericVector out_u = output["u"];
  NumericVector out_du = output["du"];
  NumericVector out_epsilon = output["epsilon"];
  NumericVector out_k = output["k"];
  NumericVector out_uw = output["uw"];
  for(int i=0;i<zmid.size();i++) {
    out_zmid[i] = zmid[i];
    out_u[i] = U1[i]*u_f;
    out_du[i] = dU1[i]*u_f;
    out_epsilon[i] = epsilon1[i]*((u_f*u_f*u_f)/hm);
    out_k[i] = k1[i]*u_f*u_f;
    out_uw[i] = uw1[i]*u_f*u_f;
  }
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
//'   \item{\code{uw1}: <uw>/(u*^2), where <uw> is the Reynolds stress and u* is friction velocity.}
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
  NumericMatrix res = windCanopyTurbulenceModel_inner(zm, Cx, hm, d0, z0, model);
  return(as<DataFrame>(res));
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
  DataFrame output = DataFrame::create(Named("zmid") = NumericVector(n, NA_REAL),
                    Named("u") = NumericVector(n, NA_REAL),
                    Named("du") = NumericVector(n, NA_REAL),
                    Named("epsilon") = NumericVector(n, NA_REAL),
                    Named("k") = NumericVector(n, NA_REAL),
                    Named("uw") = NumericVector(n, NA_REAL));
  windCanopyTurbulence_inner(output, zmid, LAD, canopyHeight,
                             u, windMeasurementHeight, model);
  return(output);
}
