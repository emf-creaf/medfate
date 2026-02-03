#include <RcppArmadillo.h>
#include "windKatul_c.h"
#include <math.h> 

Rcpp::DataFrame copyCanopyTurbulenceModelResult_c(const CanopyTurbulenceModel_RESULT& canopyTurbulenceModel) {
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("z1") = Rcpp::wrap(canopyTurbulenceModel.z1),
                                                   Rcpp::Named("U1") = Rcpp::wrap(canopyTurbulenceModel.U1),
                                                   Rcpp::Named("dU1") = Rcpp::wrap(canopyTurbulenceModel.dU1),
                                                   Rcpp::Named("epsilon1") = Rcpp::wrap(canopyTurbulenceModel.epsilon1),
                                                   Rcpp::Named("k1") = Rcpp::wrap(canopyTurbulenceModel.k1),
                                                   Rcpp::Named("uw1") = Rcpp::wrap(canopyTurbulenceModel.uw1),
                                                   Rcpp::Named("Lmix1") = Rcpp::wrap(canopyTurbulenceModel.Lmix1));
  return(output);
}

Rcpp::DataFrame copyCanopyTurbulenceResult_c(const CanopyTurbulence_RESULT& canopyTurbulence) {
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("zmid") = Rcpp::wrap(canopyTurbulence.zmid),
                                                   Rcpp::Named("u") = Rcpp::wrap(canopyTurbulence.u),
                                                   Rcpp::Named("du") = Rcpp::wrap(canopyTurbulence.du),
                                                   Rcpp::Named("epsilon") = Rcpp::wrap(canopyTurbulence.epsilon),
                                                   Rcpp::Named("k") = Rcpp::wrap(canopyTurbulence.k),
                                                   Rcpp::Named("uw") = Rcpp::wrap(canopyTurbulence.uw));
  return(output);
}
  
  
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
const double alphaCNT=0.05;

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
void thomas_c(std::vector<double>& q,
              const std::vector<double>& aa, 
              const std::vector<double>& bb, 
              const std::vector<double>& cc, 
              const std::vector<double>& dd) {
  int n=bb.size();
  std::vector<double> bet(n);
  std::vector<double> gam(n);
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
}

// Generates a vector with linear spacing
std::vector<double> linspace_c(double x1, double x2, int N) {
  std::vector<double> v(N);
  v[0] = x1;
  double s = (x2 - x1)/((double)(N - 1));
  for(int i=1;i<N;i++) {
    v[i] = v[i-1] + s;
  }
  return(v);
}
// Implementation of R which
std::vector<int> which_c(std::vector<bool>& l) {
  int c = 0;
  for(int i=0;i<l.size();i++) if(l[i]) c++;
  std::vector<int> w(c);
  int cnt=0;
  for(int i=0;i<l.size();i++) if(l[i]) {w[cnt] = i;cnt++;}
  return(w);
}

void windCanopyTurbulenceModel_inner_c(CanopyTurbulenceModel_RESULT& comm, 
                                       const std::vector<double>& zm, 
                                       const std::vector<double>& Cx, 
                                       double hm, double d0, double z0,
                                       std::string model = "k-epsilon") {
  if(zm.size() != Cx.size()) throw medfate::MedfateInternalError("Height and effective drag vectors should have the same length!");
  int N=zm.size();
  // Find zmax
  double zmax=0.0;
  for(int i=0;i<N;i++) zmax = std::max(zmax, zm[i]);
  double dz=zm[1]-zm[0];
  // ------- Define starting conditions for U/u*, k/(u*^2), epsilon/(u*3/h)
  double Ulow=0.0;
  double Uhigh=(1.0/kv)*log((zmax-d0)/z0);
  std::vector<double> U=linspace_c(Ulow,Uhigh,N);
  
  double khigh=0.5*((AAu*AAu) + (AAv*AAv) + (AAw*AAw));
  double klow=0.001*khigh;
  std::vector<double> k=linspace_c(klow,khigh,N);
  
  double epsilonhigh=(1.0/(kv*(zmax-d0)));
  double epsilonlow=0.001*epsilonhigh;
  std::vector<double> epsilon=linspace_c(epsilonlow,epsilonhigh,N);
  
  // Mixing Length model
  double alpha = kv*(hm-d0)/hm;  //fraction of mixing length (Lmixing = alpha x h) inside the canopy up to canopy top
  std::vector<double> Lmix(N, alpha*hm);
  std::vector<bool> bw(N);
  for(int i=0;i<N;i++) bw[i] = zm[i] < hm;
  std::vector<int> c1 = which_c(bw);
  int nn = 0;
  for(int i=0;i<c1.size();i++) nn = std::max(nn, c1[i]);
  for(int i=nn;i<N;i++) {
    Lmix[i] = kv*(zm[i] - d0);
  }
  Lmix[nn]=(Lmix[nn-1]+Lmix[nn+1])/2.0;
  
  double am3=-(Aq*Aq*Aq)*((AAu*AAu)-(AAw*AAw))/((AAw*AAw)-(Aq*Aq)/3.0);

  std::vector<double> y(N); 
  std::vector<double> vt(N); //viscosity
  std::vector<double> dvt(N); //derivative
  std::vector<double> dU(N); //derivative
  std::vector<double> upd(N), dia(N), lod(N);
  std::vector<double> aa(N), bb(N), cc(N), dd(N);
  std::vector<double> Su(N),Sk(N),Se1(N),Se2(N); 
  std::vector<double> a1(N), a2(N), a3(N);
  std::vector<double> uw(N);
  double eps1=0.1;
  double maxerr=9999.9;
  
  std::vector<double> Un(N);
  std::vector<double> Kn(N);
  std::vector<double> epsilonn(N);
  
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
        epsilon[i] = (3.0/3.0)*(pow(sqrt(2.0*k[i]),3.0)/am3*Lmix[i]);
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
    thomas_c(Un, aa,bb,cc,dd);
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
    thomas_c(Kn, aa,bb,cc,dd);
    maxerr=0.0;
    for(int i=0;i<N;i++) maxerr = std::max(maxerr, abs(Kn[i]-k[i]));
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
      thomas_c(epsilonn, aa,bb,cc,dd);
      //   Use successive relaxations in iterations
      for(int i=0;i<N;i++) {
        epsilon[i]=std::abs(eps1*epsilonn[i]+(1.0-eps1)*epsilon[i]);
      }
    }
    cnt++;
    if(cnt==maxcnt) Rcpp::warning("too many iterations in canopy turbulence model");
  }
  for(int i=0;i<N;i++) {
    comm.z1[i] = zm[i];
    comm.U1[i] = U[i];
    comm.dU1[i] = dU[i];
    comm.epsilon1[i] = epsilon[i];
    comm.k1[i] = k[i];
    comm.uw1[i] = uw[i]/(abs(uw[N-1]));
    comm.Lmix1[i] = Lmix[i];
  }
}

void windCanopyTurbulence_inner_c(CanopyTurbulence_RESULT& output, CanopyTurbulenceModel_RESULT& comm,
                                  const std::vector<double>& zmid, 
                                  const std::vector<double>& LAD, 
                                  double canopyHeight,
                                  double u, double windMeasurementHeight = 200, 
                                  std::string model = "k-epsilon") {
  
  //z - height vector in m
  std::vector<double> zm(zmid.size());
  for(int i=0;i<zmid.size();i++) zm[i]= zmid[i]/100.0;
  //Effective drag = Cd x leaf area density
  std::vector<double> Cx(LAD.size());
  for(int i=0;i<LAD.size();i++) Cx[i]= LAD[i]*0.2;
  //hm - canopy height (m)
  double hm = (canopyHeight/100.0);
  //d0 - displacement height (m)
  double d0 = 0.67*hm;
  //z0 - Momentum roughness height (m)
  double z0 = 0.08*hm;
  
  //u_f - Friction velocity
  double u_f = u*kv/log(((hm + windMeasurementHeight/100.0)-d0)/z0);
  
  windCanopyTurbulenceModel_inner_c(comm, zm, Cx, hm, d0, z0);
  
  for(int i=0;i<zmid.size();i++) {
    output.zmid[i] = zmid[i];
    output.u[i] = comm.U1[i]*u_f;
    output.du[i] = comm.dU1[i]*u_f;
    output.epsilon[i] = comm.epsilon1[i]*((u_f*u_f*u_f)/hm);
    output.k[i] = comm.k1[i]*u_f*u_f;
    output.uw[i] = comm.uw1[i]*u_f*u_f;
  }
}
