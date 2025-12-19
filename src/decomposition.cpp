#include <Rcpp.h>
using namespace Rcpp;


//Meentemeyer (1978)
//' @keywords internal
// [[Rcpp::export("decomposition_annualLitterDecompositionRate")]]
double annualLitterDecompositionRate(double AET, double lignin) {
  double ki = (-0.5365+0.00241*AET) - (-0.01586+0.000056*AET)*lignin;
  return(ki);
}


//' @rdname decomposition
//' @param npool number of carbon pools
//' @param sand  percent sand
//' @param clay  percent clay
//' @param O2   effect of soil anaerobic conditions on decomposition (0-1)
//' @param strlig lignin fraction: (1) SRFC and (2) SOIL structural litter (g lignin/g biomass)
//' @param cwdlig lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
//' @param rsplig fraction of carbon lost as respiration (lignin)
//' @param Kmix base mixing rate: SOM2(SRFC) -> SOM2(SOIL), 1/sec
//' @param K_s21 rate constant: total loss from SOM2(SRFC), 1/sec
//' 
//'   
//' @keywords internal
void updateCarbonTransferMatrices(List comunicationDecomp, 
                                  double sand, double clay, double O2, 
                                  NumericVector strlig, NumericVector cwdlig, 
                                  double rsplig, double Kmix, double K_s21) {
  

  int npool = 12;
  NumericMatrix A = comunicationDecomp["A"];
  NumericMatrix respf = comunicationDecomp["respf"];
  NumericMatrix pathf = comunicationDecomp["pathf"];
  
  // anaerobic factor
  double fanerb = 1.0 + 5.0 * (1.0 - O2);
  
  // Updates pathf: fractional carbon flow from pool j to pool i
  pathf(7,0) = 1.0;
  pathf(8,1) = 1.0;
  pathf(7,2) = 1.0 - strlig[0];
  pathf(9,2) = strlig[0];
  pathf(8,3) = 1.0 - strlig[1];
  pathf(10,3) = strlig[1];
  
  pathf(7,4) = 1.0 - cwdlig[0];
  pathf(9,4) = cwdlig[0];
  pathf(7,5) = 1.0 - cwdlig[1];
  pathf(9,5) = cwdlig[1];
  pathf(8,6) = 1.0 - cwdlig[2];
  pathf(10,6) = cwdlig[2];
  
  pathf(9,7) = 1.0;
  pathf(11,8) = (0.003 + 0.032 * clay/100) * fanerb;
  pathf(10,8) = 1 - pathf(12,9);
  pathf(10,9) = Kmix / K_s21;
  pathf(7,9) = 1.0 - pathf(10,9);
  pathf(11,10) = (0.003 + 0.009 * clay/100) * fanerb;
  pathf(8,10) = 1.0 - pathf(11,10);
  pathf(8,11) = 1.0;
  
  // Updates respf: fractional respiration loss for carbon flow from pool j to pool i
  
  respf(7,0) = 0.55;
  respf(8,1) = 0.55;
  respf(7,2) = 0.45;
  respf(9,2) = rsplig;
  respf(8,3) = 0.55;
  respf(10,3) = rsplig;
  
  respf(7,4) = 0.45;
  respf(9,4) = rsplig;
  respf(7,5) = 0.45;
  respf(9,5) = rsplig;
  respf(8,6) = 0.55;
  respf(10,6) = rsplig;
  
  respf(9,7) = 0.60;
  respf(11,8) = 0.0;
  respf(10,8) = (0.17 + 0.68 * sand/100.0) / pathf(10,8);
  respf(10,0) = 0.0;
  respf(7,9) = 0.55;
  respf(11,10) = 0.0;
  respf(8,10) = 0.55 / pathf(8,10);
  respf(8,11) = 0.55 * O2;
  
  
  // Update carbon transfer matrix: fractional carbon flow from pool j that enters pool i
  for(int i=0;i<npool;i++) {
    A(i,i) = -1.0; 
    for(int j=0;j<npool;j++) {
      if(j!=i) {
        A(i,j) = pathf(i,j) * (1.0 - respf(i,j));
      }
    }
  }
}

// function for pH effect
double pH_factor(double x, double a, double b, double c, double d) {
  double pi = 3.141592;
  double pHeff = b + (c / pi) * atan(d * (x - a) * pi);
  pHeff = std::max(std::min(pHeff, 1.0), 0.0);
  return(pHeff);
}


//' @param ligninPercent lignin content (% of dry)
//' @param Nmass  nitrogen content (mg N / g dry)
// [[Rcpp::export("decomposition_litterMetabolicFraction")]]
double litterMetabolicFraction(double ligninPercent, double Nmass) {
  double fnit = Nmass/1000.0; //to g N/g dry
  double flig = ligninPercent/100.0; //to g lignin /g dry
  double rlig2n = flig/fnit;
  double fmet = 0.85 - 0.013 * rlig2n;
  return(fmet);
}

// Environmental scalar for each carbon pool adjusts base
// rate for soil temperature and soil moisture scalars
// (cdi) and additionally pH, lignin, texture, anaerobic,
// and cultivation
//'  @param cdi soil temperature and moisture scalar (0-1)
//'  @param pH  soil pH
//'  @param O2 effect of soil anaerobic conditions on decomposition (0-1)
//'  @param sand percent sand
//'  @param strlig lignin fraction: (1) surface and (2) soil structural litter (g lignin/g biomass)
//'  @param cwdlig lignin fraction: (1) fine branch; (2) large wood; (3) coarse root
//'  @param cultfac effect of cultivation on decomposition (1:SOM1, 2:SOM2, 3:SOM3, 4:structural)
//'  @param K base decomposition rate (1/sec)
//'  @param Kmix base mixing rate: SOM2(surface) -> SOM2(soil), 1/sec
//' 
//' % Output
//' %   K_s21       ! rate constant: total loss from SOM2(surface), 1/sec
//' %   xi          ! environmental scalar
void updateDecompositionRateScalars(List comunicationDecomp, 
                                    double cdi, double pH, double O2, double sand, 
                                    NumericVector strlig, NumericVector cwdlig, NumericVector cultfac, 
                                    NumericMatrix K, double Kmix) {

  NumericMatrix xi = comunicationDecomp["xi"];
  double pHeff, textureEff;
  
  // ==========================================
  //  metabolic litter (surface)
  //  ==========================================
    
  pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
  xi(0,0) = cdi * pHeff;

  // ==========================================
  //  metabolic litter (soil)
  // ==========================================
    
  pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
  xi(1,1) = cdi * pHeff * O2;

  // ==========================================
  // structural litter (surface)
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(2,2) = cdi * pHeff * exp(-3.0*strlig[0]);

  // ==========================================
  //  structural litter (soil)
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(3,3) = cdi * pHeff * exp(-3.0*strlig[1]) * O2 * cultfac[3];
  
  // ==========================================
  //  coarse woody debris: fine branch
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(4,4) = cdi * pHeff * exp(-3.0*cwdlig[0]);

  // ==========================================
  //  coarse woody debris: large wood
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(5,5) = cdi * pHeff * exp(-3.0*cwdlig[1]);

  // ==========================================
  //  coarse woody debris: coarse root
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(6,6) = cdi * pHeff * exp(-3.0*cwdlig[2]) * O2;
  
  // ==========================================
  //  active soil organic matter: SOM1 (surface)
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(7,7) = cdi * pHeff;

  // ==========================================
  //  active soil organic matter: SOM1 (SOIL)
  // ==========================================
  pHeff = pH_factor(pH, 4.8, 0.5, 1.14, 0.7);
  textureEff = 0.25 + 0.75 * (sand/100);
  xi(8,8) = cdi * pHeff * O2 * textureEff * cultfac[0];

  // ==========================================
  //  slow soil organic matter: SOM2 (surface)
  // ==========================================
  // som2(surface) -> som1(surface)
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  double K_s21_to_s11 = K(9,9) * pHeff;
  // som2(surface) -> som2(soil): mixing
  double K_s21_to_s22 = Kmix;
  // total loss from som2(surface)
  double K_s21 = K_s21_to_s11 + K_s21_to_s22;
  // effective environmental scalar
  xi(9,9) = cdi * (K_s21 / K(9,9));
  
  // ==========================================
  // slow soil organic matter: SOM2 (soil)
  // ==========================================
  pHeff = pH_factor(pH, 4.0, 0.5, 1.1, 0.7);
  xi(10,10) = cdi * pHeff * O2 * cultfac[1];
  
  // ==========================================
  //  passive soil organic matter: SOM3
  // ==========================================
  pHeff = pH_factor(pH, 3.0, 0.5, 1.1, 0.7);
  xi(11,11) = cdi * pHeff * O2 * cultfac[2];
  
}

// 
// 
// function [strlig, B] = litter_partition_matrix (leaf_flig, leaf_cn, froot_flig, froot_cn)
//   
//   % Matrix to partition litter fluxes to each carbon pool
// 
// % ---------------------------------------------------------------------------------------------
// % Input
// %   leaf_flig   ! leaf litter lignin fraction
// %   leaf_cn     ! leaf litter C:N (gC/gN)
// %   froot_flig  ! fine root litter lignin fraction
// %   froot_cn    ! fine root C:N (gC/gN)
// %
// % Output
// %   strlig      ! lignin fraction: (1) SRFC and (2) SOIL structural litter (g lignin/g biomass)
// %   B           ! litter flux partitioning matrix
// % ---------------------------------------------------------------------------------------------
// 
// % --- leaf
// 
// % fraction of plant residue that is nitrogen (g N / g biomass)
// 
// fnit = 1 / (leaf_cn * 2.5);
// 
// % lignin fraction of plant residue (g lignin / g biomass)
//   
//   flig = leaf_flig;
// 
// % lignin/nitrogen ratio of plant residue
// 
// rlig2n = flig / fnit;
// 
// % fraction of plant residue that goes to metabolic litter pool
// 
// fmet = 0.85 - 0.013 * rlig2n;
// 
// % make sure the fraction of residue which is lignin is not greater than the
// % fraction which goes to structural
// 
// if (flig > (1 - fmet))
//   fmet = (1 - flig);
// end
// 
// % minimum metabolic fraction
// 
// if (fmet < 0.20) 
//   fmet = 0.20;
// end
// 
// % adjust lignin content of structural litter pool. fligst is the fraction of
// % incoming structural residue that is lignin
// 
// fligst = min(flig/(1 - fmet), 1);
// 
// % DAYCENT adjusts lignin content of structural litter pool for new material
// % that is added. That is not needed in this program because only one type
// % of leaf litter is added. If different types of leaf litter are added,
// % need to adjust structural litter lignin
// 
// strlig(1) = fligst;
// 
// % B(i,j) = litter flux partitioning matrix (litter flux j -> pool i)
//   
//   B( 1,1) = fmet;
// B( 2,1) = 0;
// B( 3,1) = 1 - fmet;
// B( 4,1) = 0;
// B( 5,1) = 0;
// B( 6,1) = 0;
// B( 7,1) = 0;
// B( 8,1) = 0;
// B( 9,1) = 0;
// B(10,1) = 0;
// B(11,1) = 0;
// B(12,1) = 0;
// 
// % --- fine root
// 
// % fraction of plant residue that is nitrogen (g N / g biomass)
// 
// fnit = 1 / (froot_cn * 2.5);
// 
// % lignin fraction of plant residue (g lignin / g biomass)
//   
//   flig = froot_flig;
// 
// % lignin/nitrogen ratio of plant residue
// 
// rlig2n = flig / fnit;
// 
// % fraction of plant residue that goes to metabolic litter pool
// 
// fmet = 0.85 - 0.013 * rlig2n;
// 
// % make sure the fraction of residue which is lignin is not greater than the
// % fraction which goes to structural
// 
// if (flig > (1 - fmet))
//   fmet = (1 - flig);
// end
// 
// % minimum metabolic fraction
// 
// if (fmet < 0.20)
//   fmet = 0.20;
// end
// 
// % adjust lignin content of structural litter pool. fligst is the fraction of
// % incoming structural residue that is lignin
// 
// fligst = min(flig/(1 - fmet), 1);
// 
// % DAYCENT adjusts lignin content of structural litter pool for new material
// % that is added. That is not needed in this program because only one type
// % of fine root litter is added. If different types of fine root litter are
// % added, need to adjust structural litter lignin
// 
// strlig(2) = fligst;
// 
// % B(i,j) = litter flux partitioning matrix (litter flux j -> pool i)
//   
//   B( 1,2) = 0;
// B( 2,2) = fmet;
// B( 3,2) = 0;
// B( 4,2) = 1 - fmet;
// B( 5,2) = 0;
// B( 6,2) = 0;
// B( 7,2) = 0;
// B( 8,2) = 0;
// B( 9,2) = 0;
// B(10,2) = 0;
// B(11,2) = 0;
// B(12,2) = 0;
// 
// % --- cwd (fine branch)
//   
//   B( 1,3) = 0;
// B( 2,3) = 0;
// B( 3,3) = 0;
// B( 4,3) = 0;
// B( 5,3) = 1;
// B( 6,3) = 0;
// B( 7,3) = 0;
// B( 8,3) = 0;
// B( 9,3) = 0;
// B(10,3) = 0;
// B(11,3) = 0;
// B(12,3) = 0;
// 
// % cwd (large wood)
//   
//   B( 1,4) = 0;
// B( 2,4) = 0;
// B( 3,4) = 0;
// B( 4,4) = 0;
// B( 5,4) = 0;
// B( 6,4) = 1;
// B( 7,4) = 0;
// B( 8,4) = 0;
// B( 9,4) = 0;
// B(10,4) = 0;
// B(11,4) = 0;
// B(12,4) = 0;
// 
// % cwd (coarse root)
//   
//   B( 1,5) = 0;
// B( 2,5) = 0;
// B( 3,5) = 0;
// B( 4,5) = 0;
// B( 5,5) = 0;
// B( 6,5) = 0;
// B( 7,5) = 1;
// B( 8,5) = 0;
// B( 9,5) = 0;
// B(10,5) = 0;
// B(11,5) = 0;
// B(12,5) = 0;