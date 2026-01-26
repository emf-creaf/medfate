#include "medfate.h"
#include "control.h"
#include "soil_c.h"
#include "Rcpp.h"

#ifndef MODELINPUT_C_H
#define MODELINPUT_C_H

struct CanopyParams {
  std::vector<double> zlow;
  std::vector<double> zmid;
  std::vector<double> zup;
  std::vector<double> LAIlive;
  std::vector<double> LAIexpanded;
  std::vector<double> LAIdead;
  std::vector<double> Tair;
  std::vector<double> Cair;
  std::vector<double> VPair;
};
struct PhenologyParams {
  std::vector<std::string> phenoType;
  std::vector<double> leafDuration;
  std::vector<double> t0gdd;
  std::vector<double> Sgdd;
  std::vector<double> Tbgdd;
  std::vector<double> Ssen;
  std::vector<double> Phsen;
  std::vector<double> Tbsen;
  std::vector<double> xsen;
  std::vector<double> ysen;
  
};

struct InterceptionParams {
  std::vector<double> LeafAngle;
  std::vector<double> LeafAngleSD;
  std::vector<double> Beta_p;
  std::vector<double> Beta_q;
  std::vector<double> ClumpingIndex;
  std::vector<double> kPAR;
  std::vector<double> kSWR;
  std::vector<double> alphaSWR;
  std::vector<double> alphaLWR;
  std::vector<double> g;
};
struct AnatomyParams {
  std::vector<double> Hmax;
  std::vector<double> Hmed;
  std::vector<double> Al2As;
  std::vector<double> Ar2Al;
  std::vector<double> SLA;
  std::vector<double> LeafWidth;
  std::vector<double> LeafDensity;
  std::vector<double> WoodDensity;
  std::vector<double> FineRootDensity;
  std::vector<double> conduit2sapwood;
  std::vector<double> SRL;
  std::vector<double> RLD;
  std::vector<double> r635;
};
struct WaterStorageParams {
  std::vector<double> maxFMC;
  std::vector<double> maxMCleaf;
  std::vector<double> maxMCstem;
  std::vector<double> LeafPI0;
  std::vector<double> LeafEPS;
  std::vector<double> LeafAF;
  std::vector<double> Vleaf;
  std::vector<double> StemPI0;
  std::vector<double> StemEPS;
  std::vector<double> StemAF;
  std::vector<double> Vsapwood;  
};

/**
 * Transpiration parameters structure, including Granier, Sperry and Sureau parameters
 */
struct TranspirationParams {
  std::vector<double> Tmax_LAI;
  std::vector<double> Tmax_LAIsq;
  std::vector<double> Psi_Extract;
  std::vector<double> Exp_Extract;
  std::vector<double> Gswmin;
  std::vector<double> Gswmax;
  std::vector<double> Gs_Toptim;
  std::vector<double> Gs_Tsens;
  std::vector<double> Gsw_AC_slope;
  std::vector<double> Gs_P50;
  std::vector<double> Gs_slope;
  std::vector<double> WUE;
  std::vector<double> WUE_par;
  std::vector<double> WUE_co2;
  std::vector<double> WUE_vpd;
  std::vector<double> Vmax298;
  std::vector<double> Jmax298;
  std::vector<double> Kmax_stemxylem;
  std::vector<double> Kmax_rootxylem;
  std::vector<double> VCleaf_kmax;
  std::vector<double> VCleafapo_kmax;
  std::vector<double> VCleaf_slope;
  std::vector<double> VCleaf_P50;
  std::vector<double> VCleaf_c;
  std::vector<double> VCleaf_d;
  std::vector<double> kleaf_symp;
  std::vector<double> VCstem_kmax;
  std::vector<double> VCstem_slope;
  std::vector<double> VCstem_P50;
  std::vector<double> VCstem_c;
  std::vector<double> VCstem_d;
  std::vector<double> kstem_xylem;
  std::vector<double> kstem_symp;
  std::vector<double> VCroottot_kmax;
  std::vector<double> VCroot_slope;
  std::vector<double> VCroot_P50;
  std::vector<double> VCroot_c;
  std::vector<double> VCroot_d;
  std::vector<double> VGrhizotot_kmax;
  std::vector<double> Plant_kmax;
  std::vector<double> FR_leaf;
  std::vector<double> FR_stem;
  std::vector<double> FR_root;
  
};

struct InternalPhenology {
  std::vector<double> gdd;
  std::vector<double> sen;
  std::vector<bool> budFormation;
  std::vector<bool> leafUnfolding;
  std::vector<bool> leafSenescence;
  std::vector<bool> leafDormancy;
  std::vector<double> phi;
};

struct InternalWater {
  std::vector<double> Einst;
  std::vector<double> Elim;
  std::vector<double> Emin_L;
  std::vector<double> Emin_S;
  std::vector<double> PlantPsi;
  std::vector<double> RootCrownPsi;
  std::vector<double> LeafPsi;
  std::vector<double> StemPsi;
  std::vector<double> LeafSympPsi;
  std::vector<double> StemSympPsi;
  std::vector<double> LeafPLC;
  std::vector<double> StemPLC;
};

struct InternalCarbon {
  std::vector<double> sugarLeaf;
  std::vector<double> starchLeaf;
  std::vector<double> sugarSapwood;
  std::vector<double> starchSapwood;
};

struct InternalMortality {
  std::vector<double> N_dead;
  std::vector<double> N_starvation;
  std::vector<double> N_dessication;
  std::vector<double> N_burnt;
  std::vector<double> N_resprouting_stumps;
  std::vector<double> Cover_dead;
  std::vector<double> Cover_starvation;
  std::vector<double> Cover_dessication;
  std::vector<double> Cover_burnt;
  std::vector<double> Cover_resprouting_stumps;
  std::vector<double> Snag_smallbranches;
  std::vector<double> Snag_largewood;
};



struct InternalAllocation {
  std::vector<double> allocationTarget;
  std::vector<double> leafAreaTarget;
  std::vector<double> sapwoodAreaTarget;
  std::vector<double> fineRootBiomassTarget;
  std::vector<double> crownBudPercent;
};

struct InternalSnags {
  std::vector<std::string> Species;
  std::vector<double> DBH;
  std::vector<double> Height;
  std::vector<double> DeathAge;
  std::vector<double> SmallBranches;
  std::vector<double> LargeWood;
};

struct InternalLitter {
  std::vector<std::string> Species;
  std::vector<double> Leaves;
  std::vector<double> Twigs;
  std::vector<double> SmallBranches;
  std::vector<double> LargeWood;
  std::vector<double> CoarseRoots;
  std::vector<double> FineRoots;
};

struct InternalSOC {
  double SurfaceMetabolic;
  double SoilMetabolic;
  double SurfaceActive;
  double SoilActive;
  double SurfaceSlow;
  double SoilSlow;
  double SoilPassive;
};

struct InternalLAIDistribution{
  std::vector<double> PrevLAIexpanded;
  std::vector<double> PrevLAIdead;
  std::vector<double> PARcohort;
  medfate::Matrix<double> live;
  medfate::Matrix<double> expanded;
  medfate::Matrix<double> dead;
};

struct InternalFCCS {
  std::vector<double> w;
  std::vector<double> cover;
  std::vector<double> hbc;
  std::vector<double> htc;
  std::vector<double> habc;
  std::vector<double> hatc;
  std::vector<double> delta;
  std::vector<double> rhob;
  std::vector<double> rhop;
  std::vector<double> PV;
  std::vector<double> beta;
  std::vector<double> betarel;
  std::vector<double> etabetarel;
  std::vector<double> sigma;
  std::vector<double> pDead;
  std::vector<double> FAI;
  std::vector<double> h;
  std::vector<double> RV;
  std::vector<double> MinFMC;
  std::vector<double> MaxFMC;
  std::vector<double> ActFMC;
};

/**
 * Internal C++ representation of the model input data structure
 */
class ModelInput {
  public:
    ControlParameters control;
    Soil soil;
    double snowpack;
    CanopyParams canopy;
    double herbLAI;
    double herbLAImax;
    PhenologyParams paramsPhenology;
    InterceptionParams paramsInterception;
    AnatomyParams paramsAnatomy;
    WaterStorageParams paramsWaterStorage;
    TranspirationParams paramsTranspiration;
    
    InternalPhenology internalPhenology;
    InternalWater internalWater;
    InternalCarbon internalCarbon;
    InternalMortality internalMortality;
    InternalAllocation internalAllocation;
    InternalSnags internalSnags;
    InternalLitter internalLitter;
    InternalSOC internalSOC;
    InternalLAIDistribution internalLAIDistribution;
    InternalFCCS internalFCCS;
    
    std::string version;
    
    ModelInput(Rcpp::List x);
};

#endif