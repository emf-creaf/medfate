#include "inner_sperry_c.h"

void initSperryNetwork_inner_c(SperryNetwork& network,
                               int c,
                               const InternalWater& internalWater, 
                               const TranspirationParams& paramsTranspiration, 
                               const WaterStorageParams& paramsWaterStorage,
                               const std::vector<double>& VCroot_kmax, 
                               const std::vector<double>& VGrhizo_kmax,
                               const std::vector<double>& psiSoil, 
                               const std::vector<double>& VG_n, 
                               const std::vector<double>& VG_alpha,
                               const ControlParameters& control,
                               double sapFluidityDay) {
  
   int nlayers = VG_n.size();
  
  network.sperryParams = control.sperry;
  network.psisoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.krhizomax = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.nsoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.alphasoil = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  network.krootmax = std::vector<double>(nlayers, medfate::NA_DOUBLE);
  for(int l=0;l<nlayers;l++) {
    network.psisoil[l] = psiSoil[l];
    network.krhizomax[l] = VGrhizo_kmax[l];
    network.krootmax[l] = sapFluidityDay*VCroot_kmax[l];
    network.alphasoil[l] = VG_alpha[l];
    network.nsoil[l] = VG_n[l];
  }
  network.kstemmax = sapFluidityDay*paramsTranspiration.VCstem_kmax[c];
  network.kleafmax = sapFluidityDay*paramsTranspiration.VCleaf_kmax[c];
  network.kleafapomax = sapFluidityDay*paramsTranspiration.VCleafapo_kmax[c];
  network.kleafsymp = sapFluidityDay*paramsTranspiration.kleaf_symp[c];
  network.stemc = paramsTranspiration.VCstem_c[c];
  network.stemd = paramsTranspiration.VCstem_d[c];
  network.leafc = paramsTranspiration.VCleaf_c[c];
  network.leafd = paramsTranspiration.VCleaf_d[c];
  network.stemd = paramsTranspiration.VCstem_d[c];
  if(control.sperry.leafCavitationEffects) {
    network.PLCleaf = internalWater.LeafPLC[c];
  } else {
    network.PLCleaf = 0.0;
  }
  network.PLCstem = internalWater.StemPLC[c];
}
