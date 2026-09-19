#include <RcppArmadillo.h>
#include "phenology_c.h"
#include "hydraulics_c.h"
#include "modelInput_c.h"
#include "carbon_c.h"
#include "decomposition_c.h"

double leafDevelopmentStatus_c(double Sgdd, double gdd, double unfoldingDD) {
  double ds = 0.0;
  if(Sgdd>0.0) {
    if(gdd>Sgdd) ds = std::min(1.0, (gdd - Sgdd)/unfoldingDD);
  } else {
    ds = 1.0;
  }
  return(ds);
}
bool leafSenescenceStatus_c(double Ssen, double sen) {
  if(sen>Ssen) return(true);
  return false;
}

void updatePhenology_c(ModelInput& x, int doy, double photoperiod, double tmean) {

  double unfoldingDD = x.control.phenology.unfoldingDD;
  
  int numCohorts = x.cohorts.CohortCode.size();
  
  for(int j=0;j<numCohorts;j++) {
    //Copy previous expansion level
    x.internalPhenology.phiPrev[j] = x.internalPhenology.phi[j];
    
    if(x.paramsPhenology.phenoType[j] == "winter-deciduous" || x.paramsPhenology.phenoType[j] == "winter-semideciduous") {
      if(doy>200) {
        x.internalPhenology.gdd[j] = 0.0;
        if(photoperiod>x.paramsPhenology.Phsen[j]) { //Primary growth still possible until decrease of photoperiod
          x.internalPhenology.sen[j] = 0.0;
          if(x.paramsPhenology.growthDeterminacy[j] == "intermediate") x.internalPhenology.leafUnfolding[j] = true;
          x.internalPhenology.leafSenescence[j] = false;
          x.internalPhenology.budFormation[j] = false;
          x.internalPhenology.leafDormancy[j] = false;
          x.internalPhenology.leafOrganogenesisDuration[j] = 0;
        } else { //Start (or continue) accumulation
          x.internalPhenology.leafUnfolding[j] = false; //Primary growth has arrested
          if(!x.internalPhenology.leafDormancy[j]) { // Check temperature accumulation until dormancy occurs
            double rsen = 0.0;
            if(tmean-x.paramsPhenology.Tbsen[j]<0.0) {
              rsen = pow(x.paramsPhenology.Tbsen[j]-tmean, x.paramsPhenology.xsen[j])*pow(photoperiod/x.paramsPhenology.Phsen[j], x.paramsPhenology.ysen[j]);
            }
            x.internalPhenology.sen[j] = x.internalPhenology.sen[j] + rsen;
            x.internalPhenology.leafSenescence[j] = leafSenescenceStatus_c(x.paramsPhenology.Ssen[j],x.internalPhenology.sen[j]);
            x.internalPhenology.leafDormancy[j] = x.internalPhenology.leafSenescence[j];
          }
          if(x.internalPhenology.leafDormancy[j]) x.internalPhenology.phi[j] = 0.0;
          x.internalPhenology.budFormation[j] = false;
          if(x.paramsPhenology.growthDeterminacy[j] == "intermediate") {
            x.internalPhenology.budFormation[j] = true;
            if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
              x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
            } else {
              //Stops bud formation after organogenesis
              x.internalPhenology.budFormation[j] = false;
            }
          } 
        }
        // Rcout << doy<< " "<< photoperiod<<" "<< rsen <<" "<< sen[j]<<" "<<  leafSenescence[j] << "\n";
      } else if (doy<=200) { //Only increase in the first part of the year and if doy > t0gdd
        x.internalPhenology.sen[j] = 0.0;
        x.internalPhenology.budFormation[j] = false;
        x.internalPhenology.leafSenescence[j] = false;
        if((tmean-x.paramsPhenology.Tbgdd[j]>0.0) && (doy >= ((int) x.paramsPhenology.t0gdd[j]))) x.internalPhenology.gdd[j] = x.internalPhenology.gdd[j] + (tmean - x.paramsPhenology.Tbgdd[j]);
        x.internalPhenology.phi[j] = leafDevelopmentStatus_c(x.paramsPhenology.Sgdd[j], x.internalPhenology.gdd[j], unfoldingDD);
        x.internalPhenology.leafUnfolding[j] = (x.internalPhenology.phi[j]>0.0);
        x.internalPhenology.leafDormancy[j] = (x.internalPhenology.phi[j]==0.0);
        if(x.paramsPhenology.growthDeterminacy[j] == "determinate") {
          if(x.internalPhenology.phi[j]==1.0) {
            x.internalPhenology.budFormation[j] = true;
            if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
              x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
            } else {
              //Stops bud formation after organogenesis
              x.internalPhenology.budFormation[j] = false;
            }
          }
        }
        // Rcout << doy<< " "<< photoperiod<<" "<< gdd[j]<<" "<<  leafUnfolding[j] << "\n";
      }
    }
    else if(x.paramsPhenology.phenoType[j] == "evergreen") {
      if(x.paramsPhenology.growthDeterminacy[j] == "indeterminate") {
        x.internalPhenology.leafSenescence[j] = true;
        x.internalPhenology.leafUnfolding[j] = true;
        x.internalPhenology.budFormation[j] = false;
        x.internalPhenology.leafDormancy[j] = false;
        x.internalPhenology.leafOrganogenesisDuration[j] = medfate::NA_INTEGER;
      } else {
        if(doy>200) {
          x.internalPhenology.leafUnfolding[j] = false; //Primary growth has arrested
          x.internalPhenology.leafSenescence[j] = false;
          x.internalPhenology.gdd[j] = 0.0;
          x.internalPhenology.phi[j] = 0.0;
          if(photoperiod>x.paramsPhenology.Phsen[j]) {
            x.internalPhenology.sen[j] = 0.0;
            if(x.paramsPhenology.growthDeterminacy[j] == "intermediate") x.internalPhenology.leafUnfolding[j] = true;
            x.internalPhenology.budFormation[j] = false;
            x.internalPhenology.leafDormancy[j] = false;
            x.internalPhenology.leafOrganogenesisDuration[j] = 0;
          } else {
            if(!x.internalPhenology.leafDormancy[j]) { // Check temperature accumulation until dormancy occurs
              double rsen = 0.0;
              if(tmean-x.paramsPhenology.Tbsen[j]<0.0) {
                rsen = pow(x.paramsPhenology.Tbsen[j]-tmean, x.paramsPhenology.xsen[j])*pow(photoperiod/x.paramsPhenology.Phsen[j], x.paramsPhenology.ysen[j]);
              }
              x.internalPhenology.sen[j] = x.internalPhenology.sen[j] + rsen;
              x.internalPhenology.leafDormancy[j] = leafSenescenceStatus_c(x.paramsPhenology.Ssen[j],x.internalPhenology.sen[j]);
            }
            x.internalPhenology.budFormation[j] = false;
            if(x.paramsPhenology.growthDeterminacy[j] == "intermediate") {
              x.internalPhenology.budFormation[j] = true;
              if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
                x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
              } else {
                //Stops bud formation after organogenesis
                x.internalPhenology.budFormation[j] = false;
              }
            } 
          }
        } else if (doy<=200) { //Only increase in the first part of the year
          x.internalPhenology.sen[j] = 0.0;
          x.internalPhenology.budFormation[j] = false;
          if(x.internalPhenology.phi[j] < 1.0) {
            if((tmean-x.paramsPhenology.Tbgdd[j]>0.0) && (doy >= ((int) x.paramsPhenology.t0gdd[j]))) x.internalPhenology.gdd[j] = x.internalPhenology.gdd[j] + (tmean - x.paramsPhenology.Tbgdd[j]);
            x.internalPhenology.phi[j] = leafDevelopmentStatus_c(x.paramsPhenology.Sgdd[j], x.internalPhenology.gdd[j],unfoldingDD);
            x.internalPhenology.leafSenescence[j] = (x.internalPhenology.phi[j]>0.0 && x.internalPhenology.phi[j]<1.0);
            x.internalPhenology.leafUnfolding[j] = (x.internalPhenology.phi[j]>0.0);
            x.internalPhenology.leafDormancy[j] = (x.internalPhenology.phi[j]==0.0);
          } else {
            x.internalPhenology.leafUnfolding[j] = false;
            x.internalPhenology.leafSenescence[j] = false;
            x.internalPhenology.leafDormancy[j] = false;
            if(x.paramsPhenology.growthDeterminacy[j] == "determinate") {
              x.internalPhenology.budFormation[j] = true;
              if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
                x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
              } else {
                //Stops bud formation after organogenesis
                x.internalPhenology.budFormation[j] = false;
              }
            }
          }
          // Rcout<<j<< " phi: "<< ph<<"\n";
        }
      } 
    }    
  }
  // Rcout<<"\n";
}

void updateLeaves_c(ModelInput& x, double wind, bool fromGrowthModel) {
  
  
  int numCohorts = x.cohorts.CohortCode.size();
  
  for(int j=0;j<numCohorts;j++) {
    //Leaf fall
    bool leafFall = true;
    if(x.paramsPhenology.phenoType[j] == "winter-semideciduous") leafFall = x.internalPhenology.leafUnfolding[j];
    if(leafFall) {
      double LAIlitter = x.above.LAI_dead[j]*(1.0 - exp(-1.0*(wind/10.0)));//Decrease dead leaf area according to wind speed
      x.above.LAI_dead[j] = x.above.LAI_dead[j] - LAIlitter;
      if(fromGrowthModel) {
          // from m2/m2 to g C/m2
          double leaf_litter = leafCperDry*1000.0*LAIlitter/x.paramsAnatomy.SLA[j];
          double twig_litter = leaf_litter/(x.paramsAnatomy.r635[j] - 1.0);
          addLeafTwigLitter_c(x.cohorts.SpeciesName[j], leaf_litter, twig_litter,
                              x.internalLitter, x.paramsLitterDecomposition,
                              x.internalSOC);
      }
    } 
    //Leaf unfolding, senescence and defoliation only dealt with if called from spwb
    if(!fromGrowthModel) {
      double psiLeafPLC = xylemPsi_c(1.0 - x.internalWater.LeafPLC[j],1.0, x.paramsTranspiration.VCleafapo_c[j], x.paramsTranspiration.VCleafapo_d[j]);
      if(x.paramsPhenology.phenoType[j] == "winter-deciduous" || x.paramsPhenology.phenoType[j] == "winter-semideciduous") {
        if((x.internalPhenology.leafSenescence[j]) && (x.above.LAI_expanded[j]>0.0)) {
          double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
          x.above.LAI_expanded[j] = 0.0; //Update expanded leaf area (will decrease if LAI_live decreases)
          x.above.LAI_dead[j] += LAI_exp_prev;//Check increase dead leaf area if expanded leaf area has decreased
          x.internalPhenology.leafSenescence[j] = false;
        } else {
          if(x.control.defoliation.cavitationInducedDefoliation) {
            double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
            x.above.LAI_expanded[j] = x.above.LAI_live[j]*std::min(x.internalPhenology.phi[j], 1.0 - proportionDefoliationWeibull_c(psiLeafPLC, x.paramsTranspiration.VCleafapo_c[j], x.paramsTranspiration.VCleafapo_d[j], x.control.defoliation.criticalLeafPLC, x.control.defoliation.cvLeafP50)); //Update expanded leaf area (will decrease if LAI_live decreases)
            x.above.LAI_dead[j] += std::max(0.0, LAI_exp_prev - x.above.LAI_expanded[j]); // Add senescence leaves to dead
          } else {
            x.above.LAI_expanded[j] = x.above.LAI_live[j]*x.internalPhenology.phi[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
          }
        }
      } else {
        //Apply defoliation effects to evergreens
        if(x.control.defoliation.cavitationInducedDefoliation) {
          double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
          x.above.LAI_expanded[j] = x.above.LAI_live[j]*(1.0 - proportionDefoliationWeibull_c(psiLeafPLC, x.paramsTranspiration.VCleafapo_c[j], x.paramsTranspiration.VCleafapo_d[j], x.control.defoliation.criticalLeafPLC, x.control.defoliation.cvLeafP50));
          x.above.LAI_dead[j] += std::max(0.0, LAI_exp_prev - x.above.LAI_expanded[j]); // Add senescence leaves to dead
        } else {
          x.above.LAI_expanded[j] = x.above.LAI_live[j];
        }
      } 
    }
  }    
}