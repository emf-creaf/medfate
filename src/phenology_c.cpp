#include <RcppArmadillo.h>
#include "phenology_c.h"
#include "hydraulics_c.h"
#include "modelInput_c.h"
#include "carbon_c.h"
#include "decomposition_c.h"

double leafDevelopmentStatus_c(double gdd, double Sgdd, double Ugdd) {
  double ds = 0.0;
  if(Sgdd>0.0) {
    if(gdd>Sgdd) ds = std::min(1.0, (gdd - Sgdd)/Ugdd);
  } else {
    ds = 1.0;
  }
  return(ds);
}


void updatePhenology_c(ModelInput& x, int doy, double photoperiod, double tmean) {


  int numCohorts = x.cohorts.CohortCode.size();
  
  for(int j=0;j<numCohorts;j++) {
    //Copy previous expansion level
    x.internalPhenology.phiPrev[j] = x.internalPhenology.phi[j];
    //Indetermined evergreen
    if(x.paramsPhenology.phenoType[j] == "evergreen" && x.paramsPhenology.growthDeterminacy[j] == "indeterminate") {
      x.internalPhenology.leafSenescence[j] = true;
      x.internalPhenology.leafUnfolding[j] = true;
      x.internalPhenology.budFormation[j] = false;
      x.internalPhenology.budDormancy[j] = false;
      x.internalPhenology.phi[j] = 1.0;
      x.internalPhenology.leafOrganogenesisDuration[j] = medfate::NA_INTEGER;
    } else { //All the rest
      if(doy>212) { // Second part of the year
        x.internalPhenology.gdd[j] = 0.0;
        x.internalPhenology.leafUnfolding[j] = false; //Primary growth has in principle arrested (but see below)
        x.internalPhenology.leafSenescence[j] = false; //Senescence has not occurred (but see below)
        if(photoperiod>x.paramsPhenology.Phsen[j]) { // Before photoperiod has decreased below threshold
          x.internalPhenology.sen[j] = 0.0; //Set to zero senescence cumulative variable
          if(x.paramsPhenology.growthDeterminacy[j] == "intermediate") { // Allow neo-growth for intermediate species
            x.internalPhenology.leafUnfolding[j] = true;
            x.internalPhenology.leafOrganogenesisDuration[j]=0;
            x.internalPhenology.budDormancy[j] = false; //No bud dormancy during neo-growth
          }
        } else { // After photoperiod has decreased below threshold
          //Address senescence model
          if(x.internalPhenology.sen[j] < (x.paramsPhenology.Ssen[j] + x.paramsPhenology.Fsen[j])) { 
            //Start bud formation for intermediate species if it has not yet occurred
            if(x.paramsPhenology.growthDeterminacy[j] == "intermediate" && !x.internalPhenology.budFormation[j] && x.internalPhenology.leafOrganogenesisDuration[j]==0) {
              x.internalPhenology.budFormation[j] = true;
            } 
            // Check temperature accumulation until dormancy occurs
            double rsen = 0.0;
            if(tmean-x.paramsPhenology.Tbsen[j]<0.0) {
              rsen = pow(x.paramsPhenology.Tbsen[j]-tmean, x.paramsPhenology.xsen[j])*pow(photoperiod/x.paramsPhenology.Phsen[j], x.paramsPhenology.ysen[j]);
            }
            x.internalPhenology.sen[j] = x.internalPhenology.sen[j] + rsen;
            bool ecoDormancy = (x.internalPhenology.sen[j] >= x.paramsPhenology.Ssen[j]) && (x.internalPhenology.sen[j] < (x.paramsPhenology.Ssen[j] + x.paramsPhenology.Fsen[j])); //Eco-dormancy
            bool completeAbscission = (x.internalPhenology.sen[j] >= (x.paramsPhenology.Ssen[j] + x.paramsPhenology.Fsen[j])); //Abscission
            if(ecoDormancy) {
              //If para-dormancy was no active then activate eco-dormancy
              if(!x.internalPhenology.budDormancy[j]) x.internalPhenology.budDormancy[j] = true;
              //Sets bud formation to false when dormancy starts
              x.internalPhenology.budFormation[j] = false;
              //Trigger senescence for winter (semi) deciduous and evergreens with autumn senescence
              if(x.paramsPhenology.senescencePeriod[j]=="autumn") {
                x.internalPhenology.phi[j] = 1.0 - (x.internalPhenology.sen[j] - x.paramsPhenology.Ssen[j])/x.paramsPhenology.Fsen[j];
                // Rcpp::Rcout<< "j" <<j << "phi " << x.internalPhenology.phi[j] <<"\n";
                x.internalPhenology.leafSenescence[j] = true;
              }
            } else if(completeAbscission) {
              if(x.paramsPhenology.senescencePeriod[j]=="autumn") {
                x.internalPhenology.phi[j] = 0.0;
                x.internalPhenology.leafSenescence[j] = false;
              }
            }
          } 
        }
        // Advance bud formation if necessary
        if(x.internalPhenology.budFormation[j]) {
          if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
            x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
          } else {
            //Stops bud formation after organogenesis
            x.internalPhenology.budFormation[j] = false;
            x.internalPhenology.budDormancy[j] = true; //Sets para-dormancy
          }
        }
      } else if (doy<=212) {
        x.internalPhenology.sen[j] = 0.0;
        x.internalPhenology.budFormation[j] = false;
        if(doy==212) x.internalPhenology.phi[j] = 1.0; //Force remaining flush
        if(x.internalPhenology.phi[j] < 1.0) {
          //Set GDD if DOY is large enough and temperature is large enough
          if(doy >= ((int) x.paramsPhenology.t0gdd[j])) x.internalPhenology.gdd[j] = x.internalPhenology.gdd[j] + std::max(0.0, tmean - x.paramsPhenology.Tbgdd[j]);
          //Update phi
          x.internalPhenology.phi[j] = leafDevelopmentStatus_c(x.internalPhenology.gdd[j], x.paramsPhenology.Sgdd[j], x.paramsPhenology.Ugdd[j]);
          // Rcpp::Rcout << " DOY: "<< doy << " GDD: " << x.internalPhenology.gdd[j] << " PHI: " <<x.internalPhenology.phi[j] <<"\n";
          //Force senescence for evergreen (determinate or intermediate) species during leaf elongation 
          if((x.paramsPhenology.phenoType[j] == "evergreen") && (x.paramsPhenology.senescencePeriod[j]=="spring")) {
            x.internalPhenology.leafSenescence[j] = (x.internalPhenology.phi[j]>0.0 && x.internalPhenology.phi[j]<1.0);
          }
          x.internalPhenology.leafUnfolding[j] = (x.internalPhenology.phi[j]>0.0);
          x.internalPhenology.budDormancy[j] = (x.internalPhenology.phi[j]==0.0);
          x.internalPhenology.leafOrganogenesisDuration[j] = 0;
        } else { // After pre-formed leaf elongation
          x.internalPhenology.leafSenescence[j] = false;
          //Determinate species conduct bud formation
          if(x.paramsPhenology.growthDeterminacy[j] == "determinate") {
            x.internalPhenology.leafUnfolding[j] = false;
            x.internalPhenology.budFormation[j] = true;
            if(x.internalPhenology.leafOrganogenesisDuration[j] < x.paramsPhenology.budFormationDays[j]) {
              x.internalPhenology.leafOrganogenesisDuration[j] = x.internalPhenology.leafOrganogenesisDuration[j] + 1;
            } else {
              //Stops bud formation after organogenesis, if finishes before DOY 200
              x.internalPhenology.budFormation[j] = false;
              x.internalPhenology.budDormancy[j] = true; //Sets para-dormancy
            }
          } else { // Plants with intermediate strategy can continue growth (bud formation delayed to autumn)
            x.internalPhenology.leafUnfolding[j] = true;
            x.internalPhenology.budFormation[j] = false;
            x.internalPhenology.budDormancy[j] = false;
          }
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
          double LAI_exp_new = x.above.LAI_live[j]*x.internalPhenology.phi[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
          if(LAI_exp_new < x.above.LAI_expanded[j]) { // Previous drought defoliation may have reduced LAI
            x.above.LAI_expanded[j] = LAI_exp_new;
            x.above.LAI_dead[j] += (LAI_exp_prev - x.above.LAI_expanded[j]); 
          }
        } else {
          if(x.control.defoliation.cavitationInducedDefoliation) {
            double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
            double LAI_exp_new = x.above.LAI_live[j]*std::min(x.internalPhenology.phi[j], 1.0 - 
                                                              proportionDefoliationWeibull_c(psiLeafPLC, x.paramsTranspiration.VCleafapo_c[j], 
                                                                                             x.paramsTranspiration.VCleafapo_d[j], x.control.defoliation.criticalLeafPLC, x.control.defoliation.cvLeafP50));
            x.above.LAI_expanded[j] = LAI_exp_new;
            if(LAI_exp_new < x.above.LAI_expanded[j]) {
              x.above.LAI_dead[j] += std::max(0.0, LAI_exp_prev - x.above.LAI_expanded[j]); // Add senescence leaves to dead
            }
          } else {
            x.above.LAI_expanded[j] = x.above.LAI_live[j]*x.internalPhenology.phi[j]; //Update expanded leaf area (will decrease if LAI_live decreases)
          }
        }
      } else if((x.paramsPhenology.phenoType[j] == "evergreen") && (x.paramsPhenology.senescencePeriod[j]=="autumn") && (x.paramsPhenology.growthPeriod[j]=="spring")) {
        if(x.internalPhenology.leafSenescence[j] && (x.above.LAI_expanded[j]>0.0)) {
          double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
          double LAI_exp_new = x.above.LAI_live[j]*((1.0 - (1.0/x.paramsPhenology.leafDuration[j])) + 
                                                    (x.internalPhenology.phi[j]/x.paramsPhenology.leafDuration[j])); //Update expanded leaf area (will decrease if LAI_live decreases)
          if(LAI_exp_new < x.above.LAI_expanded[j]) { // Previous drought defoliation may have reduced LAI
            x.above.LAI_expanded[j] = LAI_exp_new;
            x.above.LAI_dead[j] += (LAI_exp_prev - x.above.LAI_expanded[j]); 
          }
        } else {
          if(x.control.defoliation.cavitationInducedDefoliation) {
            double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
            double LAI_exp_new = x.above.LAI_live[j]*std::min(1.0 - (1.0/x.paramsPhenology.leafDuration[j]) + (x.internalPhenology.phi[j]/x.paramsPhenology.leafDuration[j]), 
                                                              1.0 - proportionDefoliationWeibull_c(psiLeafPLC, x.paramsTranspiration.VCleafapo_c[j], 
                                                                                             x.paramsTranspiration.VCleafapo_d[j], x.control.defoliation.criticalLeafPLC, x.control.defoliation.cvLeafP50));
            x.above.LAI_expanded[j] = LAI_exp_new;
            if(LAI_exp_new < x.above.LAI_expanded[j]) {
              x.above.LAI_dead[j] += std::max(0.0, LAI_exp_prev - x.above.LAI_expanded[j]); // Add senescence leaves to dead
            }
          } else {
            x.above.LAI_expanded[j] =x.above.LAI_live[j]*((1.0 - (1.0/x.paramsPhenology.leafDuration[j])) + (x.internalPhenology.phi[j]/x.paramsPhenology.leafDuration[j]));
          }
        }
      } else {
        //Apply defoliation effects to evergreens
        if(x.control.defoliation.cavitationInducedDefoliation) {
          double LAI_exp_prev= x.above.LAI_expanded[j]; //Store previous value
          double LAI_exp_new = x.above.LAI_live[j]*(1.0 - proportionDefoliationWeibull_c(psiLeafPLC, x.paramsTranspiration.VCleafapo_c[j], x.paramsTranspiration.VCleafapo_d[j], x.control.defoliation.criticalLeafPLC, x.control.defoliation.cvLeafP50));
          if(LAI_exp_new < x.above.LAI_expanded[j]) {
            x.above.LAI_expanded[j] = LAI_exp_new;
            x.above.LAI_dead[j] += (LAI_exp_prev - x.above.LAI_expanded[j]); // Add senescence leaves to dead
          }
        } else {
          x.above.LAI_expanded[j] = x.above.LAI_live[j];
        }
      }
    }
  }    
}