#ifndef GROWTH_DAY_C_H
#define GROWTH_DAY_C_H

double dailyMortalityProbability_c(double stressValue, double stressThreshold);
double phloemFlow_c(double psiUpstream, double psiDownstream,
                    double concUpstream, double concDownstream,
                    double temp, double k_f, double nonSugarConc);
double qResp_c(double Tmean);
#endif