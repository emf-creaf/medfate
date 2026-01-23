
#ifndef BIOPHYSICS_UTILS_C_H
#define BIOPHYSICS_UTILS_C_H

const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1
const double SIGMA_Wm2 = 5.67*1e-8; //Stefan-Boltzmann constant W * K^-4 * m^-2
const double defaultLambda = 546.6507; // Wavelength in nm for irradiance to photon flux conversion

double leafTemperature_c(double absRad, double airTemperature, double u, double E,  double leafWidth);
double leafTemperature2_c(double SWRabs, double LWRnet, double airTemperature, double u, double E,  double leafWidth);

double leafVapourPressure_c(double leafTemp,  double leafPsi);

double temperatureDiurnalPattern_c(double t, double tmin, double tmax, 
                                 double tminPrev, double tmaxPrev, double tminNext, double daylength);
double radiationDiurnalPattern_c(double t, double daylength);
double irradianceToPhotonFlux_c(double I, double lambda);

double waterDynamicViscosity_c(double temp);

#endif
