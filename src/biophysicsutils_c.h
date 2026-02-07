
#ifndef BIOPHYSICS_UTILS_C_H
#define BIOPHYSICS_UTILS_C_H

const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1
const double SIGMA_Wm2 = 5.67*1e-8; //Stefan-Boltzmann constant W * K^-4 * m^-2
const double defaultLambda = 546.6507; // Wavelength in nm for irradiance to photon flux conversion

double airDensity_c(double temperature, double Patm);
double atmosphericPressure_c(double elevation);
double averageDaylightTemperature_c(double Tmin, double Tmax);
double averageDailyVapourPressure_c(double Tmin, double Tmax, double RHmin, double RHmax);
double leafTemperature_c(double absRad, double airTemperature, double u, double E,  double leafWidth);
double leafTemperature2_c(double SWRabs, double LWRnet, double airTemperature, double u, double E,  double leafWidth);

double leafVapourPressure_c(double leafTemp,  double leafPsi);

double temperatureDiurnalPattern_c(double t, double tmin, double tmax, 
                                 double tminPrev, double tmaxPrev, double tminNext, double daylength);
double radiationDiurnalPattern_c(double t, double daylength);
double irradianceToPhotonFlux_c(double I, double lambda);
double saturationVapourPressure_c(double temperature);
double waterDynamicViscosity_c(double temp);
double latentHeatVaporisation_c(double temperature);
double PenmanPET_c(double latrad, double elevation, double slorad, double asprad, int J,
                 double Tmin, double Tmax, double RHmin, double RHmax, double R_s,
                 double u, double z = 10.0, double z0 = 0.001,
                 double alpha = 0.25, std::string windfun = "1956");
#endif
