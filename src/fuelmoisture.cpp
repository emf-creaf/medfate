#include <Rcpp.h>
#include "hydraulics.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include <meteoland.h>
using namespace Rcpp;



/**
 * From: Byram & Jemison (1943). See also Viney (1991).
 * airTemp - air temperature (in degrees Celsius)
 * airHumidity - air relative humidity (in percent)
 * fuelRadiation - instant solar radiation reaching fuels (in W/m2)
 * fuelWindSpeed - wind speed at fuel surface (in m/s)
 */
// [[Rcpp::export(".fuelConditions")]]
NumericVector fuelConditions(double airTemp, double airHumidity, double fuelRadiation, double fuelWindSpeed) {
  double fuelTemp = airTemp + fuelRadiation/(42.5*fuelWindSpeed + 32.7);
  double fuelHumidity = airHumidity*exp(0.059*(airTemp-fuelTemp));
  return(NumericVector::create(fuelTemp, fuelHumidity));
}

/**
 * From: Van Wagner & Pickett (1985). See also Viney (1991).
 * fuelTemperature - fuel temperature (in degrees Celsius)
 * fuelHumidity - fuel relative humidity (in percent)
 */
// [[Rcpp::export(".EMCdesorption")]]
double EMCdesorption(double fuelTemperature, double fuelHumidity) {
  return(0.942*pow(fuelHumidity, 0.679)+0.000499*exp(0.1*fuelHumidity)+0.18*(21.1-fuelTemperature)*(1.0-exp(-0.115*fuelHumidity)));
}
// [[Rcpp::export(".EMCadsorption")]]
double EMCadsorption(double fuelTemperature, double fuelHumidity) {
  return(0.618*pow(fuelHumidity, 0.753)+0.000454*exp(0.1*fuelHumidity)+0.18*(21.1-fuelTemperature)*(1.0-exp(-0.115*fuelHumidity)));
}

/**
 *  From: Simard (1968). See also Viney (1991) [metric coefficients].
 */
// [[Rcpp::export(".EMCSimard")]]
double EMCSimard(double fuelTemperature, double fuelHumidity) {
  if(fuelHumidity<10.0) {
    return(0.03+(0.2626*fuelHumidity) - (0.001040*fuelTemperature*fuelHumidity));
  }
  if(fuelHumidity < 50.0) {
    return(1.76+ (0.1601*fuelHumidity) - (0.02660*fuelTemperature));
  }
  return(21.06 - (0.4944*fuelHumidity) + (0.005565*fuelHumidity*fuelHumidity) - 0.00063*fuelTemperature*fuelHumidity);
}

/**
 * Fine Fuel Moisture Code (FFMC) of the Canadian Forestry Service (1984). See also Viney (1991).
 * 
 * m0 - initial fuel moisture content (in percent of dry weight)
 * fuelTemp - fuel temperature (in degrees Celsius)
 * fuelHumidity - fuel relative humidity (in percent)
 * fuelWindSpeed - wind speed at fuel surface (in m/s)
 * effRain - precipitation (in mm) after accounting for canopy interception
 */
double fine1hday(double m0, double fuelTemp, double fuelHumidity, double fuelWindSpeed, double effRain) {
  //apply first absorption from rain
  double delta = 0.0;
  if(m0>150.0) delta = 1.0;
  double mr =  m0 + 42.5*effRain*exp(100.0/(m0-251))*(1.0-exp(-6.93/effRain))+delta*(0.0015*sqrt(effRain)*pow(m0-150.0,2.0));
  mr = std::min(250.0, mr);
  //apply adsorption/desorption
  double Ea = EMCadsorption(fuelTemp, fuelHumidity);
  double Ed = EMCdesorption(fuelTemp, fuelHumidity);
  double EMC = 0.0;
  double eta = 0.0;
  if(mr<Ea) {
    EMC = Ea;
    eta = 1.0 - (fuelHumidity/100.0);
  } 
  else if(mr>Ed) {
    EMC = Ed;
    eta = (fuelHumidity/100.0);
  } else { //No moisture variation if Ea < mr < Ed
    return(mr);
  }
  double k0 = (0.567*(1.0-pow(eta,1.7))+0.176*sqrt(fuelWindSpeed)*(1.0 - pow(eta,8.0)))*exp(0.0365*fuelTemp);
  return(EMC  + (mr - EMC)*exp(-k0));
}

/**
 * Coarse 10-h, modified from Bradshaw et al. (1983) See Ruiz (2004).
 * 
 * m0 - initial fuel moisture content (in percent of dry weight)
 * prevFuelTempMax - Previous day maximum fuel temperature (in degrees Celsius)
 * prevFuelHumidityMin - Previous day minimum fuel relative humidity (in percent)
 * currFuelTempMin - Current day minimum fuel temperature (in degrees Celsius)
 * currFuelHumidityMax - Current day maximum fuel relative humidity (in percent)
 * rainDuration - rain duration (in hours)
 */
double coarse10hday(double m0, 
                    double prevFuelTempMax, double prevFuelHumidityMin,
                    double currFuelTempMin, double currFuelHumidityMax,
                    double rainDuration) {
  // EMC1: Equilibrium moisture content calculated from maximum temperature and maximum relative humidity
  // of previous day
  double EMC1 = EMCSimard(prevFuelTempMax, prevFuelHumidityMin);
  // EMC2: Equilibrium moisture content calculated from minimum temperature and maximum relative humidity 
  // of current day
  double EMC2 = EMCSimard(currFuelTempMin, currFuelHumidityMax);
  // Rain duration for the two periods
  double rainDuration1 = rainDuration*(16.0/24.0);
  double rainDuration2 = rainDuration*(8.0/24.0);
  // Average boundary conditions for the two periods
  double D1 = ((16.0-rainDuration1)*EMC1+(2.7*rainDuration1+76.0)*rainDuration1)/16.0;
  double D2 = ((8.0-rainDuration2)*EMC2+(2.7*rainDuration2+76.0)*rainDuration2)/8.0;
  // Moisture predictions
  double m1 = m0 + (D1-m0)*(1.0-1.11*exp(-1.6));
  return(m1  + (D2 - m1)*(1.0-0.87*exp(-0.8)));
}

/**
 * Coarse 100-h, modified from Bradshaw et al. (1983) See Ruiz (2004).
 * 
 * m0 - initial fuel moisture content (in percent of dry weight)
 * fuelTempMin - Minimum fuel temperature (in degrees Celsius)
 * fuelHumidityMax - Maximum fuel relative humidity (in percent)
 * fuelTempMax - Maximum fuel temperature (in degrees Celsius)
 * fuelHumidityMin - Minimum fuel relative humidity (in percent)
 * numSunHours - Number of daylight hours (from 0 to 24)
 * rainDuration - rain duration (in hours)
 */
double coarse100hday(double m0, 
                     double fuelTempMin, double fuelHumidityMax, 
                     double fuelTempMax, double fuelHumidityMin, 
                     double numSunHours, double rainDuration) {
  // EMC1: Equilibrium moisture content calculated from minimum temperature and maximum relative humidity
  double EMC1 = EMCSimard(fuelTempMin, fuelHumidityMax);
  // EMC2: Equilibrium moisture content calculated from maximum temperature and minimum relative humidity  
  double EMC2 = EMCSimard(fuelTempMax, fuelHumidityMin);
  // EMCmean: Mean equilibrium moisture content (gives more weight to EMC1 in winter and to EMC2 in summer)
  double EMCmean = (numSunHours*EMC2 + (24.0-numSunHours)*EMC1)/24.0;
  // D: Moisture average content (including rain effects)
  double D = ((24.0-rainDuration)*EMCmean+(0.5*rainDuration+41.0)*rainDuration)/24.0;
  return(m0  + (D - m0)*(1.0-0.87*exp(-0.24)));
}


double layerLiveFuelMoisture(double minHeight, double maxHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  double num = 0.0, den = 0.0, pin;
  int nCoh = cohortLoading.size();
  for(int i=0;i<nCoh; i++) {
    pin = (std::min(H[i], maxHeight)-std::max((1.0-CR[i])*H[i], minHeight))/(CR[i]*H[i]);
    if(pin<0.0) pin = 0.0;
    num +=(cohortFMC[i]*cohortLoading[i]*pin);
    den += (cohortLoading[i]*pin);
    // Rcout<<cohortFMC[i]<< " "<<H[i]<<" "<<maxHeight<<" "<< pBole[i]*H[i]<< " "<<minHeight<< ": "<<pin<<"\n";
  }
  if(den>0) return(num/den);
  return(NA_REAL);
}

double canopyLiveFuelMoisture(double canopyBaseHeight, double canopyTopHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  return(layerLiveFuelMoisture(canopyBaseHeight, canopyTopHeight, cohortFMC, cohortLoading, H, CR));
}

double fuelbedLiveFuelMoisture(double fuelbedHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR) {
  return(layerLiveFuelMoisture(0, fuelbedHeight, cohortFMC, cohortLoading, H, CR));
}
