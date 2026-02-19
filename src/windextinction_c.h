#ifndef WINDEXTINCTION_H
#define WINDEXTINCTION_H

double windSpeedAtCanopyHeight_c(double wind20H, double canopyHeight);
double windSpeedAtHeightOverCanopy_c(double z, double wind20H, double canopyHeight);
double windSpeedMassmanExtinction_c(double z, double wind20H, double LAIc, double canopyHeight);
double unshelteredMidflameWindSpeed_c(double wind20H, double fuelBedHeight);
double shelteredMidflameWindSpeed_c(double wind20H, double crownFillProportion, double topCanopyHeight);
double aerodynamicResistance_c(double canopyHeight, double wind);
#endif
