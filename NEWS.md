-------------------------------
 NEWS for R Package "medfate"
-------------------------------

# Version 4.0.0
- Replacement of vulnerability curve parameters
- New parameters for Jarvis-type stomatal conductance in Cochard sub-model
- Revision of Cochard sub-model

# Version 3.2.0
- Correction of a bug arisen in old Rcpp versions
- Functions for seed production and bank dynamics (new parameter SeedLongevity)
- New parameters for dispersal (SeedMass, DispersalType, DispersalDistance, DispersalShape)
- Recording of leaf PLC for Sperry and Cochard sub-models

# Version 3.1.4
- Fixing memory access errors

# Version 3.1.3
- New option 'months' in simulation summaries.
- Leaf area, foliar biomass and fuels of shrubs and herbs depend on the leaf area of plants above
- Survival model based on basal area available as an alternative to fixed mortality baseline rates
- New control thresholds 'minimumShrubCohortCover' and 'minimumTreeCohortDensity' for cohort removal
- Tree bark thickness parameters added to SpParamsMED

# Version 3.1.1
- Fire severity implemented
- Crown bud percent decreases with PLC and fire, increases with SA growth and regulates primary growth rate
- Forest objects can now have additional variables defined (LAI, foliar biomass, fuel loading), to override estimations from allometric models
- New example forest object 'exampleforestMED2'

# Version 3.1.0
- Herbaceous information in summary.forest
- Herbaceous transpiration added to the water balance
- Weather input can now have dates as Date or POSIXct objects in column 'dates' instead of row names
- Weather input can now have missing values in MinRelativeHumidity, MaxRelativeHumidity and Radiation

# Version 3.0.0
- Leaf area limited by allometries
- Functions 'growthInput()' and 'spwbInput()' no longer visible at the user level
- No calculation mode forest structure (US mode removed)
- Control flag options to restrict output in 'spwb()' 'pwb()' and 'growth()'
- Control flag to calculate fire hazard during 'spwb()' 'pwb()' and 'growth()' simulations
- Integration of SurEau-Ecos v2.0 code as new transpiration mode ('Cochard')
- Bug correction tissue moisture and water balance in granier's model
- Tree biomass allometries revised. 
- Tree foliar biomass corrected for high density. Shrub foliar biomass/fuel limited by tree basal area.

# Version 2.9.3
- Bounded leaf area index
- IFNcodes out of SpParams
- PlantPsi lower limit set to -40 MPa
- Psi_Critic no longer used in Granier submodel (Stem vulnerability curve used instead)
- New parameter Exp_Extract to model transpiration decrease in Granier submodel
- Update growth/recruitment parameters
- Removed fordyn dependency on input PET

# Version 2.9.1
- Allows species strings as input in forests objects
- Clean code in growth.cpp

# Version 2.8.3
- Replacing sprintf calls
- Update of forest_mergeShrubs and forest_mergeTrees

# Version 2.8.2
- Function 'redefineSoilLayers' moved to package 'medfateutils'
- Elements 'ID' and 'patchsize' removed from 'forest'
- Bug correction in defaultManagementFunction for 'above-systematic' and 'below-systematic'
- New management option 'targetTreeSpecies'
- Fire hazard estimation includes dead canopy fuels (Ruffault et al. 2023)

# Version 2.8.1
- New meta-modelling parameters
- New growth/senescence parameters

# Version 2.8.0
- New control parameter 'subdailyCarbonBalance', applying to growth simulations with transpiration = "Sperry"
- Hydraulic redistribution implemented for transpiration = "Granier"
- Growth with "Sperry" using canopy temperature

# Version 2.7.7
- Imputation relationship between RGRcambiummax and SRsapwood
- Dynamic modification of LeafPI0 and StemPI0 removed
- Non-stomatal limitations to photosynthesis removed

# Version 2.7.6
- MeanTemperature not longer an input to medfate. It is calculated from minimum and maximum temperature.
- MeanRelativeHumidity not longer an input to medfate. 
- PET not longer an input to medfate. It is calculated internally.
- New simulation parameter 'CO2ByYear' to specify year by year variations in atmospheric CO2
- Sensitivity of photosynthesis to VPD and CO2 concentration under Granier's model
- New species parameters 'WUE_co2', 'WUE_vpd' to regulate the effect of CO2 concentration and VPD under Granier's model
- Species parameter 'WUE_decay' renamed 'WUE_par'
- Control parameter 'Catm' renamed 'defaultCO2'
- New output data frame 'CarbonBalance' for growth 

# Version 2.7.5
- Maximum stem conductance to avoid overestimation of stem conductance in small shrubs
- Bug correction in summary.forest
- Allows filling missing Z50/Z95 values from SpParams when creating model inputs

# Version 2.7.4
- Reducing unused parameters for Sperry model
- Reducing computational time for Sperry model
- Function fireHazard now accepts objects spwb_day and growth_day

# Version 2.7.3
- Cleaning for CRAN
- New article to prepare model inputs

# Version 2.7.2
- SpParamsMED include the results of parameter estimation exercises
- RGRsapwoodmax and RGRcambiummax regulate sapwood formation for shrubs/trees, respectively 
- Simplified sapwood growth (no ring of cells)

# Version 2.7.1
- Revision of recruitment model, with the addition of a recruitment probability (ProbRecr) within the bioclimatic limits
- Temperature effects on sapwood conversion to heartwood
- Correction of estimation of Psi_Extract from turgor loss point
- Estimation of RSSG from shade tolerance
- Minimum DBH parameter in stand metrics

# Version 2.7.0
- Basic water balance model with relative water content
- 'pRootDisc' eliminated from species parameters
- Plant water balance and cuticular transpiration added to the basic water balance model
- Water pools revised for the basic and water balance models
- Shared water pools is now controlled via parameter 'rhizosphereOverlap'
- New species-specific parameter 'MortalityBaselineRate'
- Live fuel moisture content now included in simulation results
- New species-specific parameter 'RSSG' (minimum relative starch for sapwood growth)

# Version 2.6.2
- Control parameter 'modifyInput' is no longer available. Functions 'spwb', 'pwb' and 'growth' do not modify input objects and return an element 'spwbOutput' or 'growthOutput' with a copy of the final state.

# Version 2.6.1
- Tree cover (open grown assumption)
- Parameters 'ShrubCover' and 'CanopyCover' no longer required in function fuel_FCCS
- Input weather stored in output from functions 'spwb', 'pwb' and 'growth'
- New function 'fireHazard'.
- Relative bias and relative MAE as evaluation metrics
- Drought-related leaf senescence only occurring if 'StemPLC' increases
- Embolized sapwood proportion discounted from sapwood maintenance respiration 
- Evaluation of diameter increment (DI), DBH and Height series
- Output growth rates not relative to sapwood area


# Version 2.6.0
- Nitrogen content for leaves, sapwood and fine roots added. 'Nleaf' replaces 'Narea' as the latter can be calculated from 'Nleaf' using 'SLA'.
- Maintenance respiration rates based on N concentration of tissues

# Version 2.5.0
- spwb model with Granier transpiration now extracts water from soil layer according to unsaturated conductivity.
- shinyplot generic function.
- Update parameters fHDmin and fHDmax.
- New vignette 'IFNEvaluation'
- New parameter 'WUE_decay' for reduction of relative WUE in Granier's model
- Tissue construction costs are now species-specific parameters.
- Fine root growth and senescence made equal between Granier and Sperry models, inducing a new parameter 'Ar2Al' for Granier's model.

# Version 2.4.0
- Functions 'spwb_stress', 'spwb_resistances' and 'spwb_waterUseEfficiency' renamed to 'droughtStress', 'resistances' and 'waterUseEfficiency', respectively, since they can now be applied to the output of several simulation function.
- Plant biomass balance in growth. Structural changes are now daily in growth simulations. 
- Bug correction in shrub structural update. Shrub dynamics default set to TRUE.
- DBH/Height plots from growth output.
- IMPORTANT: New species parameter table.

# Version 2.3.8
- Summary functions revised, including new function 'summary.fordyn'
- Collating intra-annual 'fordyn' results for plotting.
- Forest management enabled in function 'fordyn' and default management actions defined in 'defaultManagementFunction' and 'defaultManagementArguments'
- DOY, Photoperiod and JulianDay can be taken from weather input in functions 'spwb', 'pwb' and 'growth'. 

# Version 2.3.7
- CO2 made an daily input weather variable, in addition to the default control parameter 'Catm'
- New function 'shinyplot' to create interactive graphics
- New function 'plot.fordyn' to display annual (step) summaries of forest dynamics

# Version 2.3.5
- New example vignette 'FontBlanche'
- Modification of evaluation functions to separate the evaluation of total evapotranspiration ('ETR') from the evaluation of soil evaporation + plant transpiration ('SE+TR')
- Bug correction: sub-daily stomatal conductance plots

# Version 2.3.4
- Maximum relative sapwood growth rates effective
- New option 'summary.freq' in plot.spwb and plot.growth
- New species parameters (sapwood and fine root senescence rates)
- Defaults for 'conduit2sapwood' from taxonomical family
- IMPORTANT: New species parameter tables (including estimates for conduit2sapwood)

# Version 2.3.2
- Parameter 'ParticleDensity' eliminated, as it is now calculated from 'LeafDensity', 'WoodDensity' and 'r635'

# Version 2.3.1
- IMPORTANT: New species parameter tables
- New control parameter 'fillMissingSpParams'
- Defaults for 'LeafDensity', 'WoodDensity', 'LeafPI0', 'LeafEPS' and 'LeafAF' from taxonomical family.
- New species categorical params (for inbuilt imputation) 'LeafShape' and 'LeafSize'
- Parameter 'Flammability' index eliminated (non-meansurable property of flammability).
- Parameter 'LeafLitterFuelType' eliminated because it is derived from 'LeafShape' and 'LeafSize'.
- Defaults added for 'r635', 'heatContent', 'LigninPercent' and 'SAV' according to 'LeafShape' and 'LeafSize'
- New function 'getSpParamsDefinition()' returns definition of species parameters.
- Inbuilt defaults added for shrub allometries, depending on 'LifeForm' and 'Hmax', and for tree allometries, depending on 'Group'.
- Default value added for 'pDead'.

# Version 2.2.3
- Functions soilgridsParams() and forest_map*Tables() moved to package 'medfateutils' available at GitHub (emf-creaf/medfateutils).

# Version 2.2.2
- Dependency 'spdep' removed
- Bug correction 'windKatul.cpp'

# Version 2.2.1
- Calibrated minimum bioclimatic parameters for recruitment (SpParamsMED)
- Explicit species input parameters for phenology

# Version 2.2.0
- New simulation function 'fordyn()', including recruitment process
- New function 'mergeShrubs'
- New functions 'species_parameter' and 'species_characterParameter'

# Version 2.1.4
- Dessication/defoliation homogenized across transpiration modes in growth function
- Revision of mortality (stochastic/deterministic, whole-cohort/density)

# Version 2.1.3
- Update Psi_Extract according to Psi_TLP (Hydratry)

# Version 2.1.2
- New control flags for defoliation/starvation/dessication in growth simulations
- New control flag for sink limitation in growth simulations
- Revision of phenology submodel

# Version 2.1.1
- Respiration rate for leaves made optionally species-specific.
- Update shrub allometries from De Caceres et al. (2019).
- LAI_live and LAI_expanded in growth.
- Update SpParamsUS (missing values for new parameters).

# Version 2.1.0
- Parameters of Granier's equation made species-specific if available
- Parameter modification dependencies revised
- Modification of input objects is now optional
- Update of function `transp_maximumTranspirationModel`
- Photosynthesis in Granier's model corresponds to gross photosynthesis and is proportional to transpiration

# Version 2.0.1
- Cloning initial object for optimization
- Leaf growth costs always drawn from sugar sapwood
- Bug correction: LWRnet calculation
- Bug correction: Sperry model does not crash when LAIstand = 0

# Version 2.0.0
- IMPORTANT: Soil input merged with model input. Now the calls to simulation functions (e.g., `spwb`, `growth`) do not need to include soil as input parameter.
- Bug correction: NaN values for theta > theta_sat in van Genuchten psi computation
- Function `modifyInputParams` now accepts modification of soil layer properties
- Sensitivity/Calibration vignette updated

# Version 1.1.6
- New canopy turbulence models by Katul et al (2004).
- Long-wave net radiation balance for layered canopies following Flerchinger et al. (2009).
- Multi-layer canopy energy balance as in Bonan et al. (2014)
- New option 'depthMode' to calculate fuel depth in 'fuel_FCCS'.
- Modification of 'vprofile_windExtinction' to draw turbulence models.

# Version 1.1.5
- Stomatal conductance now denoted as Gsw or GSW
- Boundary layer conductance considered
- Leaf water potential influencing leaf vapour pressure

# Version 1.1.4
- Flexible temporal resolution of model evaluation
- Basal area index evaluation

# Version 1.1.3
- New model evaluation functions
- New optimization function factories
- New function 'modifyInputParams'
- Control parameters set to a nested list
- New vignette for calibration and sensitivity analysis

# Version 1.1.2
- GW sunlit/shade minimum/maximum daily output
- FMC calculations with basic water balance output
- New output in growth simulation (biomass values)
- New function 'moisture_cohortFMCDay'
- Modification of Z50/Z95

# Version 1.1.1
- Root exudation added to carbon balance
- Revision sapwood growth
- Growth cost for fine roots in basic model
- Translocation for carbon during senescence
- Bug correction in fuel calculations with US mode

# Version 1.1.0
- Control option 'rockyLayerDrainage' instead of 'drainage' to disable macropore vertical outflow in layers with > 95% of rocks
- Soil parameter Kdrain for saturated vertical hydraulic conductivity towards groundwaters (deep drainage)
- Improved validation plots with confidence intervals

# Version 1.0.3
- Nash-Sutcliffe efficiency (NSE) implemented in spwb_validation

# Version 1.0.2
- Bug correction on the use of organic matter in Saxton (2006) equations (thanks to Milan Fischer).
- Recodification of soilgridParams due to new SoilGrids REST API (removed dependency from GSIF).
- New function 'redefineSoilLayers'

# Version 1.0.1
- New root functions
- Bulk density stored in soil object initialization
- Advanced plant water pools
- Belowground inputs restructured

# Version 1.0.0
- Reorganization of growth function
- Clarification of gross and net photosynthesis
- 'spwb_resetInputs' to 'resetInputs'
- Dependence of kmax on temperature (due to sap dynamic viscosity) incorporated
- Functions plot.spwb and plot.pwb modified to draw subdaily dynamics for a subset of dates 

# Version 0.9.1
- Small bug fixes
- update 'spwb_resetInputs'
- Update of plant water pools
- New output (annual stand summaries and aboveground structure) for function 'growth'

# Version 0.9.0
- New parametrization data set 'SpParamsUS'
- Function 'hydrology_verticalInputs' replaced by 'hydrology_soilWaterInputs' and 'hydrology_soilInfiltrationPercolation'.
- New simulation control option: 'plantWaterPools'.
- Hard (Imports) dependency from GSIF changed to soft (Suggests) one.

# Version 0.8.9
- Update supply function plot.

# Version 0.8.8
- New function 'soil_rockWeight2Volume'

# Version 0.8.7
- Corrections to energy balance for zero LAI (deciduous species)
- SFI functions moved to medfateland
- soilgridParams modified to accept a SpatialPoints object as input

# Version 0.8.6
- New function 'spwb_sensitivity' for sensitivity analyses
- New control parameter 'unlimitedSoilWater'
- Bug correction in canopy height with LAI = 0
- Modifications of spwb_ldrOptimization to work with transpirationMode = "Sperry"
- New function spwb_ldrExploration

# Version 0.8.5
- New control parameter 'fracLeafResistance'
- Different control options for parameter 'cavitationRefill'
- New control parameter 'refillMaximumRate'
- Control parameter 'hydraulicCostFunction' replaced by 'costModifier' and 'gainModifier'
- New control parameter 'cuticularTranspiration'
- Numerical controls to avoid NaN in functions 'soil_theta2psiSX' and 'soil_psi2thetaSX'
- Bug correction in estimation of root conductance proportions

# Version 0.8.4
- Water balance console output modified in spwb
- New approach to plant water compartments 'capacitance = TRUE'
- Output of plant water balance
- New option in control
- Stem segments fixed to two
- functional parameter pRootDisc removed from Sperry's advanced model
- Bug correction in fuel_cohortFineFMC
- Remove ksympver and add klatleaf/klatstem to control parameters

# Version 0.8.3
- Reference book (medfatebook) launched
- Fraction of absorbed SWR output in Granier's transpiration
- 'Stand' data frame output in spwb(), separated from 'WaterBalance'
- New function 'forest_mapTreeTable', 'forest_mapShrubTable' and 'forest_mapWoodyTables'
- 'stand_*' functions for stand-level properties
- No SWR soil absorption when snow pack is present in Sperry's model
- verticalLayerSize made a control parameter for Granier's model
- Changes in light parameters: New parameter 'alphaSWR'. 'albedo' renamed to 'gammaSWR'. 'k' renamed to 'kPAR'

# Version 0.8.2
- Modification of meteoland to better calculate direct/diffuse light on slopes
- Added Narea parameter to facilitate estimation of Vmax298
- New function 'spwb_validation'.
- Plot functions using ggplot.
- Improvement of infiltration repartition for varying macroporosity.
- Input values for latitude and topography stored in the result of simulations.
- Wind for each cohort stored in the result of transp_transpirationSperry and spwb_day.

# Version 0.8.1
- Export Ci from spwb_day.
- New function 'maximumTranspirationRatioPlot'
- Output of min/max water potential for sunlit and shade leaves.
- New function 'soil_waterRetentionPlot'
- New function 'waterUseEfficiency'
- Default Van Genuchten PTF set to 'Toth'
- New functions 'modifySpParams' and 'modifyCohortParams'.
- New function 'hydrology_interceptionPlot'

# Version 0.8.0
- Move spatial classes and methods to package 'medfateland'
- Implement underscores instead of dots to separate function groups and function names
- New function transp_Granier.
- Bug corrected in spwb.plot for snow plotting.
- New function 'pwb'.
- New function 'snowMelt'.
- Transpiration model changed to Granier and Sperry.
- PLC set to zero when DOY = 1
- New functions for leaf phenology

# Version 0.7.4
- Fixing bugs for installation in all platforms
- Update docs

# Version 0.7.3
- Shrub root system using LDR (Z50 and Z95)
- Percolation of infiltrated water consistent with layer subdivision
- New function soil.infiltrationRepartition()
- Adding interception to evapotranspiration
- Checked for CRAN


# Version 0.7.2
- Accounting of hydraulic redistribution
- Update plotting functions
- New function vprofile.RootDistribution.
- New functions for water at wilting point (-1.5 MPa)
- Print extractable water of soil
- Bug correction in photosynthesis (now done per leaf area basis)
- Leaf area distribution (and crown fine biomass distribution) following truncated normal [-1.5,1.5]

# Version 0.7.1
- Interception corrected in the complex model
- Etol set to 0.0000001
- spwb export of dEdP (equivalent to soil-plant conductance)
- Subdaily results can be stored for spwb
- New function spwb.resistances to calculate and draw segment resistances for spwb simulation results.
- Cohort parameter search by SpIndex
- Gwmin set to zero when capacitance = FALSE
- New function spwb.stress to calculate drought stress indices
- PlantStress in complex mode now is calculated as relative soil-plant conductance, for compatibility with the simple mode.
- Function name changes for interception and soil hydrology

# Version 0.7.0
- Leaf and stem water compartments added
- Update of functions 'spwb.day' and 'plot.spwb.day'
- Analytical integral of the van Genuchten function
- Analytical inverse of the incomplete gamma function

# Version 0.6.2
-  Adapt to Rcpp changes
- Stem water compartments


# Version 0.6.1
- Saturated theta in Saxton model
- Water table depth added
- Boolean option 'drainage' added to soil water balance
- Added subsurface flow processes to spwbgrid

# Version 0.6.0
- Function name 'swb' (and all related functions) changed to 'spwb'
- Added snow pack to soil state variables.
- Added new control option 'snowpack' for snowpack dynamics simulation (only when transpirationMode='Simple')
- Changed names of spatial classes from 'Forest' to 'Landscape'
- Changed name of 'exampleSPF' to 'exampleSPL'.
- Function 'spwbgrid' functional again.
- Improved documentation of function 'spwb.day'
- 'DOY' no longer needed as column in meteorological input.

# Version 0.5.6
- Growth degree days added as output of 'swb()'
- New function 'swb.resetInputs()'
- New function 'hydraulics.vulnerabilityCurvePlot()'
- Reorganization of help for hydraulics.
- Documentation of tissue moisture functions.
- Update function 'plot.growth()'

# Version 0.5.5
- New functions for tissue moisture
- New function 'fuel.cohortFineFMC'
- New option 'bySpecies' to aggregate results by species in functions 'summary.swb' and 'summary.growth'
- New option 'bySpecies' to aggregate results by species in functions 'plot.swb' and 'plot.swb.day'
- Bulk density added to soil parameters
- New set of pedotransfer functions to calculate VG parameters from texture, bulk density and organic content

# Version 0.5.4
- New functions for pressure-volume curves
- Rmarkdown vignettes
- Added Van Genuchten pedotransfer functions
- New control parameter added (for soil functions)
- New function 'forest2belowground'
- Update 'summary.swb' function

# Version 0.5.3
- Network representation of the continuum now includes a leaf segment
- Stem fraction of total resistance removed as estimation of root xylem conductance
- Default kleaf_max (8 for temperate angiosperms and 5 for gymnosperms)
- New parameter rootxylem_kmax (hydraulic conductivity of roots)

# Version 0.5.2
- Modification of default for ntrial
- Three element supply function added
- Parameters of the leaf maximum conductance and leaf vulnerability curve added

# Version 0.5.1
- Structure of swb.day output
- Increased output (leaf temperature, stomatal conductance and leaf VPD) in swb.day
- Taper modifications
- New function plot.swb.day

# Version 0.5.0
- Atmospheric CO2 is an input control parameter
- Added Hmed as species parameter (to correct reference conductivity values)
- Modified documentation

# Version 0.4.9
- Alternative way of calculating maximum root conductance
- Default value for averageFracRhizosphereResistance changed to 0.15
- Root vulnerability curve parameters taken from stem vulnerability curve when missing (d_root = d_stem/2)

# Version 0.4.8
- Leaf radiation balance with LWR from soil
- New output for swb
- Bug fix in plot.swb (match according to transpiration model)

# Version 0.4.7
- Profit maximization for sunlit and shade leaves separately
- Bug correction in swb
- Check on stomatal conductances in profit maximization

# Version 0.4.6
- Radiation absorbed by trunks
- Separation of energy balance components

# Version 0.4.5
- New control parameter: Canopy thermal capacity per LAI
- Update of calls to meteoland (diffuse radiation)
- Wind value when missing data added to control

# Version 0.4.4
- Minimum windspeed (1.0 m/s)
- Progressive leaf fall
- Bug correction (swb.plot)
- Energy balance output/plot
- Modified conductance scaling

# Version 0.4.3
- Added new parameter 'LeafWidth'

# Version 0.4.2

- Design changes for radiation balance (soil/canopy/atmosphere)
- Added temperature balance in summary.swb and plot.swb
- Added new parameter 'albedo'

# Version 0.4.1

- Documentation of soil thermodynamics
- Energy balance output in swb.day

# Version 0.4.0

- Added soil temperature state variables (for 'complex' mode)
- Added canopy temperature state variable (for 'complex' mode)
- 'gdd' in swbInput used as initial growth degree days
- Soil thermodynamics.

# Version 0.3.7

- Added vertical layer size as option in control.
- New function 'soilgridsParams' (calls functions in package 'GSIF')
- Added longwave radiation to leaf energy balance.

# Version 0.3.6

- Added new pedotransfer functions with organic matter.
- New function added 'soil.waterFC'.
- Xylem taper added as option in control.
- Hydraulic cost function added as option in control.
- New function 'light.instantaneousLightExtinctionAbsortion'.
- New functions 'transp.dayCanopyTranspiration' and 'transp.dayCanopyTranspirationPlot'.

# Version 0.3.5

- Control of numerical methods for supply function added to 'defaulControlParams'.
- Added new function 'swb.ldrCalibration' to calibrate root distribution for swb simulations (by Victor Granda)
- Added documentation for light extinction functions 'light.layerIrradianceFraction', 'light.layerSunlitFraction' and 'light.cohortSunlitShadeAbsorbedRadiation'.
