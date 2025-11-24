# Changelog

## medfate 4.8.5

- Store medfate version in initialized objects, to implement backwards
  compatibility
- Backwards compatibility of objects initialized with ver. 4.8.3 and
  4.8.4
- Checks missing key weather inputs in spwb_day()
- SLA, LeafWidth and r635 imputation for LeafShape ‘Spines’
- Age of cohorts is now allowed in forest objects and incremented in
  fordyn() simulations
- Flag ‘addunits’ in function extract() allows retrieving physical units
  for daily output variables
- Function extract() can now process fordyn() simulations

## medfate 4.8.4

CRAN release: 2025-10-21

- Summary, extract and plot functions enabled for aspwb() simulation
  results
- Bug correction: calculation of fire hazard in calls to growth() and
  fordyn()
- Fire hazard now returns canopy and understory layer loadings
- Drought-driven defoliation now occurs in spwb() simulations
- Revision of LFMC calculation
- Fine root distribution can be truncated via a control parameter,
  leading to Z100 imputation
- Sureau now updates (internally) soil water potential at subdaily steps

## medfate 4.8.3

CRAN release: 2025-08-27

- Soil pool results (REW and psi) enabled
- Areas changed by bars in PET_Precipitation plots
- New plotting function soil_conductivityCurvePlot()
- Rhizosphere overlap is a parameter in function defaultControl()
- Internal sureau structures rebuild

## medfate 4.8.2

- Bug correction detecting date vector objects (Rcpp issue)
- Bug correction in function extract()

## medfate 4.8.1

CRAN release: 2025-04-23

- New functions utils_rockOptimization, utils_ldrOptimization and
  utils_ldrExploration
- Bug correction: minimum cuticular transpiration in Granier
- Bug correction: all species deciduous in Sperry/Sureau

## medfate 4.8.0

CRAN release: 2025-01-28

- Code optimization: LAI distribution
- Communication structures to improve memory usage in medfateland
- New plot type “pointdynamics” in evaluation

## medfate 4.7.0

CRAN release: 2024-10-18

- LAI can be estimated without the competition effect of cohorts of
  larger size
- Growth rates corrected by the ratio between LAI and maximum LAI
  without competition effects
- New control flags for plant-level outputs in growth simulations
- Revision of function resistances() to produce extra results and adapt
  it to Sureau
- New parameters Dmax and SeedProductionDiameter for trees
- Code optimization: communication using data structures
- Update of SpParamsMED parameter values

## medfate 4.6.0

- ObsID column is allowed in treeData and shrubData to indicate cohort
  ID from observation
- Tree and shrub merging allow preserving cohorts with non-missing ObsID
- New control parameter to ask fordyn() to preserve (i.e. not remove or
  merge) cohorts with non-missing ObsID

## medfate 4.5.0

- Bug correction: Zero fine root proportion in Sureau
- Users can now define column Z100 in treeData or shrubData to truncate
  root distribution
- AcceptedName is not a strictly required parameter (only used in
  parameterization)

## medfate 4.4.0

CRAN release: 2024-07-23

- Multi-bucket soil model reintroduced (default) for safety
- Loading offset added as an option to fuel_FCCS

## medfate 4.3.2

- IMPORTANT: Many functions have been internalized (biophysics\_,
  fire\_,…). They are documented and exported, but not listed in
  reference list.
- Default control parameters changed
- Subdaily soil temperature change limited for stability with thin first
  layer
- Column ‘Strict’ added to SpParamsDefinition
- Detection of POSIXct date columns
- New control option ‘lfmcComponent’
- Fixed annual recovery for spwb with sureau transpiration mode
- Stem cuticular transpiration set to FALSE by default
- Imputation of Gs_P50 changed to either VCleaf_P50 (when not missing)
  or else derived from TLP

## medfate 4.3.1

CRAN release: 2024-05-25

- Functions forest2spwbInput/forest2growthInput deprecated. Preferred
  function names are spwbInput/growthInput.

## medfate 4.3.0

- IMPORTANT: soil objects are now data frames
- SWE now stored as ‘snowpack’ in spwbInput
- print.soil renamed to summary.soil
- forest_mapWoodyTables moved from medfateutils into medfate
- New function forest_reduceToDominant

## medfate 4.2.0

- Reorganization of Soil output as a list
- Snow output as independent data frame
- New function soil_redefineLayers
- Bug correction in forest_mergeShrubs
- Option segmentedXylemVulnerability set to FALSE for Sureau

## medfate 4.1.0

- Single-domain and dual-permeability soil water flows
- Infiltration excess, saturation excess and capillarity rise reported
- Lateral water flows and water table depth as inputs
- Soil output revision
- Evaluation of moisture in multiple soil layers

## medfate 4.0.0

- Options leafCavitationEffects and stemCavitationEffects for Sperry
- New taxonomic treatment (Name/AcceptedName)
- Revision of Sureau sub-model
- Replacement of vulnerability curve parameters by P12, P50 and P88
- New parameters for Jarvis-type stomatal conductance in Sureau
  sub-model
- XFT traits for vulnerability curves
- Recording of leaf PLC for all sub-models
- Richard-type soil water movement
- Soil hydrology revised
- New functions for spwb in agricultural lands
- New function ‘extract’
- New infiltration model Green-Ampt (1911)
- Interception model Liu (2001) accepted as alternative to Gash (1995)
- Optional “RainfallIntensity” (mm/h) as input in weather
- New option “defaultRainfallIntensityPerMonth”
- New functions “tree2forest”/“shrub2forest” to create forest objects
  with single cohorts
- Van Genuchten parameters and Ksat can be specified directly when
  initializing soil
- Saturation excess in soil flows
- Evaluation of GPP, H and LE
- Soil thermodynamics revised

## medfate 3.2.0

CRAN release: 2023-11-30

- Correction of a bug arisen in old Rcpp versions
- Functions for seed production and bank dynamics (new parameter
  SeedLongevity)
- New parameters for dispersal (SeedMass, DispersalType,
  DispersalDistance, DispersalShape)
- Recording of leaf PLC for Sperry and Sureau sub-models

## medfate 3.1.4

CRAN release: 2023-08-29

- Fixing memory access errors

## medfate 3.1.3

CRAN release: 2023-08-23

- New option ‘months’ in simulation summaries.
- Leaf area, foliar biomass and fuels of shrubs and herbs depend on the
  leaf area of plants above
- Survival model based on basal area available as an alternative to
  fixed mortality baseline rates
- New control thresholds ‘minimumShrubCohortCover’ and
  ‘minimumTreeCohortDensity’ for cohort removal
- Tree bark thickness parameters added to SpParamsMED

## medfate 3.1.1

- Fire severity implemented
- Crown bud percent decreases with PLC and fire, increases with SA
  growth and regulates primary growth rate
- Forest objects can now have additional variables defined (LAI, foliar
  biomass, fuel loading), to override estimations from allometric models
- New example forest object ‘exampleforest2’

## medfate 3.1.0

- Herbaceous information in summary.forest
- Herbaceous transpiration added to the water balance
- Weather input can now have dates as Date or POSIXct objects in column
  ‘dates’ instead of row names
- Weather input can now have missing values in MinRelativeHumidity,
  MaxRelativeHumidity and Radiation

## medfate 3.0.0

- Leaf area limited by allometries
- Functions ‘growthInput()’ and ‘spwbInput()’ no longer visible at the
  user level
- No calculation mode forest structure (US mode removed)
- Control flag options to restrict output in ‘spwb()’ ‘pwb()’ and
  ‘growth()’
- Control flag to calculate fire hazard during ‘spwb()’ ‘pwb()’ and
  ‘growth()’ simulations
- Integration of SurEau-Ecos v2.0 code as new transpiration mode
  (‘Sureau’)
- Bug correction tissue moisture and water balance in granier’s model
- Tree biomass allometries revised.
- Tree foliar biomass corrected for high density. Shrub foliar
  biomass/fuel limited by tree basal area.

## medfate 2.9.3

CRAN release: 2023-03-11

- Bounded leaf area index
- IFNcodes out of SpParams
- PlantPsi lower limit set to -40 MPa
- Psi_Critic no longer used in Granier submodel (Stem vulnerability
  curve used instead)
- New parameter Exp_Extract to model transpiration decrease in Granier
  submodel
- Update growth/recruitment parameters
- Removed fordyn dependency on input PET

## medfate 2.9.1

CRAN release: 2023-01-08

- Allows species strings as input in forests objects
- Clean code in growth.cpp

## medfate 2.8.3

- Replacing sprintf calls
- Update of forest_mergeShrubs and forest_mergeTrees

## medfate 2.8.2

- Function ‘redefineSoilLayers’ moved to package ‘medfateutils’
- Elements ‘ID’ and ‘patchsize’ removed from ‘forest’
- Bug correction in defaultManagementFunction for ‘above-systematic’ and
  ‘below-systematic’
- New management option ‘targetTreeSpecies’
- Fire hazard estimation includes dead canopy fuels (Ruffault et
  al. 2023)

## medfate 2.8.1

- New meta-modelling parameters
- New growth/senescence parameters

## medfate 2.8.0

CRAN release: 2022-09-14

- New control parameter ‘subdailyCarbonBalance’, applying to growth
  simulations with transpiration = “Sperry”
- Hydraulic redistribution implemented for transpiration = “Granier”
- Growth with “Sperry” using canopy temperature

## medfate 2.7.7

- Imputation relationship between RGRcambiummax and SRsapwood
- Dynamic modification of LeafPI0 and StemPI0 removed
- Non-stomatal limitations to photosynthesis removed

## medfate 2.7.6

- MeanTemperature not longer an input to medfate. It is calculated from
  minimum and maximum temperature.
- MeanRelativeHumidity not longer an input to medfate.
- PET not longer an input to medfate. It is calculated internally.
- New simulation parameter ‘CO2ByYear’ to specify year by year
  variations in atmospheric CO2
- Sensitivity of photosynthesis to VPD and CO2 concentration under
  Granier’s model
- New species parameters ‘WUE_co2’, ‘WUE_vpd’ to regulate the effect of
  CO2 concentration and VPD under Granier’s model
- Species parameter ‘WUE_decay’ renamed ‘WUE_par’
- Control parameter ‘Catm’ renamed ‘defaultCO2’
- New output data frame ‘CarbonBalance’ for growth

## medfate 2.7.5

- Maximum stem conductance to avoid overestimation of stem conductance
  in small shrubs
- Bug correction in summary.forest
- Allows filling missing Z50/Z95 values from SpParams when creating
  model inputs

## medfate 2.7.4

- Reducing unused parameters for Sperry model
- Reducing computational time for Sperry model
- Function fireHazard now accepts objects spwb_day and growth_day

## medfate 2.7.3

CRAN release: 2022-05-09

- Cleaning for CRAN
- New article to prepare model inputs

## medfate 2.7.2

- SpParamsMED include the results of parameter estimation exercises
- RGRsapwoodmax and RGRcambiummax regulate sapwood formation for
  shrubs/trees, respectively
- Simplified sapwood growth (no ring of cells)

## medfate 2.7.1

- Revision of recruitment model, with the addition of a recruitment
  probability (ProbRecr) within the bioclimatic limits
- Temperature effects on sapwood conversion to heartwood
- Correction of estimation of Psi_Extract from turgor loss point
- Estimation of RSSG from shade tolerance
- Minimum DBH parameter in stand metrics

## medfate 2.7.0

- Basic water balance model with relative water content
- ‘pRootDisc’ eliminated from species parameters
- Plant water balance and cuticular transpiration added to the basic
  water balance model
- Water pools revised for the basic and water balance models
- Shared water pools is now controlled via parameter
  ‘rhizosphereOverlap’
- New species-specific parameter ‘MortalityBaselineRate’
- Live fuel moisture content now included in simulation results
- New species-specific parameter ‘RSSG’ (minimum relative starch for
  sapwood growth)

## medfate 2.6.2

- Control parameter ‘modifyInput’ is no longer available. Functions
  ‘spwb’, ‘pwb’ and ‘growth’ do not modify input objects and return an
  element ‘spwbOutput’ or ‘growthOutput’ with a copy of the final state.

## medfate 2.6.1

- Tree cover (open grown assumption)
- Parameters ‘ShrubCover’ and ‘CanopyCover’ no longer required in
  function fuel_FCCS
- Input weather stored in output from functions ‘spwb’, ‘pwb’ and
  ‘growth’
- New function ‘fireHazard’.
- Relative bias and relative MAE as evaluation metrics
- Drought-related leaf senescence only occurring if ‘StemPLC’ increases
- Embolized sapwood proportion discounted from sapwood maintenance
  respiration
- Evaluation of diameter increment (DI), DBH and Height series
- Output growth rates not relative to sapwood area

## medfate 2.6.0

- Nitrogen content for leaves, sapwood and fine roots added. ‘Nleaf’
  replaces ‘Narea’ as the latter can be calculated from ‘Nleaf’ using
  ‘SLA’.
- Maintenance respiration rates based on N concentration of tissues

## medfate 2.5.0

- spwb model with Granier transpiration now extracts water from soil
  layer according to unsaturated conductivity.
- shinyplot generic function.
- Update parameters fHDmin and fHDmax.
- New vignette ‘IFNEvaluation’
- New parameter ‘WUE_decay’ for reduction of relative WUE in Granier’s
  model
- Tissue construction costs are now species-specific parameters.
- Fine root growth and senescence made equal between Granier and Sperry
  models, inducing a new parameter ‘Ar2Al’ for Granier’s model.

## medfate 2.4.0

- Functions ‘spwb_stress’, ‘spwb_resistances’ and
  ‘spwb_waterUseEfficiency’ renamed to ‘droughtStress’, ‘resistances’
  and ‘waterUseEfficiency’, respectively, since they can now be applied
  to the output of several simulation function.
- Plant biomass balance in growth. Structural changes are now daily in
  growth simulations.
- Bug correction in shrub structural update. Shrub dynamics default set
  to TRUE.
- DBH/Height plots from growth output.
- IMPORTANT: New species parameter table.

## medfate 2.3.8

- Summary functions revised, including new function ‘summary.fordyn’
- Collating intra-annual ‘fordyn’ results for plotting.
- Forest management enabled in function ‘fordyn’ and default management
  actions defined in ‘defaultManagementFunction’ and
  ‘defaultManagementArguments’
- DOY, Photoperiod and JulianDay can be taken from weather input in
  functions ‘spwb’, ‘pwb’ and ‘growth’.

## medfate 2.3.7

CRAN release: 2021-12-16

- CO2 made an daily input weather variable, in addition to the default
  control parameter ‘Catm’
- New function ‘shinyplot’ to create interactive graphics
- New function ‘plot.fordyn’ to display annual (step) summaries of
  forest dynamics

## medfate 2.3.5

- New example vignette ‘FontBlanche’
- Modification of evaluation functions to separate the evaluation of
  total evapotranspiration (‘ETR’) from the evaluation of soil
  evaporation + plant transpiration (‘SE+TR’)
- Bug correction: sub-daily stomatal conductance plots

## medfate 2.3.4

- Maximum relative sapwood growth rates effective
- New option ‘summary.freq’ in plot.spwb and plot.growth
- New species parameters (sapwood and fine root senescence rates)
- Defaults for ‘conduit2sapwood’ from taxonomical family
- IMPORTANT: New species parameter tables (including estimates for
  conduit2sapwood)

## medfate 2.3.2

- Parameter ‘ParticleDensity’ eliminated, as it is now calculated from
  ‘LeafDensity’, ‘WoodDensity’ and ‘r635’

## medfate 2.3.1

- IMPORTANT: New species parameter tables
- New control parameter ‘fillMissingSpParams’
- Defaults for ‘LeafDensity’, ‘WoodDensity’, ‘LeafPI0’, ‘LeafEPS’ and
  ‘LeafAF’ from taxonomical family.
- New species categorical params (for inbuilt imputation) ‘LeafShape’
  and ‘LeafSize’
- Parameter ‘Flammability’ index eliminated (non-meansurable property of
  flammability).
- Parameter ‘LeafLitterFuelType’ eliminated because it is derived from
  ‘LeafShape’ and ‘LeafSize’.
- Defaults added for ‘r635’, ‘heatContent’, ‘LigninPercent’ and ‘SAV’
  according to ‘LeafShape’ and ‘LeafSize’
- New function ‘getSpParamsDefinition()’ returns definition of species
  parameters.
- Inbuilt defaults added for shrub allometries, depending on ‘LifeForm’
  and ‘Hmax’, and for tree allometries, depending on ‘Group’.
- Default value added for ‘pDead’.

## medfate 2.2.3

CRAN release: 2021-06-18

- Functions soilgridsParams() and forest_map\*Tables() moved to package
  ‘medfateutils’ available at GitHub (emf-creaf/medfateutils).

## medfate 2.2.2

- Dependency ‘spdep’ removed
- Bug correction ‘windKatul.cpp’

## medfate 2.2.1

CRAN release: 2021-06-11

- Calibrated minimum bioclimatic parameters for recruitment
  (SpParamsMED)
- Explicit species input parameters for phenology

## medfate 2.2.0

- New simulation function ‘fordyn()’, including recruitment process
- New function ‘mergeShrubs’
- New functions ‘species_parameter’ and ‘species_characterParameter’

## medfate 2.1.4

- Dessication/defoliation homogenized across transpiration modes in
  growth function
- Revision of mortality (stochastic/deterministic, whole-cohort/density)

## medfate 2.1.3

- Update Psi_Extract according to Psi_TLP (Hydratry)

## medfate 2.1.2

- New control flags for defoliation/starvation/dessication in growth
  simulations
- New control flag for sink limitation in growth simulations
- Revision of phenology submodel

## medfate 2.1.1

- Respiration rate for leaves made optionally species-specific.
- Update shrub allometries from De Caceres et al. (2019).
- LAI_live and LAI_expanded in growth.
- Update SpParamsUS (missing values for new parameters).

## medfate 2.1.0

- Parameters of Granier’s equation made species-specific if available
- Parameter modification dependencies revised
- Modification of input objects is now optional
- Update of function `transp_maximumTranspirationModel`
- Photosynthesis in Granier’s model corresponds to gross photosynthesis
  and is proportional to transpiration

## medfate 2.0.1

- Cloning initial object for optimization
- Leaf growth costs always drawn from sugar sapwood
- Bug correction: LWRnet calculation
- Bug correction: Sperry model does not crash when LAIstand = 0

## medfate 2.0.0

- IMPORTANT: Soil input merged with model input. Now the calls to
  simulation functions (e.g., `spwb`, `growth`) do not need to include
  soil as input parameter.
- Bug correction: NaN values for theta \> theta_sat in van Genuchten psi
  computation
- Function `modifyInputParams` now accepts modification of soil layer
  properties
- Sensitivity/Calibration vignette updated

## medfate 1.1.6

- New canopy turbulence models by Katul et al (2004).
- Long-wave net radiation balance for layered canopies following
  Flerchinger et al. (2009).
- Multi-layer canopy energy balance as in Bonan et al. (2014)
- New option ‘depthMode’ to calculate fuel depth in ‘fuel_FCCS’.
- Modification of ‘vprofile_windExtinction’ to draw turbulence models.

## medfate 1.1.5

- Stomatal conductance now denoted as Gsw or GSW
- Boundary layer conductance considered
- Leaf water potential influencing leaf vapour pressure

## medfate 1.1.4

- Flexible temporal resolution of model evaluation
- Basal area index evaluation

## medfate 1.1.3

- New model evaluation functions
- New optimization function factories
- New function ‘modifyInputParams’
- Control parameters set to a nested list
- New vignette for calibration and sensitivity analysis

## medfate 1.1.2

- GW sunlit/shade minimum/maximum daily output
- FMC calculations with basic water balance output
- New output in growth simulation (biomass values)
- New function ‘moisture_cohortFMCDay’
- Modification of Z50/Z95

## medfate 1.1.1

- Root exudation added to carbon balance
- Revision sapwood growth
- Growth cost for fine roots in basic model
- Translocation for carbon during senescence
- Bug correction in fuel calculations with US mode

## medfate 1.1.0

CRAN release: 2020-11-05

- Control option ‘rockyLayerDrainage’ instead of ‘drainage’ to disable
  macropore vertical outflow in layers with \> 95% of rocks
- Soil parameter Kdrain for saturated vertical hydraulic conductivity
  towards groundwaters (deep drainage)
- Improved validation plots with confidence intervals

## medfate 1.0.3

- Nash-Sutcliffe efficiency (NSE) implemented in spwb_validation

## medfate 1.0.2

- Bug correction on the use of organic matter in Saxton (2006) equations
  (thanks to Milan Fischer).
- Recodification of soilgridParams due to new SoilGrids REST API
  (removed dependency from GSIF).
- New function ‘redefineSoilLayers’

## medfate 1.0.1

- New root functions
- Bulk density stored in soil object initialization
- Advanced plant water pools
- Belowground inputs restructured

## medfate 1.0.0

CRAN release: 2020-05-17

- Reorganization of growth function
- Clarification of gross and net photosynthesis
- ‘spwb_resetInputs’ to ‘resetInputs’
- Dependence of kmax on temperature (due to sap dynamic viscosity)
  incorporated
- Functions plot.spwb and plot.pwb modified to draw subdaily dynamics
  for a subset of dates

## medfate 0.9.1

- Small bug fixes
- update ‘spwb_resetInputs’
- Update of plant water pools
- New output (annual stand summaries and aboveground structure) for
  function ‘growth’

## medfate 0.9.0

CRAN release: 2020-03-23

- New parametrization data set ‘SpParamsUS’
- Function ‘hydrology_verticalInputs’ replaced by
  ‘hydrology_soilWaterInputs’ and
  ‘hydrology_soilInfiltrationPercolation’.
- New simulation control option: ‘plantWaterPools’.
- Hard (Imports) dependency from GSIF changed to soft (Suggests) one.

## medfate 0.8.9

- Update supply function plot.

## medfate 0.8.8

- New function ‘soil_rockWeight2Volume’

## medfate 0.8.7

- Corrections to energy balance for zero LAI (deciduous species)
- SFI functions moved to medfateland
- soilgridParams modified to accept a SpatialPoints object as input

## medfate 0.8.6

- New function ‘spwb_sensitivity’ for sensitivity analyses
- New control parameter ‘unlimitedSoilWater’
- Bug correction in canopy height with LAI = 0
- Modifications of spwb_ldrOptimization to work with transpirationMode =
  “Sperry”
- New function spwb_ldrExploration

## medfate 0.8.5

- New control parameter ‘fracLeafResistance’
- Different control options for parameter ‘cavitationRefill’
- New control parameter ‘cavitationRecoveryMaximumRate’
- Control parameter ‘hydraulicCostFunction’ replaced by ‘costModifier’
  and ‘gainModifier’
- New control parameter ‘cuticularTranspiration’
- Numerical controls to avoid NaN in functions ‘soil_theta2psiSX’ and
  ‘soil_psi2thetaSX’
- Bug correction in estimation of root conductance proportions

## medfate 0.8.4

- Water balance console output modified in spwb
- New approach to plant water compartments ‘capacitance = TRUE’
- Output of plant water balance
- New option in control
- Stem segments fixed to two
- functional parameter pRootDisc removed from Sperry’s advanced model
- Bug correction in fuel_cohortFineFMC
- Remove ksympver and add klatleaf/klatstem to control parameters

## medfate 0.8.3

- Reference book (medfatebook) launched
- Fraction of absorbed SWR output in Granier’s transpiration
- ‘Stand’ data frame output in spwb(), separated from ‘WaterBalance’
- New function ‘forest_mapTreeTable’, ‘forest_mapShrubTable’ and
  ‘forest_mapWoodyTables’
- ’stand\_\*’ functions for stand-level properties
- No SWR soil absorption when snow pack is present in Sperry’s model
- verticalLayerSize made a control parameter for Granier’s model
- Changes in light parameters: New parameter ‘alphaSWR’. ‘albedo’
  renamed to ‘gammaSWR’. ‘k’ renamed to ‘kPAR’

## medfate 0.8.2

CRAN release: 2019-05-29

- Modification of meteoland to better calculate direct/diffuse light on
  slopes
- Added Narea parameter to facilitate estimation of Vmax298
- New function ‘spwb_validation’.
- Plot functions using ggplot.
- Improvement of infiltration repartition for varying macroporosity.
- Input values for latitude and topography stored in the result of
  simulations.
- Wind for each cohort stored in the result of
  transp_transpirationSperry and spwb_day.

## medfate 0.8.1

- Export Ci from spwb_day.
- New function ‘maximumTranspirationRatioPlot’
- Output of min/max water potential for sunlit and shade leaves.
- New function ‘soil_waterRetentionPlot’
- New function ‘waterUseEfficiency’
- Default Van Genuchten PTF set to ‘Toth’
- New functions ‘modifySpParams’ and ‘modifyCohortParams’.
- New function ‘hydrology_interceptionPlot’

## medfate 0.8.0

- Move spatial classes and methods to package ‘medfateland’
- Implement underscores instead of dots to separate function groups and
  function names
- New function transp_Granier.
- Bug corrected in spwb.plot for snow plotting.
- New function ‘pwb’.
- New function ‘snowMelt’.
- Transpiration model changed to Granier and Sperry.
- PLC set to zero when DOY = 1
- New functions for leaf phenology

## medfate 0.7.4

CRAN release: 2019-03-22

- Fixing bugs for installation in all platforms
- Update docs

## medfate 0.7.3

CRAN release: 2019-03-19

- Shrub root system using LDR (Z50 and Z95)
- Percolation of infiltrated water consistent with layer subdivision
- New function soil.infiltrationRepartition()
- Adding interception to evapotranspiration
- Checked for CRAN

## medfate 0.7.2

- Accounting of hydraulic redistribution
- Update plotting functions
- New function vprofile.RootDistribution.
- New functions for water at wilting point (-1.5 MPa)
- Print extractable water of soil
- Bug correction in photosynthesis (now done per leaf area basis)
- Leaf area distribution (and crown fine biomass distribution) following
  truncated normal \[-1.5,1.5\]

## medfate 0.7.1

- Interception corrected in the complex model
- Etol set to 0.0000001
- spwb export of dEdP (equivalent to soil-plant conductance)
- Subdaily results can be stored for spwb
- New function spwb.resistances to calculate and draw segment
  resistances for spwb simulation results.
- Cohort parameter search by SpIndex
- Gwmin set to zero when capacitance = FALSE
- New function spwb.stress to calculate drought stress indices
- PlantStress in complex mode now is calculated as relative soil-plant
  conductance, for compatibility with the simple mode.
- Function name changes for interception and soil hydrology

## medfate 0.7.0

- Leaf and stem water compartments added
- Update of functions ‘spwb.day’ and ‘plot.spwb.day’
- Analytical integral of the van Genuchten function
- Analytical inverse of the incomplete gamma function

## medfate 0.6.2

- Adapt to Rcpp changes
- Stem water compartments

## medfate 0.6.1

- Saturated theta in Saxton model
- Water table depth added
- Boolean option ‘drainage’ added to soil water balance
- Added subsurface flow processes to spwbgrid

## medfate 0.6.0

- Function name ‘swb’ (and all related functions) changed to ‘spwb’
- Added snow pack to soil state variables.
- Added new control option ‘snowpack’ for snowpack dynamics simulation
  (only when transpirationMode=‘Simple’)
- Changed names of spatial classes from ‘Forest’ to ‘Landscape’
- Changed name of ‘exampleSPF’ to ‘exampleSPL’.
- Function ‘spwbgrid’ functional again.
- Improved documentation of function ‘spwb.day’
- ‘DOY’ no longer needed as column in meteorological input.

## medfate 0.5.6

- Growth degree days added as output of ‘swb()’
- New function ‘swb.resetInputs()’
- New function ‘hydraulics.vulnerabilityCurvePlot()’
- Reorganization of help for hydraulics.
- Documentation of tissue moisture functions.
- Update function ‘plot.growth()’

## medfate 0.5.5

- New functions for tissue moisture
- New function ‘fuel.cohortFineFMC’
- New option ‘bySpecies’ to aggregate results by species in functions
  ‘summary.swb’ and ‘summary.growth’
- New option ‘bySpecies’ to aggregate results by species in functions
  ‘plot.swb’ and ‘plot.swb.day’
- Bulk density added to soil parameters
- New set of pedotransfer functions to calculate VG parameters from
  texture, bulk density and organic content

## medfate 0.5.4

- New functions for pressure-volume curves
- Rmarkdown vignettes
- Added Van Genuchten pedotransfer functions
- New control parameter added (for soil functions)
- New function ‘forest2belowground’
- Update ‘summary.swb’ function

## medfate 0.5.3

- Network representation of the continuum now includes a leaf segment
- Stem fraction of total resistance removed as estimation of root xylem
  conductance
- Default kleaf_max (8 for temperate angiosperms and 5 for gymnosperms)
- New parameter rootxylem_kmax (hydraulic conductivity of roots)

## medfate 0.5.2

- Modification of default for ntrial
- Three element supply function added
- Parameters of the leaf maximum conductance and leaf vulnerability
  curve added

## medfate 0.5.1

- Structure of swb.day output
- Increased output (leaf temperature, stomatal conductance and leaf VPD)
  in swb.day
- Taper modifications
- New function plot.swb.day

## medfate 0.5.0

- Atmospheric CO2 is an input control parameter
- Added Hmed as species parameter (to correct reference conductivity
  values)
- Modified documentation

## medfate 0.4.9

- Alternative way of calculating maximum root conductance
- Default value for averageFracRhizosphereResistance changed to 0.15
- Root vulnerability curve parameters taken from stem vulnerability
  curve when missing (d_root = d_stem/2)

## medfate 0.4.8

- Leaf radiation balance with LWR from soil
- New output for swb
- Bug fix in plot.swb (match according to transpiration model)

## medfate 0.4.7

- Profit maximization for sunlit and shade leaves separately
- Bug correction in swb
- Check on stomatal conductances in profit maximization

## medfate 0.4.6

- Radiation absorbed by trunks
- Separation of energy balance components

## medfate 0.4.5

- New control parameter: Canopy thermal capacity per LAI
- Update of calls to meteoland (diffuse radiation)
- Wind value when missing data added to control

## medfate 0.4.4

- Minimum windspeed (1.0 m/s)
- Progressive leaf fall
- Bug correction (swb.plot)
- Energy balance output/plot
- Modified conductance scaling

## medfate 0.4.3

- Added new parameter ‘LeafWidth’

## medfate 0.4.2

- Design changes for radiation balance (soil/canopy/atmosphere)
- Added temperature balance in summary.swb and plot.swb
- Added new parameter ‘albedo’

## medfate 0.4.1

- Documentation of soil thermodynamics
- Energy balance output in swb.day

## medfate 0.4.0

- Added soil temperature state variables (for ‘complex’ mode)
- Added canopy temperature state variable (for ‘complex’ mode)
- ‘gdd’ in swbInput used as initial growth degree days
- Soil thermodynamics.

## medfate 0.3.7

- Added vertical layer size as option in control.
- New function ‘soilgridsParams’ (calls functions in package ‘GSIF’)
- Added longwave radiation to leaf energy balance.

## medfate 0.3.6

- Added new pedotransfer functions with organic matter.
- New function added ‘soil.waterFC’.
- Xylem taper added as option in control.
- Hydraulic cost function added as option in control.
- New function ‘light.instantaneousLightExtinctionAbsortion’.
- New functions ‘transp.dayCanopyTranspiration’ and
  ‘transp.dayCanopyTranspirationPlot’.

## medfate 0.3.5

- Control of numerical methods for supply function added to
  ‘defaulControlParams’.
- Added new function ‘swb.ldrCalibration’ to calibrate root distribution
  for swb simulations (by Victor Granda)
- Added documentation for light extinction functions
  ‘light.layerIrradianceFraction’, ‘light.layerSunlitFraction’ and
  ‘light.cohortSunlitShadeAbsorbedRadiation’.
