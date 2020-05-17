## Script to reload parameters from 'SpParams.xlsx'
SpParamsMED <-as.data.frame(readxl::read_xlsx("data-raw/SpParams.xlsx",
                                              sheet="SpParamsMED", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsMED, overwrite = T)
rm(SpParamsMED)
SpParamsUS <-as.data.frame(readxl::read_xlsx("data-raw/SpParams.xlsx",
                                              sheet="SpParamsUS", na = "NA"), stringsAsFactors=FALSE)
usethis::use_data(SpParamsUS, overwrite = T)
rm(SpParamsUS)
#Rebuild!!!