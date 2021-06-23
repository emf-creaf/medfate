getSpParamsNames<-function() {
  data("SpParamsDefinition", package = "medfate", envir = environment())
  return(SpParamsDefinition$ParameterName)
}