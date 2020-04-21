extractSubdaily<-function(x, output = "E", dates = NULL)  {
  output = match.arg(
    output, 
    c("E","Ag","An","dEdPinst","PsiRoot",
      "PsiStem","PsiLeaf","PLCstem","RWCstem","RWCleaf","PWB"))
  if(is.null(dates)) dates = as.Date(names(x$subdaily))
  numDates = length(dates)
  numSteps = x$spwbInput$control$ndailysteps
  numCohorts = nrow(x$spwbInput$above)
  m<-data.frame(matrix(nrow = numDates*numSteps, ncol = numCohorts+1))
  for(i in 1:numDates) {
    ori = x$subdaily[[as.character(dates[i])]]$PlantsInst[[output]]
    m[((i-1)*numSteps+1):(i*numSteps), 2:(numCohorts+1)] = t(ori) 
  }
  colnames(m) = c("datetime", row.names(x$spwbInput$above))
  h = 0 + (0:(numSteps-1))*(24/numSteps)
  minutes = 60*h%%1
  seconds = round(60*minutes%%1)
  minutes = floor(minutes)
  hours = floor(h)
  times = paste(hours,minutes,seconds, sep=":")
  m$datetime = as.character(as.POSIXct(paste(dates[gl(n=numDates, k=numSteps)], times)))
  return(m)
}