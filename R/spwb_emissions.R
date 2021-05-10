spwb_emissions<-function(x, verbose = FALSE) {
  data("MEGANParams")
  
  if((!x$spwbInput$control$subdailyResults) || ((x$spwbInput$control$transpirationMode) != "Sperry") || ((x$spwbInput$control$ndailysteps) != 24)) {
    stop("The emissions of gases and aerosols of MEGAN2.1 (see Guenther et al. 2012) requires subdailyResults and [Sperry] transpiration model and hourly step. Please rerun")
  } 
  NumberOfDays = length(x[["subdaily"]])
  NumberOfCohorts = nrow(x$spwbInput$cohorts) 
  NumberOfDailySteps = x$spwbInput$control$ndailysteps
  
  CompoundEmission = matrix(NA, NumberOfDays, nrow(MEGANParams))
  colnames(CompoundEmission)<- MEGANParams$Compound.Class
  rownames(CompoundEmission)<- names(x[["subdaily"]])
  cat(paste0("Processing daily emissions.\n"))
  if(!verbose) {
    pb = txtProgressBar(0, NumberOfDays, style = 3)
  }
  for(i in 1:NumberOfDays){
    if(verbose) cat(paste0("\nDay: ",i, "\n"))
    else setTxtProgressBar(pb, i)
    # # LAI values can change depending on leaf phenology
    # SunlitCohortLAI = x[["subdaily"]][[i]][["SunlitLeaves"]]$LAI
    # ShadeCohortLAI = x[["subdaily"]][[i]][["ShadeLeaves"]]$LAI
    # SunlitAndShadeCohortTotalLAI = SunlitCohortLAI + ShadeCohortLAI
    for(mm in 1:nrow(MEGANParams)){
      if(verbose) cat(paste("Compound.Class=",MEGANParams[mm,1],"and its LDF=",MEGANParams[mm,]$LDF, "\n", sep=' '))
      # CompoundEmission <- SunlitCohortLAI  #Borrow the structure
      if(mm ==1 ){
        #Megan model light response emission following the Megan paper for sunlit leaves only---start
        PPFD = x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Abs_PAR * 4.55   #Megan says it is absorbed shortwave in W.m-2, while MEGAN required ??mole.m2/s, see MEGAN paper 4.55 and also https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
        LeafTemperature = x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Temp + 273.15   #Megan says is is ShadeLeavesInst, which may not be this Abs_PAR. Will check
        PCurrentDay <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        PPreviousDay <- x[["subdaily"]][[max(1,i-1)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55  #max() is to ensure no zero
        P24 <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Abs_PAR * 4.55   #borrow the structure with value being meanglessful
        TCurrentDay <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Temp + 273.15
        TPreviousDay <- x[["subdaily"]][[max(1,i-1)]][["SunlitLeavesInst"]]$Temp + 273.15  #max() is to ensure no zero
        T24 <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Temp + 273.15   #borrow the structure with value being meanglessful    
        for(k in 1:24){
          if(k == 1){
            P24[,k] = rowMeans(PPreviousDay[,k:24])
            T24[,k] = rowMeans(TPreviousDay[,k:24])
          }else{
            P24[,k] = rowMeans(cbind2(PPreviousDay[,k:24],PCurrentDay[,1:k-1]))   #check rowMeans function
            T24[,k] = rowMeans(cbind2(TPreviousDay[,k:24],TCurrentDay[,1:k-1]))   #check rowMeans function
          }
        }
        P240 <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Abs_PAR * 4.55   #borrow the structure with value being meanglessful
        P240Negative10Day <- x[["subdaily"]][[max(1,i-10)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative9Day <- x[["subdaily"]][[max(1,i-9)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative8Day <- x[["subdaily"]][[max(1,i-8)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative7Day <- x[["subdaily"]][[max(1,i-7)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative6Day <- x[["subdaily"]][[max(1,i-6)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative5Day <- x[["subdaily"]][[max(1,i-5)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative4Day <- x[["subdaily"]][[max(1,i-4)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative3Day <- x[["subdaily"]][[max(1,i-3)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative2Day <- x[["subdaily"]][[max(1,i-2)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240Negative1Day <- x[["subdaily"]][[max(1,i-1)]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        P240CurrentDay <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Abs_PAR * 4.55
        T240 <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Temp + 273.15   #borrow the structure with value being meanglessful
        T240Negative10Day <- x[["subdaily"]][[max(1,i-10)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative9Day <- x[["subdaily"]][[max(1,i-9)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative8Day <- x[["subdaily"]][[max(1,i-8)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative7Day <- x[["subdaily"]][[max(1,i-7)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative6Day <- x[["subdaily"]][[max(1,i-6)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative5Day <- x[["subdaily"]][[max(1,i-5)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative4Day <- x[["subdaily"]][[max(1,i-4)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative3Day <- x[["subdaily"]][[max(1,i-3)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative2Day <- x[["subdaily"]][[max(1,i-2)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240Negative1Day <- x[["subdaily"]][[max(1,i-1)]][["SunlitLeavesInst"]]$Temp + 273.15
        T240CurrentDay <- x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Temp + 273.15    
        for(k in 1:24){
          First = P240Negative10Day[,k:24]
          Last = P240CurrentDay[,1:k-1]   #check rowMeans function
          P240[,k] = rowMeans(cbind2(First,P240Negative9Day,P240Negative8Day,P240Negative7Day,P240Negative6Day,P240Negative5Day,P240Negative4Day,P240Negative3Day,P240Negative2Day,P240Negative1Day,Last))
          FirstT = T240Negative10Day[,k:24]
          LastT = T240CurrentDay[,1:k-1]   #check rowMeans function
          T240[,k] = rowMeans(cbind2(FirstT,T240Negative9Day,T240Negative8Day,T240Negative7Day,T240Negative6Day,T240Negative5Day,T240Negative4Day,T240Negative3Day,T240Negative2Day,T240Negative1Day,LastT))
        }
        for(j in 1:NumberOfCohorts){
          Ps = 200 #200 is for sunlit leaves. It is 50.0 for shade leaves
          Alpha = 0.004 - 0.0005 * log(P240)
          Cp = 0.0468 * exp(0.0005 * (P24 - Ps)) * (P240 ** 0.6)   
          Yp_ldf = Cp * ((Alpha * PPFD) / (((1 + Alpha * Alpha * PPFD * PPFD) ** 0.5)))   #Yp_ldf = Cp * ((Alpha * PPFD) / (power((1 + Alpha * Alpha * PPFD * PPFD),0.5)))
        }
        #Megan model light response emission following the Megan paper for sunlit leaves only---end
        
        
        #Megan model light response emission following the Megan paper for shaded leaves only---start
        ShadePPFD = x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Abs_PAR * 4.55   #Megan says is is PPFD, which may not be this Abs_PAR. Will check
        ShadeLeafTemperature = x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Temp + 273.15   #Megan says is is PPFD, which may not be this Abs_PAR. Will check
        ShadePCurrentDay <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadePPreviousDay <- x[["subdaily"]][[max(1,i-1)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55  #max() is to ensure no zero
        ShadeP24 <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Abs_PAR * 4.55   #borrow the structure with value being meanglessful
        ShadeTCurrentDay <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeTPreviousDay <- x[["subdaily"]][[max(1,i-1)]][["ShadeLeavesInst"]]$Temp + 273.15  #max() is to ensure no zero
        ShadeT24 <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Temp + 273.15   #borrow the structure with value being meanglessful    
        for(k in 1:24){
          if(k == 1){
            ShadeP24[,k] = rowMeans(ShadePPreviousDay[,k:24])
            ShadeT24[,k] = rowMeans(ShadeTPreviousDay[,k:24])
          }else{
            ShadeP24[,k] = rowMeans(cbind2(ShadePPreviousDay[,k:24],ShadePCurrentDay[,1:k-1]))   #check rowMeans function
            ShadeT24[,k] = rowMeans(cbind2(ShadeTPreviousDay[,k:24],ShadeTCurrentDay[,1:k-1]))   #check rowMeans function
          }
        }
        ShadeP240 <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Abs_PAR * 4.55   #borrow the structure with value being meanglessful
        ShadeP240Negative10Day <- x[["subdaily"]][[max(1,i-10)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative9Day <- x[["subdaily"]][[max(1,i-9)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative8Day <- x[["subdaily"]][[max(1,i-8)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative7Day <- x[["subdaily"]][[max(1,i-7)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative6Day <- x[["subdaily"]][[max(1,i-6)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative5Day <- x[["subdaily"]][[max(1,i-5)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative4Day <- x[["subdaily"]][[max(1,i-4)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative3Day <- x[["subdaily"]][[max(1,i-3)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative2Day <- x[["subdaily"]][[max(1,i-2)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240Negative1Day <- x[["subdaily"]][[max(1,i-1)]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeP240CurrentDay <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Abs_PAR * 4.55
        ShadeT240 <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Temp + 273.15   #borrow the structure with value being meanglessful
        ShadeT240Negative10Day <- x[["subdaily"]][[max(1,i-10)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative9Day <- x[["subdaily"]][[max(1,i-9)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative8Day <- x[["subdaily"]][[max(1,i-8)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative7Day <- x[["subdaily"]][[max(1,i-7)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative6Day <- x[["subdaily"]][[max(1,i-6)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative5Day <- x[["subdaily"]][[max(1,i-5)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative4Day <- x[["subdaily"]][[max(1,i-4)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative3Day <- x[["subdaily"]][[max(1,i-3)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative2Day <- x[["subdaily"]][[max(1,i-2)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240Negative1Day <- x[["subdaily"]][[max(1,i-1)]][["ShadeLeavesInst"]]$Temp + 273.15
        ShadeT240CurrentDay <- x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Temp + 273.15    
        for(k in 1:24){
          ShadeFirst = ShadeP240Negative10Day[,k:24]
          ShadeLast = ShadeP240CurrentDay[,1:k-1]   #check rowMeans function
          ShadeP240[,k] = rowMeans(cbind2(ShadeFirst,ShadeP240Negative9Day,ShadeP240Negative8Day,ShadeP240Negative7Day,ShadeP240Negative6Day,ShadeP240Negative5Day,ShadeP240Negative4Day,ShadeP240Negative3Day,ShadeP240Negative2Day,ShadeP240Negative1Day,ShadeLast))
          ShadeFirstT = ShadeT240Negative10Day[,k:24]
          ShadeLastT = ShadeT240CurrentDay[,1:k-1]   #check rowMeans function
          ShadeT240[,k] = rowMeans(cbind2(ShadeFirstT,ShadeT240Negative9Day,ShadeT240Negative8Day,ShadeT240Negative7Day,ShadeT240Negative6Day,ShadeT240Negative5Day,ShadeT240Negative4Day,ShadeT240Negative3Day,ShadeT240Negative2Day,ShadeT240Negative1Day,ShadeLastT))
        }
        for(j in 1:NumberOfCohorts){
          Ps = 50.0 #50 is for shade leaves. It is 200 for sunlit leaves
          Alpha = 0.004 - 0.0005 * log(ShadeP240)
          Cp = 0.0468 * exp(0.0005 * (ShadeP24 - Ps)) * (ShadeP240 ** 0.6)   #continue here tomorrow
          Yp_lif = Cp * ((Alpha * ShadePPFD) / (((1 + Alpha * Alpha * ShadePPFD * ShadePPFD) ** 0.5)))   #Yp_ldf = Cp * ((Alpha * PPFD) / (power((1 + Alpha * Alpha * PPFD * PPFD),0.5)))
        }
        #Megan model light response emission following the Megan paper for shaded leaves only---end
      }      
      #Megan model temperature response emission following the Megan paper for sunlit leaves only---start
      for(j in 1:NumberOfCohorts){            
        Ts = 297
        Topt = 313 + (0.6 * (T240 - Ts))
        Ceo = MEGANParams$Ceo[mm]  #Check Megan paper Table 4 with this coefficient ranging from 1.6 to 2.37 for different emission class
        Eopt = Ceo * exp(0.05 * (T24 - Ts)) * exp(0.05 * (T240 - Ts))
        xvalue = ((1/Topt)-(1/LeafTemperature)) / 0.00831   #LeafTemperature refers to T following Megan equation 8 
        Ct2 = 230
        Ct1 = MEGANParams$Ct1[mm] #Check Megan paper Table 4 with this coefficient ranging from 60 to 130 for different emission class
        Yt_ldf = Eopt * (Ct2 * exp(Ct1 * xvalue) / (Ct2 - Ct1 * (1 - exp(Ct2 * xvalue))))
      }
      #Megan model temperature response emission following the Megan paper for sunlit leaves only---end 
      
      #Megan model temperature response emission following the Megan paper for shade leaves only---start
      for(j in 1:NumberOfCohorts){            
        Ts = 297
        Belta = MEGANParams$belta[mm]  #Check Megan paper Table 4 with this coefficient ranging from 0.08 to 0.17 for different emission class
        Yt_lif = exp(Belta *(ShadeLeafTemperature - Ts))   #note this is for ShadeLeaves, #LeafTemperature refers to T in Megan equation 11
      }
      #Megan model temperature response emission following the Megan paper for shade leaves only---end
      
      
      #This is Megan 2016 paper equation 3 (note Alex's original equation is wrong, see his email dated Mon 4/12/2021 10:29 AM) by combining Yp_lif and Yp_ldf---start
      LDF = MEGANParams$LDF[mm]   #The LDF is listed in Megan Table 4 for each compound. We will change it later
      Ypi = (1 - LDF) * Yp_lif + LDF * Yp_ldf   #The is for Megan equation 3, which is actually wrong/incomplete according to Alex email
      #This is Megan 2016 paper equation 3 (note Alex's original equation is wrong, see his email dated Mon 4/12/2021 10:29 AM) by combining Yp_lif and Yp_ldf---end
      
      #This is Megan 2016 paper equation 7 by combining Yt_lif and Yt_ldf---start
      Yti = (1 - LDF) * Yt_lif + LDF * Yt_ldf  
      #This is Megan 2016 paper equation 7 by combining Yt_lif and Yt_ldf---end
      
      
      
      
      #https://stackoverflow.com/questions/12042309/compare-two-matrices-with-criteria-in-r
      #so we will have four matrix of Fmat, Fnew, Fgro, and Fold corresponding to each cohort
      Current_DailyAtmosphereTemperature = x$Temperature$Tatm_mean[i] + 273.15
      Previous_DailyAtmosphereTemperature = x$Temperature$Tatm_mean[max(1,i-1)] + 273.15
      Current_DailyCohortLAI = x[["subdaily"]][[i]][["Plants"]]$LAI
      Previous_DailyCohortLAI = x[["subdaily"]][[max(i-1,1)]][["Plants"]]$LAI
      Fmat = Current_DailyCohortLAI   #borrow the structure, the value itself is not useful
      Fnew = Current_DailyCohortLAI   #borrow the structure, the value itself is not useful
      Fgro = Current_DailyCohortLAI   #borrow the structure, the value itself is not useful
      Fold = Current_DailyCohortLAI   #borrow the structure, the value itself is not useful
      for(p in 1:length(Current_DailyCohortLAI)){
        if(Current_DailyCohortLAI[p] == Previous_DailyCohortLAI[p]){
          Fmat[p] = 0.8
          Fnew[p] = 0.0
          Fgro[p] = 0.1
          Fold[p] = 0.1
        }else if(Current_DailyCohortLAI[p] < Previous_DailyCohortLAI[p]){
          Fnew[p] = 0.0
          Fgro[p] = 0.0
          Fold[p] = (Previous_DailyCohortLAI[p] - Current_DailyCohortLAI[p]) / Previous_DailyCohortLAI[p]    
          Fmat[p] = 1 - Fold[p]
        }else if(Current_DailyCohortLAI[p] > Previous_DailyCohortLAI[p]){
          Fold[p] = 0.0
          t = 1
          if(Previous_DailyAtmosphereTemperature <= 303){
            ti = 5 + (0.7 * (300 - Previous_DailyAtmosphereTemperature))
          }else{
            ti = 2.9
          }
          tm = 2.3 * ti
          if(t <= ti){Fnew[p] = 1 - (Previous_DailyCohortLAI[p] / Current_DailyCohortLAI[p])}
          if(t > ti){Fnew[p] = (ti / t ) * (1 - (Previous_DailyCohortLAI[p] / Current_DailyCohortLAI[p]))}
          if(t <= tm){Fmat[p] = Previous_DailyCohortLAI[p] / Current_DailyCohortLAI[p]}
          if(t > tm){Fmat[p] = (Previous_DailyCohortLAI[p] / Current_DailyCohortLAI[p]) + ((t - tm) / t) * (1 - Previous_DailyCohortLAI[p] / Current_DailyCohortLAI[p])}
          Fgro[p] = 1 - Fnew[p] - Fmat[p]
        } 
      }
      Yai = Fnew * MEGANParams$Anew[mm] + Fgro * MEGANParams$Agro[mm] + Fmat * MEGANParams$Amat[mm] + Fold * MEGANParams$Aold[mm]
      
      
      
      #This is Megan 2016 paper equation 13 for soil moisture---start
      if(MEGANParams$Compound.Class[mm] == "Isoprene"){
        THETAw = sum(soil_waterWP(x$spwbInput$soil, model = "SX"))  #soil_waterWP calculate the water volume (in mm of soil volume) of each soil layer at wilting point (-1.5 MPa),
        THETA1 = THETAw + 40   #(note 40 is defined in Megan as 0.04 m, so I converted it to 40 mm)
        THETA = x$Soil$MLTot  #MLTot": Total soil water volume (in L/m2) in Medfate. We can be converted to mm (1 mm = 1 l / m?), see https://www.researchgate.net/post/How-can-I-convert-Satellite-Soil-Moisture-data-m3-m3-into-mm
        Ysm = THETA #borrow the structure
        for(kk in 1:length(THETA)){   #I used kk instead i here because i is used for numberofay in the whole loop
          if(THETA[kk] > THETA1){Ysm[kk] = 1}
          else if ((THETA[kk] > THETAw) & (THETA[kk] < THETA1)){Ysm[kk] = (THETA[kk]-THETAw)/40.0}
          else if (THETA[kk] < THETAw){Ysm[kk] = 0}
        }
      }else{
        Ysm = 1.0
      }
      #This is Megan 2016 paper equation 13 for soil moisture---end
      
      
      
      #This is Megan 2016 paper equation 14 for internal CO2. Improvements: 1) No assumpation of 70%; 3) Separate Sunlit and Shade leaves---start
      if(MEGANParams$Compound.Class[mm] == "Isoprene"){
        CiSunlit = x[["subdaily"]][[i]][["SunlitLeavesInst"]]$Ci
        CiShade = x[["subdaily"]][[i]][["ShadeLeavesInst"]]$Ci
        Ismax = 1.344 #Check Heald et al. 2009 table 2
        h = 1.4614   #Check Heald et al. 2009 table 2
        Cstar = 585   #Check Heald et al. 2009 table 2
        Yco2Sunlit = Ismax - ((Ismax * CiSunlit ** h) / (Cstar ** h + CiSunlit ** h))  #In MEgan 2.1, Ci is assumed to be 70% of air CO2, but Medfate calculate it, which is a good choice
        Yco2Shade = Ismax - ((Ismax * CiShade ** h) / (Cstar ** h + CiShade ** h))  #In MEgan 2.1, Ci is assumed to be 70% of air CO2, but Medfate calculate it, which is a good choice
      }else{
        Yco2Sunlit = 1
        Yco2Shade = 1
      }
      #This is Megan 2016 paper equation 14 for internal CO2. Improvements: 1) No assumpation of 70%; 3) Separate Sunlit and Shade leaves---end
      
      Cce = 0.6 #0.3 and 0.57 was used in his paper, but I do not know here. Temporary value
      EF = 999.99  #We will have the EF of each vegetation cohort for each conmpound here. Tempary value here
      CompoundEmissionSunlit = Cce * EF * x[["subdaily"]][[i]][["SunlitLeaves"]]$LAI * Ypi * Yti * Ysm[i] * Yco2Sunlit
      CompoundEmissionShade = Cce * EF * x[["subdaily"]][[i]][["ShadeLeaves"]]$LAI * Ypi * Yti * Ysm[i] * Yco2Shade
      CompoundEmission[i,mm] = sum(CompoundEmissionSunlit + CompoundEmissionShade, na.rm=TRUE)
      # assign(paste("HourlyCompoundEmissionDay",mm,"Of", str_replace_all(string=MEGANParams[mm,1], pattern=" ", repl=""), sep = ""), CompoundEmission)  #https://stackoverflow.com/questions/45864305/how-can-i-include-variables-in-a-dataframe-name-in-r  and https://stackoverflow.com/questions/5992082/how-to-remove-all-whitespace-from-a-string
      #In the future, we can do more statistics such as daily sum
      
    }
    
  }
  cat(paste0("\nDone.\n"))
  return(CompoundEmission)
}