# function g.LIDS.analyse_ENMOsub
# copy of g.LIDS.analyse
# using ENMO-threshold (duration + intensity of movement) as basis for LIDS calculation
# instead of a pure binary classification (duration of movement)
# 
# lastest change by EW: 06/11/2018
g.LIDS.analyse_ENMOsub = function(acc = c(), ws3 = 5, fit.criterion.cosfit = 2, 
                                  LIDS_cosfit_periods = seq(30,180,by=5), nonstationary = FALSE) {
  # see documentation in .Rd files in the man folder, 
  # these are just notes on input / output while developing the code
  
  # Input:
  # - acc: vector of acceleration time series
  # - ws3: epoch length of acc (seconds)
  # - fit.criterion.cosfit: metric to select the best LIDS (1=cor,2=cor*range)
  # - LIDS_cosfit_periods: vector with period lengths (minutes) to consider for the cosine fitting.
  # - nonstationary: logical TRUE/FALSE if the non-stationary cosine fit should also be performed
  #
  # Output:
  # dataframe with the following columns:
  # - time_min = time (minutes)
  # - period = period of the fitted sinusoid (minutes)
  # - DC = direct current (offset) of the fitted sinusoid
  # - a = fit parameter of the fitted sinusoid
  # - b = fit parameter of the fitted sinusoid
  # - RoO = range of oscillation of the fitted sinusoid
  # - phase = phase constant of the fitted sinusoid (degrees)
  # - cor = Pearson correlation coefficient between LIDSfitted and LIDS
  # - pvalue = pvalue of the correlation between LIDSfitted and LIDS
  # - MRI = Munich Rhythmicity Index (cor*RoO)
  # - LIDS = the LIDS time series derived from acc
  # - LIDSfitted = fitted sinusoid to LIDS (ft in Roenneberg 2015)
  # - cycle = cycle number relative to start in continuous decimal numbers
  # - cycle_full = cycle number demarking full cycles (i.e. as integer)
  # - RoO_cycle = range of oscillation for each cycle
  # - Max_cycle = maximum LIDS for each cycle
  # - Min_cycle = minimum LIDS for each cycle
  # - cycle_interpol = cycle number in increments of 0.01 for interpolation
  # - LIDS_norm_interpol = interpolated LIDS for each 0.01 cycle (i.e. cycle normalized)
  # - lm_intercept = intercept of the fitted line
  # - lm_slope = slope of the fitted line
  # - lm_residuals = residuals of the fitted line
  # - lm_MeanAmplitude = mean residuals of the fitted line = mean amplitude
  # - lm_fitted = fitted line to LIDS
  #-----------------------------------------------------------------
  # Outstanding action points -VvH:
  # - Avoid hardcoded assumpution that resolution is 1 minute, and instead make epoch setting
  #  flexible, but this also means that assumption of period unit = unit of epochs needs to be account for.
  #-----------------------------------------------------------------
  # Outstanding action points - EW:
  # 
  #-----------------------------------------------------------------
  
  # Calculate Locomotor Inactivity During Sleep (LIDS):
  #
  #   Step1: 
  #   Remove noise by subtracting threshold
  #   threshold 20mg = 0.02g
  
  ENMOsub = ifelse(acc < 20, 0, acc-20)
  
  #   Step2:
  #   10-minute filter: 10-minute rolling sum
  #   10min = 600s = 120*5s
  
  rollsumwindow = 10
  anyNA = which(is.na(ENMOsub) ==  TRUE)
  if (length(anyNA) > 0) ENMOsub[anyNA] = 0
  ENMOsub_smooth = zoo::rollapply(data=ENMOsub,width=(60*rollsumwindow)/ws3, FUN=sum, partial=T) #10-min rolling sum 
  #speedier alternative if one does not want a partial rollsum:
  #binaryclassification_smooth = zoo::rollsum(x=binaryclassification, k=(60*rollsumwindow)/ws3, fill=NA) #10-min rolling sum 
  
  #   Step3:
  #   LIDS conversion + 30-min filtering
  #   30min = 1800s = 360*5s
  
  LIDS_unfiltered = 100 / (ENMOsub_smooth + 1) 
  #when using ENMOsub in mg + at 5sec resolution + 10min rollsum, values are scaled sensibly for the following transformation
  LIDS = zoo::rollapply(data=LIDS_unfiltered,width=(60*30)/ws3, FUN=function(x) mean(x,na.rm=T), partial=T)
  
  
  #   Step4:
  #   Downsample to 1 minute resolution to speed up code  
  LIDS = LIDS[seq(1,length(LIDS),by=60/ws3)]
  
  #-----------------------------------------------------------------
  
  # Derive time series
  stepsize = 1 #1 minute
  time_min = (1:length(LIDS)) * stepsize
  
  #-----------------------------------------------------------------
  
  # Define function to perform cosine fit:
  cosfit = function(time_min,LIDS,period, nonstationary = FALSE) {
    
    # Cosine fitting inspired by:
    # Roenneberg T, Keller LK, et al. Human Activity and rest in Situ, 2015
    # Methods in Enzymology, Volume 552, http://dx.doi.org/10.1016/bs.mie.2014.11.028
    # 
    #   fit(t) = DC + A * sin(wt + phi0) = DC + a*cos(wt) + b*sin(wt)
    #   w=2pi/period (angular frequency)
    #   wt:normalised time in rad
    # 
    # input:
    # - x: vector with time time series
    # - y: vector with acc time series
    # - period: cosine period in minutes
    # output:
    # - list of varibales related to the cosine fitting process
    #-----------------------------------------
    
    
    #DC
    DC = mean(LIDS,na.rm= T)
    if (is.na(DC) == T) DC = 0
    LIDS = LIDS - DC
    
    #get normalised time in radians
    norm_t <- time_min*2*pi/period
    
    #Determine a and b
    a_t = LIDS * cos(norm_t); b_t = LIDS * sin(norm_t)
    a = mean(a_t) * 2; b = mean(b_t) * 2
    
    #Determine fit
    ft = a * cos(norm_t) + b * sin(norm_t)
    
    #Parameters from fit
    cor = 0 #Pearson correlation coeff
    pvalue = 1 #pvalue of correlation
    
    if (length(ft) > 0) {
      if (length(which(is.na(ft) == FALSE)) > 0) {
        if (sd(ft,na.rm = T) != 0 & sd(LIDS) != 0) {
          ct = stats::cor.test(ft,LIDS)
          pvalue = ct$p.value
          cor = ct$estimate
        }
      }
    }
    
    RoO = 2*sqrt(a^2+b^2)  #Range of oscillation
    MRI = RoO * cor #Munich Rhythmicity Index MRI
    
    #initial phase
    #the following provides zero phase angle (phase constant) in radians for a standard sine wave
    
    if(b==0) {
      
      if(a>0) phase <- pi/2
      else phase <- -pi/2
      
    }else{
      
      if(b>0) phase <- atan(a/b)
      else phase <- atan(a/b) - pi
    }
    
    phase <- phase*360/(2*pi)
    
    
    return(data.frame(period=period, DC=DC, a=a, b=b, RoO=RoO, phase=phase, cor=cor, pvalue=pvalue, MRI=MRI))
  }
  
  #---------------------------------------------------------
  
  #Perform cosine fits
  
  #------Stationary fit
  
  #stationary results dataframe
  LIDS_S = data.frame(period=NA, DC=NA, a=NA, b=NA, RoO=NA, 
                      phase=NA, cor=NA, pvalue=NA, MRI=NA)
  
  #setup period comparison (stationary) dataframe
  PCS = data.frame(period=LIDS_cosfit_periods, DC=NA, a=NA, b=NA, RoO=NA, 
                   phase=NA, cor=NA, pvalue=NA, MRI=NA)
  
  #loop over proposed period lengths (discrete series)
  for (j in 1:length(LIDS_cosfit_periods)) { 
    period = LIDS_cosfit_periods[j]
    PCS[j,] = cosfit(time_min,LIDS = LIDS, period, nonstationary = FALSE) #Apply cosine fit
  }
  
  #select best period 
  if (length(which(is.na(PCS$cor) == FALSE)) > 2) {
    
    if (fit.criterion.cosfit == 1) { #via highest Pearson correlation coefficient
      best = which.max(PCS$cor)[1]
      LIDS_S = PCS[best,]
      
    } else if (fit.criterion.cosfit == 2) { #via highest MRI that is a peak
      peak = numeric(0)
      MRI = PCS$MRI
      
      #find MRI peak
      for(k in 2:(length(MRI)-1)) {
        if(MRI[k] > MRI[k-1] & MRI[k] > MRI[k+1]) 
          peak = c(peak,MRI[k])
      }
      
      #find max peak
      if(length(peak) >0){ #if an MRI peak was found at all
        MaxMRI = max(peak)
        LIDS_S = PCS[PCS$MRI == MaxMRI & !is.na(PCS$MRI),]
      }
    }
  }
  
  #Add LIDS and LIDSfit timeline
  row.names(LIDS_S) <- NULL
  LIDS_S = cbind(time_min, LIDS_S, LIDS)
  LIDS_S$LIDSfitted = with(LIDS_S, DC + a * cos(time_min*2*pi/period) + b * sin(time_min*2*pi/period))
  
  #add further variables
  if (nrow(LIDS_S) >= 10 & length(which(is.na(LIDS_S$LIDSfitted) == FALSE)) > 10) {
    
    #Add cycle information
    LIDS_S$cycle = with(LIDS_S, time_min/period)
    LIDS_S$cycle_full = trunc(LIDS_S$cycle)+1
    
    #Range of oscillation (RoO) per cycle
    cycle <- LIDS_S$cycle_full
    LIDS_S$RoO_cycle <- NA
    LIDS_S$max_cycle <- NA
    LIDS_S$min_cycle <- NA
    max.cyc <- max(cycle, na.rm = TRUE)
    
    for(i in 1:max.cyc) {
      
      #only calculate RoO for complete cycles
      if(all(!is.na(LIDS[cycle == i])) && #if there is no missing LIDS for part of a cycle
         length(LIDS[cycle == i]) >= LIDS_S$period[1]-(LIDS_S$period[1]/20)){ #if cycle is not truncated (at end) by more than 5%
        
        LIDS_S$max_cycle[cycle == i] <- max(LIDS[cycle == i],na.rm=FALSE)
        LIDS_S$min_cycle[cycle == i] <- min(LIDS[cycle == i],na.rm=FALSE)
        LIDS_S$RoO_cycle <- with(LIDS_S, max_cycle - min_cycle)
      }    
    }
    
    #add linear model information
    lm.LIDS = lm(LIDS_S$LIDSfitted ~ LIDS_S$cycle_full)
    coef(lm.LIDS)[1]
    LIDS_S$lm_intercept = coef(lm.LIDS)[1]
    LIDS_S$lm_slope = coef(lm.LIDS)[2]
    LIDS_S$lm_residuals = resid(lm.LIDS)
    LIDS_S$lm_MeanAmplitude = sd(LIDS_S$lm_residuals)
    LIDS_S$lm_fitted = lm.LIDS$fitted.values
    
    #Interpolate to 100th of cycle length for later averaging for plotting
    #   beware: 100th is only sensible at 1 min resolution
    maxc = max(LIDS_S$cycle, na.rm=TRUE)
    cycle_interpol = seq(from=0.01,to=maxc, by=0.01)
    A = approx(LIDS_S$cycle,LIDS_S$LIDS,xout=cycle_interpol, rule=2, method ="linear", ties=mean)
    A = data.frame(cycle_interpol=A$x, LIDS_norm_interpol=A$y)
    A$time_min = seq(1,by=1, length.out = nrow(A))
    
    #add interpolated information
    #need to use function merge for combining as both datasets can differ in length
    LIDS_S <- merge(LIDS_S, A, by="time_min", all=TRUE)
    
    LIDSan = list(LIDS_S=LIDS_S)
  }
  
  #----Non-stationary fit
  #BEWARE: 
  #The following performs a cosine fit for each epoch using the data surrounding it by +/- period/2
  #Fits based on 1 cycle may not be reliable especially for period determination
  #Large edge effects reduce length of analysed time series by the maximum period length tested
  
  if(nonstationary == TRUE){#only performed if explicitly requested
    
    #non-stationary results dataframe
    LIDS_NS = data.frame(time_min=time_min, 
                         period=NA, DC=NA, a=NA, b=NA, RoO=NA, 
                         phase=NA, cor=NA, pvalue=NA, MRI=NA)
    
    #loop over time series epoch by epoch
    for (i in 1:length(LIDS)) { 
      
      #setup period comparison (non-stationary) dataframe
      PCNS = data.frame(period=LIDS_cosfit_periods, DC=NA, a=NA, b=NA, RoO=NA, 
                        phase=NA, cor=NA, pvalue=NA, MRI=NA)
      
      #loop over proposed period lengths (discrete series)
      for (j in 1:length(LIDS_cosfit_periods)) { 
        period = LIDS_cosfit_periods[j]
        halfp = (period)/2
        halfpmax = max(LIDS_cosfit_periods)/2
        if (i > halfpmax & i < ((length(LIDS)-halfpmax)+2) ) {#only perform fit when all periods that are to be tested can be applied
          y =  LIDS[(i-halfp):(i+halfp-1)]
          x =  time_min[(i-halfp):(i+halfp-1)]
          if (sd(y, na.rm = T) > 0) {
            PCNS[j,] = cosfit(time_min = x,LIDS = y,period, nonstationary = TRUE) #Apply cosine fit
          }
        }
      }
      
      #select best period 
      if (length(which(is.na(PCNS$cor) == FALSE)) > 2) {
        
        if (fit.criterion.cosfit == 1) { #via highest Pearson correlation coefficient
          PNS = PCNS[which.max(PCNS$cor)[1],]
          
        } else if (fit.criterion.cosfit == 2) { #via highest MRI that is a peak
          peak = numeric(0)
          MRI = PCNS$MRI
          
          #find MRI peak
          for(k in 2:(length(MRI)-1)) {
            if(MRI[k] > MRI[k-1] & MRI[k] > MRI[k+1]) 
              peak = c(peak,MRI[k])
          }
          
          #find max peak
          if(length(peak) >0){ #if an MRI peak was found at all
            MaxMRI = max(peak)
            PNS = PCNS[PCNS$MRI == MaxMRI & !is.na(PCNS$MRI),]
          }
        }
        LIDS_NS[i,2:ncol(LIDS_NS)] = PNS
      }
    }
    
    #Add LIDS and LIDSfit timeline
    LIDS_NS$LIDS = LIDS
    LIDS_NS$LIDSfitted = with(LIDS_NS, DC + a * cos(time_min*2*pi/period) + b * sin(time_min*2*pi/period))
    
    if (nrow(LIDS_NS) >= 10 & length(which(is.na(LIDS_NS$LIDSfitted) == FALSE)) > 10) {
      
      #Add cycle information
      LIDS_NS$cycle = c()
      LIDS_NS$cycle = LIDS_NS$phase + abs(min(LIDS_NS$phase,na.rm = TRUE))
      dcycle = abs(diff(LIDS_NS$cycle))
      dcycle[which(dcycle > 1)] = 0
      LIDS_NS$cycle[!is.na(LIDS_NS$cycle)] = c(cumsum(dcycle[!is.na(dcycle)]) / (2*pi), NA)
      LIDS_NS$cycle_full = trunc(LIDS_NS$cycle)+1
      
      #Create regularly spaced cycle points for which we want to interpolate the above variables
      minc = min(LIDS_NS$cycle,na.rm=TRUE)
      maxc = max(LIDS_NS$cycle, na.rm=TRUE)
      LIDS_NS$cycle_interpol = ((((1:nrow(LIDS_NS))-1 ) / nrow(LIDS_NS))* (maxc-minc)) + minc
      # interpolate LIDSfitted
      A = approx(LIDS_NS$cycle,LIDS_NS$LIDS, xout = LIDS_NS$cycle_interpol, rule = 2, method = "linear", ties = mean)
      LIDS_NS$LIDS_norm_interpol = A$y
      # interpolate period
      A = approx(LIDS_NS$cycle,LIDS_NS$period, xout = LIDS_NS$cycle_interpol, rule = 2, method = "linear", ties = mean)
      LIDS_NS$LIDSperiod_interpol = A$y
      
      LIDSan = list(LIDS_S=LIDS_S, LIDS_NS=LIDS_NS)
    }
    
  }#end of if(nonstationary == TRUE) 
  
  if (length(LIDSan) > 0) {# only report LIDS if LIDS analyses were successful
    return(invisible(LIDSan))
  } else {
    return()
  }
  
}#end of function g.LIDS.analyse

