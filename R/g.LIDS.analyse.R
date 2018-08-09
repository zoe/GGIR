# By: Vincent van Hees, July 2018
g.LIDS.analyse = function(acc = c(), ws3 = 5, best.fit.criterion.cosine = 2, 
                          min_period = 30, max_period = 180, step_period = 5) {
  # Input:
  # - acc: vector of acceleration time series
  # - ws3: epoch length of acc (seconds)
  # - best.fit.criterion.cosine: metric to select the best LIDS (1=cor,2=cor*range)
  # - min_period: minimum period length to consider (minutes)
  # - max_period: maximum period length to consider (minutes)
  # - step_period: stepsize with which to interpolate min and max.
  # Output:
  # dataframe with the following columns:
  # - time = time (minutes)
  # - period = instantanious period of the fitted sinoid (minutes)
  # - phase = instantanious phase of the fitted sinoid (radians)
  # - LIDSfitted = instantanious fitted LIDS (ft in Roenneberg 2015)
  # - DC = instantanious direct current (offset) of the fitted sinoid)
  # - cor = instantanious correlation coefficients between LIDSfitted and LIDSraw
  # - pvalue = instantanious pvalue of the correlation between LIDSfitted and LIDSraw
  # - LIDSraw = inactivity score derived from acc
  # - cycle = instantanious cycle relative to start
  # - LIDSfitted_abs = LIDSfitted + DC
  # - cycle_interpol = interpolated cycle
  # - LIDSfitted_abs_norm_interpol = interpolated and cycle normalized LIDSfitted_abs
  # - LIDSfitted_norm_interpol = interpolated and cycle normalized LIDSfitted
  # - LIDSperiod_interpol = interpolated period
  #-----------------------------------------------------------------
  # Outstanding action points:
  # - Avoid hardcoded assumpution that resolution is 1 minute, and instead make epoch setting
  #  flexible, but this also means that assumption of period unit = unit of epochs needs to be account for.
  #-----------------------------------------------------------------
  # Calculate inactivity score (LIDS):
  binaryclassification = ifelse(acc < 20, 0, 1) # turn into binary classify per epoch
  rollsumwindow = 10
  anyNA = which(is.na(binaryclassification) ==  TRUE)
  if (length(anyNA) > 0) binaryclassification[anyNA] = 0
  binaryclassification_smooth = zoo::rollsum(x=binaryclassification,k=(60*rollsumwindow)/ws3) # 10 minute rolling sum
  # not using align = "center" in the previous step makes that we loose 5 minutes on each end of the time series
  # rescale to 0-100
  binaryclassification_smooth = binaryclassification_smooth / ((rollsumwindow*(60/ws3)) / 100)
  inactivity = 100 / (binaryclassification_smooth + 1) #max value devided by (score + 1)
  
  LIDSraw = zoo::rollmean(x=inactivity,k=(60*30)/ws3,fill = "extend") # 30 minute rolling average
  # Downsample to 1 minute resolution to speed up code  
  LIDSraw = LIDSraw[seq(1,length(LIDSraw),by=60/ws3)]
  #-----------------------------------------------------------------
  # Derive time series
  stepsize = 1 #1 minute
  time_min = (1:length(LIDSraw)) * stepsize
  #-----------------------------------------------------------------
  # Define function to perform cosine fit:
  cosfit = function(time_min,LIDS,period,stationary = FALSE) {
    # Cosine fitting inspired by:
    # Roenneberg T, Keller LK, et al. Human Activity and rest in Situ, 2015
    # Methods in Enzymology, Volume 552, http://dx.doi.org/10.1016/bs.mie.2014.11.028
    # input:
    # - x: vector with time series
    # - period: cosine period in minutes
    # output:
    # - list of varibales related to the cosine fitting process
    #-----------------------------------------
    halfp = period/ 2 #*NumEpochPerMin
    DC = mean(LIDS,na.rm= T)
    if (is.na(DC) == T) DC = 0
    LIDS = LIDS - DC
    if (stationary == FALSE) {
      norm_t = ((1:length(LIDS)) * 2 * pi) / length(LIDS) # 3. normalise time to period length to get t
    } else {
      norm_t <- time_min*2*pi/period
    }
    # Steps 1, 2, 3 are in the code below
    # Step 4. multiply timeseries within period (subset of LIDS2) with cos and sin of t, to get at and bt:
    a_t = LIDS * cos(norm_t); b_t = LIDS * sin(norm_t)
    # Step 5. calculate mean at and bt and multiply by two to get a and b:
    a = mean(a_t) * 2; b = mean(b_t) * 2
    # Step 6. calculate ft for each time point based on local a and b:
    ft = a * cos(norm_t) + b * sin(norm_t)
    cor = 0
    pvalue = 1
    # RoO = abs(diff(range(LIDS,na.rm = TRUE))) # I used to have x
    RoO = abs(diff(range(ft,na.rm = TRUE))) # I used to have x
    if (length(ft) > 0) {
      if (length(which(is.na(ft) == FALSE)) > 0) {
        if (sd(ft) != 0 & sd(LIDS) != 0) {
          ct = stats::cor.test(ft,LIDS)
          pvalue = ct$p.value
          cor = ct$estimate
        }
      }
    }
    MRI <- RoO * cor #Munich Rhythmicity Index MRI
    # Note: we do not multiply the correlation by the range in the data
    # because that seems to result in a very eratic score.
    phase= atan(b/a) # should always be: -pi/2 < phase < pi/2
    if (stationary == FALSE) {
      LIDSfitted = ft[halfp]
    } else {
      LIDSfitted = NA # not applicable
    }
    return(invisible(list(period=period,phase=phase,DC=DC,LIDSfitted=LIDSfitted,pvalue=pvalue,cor=cor,MRI=MRI,RoO=RoO)))
  }
  periods = seq(min_period,max_period,by=step_period)
  LIDS_S = LIDS_NS = data.frame(time=1:length(LIDSraw), period=NA, phase=NA, 
                                LIDSfitted=NA, DC=NA, pvalue=NA, cor=NA, MRI=NA, RoO=NA)
  #==========================================
  # Non stationary
  cntt = 1
  for (i in 1:length(LIDSraw)) { # Step 1. loop over SPT window epoch by epoch
    #PCNS: Period comparison non stationary
    PCNS = data.frame(period=periods, phase=NA, DC=NA, LIDSfitted=NA, pvalue=NA, cor=NA, MRI=NA,RoO=NA)
    for (j in 1:length(periods)) { # Step 2. loop over proposed period lengths (discrete series)
      period = periods[j]
      halfp = (period)/ 2
      if (i > halfp & i < ((length(LIDSraw)-halfp)+2) ) {
        y =  LIDSraw[(i-halfp):(i+halfp-1)]
        if (sd(y) > 0) {
          PCNS[j,] = cosfit(time_min,LIDS = y,period, stationary = FALSE) # Step 3-6: Apply cosine fit
        }
      }
    }
    if (length(which(is.na(PCNS$cor) == FALSE)) > 2) {
      # Step 7. select best period and from that you know instantanious phase
      if (best.fit.criterion.cosine == 1) {
        bestperiod = PCNS[which.max(PCNS$cor)[1],]
      } else if (best.fit.criterion.cosine == 2) {
        bestperiod = PCNS[which.max(PCNS$MRI)[1],]
      }
      LIDS_NS[i,2:ncol(LIDS_NS)] = bestperiod
    }
  }
  LIDS_NS$LIDSraw = LIDSraw
  #==========================================
  # Stationary
  # PCS: Period comparison stationary
  PCS = data.frame(period=periods, phase=NA, DC=NA, LIDSfitted=NA, pvalue=NA, cor=NA, MRI=NA, RoO=NA)
  for (j in 1:length(periods)) { # Step 2. loop over proposed period lengths (discrete series)
    period = periods[j]
    PCS[j,] = cosfit(time_min,LIDS = LIDSraw,period, stationary = TRUE) # Step 3-6: Apply cosine fit
  }
  if (length(which(is.na(PCS$cor) == FALSE)) > 2) {
    # Step 7. select best period and from that you know instantanious phase
    if (best.fit.criterion.cosine == 1) {
      besti = which.max(PCS$cor)[1]
      LIDS_S = PCS[besti,]
    } else if (best.fit.criterion.cosine == 2) {
      LIDS_S = PCS[which.max(PCS$MRI)[1],]
    }
  }
  LIDS_S = LIDS_S[,-c(which(colnames(LIDS_S) == "LIDSfitted"))]
  #------------------------------------------------------------
  # remove Plateau at the end if there is one
  indexLastValue = max(which(is.na(LIDS_NS$LIDSraw) == FALSE)) # find the last value that is not a NA 
  lastvalue = LIDS_NS$LIDSraw[indexLastValue] #[length(LIDS_NS$LIDSraw)]
  revLIDSraw = rev(LIDS_NS$LIDSraw)
  ReverseIndexFirstValuePlateau = which(revLIDSraw != lastvalue & is.na(revLIDSraw) == FALSE)[1]
  if (length(ReverseIndexFirstValuePlateau) > 0) {
    if (is.na(ReverseIndexFirstValuePlateau) == FALSE) {
      startPlateau = indexLastValue - ReverseIndexFirstValuePlateau
      if (length(startPlateau) > 0 & nrow(LIDS_NS) > 10) {
        if (startPlateau < length(LIDS_NS$LIDSraw)) { # kind of obvious, can be removed?
          LIDS_NS = LIDS_NS[1:startPlateau,] 
        }
      }
    }
  }
  # LIDS_NS = LIDS_NS[which(LIDS_NS$cor > 0.8)[1]:nrow(LIDS_NS),]
  if (nrow(LIDS_NS) >= 10 & length(which(is.na(LIDS_NS$LIDSfitted) == FALSE)) > 10) {
    LIDS_NS$cycle = c()
    LIDS_NS$cycle = LIDS_NS$phase + abs(min(LIDS_NS$phase,na.rm = TRUE))
    LIDS_NS$cycle[which(is.na(LIDS_NS$cycle) == TRUE)] = 0
    dcycle = abs(diff(LIDS_NS$cycle))
    dcycle[which(dcycle > 1)] = 0
    LIDS_NS$cycle[1:(nrow(LIDS_NS)-1)] = cumsum(dcycle) / (2*pi)
    LIDS_NS$cycle_full = trunc(LIDS_NS$cycle)+1
    
    # absolute fitted LIDS line, meaning: including the original offset (=DC component)
    LIDS_NS$LIDSfitted_abs = LIDS_NS$LIDSfitted+ LIDS_NS$DC
    # create regularly spaced cycle points for which we want to interpolate the above variables"
    minc = min(LIDS_NS$cycle)
    maxc = max(LIDS_NS$cycle)
    LIDS_NS$cycle_interpol = ((((1:nrow(LIDS_NS))-1 ) / nrow(LIDS_NS))* (maxc-minc)) + minc
    # interpolate LIDSfitted_abs
    A = approx(LIDS_NS$cycle,LIDS_NS$LIDSfitted_abs, xout = LIDS_NS$cycle_interpol, rule = 2, method = "linear", ties = mean)
    LIDS_NS$LIDSfitted_abs_norm_interpol = A$y
    # interpolate LIDSfitted
    A = approx(LIDS_NS$cycle,LIDS_NS$LIDSfitted, xout = LIDS_NS$cycle_interpol, rule = 2, method = "linear", ties = mean)
    LIDS_NS$LIDSfitted_norm_interpol = A$y
    # interpolate period
    A = approx(LIDS_NS$cycle,LIDS_NS$period, xout = LIDS_NS$cycle_interpol, rule = 2, method = "linear", ties = mean)
    LIDS_NS$LIDSperiod_interpol = A$y
    return(invisible(list(LIDS_NS=LIDS_NS,LIDS_S=LIDS_S)))

  } else {
    return()
  }
}