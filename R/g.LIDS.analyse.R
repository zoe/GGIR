# By: Vincent van Hees, July 2018
g.LIDS.analyse = function(acc = c(), ws3 = 5) {
  # Input:
  # - acc: vector of acceleration time series
  # - ws3: epoch length of acc (seconds)
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
  # Calculate inactivity score (LIDS):
  binaryclassification = ifelse(acc < 20, 0, 1) # turn into binary classify per epoch
  rollsumwindow = 10
  
  anyNA = which(is.na(binaryclassification) ==  TRUE)
  if (length(anyNA) > 0) binaryclassification[anyNA] = 0
  binaryclassification_smooth = zoo::rollsum(x=binaryclassification,k=(60*rollsumwindow)/ws3) # 10 minute rolling sum
  # not using align = "center" in the previous step makes that we loose 5 minutes on each end of the time series
  # Nexpandmin = Nexpandmin - (rollsumwindow/2)
  inactivity = (10*(60/ws3)) / (binaryclassification_smooth + 1) #max value devided by (score + 1)
  #rescale to 0-100
  inactivity = inactivity / (0.6*rollsumwindow/ws3)
  LIDSraw = zoo::rollmean(x=inactivity,k=(60*30)/ws3,fill = "extend") # 30 minute rolling average
  # Downsample to 1 minute resolution to speed up code  
  LIDSraw = LIDSraw[seq(1,length(LIDSraw),by=60/ws3)]
  # TODO: Make epoch setting flexible, but this also
  # means that assumption of period unit = unit of epochs needs to be account for

  # Define function to perform cosine fit:
  cosfit = function(x,period) {
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
    DC = mean(x,na.rm= T)
    if (is.na(DC) == T) DC = 0
    xbackup = x
    x = x - DC
    norm_t = ((1:length(x)) * 2 * pi) / length(x) # 3. normalise time to period length to get t
    # Steps 1, 2, 3 are in the code below
    # Step 4. multiply timeseries within period (subset of LIDS2) with cos and sin of t, to get at and bt:
    a_t = x * cos(norm_t); b_t = x * sin(norm_t)
    # Step 5. calculate mean at and bt and multiply by two to get a and b:
    a = mean(a_t) * 2; b = mean(b_t) * 2
    # Step 6. calculate ft for each time point based on local a and b:
    ft = a * cos(norm_t) + b * sin(norm_t)
    if (sd(ft) == 0 | sd(x) == 0) {
      cor = 0
      pvalue = 1
    } else {
      ct = cor.test(ft,x)
      pvalue = ct$p.value
      cor = ct$estimate
    }
    # Note: we do not multiply the correlation by the range in the data
    # because that seems to result in a very eratic score.
    phase= atan(a/b) # should always be: -pi/2 < phase < pi/2
    LIDSfitted = ft[halfp]
    return(invisible(list(period=period,phase=phase,DC=DC,LIDSfitted=LIDSfitted,pvalue=pvalue,cor=cor)))
  }
  cntt = 1
  LIDSvars = data.frame(time=1:length(LIDSraw), period=NA, phase=NA, LIDSfitted=NA, DC=NA,pvalue=NA,cor=NA)
  for (i in 1:length(LIDSraw)) { # Step 1. loop over SPT window epoch by epoch
    periods = seq(30,180,by=5)
    periodcomp = data.frame(period=periods,phase=NA,DC=NA,LIDSfitted=NA,pvalue=NA,cor=NA)
    for (j in 1:length(periods)) { # Step 2. loop over proposed period lengths (discrete series)
      period = periods[j]
      halfp = (period)/ 2
      if (i > halfp & i < ((length(LIDSraw)-halfp)+2)) {
        x = LIDSraw[(i-halfp):(i+halfp-1)]
        periodcomp[j,] = cosfit(x,period) # Step 3-6: Apply cosine fit
      }
    }
    if (length(which(is.na(periodcomp$cor) == FALSE)) > 2) {
      # Step 7. select best period and from that you know instantanious phase
      bestperiod = periodcomp[which.max(periodcomp$cor)[1],]
      LIDSvars$period[i] = bestperiod$period
      LIDSvars$phase[i] = bestperiod$phase
      LIDSvars$LIDSfitted[i] = bestperiod$LIDSfitted
      LIDSvars$DC[i] = bestperiod$DC
      LIDSvars$pvalue[i] = bestperiod$pvalue
      LIDSvars$cor[i] = bestperiod$cor
    }
  }
  LIDSvars$LIDSraw = LIDSraw
  # remove Plateau at the end
  lastvalue = LIDSvars$LIDSraw[length(LIDSvars$LIDSraw)]
  startPlateau = length(LIDSvars$LIDSraw) - which((rev(LIDSvars$LIDSraw)==lastvalue) == FALSE)[1]
  if (startPlateau < length(LIDSvars$LIDSraw)) {
    LIDSvars = LIDSvars[1:startPlateau,] 
  }
  LIDSvars = LIDSvars[which(LIDSvars$cor > 0.8)[1]:nrow(LIDSvars),]
  LIDSvars$cycle = c()
  LIDSvars$cycle = LIDSvars$phase + abs(min(LIDSvars$phase,na.rm = TRUE))
  LIDSvars$cycle[which(is.na(LIDSvars$cycle) == TRUE)] = 0
  dcycle = abs(diff(LIDSvars$cycle))
  dcycle[which(dcycle > 1)] = 0
  LIDSvars$cycle[1:(nrow(LIDSvars)-1)] = cumsum(dcycle) / (2*pi)
  # absolute fitted LIDS line, meaning: including the original offset (=DC component)
  LIDSvars$LIDSfitted_abs = LIDSvars$LIDSfitted+ LIDSvars$DC
  # create regularly spaced cycle points for which we want to interpolate the above variables"
  minc = min(LIDSvars$cycle)
  maxc = max(LIDSvars$cycle)
  LIDSvars$cycle_interpol = ((((1:nrow(LIDSvars))-1 ) / nrow(LIDSvars))* (maxc-minc)) + minc
  # interpolate LIDSfitted_abs
  A = approx(LIDSvars$cycle,LIDSvars$LIDSfitted_abs, xout = LIDSvars$cycle_interpol, rule = 2, method = "linear", ties = mean)
  LIDSvars$LIDSfitted_abs_norm_interpol = A$y
  # interpolate LIDSfitted
  A = approx(LIDSvars$cycle,LIDSvars$LIDSfitted, xout = LIDSvars$cycle_interpol, rule = 2, method = "linear", ties = mean)
  LIDSvars$LIDSfitted_norm_interpol = A$y
  # interpolate period
  A = approx(LIDSvars$cycle,LIDSvars$period, xout = LIDSvars$cycle_interpol, rule = 2, method = "linear", ties = mean)
  LIDSvars$LIDSperiod_interpol = A$y

  return(LIDSvars)
}