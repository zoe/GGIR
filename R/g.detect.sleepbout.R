g.detect.sleepbout = function(WakeBinary=c(),ws3=5,WakeBout.threshold=0.7,WakeBoutMin=30,SleepBoutMin=180) {
  # Detect sleep bouts (here defined as Bouts of at least 3 hours with < 30% wakefullness).
  # Input:
  # - WakeBinary: A binary distinction between awake and sleep (wake = 1, sleep  = 0)
  # - WakeBout.threshold: Fraction of the sleep bout that eneds to be asleep
  # - WakeBoutMin: Minimum duration of the WakeBout (minuntes)
  # - SleepBoutMin: Minimum duration of the SleepBout (minuntes)
  # Output:
  # - matrix with start and end index of sleep bouts
  #------------------------------------
  # Apply 30 minute rolling average to derive ratio of wakefulleness
  WakeBinary = zoo::rollmean(x=WakeBinary,k=(30*(60/ws3)),fill = "extend")
  # Now use threshold to make binary distinction again between sleep bout, add dummy data to the end
  WakeBinary = c(1,0,ifelse(test=WakeBinary<WakeBout.threshold,yes = 0,no = 1),rep(1,WakeBoutMin*2*(60/ws3)),0,1)
  # Same values use during testing: WakeBinary = c(1,rep(0,1200),rep(1,40),rep(0,1400),rep(1,35),rep(0,40),1)
  sleepbouts = matrix(0,20,2) # create empty matrix to be filled with sleep bout start and end times
  sleepbouti = 1
  if (length(WakeBinary == 0) > 0) { # We are interested in the Bouts of zero (sleep)
    changepoints = which(abs(diff(WakeBinary)) == 1) # where does value change?
    BoutlengthsWake = diff(changepoints)[seq(2,length(changepoints)-1,by=2)] # calculates length of Wake Bouts
    longWake = which(BoutlengthsWake > (WakeBoutMin * (60/ws3))) # identifies the Wake Bouts that meet the criteria
    for (LWi in 1:length(longWake)) { # loop over long wakeBouts
      B = changepoints[((longWake[LWi]*2)-1):((longWake[LWi]*2))] # calculates length of Sleep bouts
      # print(paste0(diff(B) / (60*(60/ws3))," hours"))
      if (diff(B) > (SleepBoutMin*(60/ws3))) { # Check whether it is at least SleepBoutMin (default 3 hours)
        sleepbouts[sleepbouti,] = B #Store these in the output matrix
        sleepbouti = sleepbouti + 1
      }
    }
  }
  return(sleepbouts)
}