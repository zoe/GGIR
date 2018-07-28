g.detect.sleepbout = function(WakeBinary=c(),WakeBout.threshold=0.5,WakeBoutMin=30,SleepBoutMin=180,ws3=5) {
  # Detect sleep bouts (here defined as Bouts of at least 3 hours with < 30% wakefullness).
  #------------------------------------
  # Apply 30 minute rolling average to derive ratio of wakefulleness
  WakeBinary = zoo::rollmean(x=WakeBinary,k=(30*(60/ws3)),fill = "extend",align="center")
  # Now use threshold to make binary distinction again between sleep bout, add dummy data to the end
  WakeBinary = c(1,0,ifelse(test=WakeBinary<WakeBout.threshold,yes = 0,no = 1),rep(1,WakeBoutMin*2*(60/ws3)),0,1)
  # Same values use during testing: WakeBinary = c(1,rep(0,1200),rep(1,40),rep(0,1400),rep(1,35),rep(0,40),1)
  sleepbouts = matrix(0,5,2) # create empty matrix to be filled with sleep bout start and end times
  sleepbouti = 1
  if (length(WakeBinary == 0) > 0) { # We are interested in the Bouts of zero (sleep)
    changepoints = which(abs(diff(WakeBinary)) == 1) # where does value change?
    BoutlengthsWake = diff(changepoints)[seq(2,length(changepoints)-1,by=2)] # calculates length of Wake Bouts
    longWake = which(BoutlengthsWake > (WakeBoutMin * (60/ws3))) # identifies the Wake Bouts that meet the criteria
    for (LWi in 1:length(longWake)) { # loop over long wakeBouts
      B = changepoints[((longWake[LWi]*2)-1):((longWake[LWi]*2))] # calculates length of Sleep bouts
      if (diff(B) > (SleepBoutMin*(60/ws3))) { # Check whether it is at least SleepBoutMin (default 3 hours)
        sleepbouts[sleepbouti,] = B #Store these in the output matrix
        sleepbouti = sleepbouti + 1
        if (sleepbouti > (nrow(sleepbouts) - 3)) sleepbouts = rbind(sleepbouts,matrix(0,5,2))
      }
    }
  }
  return(sleepbouts)
}