# function g.detect.sleepbout
# taken from https://github.com/wadpac/GGIR/blob/issue95_sleepstructure/R/g.detect.sleepbout.R
# 24th Oct 2018


g.detect.sleepbout = function(WakeBinary=c(), WakeBout.threshold=0.5, WakeBoutMin=30, SleepBoutMin=180,ws3=5) {
    # Detect sleep bouts (here defined as Bouts of at least 3 hours with < 50% wakefullness).
    # and duration of the sleep bout that follows the sleep bout within the Sleep Period Time window (SPT)
    #------------------------------------
    # Apply 30 minute rolling average to derive ratio of wakefulleness
    WakeBinary = zoo::rollmean(x=WakeBinary,k=(30*(60/ws3)),fill = "extend",align="center")
    # Now use threshold to make binary distinction again between sleep bout, add dummy data to the end
    WakeBinary = c(1,0,ifelse(test=WakeBinary<WakeBout.threshold,yes = 0,no = 1),rep(1,WakeBoutMin*10*(60/ws3)),0,1)
    # Values use during testing: WakeBinary = c(1,rep(0,1200),rep(1,40),rep(0,1400),rep(1,35),rep(0,40),1)
    sleepbouts = matrix(0,5,3) # create empty matrix to be filled with sleep bout start and end times, and wakeduration that follows it
    sleepbouti = 1
    if (length(WakeBinary == 0) > 0) { # We are interested in the Bouts of zero (sleep)
        changepoints = which(abs(diff(WakeBinary)) == 1) # where does value change?
        state = rep(c(0,1),ceiling(length(changepoints)/2)) #0 = sleep, 1 = wake
        state = state[1:length(changepoints)-1]
        statedur = diff(changepoints) #state duration
        # identify index in changepoints for which the state that folows it is 1 and the duration is less than WakeBoutMin
        tooshortwake = which(state == 1 & statedur <= (WakeBoutMin * (60/ws3)))
        tooshortwake = unique(c(tooshortwake,tooshortwake + 1)) # get start and end of those bouts
        if (length(tooshortwake) > 0) changepoints = changepoints[-tooshortwake] # remove the short wakefullness bouts
        BoutlengthsWake = diff(changepoints)[seq(2,length(changepoints)-1,by=2)] # calculates length of Wake Bouts
        longWake = which(BoutlengthsWake >= (WakeBoutMin * (60/ws3))) # identifies the Wake Bouts that meet the criteria
        for (LWi in longWake) { # loop over long wakeBouts
            B = changepoints[((LWi*2)-1):((LWi*2))] # derive start and end of Sleep bouts that precedes Wake bout
            if (length(changepoints) > ((LWi*2)+1)) {
                C = changepoints[((LWi*2)):((LWi*2)+1)] # derive start and end of Wake bout
                C[1] = C[1] + 1
            } else {
                C = c(0,0)
            }
            if (diff(B) > (SleepBoutMin*(60/ws3))) { # Check whether it is at least SleepBoutMin (default 3 hours)
                WakeDurAfter = BoutlengthsWake[LWi]/(60/ws3)
                WakeDurAfter = WakeDurAfter - (WakeBoutMin*10) # remove dummy data
                sleepbouts[sleepbouti,] = c(B, WakeDurAfter) #Store these in the output matrix
                sleepbouti = sleepbouti + 1
                if (sleepbouti > (nrow(sleepbouts) - 3)) sleepbouts = rbind(sleepbouts,matrix(0,5,3))
            }
        }
    }
    return(sleepbouts)
}