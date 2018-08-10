\name{g.detect.sleepbout}
\alias{g.detect.sleepbout}
\title{
  Detect sleep bouts within the Sleep Period Time.
}
\description{
  Detect sleep bouts from binary wake/sleep input to facilitate LIDS analyses.
}
\usage{
  g.detect.sleepbout(WakeBinary=c(), WakeBout.threshold=0.5,
            WakeBoutMin=30, SleepBoutMin=180, ws3=5)
}
\arguments{
  \item{WakeBinary}{
    A vector with binary distinction between awake and sleep (wake = 1, sleep  = 0)
  }
  \item{WakeBout.threshold}{
    Fraction of the Wakefullness bout that needs to be wakefullness.
  }
    \item{ WakeBoutMin}{
    Minimum duration of the Wakefullness Bout (minutes)
  }
    \item{SleepBoutMin}{
    Minimum duration of the Sleep Bout (minutes)
  }
    \item{ws3}{
    Integer to indicate the epoch length of the WakeBinary vector in seconds.
  }
}
\details{
  A sleepbout is defined here as a period of at least 3 hours in which the rolling
  30 minute window consists of less than 50 percent wakefullness.
}
\value{
  matrix with start and end index of sleep bouts, and as third column the duration of
  wakefullness bout that follows the slpee bout.
}

\examples{
  WakeBinary = c(rep(0,12*60*4),rep(1,12*40),rep(0,12*60,2)) # dummy input
  sleepbouts = g.detect.sleepbout(WakeBinary=WakeBinary, 
                WakeBout.threshold=0.5, WakeBoutMin=30, SleepBoutMin=180,ws3=5)
  print(sleepbouts)
}
\author{
Vincent T van Hees <vincentvanhees@gmail.com>
}
