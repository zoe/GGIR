\name{g.LIDS.analyse_ENMOsub}
\alias{g.LIDS.analyse_ENMOsub}
\title{
  Analyse Locomoter Inactivity During Sleep (LIDS)
}
\description{
  Function to perform LIDS analysis or also referred to as Ultradian structure
  analyses inspired by Winnebeck et al. 2018 https://doi.org/10.1016/j.cub.2017.11.063
  and Roenneberg et al. http://dx.doi.org/10.1016/bs.mie.2014.11.028.
}
\usage{
  g.LIDS.analyse_ENMOsub(acc = c(), ws3 = 5,fit.criterion.cosfit = 2, 
                          LIDS_cosfit_periods = seq(30,180,by=5),
                          nonstationary=FALSE)
}
\arguments{
  \item{acc}{
    Magnitude of body acceleration aggregated per epoch
  }
  \item{ws3}{
    Duration (seconds) of the epoch one value of acc corresponds to
  }
  \item{fit.criterion.cosfit}{
    Metric to select the best cosine fit (1=correlation ,2=correlation x range)
  }
  \item{LIDS_cosfit_periods}{
   vector with period lengths (minutes) to consider for the cosine fitting.
  }
  \item{nonstationary}{
   Boolean to indicate whether nonstationary analysis should be performed
  }
  
}
\value{
\item{\code{LIDSvars}}{A data.frame with the LIDS variabls, see details}
}

\details{
The value \code{LIDSvars} is a dataframe and comes with the following variables:\cr
\itemize{
\item \code{time} time (minutes)
\item \code{period} instantanious period of the fitted cosine (minutes)
\item \code{phase} instantanious phase of the fitted cosine (radians)
\item \code{LIDSfitted} instantanious fitted LIDS (ft in Roenneberg 2015)
\item \code{DC} instantanious direct current (offset) of the fitted cosine)
\item \code{cor} instantanious correlation coefficients between LIDSfitted and LIDSraw
\item \code{pvalue} instantanious pvalue of the correlation between LIDSfitted and LIDSraw
\item \code{LIDSraw} inactivity score derived from acc
\item \code{acc_aggregated} acc value aggregated: acc was first smoothed with a 30 minute
rolling mean and then aggregated as the mean per epoch, where the epoch is equal to 
the epoch length of LIDSraw (1 minute by default.
\item \code{cycle} instantanious cycle relative to start
\item \code{LIDSfitted_abs} LIDSfitted + DC
\item \code{cycle_interpol} interpolated cycle
\item \code{LIDSfitted_abs_norm_interpol} interpolated and cycle normalized LIDSfitted_abs
\item \code{LIDSfitted_norm_interpol} interpolated and cycle normalized LIDSfitted
\item \code{LIDSperiod_interpol} interpolated period
}
The value \code{LIDSstationary} is a lists with the following variables:\cr
hase=SCF$phase,period = per,amplitude = amplitude,cor = ct2$cor, pvalue = ct2$pvalue)
\itemize{
\item \code{period} period of the fitted cosine (minutes)
\item \code{phase} phase of the fitted cosine (radians)
\item \code{amplitude} amplitude of the fitted cosine
\item \code{cor} correlation coefficient between fitted and raw
\item \code{pvalue} pvalue of the correlation between fitted and raw
}
}
\examples{
  acc = abs(rnorm(n = 1200,mean = 20,sd = 5)) # generate random numbers to test
  LIDSvars = g.LIDS.analyse_ENMOsub(acc=acc,ws3=5,
              fit.criterion.cosfit = 2, 
              LIDS_cosfit_periods =seq(30,180,by=5),nonstationary=FALSE)
}
\author{
Vincent T van Hees <vincentvanhees@gmail.com>
}
