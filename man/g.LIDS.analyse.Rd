\name{g.LIDS.analyse}
\alias{g.LIDS.analyse}
\title{
  Analyse Locomoter Inactivity During Sleep (LIDS)
}
\description{
  Function to perform LIDS analysis or also referred to as Ultradian structure
  analyses inspired by Winnebeck et al. 2018 https://doi.org/10.1016/j.cub.2017.11.063
  and Roenneberg et al. http://dx.doi.org/10.1016/bs.mie.2014.11.028.
}
\usage{
  g.LIDS.analyse(acc = c(), ws3 = 5)
}
\arguments{
  \item{acc}{
    Magnitude of body acceleration aggregated per epoch
  }
  \item{ws3}{
    Duration (seconds) of the epoch one value of acc corresponds to
  }
}
\value{
\item{\code{LIDSvars}}{A data.frame with the LIDS variabls, see details}
}

\details{
The value \code{LIDSvars} is a dataframe and comes with the following variables:\cr
\itemize{
\item \code{time} time (minutes)
\item \code{period} instantanious period of the fitted sinoid (minutes)
\item \code{phase} instantanious phase of the fitted sinoid (radians)
\item \code{LIDSfitted} instantanious fitted LIDS (ft in Roenneberg 2015)
\item \code{DC} instantanious direct current (offset) of the fitted sinoid)
\item \code{cor} instantanious correlation coefficients between LIDSfitted and LIDSraw
\item \code{pvalue} instantanious pvalue of the correlation between LIDSfitted and LIDSraw
\item \code{LIDSraw} inactivity score derived from acc
\item \code{cycle} instantanious cycle relative to start
\item \code{LIDSfitted_abs} LIDSfitted + DC
\item \code{cycle_interpol} interpolated cycle
\item \code{LIDSfitted_abs_norm_interpol} interpolated and cycle normalized LIDSfitted_abs
\item \code{LIDSfitted_norm_interpol} interpolated and cycle normalized LIDSfitted
\item \code{LIDSperiod_interpol} interpolated period
}
}
\examples{
  acc = abs(rnorm(n = 1200,mean = 20,sd = 5)) # generate random numbers to test
  LIDSvars = g.LIDS.analyse(acc=acc,ws3=5)
}
\author{
Vincent T van Hees <vincentvanhees@gmail.com>
}
