% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpt_ar.R
\name{cpts_ar}
\alias{cpts_ar}
\title{Identifying changes in mean and variance of accumulation-rate records}
\usage{
cpts_ar(
  serie,
  series_name = NULL,
  proxy = NULL,
  n_rand = 1000,
  bootstrap = FALSE,
  rand_ds_distr = "uniform",
  max_cpts = 50,
  n_screen = 0.025,
  tau = NULL,
  meth_cpt = "BinSeg",
  pen_method = "Manual",
  pen_val = "4*log(n)",
  t_stat = "Normal",
  minseglen = 2
)
}
\arguments{
\item{serie}{A list returned from the \code{pretreatment_data()}
function. NB: the \code{pretreatment_data()} function must be run
with the argument \code{out = "accI"}.}

\item{series_name}{A character string giving the name of the data set,
or typically the site name (optional).}

\item{proxy}{A character string giving the name (typically the column name)
of the variable that should be analysed.}

\item{n_rand}{Positive integer giving the number of iterations used to check
for change-points in random data sets (by default
\code{n_rand = 1000}).}

\item{bootstrap}{Logical. Determines how the random data sets are generated.
If TRUE, then each of the \code{n_rand} random data sets are
generated by random sampling with replacement from the input data.
If FALSE, then the code generates in each of the \code{n_rand}
iterations a vector of random concentration values having a
user-specified distribution (based on \code{rand_ds_distr}).}

\item{rand_ds_distr}{A character string giving the distribution used
to generate the random data sets when \code{bootstrap = FALSE}.
Currently only "uniform" is supported (default). The random values
range between the minimum and the maximum concentration values of
the re-sampled series.}

\item{max_cpts}{Positive integer giving the maximum number of change points
to search for when \code{meth_cpt = "BinSeg"}. By default
\code{max_cpts = 50}).}

\item{n_screen}{Positive number giving the minimum frequency of occurrence of
change-points in the random data sets to be validated. By default,
\code{n_screen = 0.025}. Thus, with \code{n_rand = 1000} this equals
a 2.5\% chance of occurrence.}

\item{tau}{Positive integer giving the tolerance search window (in years)
for change points in the 'true proxy-accumulation rate' record. By
default \code{tau} = the number of years the raw record was re-sampled
to (as of variable \code{yr.interp} that is written in the list
generated with the \code{pretreatment_data()} function).}

\item{meth_cpt}{As of \code{changepoint::cpt.meanvar()}: choice of
"AMOC", "PELT", "BinSeg", or "SegNeigh". But note that "SegNeigh" is
computationally slow, use other methods instead.}

\item{pen_method}{Refer to \code{changepoint::cpt.meanvar()} for
argument \code{penalty}. By default, \code{pen_method = "Manual"}.}

\item{pen_val}{Refer to \code{changepoint::cpt.meanvar()} for the
argument \code{pen.value}. By default, \code{pen_val = "4*log(n)"},
where \code{n} is equal to the number of observations.}

\item{t_stat}{Refer to \code{changepoint::cpt.meanvar()} for the
argument \code{test.stat}. Currently only "Normal" is supported
(default).}

\item{minseglen}{Positive integer giving the minimum number of observations
within zones. By default \code{minseglen = 2}, see also
\code{changepoint::cpt.meanvar()}.}
}
\value{
a list with the depths and ages of the change points, the mean values
for the zones, as well as the settings of the change-point analysis.

In addition, a diagnostic plot consisting of three panels is sent to the
device.
The plot includes (from top to bottom):
\itemize{
\item the re-sampled concentration values,
\item the sediment-accumulation rate, and
\item the proxy accumulation-rate record.
}

The latter panel may show:
\itemize{
\item change points of the proxy accumulation-rate record (blue vertical lines),
\item mean proxy accumulation-rate values for the zones (red horizontal lines),
\item non-sense change points that match 'true change points' (red circles),
\item non-sense change points that do not match 'true change points' (green
circles),
}
}
\description{
Determines zone boundaries for single proxy accumulation-rate records based
on change-point analysis and checks for the influence of
sediment-accumulation rates on the detected change points.
}
\details{
This function is used to determine if any of the change points of
a proxy-accumulation record is strongly influenced by the modeled
sediment-accumulation rates. The function performs a change-point
analysis with the 'true' proxy-accumulation record (based on
the effective concentration values), and in addition it performs a
change-point analysis on a set (by default \code{n_rand = 1000}) of
non-sense proxy-accumulation records (based on randomly drawn
concentration values). The change points of the latter records
(hereafter, non-sense change points) are then first
screened based on their frequency of occurrence (\code{n_screen}, as
some change points may occur by chance). Thereafter, the change points
are screened further to retain
those that separate zones where the mean proxy-accumulation rate
changes in the same direction as the mean sediment-accumulation rate.
The rationale being that if a change point of the 'true
proxy-accumulation record' matches a non-sense change point
(with a tolerance search window equal to \code{tau}), and
if both records change in the same direction (say, both the
proxy-accumulation rate and the sediment-accumulation rate increase
across a change point), then the change point of the
proxy-accumulation rate may be strongly influenced by
sediment-accumulation rates (proxy_ar = proxy_conc * sed_ar).
}
\examples{
\dontrun{
co <- tapas::co_char_data
co_i <- tapas::pretreatment_data(co, out = "accI")

## change-point analysis with default options
co_i_cpts <- cpts_ar(co_i, proxy = "char")

## change-point analysis with higher penalty value
co_i_cpts_8logn <- cpts_ar(co_i, proxy = "char", pen_val = "8*log(n)")

## Change-point analysis with the 'red noise' record (no change point)
rdn <- tapas::red_noise
rdn_i <- tapas::pretreatment_data(rdn)
rdn_i_cpts <- tapas::cpts_ar(rdn_i, proxy = "char")


## As in the example above, but introducing a step-wise increase in the
## sediment-accumulation rate [cm/yr]
## NB: sar = 1/sediment-deposition time [yr/cm]
sdt <- c(rep_len(10, length.out = 100), rep_len(25, length.out = 360))
a_bot <- cumsum(sdt)
rdn2 <- rdn[1:length(a_bot), ]
rdn2$age_bot <- a_bot
rdn2$age_top <- a_bot - sdt
rdn2_i <- tapas::pretreatment_data(rdn2, yrInterp = 25)
rdn2_i_cpts <- tapas::cpts_ar(rdn2_i, proxy = "char")
}

}
\references{
Finsinger W, Magyari EK, Fevre J, Orban I, Pal I, Vincze I,
Hubay K, Birks HH, Braun M, Toth M (2016) Holocene fire regimes near the
treeline in the Retezat Mts. (Southern Carpathians). Quaternary
International, \emph{477}, 94-105. 10.1016/j.quaint.2016.04.029.

Killick R, Eckley IA (2014). changepoint: An R Package for Changepoint
Analysis. Journal of Statistical Software, \emph{58}(3), 1-19.
}
\seealso{
\code{\link{cpt.meanvar}}, \code{\link{pretreatment_data}}
}
\author{
Walter Finsinger
}
