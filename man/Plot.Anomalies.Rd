% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotAnomalies.R
\name{Plot.Anomalies}
\alias{Plot.Anomalies}
\title{Plots a summary figure of results obtained with 'tapas'.}
\usage{
Plot.Anomalies(
  series = NULL,
  x.lim = NULL,
  y_lim = NULL,
  plot.crosses = T,
  plot.x = F,
  plot.neg = T
)
}
\arguments{
\item{series}{The output of \code{\link{global_thresh}()},
\code{\link{local_thresh}()}, or \code{\link{peak_detection}()}.}

\item{x.lim}{A vector defining the limits for the x-axis scale (time scale).}

\item{y_lim}{A vector defining the limits for the y-axis scale.}

\item{plot.crosses}{Boolean. If \code{TRUE} (by default),
crosses are displayed to indicate
the location of anomalies.}

\item{plot.x}{Boolean. If \code{FALSE} (by default),
the x-axis labels are omitted.}

\item{plot.neg}{Boolean. If \code{FALSE} (by default),
both positive and negative anomalies
are marked with colored shaded areas.}
}
\description{
Requires output from one of the following functions:
\code{\link{global_thresh}()}, \code{\link{local_thresh}()},
or \code{\link{peak_detection}()}.
}
\author{
Walter Finsinger
}
