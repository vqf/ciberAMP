% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnaintexp.r
\name{int.plot.CiberAMP}
\alias{int.plot.CiberAMP}
\title{Interactive plot with ShinyR package}
\usage{
int.plot.CiberAMP(df, int.df)
}
\arguments{
\item{df}{First or second data frame from CNAintEXP list of results. Contains the correlation results for SCNA and mRNA differential expression for 1) queried or 2) COSMIC CGC genes.}

\item{int.df}{Third data frame from CNAintEXP list of results. Contains overlappings between SCN-altered queried and COSMIC CGC genes.}
}
\value{
It allows the user to directly interact with data using a shiny app
}
\description{
Interactive plot with ShinyR package
}
