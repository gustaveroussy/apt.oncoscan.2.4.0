\name{apt.oncoscan.process}
\alias{apt.oncoscan.process}
\title{Normalize OncoScan arrays using Affymetrix Power Tools.}
\description{
  This function performs the normalization of CEL file list for Affymetrix OncoScan microarrays using the binary apt-copynumber-onco-ssa (v2.4.0 2015-08-28) from Affymetrix Power Tools.
}
\usage{
apt.oncoscan.process(ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL,
dual.norm = TRUE, out.dir = getwd(), temp.files.keep = FALSE, force.OS = NULL,
apt.build = "na33.r2")
}
\arguments{
\item{ATChannelCel}{The first CEL file (version A).}
\item{GCChannelCel}{The second CEL file (version C).}
\item{samplename}{The name of the sample.}
\item{dual.norm}{Perform a dual-channel normalization (recommended for a better wave-effect bias normalization).}
\item{out.dir}{The path to output results.}
\item{temp.files.keep}{Keep the (numerous) intermediate files.}
\item{force.OS}{Specifiy the OS type of the execution machine (see Details).}
\item{apt.build}{Affymetrix genome build version.}
}
\details{
  Output files (all begin with the samplename as a prefix) :
  \itemize{
    \item{'*.OSCHP' : Contains the processed data. This format actually internaly is HDF5.}
  }
  \code{force.OS} is \code{NULL} by default, but can take any of these values :
  \itemize{
    \item{\code{linux}}
    \item{\code{windows}}
    \item{\code{osx}}
  }
}
\author{
Bastien Job
}

