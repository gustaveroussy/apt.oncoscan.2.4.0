\name{apt.oncoscan.process.batch}
\alias{apt.oncoscan.process.batch}
\title{Runs apt.oncoscan.process() in batch mode, with multithreading}
\description{
  Same as apt.oncoscan.process() in batch mode, with multithreading.
}

\usage{
  apt.oncoscan.process.batch(pairs.file = NULL, nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{pairs.file}{A tab-separated file giving the file names of the pair(s) of OncoScan CEL files and sample name.}
  \item{nthread}{Number of thread(s) for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any other parameter for \code{apt.oncoscan.process}.}
}
\details{
  \code{pairs.file} is a tab-separated text file containing 3 columns :
  \itemize{
    \item{ATChannelCel : Name (and path) of the A-channel CEL file(s)}
    \item{GCChannelCel : Name (and path) of the C-channel CEL file(s)}
    \item{SampleName : The output sample name(s)}
  }
}
\author{
Bastien Job
}

