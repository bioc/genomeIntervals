\name{GenomeIntervals}
\alias{GenomeIntervals}
\alias{GenomeIntervals-constructor}
\title{Constructor function for genomeIntervals objects}
\description{
  A user-friendly constructor function for creating both
  \code{Genome_intervals} and \code{Genome_intervals_stranded}
  objects.
}
\usage{
GenomeIntervals(chromosome, start, end, strand = NULL,
                inter.base = NULL, leftOpen = NULL,
                rightOpen = NULL, ...)
}
\arguments{
  \item{chromosome}{character vector of chromosome names of the
    intervals; will become the \code{seqnames} of the resulting object}
  \item{start}{numeric or integer; start (left-most) coordinate of the
    intervals}
  \item{end}{numeric or integer; end (right-most) coordinate of the
    intervals}
  \item{strand}{chacter; specifies which strand the intervals are
    located on; if specified an object of class
    \code{Genome_intervals_stranded} is created; if \code{NULL} an
    object of class \code{Genome_intervals} is created}
  \item{inter.base}{logical; if \code{TRUE} an interval is located
    between the specified coordinates, instead of spanning them; useful
    for restriction-enzym cutting sites, for example.}
  \item{leftOpen}{logical; if \code{TRUE} an interval is left-open;
    if \code{NULL} all intervals are assumed to be left-closed.}
  \item{rightOpen}{logical; if \code{TRUE} an interval is right-open;
    if \code{NULL} all intervals are assumed to be right-closed.}
  \item{\dots}{any additional annotation for supplied intervals}
}
\details{
  The arguments \code{chromosome}, \code{start}, and \code{end} need to
  be of the same length, with the first element of each vector
  corresponding to the first interval, the second element to the second
  interval, and so on.

  The same applies to \code{strand}, \code{inter.base}, \code{leftOpen},
  \code{rightOpen} and any additional vectors in '...', if they are
  specified.
}
\value{
  An object of class \code{Genome_intervals} or
  \code{Genome_intervals_stranded} depending on whether \code{strand}
  has been specified.
}
\author{J. Toedling}
\seealso{
  \code{\link{Genome_intervals-class}},
  \code{\link{Genome_intervals_stranded-class}}
}
\examples{
  ## constructing a Genome_intervals object
  G <- GenomeIntervals(start=c(1,3,4,5,8,10), end=c(5,5,6,8,9,11),
                       chromosome=rep(c("chr2","chrX","chr1"), each=2),
                       leftOpen=rep(c(FALSE, FALSE, TRUE), 2))
  show(G)

  ## constructing a Genome_intervals_stranded object with
  ##  additional interval annotation
  GS <- GenomeIntervals(start=c(1,3,4,5,8,10), end=c(5,5,6,8,9,11),
                        chromosome=rep(c("chr2","chrX","chr1"), each=2),
                        strand=c("-","-","+","+","+","+"),
                        GC.content=round(runif(6), digits=2))
  show(GS)
}
\keyword{manip}
