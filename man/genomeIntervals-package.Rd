\name{genomeIntervals-package}

\alias{genomeIntervals-package}
\alias{genomeIntervals}

\docType{package}

\title{
Operations on genomic intervals
}

\description{
Tools for operation on genomic intervals.
}

\details{
	\tabular{ll}{
	Package: \tab genomeIntervals\cr
	Version: \tab 0.9.6\cr
	Date: \tab 2009-01-15\cr
	Type: \tab Package\cr
	Depends: \tab R (>= 2.8.0), intervals (>= 0.10.3), Biobase, methods\cr
	Suggests: \tab \cr
	License: \tab Artistic 2.0\cr
	BiocViews: \tab DataImport, Infrastructure, Genetics\cr
	LazyLoad: \tab yes\cr
	}

	Index:
	\describe{
	 	\item{\code{\linkS4class{Genome_intervals}}}{Class "Genome\_intervals"}
		\item{\code{\linkS4class{Genome_intervals_stranded}}}{Class "Genome\_intervals\_stranded"}
		
		\item{\code{\link[=distance_to_nearest,Genome_intervals,Genome_intervals-method]{distance_to_nearest}}}{Distance in bases to the closest interval(s)}
		
		\item{\code{\link{gen_ints}}}{Genome Intervals examples}
		\item{\code{\link{getGffAttribute}}}{Pull one or more key/value pairs from gffAttributes strings}
		
		\item{\code{\link[=interval_overlap,Genome_intervals,Genome_intervals-method]{interval_overlap}}}{Assess overlap from one set of genomic intervals to another}
		
		\item{\code{\link[=interval_complement,Genome_intervals-method]{interval_complement}}}{Compute the complement of a set of genomic intervals}
    	
    	\item{\code{\link[=interval_intersection,Genome_intervals-method]{interval_intersection}}}{Compute the intersection of one or more sets of genomic intervals}
    	
    	\item{\code{\link[=interval_union,Genome_intervals-method]{interval_union}}}{Compute the union of genomic intervals in one or more genomic interval matrices}
    	
    	\item{\code{\link{parseGffAttributes}}}{Parse out the gffAttributes column of a \code{Genome_intervals} object}
		
		\item{\code{\link{readGff3}}}{Make a \code{Genome_intervals_stranded} object from a GFF file}
		}
}
\author{
Julien Gagneur <gagneur@embl.de>, Richard Bourgon.

Maintainer: Julien Gagneur <gagneur@embl.de>
}
\keyword{ package }
\seealso{
\code{\link[intervals:intervals-package]{intervals}}
}
