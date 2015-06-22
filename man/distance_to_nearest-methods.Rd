\name{distance_to_nearest}

\alias{distance_to_nearest}
\alias{distance_to_nearest,Genome_intervals,Genome_intervals-method}
\alias{distance_to_nearest,Genome_intervals_stranded,Genome_intervals_stranded-method}

\title{Distance in bases to the closest interval(s)}

\description{
  Given two objects, \code{from} and \code{to}, compute the distance in bases of each \code{from} interval
  to the nearest \code{to} interval(s).
  The distance between a base and the next inter-bases on either side values 0.5. Thus,
  base - base and inter-base - inter-base intervals distances are integer,
  whereas base - inter-base intervals distances are half-integers. 
}

\usage{
\S4method{distance_to_nearest}{Genome_intervals,Genome_intervals}(from, to)
\S4method{distance_to_nearest}{Genome_intervals_stranded,Genome_intervals_stranded}(from, to)
}

\arguments{
  
  \item{from}{A \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} object.}

  \item{to}{
    A \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals} } object.
  }
  
}

\details{
  A wrapper calling \code{\link[intervals:distance_to_nearest-methods]{intervals::distance_to_nearest}} by \code{seq_name} and by \code{strand}
  (if both \code{from} and \code{to} are \code{Genome_intervals_stranded} objects).
  Thus, if both are stranded, distances are computed over each strand separately.
  One object must be coerced to \code{\linkS4class{Genome_intervals}} if this is not wished. 
}

\value{
  A numeric vector of distances with one element for each row of \code{from}.
}

\seealso{\code{\link[intervals:distance_to_nearest-methods]{intervals::distance_to_nearest}}}

\examples{
## load toy examples
data(gen_ints)		 
 
## i in close_intervals notation
close_intervals(i)
 
## j in close_intervals notation
close_intervals(j)
 
## distances from i to j 
dn = distance_to_nearest(i,j)
dn

## distance == 0 if and only if the interval overlaps another one:
io = interval_overlap(i,j)
if( any( ( sapply(io, length) >0 )  != (!is.na(dn) & dn ==0) ) )
 	stop("The property 'distance == 0 if and only if the interval overlaps another one' is not followed for at least one instance.")
 
## distances without strand-specificity
distance_to_nearest( 
 		as(i,"Genome_intervals"),
 		as(j,"Genome_intervals") 
) 	
}
