\name{interval_union}

\alias{interval_union}
\alias{interval_union,Genome_intervals-method}
\alias{interval_union,Genome_intervals_stranded-method}

\alias{interval_complement}
\alias{interval_complement,Genome_intervals-method}
\alias{interval_complement,Genome_intervals_stranded-method}

\alias{interval_intersection}
\alias{interval_intersection,Genome_intervals-method}
\alias{interval_intersection,Genome_intervals_stranded-method}


\title{Genome interval set operations}

\description{
  Compute interval set operations on \code{"Genome_intervals"} or \code{"Genome_intervals_stranded"} objects.
}

\usage{
\S4method{interval_union}{Genome_intervals}(x, ...)
\S4method{interval_union}{Genome_intervals_stranded}(x, ...)

\S4method{interval_complement}{Genome_intervals}(x)
\S4method{interval_complement}{Genome_intervals_stranded}(x)

\S4method{interval_intersection}{Genome_intervals}(x,...)
\S4method{interval_intersection}{Genome_intervals_stranded}(x,...)
}

\arguments{
  \item{x}{ A \code{"Genome_intervals"} or \code{"Genome_intervals_stranded"} object.}
  \item{\dots}{Optionally, additional objects of the same class as \code{x}.}
}

\value{
   A single object of appropriate class, representing the union, complement or
   intersection of intervals computed over entries with same \code{seq_name}, \code{inter_base} and 
   also \code{strand} if all passed objects are of the class \code{"Genome_intervals_stranded"}.
}

\details{
	Wrappers calling the corresponding functions of the package \code{intervals} by 
	same \code{seq_name}, \code{inter_base} and if needed \code{strand}. Note that the   
	union of single input object \code{x} returns the reduced form of \code{x}, i.e.
	the interval representation of the covered set.
}

\seealso{
	\code{\link[intervals:interval_union-methods]{interval_union}},
	\code{\link[intervals:interval_complement-methods]{interval_complement}},
	\code{\link[intervals:interval_intersection-methods]{interval_intersection}} and
	\code{\link[intervals:reduce-methods]{reduce}} from the package \code{intervals}.
}

\examples{
## load toy examples
data(gen_ints)
## content of i object
i
	 
## complement
interval_complement(i)
 
## reduced form (non-overlapping interval representation of the covered set)
interval_union(i)

## union
interval_union(i[1:2,], i[1:4,])
 
# map to genome intervals and union again
i.nostrand = as(i,"Genome_intervals")
interval_union(i.nostrand)
 
## intersection with a second object
# print i and j in closed interval notation
close_intervals(i)
close_intervals(j)

# interval_intersection
interval_intersection(i,j)

#interval intersection non-stranded
interval_intersection(i.nostrand, as(j, "Genome_intervals"))
}
