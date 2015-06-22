\name{Genome_intervals-ordering}

\alias{Genome_intervals-ordering}
\alias{order}
\alias{rank}
\alias{sort}
\alias{xtfrm}

\alias{order,Genome_intervals-method}
\alias{order,Genome_intervals_stranded-method}

\alias{rank,Genome_intervals-method}
\alias{rank,Genome_intervals_stranded-method}

\alias{sort,Genome_intervals-method}
\alias{sort,Genome_intervals_stranded-method}

\alias{xtfrm,Genome_intervals-method}
\alias{xtfrm,Genome_intervals_stranded-method}

\title{Ordering methods for Genome intervals}

\description{
  An order is defined on genome intervals and stranded genome intervals to allow sort(), order() and rank().
 }


\usage{

\S4method{order}{Genome_intervals}(..., na.last=TRUE, decreasing=FALSE)
\S4method{order}{Genome_intervals_stranded}(..., na.last=TRUE, decreasing=FALSE)

\S4method{sort}{Genome_intervals}(x, decreasing=FALSE, ...)

\S4method{rank}{Genome_intervals}(x, na.last=TRUE,
     ties.method = c("average", "first", "random", "max", "min") )
\S4method{xtfrm}{Genome_intervals}(x)
     
}



\arguments{
	\item{x}{
    Objects of class \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}}.
  }
  
  \item{...}{
    Objects of class \code{\linkS4class{Genome_intervals}}, \code{\linkS4class{Genome_intervals_stranded}} or of any other class for \code{order}.
  }
  \item{na.last}{
    Ignored for ordering  \code{\linkS4class{Genome_intervals}} and  \code{\linkS4class{Genome_intervals_stranded}} objects 
  }
  \item{decreasing}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{ties.method}{
    A character string specifying how ties are treated. Only \code{"first"}
    is supported.
  }
}

\details{
   An order on Genome_intervals entries is defined by sorting by 
	1. seq_name
	2. start, where closed start & not inter-base < closed start & inter-base < open start & not inter-base < open start & inter-base
	3. stop, where open stop & not inter-base < open stop & inter-base < closed stop & not inter-base < closed stop & inter-base 
	4. strand (for Genome_intervals_stranded object)
  
  The factors \code{seq_name} and \code{strand} are sorted according to their levels (default R behavior).
  
  The primitive is implemented in \code{xtfrm} which is then called by the other methods.
  Hence, the \code{order}, \code{sort} and \code{rank} methods are consistent.
  
  \describe{
    \item{}{
      \code{order(..., na.last=TRUE, decreasing=TRUE)}:
      return a permutation which rearranges its first argument
     into ascending or descending order, breaking ties by further
     arguments.
        See \code{\link[base]{order}} in the base package for more details.
        \code{na.last} is ignored for \code{\linkS4class{Genome_intervals}} objects. 
    }
    \item{}{
      \code{rank(x, na.last=TRUE, ties.method=c("average", "first", "random", "max", "min"))}:
      Return the sample ranks of the (stranded) genome intervals in \code{x}.
      See \code{\link[base]{rank}} in the base package for more details.
    }
    \item{}{
      \code{sort(x)}:
      Sort \code{x}.
      See \code{\link[base]{sort}} in the base package for more details.
    }
     \item{}{
      \code{xtfrm(x)}:
      Auxiliary function that produces a numeric vector which will sort in the same order as 'x' \code{x}.
      See \code{\link[base]{xtfrm}} in the base package for more details.
      Workhorse for the other methods
    }
	
	}    
}

\seealso{
	\code{\linkS4class{Genome_intervals}}
	\code{\linkS4class{Genome_intervals_stranded}}
  \code{\link[BiocGenerics]{order}},
  \code{\link[base]{sort}},
  \code{\link[base]{rank}},
  \code{\link[base]{xtfrm}}
}

\examples{
## an example with ties
gi = GenomeIntervals(c("chr2", "chr2", "chr1", "chr1"), c(1,1,10,10), c(5,3,12,12) )

sort(gi)
rank(gi)
order(gi)

## Define order on seq_names at your convenience
## by specifying the order of the levels
## compare: 
gi = GenomeIntervals(
		c("chr2", "chr2", "chr10", "chr10"),
		c(1,1,10,10),
		c(5,3,12,12)
		)
sort(gi)

## with:
gi2 = GenomeIntervals(
	factor(c("chr2", "chr2", "chr10", "chr10"), levels=c("chr2", "chr10")),
	c(1,1,10,10),
	c(5,3,12,12)
	)
sort(gi2)

}
