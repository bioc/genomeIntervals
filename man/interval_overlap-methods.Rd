\name{interval_overlap}

\alias{interval_overlap}
\alias{interval_overlap,Genome_intervals,Genome_intervals-method}
\alias{interval_overlap,Genome_intervals_stranded,Genome_intervals_stranded-method}
\alias{interval_overlap,missing,ANY-method}
\alias{interval_overlap,ANY,missing-method}

\title{Assess overlap from one set of genomic intervals to another}

\description{
  Given two objects, a 'from' and a 'to', assess which intervals in 'to'
  overlap which of 'from'.
}

\usage{
\S4method{interval_overlap}{Genome_intervals,Genome_intervals}(
                 from, to,
                 check_valid = TRUE
                 )
\S4method{interval_overlap}{Genome_intervals_stranded,Genome_intervals_stranded}(
                 from, to,
                 check_valid = TRUE
                 )

}

\arguments{
  
  \item{from}{
    A \code{Genome_intervals} or \code{Genome_intervals_stranded} object.
    } 

  \item{to}{A \code{Genome_intervals} or \code{Genome_intervals_stranded} object.}

  
  \item{check\_valid}{
    Should \code{\link{validObject}} be called before passing to
    compiled code?
    }
  
}

\details{
  A wrapper calling \code{\link[intervals:interval_overlap-methods]{intervals:interval_overlap}} by \code{seq_name} and by \code{strand} (if both \code{to} and \code{from} are \code{"Genome_intervals_stranded"} objects).
}

\value{
  A list, with one element for each row of \code{from}. The elements
  are vectors of indices, indicating which \code{to} rows overlap
  each from. A list element of length 0 indicates a from with no
  overlapping to intervals.
}


\examples{
data(gen_ints)	 
# i as entered
i

# i in close_intervals notation
close_intervals(i)
 
# j in close_intervals notation
close_intervals(j)
 
# list of intervals of j overlapping intervals of i
interval_overlap(i,j)
	 
}
