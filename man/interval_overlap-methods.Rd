\name{interval_overlap}

\alias{interval_overlap}
\alias{interval_overlap,Genome_intervals,Genome_intervals-method}
\alias{interval_overlap,Genome_intervals_stranded,Genome_intervals_stranded-method}

\title{Assess which query intervals overlap which targets}

\description{
  Given two objects, a target and a query, assess which query intervals
  overlap which targets.
}

\usage{
\S4method{interval_overlap}{Genome_intervals,Genome_intervals}(
                 target, query,
                 check_valid = TRUE
                 )
\S4method{interval_overlap}{Genome_intervals_stranded,Genome_intervals_stranded}(
                 target, query,
                 check_valid = TRUE
                 )

}

\arguments{
  
  \item{target}{
    A \code{Genome_intervals} or \code{Genome_intervals_stranded} object.
    } 

  \item{query}{A \code{Genome_intervals} or \code{Genome_intervals_stranded} object.}

  
  \item{check\_valid}{
    Should \code{\link{validObject}} be called before passing to
    compiled code?
    }
  
}

\details{
  A wrapper calling \code{\link[intervals:interval_overlap-methods]{intervals:interval_overlap}} by \code{seq_name} and by \code{strand} (if both \code{query} and \code{target} are \code{"Genome_intervals_stranded"} objects).
}

\value{
  A list, with one element for each row of \code{target}. The elements
  are vectors of indices, indicating which \code{query} rows overlap
  each target. A list element of length 0 indicates a target with no
  overlapping query intervals.
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
