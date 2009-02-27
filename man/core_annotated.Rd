\name{core_annotated}
\alias{core_annotated}
\alias{core_annotated,Genome_intervals-method}
\alias{core_annotated,Genome_intervals_stranded-method}

\title{ Genome intervals with minimal annotation }
\description{
  returns a copy of the input (stranded) genome intervals object with annotations restricted to the minimally required ones.
}
\usage{
core_annotated(x)
}

\arguments{
  \item{x}{ A \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} object. }
}

\value{
  A copy of \code{x} with the \code{annotation} slot restricted to \code{seq_name}, \code{inter_base}
  and \code{strand} (the latter only if \code{x} is a \code{\linkS4class{Genome_intervals_stranded}} object).
}

\examples{
# load toy examples
data("gen_ints")

# add some non-core annotations to i
annotation(i)$comment = "some non-core annotation"

# i with all annotations
i

# core annotations only
core_annotated(i)

\dontrun{
# with different annotation columns, i and j cannot be combined
c( i, j )
}

# core annotated versions can
c( core_annotated(i), core_annotated(j) )

}
