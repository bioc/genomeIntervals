\name{c}

\alias{c}
\alias{c.Genome_intervals}
\alias{c.Genome_intervals_stranded}

\title{Combine genome intervals objects}

\description{
  S3 methods for combining several genome intervals into a single one.
}

\usage{
\S3method{c}{Genome_intervals}(...)
\S3method{c}{Genome_intervals_stranded}(...)
}

\arguments{
  
  \item{...}{\code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} objects.}
  
}

\details{
  If the arguments have mixed classes ( both \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}}), then
  they are coerced to \code{\linkS4class{Genome_intervals}} before combination. Otherwise, the common class is used.  
}

\value{
  A single \code{\linkS4class{Genome_intervals}} or
  \code{\linkS4class{Genome_intervals_stranded}} object. Input objects are
  combined in their order of appearance in the the argument list.
  
  If any input argument is not a \code{\linkS4class{Genome_intervals}}, \code{list(...)} is
  returned instead.
}

\note{
  These methods will be converted to S4 once the necessary dispatch on
  \code{...} is supported.
}

\examples{
# load toy examples
data("gen_ints")
    
# combine i and j returns a Genome_intervals_stranded object
c( i, j )
    
# combine a not-stranded and a stranded returns a not-stranded object
c( as(i, "Genome_intervals"), j )
}