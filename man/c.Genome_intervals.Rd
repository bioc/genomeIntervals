% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c.R
\name{c}
\alias{c}
\alias{c,Genome_intervals-method}
\alias{c,Genome_intervals_stranded-method}
\title{c extension for the \code{genomeIntervals package}}
\usage{
\S4method{c}{Genome_intervals}(x, ..., recursive = FALSE)
}
\arguments{
\item{x}{a \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} object - not mandatory.}

\item{...}{two (one if \code{x} is defined) or more \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} objects.}

\item{recursive}{inherited from the base \code{c} function definition and not used.}
}
\value{
\itemize{
  \item A single \code{\linkS4class{Genome_intervals}} or
  \code{\linkS4class{Genome_intervals_stranded}} object. Input objects are
  combined in their order of appearance in the the argument list.
  \item If any input argument is not a \code{\linkS4class{Genome_intervals}}, \code{list(...)} is
  returned instead.
}##'
}
\description{
This function combines several genome intervals (stranded or not)
objects into a single one.
}
\details{
If the arguments have mixed classes ( both \code{\linkS4class{Genome_intervals}}
or \code{\linkS4class{Genome_intervals_stranded}}), then
they are coerced to \code{\linkS4class{Genome_intervals}} before combination.
Otherwise, the common class is used. Ifa list is provided with \code{NULL}
entries, these are discarded. If a vector of object is provided with non \code{genomeIntervals}
classes, then a list, ordered as the input vector, is returned.
}
\examples{
  ##' load toy examples
  data("gen_ints")

  ##' combine i and j returns a Genome_intervals_stranded object
  c( i, j )

  ##' combine a not-stranded and a stranded returns a not-stranded object
  c( as(i, "Genome_intervals"), j )

}
