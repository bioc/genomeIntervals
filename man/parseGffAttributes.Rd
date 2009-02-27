\name{parseGffAttributes}

\alias{parseGffAttributes}

\title{Parse out the gffAttributes column of a Genome\_intervals object}

\description{
  GFF files contain a string, with key/value pairs separated by
  \dQuote{;}, and the key and value separated by \dQuote{=}. This
  function parses such strings into a list of vectors with named
  elements.
}

\usage{
parseGffAttributes(gi)
}

\arguments{
  \item{gi}{A \code{\link[=Genome_intervals-class]{Genome_intervals}} object.}
}

\value{
  A list, with one element per row of \code{gi}. Each element is a
  character vector with named components. Names correspond to keys, and
  components correspond to values.
}

\seealso{
  In many cases, \code{\link{getGffAttribute}}, in this package, is
  easier and faster. See the function \code{\link{readGff3}} for loading a GFF file. 
}

\note{
  Key/value pairs which are missing the \dQuote{=} symbol, or which have
  nothing between it and the \dQuote{;} delimiter or end of line, will
  generate a \code{NA} value, with a warning. Any key/value
  \dQuote{pairs} with more than one \dQuote{=} cause an error.
}

\examples{

# Get file path
libPath <- installed.packages()["genomeIntervals", "LibPath"]
filePath <- file.path(
                   libPath,
                   "genomeIntervals",
                   "example_files"
                   )

# Load gff and parse attributes
gff <- readGff3( file.path( filePath, "sgd_simple.gff"), isRightOpen = FALSE )
gfatt <- parseGffAttributes(gff)

head( gfatt )
}
