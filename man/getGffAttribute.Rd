\name{getGffAttribute}

\alias{getGffAttribute}

\title{Pull one or more key/value pairs from gffAttributes strings}

\description{
  GFF files contain a string, with key/value pairs separated by
  \dQuote{;}, and the key and value separated by \dQuote{=}. This
  function quickly extracts one or more key/value pairs.
}

\usage{
getGffAttribute(gi, attribute)
}

\arguments{
  \item{gi}{A \code{\link[=Genome_intervals-class]{Genome_intervals}} object.}
  \item{attribute}{A vector of key names.}
}

\value{
  A matrix with the same number of rows as \code{gi}, and one column per
  element of \code{attribute}.
}

\seealso{
  See \code{\link{parseGffAttributes}} for more complete parsing. See the function \code{\link{readGff3}} for loading a GFF file.
}

\examples{

# Get file path
libPath <- installed.packages()["genomeIntervals", "LibPath"]
filePath <- file.path(
                   libPath,
                   "genomeIntervals",
                   "example_files"
                   )

# Load gff
gff <- readGff3( file.path( filePath, "sgd_simple.gff"), isRightOpen=FALSE)

## head of full gff annotations
head(annotation(gff))

# extract ID and Parent attributes
idpa = getGffAttribute( gff, c( "ID", "Parent" ) )

head(idpa)


}
