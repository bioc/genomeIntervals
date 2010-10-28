\name{readGff3}
\alias{readGff3}

\title{Make a Genome\_intervals\_stranded object from a GFF file}

\description{
  Make a \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object from a gff file in gff3 format.
}

\usage{
readGff3(file, isRightOpen=TRUE)
}

\arguments{
  \item{file}{ The name of the gff file to read. }
  \item{isRightOpen}{
    Although a proper GFF3 file follows the convention of right-open intervals, 
    improper GFF files following the right-closed convention are frequently found.
    Set \code{isRightOpen = FALSE} in this case.
  }
}

\details{
  The file must follow gff3 format specifications as in
  \url{http://www.sequenceontology.org/gff3.shtml}.  The file is read as a
  table. Meta-information (lines starting with \#\#\#) are not parsed. A
  \dQuote{.} in, for example, the gff file's \emph{score} or
  \emph{frame} field will be converted to \code{NA}.
  When the GFF file follows the right-open interval convention (\code{isRightOpen} is \code{TRUE}), then
  GFF entries for which end base equals first base are recognized as zero-length features
  and loaded as \code{inter_base} intervals.
  Strand entries in the file are expected to be '.', '?', '+' or '-'. The two first are mapped to \code{NA}. 
  It can be that readGff3 is able to construct a \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object from the input file, although not valid.
A warning message is then generated and the constructed object is returned to allow inspection of it.
}

\value{
  A \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object image of the
  gff file.  The GFF3 fields \code{seqid}, \code{source}, \code{type}, \code{score}, \code{strand}, \code{phase} and
  \code{attributes} are stored in the \code{annotation} slot and renamed as \code{seq_name}, \code{source},
  \code{type}, \code{score}, \code{strand}, \code{phase} and \code{gffAttributes} respectively.
}


\note{
Potential FASTA entries at the end of the file are ignored.
}
\seealso{ The functions \code{\link{getGffAttribute}} and \code{\link{parseGffAttributes}} for parsing GFF attributes. }
\examples{

# Get file path
libPath <- installed.packages()["genomeIntervals", "LibPath"]
filePath <- file.path(
                      libPath,
                      "genomeIntervals",
                      "example_files"
                      )

# Load SGD gff
# SGD does not comply to the GFF3 right-open interval convention
gff <- readGff3( file.path( filePath, "sgd_simple.gff"), isRightOpen = FALSE)

head(gff,10)

head(annotation(gff),10)

}
