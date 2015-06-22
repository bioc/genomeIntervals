##' readGff3
##'
##' Read (write) a Genome_intervals_stranded object from (to) a GFF3 file
##'
##' \itemize{
##'   \item{readGff3}{
##'     Make a \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object from a gff file in gff3 format.}
##'   \item{readBasePairFeaturesGff3}{
##'     Same as \code{readGff3} assuming \code{isRightOpen='FALSE'}, i.e. no zero length intervals are created. This
##'     is the default behaviour since v1.25.1.
##'   }
##'   \item{readZeroLengthFeaturesGff3}{
##'     Same as \code{readGff3} assuming \code{isRightOpen='TRUE'}, i.e. zero length intervals are created when
##'     a feature's start is the same as its end. This was the default prior to version 1.25.1.
##'   }
##'   \item{writeGff3}{
##'     Write a \code{\link[=Genome_intervals-class]{Genome_intervals}} object to a gff file in gff3 format.
##'   }}
##'
##' The file must follow gff3 format specifications as in
##' \url{http://www.sequenceontology.org/gff3.shtml}. Due to the imprecise definition and to
##' allow for zero-length features, the default for reading a Gff3 file has been to assume
##' right open intervals (until v1.25.1). As by then, the community consensus has been to use
##' closed intervals, the default behaviour of readGff3 has been changed accordingly. The
##' \code{readGff3 file} is now a wrapper that dispatches to two sub functions - which may be used
##' directly - \code{readBasePairFeaturesGff3} and \code{readZeroLengthFeaturesGff3}. The former
##' assumes closed intervals and hence does not create zero-length intervals. The latter does the
##' opposite and uses right-open intervals!
##'
##' Some more noteworthy details:
##'
##' The file is read as a table and meta-information (lines starting with ###) are not parsed.
##'
##' A \dQuote{.} in, for example, the gff file's \emph{score} or \emph{frame}
##' field will be converted to \code{NA}.
##'
##' When the GFF file follows the right-open interval convention (\code{isRightOpen} is \code{TRUE}), then
##' GFF entries for which end base equals first base are recognized as zero-length features
##' and loaded as \code{inter_base} intervals.
##'
##' Strand entries in the file are expected to be '.', '?', '+' or '-'. The two first are mapped to \code{NA}.
##'
##' It can be that readGff3 is able to construct a \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object from the input file, although not valid.
##' A warning message is then generated and the constructed object is returned to allow inspection of it.
##'
##' Potential FASTA entries at the end of the file are ignored.
##'
##' @aliases readGff3 readGff3,character-method readZeroLengthFeaturesGff3
##' readZeroLengthFeaturesGff3,character-method readBasePairFeaturesGff3
##' readBasePairFeaturesGff3,character-method writeGff3 writeGff3,data.frame-method
##' writeGff3,Genome_intervals-method .readGff3
##' @rdname genomeIntervals-readGff3
##' @param file The name of the gff file to read/write.
##' @param isRightOpen Although it is arguable that a GFF3 file might have a right-open intervals convention - the format
##' description being at best imprecise - most GFF3 file follow a right-closed convention. Hence, as of version 1.25.1,
##' the default has been changed to \code{isRightOpen = FALSE}. See the details section on how to restore the older
##' behaviour.
##' @param object a \code{\link[=Genome_intervals-class]{Genome_intervals}} object
##' @param quiet a boolean to turn verbosity off when reading a Gff3 file
##' @return \itemize{
##'   \item{readGff3 and friends}{A \code{\link[=Genome_intervals_stranded-class]{Genome_intervals_stranded}} object image of the
##' gff file. The GFF3 fields \code{seqid}, \code{source}, \code{type}, \code{score}, \code{strand}, \code{phase} and
##' \code{attributes} are stored in the \code{annotation} slot and renamed as \code{seq_name}, \code{source},
##' \code{type}, \code{score}, \code{strand}, \code{phase} and \code{gffAttributes} respectively.}
##'   \item{writeGff3}{It dispatches to \code{write.table} and hence returns similar values.}
##' }
##' @seealso The functions \code{\link{getGffAttribute}} and \code{\link{parseGffAttributes}} for parsing GFF attributes.
##' @usage
##' readGff3(file, isRightOpen=FALSE, quiet=FALSE)
##' readBasePairFeaturesGff3(file, quiet=FALSE)
##' readZeroLengthFeaturesGff3(file, quiet=FALSE)
##' writeGff3(object, file)
##' @examples
##' # Get file path
##' libPath <- installed.packages()["genomeIntervals", "LibPath"]
##' filePath <- file.path(
##'  libPath,
##'  "genomeIntervals",
##'  "example_files"
##' )
##'
##' # Load SGD gff
##' # SGD does not comply to the GFF3 right-open interval convention
##' gff <- readGff3( file.path( filePath, "sgd_simple.gff"), isRightOpen = FALSE)
##'
##' head(gff,10)
##'
##' head(annotation(gff),10)
##'
##'\dontrun{
##' ## write the gff3 file
##' writeGff3(gff,file="sgd_simple.gff")
##' }
##'
setMethod(f = "readGff3",
          signature = "character",
          definition = function(file=character(1),
                                isRightOpen=FALSE,
                                quiet=FALSE){
            if(!quiet){
              warning(paste("'readGff3' has changed to closed interval conventions!",
                            "Use 'isRightOpen=TRUE' to restore the previous behavior",
                            "that allowed for zero-length features. Alternatively, use",
                            "the readZeroLengthFeaturesGff3 function instead.",
                            "You can turn off this warning by setting 'quiet=TRUE'",sep="\n"))
            }

            # dispatch
            return(switch(as.character(isRightOpen),
                          "TRUE"=readZeroLengthFeaturesGff3(file=file,quiet=quiet),
                          "FALSE"=readBasePairFeaturesGff3(file=file,quiet=quiet)
            ))
          })

setMethod(f="readZeroLengthFeaturesGff3",
          signature="character",
          definition=function(file=character(1),
                              quiet=FALSE){
            .readGff3(file=file,isRightOpen=TRUE,quiet=quiet)
          })

setMethod(f="readBasePairFeaturesGff3",
          signature="character",
          definition=function(file=character(1),
                              quiet=FALSE){
            .readGff3(file=file,isRightOpen=FALSE,quiet=quiet)
          })

## Internal function
.readGff3 <- function( file, isRightOpen = FALSE, quiet = FALSE ) {

    # A sanity check
    stopifnot(file.exists(file))

    # This is to ignore potential fasta sequenced entries at the end of the gff filw
    # Find total number of lines before first "##FASTA" line, for use in scan()
    # call that follows. Note that scan() is used for speed, and because nlines
    # counts both comment and non-comment lines. The nrows argument to
    # read.table() will not include comment lines, causing problems.

    l <- readLines( file )
    nlines <- -1
    if ( any(l == "##FASTA") ) nlines <- which(l == "##FASTA")[1] - 1
    rm( l ) # Free memory

    ## load data frame
    gff <- scan(
            file,
            nlines = nlines,
            sep = "\t",
            comment.char = "#",
            na.strings = ".",
            quote = "",
            what = list(
                    seq_name = character(),
                    source = character(),
                    type = character(),
                    first = numeric(),
                    end = numeric(),
                    score = numeric(),
                    strand = character(),
                    phase = integer(),
                    gffAttributes = character()
            ),
            quiet=quiet
    )

    for ( col in c( "seq_name", "source", "type" ) )
        gff[[ col ]] <- as.factor( gff[[ col ]] )

	## ? can be used for non-available strands
	gff$strand[ gff$strand =="?"] <- NA

	## every available strand should be + or -
	if( !all( gff$strand %in% c(NA, "+", "-") ) )
		stop("Incorrect GFF file. Some entries have a strand different than '+', '-', '?' or '.'")

	## we can now force strand to be a factor with levels +, -
	gff$strand <- factor(gff$strand, levels=c("+", "-"))

    gffdf <- data.frame( gff, stringsAsFactors = FALSE )

    # in case where isRightOpen == TRUE then the convention for zero_length features
    # can be trusted. In this case, set inter_base flag for all intervals for which first==end
    # otherwise, just consider all intervals as 'base' intervals.
    if(isRightOpen){
        gffdf$inter_base <- gffdf$first == gffdf$end
    }else{
        gffdf$inter_base <- FALSE
    }

    ## make Genome_intervals_stranded object
    # number of entries
    n = nrow(gffdf)
    ints <- as.matrix(gffdf[,c( "first", "end")])
    colnames(ints)  <- NULL
    rv <- new(
            "Genome_intervals_stranded",
            ints,
            closed = matrix( rep(c(TRUE, !isRightOpen), each=n), ncol =2),
            annotation = gffdf[c("seq_name", "strand", "inter_base", "source", "type", "score", "phase", "gffAttributes" )]
    )


    val = validObject(rv, test=TRUE)
    if( val!=TRUE ){
        warning("The Genome_intervals_stranded object could be constructed from the GFF3 file. However, it is not valid according to validObject(). Check whether the GFF3 file contains invalid entries.")
        warning("validObject() failure: ", val)
    }

    return(rv)

}
