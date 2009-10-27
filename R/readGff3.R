# follows GFF3 specifications as in:
# http://song.sourceforge.net/gff3.shtml

readGff3 <- function( file, isRightOpen = TRUE ) {
    
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
            )
    )
    
    for ( col in c( "seq_name", "source", "type", "strand" ) )
        gff[[ col ]] <- as.factor( gff[[ col ]] )
    
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
