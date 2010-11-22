extract <- function(pattern, string, perl = TRUE) {
    
    # Retrieve match in first parentheses pair. This function should not be exported.
    # The part (^|;) ensures that any key is at the beginning of after a ';'
	
    r <- paste(".*(^|;)", pattern, ".*", sep = "")
    matched <- grep( r, string, perl = perl )
    result <- rep( NA, length(string) )
    result[ matched ] <- sub( r, "\\2", string[matched], perl = perl )
    return( result )
    
}

parseGffAttributes <- function( gi ) {
    
    # Unlist to a vector for names processing, then reassemble as a list with
    # split. We return a list of vectors with row indexes. If only a small
    # number of key/value pairs are required, getGffAttribute is likely faster and
    # more convenient. Note that key/value pairs which split on "=" into more than
    # 2 elements generate an error; pairs which split into a single element create
    # an NA value.
    
    if( !is( gi, "Genome_intervals" ) )
        stop( "The function requires a Genome_intervals object." )
    
    ann <- annotation(gi)
    
    pairs <- strsplit( ann$gffAttributes, ";" )
    n <- sapply( pairs, length )
    
    keyVal <- strsplit( unlist( pairs ), "=" )
    nKeyVal <- sapply( keyVal, length )
    if ( any( nKeyVal > 2 ) )
        stop( "Found more than one '=' when splitting key/value pairs." )
    if ( any( nKeyVal == 1 ) )
        warning( "One or more key/value pairs with empty value found." )
    keyVal[ nKeyVal == 1 ] <- lapply(
            keyVal[ nKeyVal == 1 ],
            function(x) c( x, NA )
    )
    
    keyVal <- unlist( keyVal )
    even <- seq( 2, length(keyVal), 2 )
    odd <- even - 1
    keyValVec <- keyVal[ even ]
    names( keyValVec ) <- keyVal[ odd ]
    
    return( split( keyValVec, rep( 1:nrow(ann), n ) ) )
    
}

getGffAttribute <- function( gi, attribute ) {
    
    ann <- annotation(gi)
    sapply(
            attribute,
            function(a) {
                # Allow key/value pairs to end at a ";" or end of string.  
                re <- paste( a, "=(.+?)(;|$)", sep = "" )
                extract( re, ann$gffAttributes )
            }
    )
    
}
