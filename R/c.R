# As of 2.8.1, we are still not able to write S4 methods for rbind(). Our
# original plan was to use combine() as in biobase, but this seems to be causing
# clashes for the generic function when both packages are loaded. Improved S4
# support for functions whose first argument is ... seems to be on the way, so
# we'll just do S3 for now, and move up to S4 once it's implemented. Note that
# as of 2.8.1, rbind.data.frame() still exists and is used for data frames.

# Update! S3 method dispatch for rbind() is non-standard (see its documentation)
# and it produced unexpected dispatch to the matrix method when presented with a
# mix of Intervals and Intervals_full objects. As a consequence, we switched to
# c(), which uses standard S3 dispatch.

# see also c.Intervals_full() in package intervals

c.Genome_intervals <- function( ... ) {
    args <- list(...)
    
    ## Drop NULL arguments
    if ( any( sapply( args, is.null ) ) )
        args <- args[ !sapply( args, is.null ) ]
    
    ## If anything else than a Genome_intervals (stranded) returns a list
    classes <- sapply( args, class )
    if ( !all( classes %in% c( "Genome_intervals", "Genome_intervals_stranded" ) ) )
        return( list( ... ) )
    
    ## If mixed Genome_intervals types (stranded and not), coerce to not stranded
    same_class <- all( classes == "Genome_intervals" )
    if ( !same_class  ) {
        warning( "Coercion to 'Genome_intervals' required.", call. = FALSE )
        return( do.call( c, lapply( args, as, "Genome_intervals" ) ) )
    }
    
    # rbinds the data frame if possible
    annot =  try( do.call( rbind, lapply( args, function(x) annotation(x) ) ), silent = TRUE )
    if( class(annot) == "try-error")
        stop("impossible to combine the annotation slots:\n", annot, call.=FALSE)
    
    # actual appending
    result <- new(
            "Genome_intervals",
            do.call( c, lapply( args, as, "Intervals_full" ) ),
            annotation = annot
    )
    return( result )
}


c.Genome_intervals_stranded <- function( ... ) {
    args <- list(...)
    
    ## Drop NULL arguments
    if ( any( sapply( args, is.null ) ) )
        args <- args[ !sapply( args, is.null ) ]
    
    ## If anything else than a Genome_intervals (stranded) returns a list
    classes <- sapply( args, class )
    if ( !all( classes %in% c( "Genome_intervals", "Genome_intervals_stranded" ) ) )
        return( list( ... ) )
    
    ## If mixed Genome_intervals types (stranded and not), coerce to not stranded
    same_class <- all( classes == "Genome_intervals_stranded" )
    if ( !same_class  ) {
        warning( "Coercion to 'Genome_intervals' required.", call. = FALSE )
        return( do.call( c, lapply( args, as, "Genome_intervals" ) ) )
    }
    
    # rbinds the data frame if possible
    annot =  try( do.call( rbind, lapply( args, function(x) annotation(x) ) ), silent = TRUE )
    if( class(annot) == "try-error")
        stop("impossible to combine the annotation slots:\n", annot)
    
    # actual appending
    result <- new(
            "Genome_intervals_stranded",
            do.call( c, lapply( args, as, "Intervals_full" ) ),
            annotation = annot
    )
    return( result )
}
