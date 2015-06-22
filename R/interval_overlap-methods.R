## adapt intervals for genome interval overlap operation

intervalsForOverlap = function(x){
    rv = as(x, "Intervals_full")

    ## close the non-empty intervals and set the empty ones to NA
    e = empty(rv)
    if( any(e) ) warning("Empty intervals encountered, set to NA.", call.=FALSE)

    rv[!e,] = close_intervals(rv[!e,])

    # set NA to empty intervals
    rv[e,] = NA

    ## set intervals type to R
    type(rv) <- "R"

    ## represent inter_base intervals as right open left open intervals (requires adding one to the right bound))
    rv[inter_base(x)][,2] <- rv[inter_base(x)][,2] + 1
    closed(rv)[inter_base(x),] <- FALSE

    return(rv)
}


setMethod(
        "interval_overlap",
        signature( "Genome_intervals", "Genome_intervals" ),
        function(
                from,
                to,
                check_valid = TRUE
        ){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
                stop( "The 'to' and/or 'from' Genome_intervals are invalid." )

            ## result value, initiated as empty list
            rv = vector( mode="list", length=nrow(from) )

            ## adapt interval representation (basically moves to R with appropriate changes for inter_base)
            qints <- intervalsForOverlap( to )
            tints <- intervalsForOverlap( from)

            ## all unique seq_names
            seqlev = unique( c( levels(seqnames(to)), levels(seqnames(from)) ) )

            ## loop over seq_names and call next method
            for(s in seqlev){
                qi = which( seqnames(to) == s)
                ti = seqnames(from) == s

                rv[ti] = lapply(
                        interval_overlap(tints[ti], qints[qi], check_valid=check_valid),
                        function(v) qi[v]
                )
            }

            return(rv)
        }
)


setMethod(
        "interval_overlap",
        signature( "Genome_intervals_stranded", "Genome_intervals_stranded" ),
        function(
                from,
                to,
                check_valid = TRUE
        ){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
                stop( "The 'to' and/or 'from' Genome_intervals are invalid." )

            if( any( is.na(strand(to)) ) || any( is.na(strand(from)) ) )
                stop("NA(s) present in the strand of 'to' or 'from'.")
            ## result value, initiated as empty list
            rv = vector( mode="list", length=nrow(from) )
            strdlev = levels( strand(from) )
            ## loop over strands and call next method
            nextMethod = getMethod("interval_overlap", signature( "Genome_intervals",  "Genome_intervals" ) )
            for(s in strdlev){
                qi = which(strand(to) == s)
                ti = strand(from) == s

                rv[ti] =  lapply(
                        nextMethod( from[ti], to[qi] ),
                        function(v) qi[v]
                )
            }
            return(rv)
        }
)


argument_error <- paste(
        "The 'from' and 'to' arguments are required. Note that the",
        "  interval_overlap argument names changed at v. 1.1.1.",
        "  See documentation.",
        sep = "\n"
)

setMethod(
        "interval_overlap",
        signature( from = "missing", to = "ANY" ),
        function( from, to, check_valid, ... ) stop( argument_error )
)

setMethod(
        "interval_overlap",
        signature( from = "ANY", to = "missing" ),
        function( from, to, check_valid, ... ) stop( argument_error )
)

